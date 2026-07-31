from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .config import TierThresholds


REQUIRED_REGION_COLUMNS = {
    "Gene",
    "age_coef",
    "p_value_age",
    "BH_adjusted_age",
}


def normalize_gene_symbol(series: pd.Series) -> pd.Series:
    """Standardize gene labels before joining or grouping."""

    return series.astype("string").str.strip().str.upper()


def region_name_from_path(path: Path) -> str:
    """Convert a regional filename to a readable tissue label."""

    name = path.stem

    for suffix in ("_full", "_pvalue_results"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]

    if name.startswith("brain_"):
        name = name[len("brain_") :]

    return name.replace("_", " ")


def find_region_files(input_dir: Path) -> list[Path]:
    """Locate regional CNS CSV files in an input directory."""

    patterns = (
        "brain_*_full.csv",
        "brain_*_pvalue_results.csv",
    )
    files: list[Path] = []

    for pattern in patterns:
        files.extend(input_dir.glob(pattern))

    return sorted(set(files))


def load_region_tables(files: Iterable[Path]) -> pd.DataFrame:
    """Read and combine the regional age-regression tables."""

    frames: list[pd.DataFrame] = []

    for path in files:
        table = pd.read_csv(path)

        missing = REQUIRED_REGION_COLUMNS.difference(table.columns)
        if missing:
            raise ValueError(
                f"{path.name} is missing required columns: "
                f"{', '.join(sorted(missing))}"
            )

        table = table.copy()
        table["Gene"] = normalize_gene_symbol(table["Gene"])
        table["region"] = region_name_from_path(path)

        for column in ("age_coef", "p_value_age", "BH_adjusted_age"):
            table[column] = pd.to_numeric(table[column], errors="coerce")

        table = table.dropna(
            subset=[
                "Gene",
                "age_coef",
                "p_value_age",
                "BH_adjusted_age",
            ]
        )
        frames.append(table)

    if not frames:
        raise FileNotFoundError("No regional CNS CSV files were found.")

    return pd.concat(frames, ignore_index=True)


def summarize_age_effects(
    long_data: pd.DataFrame,
    thresholds: TierThresholds | None = None,
) -> pd.DataFrame:
    """
    Aggregate age effects for each gene across all CNS regions.

    The age coefficient is assumed to represent the annual change in
    log2(TPM + 1). The 40-year percentage change is calculated as:

        [2 ** (age_coef * 40) - 1] * 100
    """

    thresholds = thresholds or TierThresholds()
    data = long_data.copy()

    data["regional_decline_percentile"] = (
        data.groupby("region")["age_coef"]
        .rank(method="average", ascending=True, pct=True)
    )

    data["estimated_pct_change_40y"] = (
        np.power(2.0, data["age_coef"] * 40.0) - 1.0
    ) * 100.0

    strongest_index = data.groupby("Gene")["age_coef"].idxmin()
    strongest = data.loc[
        strongest_index,
        [
            "Gene",
            "region",
            "age_coef",
            "BH_adjusted_age",
        ],
    ].rename(
        columns={
            "region": "strongest_decline_region",
            "age_coef": "strongest_region_age_coef",
            "BH_adjusted_age": "strongest_region_fdr",
        }
    )

    summary = (
        data.groupby("Gene")
        .agg(
            regions_tested=("region", "nunique"),
            negative_regions=(
                "age_coef",
                lambda values: int((values < 0).sum()),
            ),
            median_age_coef=("age_coef", "median"),
            mean_age_coef=("age_coef", "mean"),
            median_pct_change_40y=(
                "estimated_pct_change_40y",
                "median",
            ),
            mean_pct_change_40y=(
                "estimated_pct_change_40y",
                "mean",
            ),
            median_regional_decline_percentile=(
                "regional_decline_percentile",
                "median",
            ),
            best_region_fdr=("BH_adjusted_age", "min"),
            median_nominal_p=("p_value_age", "median"),
        )
        .reset_index()
    )

    significant_negative = (
        data.assign(
            significant_negative=(
                (data["age_coef"] < 0)
                & (data["BH_adjusted_age"] < 0.05)
            )
        )
        .groupby("Gene")["significant_negative"]
        .sum()
        .astype(int)
        .rename("fdr_sig_negative_regions")
    )

    summary = (
        summary.merge(significant_negative, on="Gene", how="left")
        .merge(strongest, on="Gene", how="left")
    )

    summary["negative_fraction"] = (
        summary["negative_regions"] / summary["regions_tested"]
    )

    covered = summary["regions_tested"] >= thresholds.moderate_min_regions

    summary["median_effect_percentile_among_covered"] = np.nan
    summary.loc[covered, "median_effect_percentile_among_covered"] = (
        summary.loc[covered, "median_age_coef"]
        .rank(method="average", ascending=True, pct=True)
    )

    high = (
        (summary["regions_tested"] >= thresholds.high_min_regions)
        & (
            summary["negative_regions"]
            >= thresholds.high_min_negative_regions
        )
        & (
            summary["fdr_sig_negative_regions"]
            >= thresholds.high_min_fdr_regions
        )
        & (
            summary["median_effect_percentile_among_covered"]
            <= thresholds.high_max_effect_percentile
        )
    )

    moderate = (
        (
            summary["regions_tested"]
            >= thresholds.moderate_min_regions
        )
        & (
            summary["negative_fraction"]
            >= thresholds.moderate_min_negative_fraction
        )
        & (
            summary["fdr_sig_negative_regions"]
            >= thresholds.moderate_min_fdr_regions
        )
        & (
            summary["median_effect_percentile_among_covered"]
            <= thresholds.moderate_max_effect_percentile
        )
        & ~high
    )

    directional = (
        (
            summary["regions_tested"]
            >= thresholds.directional_min_regions
        )
        & (
            summary["negative_fraction"]
            >= thresholds.directional_min_negative_fraction
        )
        & (summary["median_age_coef"] < 0)
        & ~high
        & ~moderate
    )

    summary["discovery_tier"] = np.select(
        [high, moderate, directional],
        [
            "High-confidence age decline",
            "Moderate age-decline evidence",
            "Directional age decline",
        ],
        default="Not prioritized",
    )

    tier_order = {
        "High-confidence age decline": 0,
        "Moderate age-decline evidence": 1,
        "Directional age decline": 2,
        "Not prioritized": 3,
    }

    summary["_tier_order"] = summary["discovery_tier"].map(
        tier_order
    )

    summary = summary.sort_values(
        [
            "_tier_order",
            "fdr_sig_negative_regions",
            "negative_fraction",
            "median_effect_percentile_among_covered",
            "median_age_coef",
            "Gene",
        ],
        ascending=[
            True,
            False,
            False,
            True,
            True,
            True,
        ],
        na_position="last",
    ).reset_index(drop=True)

    summary.insert(
        0,
        "objective_rank",
        np.arange(1, len(summary) + 1),
    )

    return summary.drop(columns="_tier_order")


def save_age_outputs(
    long_data: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Path]:
    """Write the complete and tier-specific age-analysis tables."""

    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "regional_long": output_dir / "regional_age_effects_long.csv",
        "genomewide_summary": (
            output_dir / "genomewide_brain_aging_gene_screen.csv"
        ),
        "high_confidence": (
            output_dir / "high_confidence_age_decline_genes.csv"
        ),
        "moderate": output_dir / "moderate_age_decline_genes.csv",
        "directional": (
            output_dir / "directional_age_decline_genes.csv"
        ),
        "all_downward": output_dir / "all_downward_genes.csv",
    }

    long_data.to_csv(outputs["regional_long"], index=False)
    summary.to_csv(outputs["genomewide_summary"], index=False)

    tier_to_path = {
        "High-confidence age decline": outputs["high_confidence"],
        "Moderate age-decline evidence": outputs["moderate"],
        "Directional age decline": outputs["directional"],
    }

    for tier, path in tier_to_path.items():
        summary.loc[summary["discovery_tier"] == tier].to_csv(
            path,
            index=False,
        )

    summary.loc[
        summary["discovery_tier"] != "Not prioritized"
    ].to_csv(
        outputs["all_downward"],
        index=False,
    )

    return outputs
