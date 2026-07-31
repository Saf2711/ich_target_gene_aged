from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_age_tier_counts(
    summary: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the number of genes in each age-evidence tier."""

    order = [
        "High-confidence age decline",
        "Moderate age-decline evidence",
        "Directional age decline",
        "Not prioritized",
    ]

    counts = (
        summary["discovery_tier"]
        .value_counts()
        .reindex(order, fill_value=0)
    )

    plt.figure(figsize=(9, 5))
    plt.bar(
        counts.index,
        counts.values,
    )
    plt.ylabel("Number of genes")
    plt.xlabel("Age-evidence tier")
    plt.title("Genome-wide CNS age-decline classification")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def plot_gene_regional_age_effect(
    long_data: pd.DataFrame,
    gene: str,
    output_path: Path,
) -> None:
    """Plot one gene's estimated 40-year effect across CNS regions."""

    selected = long_data.loc[
        long_data["Gene"].eq(gene.upper())
    ].copy()

    if selected.empty:
        return

    selected["estimated_pct_change_40y"] = (
        np.power(
            2.0,
            selected["age_coef"] * 40.0,
        )
        - 1.0
    ) * 100.0

    selected = selected.sort_values(
        "estimated_pct_change_40y"
    )

    plt.figure(figsize=(10, 6))
    plt.barh(
        selected["region"],
        selected["estimated_pct_change_40y"],
    )
    plt.axvline(0, linewidth=1)
    plt.xlabel(
        "Estimated expression change over 40 years (%)"
    )
    plt.ylabel("GTEx CNS region")
    plt.title(
        f"{gene.upper()} age association across GTEx CNS regions"
    )
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def plot_ich_evidence_distribution(
    ranking: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the number of genes assigned to each automated ICH score."""

    counts = (
        ranking["auto_ich_score"]
        .value_counts()
        .reindex(
            range(6),
            fill_value=0,
        )
        .sort_index()
    )

    plt.figure(figsize=(8, 5))
    plt.bar(
        counts.index.astype(str),
        counts.values,
    )
    plt.xlabel("Automated ICH evidence score")
    plt.ylabel("Number of age-declining genes")
    plt.title(
        "ICH literature evidence across age-declining genes"
    )
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def plot_top_priority_targets(
    ranking: pd.DataFrame,
    output_path: Path,
    top_n: int = 20,
) -> None:
    """Plot the highest-ranked restoration-compatible targets."""

    plotted = (
        ranking.loc[
            ~ranking[
                "auto_direction_gate"
            ].eq(
                "restoration_contraindicated"
            )
        ]
        .head(top_n)
        .sort_values("priority_score")
    )

    if plotted.empty:
        return

    plt.figure(
        figsize=(
            10,
            max(
                6,
                top_n * 0.30,
            ),
        )
    )

    plt.barh(
        plotted["Gene"],
        plotted["priority_score"],
    )

    plt.xlabel("Priority score (0-100)")
    plt.ylabel("Gene")
    plt.title(
        "Top age-declining ICH gene-delivery targets"
    )
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
