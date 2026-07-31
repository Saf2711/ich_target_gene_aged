from __future__ import annotations

from pathlib import Path

import pandas as pd

from .aging import (
    find_region_files,
    load_region_tables,
    save_age_outputs,
    summarize_age_effects,
)
from .config import (
    DEFAULT_ICH_QUERY,
    GENCODE_GTF_URL,
    NCBI_GENE_INFO_URL,
    ScoreWeights,
    TierThresholds,
)
from .download import download_cns_age_tables
from .figures import (
    plot_age_tier_counts,
    plot_gene_regional_age_effect,
    plot_ich_evidence_distribution,
    plot_top_priority_targets,
)
from .literature import (
    ApiSettings,
    NcbiClient,
    run_literature_screen,
)
from .payload import (
    add_payload_annotations,
    parse_gencode_cds_lengths,
)
from .scoring import assign_priority_scores


def run_aging_analysis(
    input_dir: Path,
    output_dir: Path,
    thresholds: TierThresholds | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the complete genome-wide CNS age analysis."""

    files = find_region_files(input_dir)

    if len(files) != 13:
        raise RuntimeError(
            f"Expected 13 CNS-region files; found {len(files)} in {input_dir}."
        )

    long_data = load_region_tables(files)
    summary = summarize_age_effects(
        long_data,
        thresholds=thresholds,
    )

    save_age_outputs(
        long_data=long_data,
        summary=summary,
        output_dir=output_dir,
    )

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    plot_age_tier_counts(
        summary,
        figures_dir / "age_tier_counts.png",
    )

    plot_gene_regional_age_effect(
        long_data=long_data,
        gene="GPX4",
        output_path=(
            figures_dir
            / "gpx4_age_effect_by_region.png"
        ),
    )

    return long_data, summary


def run_target_prioritization(
    downward_genes_path: Path,
    output_dir: Path,
    email: str,
    api_key: str | None = None,
    query: str = DEFAULT_ICH_QUERY,
    gene_info_url: str = NCBI_GENE_INFO_URL,
    gencode_url: str = GENCODE_GTF_URL,
    weights: ScoreWeights | None = None,
    top_n: int = 20,
) -> pd.DataFrame:
    """Run literature screening, payload annotation, and final ranking."""

    output_dir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(downward_genes_path)
    genes["Gene"] = (
        genes["Gene"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    screened = run_literature_screen(
        downward_genes=genes,
        output_dir=output_dir / "literature",
        email=email,
        api_key=api_key,
        query=query,
        gene_info_url=gene_info_url,
    )

    payload_dir = output_dir / "payload"
    payload_dir.mkdir(parents=True, exist_ok=True)

    settings = ApiSettings(
        email=email,
        api_key=api_key,
    )

    client = NcbiClient(
        settings=settings,
        cache_dir=output_dir / "cache",
    )

    gtf_path = (
        payload_dir
        / Path(gencode_url).name
    )

    client.download(
        gencode_url,
        gtf_path,
    )

    gencode = parse_gencode_cds_lengths(
        gtf_path
    )

    gencode.to_csv(
        payload_dir
        / "gencode_selected_transcripts.csv",
        index=False,
    )

    annotated = add_payload_annotations(
        screened,
        gencode,
    )

    ranking = assign_priority_scores(
        annotated,
        weights=weights,
    )

    ranking.to_csv(
        output_dir
        / "all_age_declining_genes_ranked.csv",
        index=False,
    )

    ranking.loc[
        ranking[
            "auto_direction_gate"
        ].eq(
            "restoration_supported"
        )
    ].to_csv(
        output_dir
        / "restoration_supported_targets.csv",
        index=False,
    )

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    plot_ich_evidence_distribution(
        ranking,
        figures_dir
        / "ich_evidence_distribution.png",
    )

    plot_top_priority_targets(
        ranking,
        figures_dir
        / "top_priority_targets.png",
        top_n=top_n,
    )

    return ranking


def run_all(
    project_dir: Path,
    email: str,
    api_key: str | None = None,
    query: str = DEFAULT_ICH_QUERY,
    weights: ScoreWeights | None = None,
    top_n: int = 20,
) -> pd.DataFrame:
    """Download the aging tables and run the complete pipeline."""

    raw_age_dir = (
        project_dir
        / "data"
        / "raw"
        / "age_effect_explorer"
    )

    age_results_dir = (
        project_dir
        / "results"
        / "aging"
    )

    target_results_dir = (
        project_dir
        / "results"
        / "target_prioritization"
    )

    download_cns_age_tables(
        output_dir=raw_age_dir,
    )

    run_aging_analysis(
        input_dir=raw_age_dir,
        output_dir=age_results_dir,
    )

    return run_target_prioritization(
        downward_genes_path=(
            age_results_dir
            / "all_downward_genes.csv"
        ),
        output_dir=target_results_dir,
        email=email,
        api_key=api_key,
        query=query,
        weights=weights,
        top_n=top_n,
    )
