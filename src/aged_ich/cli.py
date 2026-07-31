from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests

from .aging import (
    find_region_files,
    load_region_tables,
    save_age_outputs,
    summarize_age_effects,
)
from .config import (
    DEFAULT_ICH_QUERY,
    ScoreWeights,
    TierThresholds,
)
from .download import download_cns_age_tables
from .pipeline import (
    run_aging_analysis,
    run_all,
    run_target_prioritization,
)


def add_weight_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    parser.add_argument(
        "--age-weight",
        type=float,
        default=0.20,
    )
    parser.add_argument(
        "--ich-weight",
        type=float,
        default=0.60,
    )
    parser.add_argument(
        "--payload-weight",
        type=float,
        default=0.20,
    )


def weights_from_args(
    args: argparse.Namespace,
) -> ScoreWeights:
    weights = ScoreWeights(
        age=args.age_weight,
        ich=args.ich_weight,
        payload=args.payload_weight,
    )
    weights.validate()
    return weights


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aged-ich",
        description=(
            "Genome-wide aging analysis and ICH "
            "gene-delivery target prioritization."
        ),
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    download = subparsers.add_parser(
        "download-aging-data",
        help=(
            "Download the 13 Age Effect Explorer "
            "CNS age-regression tables."
        ),
    )
    download.add_argument(
        "--output-dir",
        type=Path,
        default=(
            Path("data")
            / "raw"
            / "age_effect_explorer"
        ),
    )

    analyze = subparsers.add_parser(
        "analyze-aging",
        help=(
            "Build the genome-wide CNS aging summary."
        ),
    )
    analyze.add_argument(
        "--input-dir",
        type=Path,
        required=True,
    )
    analyze.add_argument(
        "--output-dir",
        type=Path,
        default=(
            Path("results")
            / "aging"
        ),
    )

    prioritize = subparsers.add_parser(
        "prioritize-targets",
        help=(
            "Run the ICH literature, payload, "
            "and ranking stages."
        ),
    )
    prioritize.add_argument(
        "--downward-genes",
        type=Path,
        required=True,
    )
    prioritize.add_argument(
        "--output-dir",
        type=Path,
        default=(
            Path("results")
            / "target_prioritization"
        ),
    )
    prioritize.add_argument(
        "--email",
        required=True,
    )
    prioritize.add_argument(
        "--api-key",
        default=None,
    )
    prioritize.add_argument(
        "--ich-query",
        default=DEFAULT_ICH_QUERY,
    )
    prioritize.add_argument(
        "--top-n",
        type=int,
        default=20,
    )
    add_weight_arguments(prioritize)

    run = subparsers.add_parser(
        "run-all",
        help=(
            "Download data and run the complete "
            "pipeline from start to finish."
        ),
    )
    run.add_argument(
        "--project-dir",
        type=Path,
        default=Path("."),
    )
    run.add_argument(
        "--email",
        required=True,
    )
    run.add_argument(
        "--api-key",
        default=None,
    )
    run.add_argument(
        "--ich-query",
        default=DEFAULT_ICH_QUERY,
    )
    run.add_argument(
        "--top-n",
        type=int,
        default=20,
    )
    add_weight_arguments(run)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "download-aging-data":
        files = download_cns_age_tables(
            output_dir=args.output_dir,
        )
        print(
            f"Downloaded {len(files)} CNS tables "
            f"to {args.output_dir}."
        )
        return 0

    if args.command == "analyze-aging":
        _, summary = run_aging_analysis(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            thresholds=TierThresholds(),
        )

        print(
            summary[
                "discovery_tier"
            ].value_counts().to_string()
        )
        return 0

    if args.command == "prioritize-targets":
        ranking = run_target_prioritization(
            downward_genes_path=(
                args.downward_genes
            ),
            output_dir=args.output_dir,
            email=args.email,
            api_key=args.api_key,
            query=args.ich_query,
            weights=weights_from_args(args),
            top_n=args.top_n,
        )

        print(
            ranking[
                [
                    "screen_rank",
                    "Gene",
                    "age_score",
                    "auto_ich_score",
                    "payload_score",
                    "priority_score",
                ]
            ]
            .head(args.top_n)
            .to_string(index=False)
        )
        return 0

    if args.command == "run-all":
        ranking = run_all(
            project_dir=args.project_dir,
            email=args.email,
            api_key=args.api_key,
            query=args.ich_query,
            weights=weights_from_args(args),
            top_n=args.top_n,
        )

        print(
            ranking[
                [
                    "screen_rank",
                    "Gene",
                    "age_score",
                    "auto_ich_score",
                    "payload_score",
                    "priority_score",
                ]
            ]
            .head(args.top_n)
            .to_string(index=False)
        )
        return 0

    parser.error(
        f"Unknown command: {args.command}"
    )
    return 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())

    except (
        FileNotFoundError,
        RuntimeError,
        ValueError,
        requests.RequestException,
    ) as error:
        print(
            f"Error: {error}",
            file=sys.stderr,
        )
        raise SystemExit(1)
