from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def parse_gtf_attributes(
    value: str,
) -> dict[str, list[str]]:
    """Parse the ninth column of a GENCODE GTF file."""

    attributes: dict[str, list[str]] = defaultdict(list)

    for item in value.strip().strip(";").split(";"):
        item = item.strip()

        if not item:
            continue

        key, _, raw_value = item.partition(" ")

        attributes[key].append(
            raw_value.strip().strip('"')
        )

    return dict(attributes)


def parse_gencode_cds_lengths(
    path: Path,
) -> pd.DataFrame:
    """
    Calculate transcript CDS lengths from a GENCODE GTF.

    MANE Select transcripts are preferred. When no MANE Select transcript is
    available, the longest basic protein-coding transcript is selected.
    """

    opener = gzip.open if path.suffix == ".gz" else open

    transcript_lengths: dict[
        tuple[str, str],
        int,
    ] = defaultdict(int)

    transcript_metadata: dict[
        tuple[str, str],
        dict[str, object],
    ] = {}

    with opener(
        path,
        "rt",
        encoding="utf-8",
    ) as handle:
        for line in tqdm(
            handle,
            desc="Parsing GENCODE CDS features",
        ):
            if (
                not line
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip("\n").split("\t")

            if (
                len(fields) != 9
                or fields[2] != "CDS"
            ):
                continue

            start = int(fields[3])
            end = int(fields[4])

            attrs = parse_gtf_attributes(
                fields[8]
            )

            gene_name = (
                attrs.get("gene_name")
                or [""]
            )[0].upper()

            transcript_id = (
                attrs.get("transcript_id")
                or [""]
            )[0]

            gene_type = (
                attrs.get("gene_type")
                or [""]
            )[0]

            transcript_type = (
                attrs.get("transcript_type")
                or [""]
            )[0]

            tags = set(
                attrs.get("tag", [])
            )

            if (
                not gene_name
                or not transcript_id
            ):
                continue

            key = (
                gene_name,
                transcript_id,
            )

            transcript_lengths[key] += (
                end - start + 1
            )

            transcript_metadata[key] = {
                "gene": gene_name,
                "transcript_id": transcript_id,
                "gene_type_gencode": gene_type,
                "transcript_type": transcript_type,
                "is_mane_select": (
                    "MANE_Select" in tags
                ),
                "is_basic": "basic" in tags,
            }

    rows: list[dict[str, object]] = []

    for key, cds_length in (
        transcript_lengths.items()
    ):
        metadata = transcript_metadata[key]

        rows.append(
            {
                **metadata,
                "cds_bp": cds_length,
            }
        )

    transcripts = pd.DataFrame(rows)

    if transcripts.empty:
        raise ValueError(
            "No CDS features were parsed from the GENCODE GTF."
        )

    transcripts["_priority"] = (
        transcripts[
            "is_mane_select"
        ].astype(int) * 2
        + transcripts[
            "is_basic"
        ].astype(int)
    )

    selected = (
        transcripts.sort_values(
            [
                "gene",
                "_priority",
                "cds_bp",
            ],
            ascending=[
                True,
                False,
                False,
            ],
        )
        .drop_duplicates("gene")
        .drop(columns="_priority")
    )

    return selected


def assign_payload_score(
    cds_bp: float,
) -> int:
    """Assign the coding-sequence payload score."""

    if pd.isna(cds_bp):
        return 0

    if cds_bp <= 800:
        return 5

    if cds_bp <= 1500:
        return 4

    if cds_bp <= 2250:
        return 3

    if cds_bp <= 3000:
        return 2

    if cds_bp <= 4000:
        return 1

    return 0


def add_payload_annotations(
    genes: pd.DataFrame,
    gencode: pd.DataFrame,
) -> pd.DataFrame:
    """Join the selected GENCODE transcript and payload score."""

    result = genes.merge(
        gencode,
        left_on="Gene",
        right_on="gene",
        how="left",
    )

    result["payload_score"] = (
        result["cds_bp"]
        .map(assign_payload_score)
        .astype(int)
    )

    result["protein_coding"] = (
        result["gene_type_gencode"]
        .fillna("")
        .eq("protein_coding")
    )

    return result.drop(
        columns=["gene"],
        errors="ignore",
    )
