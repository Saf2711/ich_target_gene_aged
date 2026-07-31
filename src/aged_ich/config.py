from __future__ import annotations

from dataclasses import dataclass


AGE_EFFECT_GITHUB_REPOSITORIES = (
    "ML198/GTEx-Explorer",
    "ML198/Age-Effect-Explorer",
)

AGE_EFFECT_DATA_DIRECTORY = "data/p_value"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PUBTATOR_EXPORT = (
    "https://www.ncbi.nlm.nih.gov/research/"
    "pubtator-api/publications/export/biocjson"
)

NCBI_GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
    "Mammalia/Homo_sapiens.gene_info.gz"
)

GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_50/gencode.v50.annotation.gtf.gz"
)

DEFAULT_ICH_QUERY = (
    '("intracerebral hemorrhage"[Title/Abstract] '
    'OR "intracerebral haemorrhage"[Title/Abstract] '
    'OR "hemorrhagic stroke"[Title/Abstract] '
    'OR "haemorrhagic stroke"[Title/Abstract] '
    'OR "brain hemorrhage"[Title/Abstract] '
    'OR "brain haemorrhage"[Title/Abstract])'
)

AGE_SCORE_MAP = {
    "High-confidence age decline": 5,
    "Moderate age-decline evidence": 4,
    "Directional age decline": 3,
    "Not prioritized": 0,
}


@dataclass(frozen=True)
class TierThresholds:
    """Rules for assigning genome-wide age-decline evidence tiers."""

    high_min_regions: int = 13
    high_min_negative_regions: int = 12
    high_min_fdr_regions: int = 2
    high_max_effect_percentile: float = 0.10

    moderate_min_regions: int = 10
    moderate_min_negative_fraction: float = 0.80
    moderate_min_fdr_regions: int = 1
    moderate_max_effect_percentile: float = 0.25

    directional_min_regions: int = 10
    directional_min_negative_fraction: float = 0.70


@dataclass(frozen=True)
class ScoreWeights:
    """Weights used for the final 0-100 priority score."""

    age: float = 0.20
    ich: float = 0.60
    payload: float = 0.20

    def validate(self) -> None:
        total = self.age + self.ich + self.payload
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                f"Age, ICH, and payload weights must sum to 1.0; received {total:.3f}."
            )
