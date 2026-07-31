# Genome-wide aging analysis and ICH gene-delivery target prioritization

This repository contains the complete Python workflow used to:

1. download the Age Effect Explorer CNS age-regression tables;
2. perform a genome-wide aging analysis across 13 human CNS regions;
3. identify genes with reproducible age-related decline;
4. screen the complete age-declining gene set against the intracerebral
   hemorrhage literature;
5. annotate coding-sequence size using GENCODE;
6. calculate age, ICH-evidence, payload, and final priority scores;
7. export ranked tables and publication-quality figures.

## Analysis overview

The aging analysis starts with all genes represented in the 13 CNS tables. It
does not begin with GPX4 or a predefined candidate list.

The processed dataset in this repository contains:

| Metric | Count |
|---|---:|
| CNS regions | 13 |
| Unique gene labels | 28,250 |
| High-confidence age-declining genes | 153 |
| Moderate age-decline genes | 1,927 |
| Directional age-decline genes | 7,738 |
| Total downward-gene universe | 9,818 |

The target-prioritization stage intersects the full downward-gene universe with
an ICH PubMed corpus, adds PubTator gene annotations, assigns article-derived
ICH evidence scores, adds GENCODE coding-sequence lengths, and calculates the
final priority score.

## Repository structure

```text
.
├── data
│   ├── processed
│   │   ├── genomewide_brain_aging_gene_screen.csv
│   │   └── all_downward_genes.csv
│   └── raw
├── docs
│   ├── DATA_DICTIONARY.md
│   ├── METHODS.md
│   ├── REPOSITORY_STRUCTURE.md
│   └── SCORING.md
├── results
├── src
│   └── aged_ich
│       ├── aging.py
│       ├── cli.py
│       ├── config.py
│       ├── download.py
│       ├── figures.py
│       ├── literature.py
│       ├── payload.py
│       ├── pipeline.py
│       └── scoring.py
├── tests
├── pyproject.toml
├── requirements.txt
└── run_pipeline.py
```

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
```

macOS/Linux:

```bash
source .venv/bin/activate
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the package:

```bash
pip install -e .
```

## Complete start-to-finish run

NCBI E-utilities requires an email address in API requests.

```bash
aged-ich run-all \
  --project-dir . \
  --email your_email@example.com
```

With an NCBI API key:

```bash
aged-ich run-all \
  --project-dir . \
  --email your_email@example.com \
  --api-key "$NCBI_API_KEY"
```

## Run each stage separately

### 1. Download the 13 CNS age tables

```bash
aged-ich download-aging-data \
  --output-dir data/raw/age_effect_explorer
```

### 2. Analyze the aging dataset

```bash
aged-ich analyze-aging \
  --input-dir data/raw/age_effect_explorer \
  --output-dir results/aging
```

### 3. Prioritize ICH gene-delivery targets

```bash
aged-ich prioritize-targets \
  --downward-genes results/aging/all_downward_genes.csv \
  --email your_email@example.com \
  --output-dir results/target_prioritization
```

## Final outputs

The aging stage writes:

```text
regional_age_effects_long.csv
genomewide_brain_aging_gene_screen.csv
high_confidence_age_decline_genes.csv
moderate_age_decline_genes.csv
directional_age_decline_genes.csv
all_downward_genes.csv
figures/age_tier_counts.png
figures/gpx4_age_effect_by_region.png
```

The target-prioritization stage writes:

```text
literature/gene_article_evidence.csv
literature/all_downward_genes_literature_screened.csv
payload/gencode_selected_transcripts.csv
all_age_declining_genes_ranked.csv
restoration_supported_targets.csv
figures/ich_evidence_distribution.png
figures/top_priority_targets.png
```

## Default scoring weights

```text
Age-decline evidence: 20%
ICH literature evidence: 60%
Payload suitability: 20%
```

The weights can be changed from the command line:

```bash
aged-ich prioritize-targets \
  --downward-genes results/aging/all_downward_genes.csv \
  --email your_email@example.com \
  --age-weight 0.20 \
  --ich-weight 0.60 \
  --payload-weight 0.20
```

The three values must sum to 1.0.

## Processed input tables

The repository includes the previously generated genome-wide summary and the
complete 9,818-gene downward set under `data/processed/`. These files allow the
literature and payload stages to be run without repeating the aging analysis.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

## Documentation

- [`docs/METHODS.md`](docs/METHODS.md): full computational method
- [`docs/SCORING.md`](docs/SCORING.md): tier and score definitions
- [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md): output columns
- [`docs/REPOSITORY_STRUCTURE.md`](docs/REPOSITORY_STRUCTURE.md): file guide
