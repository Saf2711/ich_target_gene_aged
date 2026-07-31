# Repository structure

## `src/aged_ich/download.py`

Downloads the 13 CNS age-regression tables from the Age Effect Explorer GitHub
repository.

## `src/aged_ich/aging.py`

Loads the regional tables, calculates regional and cross-region metrics, assigns
age-evidence tiers, and writes the genome-wide outputs.

## `src/aged_ich/literature.py`

Builds the ICH PubMed corpus, retrieves PubTator annotations, maps NCBI Gene IDs
to human symbols, classifies title/abstract evidence, and aggregates gene-level
ICH scores.

## `src/aged_ich/payload.py`

Parses the GENCODE GTF file, selects one transcript per gene, calculates CDS
length, and assigns payload scores.

## `src/aged_ich/scoring.py`

Combines the age, ICH, and payload scores and generates the final rank.

## `src/aged_ich/figures.py`

Creates the age-tier, regional age-effect, evidence-distribution, and priority
figures.

## `src/aged_ich/pipeline.py`

Provides the reusable start-to-finish workflow functions.

## `src/aged_ich/cli.py`

Implements the `aged-ich` command-line interface.

## `data/processed`

Contains the processed genome-wide age summary and the full downward-gene
universe generated during the original analysis.

## `results`

Default location for newly generated output tables and figures.
