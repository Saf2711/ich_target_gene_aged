# Computational methods

## 1. Aging dataset

The aging analysis uses tissue-specific age-regression tables from Age Effect
Explorer, which were generated from GTEx v10 gene-expression data. Thirteen
human CNS regions are included:

1. amygdala;
2. anterior cingulate cortex, BA24;
3. caudate, basal ganglia;
4. cerebellar hemisphere;
5. cerebellum;
6. cortex;
7. frontal cortex, BA9;
8. hippocampus;
9. hypothalamus;
10. nucleus accumbens, basal ganglia;
11. putamen, basal ganglia;
12. spinal cord, cervical C1;
13. substantia nigra.

Each regional file contains a gene symbol, age coefficient, nominal age
P-value, and Benjamini-Hochberg-adjusted age P-value.

The age coefficient represents annual change in log2(TPM + 1).

## 2. Regional table processing

The script:

1. loads all 13 regional CSV files;
2. standardizes gene symbols to uppercase;
3. converts the statistical columns to numeric values;
4. removes rows missing the gene or age statistics;
5. appends the CNS region name;
6. concatenates the regional tables into one long dataframe.

## 3. Effect-size transformation

For each regional gene result, the age coefficient is converted to an
estimated expression difference over 40 years:

```text
estimated 40-year change (%) =
[2^(age coefficient × 40) − 1] × 100
```

A negative result indicates lower estimated expression at the older end of the
40-year interval.

## 4. Genome-wide aggregation

For every gene, the script calculates:

- number of CNS regions tested;
- number of regions with a negative age coefficient;
- fraction of tested regions showing decline;
- number of FDR-significant negative regions;
- median and mean age coefficients;
- median and mean estimated 40-year changes;
- minimum regional FDR;
- median nominal P-value;
- region with the strongest negative coefficient;
- median within-region decline percentile.

The within-region percentile ranks all genes from the strongest decline to the
strongest increase within each CNS tissue.

## 5. Age-decline evidence tiers

### High-confidence age decline

A gene must:

- be measured in all 13 regions;
- decline in at least 12 regions;
- have at least two FDR-significant regional declines;
- have a median age coefficient within the strongest 10% of genes measured in
  at least 10 regions.

### Moderate age-decline evidence

A gene must:

- be measured in at least 10 regions;
- decline in at least 80% of measured regions;
- have at least one FDR-significant regional decline;
- have a median age coefficient within the strongest 25% of genes measured in
  at least 10 regions.

### Directional age decline

A gene must:

- be measured in at least 10 regions;
- decline in at least 70% of measured regions;
- have a negative median age coefficient;
- not meet the high-confidence or moderate criteria.

The complete downward-gene universe is the union of the high-confidence,
moderate, and directional groups.

## 6. ICH publication corpus

The literature stage runs one PubMed search using the following title/abstract
terms:

```text
intracerebral hemorrhage
intracerebral haemorrhage
hemorrhagic stroke
haemorrhagic stroke
brain hemorrhage
brain haemorrhage
```

All matching PubMed identifiers are retrieved through NCBI E-utilities.

The associated titles, abstracts, publication types, years, and DOIs are
downloaded in batches and stored locally.

## 7. Gene normalization

PubTator Central annotations are downloaded for the ICH publication corpus.

PubTator gene mentions are normalized to NCBI Gene identifiers. Human NCBI Gene
identifiers are converted to official gene symbols using
`Homo_sapiens.gene_info.gz`.

The normalized ICH gene symbols are intersected with the full age-declining gene
universe.

## 8. Automated ICH evidence classification

Titles and abstracts are evaluated for:

- increased expression or protein delivery;
- decreased expression or inhibition;
- beneficial outcomes;
- harmful outcomes;
- mechanistic evidence;
- primary experimental models or human studies;
- publication type.

Each gene-publication pair is assigned one of the following classes:

```text
increase_beneficial
decrease_harmful
decrease_beneficial
increase_harmful
mechanistic_primary
mention_only
review
```

Gene-level counts are calculated for each class.

## 9. ICH evidence score

The automated ICH score is assigned from the gene-level publication counts:

| Score | Rule |
|---:|---|
| 5 | At least two restorative-direction primary articles |
| 4 | One restorative-direction primary article |
| 3 | Mechanistic primary evidence or a contraindicating-direction intervention |
| 2 | Primary article mention without a stronger intervention class |
| 1 | Review-only evidence |
| 0 | No annotated ICH article detected |

The direction gate is:

```text
restoration_supported
restoration_contraindicated
conflicting_direction
direction_unclear
no_ich_literature_detected
```

## 10. GENCODE payload annotation

The script downloads the GENCODE v50 human GTF file.

Coding-sequence features are summed by transcript. The selected transcript is:

1. MANE Select, when present;
2. otherwise, the longest basic protein-coding transcript.

The CDS length is reported in base pairs.

## 11. Payload score

| Selected CDS length | Score |
|---:|---:|
| ≤800 bp | 5 |
| 801–1,500 bp | 4 |
| 1,501–2,250 bp | 3 |
| 2,251–3,000 bp | 2 |
| 3,001–4,000 bp | 1 |
| >4,000 bp or unavailable | 0 |

## 12. Final priority score

The default final score is:

```text
priority score =
20 × [
    age score × 0.20
    + ICH score × 0.60
    + payload score × 0.20
]
```

Age scores are:

| Age tier | Score |
|---|---:|
| High-confidence | 5 |
| Moderate | 4 |
| Directional | 3 |
| Not prioritized | 0 |

Genes assigned `restoration_contraindicated` receive a final priority score of
zero.

The final ranked table is sorted by:

1. priority score;
2. ICH evidence score;
3. ICH article count;
4. age score;
5. payload score;
6. genome-wide age rank;
7. gene symbol.

## 13. Generated figures

The code generates:

- number of genes in each age-evidence tier;
- GPX4 age effect across the 13 CNS regions;
- distribution of automated ICH evidence scores;
- top final priority scores.
