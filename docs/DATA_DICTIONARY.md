# Data dictionary

## `genomewide_brain_aging_gene_screen.csv`

| Column | Description |
|---|---|
| objective_rank | Deterministic genome-wide age rank |
| Gene | Uppercase gene symbol |
| discovery_tier | Age-evidence category |
| regions_tested | Number of CNS regions with a result |
| negative_regions | Number of regions with age coefficient <0 |
| negative_fraction | Negative regions divided by tested regions |
| fdr_sig_negative_regions | Negative regions with BH-adjusted P<0.05 |
| median_age_coef | Median regional age coefficient |
| median_pct_change_40y | Median estimated 40-year change |
| mean_age_coef | Mean regional age coefficient |
| mean_pct_change_40y | Mean estimated 40-year change |
| median_effect_percentile_among_covered | Genome-wide percentile of median effect |
| median_regional_decline_percentile | Median within-region percentile |
| best_region_fdr | Lowest regional age FDR |
| strongest_decline_region | Region with most negative age coefficient |
| strongest_region_age_coef | Most negative regional coefficient |
| strongest_region_fdr | FDR for the strongest negative region |
| median_nominal_p | Median regional nominal age P-value |

## `all_downward_genes_literature_screened.csv`

Contains all columns above plus:

| Column | Description |
|---|---|
| ich_article_count | Number of annotated ICH publications |
| primary_article_count | Number of non-review publications |
| review_article_count | Number of review publications |
| increase_beneficial_count | Increased/delivered gene plus beneficial outcome |
| decrease_beneficial_count | Decreased/inhibited gene plus beneficial outcome |
| increase_harmful_count | Increased gene plus harmful outcome |
| decrease_harmful_count | Decreased gene plus harmful outcome |
| mechanistic_primary_count | Mechanistic primary articles |
| auto_ich_score | Automated ICH evidence score |
| auto_direction_gate | Restorative direction classification |
| auto_evidence_rationale | Summary of article counts |
| top_pmids | Up to five prioritized PubMed identifiers |
| top_article_titles | Titles corresponding to the prioritized records |

## `all_age_declining_genes_ranked.csv`

Contains all age, literature, and payload fields plus:

| Column | Description |
|---|---|
| transcript_id | Selected GENCODE transcript |
| gene_type_gencode | GENCODE gene type |
| transcript_type | GENCODE transcript type |
| is_mane_select | Whether the selected transcript is MANE Select |
| is_basic | Whether the selected transcript carries the basic tag |
| cds_bp | Selected coding-sequence length |
| payload_score | Coding-sequence payload score |
| protein_coding | Protein-coding flag |
| age_score | Numeric age-evidence score |
| priority_score | Final weighted score |
| screen_rank | Final rank |
