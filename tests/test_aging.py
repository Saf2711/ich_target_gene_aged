import numpy as np
import pandas as pd

from aged_ich.aging import summarize_age_effects


def synthetic_long_table() -> pd.DataFrame:
    rows = []

    for region_index in range(13):
        region = f"region_{region_index + 1}"

        rows.extend(
            [
                {
                    "Gene": "DECLINE",
                    "region": region,
                    "age_coef": -0.010,
                    "p_value_age": 0.001,
                    "BH_adjusted_age": 0.010,
                },
                {
                    "Gene": "NEUTRAL",
                    "region": region,
                    "age_coef": 0.000,
                    "p_value_age": 0.900,
                    "BH_adjusted_age": 1.000,
                },
                {
                    "Gene": "INCREASE",
                    "region": region,
                    "age_coef": 0.010,
                    "p_value_age": 0.001,
                    "BH_adjusted_age": 0.010,
                },
            ]
        )

    return pd.DataFrame(rows)


def test_summary_contains_all_genes():
    summary = summarize_age_effects(synthetic_long_table())
    assert set(summary["Gene"]) == {"DECLINE", "NEUTRAL", "INCREASE"}


def test_estimated_40_year_change_is_negative_for_declining_gene():
    summary = summarize_age_effects(synthetic_long_table())
    decline = summary.loc[summary["Gene"] == "DECLINE"].iloc[0]
    assert decline["median_pct_change_40y"] < 0


def test_declining_gene_is_ranked_before_increasing_gene():
    summary = summarize_age_effects(synthetic_long_table())
    ranks = summary.set_index("Gene")["objective_rank"]
    assert ranks["DECLINE"] < ranks["INCREASE"]
