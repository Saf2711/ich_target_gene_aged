from __future__ import annotations

import numpy as np
import pandas as pd

from .config import AGE_SCORE_MAP, ScoreWeights


def assign_priority_scores(
    genes: pd.DataFrame,
    weights: ScoreWeights | None = None,
) -> pd.DataFrame:
    """Calculate the final age-ICH-payload priority score."""

    weights = weights or ScoreWeights()
    weights.validate()

    result = genes.copy()

    result["age_score"] = (
        result["discovery_tier"]
        .map(AGE_SCORE_MAP)
        .fillna(0)
        .astype(int)
    )

    result["priority_score"] = 20.0 * (
        result["age_score"] * weights.age
        + result["auto_ich_score"] * weights.ich
        + result["payload_score"] * weights.payload
    )

    contraindicated = result[
        "auto_direction_gate"
    ].eq(
        "restoration_contraindicated"
    )

    result.loc[
        contraindicated,
        "priority_score",
    ] = 0.0

    result = result.sort_values(
        [
            "priority_score",
            "auto_ich_score",
            "ich_article_count",
            "age_score",
            "payload_score",
            "objective_rank",
            "Gene",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            False,
            True,
            True,
        ],
    ).reset_index(drop=True)

    result.insert(
        0,
        "screen_rank",
        np.arange(1, len(result) + 1),
    )

    return result
