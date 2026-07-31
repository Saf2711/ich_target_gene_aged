import pandas as pd

from aged_ich.config import ScoreWeights
from aged_ich.scoring import assign_priority_scores


def test_priority_score_formula():
    data = pd.DataFrame(
        {
            "Gene": ["GPX4"],
            "discovery_tier": ["Directional age decline"],
            "auto_ich_score": [5],
            "payload_score": [5],
            "auto_direction_gate": ["restoration_supported"],
            "ich_article_count": [2],
            "objective_rank": [100],
        }
    )

    result = assign_priority_scores(
        data,
        weights=ScoreWeights(
            age=0.20,
            ich=0.60,
            payload=0.20,
        ),
    )

    assert result.loc[0, "priority_score"] == 92.0


def test_contraindicated_gene_receives_zero():
    data = pd.DataFrame(
        {
            "Gene": ["TEST"],
            "discovery_tier": ["High-confidence age decline"],
            "auto_ich_score": [5],
            "payload_score": [5],
            "auto_direction_gate": ["restoration_contraindicated"],
            "ich_article_count": [10],
            "objective_rank": [1],
        }
    )

    result = assign_priority_scores(data)

    assert result.loc[0, "priority_score"] == 0.0
