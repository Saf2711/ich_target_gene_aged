import pandas as pd

from aged_ich.literature import auto_ich_score


def test_two_restorative_articles_receive_score_five():
    group = pd.DataFrame(
        {
            "evidence_class": [
                "increase_beneficial",
                "increase_beneficial",
            ],
            "is_review": [False, False],
            "pmid": ["1", "2"],
        }
    )

    score, gate, _ = auto_ich_score(group)

    assert score == 5
    assert gate == "restoration_supported"


def test_beneficial_inhibition_is_contraindicated_for_restoration():
    group = pd.DataFrame(
        {
            "evidence_class": [
                "decrease_beneficial",
            ],
            "is_review": [False],
            "pmid": ["1"],
        }
    )

    score, gate, _ = auto_ich_score(group)

    assert score == 3
    assert gate == "restoration_contraindicated"
