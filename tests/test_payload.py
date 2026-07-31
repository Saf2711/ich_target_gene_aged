from aged_ich.payload import assign_payload_score


def test_payload_score_boundaries():
    assert assign_payload_score(800) == 5
    assert assign_payload_score(801) == 4
    assert assign_payload_score(1500) == 4
    assert assign_payload_score(1501) == 3
    assert assign_payload_score(2250) == 3
    assert assign_payload_score(2251) == 2
    assert assign_payload_score(3000) == 2
    assert assign_payload_score(3001) == 1
    assert assign_payload_score(4000) == 1
    assert assign_payload_score(4001) == 0
