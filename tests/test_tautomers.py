from dsvr.chemistry.tautomers import enumerate_tautomers_placeholder


def test_tautomer_placeholder_returns_input() -> None:
    assert enumerate_tautomers_placeholder("CCO") == ["CCO"]

