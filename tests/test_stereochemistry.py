from dsvr.chemistry.stereochemistry import enumerate_stereoisomers_placeholder


def test_stereochemistry_placeholder_returns_input() -> None:
    assert enumerate_stereoisomers_placeholder("CCO") == ["CCO"]

