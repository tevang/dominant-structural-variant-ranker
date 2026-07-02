from dsvr.chemistry.standardize import normalize_smiles_text


def test_normalize_smiles_text_strips_whitespace() -> None:
    assert normalize_smiles_text(" CCO ") == "CCO"

