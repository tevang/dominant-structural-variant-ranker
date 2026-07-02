from dsvr.chemistry.identifiers import stable_molecule_id


def test_stable_molecule_id_is_stable() -> None:
    assert stable_molecule_id("CCO") == stable_molecule_id("CCO")
    assert len(stable_molecule_id("CCO")) == 16

