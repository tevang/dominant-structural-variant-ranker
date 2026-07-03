from dsvr.models import (
    make_crest_conformer_id,
    make_input_id,
    make_protomer_id,
    make_seed_id,
    make_stereo_id,
    make_tautomer_id,
)


def test_lineage_ids_are_deterministic_and_nested() -> None:
    input_id = make_input_id("65986444", "CCC(C)CC(C)(CCN)O")
    protomer_id = make_protomer_id(input_id, 1, "CCC(C)CC(C)(CCN)O", "CCC(C)CC(C)(CCN)O")
    tautomer_id = make_tautomer_id(
        protomer_id,
        1,
        "CCC(C)CC(C)(CCN)O",
        "CCC(C)CC(C)(CCN)O",
    )
    stereo_id = make_stereo_id(
        tautomer_id,
        1,
        "CCC(C)CC(C)(CCN)O",
        "CCC(C)CC(C)(CCN)O",
    )
    seed_id = make_seed_id(
        stereo_id,
        1,
        "CCC(C)CC(C)(CCN)O",
        "CCC(C)CC(C)(CCN)O",
    )
    crest_id = make_crest_conformer_id(
        stereo_id,
        1,
        "CCC(C)CC(C)(CCN)O",
        "CCC(C)CC(C)(CCN)O",
    )

    assert protomer_id.startswith(input_id + "_p01_")
    assert tautomer_id.startswith(protomer_id + "_t01_")
    assert stereo_id.startswith(tautomer_id + "_s01_")
    assert seed_id.startswith(stereo_id + "_c01_")
    assert crest_id.startswith(stereo_id + "_crest01_")
    assert make_input_id("65986444", "CCC(C)CC(C)(CCN)O") == input_id
