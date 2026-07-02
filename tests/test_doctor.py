from dsvr.utils.tool_check import check_tools


def test_doctor_returns_optional_tool_statuses() -> None:
    statuses = check_tools()

    names = {status.name for status in statuses}
    assert {"rdkit", "molscrub", "Auto3D", "xtb", "crest"}.issubset(names)
    assert all(status.required is False for status in statuses)

