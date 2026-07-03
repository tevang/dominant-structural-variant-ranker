# ruff: noqa: I001
import os
import shutil
from pathlib import Path

import pytest

from dsvr.config import RunConfig
from dsvr.runners.crest_runner import run_crest_for_seed
from tests.test_crest_runner import _seed

pytestmark = pytest.mark.external


@pytest.mark.skipif(
    os.environ.get("DSVR_RUN_EXTERNAL") != "1",
    reason="External CREST/xTB tests require DSVR_RUN_EXTERNAL=1.",
)
@pytest.mark.skipif(
    shutil.which("crest") is None or shutil.which("xtb") is None,
    reason="External CREST/xTB tests require crest and xtb on PATH.",
)
def test_external_crest_xtb_smoke(tmp_path: Path) -> None:
    seed = _seed("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "mols.smi",
        output_dir=tmp_path / "run",
        crest={"enabled": True, "nproc": 1, "ewin_kcal_mol": 2.0},
        thermo={"enabled": False, "xtb_hessian": False, "xtb_thermo": False},
    )

    records = run_crest_for_seed(seed, config)

    assert records
