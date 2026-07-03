import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from dsvr import cli
from dsvr.config import RunConfig
from dsvr.models import RankedVariantRecord
from dsvr.runners import psi4_runner, pyscf_runner
from dsvr.runners.psi4_runner import parse_psi4_energy, rescore_top_ranked_with_psi4
from dsvr.runners.pyscf_runner import parse_pyscf_energy, rescore_top_ranked_with_pyscf
from dsvr.utils.units import hartree_to_kcal_mol


def test_qm_modules_import_without_optional_backends() -> None:
    assert isinstance(psi4_runner.psi4_available(), bool)
    assert isinstance(pyscf_runner.pyscf_available(), bool)


def test_parse_mocked_psi4_output(tmp_path: Path) -> None:
    log = tmp_path / "psi4_output.log"
    log.write_text("@DF-RKS Final Energy: -40.123456\n", encoding="utf-8")

    assert parse_psi4_energy(log) == -40.123456


def test_parse_mocked_pyscf_output(tmp_path: Path) -> None:
    log = tmp_path / "pyscf_output.log"
    log.write_text("DSVR_PYSCF_ENERGY_HARTREE = -41.0\n", encoding="utf-8")

    assert parse_pyscf_energy(log) == -41.0


def test_psi4_rescoring_writes_separate_qm_ranking(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = RunConfig(
        input_path=tmp_path / "input.smi",
        output_dir=tmp_path / "run",
        refinement={"qm_backend": "psi4", "psi4_enabled": True, "max_candidates_for_refinement": 1},
    )
    candidates = [_ranked(tmp_path, "a", 0.0), _ranked(tmp_path, "b", 10.0)]
    monkeypatch.setattr(psi4_runner.shutil, "which", lambda executable: "/mock/psi4")

    def fake_run_command(command, cwd=None, **kwargs):
        assert cwd is not None
        Path(cwd, "psi4_output.log").write_text(
            "Total Energy = -40.000000\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="psi4 ok\n", stderr="")

    monkeypatch.setattr(psi4_runner, "run_command", fake_run_command)

    rescored = rescore_top_ranked_with_psi4(candidates, config)

    assert len(rescored) == 1
    assert rescored[0].source_software == "psi4"
    assert rescored[0].score_kcal_mol == pytest.approx(hartree_to_kcal_mol(-40.0))
    assert any("electronic energies only" in warning for warning in rescored[0].warnings)
    assert (config.output_dir / "qm" / "psi4" / "ranked_variants_qm.csv").exists()


def test_pyscf_rescoring_writes_separate_qm_ranking(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = RunConfig(
        input_path=tmp_path / "input.smi",
        output_dir=tmp_path / "run",
        refinement={"qm_backend": "pyscf", "pyscf_enabled": True},
    )
    candidate = _ranked(tmp_path, "a", 0.0)
    monkeypatch.setattr(pyscf_runner, "pyscf_available", lambda: True)

    def fake_run_command(command, cwd=None, **kwargs):
        assert cwd is not None
        Path(cwd, "pyscf_output.log").write_text(
            "DSVR_PYSCF_ENERGY_HARTREE = -39.5\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="pyscf ok\n", stderr="")

    monkeypatch.setattr(pyscf_runner, "run_command", fake_run_command)

    rescored = rescore_top_ranked_with_pyscf([candidate], config)

    assert len(rescored) == 1
    assert rescored[0].source_software == "pyscf"
    assert rescored[0].score_kcal_mol == pytest.approx(hartree_to_kcal_mol(-39.5))
    assert (config.output_dir / "qm" / "pyscf" / "ranked_variants_qm.csv").exists()


def test_run_qm_cli_explains_optional_nature(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    ranking_dir = run_dir / "ranking"
    ranking_dir.mkdir(parents=True)
    ranked = [_ranked(tmp_path, "a", 0.0)]
    (ranking_dir / "ranked_variants.json").write_text(
        json.dumps([record.model_dump(mode="json") for record in ranked]) + "\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(cli.app, ["run-qm", str(run_dir), "--backend", "none"])

    assert result.exit_code == 0, result.output
    assert "optional" in result.output.lower()
    assert "disabled" in result.output.lower()


def _ranked(tmp_path: Path, suffix: str, relative: float) -> RankedVariantRecord:
    source_workdir = tmp_path / f"crest_{suffix}"
    source_workdir.mkdir(exist_ok=True)
    (source_workdir / "crest_conformer_0001.xyz").write_text(
        "1\nmock\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    return RankedVariantRecord(
        id=f"ranked_{suffix}",
        parent_id=f"thermo_{suffix}",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles="[H]",
        isomeric_smiles="[H]",
        molecular_formula="H",
        formal_charge=0,
        explicit_proton_count=1,
        source_software="dsvr-ranking",
        source_python_function="test",
        metadata={"ranking": {"source_workdir": str(source_workdir)}},
        rank=1,
        score_kcal_mol=relative,
        relative_free_energy_kcal_mol=relative,
        boltzmann_population=1.0,
    )
