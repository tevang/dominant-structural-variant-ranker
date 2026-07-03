import subprocess
from pathlib import Path

import pytest

from dsvr.config import RunConfig
from dsvr.models import RankedVariantRecord
from dsvr.parsing.censo_outputs import parse_censo_output
from dsvr.runners import censo_runner
from dsvr.runners.censo_runner import (
    CensoUnavailableError,
    censo_workdir,
    refine_top_ranked_with_censo,
)


def test_parse_mocked_censo_output(tmp_path: Path) -> None:
    output = tmp_path / "censo.out"
    output.write_text(
        "\n".join(
            [
                "CENSO final ranking",
                "CONF 1 -100.000 kcal/mol 0.000 kcal/mol pop=0.850",
                "CONF 2 -99.000 kcal/mol 1.000 kcal/mol pop=0.150",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    parsed = parse_censo_output(output)

    assert len(parsed.candidates) == 2
    assert parsed.candidates[0].free_energy_kcal_mol == -100.0
    assert parsed.candidates[0].relative_free_energy_kcal_mol == 0.0
    assert parsed.candidates[0].population == 0.85


def test_censo_disabled_is_noop(tmp_path: Path) -> None:
    config = RunConfig(input_path=tmp_path / "input.smi", output_dir=tmp_path / "run")

    assert refine_top_ranked_with_censo([_ranked(tmp_path, "a", 0.0)], config) == []


def test_censo_requested_but_unavailable_fails_early(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = RunConfig(
        input_path=tmp_path / "input.smi",
        output_dir=tmp_path / "run",
        refinement={"censo_enabled": True},
    )
    monkeypatch.setattr(censo_runner.shutil, "which", lambda executable: None)

    with pytest.raises(CensoUnavailableError, match="Install CENSO"):
        refine_top_ranked_with_censo([_ranked(tmp_path, "a", 0.0)], config)


def test_censo_refines_top_n_and_preserves_preliminary_ranking(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = RunConfig(
        input_path=tmp_path / "input.smi",
        output_dir=tmp_path / "run",
        refinement={"censo_enabled": True, "max_candidates_for_refinement": 1},
    )
    preliminary_dir = config.output_dir / "ranking"
    preliminary_dir.mkdir(parents=True)
    preliminary_csv = preliminary_dir / "ranked_variants.csv"
    preliminary_csv.write_text("preliminary\n", encoding="utf-8")
    candidates = [_ranked(tmp_path, "a", 0.0), _ranked(tmp_path, "b", 5.0)]
    monkeypatch.setattr(censo_runner.shutil, "which", lambda executable: "/mock/censo")

    def fake_run_command(command, cwd=None, **kwargs):
        assert cwd is not None
        Path(cwd, "censo.out").write_text(
            "CONF 1 -101.000 kcal/mol 0.000 kcal/mol pop=1.000\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="censo ok\n", stderr="")

    monkeypatch.setattr(censo_runner, "run_command", fake_run_command)

    refined = refine_top_ranked_with_censo(candidates, config)

    assert len(refined) == 1
    assert refined[0].source_software == "censo"
    assert refined[0].score_kcal_mol == -101.0
    assert (config.output_dir / "censo" / "ranked_variants_refined.csv").exists()
    assert preliminary_csv.read_text(encoding="utf-8") == "preliminary\n"
    assert censo_workdir(candidates[0], config).exists()
    assert not censo_workdir(candidates[1], config).exists()


def test_censo_command_template_is_used(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = RunConfig(
        input_path=tmp_path / "input.smi",
        output_dir=tmp_path / "run",
        refinement={
            "censo_enabled": True,
            "censo_command_template": "custom-censo --input {input_path} --nproc {nproc}",
        },
    )
    candidate = _ranked(tmp_path, "a", 0.0)
    commands = []

    def fake_run_command(command, cwd=None, **kwargs):
        commands.append(command)
        Path(cwd, "censo.out").write_text(
            "final Gibbs free energy = -10.0 kcal/mol\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="ok\n", stderr="")

    monkeypatch.setattr(censo_runner, "run_command", fake_run_command)

    refine_top_ranked_with_censo([candidate], config)

    assert commands[0][0] == "custom-censo"
    assert "--input" in commands[0]


def _ranked(tmp_path: Path, suffix: str, relative: float) -> RankedVariantRecord:
    source_workdir = tmp_path / f"crest_{suffix}"
    source_workdir.mkdir()
    (source_workdir / "crest_conformers.xyz").write_text(
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
