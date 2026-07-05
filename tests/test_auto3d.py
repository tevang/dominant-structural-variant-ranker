import csv
from pathlib import Path

from rdkit import Chem
from typer.testing import CliRunner

from dsvr.chemistry import conformers_auto3d
from dsvr.chemistry.conformers_auto3d import (
    generate_auto3d_seeds,
    generate_auto3d_seeds_from_protomers,
    reduce_auto3d_entropy_ensemble,
    score_auto3d_representative_variants,
)
from dsvr.cli import app
from dsvr.config import RunConfig
from dsvr.models import (
    ProtomerRecord,
    SeedConformerRecord,
    StereoRecord,
    make_input_id,
    make_protomer_id,
    make_stereo_id,
    make_tautomer_id,
)
from dsvr.runners.auto3d_runner import Auto3DExecutionError, Auto3DUnavailableError


def test_generate_auto3d_seeds_from_mock_output_preserves_parent(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={"method": "auto3d", "auto3d_k": 2, "auto3d_model": "AIMNet2"},
    )

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        assert input_path.exists()
        assert k == 2
        assert model == "AIMNet2"
        assert internal_tautomer_stereo_enum is False
        output_sdf = output_dir / "mock_auto3d.sdf"
        _write_auto3d_output(output_sdf, "CCO", stereo.id, energy="-12.34", prop="E_kcal_mol")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds([stereo], config)

    assert len(records) == 1
    assert records[0].parent_id == stereo.id
    assert records[0].metadata["auto3d"]["lineage_mode"] == "post_stereo_seed"
    assert records[0].energy_kcal_mol == -12.34
    assert records[0].forcefield_status == "auto3d_optimized"
    assert any("disabled to avoid double enumeration" in w for w in records[0].warnings)
    seed_dir = tmp_path / "run" / "seeding" / "auto3d"
    assert (seed_dir / "auto3d_input.sdf").exists()
    assert (seed_dir / "auto3d_seeds.sdf").exists()
    assert (seed_dir / "auto3d_seeds.csv").exists()


def test_auto3d_internal_enum_marks_less_controlled_lineage(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={
            "method": "auto3d",
            "auto3d_internal_tautomer_stereo_enum": True,
        },
    )

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        assert internal_tautomer_stereo_enum is True
        output_sdf = output_dir / "mock_auto3d_internal.sdf"
        _write_auto3d_output(output_sdf, "CCO", None, energy="-1.0", prop="E_kcal_mol")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds([stereo], config)

    assert len(records) == 1
    assert records[0].parent_id == stereo.id
    assert records[0].metadata["auto3d"]["lineage_mode"] == "auto3d_internal_enum"
    assert any("less controlled" in warning for warning in records[0].warnings)


def test_cli_seed_auto3d_reports_missing_auto3d(tmp_path: Path, monkeypatch) -> None:
    stereo = _stereo("ethanol", "CCO")
    stereo_sdf = tmp_path / "stereo.sdf"
    _write_stereo_sdf(stereo_sdf, stereo)

    def missing_auto3d(*args, **kwargs):
        raise Auto3DUnavailableError("Install Auto3D before running seed-auto3d")

    monkeypatch.setattr("dsvr.cli.generate_auto3d_seeds", missing_auto3d)

    result = CliRunner().invoke(
        app,
        ["seed-auto3d", str(stereo_sdf), "--out", str(tmp_path / "out"), "--k", "5"],
    )

    assert result.exit_code != 0
    assert "Install Auto3D" in result.output


def test_cli_seed_auto3d_success_with_mocked_runner(tmp_path: Path, monkeypatch) -> None:
    stereo = _stereo("ethanol", "CCO")
    stereo_sdf = tmp_path / "stereo.sdf"
    _write_stereo_sdf(stereo_sdf, stereo)

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        assert internal_tautomer_stereo_enum is False
        output_sdf = output_dir / "mock_auto3d.sdf"
        _write_auto3d_output(output_sdf, "CCO", stereo.id, energy="-12.34", prop="E_kcal_mol")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    result = CliRunner().invoke(
        app,
        ["seed-auto3d", str(stereo_sdf), "--out", str(tmp_path / "out"), "--k", "5"],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "seeding" / "auto3d" / "auto3d_report.json").exists()
    assert "internal tautomer/stereo enumeration disabled" in result.output


def test_generate_auto3d_protocol_seeds_from_protomers(tmp_path: Path, monkeypatch) -> None:
    protomer = _protomer("ethanol", "CCO")
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
        seeding={
            "method": "auto3d",
            "auto3d_k": 3,
            "auto3d_model": "ANI2xt",
            "auto3d_max_confs": 10,
            "auto3d_patience": 200,
            "auto3d_threshold": 0.3,
            "auto3d_opt_steps": 2000,
        },
    )

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        mpi_np: int | None,
        cpu_workers: int | None,
        memory_gb: int | None,
        capacity: int | None,
        max_confs: int | None,
        patience: int | None,
        threshold: float | None,
        opt_steps: int | None,
        use_gpu: bool,
        stream_output: bool,
    ) -> tuple[Path, list[str]]:
        assert input_path.exists()
        assert k == 3
        assert model == "ANI2xt"
        assert internal_tautomer_stereo_enum is True
        assert mpi_np == 4
        assert cpu_workers is None
        assert memory_gb is None
        assert capacity is None
        assert max_confs == 10
        assert patience == 200
        assert threshold == 0.3
        assert opt_steps == 2000
        assert use_gpu is True
        assert stream_output is True
        output_sdf = output_dir / "mock_auto3d_protocol.sdf"
        _write_auto3d_output(output_sdf, "CCO", protomer.id, energy="-1.0", prop="E_kcal_mol")
        return output_sdf, ["auto3d", "run", str(input_path), "--enumerate-tautomer"]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds_from_protomers([protomer], config)

    assert len(records) == 1
    assert records[0].parent_id == protomer.id
    assert records[0].energy_kcal_mol == -1.0
    assert records[0].metadata["auto3d"]["lineage_mode"] == (
        "protomer_to_auto3d_internal_tautomer_stereo_enum"
    )
    assert (tmp_path / "run" / "seeding" / "auto3d_protocol" / "auto3d_protocol_seeds.sdf").exists()


def test_generate_auto3d_protocol_adapts_large_protomer_settings(
    tmp_path: Path,
    monkeypatch,
) -> None:
    small = _protomer("ethanol", "CCO")
    large = _protomer("flexible_large", "C" * 50)
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
        seeding={
            "method": "auto3d",
            "auto3d_mpi_np": 28,
            "auto3d_cpu_workers": 28,
            "auto3d_memory_gb": 1,
            "auto3d_capacity": 8,
            "auto3d_max_confs": 8,
            "auto3d_patience": 80,
            "auto3d_opt_steps": 1000,
            "auto3d_use_gpu": True,
        },
    )
    calls = []

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        mpi_np: int | None,
        cpu_workers: int | None,
        memory_gb: int | None,
        capacity: int | None,
        max_confs: int | None,
        patience: int | None,
        threshold: float | None,
        opt_steps: int | None,
        use_gpu: bool,
        stream_output: bool,
    ) -> tuple[Path, list[str]]:
        del k, model, threshold, use_gpu, stream_output
        ids_and_smiles = [line.split(maxsplit=1) for line in input_path.read_text().splitlines()]
        calls.append(
            {
                "ids": [parts[1] for parts in ids_and_smiles],
                "internal_enum": internal_tautomer_stereo_enum,
                "mpi_np": mpi_np,
                "cpu_workers": cpu_workers,
                "memory_gb": memory_gb,
                "capacity": capacity,
                "max_confs": max_confs,
                "patience": patience,
                "opt_steps": opt_steps,
            }
        )
        output_sdf = output_dir / "mock_auto3d_protocol.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for smiles, protomer_id in ids_and_smiles:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            mol.SetProp("_Name", protomer_id)
            mol.SetProp("E_kcal_mol", "-1.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds_from_protomers([small, large], config)

    assert len(records) == 2
    small_call = next(call for call in calls if small.id in call["ids"])
    large_call = next(call for call in calls if large.id in call["ids"])
    assert small_call["internal_enum"] is True
    assert large_call["internal_enum"] is False
    assert large_call["max_confs"] == 1
    assert large_call["capacity"] == 1
    assert large_call["memory_gb"] == 2
    assert large_call["patience"] == 30
    assert large_call["opt_steps"] == 250
    assert large_call["mpi_np"] < 28
    assert large_call["cpu_workers"] < 28
    by_parent = {record.parent_id: record for record in records}
    assert by_parent[large.id].metadata["auto3d"]["internal_tautomer_stereo_enum"] is False
    assert by_parent[large.id].metadata["auto3d"]["max_confs"] == 1
    assert any(
        "excessive conformer expansion" in warning
        for warning in by_parent[large.id].warnings
    )

    plan_path = tmp_path / "run" / "seeding" / "auto3d_protocol" / "auto3d_adaptive_plan.csv"
    plan_rows = list(csv.DictReader(plan_path.open()))
    large_plan = next(row for row in plan_rows if row["protomer_id"] == large.id)
    assert large_plan["large_molecule"] == "True"
    assert large_plan["internal_tautomer_stereo_enum"] == "False"
    assert large_plan["max_confs"] == "1"


def test_generate_auto3d_protocol_seeds_from_protomers_falls_back_for_missing_outputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    protomer_a = _protomer("ethanol", "CCO")
    protomer_b = _protomer("ethylamine", "CCN")
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
    )
    config = config.model_copy(
        update={
            "seeding": config.seeding.model_copy(
                update={"auto3d_allow_rdkit_fallback": True}
            )
        }
    )

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        output_sdf = output_dir / "mock_auto3d_protocol_partial.sdf"
        _write_auto3d_output(output_sdf, "CCO", protomer_a.id, energy="-1.0", prop="E_kcal_mol")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds_from_protomers([protomer_a, protomer_b], config)

    assert len(records) == 2
    by_parent = {record.parent_id: record for record in records}
    assert by_parent[protomer_a.id].energy_kcal_mol == -1.0
    assert by_parent[protomer_b.id].source_software == "rdkit"
    assert by_parent[protomer_b.id].energy_kcal_mol is None
    assert by_parent[protomer_b.id].embedding_status == "success"
    assert "auto3d_fallback" in by_parent[protomer_b.id].metadata
    assert any(
        "falling back to a single RDKit ETKDG seed" in warning
        for warning in by_parent[protomer_b.id].warnings
    )


def test_generate_auto3d_protocol_seeds_recovers_even_when_fallback_config_disabled(
    tmp_path: Path,
    monkeypatch,
) -> None:
    protomer = _protomer("ethylamine", "CCN")
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
    )
    config = config.model_copy(
        update={
            "seeding": config.seeding.model_copy(
                update={"auto3d_allow_rdkit_fallback": False}
            )
        }
    )

    def crashing_run_auto3d(
        input_path: Path,
        output_dir: Path,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        del input_path, kwargs
        log_dir = output_dir / "logs" / "20260705T000000Z_auto3d"
        log_dir.mkdir(parents=True)
        (log_dir / "combined.log").write_text(
            "Optimization finished: Dropped(Oscillating): 1\n"
            "OSError: File error: Invalid input file output_out.sdf\n",
            encoding="utf-8",
        )
        raise Auto3DExecutionError("Auto3D failed. Tried commands: mock")

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", crashing_run_auto3d)

    records = generate_auto3d_seeds_from_protomers([protomer], config)

    assert len(records) == 1
    assert records[0].source_software == "rdkit"
    assert "dropped all candidate structures as oscillating" in (
        records[0].metadata["auto3d_fallback"]["reason"]
    )


def test_generate_auto3d_protocol_seeds_falls_back_when_batch_crashes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    protomer = _protomer("ethylamine", "CCN")
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
    )
    config = config.model_copy(
        update={
            "seeding": config.seeding.model_copy(
                update={"auto3d_allow_rdkit_fallback": True}
            )
        }
    )

    def crashing_run_auto3d(
        input_path: Path,
        output_dir: Path,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        del input_path, kwargs
        log_dir = output_dir / "logs" / "20260705T000000Z_auto3d"
        log_dir.mkdir(parents=True)
        (log_dir / "combined.log").write_text(
            "Optimization finished: Dropped(Oscillating): 1\n"
            "OSError: File error: Invalid input file output_out.sdf\n",
            encoding="utf-8",
        )
        raise Auto3DExecutionError("Auto3D failed. Tried commands: mock")

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", crashing_run_auto3d)

    records = generate_auto3d_seeds_from_protomers([protomer], config)

    assert len(records) == 1
    record = records[0]
    assert record.parent_id == protomer.id
    assert record.source_software == "rdkit"
    assert record.energy_kcal_mol is None
    assert record.metadata["auto3d_fallback"]["reason"] == (
        "Auto3D dropped all candidate structures as oscillating during optimization"
    )
    assert (
        tmp_path
        / "run"
        / "seeding"
        / "auto3d_protocol"
        / "batch_001"
        / "auto3d_batch_failure.txt"
    ).exists()
    assert (
        tmp_path / "run" / "seeding" / "auto3d_protocol" / "auto3d_protocol_seeds.sdf"
    ).exists()



def test_reduce_auto3d_entropy_ensemble_includes_configurational_entropy(
    tmp_path: Path,
) -> None:
    protomer = _protomer("ethanol", "CCO")
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
    )
    seeds = [
        _seed_from_protomer(protomer, 1, -10.0),
        _seed_from_protomer(protomer, 2, -9.5),
    ]

    records = reduce_auto3d_entropy_ensemble(seeds, config)

    assert len(records) == 1
    assert records[0].free_energy_kcal_mol is not None
    assert records[0].free_energy_kcal_mol < -10.0
    assert records[0].entropy_cal_mol_k is not None
    assert records[0].entropy_cal_mol_k > 0.0
    assert (tmp_path / "run" / "auto3d_entropy" / "auto3d_entropy_records.csv").exists()


def test_score_auto3d_representative_variants_uses_svp_score(
    tmp_path: Path,
) -> None:
    protomer = _protomer("ethanol", "CCO")
    config = RunConfig(
        protocol="auto3d_entropy",
        input_path=tmp_path / "protomers.sdf",
        output_dir=tmp_path / "run",
        thermo={"enabled": False, "population_scope": "all_approximate"},
    )
    seeds = [
        _seed_from_protomer(protomer, 1, -10.0),
        _seed_from_protomer(protomer, 2, -9.5),
    ]

    records = score_auto3d_representative_variants(seeds, config)

    assert len(records) == 1
    assert records[0].free_energy_kcal_mol is not None
    assert records[0].free_energy_kcal_mol >= 0.0
    assert records[0].entropy_cal_mol_k is None
    assert records[0].metadata["auto3d_representative"]["candidate_conformer_count"] == 2
    assert "svp_score" in records[0].metadata
    assert (
        tmp_path / "run" / "auto3d_representatives" / "auto3d_representative_scores.csv"
    ).exists()


def _stereo(molname: str, smiles: str) -> StereoRecord:
    molecule = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    tautomer_id = make_tautomer_id(protomer_id, 1, canonical, isomeric)
    stereo_id = make_stereo_id(tautomer_id, 1, canonical, isomeric)
    return StereoRecord(
        id=stereo_id,
        parent_id=tautomer_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="",
        formal_charge=Chem.GetFormalCharge(molecule),
        explicit_proton_count=None,
        source_software="test",
        source_python_function="test",
        stereo_index=1,
        rdkit_mol=molecule,
    )


def _protomer(molname: str, smiles: str) -> ProtomerRecord:
    molecule = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    return ProtomerRecord(
        id=protomer_id,
        parent_id=input_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="C2H6O",
        formal_charge=Chem.GetFormalCharge(molecule),
        explicit_proton_count=6,
        source_software="test",
        source_python_function="test",
        protomer_index=1,
        rdkit_mol=molecule,
    )


def _seed_from_protomer(
    protomer: ProtomerRecord,
    index: int,
    energy: float,
) -> SeedConformerRecord:
    return SeedConformerRecord(
        id=f"{protomer.id}_seed_{index}",
        parent_id=protomer.id,
        input_molecule_id=protomer.input_molecule_id,
        molname=protomer.molname,
        canonical_smiles=protomer.canonical_smiles,
        isomeric_smiles=protomer.isomeric_smiles,
        molecular_formula=protomer.molecular_formula,
        formal_charge=protomer.formal_charge,
        explicit_proton_count=protomer.explicit_proton_count,
        source_software="auto3d",
        energy_kcal_mol=energy,
        conformer_index=index,
    )


def _write_auto3d_output(
    path: Path,
    smiles: str,
    stereo_id: str | None,
    *,
    energy: str,
    prop: str = "E_tot",
) -> None:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", stereo_id or "auto3d_internal_output")
    if stereo_id is not None:
        mol.SetProp("DSVR_STEREO_ID", stereo_id)
    mol.SetProp(prop, energy)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def _write_stereo_sdf(path: Path, stereo: StereoRecord) -> None:
    mol = Chem.Mol(stereo.rdkit_mol)
    mol.SetProp("_Name", stereo.id)
    mol.SetProp("DSVR_STEREO_ID", stereo.id)
    mol.SetProp("DSVR_INPUT_ID", stereo.input_molecule_id)
    mol.SetProp("DSVR_PARENT_TAUTOMER_ID", stereo.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", stereo.molname)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
