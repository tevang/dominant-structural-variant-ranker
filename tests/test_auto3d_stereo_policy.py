"""fix-auto3d-integration §4: unspecified-stereochemistry policy tests.

Verifies both ``on_unspecified_stereo`` policies (``enumerate`` up-front and
``auto3d_enumerate`` via ``--enumerate-isomer`` sub-batches), the energy
aggregation of policy-enumerated isomers, and the provenance recorded on
treated variants.
"""

from pathlib import Path

from rdkit import Chem

from dsvr.chemistry import tautomer_auto3d_filter as tautomer_filter
from dsvr.chemistry.auto3d_stereo_policy import apply_stereo_policy, has_unspecified_stereo
from dsvr.chemistry.final3d import generate_final_3d_variants
from dsvr.chemistry.tautomer_auto3d_filter import _Candidate, filter_tautomers_with_auto3d
from dsvr.config import RunConfig
from dsvr.models import (
    ProtomerRecord,
    StereoRecord,
    make_stereo_id,
    make_tautomer_id,
)

UNSPECIFIED_CHIRAL = "FC(Cl)Br"  # one undefined chiral center


def test_has_unspecified_stereo_detection():
    assert has_unspecified_stereo(Chem.MolFromSmiles(UNSPECIFIED_CHIRAL))
    assert not has_unspecified_stereo(Chem.MolFromSmiles("F[C@H](Cl)Br"))
    assert not has_unspecified_stereo(Chem.MolFromSmiles("CCO"))


def test_enumerate_policy_expands_lines_and_restores_none_needed():
    config = RunConfig()
    mol = Chem.MolFromSmiles(UNSPECIFIED_CHIRAL)
    plan = apply_stereo_policy([("x", mol)], config)
    assert plan.treatment == "enumerated_upfront"
    assert len(plan.expanded) == 2  # R/S pair
    assert all(has_unspecified_stereo(isomer) is False for _b, isomer, _l in plan.expanded)
    assert len({line_id for _b, _m, line_id in plan.expanded}) == 2

    plan_specified = apply_stereo_policy([("y", Chem.MolFromSmiles("CCO"))], config)
    assert plan_specified.treatment == "none_needed"
    assert [line for _b, _m, line in plan_specified.expanded] == ["y"]


def _protomer(smiles: str = "CC(=O)C") -> ProtomerRecord:
    mol = Chem.MolFromSmiles(smiles)
    return ProtomerRecord(
        id="mol_p01",
        parent_id="mol",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        molecular_formula="C3H6O",
        formal_charge=0,
        explicit_proton_count=6,
        source_software="test",
        protomer_index=1,
        rdkit_mol=mol,
    )


def _candidate(protomer: ProtomerRecord, index: int, smiles: str) -> _Candidate:
    mol = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    metadata = {"auto3d_tautomer_filtering": True, "tautomer_filtering_stage": "pre_stereo"}
    return _Candidate(
        index=index,
        tautomer_id=make_tautomer_id(protomer.id, index, canonical, isomeric, metadata),
        molecule=mol,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        rdkit_score=0.0,
        is_input_tautomer=index == 1,
        is_canonical_tautomer=index == 1,
    )


def _fake_energy_auto3d(lines_energies: dict[str, float], calls: list[dict]):
    def fake(input_path: Path, output_dir: Path, **kwargs):
        calls.append(kwargs)
        lines = input_path.read_text(encoding="utf-8").splitlines()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "mock.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for line in lines:
            smiles, line_id = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", line_id)
            mol.SetProp("E_kcal_mol", str(lines_energies.get(line_id, 0.0)))
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "run"]

    return fake


def test_enumerate_policy_sends_specified_isomers_and_records_provenance(tmp_path, monkeypatch) -> None:
    protomer = _protomer()
    specified = _candidate(protomer, 1, "CC(=O)C")
    chiral = _candidate(protomer, 2, UNSPECIFIED_CHIRAL)
    monkeypatch.setattr(
        tautomer_filter, "_enumerate_candidates", lambda *a, **k: ([specified, chiral], None)
    )
    seen_inputs: list[str] = []
    calls: list[dict] = []

    def fake(input_path: Path, output_dir: Path, **kwargs):
        seen_inputs.append(input_path.read_text(encoding="utf-8"))
        return _fake_energy_auto3d({}, calls)(input_path, output_dir, **kwargs)

    monkeypatch.setattr(tautomer_filter, "run_auto3d", fake)
    config = RunConfig(output_dir=tmp_path / "run", tautomer_filtering={"tauto_k": 2})

    records = filter_tautomers_with_auto3d([protomer], config)

    assert all(kw.get("isomer_enum_only") is False for kw in calls)
    all_lines = "\n".join(seen_inputs)
    # The unspecified center must never reach Auto3D: input contains the
    # enumerated R/S forms instead of the bare unspecified SMILES.
    assert "@" in all_lines
    assert f"{UNSPECIFIED_CHIRAL} " not in all_lines
    assert f"{chiral.tautomer_id}__st1" in all_lines
    treated = next(record for record in records if record.isomeric_smiles == chiral.isomeric_smiles)
    policy = treated.metadata["auto3d_tautomer_filtering"]["stereo_policy"]
    assert policy["policy"] == "enumerate"
    assert policy["treatment"] == "enumerated_upfront"
    assert policy["unspecified_count"] == 1
    assert any("unspecified stereochemistry treated" in warning for warning in treated.warnings)


def test_auto3d_enumerate_policy_enables_isomer_enumeration_for_subset(tmp_path, monkeypatch) -> None:
    protomer = _protomer()
    specified = _candidate(protomer, 1, "CC(=O)C")
    chiral = _candidate(protomer, 2, UNSPECIFIED_CHIRAL)
    monkeypatch.setattr(
        tautomer_filter, "_enumerate_candidates", lambda *a, **k: ([specified, chiral], None)
    )
    seen: dict[bool, str] = {}

    def fake(input_path: Path, output_dir: Path, **kwargs):
        seen[bool(kwargs.get("isomer_enum_only"))] = input_path.read_text(encoding="utf-8")
        return _fake_energy_auto3d({}, [])(input_path, output_dir, **kwargs)

    monkeypatch.setattr(tautomer_filter, "run_auto3d", fake)
    config = RunConfig(
        output_dir=tmp_path / "run",
        tautomer_filtering={"tauto_k": 2},
        auto3d={"on_unspecified_stereo": "auto3d_enumerate"},
    )

    records = filter_tautomers_with_auto3d([protomer], config)

    assert set(seen) == {False, True}, seen.keys()
    assert f"{chiral.tautomer_id}\n" in seen[True]  # original id, not expanded
    assert specified.tautomer_id in seen[False]
    treated = next(record for record in records if record.isomeric_smiles == chiral.isomeric_smiles)
    policy = treated.metadata["auto3d_tautomer_filtering"]["stereo_policy"]
    assert policy["policy"] == "auto3d_enumerate"
    assert policy["treatment"] == "auto3d_enumerate_isomer"


def _chiral_stereo_record() -> StereoRecord:
    mol = Chem.MolFromSmiles(UNSPECIFIED_CHIRAL)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    metadata = {"mock": True}
    tautomer_id = make_tautomer_id("mol_p01", 1, canonical, isomeric, metadata)
    return StereoRecord(
        id=make_stereo_id(tautomer_id, 1, canonical, isomeric, metadata),
        parent_id=tautomer_id,
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="CHBrClF",
        formal_charge=0,
        explicit_proton_count=1,
        source_software="test",
        stereo_index=1,
        rdkit_mol=mol,
        metadata=metadata,
    )


def test_final3d_enumerate_policy_expands_input_and_records_provenance(tmp_path, monkeypatch) -> None:
    record = _chiral_stereo_record()
    inputs: list[Path] = []

    def fake(input_path: Path, output_dir: Path, **kwargs):
        inputs.append(input_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "mock.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for mol in Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False):
            if mol is None:
                continue
            mol.SetProp("E_kcal_mol", "-1.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "run"]

    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", fake)
    config = RunConfig(output_dir=tmp_path / "run", final_3d={"use_gpu": False})

    result = generate_final_3d_variants([record], config)

    assert len(inputs) == 1
    input_text = inputs[0].read_text(encoding="utf-8")
    # Two enumerated isomers, both mapped to the base stereo id.
    assert input_text.count("$$$$\n") == 2
    assert input_text.count(record.id) >= 2
    assert result.records[0].metadata["final_3d"]["stereo_policy"]["treatment"] == "enumerated_upfront"
    assert any("unspecified stereochemistry treated" in w for w in result.records[0].warnings)
