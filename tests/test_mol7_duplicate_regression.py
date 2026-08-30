"""Regression test for the mol_000007 duplicate case: an input molecule whose
protomer candidates are tautomers of one another (neutral multi-azole), so
every protomer branch enumerates the same tautomer space and branch-local
selection produces exact duplicates downstream.

Runs both tautomer generation paths with Auto3D mocked out and asserts that
no exact duplicates survive at the tautomer, stereoisomer, and final-3D
stages, and that branches are refilled to ``tauto_k`` after deduplication.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest
from rdkit import Chem

from dsvr.chemistry.final3d import generate_final_3d_variants
from dsvr.chemistry.identity import exact_duplicate_key
from dsvr.chemistry.protonation import _records_from_candidates
from dsvr.chemistry.stereo_auto3d_filter import filter_stereoisomers_with_auto3d
from dsvr.chemistry.stereochemistry import enumerate_stereoisomers
from dsvr.chemistry.tautomer_auto3d_filter import filter_tautomers_with_auto3d
from dsvr.chemistry.tautomers import enumerate_tautomers
from dsvr.config import RunConfig
from dsvr.models import MoleculeInput, TautomerRecord
from dsvr.workflow import engine as engine_module
from dsvr.workflow.engine import dedupe_and_refill_tautomers, run_workflow

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

MOL7_SMILES = "O=C(NC1CCCCC1)c1cccc(-c2n[nH]c3ccc(-c4nc[nH]n4)cc23)c1"
# Protomer candidate that is a tautomer of MOL7_SMILES (proton on different
# ring nitrogens of the azole system); both enumerate the same 16-member
# tautomer space, which is exactly the mol_000007 duplicate pathology.
MOL7_ALT_TAUTOMERIC_PROTOMER = (
    "O=C(NC1CCCCC1)c1cccc(-c2[nH]nc3ccc(-c4ncn[nH]4)cc23)c1"
)


def _mol7_protomers(config: RunConfig) -> list:
    mol = Chem.MolFromSmiles(MOL7_SMILES)
    alt = Chem.MolFromSmiles(MOL7_ALT_TAUTOMERIC_PROTOMER)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    molecule = MoleculeInput(
        input_id="mol_000007",
        molname="CHEMBL3957293",
        source_format="smiles",
        original_smiles=MOL7_SMILES,
        canonical_smiles=canonical,
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        rdkit_mol=mol,
    )
    output_dir = config.output_dir / "enumeration" / "protomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    return _records_from_candidates(
        molecule,
        [mol, alt],
        config=config,
        source_software="mock-molscrub",
        source_command="mock molscrub emitted tautomeric protomers",
        output_dir=output_dir,
    )


def _config(tmp_path: Path, **overrides) -> RunConfig:
    tautomer_filtering = {
        "tauto_k": 2,
        "keep_input_tautomer": False,
        "tauto_window_kcal_mol": 100000.0,
        "use_gpu": False,
    }
    tautomer_filtering.update(overrides.pop("tautomer_filtering", {}))
    enumeration = {"max_tautomers_per_protomer": 8}
    enumeration.update(overrides.pop("enumeration", {}))
    return RunConfig(
        input_path=tmp_path / "input.smi",
        output_dir=tmp_path / "run",
        overwrite=True,
        resume=False,
        tautomer_filtering=tautomer_filtering,
        enumeration=enumeration,
        **overrides,
    )


def _hash_energy(smiles: str) -> float:
    digest = int.from_bytes(hashlib.sha256(smiles.encode()).digest()[:4], "big")
    return (digest % 100000) / 100.0


def _fake_auto3d_smi(input_path: Path, output_dir: Path, **kwargs):
    output_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = output_dir / "mock_auto3d.sdf"
    writer = Chem.SDWriter(str(output_sdf))
    for line in input_path.read_text(encoding="utf-8").splitlines():
        smiles, line_id = line.split(maxsplit=1)
        mol = Chem.MolFromSmiles(smiles)
        isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        mol.SetProp("_Name", line_id)
        mol.SetProp("E_kcal_mol", str(_hash_energy(isomeric)))
        writer.write(mol)
    writer.close()
    return output_sdf, ["auto3d", "mock", str(input_path)]


def _fake_auto3d_sdf(input_path: Path, output_dir: Path, **kwargs):
    if str(input_path).endswith(".smi"):
        return _fake_auto3d_smi(input_path, output_dir, **kwargs)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = output_dir / "mock_auto3d.sdf"
    writer = Chem.SDWriter(str(output_sdf))
    supplier = Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False)
    for mol in supplier:
        if mol is None:
            continue
        mol.SetProp("E_kcal_mol", "0.0")
        writer.write(mol)
    writer.close()
    return output_sdf, ["auto3d", "mock", str(input_path)]


def _exact_keys(records) -> list:
    return [exact_duplicate_key(record.rdkit_mol) for record in records]


def _assert_no_exact_duplicates(records) -> None:
    keys = _exact_keys(records)
    assert len(keys) == len(set(keys)), "exact duplicates survived"


def test_auto3d_filter_path_duplicate_branches_dedupe_and_refill(
    tmp_path: Path, monkeypatch
) -> None:
    config = _config(tmp_path)
    protomers = _mol7_protomers(config)
    assert len(protomers) == 2
    monkeypatch.setattr(
        "dsvr.chemistry.tautomer_auto3d_filter.run_auto3d", _fake_auto3d_smi
    )

    result = filter_tautomers_with_auto3d(protomers, config)
    # Both branches enumerate the same 16-member tautomer space and select the
    # same top-2, so four branch-local records cover only two structures.
    assert len(result.selected_records) == 4
    assert len(set(_exact_keys(result.selected_records))) == 2

    deduped = dedupe_and_refill_tautomers(
        result.selected_records, result.pool_records, config
    )

    _assert_no_exact_duplicates(deduped)
    assert len(deduped) == 4  # tauto_k per branch * 2 branches
    by_branch: dict[str, list[TautomerRecord]] = {}
    for record in deduped:
        by_branch.setdefault(record.parent_id or "", []).append(record)
    assert sorted(len(branch) for branch in by_branch.values()) == [2, 2]

    audit = (
        tmp_path / "run" / "enumeration" / "tautomers" / "tautomer_dedupe.csv"
    ).read_text(encoding="utf-8")
    merge_rows = [row for row in csv.DictReader(audit.splitlines()) if row["action"] == "merge"]
    refill_rows = [row for row in csv.DictReader(audit.splitlines()) if row["action"] == "refill"]
    shortfall_rows = [
        row for row in csv.DictReader(audit.splitlines()) if row["action"] == "shortfall"
    ]
    assert len(merge_rows) == 2
    assert len(refill_rows) == 2
    assert shortfall_rows == []
    eliminated_ids = {row["eliminated_tautomer_ids"] for row in merge_rows}
    retained_ids = {record.id for record in deduped}
    assert eliminated_ids.isdisjoint(retained_ids)


def test_rdkit_enumeration_path_duplicate_branches_dedupe(tmp_path: Path) -> None:
    config = _config(tmp_path)
    protomers = _mol7_protomers(config)

    first = enumerate_tautomers(protomers[0], config)
    second = enumerate_tautomers(protomers[1], config)
    selected = [*first.selected_records, *second.selected_records]
    pool = [*first.pool_records, *second.pool_records]
    assert len(selected) == 16  # 8 per branch, cap 8 below the 16-member space
    assert pool  # over-cap candidates are exposed for refill

    deduped = dedupe_and_refill_tautomers(selected, pool, config)

    _assert_no_exact_duplicates(deduped)
    retained_smiles = {record.isomeric_smiles for record in deduped}
    assert len(retained_smiles) == len(deduped)
    per_branch: dict[str, int] = {}
    for record in deduped:
        per_branch[record.parent_id or ""] = per_branch.get(record.parent_id or "", 0) + 1
    # The branch winning representatives keeps its 8; the other branch is
    # refilled to tauto_k=2 from its own unused unique pool.
    assert sorted(per_branch.values()) == [2, 8]

    audit_path = tmp_path / "run" / "enumeration" / "tautomers" / "tautomer_dedupe.csv"
    rows = list(csv.DictReader(audit_path.open(encoding="utf-8")))
    assert len([row for row in rows if row["action"] == "merge"]) == 8
    assert len([row for row in rows if row["action"] == "refill"]) == 2


def test_end_to_end_no_duplicate_structures_anywhere(tmp_path: Path, monkeypatch) -> None:
    input_path = tmp_path / "input.smi"
    input_path.write_text(
        f"SMILES  molname\n{MOL7_SMILES} CHEMBL3957293\n", encoding="utf-8"
    )
    config = _config(tmp_path, variant_filtering={"enabled": False}, final_3d={"use_gpu": False})

    def fake_protomers(molecule: MoleculeInput, run_config: RunConfig):
        (run_config.output_dir / "enumeration" / "protomers").mkdir(
            parents=True, exist_ok=True
        )
        return _records_from_candidates(
            molecule,
            [
                Chem.MolFromSmiles(MOL7_SMILES),
                Chem.MolFromSmiles(MOL7_ALT_TAUTOMERIC_PROTOMER),
            ],
            config=run_config,
            source_software="mock-molscrub",
            source_command="mock molscrub emitted tautomeric protomers",
            output_dir=run_config.output_dir / "enumeration" / "protomers",
        )

    monkeypatch.setattr(engine_module, "generate_protomer_candidates", fake_protomers)
    monkeypatch.setattr(
        "dsvr.chemistry.tautomer_auto3d_filter.run_auto3d", _fake_auto3d_smi
    )
    monkeypatch.setattr(
        "dsvr.chemistry.stereo_auto3d_filter.run_auto3d", _fake_auto3d_smi
    )
    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", _fake_auto3d_sdf)

    result = run_workflow(config)
    assert result.molecule_count == 1
    outdir = config.output_dir

    tautomers_sdf = list((outdir / "enumeration" / "tautomers").glob("*_tautomers.sdf"))
    assert tautomers_sdf, "expected per-branch tautomer SDFs"
    # Per-branch SDFs are branch-local; the engine-level dedupe makes the
    # downstream set unique, which is what all_tautomers.sdf reflects.
    engine_tautomers = [
        mol
        for mol in Chem.SDMolSupplier(str(outdir / "all_tautomers.sdf"), sanitize=True, removeHs=False)
        if mol is not None
    ]
    assert len(engine_tautomers) == 4
    _assert_no_exact_duplicates_tautomers_sdf(engine_tautomers)

    stereo_selected = outdir / "stereoisomers_selected.csv"
    assert stereo_selected.exists()
    with stereo_selected.open(encoding="utf-8", newline="") as handle:
        stereo_rows = list(csv.DictReader(handle))
    assert len(stereo_rows) == 4
    assert len({row["isomeric_smiles"] for row in stereo_rows}) == 4

    final_supplier = Chem.SDMolSupplier(
        str(outdir / "final_variants.sdf"), sanitize=True, removeHs=False
    )
    final_mols = [mol for mol in final_supplier if mol is not None]
    assert len(final_mols) == 4
    final_keys = [exact_duplicate_key(mol) for mol in final_mols]
    assert len(set(final_keys)) == len(final_keys)
    assert (outdir / "final_dedupe_audit.csv").exists()
    assert (outdir / "stereoisomer_filtering" / "stereo_dedupe.csv").exists()
    assert (outdir / "enumeration" / "tautomers" / "tautomer_dedupe.csv").exists()


def _assert_no_exact_duplicates_tautomers_sdf(mols) -> None:
    keys = [exact_duplicate_key(mol) for mol in mols]
    assert len(keys) == len(set(keys)), "duplicate tautomers in engine output"


def test_stereo_and_final_stages_stay_duplicate_free_for_mol7(
    tmp_path: Path, monkeypatch
) -> None:
    """Stereoisomer and final-3D stages inherit the deduped tautomer set and
    remain free of exact duplicates."""

    config = _config(tmp_path)
    protomers = _mol7_protomers(config)
    monkeypatch.setattr(
        "dsvr.chemistry.tautomer_auto3d_filter.run_auto3d", _fake_auto3d_smi
    )
    monkeypatch.setattr(
        "dsvr.chemistry.stereo_auto3d_filter.run_auto3d", _fake_auto3d_smi
    )
    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", _fake_auto3d_sdf)

    tauto_result = filter_tautomers_with_auto3d(protomers, config)
    deduped = dedupe_and_refill_tautomers(
        tauto_result.selected_records, tauto_result.pool_records, config
    )
    stereos = [
        record
        for tautomer in deduped
        for record in enumerate_stereoisomers(tautomer, config)
    ]
    stereo_result = filter_stereoisomers_with_auto3d(stereos, config)
    _assert_no_exact_duplicates(stereo_result.selected_records)

    final_result = generate_final_3d_variants(stereo_result.selected_records, config)
    _assert_no_exact_duplicates(final_result.records)
    assert len(final_result.records) == len(stereo_result.selected_records)
