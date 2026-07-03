from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.chemistry.conformers_auto3d import generate_auto3d_seeds
from dsvr.chemistry.conformers_rdkit import generate_rdkit_seeds
from dsvr.chemistry.protonation import generate_protomer_candidates
from dsvr.chemistry.stereochemistry import enumerate_stereoisomers
from dsvr.chemistry.tautomers import enumerate_tautomers
from dsvr.config import RunConfig, write_resolved_config
from dsvr.io.read_inputs import read_molecules
from dsvr.io.write_outputs import write_final_ranked_outputs, write_json, write_ranked_csv
from dsvr.models import (
    AnyLineageRecord,
    CrestConformerRecord,
    MoleculeInput,
    ProtomerRecord,
    RankedVariantRecord,
    SeedConformerRecord,
    VariantRecord,
    WorkflowResult,
    make_protomer_id,
)
from dsvr.ranking.population import compute_delta_g_and_populations, write_ranked_outputs
from dsvr.reporting.markdown import write_run_report, write_summary_markdown
from dsvr.runners.auto3d_runner import Auto3DUnavailableError
from dsvr.runners.censo_runner import refine_top_ranked_with_censo
from dsvr.runners.crest_runner import run_crest_for_seed
from dsvr.runners.molscrub_runner import MolscrubUnavailableError
from dsvr.runners.psi4_runner import rescore_top_ranked_with_psi4
from dsvr.runners.pyscf_runner import rescore_top_ranked_with_pyscf
from dsvr.runners.xtb_runner import run_xtb_thermo
from dsvr.utils.logging import configure_logging
from dsvr.workflow.provenance import build_provenance, write_all_provenance_outputs
from dsvr.workflow.steps import (
    StepState,
    file_hash,
    mark_done,
    planned_steps,
    records_hash,
    should_skip_step,
    skipped_state,
    write_dry_run_plan,
)


def run_workflow(config: RunConfig) -> WorkflowResult:
    outdir = config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    configure_logging(level=config.logging.level, log_file=outdir / "logs" / "workflow.log")
    write_resolved_config(config, outdir)

    if config.dry_run:
        plan_path = write_dry_run_plan(config)
        manifest = _manifest(config, [], [], dry_run=True, extra={"dry_run_plan": str(plan_path)})
        write_json(outdir / "manifest.json", manifest)
        write_summary_markdown(outdir / "summary.md", molecule_count=0, variant_count=0)
        return WorkflowResult(outdir=outdir, molecule_count=0, ranked_records=[])

    steps = {step.name: step for step in planned_steps(config)}
    states: list[StepState] = []
    records: list[AnyLineageRecord] = []
    warnings: list[str] = []
    input_hash = file_hash(config.input_path)

    molecules = read_molecules(
        config.input_path,
        input_format=config.input_format,
        invalid_output_path=outdir / "invalid_inputs.csv",
    )
    _ensure_invalid_inputs_csv(outdir / "invalid_inputs.csv")
    _write_input_table(outdir / "input" / "inputs.csv", molecules)
    states.append(mark_done(steps["input"], input_hash, config, details={"count": len(molecules)}))

    if should_skip_step(steps["standardize"], records_hash(molecules), config):
        states.append(skipped_state(steps["standardize"], records_hash(molecules), config))
    else:
        molecules = _standardize_molecules(molecules, config)
        _write_input_table(outdir / "input" / "standardized_inputs.csv", molecules)
        states.append(
            mark_done(
                steps["standardize"],
                records_hash(molecules),
                config,
                details={"enabled": config.chemistry.standardize},
            )
        )

    protomers: list[ProtomerRecord] = []
    protomer_hash = records_hash(molecules)
    if should_skip_step(steps["protonation"], protomer_hash, config):
        states.append(skipped_state(steps["protonation"], protomer_hash, config))
        protomers = _load_protomers(outdir)
    else:
        for molecule in molecules:
            try:
                protomers.extend(generate_protomer_candidates(molecule, config))
            except MolscrubUnavailableError as exc:
                warnings.append(str(exc))
                protomers.extend(_fallback_protomer(molecule, config))
        states.append(
            mark_done(
                steps["protonation"],
                protomer_hash,
                config,
                details={"count": len(protomers), "fallback_warnings": warnings},
            )
        )
    records.extend(protomers)

    tautomers = _run_step_list(
        "tautomers",
        protomers,
        config,
        states,
        lambda item: enumerate_tautomers(item, config),
    )
    records.extend(tautomers)

    stereos = _run_step_list(
        "stereochemistry",
        tautomers,
        config,
        states,
        lambda item: enumerate_stereoisomers(item, config),
    )
    records.extend(stereos)

    seeds: list[SeedConformerRecord] = []
    seed_input_hash = records_hash(stereos)
    if should_skip_step(steps["seeding"], seed_input_hash, config):
        states.append(skipped_state(steps["seeding"], seed_input_hash, config))
        seeds = []
    else:
        if config.seeding.method in {"etkdg", "both"}:
            for stereo in stereos:
                seeds.extend(generate_rdkit_seeds(stereo, config))
        if config.seeding.method in {"auto3d", "both"}:
            try:
                seeds.extend(generate_auto3d_seeds(stereos, config))
            except Auto3DUnavailableError as exc:
                warnings.append(str(exc))
                if config.seeding.method == "auto3d":
                    warnings.append(
                        "Auto3D unavailable in integrated workflow; falling back to RDKit ETKDG."
                    )
                    for stereo in stereos:
                        seeds.extend(generate_rdkit_seeds(stereo, config))
        states.append(
            mark_done(
                steps["seeding"],
                seed_input_hash,
                config,
                details={"count": len(seeds)},
            )
        )
    records.extend(seeds)

    crest_records = _run_crest_or_seed_ranking(seeds, config, states, warnings)
    records.extend(crest_records)

    thermo_records = []
    if config.thermo.enabled and crest_records and _tool_available(config.crest.xtb_executable):
        for conformer in crest_records:
            thermo_records.append(run_xtb_thermo(conformer, config))
    records.extend(thermo_records)
    states.append(
        mark_done(
            steps["xtb_thermo"],
            records_hash(crest_records),
            config,
            details={"count": len(thermo_records), "enabled": config.thermo.enabled},
        )
    )

    rank_source = thermo_records if thermo_records else crest_records
    ranked = compute_delta_g_and_populations(rank_source, config)
    write_ranked_outputs(ranked, outdir / "ranking")
    states.append(
        mark_done(
            steps["ranking"],
            records_hash(rank_source),
            config,
            details={"count": len(ranked)},
        )
    )
    records.extend(ranked)

    censo_records = []
    if config.refinement.censo_enabled and _tool_available(config.refinement.censo_executable):
        censo_records = refine_top_ranked_with_censo(ranked, config)
    elif config.refinement.censo_enabled:
        warnings.append("CENSO requested but unavailable; integrated workflow skipped CENSO.")
    if censo_records:
        ranked_for_qm = censo_records
        records.extend(censo_records)
    else:
        ranked_for_qm = ranked
    states.append(
        mark_done(
            steps["censo"],
            records_hash(ranked),
            config,
            details={"count": len(censo_records), "enabled": config.refinement.censo_enabled},
        )
    )

    qm_records = _run_optional_qm(ranked_for_qm, config)
    records.extend(qm_records)
    states.append(
        mark_done(
            steps["qm"],
            records_hash(ranked_for_qm),
            config,
            details={"count": len(qm_records), "backend": config.refinement.qm_backend},
        )
    )

    write_all_provenance_outputs(records, outdir)
    write_final_ranked_outputs(outdir, ranked, config)
    write_summary_markdown(
        outdir / "summary.md",
        molecule_count=len(molecules),
        variant_count=len(ranked),
    )
    manifest = _manifest(config, states, warnings, dry_run=False)
    write_json(outdir / "manifest.json", manifest)
    write_run_report(
        outdir / "report.md",
        config=config,
        records=records,
        ranked_records=ranked,
        manifest=manifest,
        output_files=_final_output_files(outdir),
    )
    states.append(
        mark_done(
            steps["reports"],
            records_hash(records),
            config,
            details={"manifest": "manifest.json"},
        )
    )
    legacy_records = [
        VariantRecord(
            variant_id=record.id,
            parent_name=record.molname,
            smiles=record.isomeric_smiles,
            relative_energy_kcal_mol=record.relative_free_energy_kcal_mol,
            approximate_population=record.boltzmann_population,
            status="ranked",
        )
        for record in ranked
    ]
    write_ranked_csv(outdir / "ranked.csv", legacy_records)
    write_json(outdir / "provenance.json", build_provenance(config.input_path, config, molecules))
    return WorkflowResult(
        outdir=outdir,
        molecule_count=len(molecules),
        ranked_records=legacy_records,
    )


def run_smoke_workflow(config: RunConfig) -> WorkflowResult:
    return run_workflow(config)


def _run_step_list(
    step_name: str,
    inputs: list[Any],
    config: RunConfig,
    states: list[StepState],
    fn,
) -> list[Any]:
    step = {item.name: item for item in planned_steps(config)}[step_name]
    input_hash = records_hash(inputs)
    if should_skip_step(step, input_hash, config):
        states.append(skipped_state(step, input_hash, config))
        return []
    outputs = []
    for item in inputs:
        outputs.extend(fn(item))
    states.append(mark_done(step, input_hash, config, details={"count": len(outputs)}))
    return outputs


def _standardize_molecules(
    molecules: list[MoleculeInput],
    config: RunConfig,
) -> list[MoleculeInput]:
    if not config.chemistry.standardize:
        return molecules
    standardized = []
    for molecule in molecules:
        mol = Chem.Mol(molecule.rdkit_mol)
        Chem.SanitizeMol(mol)
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        standardized.append(
            molecule.model_copy(
                update={
                    "rdkit_mol": mol,
                    "canonical_smiles": canonical,
                    "isomeric_smiles": isomeric,
                }
            )
        )
    return standardized


def _fallback_protomer(molecule: MoleculeInput, config: RunConfig) -> list[ProtomerRecord]:
    output_dir = config.output_dir / "enumeration" / "protomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    mol = Chem.Mol(molecule.rdkit_mol)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    metadata = {
        "fallback": True,
        "reason": "molscrub unavailable; original molecule used as single protomer candidate",
    }
    record = ProtomerRecord(
        id=make_protomer_id(molecule.input_id, 1, canonical, isomeric, metadata),
        parent_id=molecule.input_id,
        input_molecule_id=molecule.input_id,
        molname=molecule.molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula=formula,
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=_explicit_proton_count(mol),
        source_software="dsvr-fallback",
        source_python_function="dsvr.workflow.engine._fallback_protomer",
        output_paths=[
            output_dir / f"{molecule.input_id}_protomers.sdf",
            output_dir / f"{molecule.input_id}_protomers.csv",
        ],
        warnings=[
            "molscrub unavailable; original molecule used as a single protomer candidate. "
            "This is suitable for smoke tests only."
        ],
        metadata=metadata,
        protomer_index=1,
        rdkit_mol=mol,
    )
    _write_protomer_outputs(record)
    return [record]


def _run_crest_or_seed_ranking(
    seeds: list[SeedConformerRecord],
    config: RunConfig,
    states: list[StepState],
    warnings: list[str],
) -> list[CrestConformerRecord]:
    step = {item.name: item for item in planned_steps(config)}["crest"]
    seed_hash = records_hash(seeds)
    if should_skip_step(step, seed_hash, config):
        states.append(skipped_state(step, seed_hash, config))
        return []
    records = []
    crest_tools_available = _tool_available(config.crest.executable) and _tool_available(
        config.crest.xtb_executable
    )
    if config.crest.enabled and crest_tools_available:
        for seed in seeds:
            records.extend(run_crest_for_seed(seed, config))
    else:
        if config.crest.enabled:
            warnings.append(
                "CREST/xTB requested but unavailable; using seed energies for smoke ranking."
            )
        records.extend(_crest_like_records_from_seeds(seeds, config))
    states.append(
        mark_done(
            step,
            seed_hash,
            config,
            details={"count": len(records), "enabled": config.crest.enabled},
        )
    )
    return records


def _crest_like_records_from_seeds(
    seeds: list[SeedConformerRecord],
    config: RunConfig,
) -> list[CrestConformerRecord]:
    records = []
    for index, seed in enumerate(seeds, start=1):
        metadata = {
            "ranking": {"source_workdir": None},
            "smoke_mode": True,
            "source_seed_id": seed.id,
        }
        records.append(
            CrestConformerRecord(
                id=f"{seed.id}_crest_smoke",
                parent_id=seed.id,
                input_molecule_id=seed.input_molecule_id,
                molname=seed.molname,
                canonical_smiles=seed.canonical_smiles,
                isomeric_smiles=seed.isomeric_smiles,
                molecular_formula=seed.molecular_formula,
                formal_charge=seed.formal_charge,
                explicit_proton_count=seed.explicit_proton_count,
                source_software="rdkit-seed-smoke-ranking",
                source_python_function="dsvr.workflow.engine._crest_like_records_from_seeds",
                warnings=["CREST/xTB skipped; RDKit seed energy used for smoke-mode ranking."],
                metadata=metadata,
                crest_index=index,
                energy_kcal_mol=seed.energy_kcal_mol if seed.energy_kcal_mol is not None else 0.0,
                relative_energy_kcal_mol=None,
            )
        )
    return records


def _run_optional_qm(
    ranked: list[RankedVariantRecord],
    config: RunConfig,
) -> list[RankedVariantRecord]:
    if config.refinement.qm_backend == "psi4" or config.refinement.psi4_enabled:
        return rescore_top_ranked_with_psi4(ranked, config)
    if config.refinement.qm_backend == "pyscf" or config.refinement.pyscf_enabled:
        return rescore_top_ranked_with_pyscf(ranked, config)
    return []


def _load_protomers(outdir: Path) -> list[ProtomerRecord]:
    from dsvr.chemistry.tautomers import read_protomers_sdf

    records = []
    for path in sorted((outdir / "enumeration" / "protomers").glob("*_protomers.sdf")):
        records.extend(read_protomers_sdf(path))
    return records


def _write_input_table(path: Path, molecules: list[MoleculeInput]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "input_id",
                "molname",
                "source_format",
                "canonical_smiles",
                "isomeric_smiles",
            ],
        )
        writer.writeheader()
        for molecule in molecules:
            writer.writerow(
                {
                    "input_id": molecule.input_id,
                    "molname": molecule.molname,
                    "source_format": molecule.source_format,
                    "canonical_smiles": molecule.canonical_smiles,
                    "isomeric_smiles": molecule.isomeric_smiles,
                }
            )


def _ensure_invalid_inputs_csv(path: Path) -> None:
    if path.exists():
        return
    path.write_text("record_index,raw_record,error\n", encoding="utf-8")


def _write_protomer_outputs(record: ProtomerRecord) -> None:
    sdf_path, csv_path = record.output_paths
    writer = Chem.SDWriter(str(sdf_path))
    mol = Chem.Mol(record.rdkit_mol)
    mol.SetProp("_Name", record.id)
    mol.SetProp("DSVR_PROTOMER_ID", record.id)
    mol.SetProp("DSVR_INPUT_ID", record.input_molecule_id)
    mol.SetProp("DSVR_PARENT_ID", record.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", record.molname)
    writer.write(mol)
    writer.close()
    with Path(csv_path).open("w", encoding="utf-8", newline="") as handle:
        writer_csv = csv.DictWriter(handle, fieldnames=["id", "parent_id", "warnings"])
        writer_csv.writeheader()
        writer_csv.writerow(
            {
                "id": record.id,
                "parent_id": record.parent_id,
                "warnings": " | ".join(record.warnings),
            }
        )


def _manifest(
    config: RunConfig,
    states: list[StepState],
    warnings: list[str],
    *,
    dry_run: bool,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "run_name": config.run_name,
        "input_path": str(config.input_path),
        "output_dir": str(config.output_dir),
        "dry_run": dry_run,
        "max_workers": config.max_workers,
        "tool_nproc": {"crest": config.crest.nproc},
        "steps": [
            {
                "name": state.step.name,
                "status": state.status,
                "skipped": state.skipped,
                "output_dir": str(state.output_dir),
                "details": state.details,
            }
            for state in states
        ],
        "warnings": warnings,
    } | (extra or {})


def _final_output_files(outdir: Path) -> list[Path]:
    names = [
        "manifest.json",
        "resolved_config.yaml",
        "invalid_inputs.csv",
        "inputs.csv",
        "protomers.csv",
        "tautomers.csv",
        "stereoisomers.csv",
        "seeds.csv",
        "crest_conformers.csv",
        "thermo.csv",
        "ranked_variants.csv",
        "ranked_variants.json",
        "ranked_variants.sdf",
        "report.md",
    ]
    return [outdir / name for name in names]


def _tool_available(executable: str) -> bool:
    return shutil.which(executable) is not None


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
