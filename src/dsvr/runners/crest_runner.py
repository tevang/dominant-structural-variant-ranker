from __future__ import annotations

import csv
import gzip
import json
import re
import shlex
import shutil
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.filtering.disk_guard import DiskLimitError, directory_size_gb, enforce_run_disk_limit
from dsvr.models import CrestConformerRecord, SeedConformerRecord, make_crest_conformer_id
from dsvr.parsing.crest_outputs import ParsedCrestOutput, parse_crest_outputs
from dsvr.runners.subprocess_utils import ExternalToolError, run_command
from dsvr.workflow.provenance import write_jsonl


class CrestUnavailableError(RuntimeError):
    """Raised when CREST is required but unavailable."""


class CrestExecutionError(RuntimeError):
    """Raised when CREST execution fails unexpectedly."""


def run_crest_for_seed(
    seed_record: SeedConformerRecord,
    config: RunConfig,
) -> list[CrestConformerRecord]:
    workdir = crest_workdir(seed_record, config)
    workdir.mkdir(parents=True, exist_ok=True)
    _enforce_disk_guard(config.output_dir, config, seed_record)
    input_xyz = workdir / "input.xyz"
    _write_seed_xyz(seed_record, input_xyz)

    warnings: list[str] = []
    charge = _formal_charge(seed_record)
    uhf = _uhf(seed_record, charge, warnings)
    command = build_crest_command(
        input_xyz=input_xyz,
        workdir=workdir,
        charge=charge,
        uhf=uhf,
        config=config,
        warnings=warnings,
    )

    if config.resume and _has_crest_outputs(workdir):
        parsed = parse_crest_outputs(workdir)
        records = _records_from_parsed(seed_record, config, workdir, command, parsed, warnings)
        _write_outputs(workdir, records)
        return records

    try:
        run_command(
            command,
            cwd=workdir,
            timeout_s=None
            if config.crest.walltime_minutes is None
            else config.crest.walltime_minutes * 60,
            log_dir=workdir / "logs",
            command_name="crest",
            check=True,
            show_progress=config.logging.tail_subprocess_logs,
        )
    except ExternalToolError as exc:
        record = _failure_record(
            seed_record,
            config,
            workdir,
            command,
            warnings + [f"CREST failed: {exc}"],
            exc.metadata,
        )
        _write_outputs(workdir, [record])
        return [record]
    except DiskLimitError as exc:
        record = _failure_record(
            seed_record,
            config,
            workdir,
            command,
            warnings + [f"CREST disk guard stopped job: {exc}"],
            {"disk_guard": str(exc)},
        )
        _write_outputs(workdir, [record])
        return [record]

    parsed = parse_crest_outputs(workdir)
    if not parsed.conformers:
        record = _failure_record(
            seed_record,
            config,
            workdir,
            command,
            warnings + parsed.warnings + ["CREST completed but no conformers were parsed."],
            {},
        )
        _write_outputs(workdir, [record])
        return [record]
    records = _records_from_parsed(seed_record, config, workdir, command, parsed, warnings)
    _write_outputs(workdir, records)
    _cleanup_crest_workdir(workdir, config)
    _enforce_disk_guard(config.output_dir, config, seed_record)
    return records


def build_crest_command(
    *,
    input_xyz: Path,
    workdir: Path,
    charge: int,
    uhf: int,
    config: RunConfig,
    warnings: list[str] | None = None,
) -> list[str]:
    if config.crest.command_template:
        return _command_from_template(
            config.crest.command_template,
            input_xyz=input_xyz,
            charge=charge,
            uhf=uhf,
            config=config,
        )

    executable = shutil.which(config.crest.executable)
    if executable is None:
        raise CrestUnavailableError(
            f"CREST executable '{config.crest.executable}' was not found on PATH. "
            "Install CREST or set crest.executable/crest.command_template in the config."
        )

    help_text = _crest_help(executable, workdir)
    if not help_text:
        raise CrestUnavailableError(
            "CREST help/version could not be inspected. Set crest.command_template for "
            "this CREST installation to avoid unvalidated default flags."
        )

    command = [executable, str(input_xyz.resolve())]
    _append_if_supported(command, help_text, f"--gfn{config.crest.gfn}", warnings)
    _append_pair_if_supported(command, help_text, "--chrg", str(charge), warnings)
    if uhf > 0:
        _append_pair_if_supported(command, help_text, "--uhf", str(uhf), warnings)
    if config.chemistry.solvent_model != "none":
        solvent_flag = f"--{config.chemistry.solvent_model}"
        _append_pair_if_supported(
            command,
            help_text,
            solvent_flag,
            config.chemistry.solvent,
            warnings,
        )
    _append_pair_if_supported(
        command,
        help_text,
        "--ewin",
        str(config.crest.ewin_kcal_mol),
        warnings,
    )
    _append_pair_if_supported(command, help_text, "-T", str(config.crest.nproc), warnings)
    command.extend(config.crest.extra_args)
    return command


def crest_workdir(seed_record: SeedConformerRecord, config: RunConfig) -> Path:
    stereo_id = seed_record.parent_id or "unknown_stereo"
    return (
        config.output_dir
        / "crest"
        / seed_record.input_molecule_id
        / stereo_id
        / seed_record.id
    )


def read_seed_sdf(path: Path) -> list[SeedConformerRecord]:
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
    records: list[SeedConformerRecord] = []
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        seed_id = _prop_or_default(molecule, "DSVR_SEED_ID", molecule.GetProp("_Name"))
        parent_id = _prop_or_default(molecule, "DSVR_PARENT_STEREO_ID", "unknown_stereo")
        input_id = _prop_or_default(molecule, "DSVR_INPUT_ID", "unknown_input")
        canonical_smiles = _prop_or_default(
            molecule,
            "DSVR_CANONICAL_SMILES",
            Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False),
        )
        isomeric_smiles = _prop_or_default(
            molecule,
            "DSVR_ISOMERIC_SMILES",
            Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
        )
        energy = _optional_float_prop(molecule, "DSVR_ENERGY_KCAL_MOL")
        records.append(
            SeedConformerRecord(
                id=seed_id,
                parent_id=parent_id,
                input_molecule_id=input_id,
                molname=_prop_or_default(molecule, "DSVR_MOLNAME", molecule.GetProp("_Name")),
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=_prop_or_default(
                    molecule,
                    "DSVR_FORMULA",
                    rdMolDescriptors.CalcMolFormula(molecule),
                ),
                formal_charge=int(
                    _prop_or_default(
                        molecule,
                        "DSVR_FORMAL_CHARGE",
                        str(Chem.GetFormalCharge(molecule)),
                    )
                ),
                explicit_proton_count=_optional_int_prop(molecule, "DSVR_EXPLICIT_PROTON_COUNT"),
                source_software=_prop_or_default(molecule, "DSVR_SOURCE_SOFTWARE", "sdf"),
                source_python_function="dsvr.runners.crest_runner.read_seed_sdf",
                metadata={"source_seed_sdf": str(path)},
                conformer_index=index,
                energy_kcal_mol=energy,
                rdkit_mol=molecule,
                rdkit_conformer_id=0 if molecule.GetNumConformers() else None,
                forcefield=_prop_or_default(molecule, "DSVR_FORCEFIELD", ""),
                forcefield_status=_prop_or_default(molecule, "DSVR_FORCEFIELD_STATUS", "unknown"),
                embedding_status="success",
            )
        )
    return records


def _records_from_parsed(
    seed_record: SeedConformerRecord,
    config: RunConfig,
    workdir: Path,
    command: list[str],
    parsed: ParsedCrestOutput,
    warnings: list[str],
) -> list[CrestConformerRecord]:
    records: list[CrestConformerRecord] = []
    kept_conformers = parsed.conformers[: config.crest.max_conformers_to_parse]
    for conformer in kept_conformers:
        metadata = {
            "crest": {
                "workdir": str(workdir),
                "gfn": config.crest.gfn,
                "solvent": config.chemistry.solvent,
                "solvent_model": config.chemistry.solvent_model,
                "energy_source": str(parsed.energy_source) if parsed.energy_source else None,
                "comment": conformer.comment,
            },
            "lineage": _lineage_metadata(seed_record),
        }
        records.append(
            CrestConformerRecord(
                id=make_crest_conformer_id(
                    seed_record.parent_id or seed_record.id,
                    conformer.index,
                    seed_record.canonical_smiles,
                    seed_record.isomeric_smiles,
                    metadata,
                ),
                parent_id=seed_record.id,
                input_molecule_id=seed_record.input_molecule_id,
                molname=seed_record.molname,
                canonical_smiles=seed_record.canonical_smiles,
                isomeric_smiles=seed_record.isomeric_smiles,
                molecular_formula=seed_record.molecular_formula,
                formal_charge=seed_record.formal_charge,
                explicit_proton_count=seed_record.explicit_proton_count,
                source_software="crest",
                source_command=" ".join(command),
                source_python_function="dsvr.runners.crest_runner.run_crest_for_seed",
                output_paths=[
                    workdir / "crest_conformers.xyz",
                    workdir / "crest_conformers.csv",
                    workdir / "crest_provenance.jsonl",
                ],
                warnings=[*warnings, *parsed.warnings],
                metadata=metadata,
                crest_index=conformer.index,
                energy_kcal_mol=conformer.energy_kcal_mol,
                relative_energy_kcal_mol=conformer.relative_energy_kcal_mol,
            )
        )
        if conformer.index <= config.crest.max_conformers_to_keep or config.crest.keep_raw_xyz:
            (workdir / f"crest_conformer_{conformer.index:04d}.xyz").write_text(
                conformer.xyz,
                encoding="utf-8",
            )
    return records


def _failure_record(
    seed_record: SeedConformerRecord,
    config: RunConfig,
    workdir: Path,
    command: list[str],
    warnings: list[str],
    failure_metadata: dict[str, Any],
) -> CrestConformerRecord:
    metadata = {
        "crest": {
            "workdir": str(workdir),
            "gfn": config.crest.gfn,
            "solvent": config.chemistry.solvent,
            "solvent_model": config.chemistry.solvent_model,
            "status": "failed",
        },
        "lineage": _lineage_metadata(seed_record),
        "failure": failure_metadata,
    }
    return CrestConformerRecord(
        id=make_crest_conformer_id(
            seed_record.parent_id or seed_record.id,
            0,
            seed_record.canonical_smiles,
            seed_record.isomeric_smiles,
            metadata,
        ),
        parent_id=seed_record.id,
        input_molecule_id=seed_record.input_molecule_id,
        molname=seed_record.molname,
        canonical_smiles=seed_record.canonical_smiles,
        isomeric_smiles=seed_record.isomeric_smiles,
        molecular_formula=seed_record.molecular_formula,
        formal_charge=seed_record.formal_charge,
        explicit_proton_count=seed_record.explicit_proton_count,
        source_software="crest",
        source_command=" ".join(command),
        source_python_function="dsvr.runners.crest_runner.run_crest_for_seed",
        output_paths=[workdir / "crest_failures.json", workdir / "crest_provenance.jsonl"],
        warnings=warnings,
        metadata=metadata,
        crest_index=0,
    )


def _write_outputs(workdir: Path, records: list[CrestConformerRecord]) -> None:
    write_jsonl(workdir / "crest_provenance.jsonl", records)
    _write_csv(workdir / "crest_conformers.csv", records)
    failures = [record for record in records if record.crest_index == 0 or record.warnings]
    if failures:
        (workdir / "crest_failures.json").write_text(
            json.dumps([record.model_dump(mode="json") for record in failures], indent=2) + "\n",
            encoding="utf-8",
        )
    _write_energy_table(workdir / "crest_energy_table.csv", records)


def _write_csv(path: Path, records: list[CrestConformerRecord]) -> None:
    columns = [
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "canonical_smiles",
        "isomeric_smiles",
        "molecular_formula",
        "formal_charge",
        "explicit_proton_count",
        "crest_index",
        "energy_kcal_mol",
        "relative_energy_kcal_mol",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})


def _write_energy_table(path: Path, records: list[CrestConformerRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "conformer_id",
                "crest_index",
                "energy_kcal_mol",
                "relative_energy_kcal_mol",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "conformer_id": record.id,
                    "crest_index": record.crest_index,
                    "energy_kcal_mol": record.energy_kcal_mol,
                    "relative_energy_kcal_mol": record.relative_energy_kcal_mol,
                }
            )


def _cleanup_crest_workdir(workdir: Path, config: RunConfig) -> None:
    if config.crest.delete_intermediate_xyz or config.disk.delete_intermediate_xyz:
        for pattern in config.crest.cleanup_patterns:
            for path in workdir.glob(pattern):
                if path.is_file():
                    path.unlink(missing_ok=True)
    if not config.crest.keep_raw_xyz and not config.disk.keep_raw_xyz:
        keep_names = {"input.xyz", "crest_conformers.xyz", "crest_ensemble.xyz", "crest_best.xyz"}
        keep_names.update(
            f"crest_conformer_{index:04d}.xyz"
            for index in range(1, config.crest.max_conformers_to_keep + 1)
        )
        for path in workdir.glob("*.xyz"):
            if path.name not in keep_names:
                path.unlink(missing_ok=True)
    if config.crest.compress_raw_outputs or config.disk.compress_raw_outputs:
        for path in list(workdir.glob("*.log")) + list(workdir.glob("*.out")):
            _gzip_file(path)


def _gzip_file(path: Path) -> None:
    if not path.exists() or path.suffix == ".gz":
        return
    gz_path = path.with_suffix(path.suffix + ".gz")
    with path.open("rb") as source, gzip.open(gz_path, "wb") as target:
        shutil.copyfileobj(source, target)
    path.unlink(missing_ok=True)


def _enforce_disk_guard(run_dir: Path, config: RunConfig, seed_record: SeedConformerRecord) -> None:
    enforce_run_disk_limit(run_dir, config)
    molecule_dir = config.output_dir / "crest" / seed_record.input_molecule_id
    size_gb = directory_size_gb(molecule_dir)
    if size_gb > config.disk.max_molecule_dir_gb:
        message = (
            f"CREST molecule directory {molecule_dir} is {size_gb:.2f} GiB, "
            f"above {config.disk.max_molecule_dir_gb:.2f} GiB"
        )
        if config.disk.fail_on_disk_limit:
            raise DiskLimitError(message)


def _write_seed_xyz(seed_record: SeedConformerRecord, path: Path) -> None:
    if seed_record.rdkit_mol is None:
        raise CrestExecutionError("Seed record does not contain an RDKit molecule for XYZ export.")
    molecule = Chem.Mol(seed_record.rdkit_mol)
    if molecule.GetNumConformers() == 0:
        raise CrestExecutionError("Seed record RDKit molecule has no conformer coordinates.")
    conformer = molecule.GetConformer()
    lines = [str(molecule.GetNumAtoms()), seed_record.id]
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol()} {position.x:.10f} {position.y:.10f} {position.z:.10f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _formal_charge(seed_record: SeedConformerRecord) -> int:
    if seed_record.rdkit_mol is not None:
        return Chem.GetFormalCharge(seed_record.rdkit_mol)
    return seed_record.formal_charge or 0


def _uhf(seed_record: SeedConformerRecord, charge: int, warnings: list[str]) -> int:
    if seed_record.rdkit_mol is None:
        return 0
    electron_count = sum(atom.GetAtomicNum() for atom in seed_record.rdkit_mol.GetAtoms()) - charge
    if electron_count % 2:
        warnings.append(
            "Odd electron count suspected; defaulting CREST --uhf to 1. "
            "Override with crest.command_template if a different spin treatment is required."
        )
        return 1
    return 0


def _crest_help(executable: str, workdir: Path) -> str:
    for args in ([executable, "--help"], [executable, "-h"]):
        try:
            completed = run_command(
                args,
                cwd=workdir,
                timeout_s=20,
                log_dir=workdir / "logs",
                command_name="crest_help",
                check=False,
            )
        except ExternalToolError:
            continue
        text = f"{completed.stdout}\n{completed.stderr}"
        if text.strip():
            return text
    return ""


def _append_if_supported(
    command: list[str],
    help_text: str,
    flag: str,
    warnings: list[str] | None,
) -> None:
    if flag in help_text:
        command.append(flag)
    elif warnings is not None:
        warnings.append(f"CREST help did not advertise {flag}; flag omitted.")


def _append_pair_if_supported(
    command: list[str],
    help_text: str,
    flag: str,
    value: str,
    warnings: list[str] | None,
) -> None:
    if flag in help_text:
        command.extend([flag, value])
    elif warnings is not None:
        warnings.append(f"CREST help did not advertise {flag}; flag omitted.")


def _command_from_template(
    template: str,
    *,
    input_xyz: Path,
    charge: int,
    uhf: int,
    config: RunConfig,
) -> list[str]:
    rendered = template.format(
        input_xyz=input_xyz.resolve(),
        crest_executable=config.crest.executable,
        xtb_executable=config.crest.xtb_executable,
        gfn=config.crest.gfn,
        charge=charge,
        uhf=uhf,
        solvent=config.chemistry.solvent,
        solvent_model=config.chemistry.solvent_model,
        ewin=config.crest.ewin_kcal_mol,
        nproc=config.crest.nproc,
        extra_args=" ".join(config.crest.extra_args),
    )
    return shlex.split(rendered)


def _has_crest_outputs(workdir: Path) -> bool:
    return any(
        path.exists()
        for path in (
            workdir / "crest_conformers.xyz",
            workdir / "crest_ensemble.xyz",
            workdir / "crest_best.xyz",
        )
    )


def _lineage_metadata(seed_record: SeedConformerRecord) -> dict[str, str | None]:
    stereo_id = seed_record.parent_id
    return {
        "input_id": seed_record.input_molecule_id,
        "protomer_id": _ancestor_id(stereo_id, "p"),
        "tautomer_id": _ancestor_id(stereo_id, "t"),
        "stereo_id": stereo_id,
        "seed_id": seed_record.id,
    }


def _ancestor_id(stereo_id: str | None, marker: str) -> str | None:
    if stereo_id is None:
        return None
    pattern = rf"^(.+_{marker}\d{{2}}_[0-9a-f]{{10}})"
    match = re.match(pattern, stereo_id)
    return match.group(1) if match else None


def _prop_or_default(molecule: Chem.Mol, key: str, default: str) -> str:
    value = molecule.GetProp(key).strip() if molecule.HasProp(key) else ""
    return value or default


def _optional_float_prop(molecule: Chem.Mol, key: str) -> float | None:
    if not molecule.HasProp(key):
        return None
    try:
        return float(molecule.GetProp(key))
    except ValueError:
        return None


def _optional_int_prop(molecule: Chem.Mol, key: str) -> int | None:
    if not molecule.HasProp(key):
        return None
    try:
        return int(molecule.GetProp(key))
    except ValueError:
        return None
