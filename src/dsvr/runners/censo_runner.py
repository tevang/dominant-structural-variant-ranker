from __future__ import annotations

import csv
import json
import shlex
import shutil
from pathlib import Path

from dsvr.config import RunConfig
from dsvr.models import RankedVariantRecord
from dsvr.parsing.censo_outputs import CensoCandidateResult, parse_censo_output
from dsvr.runners.subprocess_utils import ExternalToolError, run_command


class CensoUnavailableError(RuntimeError):
    """Raised when CENSO refinement is requested but unavailable."""


class CensoExecutionError(RuntimeError):
    """Raised when CENSO execution cannot be prepared."""


def censo_available(config: RunConfig | None = None) -> bool:
    executable = config.refinement.censo_executable if config is not None else "censo"
    return shutil.which(executable) is not None


def refine_top_ranked_with_censo(
    ranked_records: list[RankedVariantRecord],
    config: RunConfig,
) -> list[RankedVariantRecord]:
    if not config.refinement.censo_enabled:
        return []
    executable = shutil.which(config.refinement.censo_executable)
    if executable is None and config.refinement.censo_command_template is None:
        raise CensoUnavailableError(
            f"CENSO was requested but executable '{config.refinement.censo_executable}' "
            "was not found on PATH. Install CENSO, load the relevant module, or set "
            "refinement.censo_command_template for your environment."
        )

    selected = sorted(
        ranked_records,
        key=lambda record: (
            float("inf")
            if record.relative_free_energy_kcal_mol is None
            else record.relative_free_energy_kcal_mol,
            record.rank,
            record.id,
        ),
    )[: config.refinement.max_candidates_for_refinement]

    refined: list[RankedVariantRecord] = []
    for candidate in selected:
        workdir = censo_workdir(candidate, config)
        workdir.mkdir(parents=True, exist_ok=True)
        input_path = prepare_censo_input(candidate, workdir)
        command = build_censo_command(candidate, config, input_path, workdir, executable)
        warnings = list(candidate.warnings)
        try:
            run_command(
                command,
                cwd=workdir,
                timeout_s=None
                if config.crest.walltime_minutes is None
                else config.crest.walltime_minutes * 60,
                log_dir=workdir / "logs",
                command_name="censo",
                check=True,
                show_progress=config.logging.tail_subprocess_logs,
            )
        except ExternalToolError as exc:
            warnings.append(f"CENSO failed: {exc}")
        result = parse_censo_output(_find_censo_output(workdir))
        warnings.extend(result.warnings)
        refined.append(_refined_record(candidate, result.candidates, command, workdir, warnings))

    refined = _rerank_refined(refined)
    write_censo_refined_outputs(refined, config.output_dir / "censo")
    return refined


def build_censo_command(
    candidate: RankedVariantRecord,
    config: RunConfig,
    input_path: Path,
    workdir: Path,
    executable: str | None = None,
) -> list[str]:
    if config.refinement.censo_command_template:
        rendered = config.refinement.censo_command_template.format(
            censo_executable=config.refinement.censo_executable,
            input_path=input_path,
            candidate_id=candidate.id,
            workdir=workdir,
            nproc=config.crest.nproc,
            extra_args=" ".join(config.refinement.censo_extra_args),
        )
        return shlex.split(rendered)
    return [
        executable or config.refinement.censo_executable,
        str(input_path),
        *config.refinement.censo_extra_args,
    ]


def prepare_censo_input(candidate: RankedVariantRecord, workdir: Path) -> Path:
    source_workdir = candidate.metadata.get("ranking", {}).get("source_workdir")
    if source_workdir is None:
        raise CensoExecutionError(
            f"Ranked candidate {candidate.id} lacks source CREST workdir metadata."
        )
    source = _select_crest_input(Path(source_workdir), candidate)
    target = workdir / source.name
    target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    return target


def censo_workdir(candidate: RankedVariantRecord, config: RunConfig) -> Path:
    return config.output_dir / "censo" / candidate.id


def write_censo_refined_outputs(records: list[RankedVariantRecord], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_refined_csv(output_dir / "ranked_variants_refined.csv", records)
    (output_dir / "ranked_variants_refined.json").write_text(
        json.dumps([record.model_dump(mode="json") for record in records], indent=2) + "\n",
        encoding="utf-8",
    )


def _refined_record(
    candidate: RankedVariantRecord,
    censo_candidates: list[CensoCandidateResult],
    command: list[str],
    workdir: Path,
    warnings: list[str],
) -> RankedVariantRecord:
    best = censo_candidates[0] if censo_candidates else None
    score = best.free_energy_kcal_mol if best is not None else candidate.score_kcal_mol
    relative = (
        best.relative_free_energy_kcal_mol
        if best is not None and best.relative_free_energy_kcal_mol is not None
        else candidate.relative_free_energy_kcal_mol
    )
    population = (
        best.population
        if best is not None and best.population is not None
        else candidate.boltzmann_population
    )
    data = candidate.model_dump(mode="python")
    data.update(
        {
            "source_software": "censo",
            "source_command": " ".join(command),
            "source_python_function": "dsvr.runners.censo_runner.refine_top_ranked_with_censo",
            "score_kcal_mol": score,
            "relative_free_energy_kcal_mol": relative,
            "boltzmann_population": population,
            "warnings": warnings,
            "metadata": candidate.metadata
            | {
                "censo": {
                    "workdir": str(workdir),
                    "candidate_count": len(censo_candidates),
                    "refined": bool(censo_candidates),
                    "preserves_preliminary_ranking_id": candidate.id,
                }
            },
        }
    )
    return RankedVariantRecord.model_validate(data)


def _rerank_refined(records: list[RankedVariantRecord]) -> list[RankedVariantRecord]:
    sorted_records = sorted(
        records,
        key=lambda record: (
            float("inf")
            if record.relative_free_energy_kcal_mol is None
            else record.relative_free_energy_kcal_mol,
            record.id,
        ),
    )
    reranked: list[RankedVariantRecord] = []
    for rank, record in enumerate(sorted_records, start=1):
        data = record.model_dump(mode="python")
        data["rank"] = rank
        reranked.append(RankedVariantRecord.model_validate(data))
    return reranked


def _select_crest_input(source_workdir: Path, candidate: RankedVariantRecord) -> Path:
    candidates = [
        source_workdir / "crest_conformers.xyz",
        source_workdir / "crest_ensemble.xyz",
        source_workdir / "crest_best.xyz",
    ]
    candidates.extend(sorted(source_workdir.glob("crest_conformer_*.xyz")))
    for path in candidates:
        if path.exists():
            return path
    raise CensoExecutionError(
        f"No CREST conformer ensemble found for {candidate.id} in {source_workdir}."
    )


def _find_censo_output(workdir: Path) -> Path:
    for name in ("censo.out", "censo.log", "output.log"):
        path = workdir / name
        if path.exists():
            return path
    return workdir / "censo.out"


def _write_refined_csv(path: Path, records: list[RankedVariantRecord]) -> None:
    columns = [
        "rank",
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "score_kcal_mol",
        "relative_free_energy_kcal_mol",
        "boltzmann_population",
        "population_scope",
        "approximate_population",
        "source_software",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})
