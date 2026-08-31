"""Uni-Pka protomer generation through the EasyDock container implementation.

The Uni-Pka container (``ci-lab-cz/easydock`` ``containers/unipka``) reads
tab-separated ``SMILES<TAB>name`` records and writes one line per protonation
form: ``form_smi<TAB>name<TAB>occupancy`` (``NA`` occupancy with the input
SMILES echoed when nothing could be predicted). With ``--distribution-file``
it additionally writes the microspecies free energies and occupancies over a
pH grid. DSVR calls it once per run for the whole protonation batch.
"""

from __future__ import annotations

import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from dsvr.config import UnipkaConfig
from dsvr.models import MoleculeInput
from dsvr.runners.subprocess_utils import ExternalToolError, run_command


class UnipkaUnavailableError(RuntimeError):
    """Raised when the container runtime or Uni-Pka image is not available."""


class UnipkaExecutionError(RuntimeError):
    """Raised when the Uni-Pka container runs but fails as a whole."""


@dataclass(frozen=True)
class UnipkaForm:
    form_smiles: str
    occupancy: float


@dataclass(frozen=True)
class UnipkaEnvelope:
    """Per-microstate free energies (dG), pH-independent, for one molecule."""

    microstates: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class UnipkaMoleculeResult:
    """Container outputs for one input molecule.

    ``forms`` are sorted by decreasing occupancy; empty when the molecule
    could not be predicted (``failed`` is then True and ``failed_warning``
    explains).
    """

    input_id: str
    forms: list[UnipkaForm]
    envelope: UnipkaEnvelope
    failed: bool = False
    failed_warning: str | None = None
    source_command: str = ""


@dataclass(frozen=True)
class UnipkaBatchResult:
    results: dict[str, UnipkaMoleculeResult]
    distribution_path: Path
    input_path: Path
    command: list[str]
    runtime: str
    container: str


@dataclass(frozen=True)
class _ResolvedRuntime:
    name: str
    executable: str


#: Container-internal path of the Uni-Pka script inside the EasyDock image.
UNIPKA_IMAGE_SCRIPT = "/unipka/unipka.py"

#: Vendored current EasyDock script used when ``script_path`` is unset
#: (see https://github.com/ci-lab-cz/easydock/tree/master/containers/unipka, BSD-3-Clause).
VENDORED_SCRIPT_RELATIVE = "containers/unipka.py"

#: Every published Zenodo ``unipka.sif`` build (through record 19627026,
#: 2026-04-19) bakes the pre-occupancy ``unipka.py`` that rejects ``-n``,
#: ``--occupancy`` and ``--distribution-file``; the vendored script overrides it.
#: Set ``protonation.unipka.script_path`` to another copy, or to ``""`` to run
#: the image's own script (only correct for freshly built images).
STALE_IMAGE_HINT = (
    "The Uni-Pka image ships an outdated unipka.py without the occupancy flags "
    "(-n/--occupancy/--distribution-file). Set protonation.unipka.script_path to the "
    "current EasyDock containers/unipka/unipka.py (the repo vendors containers/unipka.py), "
    "or rebuild the image from the EasyDock recipe (docs/external_tools.md)."
)


def _package_root() -> Path:
    return Path(__file__).resolve().parents[3]


def resolve_unipka_path(value: str) -> Path:
    """Resolve a configured file path; relative values anchor at the repo root.

    Relative so that checked-in configs and ``resolved_config.yaml`` copies stay
    portable across machines and working directories.
    """

    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    vendored = _package_root() / path
    return vendored if vendored.exists() else Path.cwd() / path


def resolve_unipka_script(config: UnipkaConfig) -> Path | None:
    """Absolute path of the host ``unipka.py`` that overrides the baked image copy.

    ``script_path`` unset uses the vendored ``containers/unipka.py`` when present
    (all published Zenodo SIF builds are stale); an empty string disables the
    override; any other value must exist.
    """

    if config.script_path is None:
        vendored = _package_root() / VENDORED_SCRIPT_RELATIVE
        return vendored if vendored.exists() else None
    if not config.script_path:
        return None
    path = resolve_unipka_path(config.script_path)
    if not path.exists():
        raise UnipkaUnavailableError(
            f"Uni-Pka script override not found: {config.script_path} "
            f"(resolved to {path}). Set protonation.unipka.script_path to the current "
            "EasyDock containers/unipka/unipka.py or remove the override."
        )
    return path


def _looks_like_local_path(value: str) -> bool:
    """True for explicit file paths (.sif or absolute/./~/../ rooted).

    Everything else — including namespaced or registry-qualified image names
    like ``easydock/unipka:latest`` — is a container image reference and must
    never be anchored to the filesystem.
    """

    return value.endswith(".sif") or value.startswith(("/", "./", "../", "~"))


def resolve_unipka_container(config: UnipkaConfig) -> str:
    """Container reference with file paths resolved; the bare default may map to a SIF."""

    container = config.container
    if not container:
        raise UnipkaUnavailableError(
            "protonation.unipka.container is not set. Configure a Docker image name or an "
            "Apptainer .sif path (see docs/external_tools.md, Uni-Pka entry)."
        )
    if container == "unipka":
        vendored = _package_root() / "containers" / "unipka.sif"
        if vendored.exists():
            return str(vendored)
        return container
    if _looks_like_local_path(container):
        return str(resolve_unipka_path(container))
    return container


def resolve_unipka_runtime(config: UnipkaConfig) -> _ResolvedRuntime:
    """Pick the container runtime and executable for the configured Uni-Pka image."""

    is_sif = resolve_unipka_container(config).endswith(".sif")
    if config.runtime == "auto":
        preferred = ["apptainer", "docker"] if is_sif else ["docker", "apptainer"]
    else:
        preferred = [config.runtime]
    for runtime in preferred:
        names = ("apptainer", "singularity") if runtime == "apptainer" else (runtime,)
        for name in names:
            path = shutil.which(name)
            if path:
                return _ResolvedRuntime(name=runtime, executable=path)
    if config.runtime != "auto":
        names = ("apptainer", "singularity") if config.runtime == "apptainer" else (config.runtime,)
        raise UnipkaUnavailableError(
            f"Container runtime {config.runtime!r} was configured for Uni-Pka but no "
            f"{' or '.join(names)!r} executable is on PATH. Install Apptainer/Singularity "
            "or Docker (see docs/installation.md)."
        )
    raise UnipkaUnavailableError(
        "No container runtime found for Uni-Pka (looked for apptainer, singularity, docker). "
        "Install Apptainer or Docker, then acquire the Uni-Pka image "
        "(see docs/external_tools.md, Uni-Pka entry)."
    )


def check_unipka_image(config: UnipkaConfig) -> str | None:
    """Return an availability detail string, or None when the image is missing.

    Abspath-like or .sif values must exist on disk; bare image names are
    assumed pullable by the runtime and are checked at run time.
    """

    try:
        container = resolve_unipka_container(config)
    except UnipkaUnavailableError:
        return None
    if _looks_like_local_path(container):
        return container if Path(container).exists() else None
    return container


def inspect_unipka(config: UnipkaConfig) -> dict[str, object]:
    """Doctor-facing probe of Uni-Pka availability."""

    try:
        script = resolve_unipka_script(config)
    except UnipkaUnavailableError as exc:
        return {
            "container": config.container,
            "runtime": None,
            "runtime_error": str(exc),
            "image": None,
        }
    try:
        runtime = resolve_unipka_runtime(config)
    except UnipkaUnavailableError as exc:
        return {
            "container": config.container,
            "runtime": None,
            "runtime_error": str(exc),
            "image": None,
        }
    return {
        "container": config.container,
        "runtime": runtime.name,
        "runtime_error": None,
        "image": check_unipka_image(config),
        "script_override": str(script) if script else None,
    }


def build_unipka_command(
    resolved: _ResolvedRuntime,
    config: UnipkaConfig,
    *,
    input_path: Path,
    output_path: Path,
    distribution_path: Path,
    ph: float,
) -> list[str]:
    args = [
        "-i",
        str(input_path),
        "-o",
        str(output_path),
        "--pH",
        f"{ph:g}",
        "-n",
        str(config.max_forms),
        "--occupancy",
        f"{config.min_occupancy:g}",
        "--distribution-file",
        str(distribution_path),
        "--ph-range",
        f"{config.ph_range_low:g}",
        f"{config.ph_range_high:g}",
        "--ph-step",
        f"{config.ph_step:g}",
        "--distribution-min-occupancy",
        f"{config.distribution_min_occupancy:g}",
    ]
    workdir = input_path.parent.resolve()
    container = resolve_unipka_container(config)
    script = resolve_unipka_script(config)
    if resolved.name == "docker":
        command = [
            resolved.executable,
            "run",
            "--rm",
            "-v",
            f"{workdir}:{workdir}",
        ]
        if script is not None:
            command += ["-v", f"{script}:{UNIPKA_IMAGE_SCRIPT}:ro"]
        return [*command, container, "protonate", *args]
    command = [resolved.executable, "run", "--bind", f"{workdir}:{workdir}"]
    if script is not None:
        command += ["--bind", f"{script}:{UNIPKA_IMAGE_SCRIPT}"]
    return [*command, container, "protonate", *args]


def generate_unipka_batch(
    molecules: list[MoleculeInput],
    *,
    config: UnipkaConfig,
    ph: float,
    workdir: Path,
    log_dir: Path | None = None,
) -> UnipkaBatchResult:
    """Run Uni-Pka once for all given molecules and return per-molecule results."""

    if not molecules:
        raise ValueError("generate_unipka_batch requires at least one molecule")
    resolved = resolve_unipka_runtime(config)
    container = resolve_unipka_container(config)
    if container.endswith(".sif") and not Path(container).exists():
        raise UnipkaUnavailableError(
            f"Uni-Pka image file not found: {container}. Download the pre-built unipka.sif "
            "from Zenodo or build it from the EasyDock recipe (docs/external_tools.md)."
        )

    workdir = workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    input_path = workdir / "unipka_input.tsv"
    output_path = workdir / "unipka_output.tsv"
    distribution_path = workdir / "unipka_distribution.tsv"
    writer = input_path.open("w", encoding="utf-8")
    try:
        for molecule in molecules:
            smiles = molecule.isomeric_smiles or molecule.canonical_smiles
            writer.write(f"{smiles}\t{molecule.input_id}\n")
    finally:
        writer.close()

    command = build_unipka_command(
        resolved,
        config,
        input_path=input_path,
        output_path=output_path,
        distribution_path=distribution_path,
        ph=ph,
    )
    try:
        completed = run_command(
            command,
            cwd=input_path.parent,
            timeout_s=config.timeout_seconds,
            log_dir=log_dir,
            command_name="unipka",
            check=False,
        )
    except ExternalToolError as exc:
        raise UnipkaExecutionError(f"Uni-Pka container failed to run: {exc}") from exc
    if completed.returncode != 0:
        stderr_tail = completed.stderr.strip()[-500:]
        hint = ""
        if "unrecognized arguments" in completed.stderr and resolve_unipka_script(config) is None:
            hint = f" {STALE_IMAGE_HINT}"
        raise UnipkaExecutionError(
            f"Uni-Pka container exited with code {completed.returncode}: {stderr_tail}{hint}"
        )
    if not output_path.exists() or not output_path.read_text(encoding="utf-8").strip():
        raise UnipkaExecutionError(
            "Uni-Pka container produced no output rows; see subprocess logs and stderr."
        )

    results = parse_unipka_outputs(
        molecules,
        output_path=output_path,
        distribution_path=distribution_path,
    )
    command_text = " ".join(command)
    results = {
        input_id: UnipkaMoleculeResult(
            input_id=result.input_id,
            forms=result.forms,
            envelope=result.envelope,
            failed=result.failed,
            failed_warning=result.failed_warning,
            source_command=command_text,
        )
        for input_id, result in results.items()
    }
    return UnipkaBatchResult(
        results=results,
        distribution_path=distribution_path,
        input_path=input_path,
        command=command,
        runtime=resolved.name,
        container=container,
    )


def parse_unipka_outputs(
    molecules: list[MoleculeInput],
    *,
    output_path: Path,
    distribution_path: Path | None,
) -> dict[str, UnipkaMoleculeResult]:
    """Parse the container's occupancy output and distribution file by input id."""

    forms_by_id: dict[str, list[UnipkaForm]] = {mol.input_id: [] for mol in molecules}
    failed_ids: set[str] = set()
    for line in output_path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        form_smi, name, occupancy_raw = parts[0], parts[1], parts[2]
        if name not in forms_by_id:
            continue
        if occupancy_raw.strip() == "NA":
            failed_ids.add(name)
            continue
        try:
            occupancy = float(occupancy_raw)
        except ValueError:
            continue
        if not math.isfinite(occupancy):
            continue
        forms_by_id[name].append(UnipkaForm(form_smiles=form_smi, occupancy=occupancy))

    envelopes: dict[str, UnipkaEnvelope] = {}
    if distribution_path is not None and distribution_path.exists():
        envelopes = parse_distribution_file(
            distribution_path, {mol.input_id for mol in molecules}
        )

    results: dict[str, UnipkaMoleculeResult] = {}
    for molecule in molecules:
        forms = sorted(
            forms_by_id[molecule.input_id], key=lambda form: (-form.occupancy, form.form_smiles)
        )
        failed = molecule.input_id in failed_ids or not forms
        warning = (
            "Uni-Pka returned no protonation form for this molecule; retained input state"
            if failed
            else None
        )
        results[molecule.input_id] = UnipkaMoleculeResult(
            input_id=molecule.input_id,
            forms=forms,
            envelope=envelopes.get(molecule.input_id, UnipkaEnvelope()),
            failed=failed,
            failed_warning=warning,
        )
    return results


def parse_distribution_file(path: Path, wanted_ids: set[str]) -> dict[str, UnipkaEnvelope]:
    """Collect per-microstate dG values from the distribution TSV.

    Columns: name, input_smi, microstate_smi, dG, occupancy, pH. dG is
    pH-independent, so the first row seen per microstate is authoritative.
    """

    microstates: dict[str, dict[str, float]] = {}
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        if lineno == 0 and line.lower().startswith("name"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        name, _input_smi, microstate_smi, dg_raw = parts[0], parts[1], parts[2], parts[3]
        if name not in wanted_ids:
            continue
        bucket = microstates.setdefault(name, {})
        if microstate_smi in bucket:
            continue
        try:
            dg = float(dg_raw)
        except ValueError:
            continue
        if math.isfinite(dg):
            bucket[microstate_smi] = dg
    return {name: UnipkaEnvelope(microstates=states) for name, states in microstates.items()}
