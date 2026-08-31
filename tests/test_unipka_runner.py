"""Unit tests for the Uni-Pka runner using a fake container executable and fixtures."""

from __future__ import annotations

import stat
from pathlib import Path

import pytest
from rdkit import Chem

from dsvr.config import UnipkaConfig
from dsvr.models import MoleculeInput
from dsvr.runners import unipka_runner
from dsvr.runners.unipka_runner import (
    UnipkaExecutionError,
    UnipkaUnavailableError,
    generate_unipka_batch,
    inspect_unipka,
    parse_unipka_outputs,
    resolve_unipka_container,
    resolve_unipka_runtime,
    resolve_unipka_script,
)


@pytest.fixture(autouse=True)
def _isolate_package_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep tests hermetic: pretend no vendored containers/unipka.{sif,py} is present.

    The repo vendors containers/unipka.py (and may hold a local unipka.sif); without
    this, bare container="unipka" would resolve to those and leak bind mounts.
    """

    monkeypatch.setattr(unipka_runner, "_package_root", lambda: tmp_path / "empty_repo")


def _molecule(smiles: str, input_id: str) -> MoleculeInput:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return MoleculeInput(
        input_id=input_id,
        molname=input_id,
        source_format="smi",
        original_smiles=smiles,
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        rdkit_mol=mol,
    )


MAIN_OUTPUT_FIXTURE = """\
CC(=O)[O-]\tmol-1\t0.9750
CC(=O)O\tmol-1\t0.0250
c1c[nH]cn1\tmol-2\t0.6254
c1c[nH+]c[nH]1\tmol-2\t0.3746
CCO\tmol-3\tNA
"""

DISTRIBUTION_FIXTURE = """\
name\tinput_smi\tmicrostate_smi\tdG\toccupancy\tpH
mol-1\tCC(=O)O\tCC(=O)[O-]\t-6.7938\t0.999\t2
mol-1\tCC(=O)O\tCC(=O)O\t-5.24537\t0.001\t2
mol-1\tCC(=O)O\tCC(=O)O\t-5.24537\t0.5\t4.76
mol-1\tCC(=O)O\tCC(=O)[O-]\t-6.7938\t0.5\t4.76
mol-2\tc1c[nH]cn1\tc1c[nH]cn1\t-5.0\t0.5\t7.0
mol-2\tc1c[nH]cn1\tc1c[nH+]c[nH]1\t-6.0\t0.5\t7.0
"""


def _write_fixtures(workdir: Path) -> tuple[Path, Path]:
    output = workdir / "unipka_output.tsv"
    distribution = workdir / "unipka_distribution.tsv"
    output.write_text(MAIN_OUTPUT_FIXTURE, encoding="utf-8")
    distribution.write_text(DISTRIBUTION_FIXTURE, encoding="utf-8")
    return output, distribution


def test_parse_outputs_sorts_forms_and_maps_ids(tmp_path: Path) -> None:
    molecules = [
        _molecule("CC(=O)O", "mol-1"),
        _molecule("c1c[nH]cn1", "mol-2"),
        _molecule("CCO", "mol-3"),
    ]
    output, distribution = _write_fixtures(tmp_path)

    results = parse_unipka_outputs(molecules, output_path=output, distribution_path=distribution)

    assert len(results) == 3
    assert [f.occupancy for f in results["mol-1"].forms] == [0.9750, 0.0250]
    assert results["mol-1"].failed is False
    assert results["mol-1"].envelope.microstates == {"CC(=O)[O-]": -6.7938, "CC(=O)O": -5.24537}
    assert results["mol-2"].envelope.microstates["c1c[nH+]c[nH]1"] == -6.0
    assert results["mol-3"].failed is True
    assert "retained input state" in (results["mol-3"].failed_warning or "")


def test_parse_outputs_missing_envelope_is_empty(tmp_path: Path) -> None:
    molecules = [_molecule("CCO", "mol-x")]
    output = tmp_path / "out.tsv"
    output.write_text(MAIN_OUTPUT_FIXTURE.split("\n")[0] + "\n", encoding="utf-8")

    results = parse_unipka_outputs([molecules[0]], output_path=output, distribution_path=None)

    assert results["mol-x"].envelope.microstates == {}


def test_resolve_runtime_missing_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda _name: None)

    with pytest.raises(UnipkaUnavailableError, match="No container runtime"):
        resolve_unipka_runtime(UnipkaConfig(container="unipka"))

    probe = inspect_unipka(UnipkaConfig(container="unipka"))
    assert probe["runtime"] is None
    assert probe["runtime_error"]


def test_resolve_runtime_prefers_configured(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: f"/bin/{name}")
    config = UnipkaConfig(container="/data/unipka.sif", runtime="apptainer")

    resolved = resolve_unipka_runtime(config)

    assert resolved.name == "apptainer"
    assert resolved.executable == "/bin/apptainer"


def test_bare_container_default_maps_to_vendored_sif(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = tmp_path / "repo"
    (repo / "containers").mkdir(parents=True)
    vendored = repo / "containers" / "unipka.sif"
    vendored.write_bytes(b"sif")
    monkeypatch.setattr(unipka_runner, "_package_root", lambda: repo)

    assert resolve_unipka_container(UnipkaConfig(container="unipka")) == str(vendored)
    # absent vendored image keeps the bare name (assumed pullable)
    vendored.unlink()
    assert resolve_unipka_container(UnipkaConfig(container="unipka")) == "unipka"


@pytest.mark.parametrize(
    "name",
    ["easydock/unipka:latest", "registry.example.org:5000/unipka", "unipka:2.0"],
)
def test_namespaced_image_names_stay_image_references(name: str) -> None:
    # registry-qualified image names must not be anchored to the filesystem
    config = UnipkaConfig(container=name)
    assert resolve_unipka_container(config) == name
    assert unipka_runner.check_unipka_image(config) == name


def test_file_container_paths_still_resolve(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert resolve_unipka_container(UnipkaConfig(container="/data/unipka.sif")) == (
        "/data/unipka.sif"
    )
    (tmp_path / "local.sif").write_bytes(b"sif")
    monkeypatch.chdir(tmp_path)
    assert resolve_unipka_container(UnipkaConfig(container="./local.sif")) == str(
        tmp_path / "local.sif"
    )


def test_distribution_min_occupancy_passed_to_command(tmp_path: Path) -> None:
    host_script = tmp_path / "unipka.py"
    host_script.write_text("# x\n", encoding="utf-8")
    workdir = tmp_path / "w"
    workdir.mkdir()
    config = UnipkaConfig(
        container="img",
        runtime="docker",
        script_path=str(host_script),
        distribution_min_occupancy=0.001,
    )
    resolved = unipka_runner._ResolvedRuntime(name="docker", executable="/bin/docker")

    command = unipka_runner.build_unipka_command(
        resolved,
        config,
        input_path=workdir / "in.tsv",
        output_path=workdir / "out.tsv",
        distribution_path=workdir / "dist.tsv",
        ph=7.4,
    )

    idx = command.index("--distribution-min-occupancy")
    assert command[idx + 1] == "0.001"


def test_script_path_resolution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = tmp_path / "repo"
    (repo / "containers").mkdir(parents=True)
    vendored_script = repo / "containers" / "unipka.py"
    vendored_script.write_text("# vendored\n", encoding="utf-8")
    monkeypatch.setattr(unipka_runner, "_package_root", lambda: repo)

    # unset → vendored copy when present
    assert resolve_unipka_script(UnipkaConfig()) == vendored_script
    # empty string disables the override
    assert resolve_unipka_script(UnipkaConfig(script_path="")) is None
    # explicit relative path anchors at the repo root
    assert resolve_unipka_script(UnipkaConfig(script_path="containers/unipka.py")) == (
        vendored_script
    )
    # missing explicit path raises
    with pytest.raises(UnipkaUnavailableError, match="script override not found"):
        resolve_unipka_script(UnipkaConfig(script_path="nope/unipka.py"))
    # unset with no vendored copy → no override
    (repo / "containers" / "unipka.py").unlink()
    assert resolve_unipka_script(UnipkaConfig()) is None


def test_command_includes_script_bind_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("shutil.which", lambda name: f"/bin/{name}")
    host_script = tmp_path / "unipka.py"
    host_script.write_text("# host\n", encoding="utf-8")
    workdir = tmp_path / "w"
    workdir.mkdir()
    config = UnipkaConfig(
        container="/data/unipka.sif", runtime="apptainer", script_path=str(host_script)
    )
    resolved = resolve_unipka_runtime(config)

    command = unipka_runner.build_unipka_command(
        resolved,
        config,
        input_path=workdir / "in.tsv",
        output_path=workdir / "out.tsv",
        distribution_path=workdir / "dist.tsv",
        ph=7.4,
    )

    assert f"{host_script}:/unipka/unipka.py" in command
    docker_config = UnipkaConfig(container="img", runtime="docker", script_path=str(host_script))
    docker_cmd = unipka_runner.build_unipka_command(
        unipka_runner._ResolvedRuntime(name="docker", executable="/bin/docker"),
        docker_config,
        input_path=workdir / "in.tsv",
        output_path=workdir / "out.tsv",
        distribution_path=workdir / "dist.tsv",
        ph=7.4,
    )
    assert f"{host_script}:/unipka/unipka.py:ro" in docker_cmd


def test_stale_image_error_carries_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stale = tmp_path / "fake_docker"
    stale.write_text(
        "#!/bin/bash\necho 'unipka.py: error: unrecognized arguments: -n 4 --occupancy 0.05' "
        ">&2\nexit 2\n",
        encoding="utf-8",
    )
    stale.chmod(stale.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(
        "shutil.which", lambda name: str(stale) if name == "docker" else None
    )
    config = UnipkaConfig(container="img", runtime="docker", timeout_seconds=60)

    with pytest.raises(UnipkaExecutionError, match="outdated unipka.py"):
        generate_unipka_batch(
            [_molecule("CCO", "mol-1")], config=config, ph=7.0, workdir=tmp_path / "w"
        )


def test_batch_missing_sif_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr("shutil.which", lambda name: f"/bin/{name}")
    config = UnipkaConfig(container=str(tmp_path / "missing.sif"))

    with pytest.raises(UnipkaUnavailableError, match="image file not found"):
        generate_unipka_batch(
            [_molecule("CCO", "mol-1")], config=config, ph=7.0, workdir=tmp_path / "w"
        )


def _make_fake_container(tmp_path: Path) -> Path:
    """Docker stand-in: copies fixture outputs into the argv-declared paths."""

    script = tmp_path / "fake_docker"
    script.write_text(
        """#!/bin/bash
set -euo pipefail
# args: <exe> run --rm -v DIR:DIR <image> protonate -i IN -o OUT ... --distribution-file DIST ...
prev=""
for arg in "$@"; do
  case "$prev" in
    -i) INPUT="$arg" ;;
    -o) OUTPUT="$arg" ;;
    --distribution-file) DIST="$arg" ;;
  esac
  prev="$arg"
done
"${FAKE_UNIPKA_MODE:-ok}" "$INPUT" "$OUTPUT" "$DIST"
""",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)

    for mode_name in [tmp_path / "unipka_mode_ok"]:
        mode_name.write_text(
            """#!/bin/bash
set -euo pipefail
INPUT="$1"; OUTPUT="$2"; DIST="$3"
catalog="${FAKE_UNIPKA_CATALOG:?catalog file required}"
: > "$OUTPUT"
header_written=0
echo -e "name\\tinput_smi\\tmicrostate_smi\\tdG\\toccupancy\\tpH" > "$DIST"
while IFS=$'\\t' read -r smi name; do
  while IFS=$'\\t' read -r csmi cname cocc cdg; do
    if [ "$cname" = "$name" ]; then
      if [ "$cocc" = "NA" ]; then
        echo -e "${smi}\\t${name}\\tNA" >> "$OUTPUT"
      else
        echo -e "${csmi}\\t${cname}\\t${cocc}" >> "$OUTPUT"
        echo -e "${cname}\\t${smi}\\t${csmi}\\t${cdg}\\t${cocc}\\t7.0" >> "$DIST"
      fi
    fi
  done < "$catalog"
done < "$INPUT"
""",
            encoding="utf-8",
        )
        mode_name.chmod(mode_name.stat().st_mode | stat.S_IEXEC)
    failing = tmp_path / "unipka_mode_fail"
    failing.write_text("#!/bin/bash\nexit 1\n", encoding="utf-8")
    failing.chmod(failing.stat().st_mode | stat.S_IEXEC)
    empty = tmp_path / "unipka_mode_empty"
    empty.write_text("#!/bin/bash\n: > \"$2\"\nexit 0\n", encoding="utf-8")
    empty.chmod(empty.stat().st_mode | stat.S_IEXEC)
    return script


def test_batch_runs_fake_container_and_parses(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    docker = _make_fake_container(tmp_path)
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "\n".join(
            [
                "CC(=O)[O-]\tmol-1\t0.9750\t-6.7938",
                "CC(=O)O\tmol-1\t0.0250\t-5.24537",
                "X_HASH_NA_SENTINEL\tmol-2\tNA\t0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("FAKE_UNIPKA_MODE", str(tmp_path / "unipka_mode_ok"))
    monkeypatch.setenv("FAKE_UNIPKA_CATALOG", str(catalog))
    monkeypatch.setattr(
        "shutil.which",
        lambda name: str(docker) if name == "docker" else None,
    )
    config = UnipkaConfig(container="unipka", runtime="docker", timeout_seconds=120)
    molecules = [_molecule("CC(=O)O", "mol-1"), _molecule("C1CCCCC1", "mol-2")]
    workdir = tmp_path / "work"

    result = generate_unipka_batch(molecules, config=config, ph=7.0, workdir=workdir)

    assert result.runtime == "docker"
    assert result.container == "unipka"
    assert "protonate" in result.command
    assert [f.occupancy for f in result.results["mol-1"].forms] == [0.9750, 0.0250]
    assert result.results["mol-1"].envelope.microstates["CC(=O)[O-]"] == -6.7938
    # NA rows: the container echoes the input SMILES with NA occupancy
    assert result.results["mol-2"].failed is True
    assert result.distribution_path.exists()
    assert result.input_path.read_text(encoding="utf-8").strip() != ""


def test_batch_failing_container_raises_execution_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    docker = _make_fake_container(tmp_path)
    monkeypatch.setenv("FAKE_UNIPKA_MODE", str(tmp_path / "unipka_mode_fail"))
    monkeypatch.setenv("FAKE_UNIPKA_CATALOG", str(tmp_path / "cat.tsv"))
    monkeypatch.setattr(
        "shutil.which",
        lambda name: str(docker) if name == "docker" else None,
    )
    config = UnipkaConfig(container="unipka", runtime="docker", timeout_seconds=120)

    with pytest.raises(UnipkaExecutionError, match="exited with code"):
        generate_unipka_batch(
            [_molecule("CC(=O)O", "mol-1")], config=config, ph=7.0, workdir=tmp_path / "w"
        )


def test_batch_empty_output_raises_execution_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    docker = _make_fake_container(tmp_path)
    monkeypatch.setenv("FAKE_UNIPKA_MODE", str(tmp_path / "unipka_mode_empty"))
    monkeypatch.setenv("FAKE_UNIPKA_CATALOG", str(tmp_path / "cat.tsv"))
    monkeypatch.setattr(
        "shutil.which",
        lambda name: str(docker) if name == "docker" else None,
    )
    config = UnipkaConfig(container="unipka", runtime="docker", timeout_seconds=120)

    with pytest.raises(UnipkaExecutionError, match="no output"):
        generate_unipka_batch(
            [_molecule("CC(=O)O", "mol-1")], config=config, ph=7.0, workdir=tmp_path / "w"
        )
