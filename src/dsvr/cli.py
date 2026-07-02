from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from dsvr.config import DsvrConfig, load_config
from dsvr.io.read_inputs import read_molecules
from dsvr.utils.tool_check import check_tools
from dsvr.workflow.engine import run_smoke_workflow

app = typer.Typer(help="Rank pH- and solvent-dependent structural variants of small molecules.")
console = Console()


@app.callback()
def main() -> None:
    """Dominant Structural Variant Ranker CLI."""


@app.command()
def version() -> None:
    """Print package version."""
    from dsvr import __version__

    console.print(__version__)


@app.command()
def doctor() -> None:
    """Check optional external tools without importing them at package import time."""
    table = Table(title="DSVR external tool check")
    table.add_column("Tool")
    table.add_column("Kind")
    table.add_column("Required")
    table.add_column("Status")
    table.add_column("Detail")
    for item in check_tools():
        table.add_row(
            item.name,
            item.kind,
            "yes" if item.required else "optional",
            "ok" if item.available else "missing",
            item.detail,
        )
    console.print(table)


@app.command()
def inspect(
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
) -> None:
    """Inspect input molecules without running external tools."""
    molecules = read_molecules(input_path)
    table = Table(title=f"Input molecules: {input_path}")
    table.add_column("Index")
    table.add_column("Name")
    table.add_column("Format")
    table.add_column("SMILES")
    for mol in molecules:
        table.add_row(str(mol.index), mol.name, mol.input_format, mol.smiles or "")
    console.print(table)


@app.command()
def run(
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    config_path: Annotated[
        Path | None,
        typer.Option("--config", "-c", exists=True, dir_okay=False, help="YAML workflow config."),
    ] = None,
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Run output directory.")] = Path(
        "runs/dsvr"
    ),
) -> None:
    """Run the current smoke workflow scaffold."""
    config = load_config(config_path) if config_path else DsvrConfig()
    result = run_smoke_workflow(input_path=input_path, outdir=outdir, config=config)
    console.print(f"Wrote smoke workflow outputs to [bold]{result.outdir}[/bold]")
    console.print(f"Molecules: {result.molecule_count}")


if __name__ == "__main__":
    app()
