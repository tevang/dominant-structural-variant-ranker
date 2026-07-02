from __future__ import annotations

from pathlib import Path

from dsvr.chemistry.identifiers import stable_molecule_id
from dsvr.config import DsvrConfig
from dsvr.io.read_inputs import read_molecules
from dsvr.io.write_outputs import write_json, write_ranked_csv
from dsvr.models import VariantRecord, WorkflowResult
from dsvr.reporting.markdown import write_summary_markdown
from dsvr.workflow.provenance import build_provenance


def run_smoke_workflow(input_path: Path, outdir: Path, config: DsvrConfig) -> WorkflowResult:
    outdir.mkdir(parents=True, exist_ok=True)
    molecules = read_molecules(input_path)
    records = [
        VariantRecord(
            variant_id=stable_molecule_id(mol.smiles or mol.input_hash),
            parent_name=mol.name,
            smiles=mol.smiles,
            relative_energy_kcal_mol=0.0,
            approximate_population=1.0,
            status="smoke-placeholder",
        )
        for mol in molecules
    ]
    write_ranked_csv(outdir / "ranked.csv", records)
    write_json(outdir / "provenance.json", build_provenance(input_path, config, molecules))
    write_summary_markdown(
        outdir / "summary.md",
        molecule_count=len(molecules),
        variant_count=len(records),
    )
    return WorkflowResult(outdir=outdir, molecule_count=len(molecules), ranked_records=records)
