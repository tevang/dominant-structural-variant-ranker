# Workflow

The default DSVR workflow is physics-heavy and maintenance-oriented. It is
designed to generate chemically plausible candidate structural variants, search
their conformational ensembles, and rank generated candidates by parsed
free-energy estimates.

## Defaults

| Setting | Default |
| --- | --- |
| pH | `7.0` |
| Solvent | `water` |
| Temperature | `298.15 K` |
| Initial seeder | RDKit ETKDG |
| Optional seeder/prefilter | Auto3D |
| Main decision engine | CREST/xTB |
| High-confidence refinement | Optional CENSO |
| Final QM rescoring | Optional Psi4 or PySCF |

## Step-by-Step Table

| Step | Software | Function | Status |
| --- | --- | --- | --- |
| 1 | molscrub | Generate practical pH/protomer/protonation candidates at the target pH or pH window. | Compulsory for pH/protomer workflow |
| 2 | RDKit | Enumerate tautomers from generated candidates. | Compulsory by default |
| 3 | RDKit | Enumerate stereoisomers explicitly under configured limits. | Compulsory by default |
| 4 | RDKit ETKDG | Generate initial 3D conformer seeds. | Compulsory default seeder |
| 5 | Auto3D | Optional seed generation or prefiltering with neural-network potentials. | Optional |
| 6 | CREST/xTB | Perform conformer search and ensemble reduction. | Compulsory for physics-heavy ranking |
| 7 | xTB thermo / CREST entropy | Extract free-energy terms and rank by relative Delta G. | Compulsory for population-oriented ranking |
| 8 | CENSO | Refine conformer ensemble ranking with higher-confidence settings. | Optional |
| 9 | Psi4 or PySCF | Final quantum-chemistry rescoring of selected structures. | Optional |
| 10 | DSVR reporting | Write ranked tables, provenance, logs, summaries, and scope labels. | Compulsory |

## Operational Sequence

1. Read SMILES or SDF inputs and assign stable input hashes.
2. Standardize molecules enough to support deterministic downstream records.
3. Use molscrub for practical pH/protomer candidate generation.
4. Enumerate RDKit tautomers within configured caps.
5. Enumerate RDKit stereoisomers within configured caps.
6. Deduplicate generated candidates by configured identifiers.
7. Generate seed 3D geometries with RDKit ETKDG by default.
8. Optionally use Auto3D for seed generation or prefiltering.
9. Run CREST/xTB to search conformer space and reduce ensembles.
10. Parse xTB thermo and CREST entropy/free-energy outputs.
11. Optionally run CENSO for higher-confidence ensemble refinement.
12. Optionally rescore selected structures with Psi4 or PySCF.
13. Compute relative free energies and scoped Boltzmann weights.
14. Emit ranked outputs, logs, provenance, and reports.

## Enumeration Controls

RDKit tautomer canonicalization is not a stability ranking. Tautomer
canonicalization can be useful for deduplication or representative naming, but
the workflow must not treat the canonical tautomer as the physically dominant
tautomer.

RDKit stereoisomer enumeration is explicit and controlled. The workflow should
record stereochemistry assumptions, maximum enumeration counts, skipped states,
and whether undefined stereocenters were expanded.

Auto3D can internally enumerate tautomers and stereoisomers. DSVR should disable
that behavior by default when RDKit already performed enumeration. Enable Auto3D
internal enumeration only when the user explicitly selects that mode.

## Ranking and Population Scope

CREST/xTB is the main physics-based ranking layer. DSVR ranks generated
candidates by parsed relative free energies when available.

Boltzmann populations are derived from relative free energies and must be
labeled with their scope:

- Comparable within same formula/proton count.
- Approximate across different protonation/protomer states unless micro-pKa or
  proton chemical-potential corrections are available.

Without those corrections, DSVR reports approximate ranking of generated
candidates under the configured preparation assumptions, not rigorous
pH-dependent solution populations.
