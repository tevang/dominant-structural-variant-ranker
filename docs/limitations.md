# Limitations

DSVR is an orchestration layer for practical structural-variant ranking. It can
coordinate high-value open-source tools, but it does not make their results more
rigorous than the underlying physics, enumeration coverage, and thermodynamic
corrections support.

## pH and Protonation

Default pH handling is candidate generation, not rigorous constant-pH
thermodynamics. molscrub is used to generate practical pH/protomer/protonation
candidates at the requested pH or pH window.

DSVR does not perform rigorous pH-dependent population calculations unless a
micro-pKa/proton chemical potential correction plugin is added. Without that
plugin, free energies from CREST/xTB rank generated candidates under the chosen
solvent model, but cross-protomer and cross-charge populations are approximate.

## Tautomers

RDKit tautomer canonicalization is not stability ranking. It should not be used
to assert the dominant tautomer in solution. Stability ranking requires a
physics-based or empirically calibrated scoring layer and careful scope labels.

## Stereoisomers

RDKit stereoisomer enumeration is explicit and controlled. The workflow can
expand undefined stereocenters, but the resulting candidates are only as
complete as the configured enumeration caps and input chemistry allow.

## Auto3D Enumeration

Auto3D is useful for seed conformer generation and optional prefiltering. It can
also perform internal tautomer and stereoisomer enumeration. DSVR must not let
Auto3D double-enumerate tautomers or stereoisomers after RDKit has already done
that work unless the user explicitly enables Auto3D internal enumeration.

## CREST/xTB Ranking

CREST/xTB is the main physics-based ranking layer. Its rankings depend on
conformer search completeness, xTB method choice, solvation model, charge,
multiplicity, temperature treatment, and parser correctness.

Implicit solvent models and semiempirical methods are approximations. They are
often useful for triage and ensemble ranking, but they are not universal
substitutes for experimental data or high-level thermodynamic cycles.

## Boltzmann Populations

Boltzmann populations are derived from relative free energies. They must be
labeled with their comparison scope:

- Comparable within same formula/proton count.
- Approximate across different protonation/protomer states unless micro-pKa or
  proton chemical-potential corrections are available.

If a ranked output mixes protonation states, protomers, charges, or formulas,
DSVR should report those populations as approximate over generated candidates,
not as rigorous solution speciation.

## Enumeration Bias

Missing candidate states cannot be recovered by downstream ranking. If molscrub,
RDKit, or Auto3D omits a relevant protomer, tautomer, stereoisomer, or conformer,
CREST/xTB can only rank the candidates it receives.

## Dependency and Version Sensitivity

Results may change with tool versions, default parameters, CPU/GPU settings,
threading, random seeds, and parser behavior. DSVR should record versions,
command lines, configuration, input hashes, and logs for every run.
