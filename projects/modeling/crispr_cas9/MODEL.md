# Model and evidence

## Intracellular system

The simulator tracks free cytoplasmic and nuclear Cas9, two nuclear guide RNAs and their active nuclear complexes as integer counts. Production is an effective protein/RNA source; transcription, translation and guide processing are lumped. Import/export, degradation, reversible RNP assembly and target-bound complexes are explicit.

For nuclear volume V in fL, one nM corresponds to `6.02214076e23 × 1e-24 × V` molecules. RNP assembly has propensity `k × N_Cas9 × N_gRNA / molecules_per_nM`. DNA binding consumes a free RNP and product release returns it. Complex degradation removes its Cas9 and guide together. Births, losses, free pools and target-bound pools reconcile exactly across the founder lineage.

Each locus has two alleles. Their states are:

```text
intact <-> PAM-bound <-> R-loop -> cleaved/Cas9-bound -> released DSB
  ^                                                       |
  +-------------------- precise repair ------------------+
                                                          |
                         three resistant mutation categories
```

Two categories represent hypothetical disruptive alleles and one other resistant allele. Their frequencies are scenario assumptions, not sequence-resolved indel predictions. Disruptive labels do not establish functional knockout. Intended loci share Cas9, while one illustrative off-target uses guide 1. Genome-wide off-target search, inter-locus deletions, chromatin and cell death are absent.

## Expression and lineage

Default nucleus/cytoplasm birth volumes are 30/150 fL. Volume grows exponentially to twice birth volume over an assumed eight-hour cycle. Production is piecewise constant per volume; first-order reactions depend on counts and bimolecular reactions on concentration. Exact stochastic thinning uses endpoint bounds over each expression/growth segment. There is no simulation time-step approximation. Growth-dependent births are tested against an independent nonhomogeneous Poisson reference.

Division halves volume and partitions free molecules binomially. Pre-cleavage complexes dissociate before DNA replication. Unrepaired breaks block division and pause growth until repaired. Both daughters inherit the diploid mutation states. This strong, uncalibrated checkpoint can dominate lineage size and mosaicism. The 24-hour study permits at most three divisions and does not model whole-plant development, germline transmission or seed segregation.

Expression windows represent developmental timing or controlled exposure. Named promoters are not assigned invented production rates. Early and late cases both use a 12-hour production window; this does not equalize total protein expenditure across growing lineages. Promoter/terminator context requires independent evidence: [Wang et al., 2015](https://doi.org/10.1186/s13059-015-0715-0).

## Published-assay fit

The [source manifest](src/crispr_kinetics/data/plant_repair_sources.json) identifies the paper, pinned repository commit, checksums and Table 1 repair references from [Ben Tov, Mafessoni et al., 2024](https://doi.org/10.1038/s41467-024-49410-x). The count table is the authors' public Psy1 example, not asserted identical to the final paper's Source Data workbook. No author R code is redistributed.

The repository README conflicts with its preparation code. We follow `csv2timecourse.R` and the plotting function in `likelihood_functions.R`: **y1 = intact, y2 = indels, y3 = direct DSB, y4 = processed DSB**. Counts and the independent control-error matrix remain unchanged.

Our reduced latent model is `I <-> B -> M`, with constant cutting, precise repair and mutation repair. A fitted observation fraction splits B into direct/processed DSB categories before applying the error matrix. Initial DSB fraction is fitted from training observations. This simplification omits the publication's time-dependent induction model and does not reproduce its full parameterization.

Eight fixed-seed optimization starts fit observations through 48 h. Both 72 h replicates are reserved for extrapolation evaluation. Rate bounds are cut/mutation `exp(-12)..1 h^-1` and precise repair `exp(-12)..8 h^-1`; an ablation fixes precise repair to zero. Observation fractions lie inside `(0,1)` and initial DSB inside `(0,0.2]`. A training-only likelihood profile assesses precise-repair sensitivity. Conditional read-count likelihoods do not fully account for replicate overdispersion and are not biological confidence intervals. Finite multistart optimization cannot certify a global optimum.

Example-fitted repair parameters and published tomato point estimates enter separate cellular scenarios. No expression, embryo, Arabidopsis or ABCG18/19 measurements were fitted. Public kinetic studies provide context for the target-engagement layer, whose current rates remain illustrative: [Gong et al., 2018](https://doi.org/10.1016/j.celrep.2017.12.041), [Eslami-Mossallam et al., 2022](https://doi.org/10.1038/s41467-022-28994-2).

## Readouts and uncertainty

Readouts separate free nuclear RNP, mutant alleles, a double-biallelic disruption proxy, selected-off-target mutations and unresolved breaks. Cell fractions first average within each founder, then across independent founders. Descendants are correlated and never treated as independent replicates.

The mosaicism lower event requires two fully resolved terminal cells with different synthetic genotypes. The upper event also includes multi-cell lineages with unresolved breaks. Single-cell or incompletely repaired lineages cannot establish non-mosaic plants. Per-quantity 95% Hoeffding bounds describe Monte Carlo sampling only; they remain nonzero after observing no events and are not simultaneous across quantities. Expression/repair cases are sensitivity comparisons, not a biological posterior or a validated ranking of constructs.

For custom inputs, use `dataclasses.replace(CellParameters(), ...)` and `simulate_founder(parameters, loci=..., seed=...)`. Supplied locus rates and expression profiles require validation in their own experimental context before interpreting absolute outcomes.
