# Biological modeling specification

This specification combines the three resource proposals with the biological-modeling scope. Bioinformatics is outside this implementation. Work proceeds one model at a time, starting with the CRISPR-Cas9 kinetic baseline.

## Three themes and four models

| Theme | Individual model | Primary observables |
|---|---|---|
| Molecular and cellular dynamics | CRISPR-Cas9 editing | Binding/cleavage time courses, first-passage times and conditional repair outcomes. |
| Molecular and cellular dynamics | Dual-target CAR-T signaling | Antigen-dependent signaling, tumor and T-cell dynamics, and cytokine proxies. |
| Anesthesia and awareness | Anesthetic-modulated neural dynamics | Oscillation spectra, synchrony, perturbation recovery and independently defined state labels. |
| Mitochondrial electron transfer | Complex III -> cytochrome c -> complex IV | Redox populations, docking/transport contributions, transfer rates and flux. |

Resource mapping: `Biomolecular_Simulation.md` specifies Cas9; `Molecular_Signaling_Cascades.md` specifies CAR-T; `Anesthesia_Awareness.md` specifies the neural study. The electron-transfer scope is mitochondrial complex III → cytochrome c → complex IV, distinct from a multi-heme protein wire.

## First implementation and acceptance

The first Cas9 milestone is a transparent probability model: unbound <-> PAM-bound <-> R-loop -> cleaved, followed by precise repair, indel or HDR branches. Free RNP is initially a fixed reservoir. Terminal repair makes this a single-cleavage model; delivery, re-cutting and cellular context must not be implied by this baseline.

Acceptance requires nonnegative transition rates, conserved target probability, explicit units, reproducible local random streams, matrix-exponential/ODE agreement and Gillespie agreement within a declared statistical bound. Binding-only, irreversible-cascade and repair-only limits must match closed-form formulas. Failed convergence or exhausted sampling budgets must not produce a successful report.

Synthetic parameter examples establish numerical behavior only. Assay calibration is a separate milestone: bind each dataset to its nuclease, guide/target sequences, concentrations, time units, observation definition, source and licence. Separate binding, cleavage and repair readouts. Assess parameter identifiability and hold out guides/assays rather than randomly separating correlated measurements. Do not infer indel spectra or T-cell editing accuracy from a cleavage-only fit.

## Subsequent classical models

CAR-T follows after Cas9. Start with receptor/antigen kinetics and an explicit AND/OR dual-target rule, then connect signaling to population dynamics. Validate zero-antigen controls, state positivity, units and parameter sensitivity. Cytokine output is a model proxy until clinical-outcome validation exists. Connecting edit outcomes to signaling requires an independently supported genotype-to-phenotype mapping.

For anesthesia, a propofol-inspired thalamocortical model is a proposed tractable starting point, subject to approval. Reproduce selected dynamical references before fitting EEG features. Use subject/session-separated evaluation and distinguish observable response labels from awareness itself. A simulated bifurcation or classifier score does not establish a theory or measurement of consciousness.

For mitochondrial transport, separate cytochrome c diffusion/docking from intracomplex electron transfer. Resolve the complex-III donor interface and cytochrome-c donation to complex-IV CuA. Begin with thermodynamically consistent rate/population models and published kinetic references. Do not model mobile cytochrome c as a permanently coherent wire or imply that a reduced transfer model reproduces oxygen reduction/proton pumping.

## DFT and FeNNix

DFT/QM-MM is optional where electronic energetics are needed, especially redox states, reorganization and donor/acceptor coupling in the electron-transfer study. State geometry, charge, spin, environment and functional/basis assumptions. Whole-cell and neural dynamics do not require DFT by default.

FeNNol is the molecular-model engine; FeNNix-Bio1 is a particular set of interatomic-potential weights. Its model card identifies molecular dynamics and geometry/binding applications and an Academic Software License. The locally inspected Bio1 reference-energy table has neither Fe nor Cu entries. Therefore the model is not assumed suitable for heme/CuA redox chemistry, reactive Cas9 catalysis or electronic couplings. Any use needs explicit chemical-domain, property and licence checks. It is not a substitute for validated redox electronic-structure calculations.

## Reuse and quantum extensions

Reuse classical optimization, stochastic sampling, leakage-resistant evaluation, Hamiltonian/reference comparisons and resource-accounting patterns from the existing repositories. Preserve attribution and scientific assumptions; do not copy chemical scaffold rules into guide-RNA splits or treat nonlinear biological ODEs as direct HHL inputs.

Possible later experiments include uncertainty-estimation or optimization circuits, small electronic-structure/dynamics benchmarks, and quantum feature maps with matched classical controls. Report classical cost, encoding/preparation, condition numbers, truncation, measurements and readout cost. Physical electron tunneling is not evidence of a quantum-computing speedup. Advantages are hypotheses to test, not acceptance criteria to assume.

## Primary starting references

- [Cas9 thermodynamic/kinetic model, 2016](https://doi.org/10.1371/journal.pcbi.1004724).
- [Cas9 unwinding/cleavage kinetics, 2018](https://doi.org/10.1016/j.celrep.2017.12.041).
- [SpCas9 kinetic model, 2022](https://doi.org/10.1038/s41467-022-28994-2).
- [Mechanistic CAR-T signaling and heterogeneity](https://pmc.ncbi.nlm.nih.gov/articles/PMC7467855/).
- [Propofol thalamocortical alpha model](https://doi.org/10.1073/pnas.1017069108).
- [Respiratory complex III/IV transport kinetics](https://pubmed.ncbi.nlm.nih.gov/32201307/).
- [Cytochrome c/CuA electron transfer](https://www.pnas.org/content/pnas/109/43/17348.full.pdf).
- [FeNNix-Bio1 model card and licence](https://huggingface.co/fennol-tools/FeNNix-Bio1).
