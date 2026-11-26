# Model and verification

## Receptor component

The finite-state process is `U ⇄ C0 → C1 → … → Cn`. Association is `kon × a`, where `a` is a dimensionless effective antigen-contact input. Each bound state dissociates to `U` at `koff`; progression occurs at `kp`. Only `Cn` is active. Rates are in **min⁻¹**. These synthetic effective rates do not specify a measured CD19/CD22 affinity or antigen count per cell. Proofreading is a coarse-grained delay, not an explicit CD3ζ/ZAP70/MAPK network.

For row generator `Q`, probabilities satisfy `dp/dt = Qᵀp`. The active-state equilibrium is

`f(a) = kon·a / (kon·a + koff) × [kp / (kp + koff)]ⁿ`.

Matrix exponentiation and an exact finite-receptor Gillespie process give independent deterministic/stochastic references. The stochastic check uses 4,000 independent receptors initially unbound. Its simultaneous absolute-error bound is `sqrt(log(2m/δ)/(2N))`, with `m` states and `δ = 10⁻⁶`. Event-budget exhaustion raises an error. This check concerns receptor sampling, not biological uncertainty.

The population model uses receptor equilibrium. With the nominal `koff=0.2/min`, binding relaxation is on a minutes scale, much faster than the synthetic day-scale population dynamics. Fast equilibration and constant per-phenotype contact inputs are assumptions; ligand depletion, changing receptor density and spatial synapses are absent.

For arm activities `f19,f22`, CD19-only and CD22-only use the corresponding arm. Idealized independent-opportunity OR is `f19 + f22 − f19·f22`; AND is `f19·f22`. The four target phenotypes are `(++,+−,−+,−−)`, with present antigen assigned `a=1` and absent antigen `a=0`. Double-negative targets cannot activate or be killed. OR/AND are architecture assumptions, not inferred CTA101 logic. Independent-arm activities do not establish equal expression, costimulation or potency between constructs. A structural control replaces OR with `max(f19,f22)` to remove its double-positive bonus.

## Conditional population extension

Time is in **days**; `E`, `M` and `Tᵢ` are in **millions of cells**. Let `N=ΣTᵢ`, `A=ΣsᵢTᵢ` and `h=A/(H+A)`.

```
dE/dt  = (φh − ρ − αN)E + θAM
dM/dt  = εE − θAM − μM
dTᵢ/dt = rTᵢ(1 − N/K) − γsᵢETᵢ + (LᵀT)ᵢ
dC/dt  = qEh − dC·C
dJ/dt  = C
```

`C` is a cytokine proxy in arbitrary units; `J` is its time integral. It is not a measured cytokine concentration, CRS grade, toxicity threshold or safety endpoint. Proxy peaks are evaluated on the declared 0.1-day output grid. No cytokine-feedback effects are asserted.

| Parameter | Unit | Role |
|---|---|---|
| `r, φ, ρ, ε, μ` | day⁻¹ | Tumor growth; activated effector proliferation; total effector reduction; memory formation; memory death |
| `K, H` | million cells | Tumor capacity and activation half-saturation |
| `α, θ, γ` | (million cells·day)⁻¹ | Tumor-dependent suppression, memory recall and killing |
| `loss19, loss22` | day⁻¹ | Effective irreversible antigen-phenotype loss |
| `q` | arbitrary units/(million cells·day) | Proxy production |
| `dC` | day⁻¹ | Proxy clearance |

Memory formation is included in `ρ`, so `0 ≤ ε ≤ ρ`; natural effector death is `ρ−ε`. No-antigen activation is zero. `L` is a conservative row generator: loss of CD19 maps `++ → −+` and `+− → −−`; loss of CD22 maps `++ → +−` and `−+ → −−`. These transitions are effective phenotype changes, not an inference about genetic mutations, lineage switches or measured escape rates. Shared capacity prevents each phenotype from independently growing to `K`.

Two integrated counters record `Ḃ = ΣrTᵢ(1−N/K)` and `Ḋ = ΣγsᵢETᵢ`. The ledger `N(t)−N(0)−B(t)+D(t)=0` must close. With nonnegative states, all zero-state derivatives are nonnegative. Starting at `N≤K`, total tumor cannot cross above `K`, since loss transitions conserve cells and killing is nonnegative. The ledger alone does not validate the biological rate laws; analytical no-treatment and phenotype-switching references test the underlying dynamics separately.

No state is truncated to zero and no fractional population is called extinction. DOP853 is compared with tighter Radau integration. Nonfinite results, substantial negative populations, accounting errors and solver failure raise errors. Signed roundoff remains visible rather than being silently repaired.

## Published population reference

The independent `reference.py` implements [Barros et al. (2021), equations 1–3](https://doi.org/10.3390/cancers13122941), in cells and days:

```
dE/dt = (φ − ρ − αT)E + θTM
dM/dt = εE − θTM − μM
dT/dt = rT(1−bT) − γET
```

The **HDLM-2/CD123 mouse** column of Table 2 supplies all reference parameters, with attribution and units in `data/sources.json`. This model has antigen-independent effector proliferation, unlike the gated extension. It is a separate reference, not a direct calibration of the B-ALL model. No author source code is copied. The reported RAJI parameters are not used.

At coexistence, `E=r(1−bT)/γ` and `M=εE/(θT+μ)`. Substitution gives

`−αθT² + [(φ−ρ+ε)θ−αμ]T + (φ−ρ)μ = 0`.

Both positive roots and all three populations are compared with section 2.3. Table values give a lower tumor root of **1,091.197**, while the text prints **1,089**, a **0.202% discrepancy**. We retain the table parameters, report the printed-state residuals and use a 0.3% comparison tolerance. The computed roots themselves must satisfy the original ODEs to numerical precision. Untreated logistic growth and the tumor-free effector/memory solution provide additional independent references.

This verifies an implementation of published equations. It does **not** refit experimental curves or independently reproduce the paper's biological calibration. The original study assumed BLI units corresponded to cells without an established conversion and noted non-unique parameter choices.

## Separate adapter-binding component

[Lee et al. (2019)](https://doi.org/10.1038/s41467-019-10565-7) studied anti-fluorescein CAR-T cells, FITC-folate and folate-receptor-positive targets. Their reported cell-binding affinities are `KR=0.03 nM` and `KF=1 nM`; the isolated-antibody affinity is a different quantity. This architecture is not the clinical CD19/CD22 construct.

For specified **free** adapter `a`, CAR sites `R`, target sites `F` and ternary bridge `B`, assume equilibrium without cooperativity:

`RA=Ra/KR`, `FA=Fa/KF`, `B=RFa/(KR·KF)`.

Conserve `Rtot=R+RA+B`, `Ftot=F+FA+B`. Eliminating free sites yields `B=c(Rtot−B)(Ftot−B)`, where `c=a/[(KR+a)(KF+a)]`. The smaller physical root is evaluated in a numerically stable form. The analytic maximum occurs at `a=√(KR·KF)`. Zero adapter gives zero bridges; very high adapter suppresses bridging through separate saturation of the two sites (hook effect).

Synthetic receptor totals of `0.01 nM` each describe an illustrative contact compartment; their conversion from surface densities is unknown. Synthetic free-concentration jumps and a 60-minute free-exposure half-life illustrate exposure changes only. **No administered-dose conversion, full pharmacokinetics, competition intervention or receptor-residence kinetics is implemented.** Rapid equilibrium may fail for slow unbinding; declining free exposure does not establish rapid bridge rupture or toxicity reversal. This component is not coupled to the CD19/CD22 population study.

## Comparison and limits

All policies share the same initial tumor mixture, 30-day horizon and nominal population parameters. The exploratory sensitivity study uses 64 paired Latin-hypercube samples, independently varying growth, killing, recall, suppression and both antigen-loss rates by log-uniform factors from 0.5 to 2. These are declared design bounds, not fitted distributions, confidence intervals or virtual patients. Every parameter draw and result is retained. Spearman coefficients summarize these chosen scenarios without causal or clinical interpretation.

## Modest PK/PD extension: exposure and response tradeoffs

The question motivated by the review is **whether timing changes tumor suppression and inflammatory activation when total adapter input is held fixed, and whether that comparison survives unmeasured kinetic assumptions**. This is a scenario-comparison contribution built around published components, rather than a new biological mechanism or a refit of an existing population model.

The separate FR-positive adapter arm uses one specified systemic concentration `a`:

`da/dt = u(t) − ke·a`, with `ke = log(2)·1440 / half_life_minutes`.

Boluses increment `a` exactly at their declared times; continuous infusions specify `u` in nM/day. All active schedules deliver **0.7 nM of concentration-equivalent input** over seven days and are followed to day 14. This is an assumed fixed-volume exposure input, not nmol/kg or a clinical administered dose. The nominal half-life is 60 minutes. The total-input/elimination ledger must satisfy `a + eliminated = input`; `AUC = eliminated/ke`. Withdrawal permits nearly complete washout, giving equal free-adapter AUC across the active schedules. Equal free AUC need not give equal engagement or PD response.

Binding now uses finite mass-action **kinetics**, not instantaneous equilibrium. Let `x=RA`, `y=FA`, `z=B`, `R=Rtot−x−z` and `F=Ftot−y−z`. With rates in day units:

```
vR = konR·a·R − koffR·x
vF = konF·a·F − koffF·y
wF = konF·x·F − koffF·z
wR = konR·y·R − koffR·z
dx/dt = vR − wF
dy/dt = vF − wR
dz/dt = wF + wR
```

`kon = koff/Kd` preserves each measured equilibrium affinity while allowing unknown residence times to vary. Both binding orders agree with the equilibrium component and conserve local receptor pools. The representative contact compartment is assumed small relative to the systemic reservoir, so its binding does not deplete plasma adapter. Plasma and local free ligand are assumed equal; tissue distribution and transport are absent. These assumptions need testing before interpreting real exposure.

Relative engagement is `g=z/(Bhalf+z)`. The cellular model uses `sFR+=g` and `sFR−=0`, with only FR-positive and FR-negative tumor pools. CD19/CD22 loss transitions are disabled. The initial mixture is 99% FR-positive and 1% FR-negative. Other population rates are synthetic, matched across schedules; they are not transferred as measured rates from a B-ALL trial. `Bhalf` and both dissociation rates are unmeasured. Cellular expansion/persistence and engagement-dependent tumor/proxy responses provide a limited PK/PD connection, without an ADMET claim.

Compare daily boluses, four split boluses/day, ramped boluses and continuous input. All input boundaries split Radau integration exactly; sampled bolus states are post-input. The total tumor, systemic adapter ledger and finite receptor bounds are checked at adaptive solver steps and output times. Tightened tolerances check PD endpoint convergence. Cytokine peaks use a 0.02-day observation grid.

Sixteen paired sensitivity scenarios vary exposure half-life (30–89.9 min), CAR and target off-rates (0.002–0.2/min and 0.02–2/min), and bridge half-activation, killing and proliferation (0.5–2× nominal). Kds stay fixed. The half-life interval is illustrative and does not calibrate mouse tumor kinetics from a human circulation bound. A schedule is nondominated if no other tested schedule has both lower-or-equal final tumor and lower-or-equal proxy, with at least one strict improvement. Separate comparisons use peak and integrated proxy. Frontier counts describe these selected scenarios, not clinical success probabilities or a globally optimized regimen.

The model may expose a tradeoff but cannot explain or establish the review's observed improvement in antitumor potency alongside reduced toxicity. Exhaustion, host cytokine amplification and tissue injury are not represented. Unresolved mechanisms are not added simply to force agreement with the paper. Raw assay observations would be needed to calibrate these responses and discriminate mechanisms; published affinities and a numerical reference are not enough.

The baseline excludes normal B cells, host immunity, allogeneic rejection, exhaustion states, spatial trafficking, manufacturing and stochastic tumor extinction. Normal B cells also express the target antigens, so tumor selectivity or clinical safety cannot be inferred. TRAC/CD52 edits and their effects are not simulated; no genotype-to-phenotype mapping is borrowed from the Cas9 prototype. DFT and quantum algorithms are not required to solve this classical baseline.
