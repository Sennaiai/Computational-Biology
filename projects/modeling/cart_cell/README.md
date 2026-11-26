# CAR-T: tumor control and activation tradeoffs

A mechanistic simulation motivated by my 2023 CAR-T analytical review, coauthored with Ryan Owchar: **how do recognition strategy and controlled exposure affect tumor suppression, antigen escape and inflammatory activation?** The contribution is a reproducible comparison tool with explicit uncertainty, rather than a claim of a new treatment mechanism.

Two complementary experiments keep their biological architectures separate:

- **CD19/CD22 recognition:** single-target, idealized OR and AND rules, finite receptor proofreading, effector/memory dynamics and four antigen phenotypes.
- **Adapter exposure:** a separate anti-fluorescein/FITC-folate/folate-receptor model compares bolus, split, ramped and continuous input. A modest PK/PD layer connects exposure and binding kinetics to cellular response; active schedules have equal total input and free-adapter AUC.

## Evidence and results

Published mouse-model equations and reported equilibria provide an independent numerical reference. Published adapter-binding affinities anchor a separate component. **Population-extension and exposure-response parameters remain illustrative; this is not assay-calibrated B-ALL or clinical risk prediction.**

The nominal antigen study gives day-30 tumor burdens of **12.29 million cells (CD19), 12.97 (CD22), 1.11 (OR) and 22.12 (AND)** from one million initial tumor cells. OR reaches a transient minimum of 0.072 million, then double-negative escape drives regrowth. These are conditional simulations, not treatment success rates.

The exposure-matched adapter study shows a tradeoff:

| Input schedule | Day-14 tumor, million cells | Peak activation proxy, arbitrary units |
|---|---:|---:|
| No adapter | 7.62 | 0 |
| Daily bolus | 6.33 | 0.0266 |
| Split bolus | 3.48 | 0.1014 |
| Ramped bolus | 6.71 | 0.0195 |
| Continuous | 2.64 | 0.1303 |

Lower tumor burden does not automatically mean a better overall tradeoff. Paired sensitivity analyses retain all parameter draws and assess dominance under both peak and integrated proxy. Unknown binding residence times are varied while measured affinities stay fixed. No clinical risk probability, safe dose or preferred patient regimen is inferred.

![Matched-input comparison](results/pkpd_tradeoffs.svg)

## Run and inspect

Python 3.11; from the repository root:

```bash
python -m pip install -e 'projects/modeling/cart_cell[notebook]'
python -m unittest discover -s tests/modeling/cart_cell -v
python projects/modeling/cart_cell/scripts/execute_notebook.py --output-dir artifacts/cart_cell
```

[Executed notebook](CAR_T_Tradeoffs.ipynb) · [Equations, units and checks](MODEL.md) · [Sources](src/cart_dynamics/data/sources.json) · [Recognition results](results/cart_results.json) · [Exposure comparison](results/pkpd_results.json)

Analytical controls, exact receptor sampling, independent solvers, input boundaries and receptor/tumor/exposure accounting are tested. Scoped CI executes this study without rerunning completed projects. Results record configurations, seeds, source hashes, numerical residuals and runtime. The published table/text equilibrium discrepancy (0.20%) is retained and documented.

Normal B cells, host immunity, exhaustion, tissue distribution and clinical cytokine feedback are absent. The proxy is not CRS; universal-cell edits and combined CD19/CD22 adapter constructs are not simulated. The review supplies the motivation; this prototype does not establish clinical efficacy or safety.
