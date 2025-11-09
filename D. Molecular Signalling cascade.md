### Project 3 — CAR‑T Cell Signaling Cascade: Efficacy, Safety, and System Dynamics

### Overview

Model and simulate the intracellular signaling cascade of dual‑target CAR‑T cells (CD19/CD22) to analyze activation, proliferation, cytotoxicity, and safety (cytokine release). The project integrates mechanistic ODE modeling with stochastic simulation and ML overlays to predict therapy efficacy and toxicity, linking directly to the T‑cell engineering data from Project 2.

## Project Steps

1. **Model Definition**: Construct ODE‑based signaling and interaction network for dual‑antigen (CD19/CD22) CAR‑T activation.
2. **Parameterization**: Set kinetic parameters using literature and priors from CRISPR editing outcomes (TRAC knock‑in, PD‑1 KO, off‑target burden).
3. **Simulation & Sensitivity**: Run deterministic/stochastic simulations; analyze sensitivity and stability.
4. **Efficacy & Safety Metrics**: Quantify tumor clearance, CAR‑T expansion, cytokine surges (CRS proxy).
5. **Visualization & Reporting**: Generate activation maps, cytokine timecourses, and Pareto plots (efficacy vs toxicity).

## Methods (Proposed)

### Step 1 — Model Structure & Equations

* **Modules**:

  * **Receptor Binding**: CAR engagement with CD19/CD22 → signal initiation.
  * **Signal Transduction**: Lck/ZAP70 → PI3K/AKT, MAPK, NF‑κB pathways.
  * **Cellular Responses**: cytokine production, proliferation, cytotoxic activity, exhaustion.
  * **Tumor Dynamics**: predator–prey equations; tumor antigen loss (escape).
* **Equations**: System of ODEs (mass‑action or Michaelis–Menten forms) representing activation and feedback loops.

### Step 2 — Parameterization

* **Sources**: Literature values, in‑vitro CAR‑T assays, and JSON priors from Project 2 (TRAC success, PD‑1 KO fraction, off‑target burden).
* **Strategy**:

  * Baseline parameters for Lck/ZAP70 and downstream cascades from existing models.
  * Scale activation/persistence based on CRISPR priors: PD‑1 KO → ↑activation/persistence; off‑target burden → ↓fitness.
  * Dual‑antigen logic: OR/AND gating (CD19/CD22 density × CAR affinity).

### Step 3 — Simulation & Sensitivity Analysis

* **Simulation**: deterministic ODE solver (SciPy.integrate) + stochastic Gillespie variant for noise.
* **Sensitivity**: Sobol/FAST methods for parameter influence.
* **Outputs**: activation curves, tumor burden over time, cytokine levels, CAR‑T counts.

### Step 4 — Efficacy & Safety Quantification

* **Metrics**:

  * Efficacy → tumor burden AUC↓, time‑to‑clearance, CAR‑T expansion/persistence.
  * Safety → cytokine peak (IL‑6, IFN‑γ), activation fraction > threshold.
* **Evaluation**: trade‑off frontiers (Pareto: efficacy ↑ vs toxicity ↓), compare single‑ vs dual‑target.

### Step 5 — Visualization & Reporting

* **Plots**: timecourses (activation, cytokines, tumor burden), heatmaps (affinity/density sweeps), phase diagrams.
* **Reports**: `reports/cascade_summary.html` with key metrics and sensitivity tables.
* **Optional App**: sliders for antigen density, affinity, or CRISPR priors → interactive visualization.

## Datasets & Sources

* **Biochemical Parameters**: published CAR‑T signaling models (e.g., PI3K/AKT, NF‑κB kinetics).
* **Experimental References**: cytokine release assays, cell killing curves (B‑ALL context).
* **Inputs from Project 2**: JSON of TRAC knock‑in, PD‑1 KO, off‑target burden indices.

## Expected Deliverables

* `models/cascade_equations.py` (ODE implementation)
* `datasets/parameters.yaml` (kinetic constants, priors)
* `reports/cascade_summary.html`
* `notebooks/03_sensitivity_and_tradeoffs.ipynb`
* `app/` (optional visualization dashboard)

## Future Directions

* **Quantum‑Enhanced Optimization**: Use QAOA or hybrid quantum kernels to optimize CAR design (affinity, costimulatory domains) under safety constraints.
* **Integration with Clinical Data**: fit to patient cytokine/time‑course datasets for transfer learning.
* **Multi‑Scale Modeling**: extend from intracellular signaling → tissue‑level tumor–immune dynamics.

## Repo Skeleton

```
car-t-simulation/
  README.md
  config.yaml
  data/  datasets/
  models/  notebooks/
    01_model_definition.ipynb
    02_parameterization_and_simulation.ipynb
    03_sensitivity_and_tradeoffs.ipynb
  reports/  app/
```

## Milestones (MVP → V1)

* **MVP**: Build ODE model + baseline simulation + efficacy/safety metrics for CD19/CD22 dual targeting.
* **V1**: Add stochastic simulation, sensitivity analysis, and interactive plots.
* **Stretch**: quantum‑enhanced optimization; fit to experimental data; tissue‑level coupling.

