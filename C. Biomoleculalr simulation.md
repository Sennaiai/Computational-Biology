# Project 2 — CRISPR‑Cas9 in T Cells: Mechanism, Off‑Target Risk, and Probabilistic Outcomes

      ### Overview

Design or import gRNAs for a (adaptable) target locus (CD19/CD22), predict calibrated on/off‑target cut probabilities with T‑cell epigenetic features, estimate indel spectra, aggregate wet‑lab‑like scenario probabilities, and export a JSON handoff for (Project 3) CAR‑T signaling simulations

      ### Project Steps

1. **Design/Import & Rank**: Enumerate or ingest gRNAs -> target genome scan for off‑targets -> baseline scoring -> rank top‑k

2. **Cut Probability Model**: Mechanistic prior (penalties + simple R‑loop barrier) -> ML head -> probability calibration

3. **Edit Outcome Model**: Microhomology‑aware indel spectrum -> frameshift vs in‑frame distribution

4. **Scenario Aggregation**: Monte‑Carlo sampling -> natural‑language summaries per gRNA (on/off‑target outcome mix)

5. **Visualization & Export**: Design summary, per‑site probabilities, repair outcome plots -> JSON handoff for Project- 3

        ### Methods

### Step 1 — Data Access -> Preprocessing -> Feature Extraction

* **Access**: Cas‑OFFinder/CRISPRitz for off‑target enumeration (≤3–4 mismatches/bulges). Reference: hg38
  
* **Organization**: `datasets/{raw,features}/`, `designs/`, `reports/`, `models/`

* **Feature Engineering**
      
      * *Sequence/Context*: PAM type, seed mismatch pattern, GC%, homopolymers, RNAfold hairpin flag
  
      * *Biophysical*: position‑dependent mismatch penalties; simple R‑loop energy proxy
  
      * *T‑cell Epigenetics*: ATAC/DNase signal, H3K27ac/H3K4me3/H3K27me3, CpG methylation
  
      * *Repair Bias*: microhomology length/content
 
  
              * **Outputs**: `designs/candidates.csv`, `datasets/features/features.parquet`


### Step 2 — Cut Probability (Binding -> Cleavage) & Calibration

* **Baseline**: logistic/GBM on engineered features
  
* **Calibration**: temperature scaling (ECE/NLL check)/conformal risk tiers
  
* **Evaluation**: PR‑AUC, top‑k recall, calibration curves; ablation ± epigenetics
  
      * **Outputs**: `models/cut_probability.pkl`, per‑site probability tables


### Step 3 — Edit Outcome Distribution (Repair)

* **Method**: inDelphi/FORECasT‑style heuristic using microhomology + context
  
* **Outputs**: frameshift %, in‑frame %, top indel probabilities per site
  

### Step 4 — Scenario Aggregation & Reporting

* **Monte‑Carlo** over on‑target + top‑K off‑targets -> probability of practical outcomes
  
* **Reports**: `reports/top5_scenarios.html` with per‑gRNA summary (e.g. “≥80% frameshift with 64% probability; off‑target chr7 >5% in 7% of runs”)

### Step 5 — Visualization & Export

* **Visuals**: locus lollipop map of predicted cuts; seed mismatch heatmaps; indel bars & frameshift pies; simple R‑loop schematic
  
* **Handoff (to Project 3)**: JSON

```json
{
  "t_cell_engineering": {
    "trac_knockin_success": 0.68,
    "pd1_knockout_fraction": 0.74,
    "offtarget_burden_index": 0.12
  },
  "notes": "T‑cell epi‑aware; calibrated"
}
```

            ## Datasets & Sources (targeted)

* **Off‑target assays**: GUIDE‑seq, CIRCLE‑seq, SITE‑seq, DISCOVER‑seq (for training/validation)

  
* **Genomes**: hg38; gene annotation for CD19/CD22

  
* **T‑cell epigenetics**: ENCODE/Roadmap ATAC/DNase, histone marks; WGBS methylation

  
* **Baselines**: CFD/MIT scores for benchmarking

## Expected Deliverables

* `designs/candidates.csv`, `datasets/features/features.parquet`
* `models/cut_probability.pkl`
* `reports/design_summary.html`, `reports/top5_scenarios.html`
* `app/` (Streamlit/Dash)
* **Handoff JSON** for CAR‑T model



      ### Future Directions (T‑cell focused)

* **Quantum‑Enhanced Prediction**: hybrid QML kernel or CV‑QRC re‑ranking; optional VQE/QAOA for constrained gRNA selection
  
* **Editing Persistence**: stochastic division model for T‑cell expansion (Cas9 decay, re‑cut, mosaicism) → updated priors for CAR‑T


### Repo Skeleton

```
crispr-design/
  README.md
  config.yaml
  data/  designs/
  datasets/raw/  datasets/features/
  models/  mechanism/
  notebooks/
    01_design.ipynb
    02_cut_probability_and_calibration.ipynb
    03_scenarios_and_reports.ipynb
  reports/  app/
```

## Milestones (MVP -> V1)

* **MVP**: Step 1 + Step 2 (logistic + calibration) + Step 4 (scenarios) → design & scenario reports for top‑5 gRNAs
  
* **V1**: add epigenetic features + ablation; repair outcome plots; optional UI and JSON handoff

* **Stretch**: deep model + conformal intervals; optimizer; lightweight animations; alternate nucleases

