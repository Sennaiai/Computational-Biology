# Overview

This project investigates how **self-awareness** emerges from the **physical and informational organization** of living 
systems.
It begins from the premise that **consciousness** is matter’s capacity to process information, while **awareness** is 
the system’s recursive modeling of itself.

I hypothesize that **self-awareness arises near a critical point between order and chaos**, where information flow, 
memory, and stability are co-optimized. **Anesthesia** provides a natural testbed: awareness fades as the brain departs 
from this critical balance and returns as it is restored.

Through **computational modeling** and **multimodal data analysis**, this project aims to identify **quantitative 
markers of transitions** between awareness and unawareness—reframing consciousness as a phenomenon of **self-organized
critical homeostasis** that bridges **molecular biology, neuroscience, and quantum information**.

**– Link to Blog**

---

## Project Steps

A multi-stage investigation combining:

1. **Data Analysis:** Statistical and information-theoretic characterization of EEG and fMRI data during anesthesia.
2. **Simulation:** Network-based and dynamical-systems modeling of brain criticality (e.g., coupled oscillators, Ising,
or neural field models).
3. **Machine Learning:** Classification and prediction of awareness states using entropy, mutual information, and 
criticality metrics.
4. **Quantum/Information Integration:** Exploratory modeling of integrated information and coherence using **Integrated 
Information Theory (IIT)** and **quantum-inspired frameworks**.

---

## Methods (Proposed)

### Step 1 — Data Access → Preprocessing → Feature Extraction

#### a. Access

**i. EEG Pipeline**

* **Source:** PhysioNet anesthesia EEG datasets (MNE-compatible: `.edf`, `.set`, or `.fif`)
* **Organization:** `datasets/eeg/raw/`, `datasets/eeg/derivatives/`, `notebooks/eeg/`

**ii. fMRI Pipeline**

* **Source:** OpenNeuro propofol sedation datasets
* **Organization:** BIDS layout in `datasets/fmri/raw/` and derivatives in `datasets/fmri/derivatives/`

---

#### b. Preprocessing

**i. EEG Pipeline**

* **Re-reference:** Average mastoids or common average
* **Filters:** 0.3–45 Hz band-pass; 50/60 Hz notch filter
* **Downsampling:** e.g., 250 Hz (if needed)
* **Artifact Removal:** Blink and muscle artifact correction via ICA; reject and interpolate bad channels
* **Epoching:** Fixed windows (e.g., 10 s) with overlap for stable estimates
* **Labeling:** Derived from metadata or provisional labels based on depth-of-anesthesia notes

**ii. fMRI Pipeline**

* **Preprocessing:** fMRIPrep-style (motion correction, slice timing, susceptibility distortion correction, 
normalization to MNI)
* **Nuisance Regression:** Motion, WM, CSF; band-pass filter 0.01–0.1 Hz; scrub high-motion volumes
* **Parcellation:** AAL, Schaefer-200, or Glasser atlas; extract regional time series

---

#### c. Feature Extraction

**i. EEG Per Epoch**

**Spectral and Aperiodic Features**

* Absolute and relative power in delta (0.5–4 Hz), theta (4–8 Hz), alpha (8–14 Hz), beta (14–30 Hz) bands
* Frontal alpha peak frequency and power
* Aperiodic 1/f exponent and offset (FOOOF)

**Complexity and Entropy**

* Multiscale entropy
* Permutation entropy
* Lempel–Ziv complexity (binary coarse graining)
* Sample entropy

**Criticality Proxies**

* Burst-suppression ratio (time below amplitude threshold)
* Detrended fluctuation analysis (DFA α)
* Autocorrelation at lag-1, rolling variance (critical slowing down)
* Neuronal “avalanches”: event detection via z-threshold, size/duration distributions, power-law exponents, branching
ratio ( m )

**Connectivity Features (Optional)**

* Phase-locking value (PLV) or weighted phase-lag index (wPLI) in alpha and delta bands
* Transfer entropy or phase-slope index for directionality

**Output**

* Tidy feature table per epoch (one row per epoch with features and state label)
* Save to: `datasets/eeg/derivatives/features.parquet`

---

**ii. fMRI**

**Functional Connectivity**

* Static FC: Pearson and partial correlations (Fisher-z transformed)
* Spectral FC: Coherence in low-frequency bands
* Dynamic FC: Sliding windows and state clustering

**Effective Connectivity**

* Granger causality on prewhitened time series
* Transfer entropy on downsampled time series

**Graph Measures**

* Global efficiency
* Modularity
* Participation coefficient
* Clustering coefficient
* Hubness measures

**Complexity and Criticality Proxies**

* Lempel–Ziv complexity on regional or principal component time series
* Multiscale entropy
* Variance and lag-1 autocorrelation as markers of critical slowing

**Output**

* Per-scan feature tables and FC matrices saved to `datasets/fmri/derivatives/`

---

### Step 2 — Criticality and Complexity Analysis → State Mapping → Initial Model

#### A. States

1. Deep Anesthesia
2. Intermediate Anesthesia
3. Sedated
4. Intermediate Aware
5. Aware

---

#### B. Expected Signatures per State (Guide Rails)

| Feature Family                           | Deep Anesthesia   | Intermediate Anesthesia | Sedated   | Intermediate Aware | Aware         |
| ---------------------------------------- | ----------------- | ----------------------- | --------- | ------------------ | ------------- |
| **Burst-suppression ratio**              | High              | Low–moderate            | Near zero | Zero               | Zero          |
| **Frontal alpha power**                  | Fragmented        | Low/slow                | Weak      | Moderate           | Stronger      |
| **Delta power**                          | Very high         | High                    | Elevated  | Mild               | Baseline      |
| **1/f aperiodic exponent**               | Steeper           | Steep                   | Moderate  | Mild               | Flatter       |
| **Lempel–Ziv complexity**                | Lowest            | Low                     | Moderate  | Higher             | Highest       |
| **Multiscale entropy**                   | Low               | Low–moderate            | Moderate  | Higher             | High          |
| **DFA α**                                | Far from critical | Closer                  | Closer    | Near critical      | Near critical |
| **Autocorr (lag-1, variance)**           | Transition peak   | Elevated                | Moderate  | Mild               | Baseline      |
| **Avalanche exponents / ( m )**          | ( m < 1 )         | Approaching 1           | Near 1    | ≈1                 | ≈1            |
| **Functional connectivity (efficiency)** | Low               | Low–moderate            | Moderate  | Higher             | Highest       |

**Goal:** Learn quantitative **boundaries** between these states empirically from data.

---

#### C. Criticality Tests

* **Power-law Fits:** Apply Clauset method on avalanche size/duration distributions; compare exponents to theoretical 
near-critical values and log-normal alternatives.
* **Branching Ratio ( m ):** Estimate from avalanche trees or point processes; near 1 indicates criticality.
* **Critical Slowing Down:** Identify rising variance and lag-1 autocorrelation near loss or recovery of responsiveness.
* **DFA:** Apply to amplitude envelopes and principal components (values near 0.5 = uncorrelated; higher = long-range 
dependence).
* **Aperiodic Exponent:** Derive from FOOOF as proxy for excitation-inhibition balance.

---

#### D. State Model

* Create feature vector per epoch or scan:
  `[spectral, aperiodic, entropy, LZC, DFA, autocorr, avalanche, FC graph metrics]`
* Fit a **Hidden Markov Model (HMM)** or **Bayesian state-space model** with 5 latent states.
* Initialize with k-means centroids for stability.
* Constrain by anchoring one state to **burst suppression** and another to **fully awake**; allow remaining three to be 
learned.
* **Output:** State trajectories over time and centroids defining the five levels.

---

#### E. Evaluation

* **Internal:** Cluster silhouette, temporal dwell times, and transition probabilities vs. clinical phases.
* **External:** Accuracy vs. behavioral responsiveness (if available); replication across datasets.
* **Interpretability:** SHAP values or permutation importance to reveal defining features per state.

---

### Step 3 — Modeling, Simulation, and Machine Learning

* Simulate critical transitions with **Kuramoto** or **Wilson–Cowan** networks and fit to observed features.
* Train classifiers to predict the five states from extracted features (logistic regression, SVM, GBM, or neural 
networks).

---

### Step 4 — IIT Integration

* Compute toy **Φ** on reduced networks (using pyPhi or custom implementation).
* Relate integrated information to state indices and include Φ as a feature in the overall model.

---

### Step 5 — Quantum Extensions

* **CV-compatible encodings** for spectral and aperiodic features.
* **Quantum Machine Learning (QML)** hybrid kernels as baselines alongside classical models.
* Long-term: simulate **drug mechanism effects** as parameter shifts in network models to explore coherence and 
decoherence processes.

---

## References

* Brown et al., *PNAS* (2010)
* Purdon et al., *Science* (2013)
* Schartner et al., *Neuroscience of Consciousness* (2015)
* Casali et al., *PNAS* (2013) — PCI Study
* Shew & Plenz, *J. Neurosci.* (2013)
* Tagliazucchi et al., *Frontiers in Systems Neuroscience* (2016)
* Barttfeld et al., *PNAS* (2015)
* Demertzi et al., *PNAS* (2019)
