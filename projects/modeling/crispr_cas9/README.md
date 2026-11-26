# CRISPR-Cas9 kinetics

A classical numerical baseline for molecular editing. It models target binding, reversible R-loop formation, cleavage and a separate repair layer. **Parameters are synthetic; biological calibration has not been performed.**

Seven target states form a continuous-time Markov model. A fixed free-RNP reservoir gives pseudo-first-order binding. Repair outcomes are terminal: exact repair is not re-cut, indels are not automatically functional knockouts, and HDR is a labeled branch rather than a sequence-resolved repair model. Chromatin, delivery and enzyme depletion are outside this first baseline.

The matrix exponential and an independent ODE integration solve the same generator. Exact Gillespie trajectories are compared under a predeclared simultaneous sampling bound. Closed-form binding, irreversible-cleavage and branching-repair cases provide independent tests. The model reports first-cleavage probability separately from current unrepaired breaks and altered-repair outcomes.

From the repository root, using Python 3.11:

```bash
python -m pip install -e 'projects/modeling/crispr_cas9[notebook]'
python -m unittest discover -s tests/modeling/crispr_cas9 -v
crispr-kinetics --output artifacts/crispr_cas9
python projects/modeling/crispr_cas9/scripts/execute_notebook.py --output-dir artifacts/crispr_cas9
```

The [executed notebook](CRISPR_Cas9_Kinetics.ipynb) shows the generator and distinct cleavage/repair readouts. The exported JSON records all state probabilities, stochastic counts, seeds, source hashes, dependency versions and numerical runtime. CI executes this study only when its files change; manual execution remains available.

The [reference run](results/kinetics.json) gives a maximum ODE/matrix-exponential difference of **2.65 × 10⁻¹¹**. Three ensembles of 4,000 trajectories have a maximum state-probability gap of **0.0223**, below the predeclared simultaneous bound **0.0520**. These measure numerical consistency, not biological predictive accuracy.

Mechanistic references: [Farasat and Salis, 2016](https://doi.org/10.1371/journal.pcbi.1004724), [Gong et al., 2018](https://doi.org/10.1016/j.celrep.2017.12.041), and [Eslami-Mossallam et al., 2022](https://doi.org/10.1038/s41467-022-28994-2). This coarse-grained baseline is not a reproduction of their full parameterized models. Published-assay calibration, sequence-level learning and any quantum extension require separate validation.
