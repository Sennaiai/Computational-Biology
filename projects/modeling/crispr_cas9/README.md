# CRISPR-Cas9 cell and lineage simulation

A mechanistic prototype for exploring how expression timing, Cas9/gRNA availability, binding, cleavage, repair and cell divisions influence editing. The case study is the Arabidopsis ABCG18/ABCG19 experiment. **Cellular results are conditional simulations, not validated predictions for those genes or a named promoter.**

The [cellular notebook](Cellular_Editing_Simulation.ipynb) compares seven scenarios, including guide imbalance, alternative repair rates and a no-Cas9 control. Finite molecule counts, shared Cas9 competition, nuclear transport, precise repair/re-cutting and allele inheritance are explicit. A simplified break-dependent division checkpoint is a consequential assumption. Mosaicism is a synthetic-genotype proxy; unresolved breaks remain separate.

A reduced repair model is fitted to the authors' public [Psy1 example assay](https://github.com/fabrimafe/DSBtimecourse), with both 72-hour replicates held out. Held-out mean absolute fraction error is **0.0398**, versus **0.0450** without precise repair. These are tomato example-assay results. Expression, development and Arabidopsis transfer remain uncalibrated. [Model and data definitions](MODEL.md).

From the repository root with Python 3.11:

```bash
python -m pip install -e 'projects/modeling/crispr_cas9[notebook]'
python -m unittest discover -s tests/modeling/crispr_cas9 -v
python -m crispr_kinetics.cell_study --output artifacts/cellular --founders 64
python projects/modeling/crispr_cas9/scripts/execute_notebook.py --notebook Cellular_Editing_Simulation.ipynb --output-dir artifacts/cellular
```

[Reference results](results/cellular_results.json) record seeds, parameters, molecular ledgers, founder-level uncertainty and source hashes. The [single-target notebook](CRISPR_Cas9_Kinetics.ipynb) retains the matrix-exponential/ODE/Gillespie baseline. CI executes this project only for its own changes.

## Future directions

The author's 2022 pre-analysis cuticular-transport study sought ABCG18/ABCG19 double mutants by editing one locus in an existing mutant background. It describes Gibson-assembled tRNA–gRNA constructs and assembly ambiguities, without establishing editing or phenotype outcomes.

- Benchmark matched components against public Arabidopsis [promoter/genotype observations](https://doi.org/10.1186/s13059-015-0715-0) and [repair-outcome data](https://doi.org/10.1093/pnasnexus/pgae094).
- Represent existing mutant backgrounds and verified vector/guide architectures.
- If recovered, use T1 HRM screens with controls and sample metadata for screening-level comparisons. Whole-genome sequencing is unavailable; HRM alone cannot establish exact genotypes or repair kinetics.

Remaining predictions stay explicitly conditional. The prototype does not determine whether repair or mosaicism explained the original experiment.
