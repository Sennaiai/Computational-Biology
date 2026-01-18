## Genomic Ancestry Inference Pipeline

This repository implements a fully automated, end-to-end genomic ancestry inference pipeline built using **Snakemake**. The workflow is designed for reproducible and scalable population assignment in human genomic datasets, particularly for low-coverage sequencing data.

The pipeline integrates contamination assessment, principal component analysis (PCA) extraction, and ancestry classification using **1000 Genomes reference populations**, providing a unified framework for ancestry inference suitable for large-scale genomic studies.

---

## Workflow Overview

The pipeline executes the following major stages:

1. **Reference Genome Preparation**
   Indexing the GRCh38 reference genome to support downstream alignment-based analyses.

2. **Contamination Assessment and PCA Estimation**
   Running **VerifyBamID2** on each input CRAM file to estimate sample contamination and extract ancestry-informative principal components.

3. **PCA Standardization**
   Parsing VerifyBamID2 outputs and harmonizing PCA coordinates across samples.

4. **Reference Population Integration**
   Extracting and harmonizing PCA coordinates from the 1000 Genomes reference panel.

5. **Ancestry Assignment**
   Training a **k-nearest neighbors (KNN)** classifier on reference population PCs and assigning ancestry labels to study samples.

6. **Visualization and Reporting**
   Generating PCA plots and summary tables for downstream interpretation and quality control.

All steps are orchestrated through **Snakemake**, ensuring reproducibility, modular execution, and transparent dependency tracking.

---

## Inputs

The workflow requires the following inputs:

* Low-coverage CRAM files (e.g., HGDP or cohort samples)
* GRCh38 reference genome with associated index files
* 1000 Genomes population metadata table
* 1000 Genomes PCA matrix (`.V`) compatible with VerifyBamID2

Sample configuration is managed through:

```
config/samples.tsv
```

---

## Outputs

The pipeline produces a comprehensive set of analysis-ready outputs, including:

* VerifyBamID2 contamination metrics (`*.selfSM` and summary tables)
* Per-sample PCA coordinates extracted from ancestry outputs
* Harmonized reference PCA table with population labels
* Predicted ancestry labels for study samples
* PCA visualization plots (PC1–PC4 combinations)
* Final summary tables integrating contamination and ancestry results

These outputs support both quality control and downstream population genetics analysis.

---

## Reproducibility

All software dependencies are defined in `environment.yml`. The pipeline can be executed in a fully reproducible environment using:

```bash
conda env create -f environment.yml
conda activate genomics-env
snakemake --cores N
```

All core data-processing logic is encapsulated in the `scripts/` directory, enabling transparent inspection and modification of each analytical step.

---

## Adaptability

The workflow is modular by design and can be readily extended to:

* Alternative reference panels (e.g., gnomAD, TOPMed)
* Additional principal components
* Different ancestry classifiers
* Expanded quality control metrics

This makes the pipeline suitable for both research and clinical genomics applications.
