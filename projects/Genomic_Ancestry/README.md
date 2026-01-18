## Genomic Ancestry Pipeline

A Snakemake workflow for processing low-coverage CRAM files, running VerifyBamID2, extracting population PCs, and assigning ancestry labels using 1000 Genomes reference data.

This repository contains a fully automated end-to-end pipeline designed for reproducible ancestry inference in human genomic datasets. It integrates reference resources, contamination assessment, PCA extraction, population assignment, and visualization in a single, modular workflow.

## Workflow Overview

The pipeline performs the following tasks:

# Indexing reference genome (GRCh38).

# Running VerifyBamID2 on each input CRAM to generate contamination metrics and PC estimates.

# Parsing VerifyBamID2 outputs and standardizing the PCA representation.

# Extracting reference PCs from 1000 Genomes data and harmonizing sample IDs.

# Assigning ancestry to study samples using a K-nearest neighbors classifier trained on reference PCs.

# Producing summary tables and PCA plots for downstream interpretation.

All steps are managed through Snakemake, ensuring reproducibility, modularity, and transparent dependency handling.

## Inputs

The workflow expects:

# Low-coverage CRAM files (e.g., HGDP samples).

# GRCh38 full reference genome with index files.

# 1000 Genomes population reference table.

# 1000 Genomes PCA matrix (.V) compatible with VerifyBamID2.

Configuration of input samples is handled through config/samples.tsv.

## Outputs

The pipeline produces:

# VerifyBamID2 contamination estimates (*.selfSM and summary tables).

# Per-sample PCA coordinates extracted from .Ancestry outputs.

# Harmonized reference PCA table with population labels.

# Assigned population labels for study samples.

# PCA visualization plots (PC1–PC4 combinations).

# Final summary tables collating contamination and ancestry calls.


## Reproducibility

All software dependencies are specified in environment.yml.
The workflow is portable and can be executed with:

            conda env create -f environment.yml
            conda activate genomics-env
            snakemake --cores N


Scripts under scripts/ encapsulate all core processing steps and can be adapted independently if required for alternative population panels or extended PCA dimensions.

## Adaptability

The pipeline is modular and can be repurposed.