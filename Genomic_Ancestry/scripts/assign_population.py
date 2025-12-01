#!/usr/bin/env python3
"""
assign_population.py

Assign populations to study samples using K-nearest neighbors (KNN) in PC space.

Two modes of use:

1) Snakemake mode (recommended in the pipeline):
   - Expects:
       snakemake.input.ref   -> reference_pcs.tsv
       snakemake.input.pcs   -> study_pcs.tsv
       snakemake.output[0]   -> Populations.txt (or similar)

2) CLI mode:
   Usage:
       python scripts/assign_population.py REF_PCS STU_PCS OUT_TSV

   Where:
       REF_PCS = output/ancestry/reference_pcs.tsv
                 (columns: ID, POPULATION, PC1, PC2, PC3, PC4)
       STU_PCS = output/ancestry/study_pcs.tsv
                 (columns: SAMPLE, PC1, PC2, PC3, PC4)
       OUT_TSV = path to write the population assignments
"""

import sys
import pathlib
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier


PC_COLS = ["PC1", "PC2", "PC3", "PC4"]
REQUIRED_REF_COLS = {"ID", "POPULATION", *PC_COLS}
REQUIRED_STU_COLS = {"SAMPLE", *PC_COLS}


def assign_populations(ref_file: str, stu_file: str, out_file: str) -> None:
    """Core logic: load reference + study PCs, fit KNN, write assignments."""
    ref = pd.read_csv(ref_file, sep="\t")
    stu = pd.read_csv(stu_file, sep="\t")

    # Check required columns
    if not REQUIRED_REF_COLS.issubset(ref.columns):
        missing = REQUIRED_REF_COLS - set(ref.columns)
        raise ValueError(f"Reference file missing required columns: {missing}")

    if not REQUIRED_STU_COLS.issubset(stu.columns):
        missing = REQUIRED_STU_COLS - set(stu.columns)
        raise ValueError(f"Study file missing required columns: {missing}")

    # Ensure PCs are numeric
    for c in PC_COLS:
        ref[c] = pd.to_numeric(ref[c], errors="coerce")
        stu[c] = pd.to_numeric(stu[c], errors="coerce")

    before_ref = len(ref)
    ref = ref.dropna(subset=PC_COLS + ["POPULATION"])
    after_ref = len(ref)
    print(f"[INFO] Reference rows: {before_ref} → {after_ref} after dropping NaNs")

    if after_ref == 0:
        raise ValueError("No valid reference rows left after dropping NaNs.")

    before_stu = len(stu)
    stu = stu.dropna(subset=PC_COLS)
    after_stu = len(stu)
    print(f"[INFO] Study rows: {before_stu} → {after_stu} after dropping NaNs")

    if after_stu == 0:
        raise ValueError("No valid study rows left after dropping NaNs.")

    # Train KNN on reference PCs
    X_train = ref[PC_COLS].values
    y_train = ref["POPULATION"].values

    X_test = stu[PC_COLS].values

    knn = KNeighborsClassifier(n_neighbors=5)
    knn.fit(X_train, y_train)

    preds = knn.predict(X_test)

    out_df = stu.copy()
    # Use POPULATION as final column name (keep things clean)
    out_df["POPULATION"] = preds

    pathlib.Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_file, sep="\t", index=False)

    print(f"[DONE] Saved population assignments → {out_file}")


def main_cli():
    """Entry point for CLI usage."""
    if len(sys.argv) != 4:
        print(
            "Usage: assign_population.py REF_PCS STU_PCS OUT_TSV",
            file=sys.stderr,
        )
        sys.exit(1)

    ref_file = sys.argv[1]
    stu_file = sys.argv[2]
    out_file = sys.argv[3]

    assign_populations(ref_file, stu_file, out_file)


# Dispatch between Snakemake mode and CLI mode
if __name__ == "__main__":
    if "snakemake" in globals():
        # Snakemake provides the input/output paths
        ref_file = snakemake.input.ref       # noqa: F821
        stu_file = snakemake.input.pcs       # noqa: F821
        out_file = snakemake.output[0]       # noqa: F821

        assign_populations(ref_file, stu_file, out_file)
    else:
        main_cli()
