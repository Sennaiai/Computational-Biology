#!/usr/bin/env python3
"""
extract_pcs.py (CLI version)

Usage:
    python scripts/extract_pcs.py \
        V_FILE POP_FILE OUT_REF OUT_STUDY ANCESTRY1 [ANCESTRY2 ...]

Where:
- V_FILE     = data/1000G/1000g.phase3.100k.b38.vcf.gz.dat.V
- POP_FILE   = data/1000G/1000G_reference_populations.txt
- OUT_REF    = output/ancestry/reference_pcs.tsv
- OUT_STUDY  = output/ancestry/study_pcs.tsv
- ANCESTRY*  = one or more *.Ancestry files from VerifyBamID2
"""

import sys
import pathlib
import pandas as pd


# 1. Load 1000G reference PCs + populations

def load_reference_pcs(v_file: str, pop_file: str) -> pd.DataFrame:
    """
    Load 1000G reference PCs (from .V) and population labels
    (from 1000G_reference_populations.txt), and join them horizontally.

    Handles two common .V layouts:

    A) PC1..PC4 only:
           PC1 PC2 PC3 PC4
    B) sampleID + PCs:
           ID  PC1 PC2 PC3 PC4

    Returns: ID, POPULATION, PC1, PC2, PC3, PC4
    """

    # Raw V matrix (any whitespace)
    v_raw = pd.read_csv(
        v_file,
        delim_whitespace=True,
        header=None,
        comment="#"
    )

    # Population labels
    pops = pd.read_csv(
        pop_file,
        delim_whitespace=True,
        header=None,
        names=["ID", "POPULATION"],
    )

    # Align row counts
    if len(v_raw) != len(pops):
        print(
            f"[WARN] V has {len(v_raw)} rows, population file has {len(pops)} rows. "
            "Truncating both to the minimum length.",
            file=sys.stderr,
        )
        n = min(len(v_raw), len(pops))
        v_raw = v_raw.iloc[:n].reset_index(drop=True)
        pops = pops.iloc[:n].reset_index(drop=True)

    # Decide where PCs start:
    # If first column of v_raw matches ID (e.g. HG00096), treat col 0 as ID,
    # and use the next 4 columns as PC1..PC4.
    pc_start = 0

    if v_raw.shape[1] >= 5:
        try:
            if str(v_raw.iloc[0, 0]) == str(pops.iloc[0, 0]):
                pc_start = 1
        except Exception:
            pc_start = 0

    # Take 4 PC columns from v_raw starting at pc_start
    if v_raw.shape[1] < pc_start + 4:
        raise ValueError(
            f"Not enough columns in V file {v_file} to extract 4 PCs "
            f"starting at column {pc_start} (total columns: {v_raw.shape[1]})."
        )

    V = v_raw.iloc[:, pc_start:pc_start + 4].copy()
    V.columns = ["PC1", "PC2", "PC3", "PC4"]

    # Ensure PCs are numeric
    for c in ["PC1", "PC2", "PC3", "PC4"]:
        V[c] = pd.to_numeric(V[c], errors="coerce")

    ref = pd.concat([pops, V], axis=1)
    return ref



# 2. Helper: infer sample ID from row or filename

def infer_sample_id_from_row_or_path(row: pd.Series, path: pathlib.Path) -> str:
    """
    Try to infer sample ID from row columns; if that fails, fall back to filename.
    """
    for key in ["SampleID", "SEQ_ID", "#SEQ_ID", "Sample", "SAMPLE", "ID"]:
        if key in row.index and pd.notna(row[key]):
            return str(row[key])

    stem = path.stem
    return stem.split(".")[0]



# 3. Load study PCs from *.Ancestry files

def load_study_pcs(ancestry_paths) -> pd.DataFrame:
    """
    Load PCs for study samples from VerifyBamID2 *.Ancestry files.

    Handles two formats:

    1) Vertical format:
       columns: ['PC', 'ContaminatingSample', 'IntendedSample']
       - One row per PC
       - Take first 4 PCs (PC=1..4) from 'IntendedSample'

    2) Wide format:
       - Row(s) with PC1, PC2, PC3, PC4, etc.
    """
    rows = []

    for p in map(pathlib.Path, ancestry_paths):
        df = pd.read_csv(p, sep=r"\s+", comment="#", engine="python")

        if df.empty:
            print(f"[WARN] Ancestry file {p} is empty, skipping.", file=sys.stderr)
            continue

        cols = list(df.columns)

        # Case 1: vertical PC format
        if set(cols) == {"PC", "ContaminatingSample", "IntendedSample"}:
            df_sorted = df.sort_values("PC")
            top = df_sorted.head(4).reset_index(drop=True)

            sample_id = p.stem.split(".")[0]

            pcs = {"SAMPLE": sample_id}
            for i in range(4):
                if i < len(top):
                    pcs[f"PC{i+1}"] = float(top.loc[i, "IntendedSample"])
                else:
                    pcs[f"PC{i+1}"] = 0.0

            rows.append(pcs)
            continue

        # Case 2: wide format
        if "SampleType" in df.columns:
            df_intended = df.loc[df["SampleType"] == "IntendedSample"]
            if df_intended.empty:
                row = df.iloc[0]
            else:
                row = df_intended.iloc[0]
        else:
            row = df.iloc[0]

        pc_like = [
            c for c in df.columns
            if isinstance(c, str) and c.upper().startswith("PC")
        ]
        pc_like_sorted = sorted(pc_like)
        if len(pc_like_sorted) < 4:
            raise ValueError(
                f"Could not find 4 PC columns in Ancestry file {p}. "
                f"Available columns: {list(df.columns)}"
            )

        pc_cols = pc_like_sorted[:4]
        sample_id = infer_sample_id_from_row_or_path(row, p)

        pcs = {"SAMPLE": sample_id}
        for i, col in enumerate(pc_cols, start=1):
            pcs[f"PC{i}"] = float(row[col])

        rows.append(pcs)

    if not rows:
        raise ValueError("No valid PCs found in any Ancestry files.")

    return pd.DataFrame(rows)


# 4. Main CLI

def main():
    if len(sys.argv) < 6:
        print(
            "Usage: extract_pcs.py V_FILE POP_FILE OUT_REF OUT_STUDY "
            "ANCESTRY1 [ANCESTRY2 ...]",
            file=sys.stderr,
        )
        sys.exit(1)

    v_file = sys.argv[1]
    pop_file = sys.argv[2]
    out_ref = sys.argv[3]
    out_study = sys.argv[4]
    ancestry_files = sys.argv[5:]

    ref_df = load_reference_pcs(v_file, pop_file)
    study_df = load_study_pcs(ancestry_files)

    pathlib.Path(out_ref).parent.mkdir(parents=True, exist_ok=True)
    ref_df.to_csv(out_ref, sep="\t", index=False)

    pathlib.Path(out_study).parent.mkdir(parents=True, exist_ok=True)
    study_df.to_csv(out_study, sep="\t", index=False)

    print("Extracted PCs!")


if __name__ == "__main__":
    main()
