#!/usr/bin/env python3
"""
parse_contamination.py

Combine multiple VerifyBamID2 *.selfSM files into a single contamination summary.

Used in Snakemake:
    snakemake.input   -> list of .selfSM paths
    snakemake.output[0] -> Contamination.txt

Output columns:
    Sample, SNPS, READS, AVG_DP, FREEMIX
"""

import pandas as pd

# Snakemake provides:
input_files = list(snakemake.input)       # noqa: F821
out_file = str(snakemake.output[0])       # noqa: F821

rows = []

for path in input_files:
    # Load .selfSM file (whitespace-separated)
    df = pd.read_csv(path, sep=r"\s+")

    if df.empty:
        continue

    # Normalize SEQ_ID field if header contains '#SEQ_ID'
    if "#SEQ_ID" in df.columns:
        df = df.rename(columns={"#SEQ_ID": "SEQ_ID"})

    row = df.iloc[0]

    rows.append({
        "Sample": row.get("SEQ_ID", "NA"),
        "SNPS": row.get("#SNPS", row.get("SNPS", pd.NA)),
        "READS": row.get("#READS", row.get("READS", pd.NA)),
        "AVG_DP": row.get("AVG_DP", pd.NA),
        "FREEMIX": row.get("FREEMIX", pd.NA),
    })

# Assemble into a clean table
out_df = pd.DataFrame(rows)

# Guarantee header even if empty
if out_df.empty:
    out_df = pd.DataFrame(columns=["Sample", "SNPS", "READS", "AVG_DP", "FREEMIX"])

# Save table
out_df.to_csv(out_file, sep="\t", index=False)
