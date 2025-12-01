#!/usr/bin/env python3
"""
plot_pcs.py

Reads:
  - output/ancestry/reference_pcs.tsv (ID, POPULATION, PC1–PC4)
  - output/ancestry/study_pcs.tsv     (SAMPLE, PC1–PC4)

Produces:
  - output/plots/PC1_PC2.png
  - output/plots/PC2_PC3.png
  - output/plots/PC3_PC4.png
  - output/plots/PC1_PC2_PC3.png  (3D scatter plot)
"""

import pathlib
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


POP_COLORS = {
    "EUR": "tab:blue",
    "EAS": "tab:orange",
    "AFR": "tab:green",
    "SAS": "tab:red",
    "AMR": "tab:purple",
}


def make_2d_plot(ref, study, x, y, out_path):
    plt.figure(figsize=(6, 5))

    # Reference samples colored by population
    for pop, group in ref.groupby("POPULATION"):
        c = POP_COLORS.get(pop, "gray")
        plt.scatter(
            group[x],
            group[y],
            s=10,
            alpha=0.5,
            label=pop,
            color=c,
        )

    # Study samples as yellow stars
    plt.scatter(
        study[x],
        study[y],
        s=40,
        marker="*",
        edgecolor="k",
        facecolor="yellow",
        label="Study samples",
        zorder=5,
    )

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{x} vs {y}")
    plt.legend(frameon=False, fontsize=8)

    plt.tight_layout()
    pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.close()


def make_3d_plot(ref, study, out_path):
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Reference samples
    for pop, group in ref.groupby("POPULATION"):
        c = POP_COLORS.get(pop, "gray")
        ax.scatter(
            group["PC1"],
            group["PC2"],
            group["PC3"],
            s=10,
            alpha=0.5,
            label=pop,
            color=c,
        )

    # Study samples
    ax.scatter(
        study["PC1"],
        study["PC2"],
        study["PC3"],
        s=50,
        marker="*",
        edgecolor="k",
        facecolor="yellow",
        label="Study samples",
    )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title("PC1 vs PC2 vs PC3")
    ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0.0)

    plt.tight_layout()
    pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def main(snakemake=None):
    if snakemake is not None:
        pcs_file = snakemake.input.pcs
        ref_file = snakemake.input.ref
        (
            out_pc1_pc2,
            out_pc2_pc3,
            out_pc3_pc4,
            out_pc1_pc2_pc3,
        ) = snakemake.output
    else:
        # Manual test mode
        ref_file = "output/ancestry/reference_pcs.tsv"
        pcs_file = "output/ancestry/study_pcs.tsv"
        out_pc1_pc2 = "output/plots/PC1_PC2.png"
        out_pc2_pc3 = "output/plots/PC2_PC3.png"
        out_pc3_pc4 = "output/plots/PC3_PC4.png"
        out_pc1_pc2_pc3 = "output/plots/PC1_PC2_PC3.png"

    ref = pd.read_csv(ref_file, sep="\t")
    study = pd.read_csv(pcs_file, sep="\t")

    make_2d_plot(ref, study, "PC1", "PC2", out_pc1_pc2)
    make_2d_plot(ref, study, "PC2", "PC3", out_pc2_pc3)
    make_2d_plot(ref, study, "PC3", "PC4", out_pc3_pc4)
    make_3d_plot(ref, study, out_pc1_pc2_pc3)

    print("All PCA plots generated.")


if "snakemake" in globals():
    main(snakemake)
else:
    main()
