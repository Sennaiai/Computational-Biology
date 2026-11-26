"""Run the simulation-first Cas9 prototype and its published-assay reference."""

import argparse
from dataclasses import asdict, replace
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path
import platform
from time import perf_counter

import numpy as np

from .calibration import run_calibration
from .cell import (
    CellParameters,
    DEFAULT_LOCI,
    DNA_STATES,
    MOLECULES_PER_NM_FL,
    POOLS,
    simulate_founder,
)


def scenarios(calibration):
    p = CellParameters()
    fitted = calibration["full"]["parameters"]
    return {
        "early_expression": p,
        "late_expression": replace(p, expression_on_h=8, expression_off_h=20),
        "higher_Cas9_production": replace(
            p, cas9_production_nM_h=2 * p.cas9_production_nM_h
        ),
        "guide2_production_imbalance": replace(p, guide_production_nM_h=(1, 0.1)),
        "example_fitted_repair": replace(
            p, precise_repair_per_h=fitted[1], mutation_repair_per_h=fitted[2]
        ),
        "PhyB2_repair_reference": replace(
            p, precise_repair_per_h=0.0583, mutation_repair_per_h=0.0864
        ),
        "no_Cas9_control": replace(p, cas9_production_nM_h=0),
    }


def observables(founder):
    a = np.asarray(founder["observations"])
    cells = a[:, 0]
    return {
        "modeled_cells": cells,
        "double_biallelic_disruption_proxy": a[:, 1] / cells,
        "selected_off_target_mutant_cell_fraction": a[:, 2] / cells,
        "unresolved_break_cell_fraction": a[:, 3] / cells,
        "free_RNP1_nM": a[:, 4] / a[:, 6] / MOLECULES_PER_NM_FL,
        "free_RNP2_nM": a[:, 5] / a[:, 6] / MOLECULES_PER_NM_FL,
        "target1_mutant_allele_fraction": a[:, 7] / (2 * cells),
        "target2_mutant_allele_fraction": a[:, 8] / (2 * cells),
    }


def bounded_mean(values, delta=0.05):
    if (
        not isinstance(delta, (int, float))
        or isinstance(delta, bool)
        or not 0 < delta < 1
    ):
        raise ValueError("Invalid sampling-error probability")
    values = np.asarray(values, dtype=float)
    if (
        values.ndim != 1
        or len(values) == 0
        or not np.isfinite(values).all()
        or np.any((values < 0) | (values > 1))
    ):
        raise ValueError("Expected bounded independent-founder outcomes")
    mean = float(values.mean())
    half = math.sqrt(math.log(2 / delta) / (2 * len(values)))
    return {
        "mean": mean,
        "MC_lower": max(0, mean - half),
        "MC_upper": min(1, mean + half),
    }


def summarize(founders, seeds):
    values = [observables(r) for r in founders]
    final = {}
    for name in (
        "double_biallelic_disruption_proxy",
        "selected_off_target_mutant_cell_fraction",
        "unresolved_break_cell_fraction",
        "target1_mutant_allele_fraction",
        "target2_mutant_allele_fraction",
    ):
        final[name] = bounded_mean([r[name][-1] for r in values])
    lower = np.array([r["mosaicism_proxy"] for r in founders], dtype=float)
    upper = np.array(
        [
            r["mosaicism_proxy"]
            or (
                r["terminal_cells_count"] >= 2
                and r["resolved_terminal_cells_count"] < r["terminal_cells_count"]
            )
            for r in founders
        ],
        dtype=float,
    )
    final["mosaicism_proxy_lower_event"] = bounded_mean(lower)
    final["mosaicism_proxy_upper_event"] = bounded_mean(upper)
    return {
        "founder_count": len(founders),
        "seeds": seeds,
        "mean_timecourses": {
            name: np.mean([v[name] for v in values], axis=0).tolist()
            for name in values[0]
        },
        "final": final,
        "mosaicism_evaluable_founders": sum(r["mosaicism_evaluable"] for r in founders),
        "terminal_cell_counts": [r["terminal_cells_count"] for r in founders],
        "founder_outcomes": [
            {
                "seed": seed,
                "terminal_cells": r["terminal_cells_count"],
                "resolved_terminal_cells": r["resolved_terminal_cells_count"],
                "mosaicism_proxy": r["mosaicism_proxy"],
                "molecular_ledger_error": r["molecular_ledger_error"],
                "cut_events": r["cut_events"],
                "precise_repair_events": r["precise_repair_events"],
                "mutation_repair_events": r["mutation_repair_events"],
                "terminal_observables": {name: float(v[name][-1]) for name in v},
            }
            for r, v, seed in zip(founders, values, seeds)
        ],
        "first_founder_trace": founders[0],
        "total_accepted_reactions": sum(r["accepted_reactions"] for r in founders),
        "total_thinning_proposals": sum(r["thinning_proposals"] for r in founders),
        "all_molecular_ledgers_exact": all(
            r["molecular_ledger_error"] == [0, 0, 0] for r in founders
        ),
    }


def run(output, *, founders=64, seed=2026):
    if (
        isinstance(founders, bool)
        or not isinstance(founders, int)
        or not 2 <= founders <= 256
    ):
        raise ValueError("Choose 2..256 independent founder lineages")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**32:
        raise ValueError("Invalid seed")
    start = perf_counter()
    calibration = run_calibration()
    cases = scenarios(calibration)
    streams = np.random.SeedSequence(seed).spawn(len(cases))
    times = np.linspace(0, 24, 13)
    records = {}
    for (name, parameters), stream in zip(cases.items(), streams):
        seeds = [int(s.generate_state(1)[0]) for s in stream.spawn(founders)]
        data = [simulate_founder(parameters, times_h=times, seed=s) for s in seeds]
        records[name] = {
            "parameters": {
                **asdict(parameters),
                "guide_production_nM_h": list(parameters.guide_production_nM_h),
            },
            **summarize(data, seeds),
        }
    if records["no_Cas9_control"]["final"]["target1_mutant_allele_fraction"][
        "mean"
    ] != 0 or any(not r["all_molecular_ledgers_exact"] for r in records.values()):
        raise ArithmeticError(
            "Cell study failed a negative control or conservation check"
        )
    package = Path(__file__).parent
    report = {
        "status": "simulation_prototype_executed",
        "biological_status": "Conditional scenario simulation; not validated Arabidopsis editing or mosaicism predictions",
        "case_study": "Two diploid target loci labeled ABCG18/ABCG19; no actual guide sequences or genome coordinates supplied",
        "units": {
            "time": "h",
            "volume": "fL",
            "production": "nM/h in the source compartment",
            "counts": "molecules",
        },
        "molecule_pools": list(POOLS),
        "DNA_states": list(DNA_STATES),
        "loci": [asdict(x) for x in DEFAULT_LOCI],
        "times_h": times.tolist(),
        "base_seed": seed,
        "founders_per_scenario": founders,
        "scenarios": records,
        "assay_calibration": calibration,
        "assumptions": [
            "Effective Cas9 protein production in cytoplasm; effective guide production in nucleus. Transcription/translation and guide processing are lumped.",
            "Expression scenarios have a 12-hour production window. Promoter names are not mapped to unmeasured concentrations.",
            "Cells grow to twice birth volume; growth pauses at scheduled division. Unrepaired DSBs delay division; pre-cleavage complexes dissociate before DNA replication. This checkpoint is an uncalibrated modeling choice.",
            "Division partitions free molecules binomially and copies diploid DNA states to both daughters. At most three divisions per lineage; not a whole-plant developmental model.",
            "Three synthetic resistant mutation categories: two disruption proxies, one other category. No sequence-resolved indel or functional-knockout prediction.",
            "One illustrative off-target locus; no genome-wide search, chromatin, inter-locus deletions, viability or tissue-specific developmental model.",
            "Repair parameters transferred from tomato protoplast references/example calibration; expression, target rates, cell-cycle settings and repair-category probabilities are assumptions.",
        ],
        "uncertainty": {
            "MC_bounds": "Per-quantity 95% Hoeffding bounds across independent founders; descendants within a founder are correlated and are not independent replicates. Bounds are not simultaneous across quantities or scenarios.",
            "parameter_uncertainty": "Expression/repair scenarios are sensitivity comparisons, not samples from a fitted biological posterior.",
            "mosaicism": "Lower event: at least two resolved terminal cells have different synthetic genotypes. Upper event also includes unresolved multi-cell lineages; incomplete repair leaves mosaicism undetermined. Proxy uses synthetic allele identities.",
        },
        "numerical_seconds": perf_counter() - start,
        "provenance": {
            "python": platform.python_version(),
            "packages": {
                name: metadata.version(name)
                for name in ("numpy", "scipy", "matplotlib")
            },
            "source_sha256": {
                str(p.relative_to(package)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (
                    *sorted(package.glob("*.py")),
                    *sorted((package / "data").glob("*.*")),
                )
            },
        },
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plot(report, output / "cellular_simulation.svg")
    (output / "cellular_results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(
        f"Executed {len(cases)} scenarios with {founders} independent founders each; all molecule ledgers conserved."
    )
    return report


def plot(report, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    t = report["times_h"]
    records = report["scenarios"]
    for name, color in (
        ("early_expression", "tab:blue"),
        ("guide2_production_imbalance", "tab:orange"),
    ):
        row = records[name]["mean_timecourses"]
        for guide, style in ((1, "-"), (2, "--")):
            axes[0, 0].plot(
                t,
                row[f"free_RNP{guide}_nM"],
                style,
                color=color,
                label=f"{name.replace('_', ' ')}; guide {guide}",
            )
    axes[0, 0].set(
        xlabel="Time (h)",
        ylabel="Free nuclear RNP (nM)",
        title="Expression and guide competition",
    )
    for name in list(records)[:4]:
        (line,) = axes[0, 1].plot(
            t,
            records[name]["mean_timecourses"]["target2_mutant_allele_fraction"],
            marker=".",
            label=name.replace("_", " "),
        )
        mean = np.asarray(
            records[name]["mean_timecourses"]["target2_mutant_allele_fraction"]
        )
        half = math.sqrt(math.log(40) / (2 * records[name]["founder_count"]))
        axes[0, 1].fill_between(
            t,
            np.maximum(0, mean - half),
            np.minimum(1, mean + half),
            color=line.get_color(),
            alpha=0.07,
        )
    axes[0, 1].set(
        xlabel="Time (h)",
        ylabel="Mutant-allele fraction at target 2",
        ylim=(0, 1),
        title="Conditional editing; pointwise 95% MC bounds",
    )
    names = list(records)[:-1]
    y = np.arange(len(names))
    lower = [records[n]["final"]["mosaicism_proxy_lower_event"]["mean"] for n in names]
    upper = [records[n]["final"]["mosaicism_proxy_upper_event"]["mean"] for n in names]
    axes[1, 0].hlines(y, lower, upper, color="tab:purple", linewidth=3)
    axes[1, 0].scatter(
        lower, y, label="Resolved heterogeneity", marker="o", color="tab:purple"
    )
    axes[1, 0].scatter(
        upper, y, label="Including unresolved lineages", marker="|", s=70, color="black"
    )
    axes[1, 0].set(
        yticks=y,
        yticklabels=[n.replace("_", " ") for n in names],
        xlim=(-0.02, 1.02),
        xlabel="Fraction of simulated founders",
        title="Mosaicism proxy: incomplete repair limits inference",
    )
    calibration = report["assay_calibration"]
    times = np.asarray(calibration["times_h"])
    counts = np.asarray(calibration["counts_in_source_order"])
    observed = counts / counts.sum(axis=1, keepdims=True)
    predicted = np.asarray(calibration["predicted_observed_fractions"])
    for column, label, color in (
        (1, "Indels", "tab:green"),
        (2, "Direct DSB", "tab:red"),
    ):
        axes[1, 1].scatter(
            times,
            observed[:, column],
            label=f"{label}: example observations",
            color=color,
            s=18,
        )
        axes[1, 1].plot(
            times[::2],
            predicted[::2, column],
            color=color,
            label=f"{label}: reduced model",
        )
    axes[1, 1].axvspan(60, 74, alpha=0.08, color="black")
    axes[1, 1].set(
        xlabel="Time (h); 72 h held out",
        ylabel="Observed molecule fraction",
        title="Published Psy1 example: assay reference",
    )
    for ax in axes.flat:
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)
    fig.suptitle(
        "Conditional cell-lineage prototype — cellular parameters are not Arabidopsis-calibrated",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
    path = Path(output)
    if path.suffix == ".svg":
        path.write_text(
            "\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--founders", type=int, default=64)
    args = parser.parse_args()
    run(args.output, founders=args.founders)


if __name__ == "__main__":
    main()
