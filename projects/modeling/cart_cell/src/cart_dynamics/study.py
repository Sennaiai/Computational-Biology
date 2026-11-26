"""Reproducible CAR-T component checks and conditional comparative experiment."""

import argparse
from dataclasses import asdict, replace
import hashlib
from importlib.metadata import version
from importlib.resources import files
import json
from pathlib import Path
import platform
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc, spearmanr

from .adapter import AdapterBinding, free_adapter_profile
from . import population, reference
from .signaling import Receptor, phenotype_signals
from .validation import integer


POLICIES = ("cd19", "cd22", "or", "and")
INITIAL = [0.2, 0.0, 0.7, 0.15, 0.14, 0.01]
VARIED = ("growth", "killing", "recall", "suppression", "loss19", "loss22")


def source_digest():
    root = Path(__file__).parent
    digest = hashlib.sha256()
    paths = sorted(root.glob("*.py")) + sorted((root / "data").glob("*.json"))
    for path in paths:
        digest.update(str(path.relative_to(root)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


def run_study(output, *, samples=64, seed=20261126):
    integer(samples, "samples", 4, 1024)
    integer(seed, "seed", 0, 2**32 - 1)
    start = time.perf_counter()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    receptor = Receptor()
    times = np.linspace(0, 30, 301)
    p = population.Parameters()
    trajectories = {
        policy: population.simulate(
            INITIAL, times, phenotype_signals(policy, receptor), p
        )
        for policy in POLICIES
    }
    no_car = population.simulate(
        [0, 0, *INITIAL[2:]], times, phenotype_signals("none"), p
    )
    # Structural control: remove the independent-arms bonus on double-positive targets.
    no_joint_bonus = np.array([receptor.equilibrium_active(1)] * 3 + [0])
    structural_control = population.simulate(INITIAL, times, no_joint_bonus, p)
    independent = population.simulate(
        INITIAL,
        times,
        phenotype_signals("or", receptor),
        p,
        method="Radau",
        rtol=1e-10,
        atol=1e-12,
    )
    solver_error = float(np.max(np.abs(trajectories["or"].values - independent.values)))
    if solver_error > 2e-6:
        raise RuntimeError("independent population solvers disagree")

    # Paired parameter draws, identical initial mixtures and horizons for every policy.
    draws = qmc.LatinHypercube(d=len(VARIED), seed=seed).random(samples)
    factors = np.exp(np.log(0.5) + draws * np.log(4))
    rows = []
    for index, factor in enumerate(factors):
        parameters = replace(
            p,
            **{
                name: getattr(p, name) * f
                for name, f in zip(VARIED, factor, strict=True)
            },
        )
        results = {}
        for policy in ("cd19", "or", "and"):
            trajectory = population.simulate(
                INITIAL, times, phenotype_signals(policy, receptor), parameters
            )
            results[policy] = trajectory.summary()
        rows.append(
            {"sample": index, "parameters": asdict(parameters), "results": results}
        )
    quantiles = {}
    for policy in ("cd19", "or", "and"):
        values = [r["results"][policy]["final_tumor_million_cells"] for r in rows]
        quantiles[policy] = dict(
            zip(
                ("q10", "median", "q90"),
                np.quantile(values, [0.1, 0.5, 0.9]).tolist(),
                strict=True,
            )
        )
    correlations = {}
    final_or = [row["results"]["or"]["final_tumor_million_cells"] for row in rows]
    for name in VARIED:
        coefficient = spearmanr(
            [row["parameters"][name] for row in rows], final_or
        ).statistic
        if not np.isfinite(coefficient):
            raise RuntimeError("undefined sensitivity correlation")
        correlations[name] = float(coefficient)

    # Reference roots are compared with rounded values independently printed in the article.
    sources = json.loads(
        files("cart_dynamics").joinpath("data/sources.json").read_text()
    )
    expected = np.asarray(sources["barros2021"]["reported_equilibria_cells"])
    equilibria = reference.coexistence_equilibria()
    equilibrium_error = float(np.max(np.abs(equilibria / expected - 1)))
    # Table-derived lower tumor root is 1091.197; printed text says 1089.
    # Preserve this source discrepancy; do not retune the table to match it.
    if equilibrium_error > 0.003:
        raise RuntimeError("published equilibrium comparison failed")
    counts = receptor.sample_counts(0.6, 2, 4000, seed)
    probability = receptor.probabilities(0.6, 2)
    mc_error = float(np.max(np.abs(counts / counts.sum() - probability)))
    mc_bound = float(np.sqrt(np.log(2 * len(counts) / 1e-6) / (2 * counts.sum())))
    if mc_error > mc_bound:
        raise RuntimeError("stochastic receptor reference exceeded its declared bound")

    adapter = AdapterBinding.published_affinities()
    report = {
        "schema_version": 1,
        "evidence_boundary": "Numerically verified, literature-referenced mechanistic prototype. No B-ALL assay calibration, clinical response, CRS prediction or quantum-advantage result.",
        "seed": seed,
        "environment": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "scipy", "matplotlib")},
        },
        "source_sha256": source_digest(),
        "configuration": {
            "horizon_days": 30,
            "output_interval_days": 0.1,
            "initial_million_cells": dict(
                zip(population.STATES[:6], INITIAL, strict=True)
            ),
            "parameters": asdict(p),
            "receptor": asdict(receptor),
            "solver": "DOP853",
            "rtol": 1e-8,
            "atol": 1e-10,
            "max_step_days": 0.5,
        },
        "nominal": {
            key: trajectory.summary() for key, trajectory in trajectories.items()
        },
        "no_car_control": no_car.summary(),
        "structural_control": {
            "description": "Dual recognition with no double-positive joint-probability bonus: max(s19,s22), instead of independent-opportunity OR. This is a dependence sensitivity check, not a receptor-expression calibration.",
            "signals": no_joint_bonus.tolist(),
            "result": structural_control.summary(),
        },
        "numerical_checks": {
            "dop853_radau_max_absolute_difference": solver_error,
            "published_equilibria_cells_E_M_T": equilibria.tolist(),
            "published_equilibrium_max_relative_difference": equilibrium_error,
            "printed_equilibrium_residual_cells_per_day": [
                reference.rhs(0, state, reference.published_parameters()).tolist()
                for state in expected
            ],
            "receptor_mc_max_absolute_difference": mc_error,
            "receptor_mc_simultaneous_hoeffding_bound_delta_1e_6": mc_bound,
            "receptor_mc_counts": counts.tolist(),
            "receptor_matrix_exponential_probabilities": probability.tolist(),
        },
        "recognition": {
            policy: phenotype_signals(policy, receptor).tolist() for policy in POLICIES
        },
        "uncertainty": {
            "design": "Paired Latin-hypercube draws; independent log-uniform factors from 0.5 to 2 of six synthetic nominal rates. Exploratory bounds, not a fitted posterior or population distribution.",
            "samples": samples,
            "varied": list(VARIED),
            "final_tumor_million_cells_quantiles": quantiles,
            "or_final_tumor_spearman_descriptive_only": correlations,
            "paired_or_below_cd19_count": sum(
                r["results"]["or"]["final_tumor_million_cells"]
                < r["results"]["cd19"]["final_tumor_million_cells"]
                for r in rows
            ),
            "runs": rows,
        },
        "adapter_component": {
            "parameters_nM": asdict(adapter),
            "peak_free_adapter_nM": adapter.peak_free_adapter,
            "maximum_bridge_fraction_of_limiting_receptor": adapter.bridge_fraction(
                adapter.peak_free_adapter
            ),
            "exposure_half_life_minutes_synthetic": 60,
            "boundary": "Separate free-ligand equilibrium demonstration for FITC-folate; no coupling to CD19/CD22 populations, dose conversion or residence-time prediction.",
        },
    }
    _figures(output, receptor, adapter, times, trajectories, no_car, rows)
    report["end_to_end_seconds_including_figures"] = time.perf_counter() - start
    (output / "cart_results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


def _figures(output, receptor, adapter, times, trajectories, no_car, rows):
    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "cart-study-v1"})
    colors = {"cd19": "#285f95", "cd22": "#6886a7", "or": "#147a63", "and": "#b96836"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.7), layout="constrained")
    for name, trajectory in trajectories.items():
        axes[0].plot(times, trajectory.tumor, label=name.upper(), color=colors[name])
        axes[1].plot(times, trajectory.values[:, 6], color=colors[name])
    axes[0].plot(times, no_car.tumor, "--", color="#555555", label="No CAR-T")
    axes[0].set(
        yscale="log",
        xlabel="Time (days)",
        ylabel="Tumor (million cells)",
        title="Conditional antigen-escape dynamics",
    )
    axes[0].legend(fontsize=8)
    axes[1].set(
        xlabel="Time (days)",
        ylabel="Cytokine proxy (arbitrary units)",
        title="Activation-dependent model output",
    )
    boxdata = [
        [r["results"][name]["final_tumor_million_cells"] for r in rows]
        for name in ("cd19", "or", "and")
    ]
    axes[2].boxplot(boxdata, tick_labels=["CD19", "OR", "AND"], showfliers=True)
    axes[2].set(
        yscale="log",
        ylabel="Day-30 tumor (million cells)",
        title=f"{len(rows)} paired parameter scenarios",
    )
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "cart_population.svg", metadata={"Date": None})
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.7), layout="constrained")
    minutes = np.linspace(0, 30, 151)
    for off in (0.2, 1.0):
        model = replace(receptor, dissociation=off)
        axes[0].plot(
            minutes,
            [model.probabilities(1, t)[-1] for t in minutes],
            label=f"koff={off}/min",
        )
    axes[0].set(
        xlabel="Time (minutes)",
        ylabel="Active receptor probability",
        title="Finite-state proofreading component",
    )
    axes[0].legend(fontsize=8)
    a = np.geomspace(1e-5, 1e4, 301)
    axes[1].semilogx(
        a, [adapter.bridge_fraction(value) for value in a], color="#147a63"
    )
    axes[1].axvline(adapter.peak_free_adapter, linestyle="--", color="#777777")
    axes[1].set(
        xlabel="Free FITC-folate (nM)",
        ylabel="Bridge / limiting receptor total",
        title="Separate equilibrium adapter component",
    )
    minutes = np.linspace(0, 480, 481)
    for label, pulses in (
        ("One specified pulse", [(0, 0.1)]),
        ("Three specified pulses", [(0, 0.1), (120, 0.1), (240, 0.1)]),
    ):
        free = free_adapter_profile(minutes, pulses, 60)
        axes[2].plot(
            minutes, [adapter.bridge_fraction(value) for value in free], label=label
        )
    axes[2].set(
        xlabel="Time (minutes)",
        ylabel="Equilibrium bridge fraction",
        title="Synthetic free-adapter exposure",
    )
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "cart_components.svg", metadata={"Date": None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20261126)
    parser.add_argument("--pkpd-samples", type=int, default=16)
    args = parser.parse_args()
    report = run_study(args.output_dir, samples=args.samples, seed=args.seed)
    for name, result in report["nominal"].items():
        print(
            f"{name:5s}: day-30 tumor {result['final_tumor_million_cells']:.6f} million cells; proxy AUC {result['cytokine_proxy_auc_arbitrary_unit_days']:.4f}"
        )
    from .pkpd_study import run_study as run_pkpd

    pkpd_report = run_pkpd(args.output_dir, samples=args.pkpd_samples, seed=args.seed)
    print("Separate adapter architecture (day 14; proxies in arbitrary units):")
    for name, result in pkpd_report["nominal"].items():
        print(
            f"{name:14s}: tumor {result['final_tumor_million_cells']:.4f} million; peak proxy {result['peak_cytokine_proxy_arbitrary_units']:.4f}"
        )
    print(f"Saved evidence to {args.output_dir}")


if __name__ == "__main__":
    main()
