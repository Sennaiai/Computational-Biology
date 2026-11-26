"""Exposure-matched comparison motivated by the efficacy/toxicity review thesis."""

from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc

from . import pkpd, population
from .adapter import AdapterBinding
from .study import source_digest
from .validation import integer


def run_study(output, *, samples=16, seed=20261126):
    integer(samples, "samples", 4, 256)
    integer(seed, "seed", 0, 2**32 - 1)
    start = time.perf_counter()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    times = np.linspace(0, 14, 701)
    protocols = pkpd.protocols()
    trajectories = {
        protocol.name: pkpd.simulate(protocol, times) for protocol in protocols
    }
    nominal = {name: result.summary() for name, result in trajectories.items()}
    binding = AdapterBinding.published_affinities()
    base_p = replace(population.Parameters(), loss19=0, loss22=0)
    base_k = pkpd.Kinetics()
    tighter = pkpd.simulate(protocols[1], times, rtol=1e-9, atol=1e-13)
    endpoint_error = abs(
        tighter.summary()["final_tumor_million_cells"]
        - nominal["daily_bolus"]["final_tumor_million_cells"]
    )
    if endpoint_error > 1e-5:
        raise RuntimeError(
            "PKPD endpoints changed under tighter integration tolerances"
        )
    # Vary kinetics at fixed published equilibrium affinities. This deliberately
    # tests what Kd values alone cannot identify about residence and scheduling.
    design = qmc.LatinHypercube(d=6, seed=seed).random(samples)
    rows = []
    for i, draw in enumerate(design):
        k = replace(
            base_k,
            half_life_minutes=30 + draw[0] * 59.9,
            car_off_per_minute=0.002 * 100 ** draw[1],
            target_off_per_minute=0.02 * 100 ** draw[2],
            bridge_half_nM=base_k.bridge_half_nM * 0.5 * 4 ** draw[3],
        )
        p = replace(
            base_p,
            killing=base_p.killing * 0.5 * 4 ** draw[4],
            proliferation=base_p.proliferation * 0.5 * 4 ** draw[5],
        )
        outcomes = {
            protocol.name: pkpd.simulate(
                protocol, times, kinetics=k, parameters=p
            ).summary()
            for protocol in protocols[1:]
        }
        rows.append(
            {
                "sample": i,
                "kinetics": asdict(k),
                "population_parameters": asdict(p),
                "outcomes": outcomes,
                "peak_proxy_frontier": pkpd.nondominated(outcomes),
                "auc_proxy_frontier": pkpd.nondominated(
                    outcomes, "cytokine_proxy_auc_arbitrary_unit_days"
                ),
            }
        )
    report = {
        "schema_version": 1,
        "source_sha256": source_digest(),
        "seed": seed,
        "question": "At matched cumulative adapter input, how do exposure timing and unmeasured binding kinetics change tumor-control/activation-proxy tradeoffs?",
        "evidence_boundary": "Separate FR-positive adapter architecture; synthetic PK and cellular PD, with published binding Kds. No clinical dose, CRS risk, patient-response calibration or validated combined CD19/CD22 adapter product.",
        "configuration": {
            "horizon_days": 14,
            "input_window_days": 7,
            "output_interval_days": 0.02,
            "initial_million_cells": {
                "effector": 0.2,
                "memory": 0,
                "FR_positive_tumor": 0.99,
                "FR_negative_tumor": 0.01,
            },
            "population_parameters": asdict(base_p),
            "binding_nM": asdict(binding),
            "kinetics": asdict(base_k),
            "protocols": [asdict(protocol) for protocol in protocols],
            "integration": {
                "solver": "Radau",
                "rtol": 2e-7,
                "atol": 1e-11,
                "max_step_days": 0.1,
                "input_boundaries": "all boluses and infusion boundaries split exactly; sampled boundary states are post-input",
            },
        },
        "nominal": nominal,
        "nominal_peak_proxy_frontier_active_protocols": pkpd.nondominated(
            {name: value for name, value in nominal.items() if name != "no_adapter"}
        ),
        "nominal_auc_proxy_frontier_active_protocols": pkpd.nondominated(
            {name: value for name, value in nominal.items() if name != "no_adapter"},
            "cytokine_proxy_auc_arbitrary_unit_days",
        ),
        "numerical_checks": {
            "daily_bolus_final_tumor_difference_tighter_tolerances": endpoint_error,
            "matched_active_input_nM": 0.7,
            "nominal_free_adapter_auc_nM_days": 0.7 / base_k.elimination_per_day,
        },
        "sensitivity": {
            "design": "Paired Latin-hypercube scenarios, not a posterior or clinical population. Equilibrium Kds fixed; kon varies consistently with each sampled koff.",
            "samples": samples,
            "bounds": {
                "half_life_minutes": [30, 89.9],
                "car_off_per_minute_log_uniform": [0.002, 0.2],
                "target_off_per_minute_log_uniform": [0.02, 2],
                "bridge_half_factor_log_uniform": [0.5, 2],
                "killing_factor_log_uniform": [0.5, 2],
                "proliferation_factor_log_uniform": [0.5, 2],
            },
            "peak_proxy_frontier_counts": {
                p.name: sum(p.name in row["peak_proxy_frontier"] for row in rows)
                for p in protocols[1:]
            },
            "auc_proxy_frontier_counts": {
                p.name: sum(p.name in row["auc_proxy_frontier"] for row in rows)
                for p in protocols[1:]
            },
            "lowest_final_tumor_counts": {
                p.name: sum(
                    p.name
                    == min(
                        row["outcomes"],
                        key=lambda name: row["outcomes"][name][
                            "final_tumor_million_cells"
                        ],
                    )
                    for row in rows
                )
                for p in protocols[1:]
            },
            "runs": rows,
        },
        "interpretation": "Nondominance applies only to these schedules and declared endpoints. A high or low arbitrary-unit cytokine proxy is not a probability of toxicity. The minimal model does not establish the clinical decoupling of potency and toxicity discussed in the review; host cytokine feedback and exhaustion are not modeled.",
    }
    _figures(output, trajectories, nominal, rows)
    report["end_to_end_seconds_including_figures"] = time.perf_counter() - start
    (output / "pkpd_results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


def _figures(output, trajectories, nominal, rows):
    colors = {
        "daily_bolus": "#285f95",
        "split_bolus": "#147a63",
        "ramped_bolus": "#b96836",
        "continuous": "#70549b",
        "no_adapter": "#666666",
    }
    labels = {name: name.replace("_", " ").capitalize() for name in trajectories}
    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "cart-pkpd-v1"})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for name, trajectory in trajectories.items():
        t, y = trajectory.times, trajectory.values
        axes[0, 0].plot(
            t, y[:, 10], color=colors[name], label=labels[name], linewidth=1
        )
        axes[0, 1].plot(t, y[:, :2].sum(axis=1), color=colors[name])
        axes[1, 0].plot(t, y[:, 2:6].sum(axis=1), color=colors[name])
        axes[1, 1].plot(t, y[:, 6], color=colors[name])
    axes[0, 0].set(ylabel="Free adapter (nM)", title="Specified systemic exposure")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(
        ylabel="Effector + memory (million cells)",
        title="Conditional cellular expansion / persistence",
    )
    axes[1, 0].set(ylabel="Tumor (million cells)", title="Conditional tumor response")
    axes[1, 1].set(
        ylabel="Cytokine proxy (arbitrary units)", title="Inflammatory-activation proxy"
    )
    for ax in axes.flat:
        ax.set_xlabel("Time (days)")
        ax.axvline(7, color="#aaaaaa", linestyle="--", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "pkpd_trajectories.svg", metadata={"Date": None})
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for ax, endpoint, label in zip(
        axes,
        (
            "peak_cytokine_proxy_arbitrary_units",
            "cytokine_proxy_auc_arbitrary_unit_days",
        ),
        ("Peak proxy (arbitrary units)", "Proxy AUC (arbitrary unit·days)"),
        strict=True,
    ):
        for name in colors:
            if name == "no_adapter":
                continue
            cloud = np.array(
                [
                    [
                        r["outcomes"][name][endpoint],
                        r["outcomes"][name]["final_tumor_million_cells"],
                    ]
                    for r in rows
                ]
            )
            ax.scatter(cloud[:, 0], cloud[:, 1], color=colors[name], alpha=0.18, s=15)
            ax.scatter(
                nominal[name][endpoint],
                nominal[name]["final_tumor_million_cells"],
                color=colors[name],
                label=labels[name],
                s=65,
                edgecolors="white",
            )
        ax.set(
            xlabel=label,
            ylabel="Day-14 tumor (million cells)",
            title="Matched-input tradeoff: nominal + sensitivity",
        )
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=8)
    fig.savefig(output / "pkpd_tradeoffs.svg", metadata={"Date": None})
    plt.close(fig)
