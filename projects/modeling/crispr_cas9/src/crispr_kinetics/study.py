"""Execute a bounded synthetic kinetic benchmark and export its evidence."""

import argparse
from dataclasses import asdict, replace
import hashlib
from importlib import metadata, resources
import json
import math
import os
from pathlib import Path
import platform
from time import perf_counter

import numpy as np

from .kinetics import (
    STATES,
    Rates,
    diagnostics,
    first_cleavage_time,
    generator,
    readouts,
    simulate,
    solve,
)


def demonstration():
    data = json.loads(
        resources.files("crispr_kinetics").joinpath("data/demo_rates.json").read_text()
    )
    if data["provenance"] != "synthetic_demonstration":
        raise ValueError(
            "This experiment requires the disclosed synthetic configuration"
        )
    return Rates(**data["rates"]), data


def run(output):
    start = perf_counter()
    rates, source = demonstration()
    times = np.linspace(0, 600, 61)
    exact = solve(rates, times)
    ode = solve(rates, times, method="ode")
    gap = float(abs(exact - ode).max())
    if gap > 2e-9:
        raise ArithmeticError("ODE and matrix-exponential references disagree")
    paths, seeds, delta = 4000, (7, 19, 41), 1e-6
    epsilon = math.sqrt(math.log(2 * exact.size * len(seeds) / delta) / (2 * paths))
    stochastic = []
    for seed in seeds:
        counts = simulate(rates, times, paths=paths, seed=seed)
        error = float(abs(counts / paths - exact).max())
        if error > epsilon:
            raise ArithmeticError(
                "SSA/reference gap exceeds the predeclared simultaneous bound"
            )
        stochastic.append(
            {
                "seed": seed,
                "paths": paths,
                "max_absolute_probability_gap": error,
                "counts": counts.tolist(),
            }
        )
    sensitivity = []
    for concentration in (0.0, 1.0, 5.0, 10.0, 25.0, 50.0):
        config = replace(rates, free_rnp_nM=concentration)
        p = solve(config, [120.0])
        sensitivity.append(
            {
                "free_rnp_nM": concentration,
                "cleaved_ever_at_120_s": float(readouts(p)["cleaved_ever"][0]),
                "first_cleavage": first_cleavage_time(config),
            }
        )
    package = Path(__file__).parent
    report = {
        "status": "numerical_baseline_passed",
        "biological_calibration": "not_performed",
        "parameter_source": source,
        "scope": "One target molecule, fixed free-RNP reservoir, one cleavage followed by terminal repair outcomes; no re-cutting.",
        "interpretation": "Cleavage and repair probabilities under synthetic rates; not measured T-cell editing, sequence-specific off-target predictions, indel spectra or functional editing probabilities.",
        "states": list(STATES),
        "generator_convention": "dp/dt = Q @ p; columns sum to zero; rates in s^-1",
        "rates": asdict(rates),
        "generator_per_s": generator(rates).tolist(),
        "times_s": times.tolist(),
        "reference_probabilities": exact.tolist(),
        "ode_expm_max_absolute_gap": gap,
        "reference_invariants": diagnostics(exact),
        "ode_invariants": diagnostics(ode),
        "stochastic_validation": stochastic,
        "simultaneous_sampling_bound": {
            "inequality": "Hoeffding union bound over states, times and all three seed ensembles",
            "global_delta": delta,
            "epsilon": epsilon,
            "meaning": "Monte Carlo sampling validation at fixed model parameters; not biological or parameter uncertainty",
        },
        "mean_first_cleavage": first_cleavage_time(rates),
        "concentration_sensitivity": sensitivity,
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
                    package / "data/demo_rates.json",
                )
            },
            "thread_settings": {
                k: os.environ.get(k)
                for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            },
        },
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plot(report, output / "kinetics.svg")
    (output / "kinetics.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(summary(report), allow_nan=False), flush=True)
    return report


def summary(report):
    return {
        "status": report["status"],
        "parameters": report["parameter_source"]["provenance"],
        "biological_calibration": report["biological_calibration"],
        "ode_expm_max_absolute_gap": report["ode_expm_max_absolute_gap"],
        "max_SSA_probability_gap": max(
            v["max_absolute_probability_gap"] for v in report["stochastic_validation"]
        ),
        "simultaneous_sampling_bound": report["simultaneous_sampling_bound"]["epsilon"],
        "numerical_seconds": report["numerical_seconds"],
    }


def plot(report, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    times = np.asarray(report["times_s"])
    exact = np.asarray(report["reference_probabilities"])
    first = report["stochastic_validation"][0]
    sampled = np.asarray(first["counts"]) / first["paths"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].plot(times, readouts(exact)["cleaved_ever"], label="Matrix exponential")
    axes[0].scatter(
        times[::4],
        readouts(sampled)["cleaved_ever"][::4],
        s=18,
        label="Gillespie; seed 7",
    )
    axes[0].set(xlabel="Time (s)", ylabel="Probability", title="First cleavage")
    for state in (3, 4, 5, 6):
        axes[1].plot(times, exact[:, state], label=STATES[state].replace("_", " "))
    axes[1].set(xlabel="Time (s)", ylabel="Probability", title="Separate repair layer")
    rows = report["concentration_sensitivity"]
    axes[2].plot(
        [r["free_rnp_nM"] for r in rows],
        [r["cleaved_ever_at_120_s"] for r in rows],
        marker="o",
    )
    axes[2].set(
        xlabel="Fixed free-RNP concentration (nM)",
        ylabel="P(cleaved by 120 s)",
        title="Illustrative sensitivity",
    )
    for ax in axes:
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    fig.suptitle(
        "Synthetic kinetic baseline - not calibrated biological predictions",
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
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
