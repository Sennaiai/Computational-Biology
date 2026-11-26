"""A reduced constant-cutting repair fit to the authors' public Psy1 example.

This is an independent, smaller model, not a reproduction of the publication's
time-dependent induction fit. Rates are effective assay parameters. Neither
read-count likelihood nor its profile represents Arabidopsis predictive accuracy.
"""

import hashlib
from importlib import resources
import io
import json

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from scipy.special import xlogy


def assay_data():
    root = resources.files("crispr_kinetics").joinpath("data")
    source = json.loads(root.joinpath("plant_repair_sources.json").read_text())
    for name, entry in source["files"].items():
        if (
            hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest()
            != entry["sha256"]
        ):
            raise ValueError("Published input checksum mismatch")
    raw = np.loadtxt(
        io.StringIO(root.joinpath("psy1_timecourse.tsv").read_text()), skiprows=1
    )
    error = np.loadtxt(io.StringIO(root.joinpath("psy1_error.tsv").read_text()))
    if raw.shape != (14, 5) or np.any(raw < 0) or np.any(raw != np.floor(raw)):
        raise ValueError("Unexpected assay counts")
    if (
        error.shape != (4, 4)
        or np.any(error < 0)
        or not np.allclose(error.sum(axis=1), 1)
    ):
        raise ValueError("Invalid observation-error matrix")
    times, counts = raw[:, 0], raw[:, 1:].astype(np.int64)
    if np.any(counts.sum(axis=1) <= 0):
        raise ValueError("Empty assay sample")
    return times, counts, error, source


def predict(theta, times, error):
    """theta: cut/h, precise/h, mutation/h, processed-DSB fraction, initial DSB.

    Latent states I, B, M use a column generator. Observed columns follow the
    authors' preparation code: intact, indels, direct DSB, processed DSB.
    """
    theta, times, error = np.asarray(theta), np.asarray(times), np.asarray(error)
    if any(x.dtype.kind not in "iuf" for x in (theta, times, error)):
        raise ValueError("Parameters, times and observation probabilities must be real")
    if (
        theta.shape != (5,)
        or not np.isfinite(theta).all()
        or np.any(theta < 0)
        or np.any(theta[3:] > 1)
    ):
        raise ValueError("Invalid repair parameters")
    if times.ndim != 1 or not np.isfinite(times).all() or np.any(times < 0):
        raise ValueError("Invalid times")
    if (
        error.shape != (4, 4)
        or np.any(error < 0)
        or not np.isfinite(error).all()
        or not np.allclose(error.sum(axis=1), 1, atol=1e-12)
    ):
        raise ValueError("Invalid observation matrix")
    cut, precise, mutation, processed, initial = theta
    q = np.array([[-cut, precise, 0], [cut, -precise - mutation, 0], [0, mutation, 0]])
    p0 = np.array([1 - initial, initial, 0])
    mapping = np.array(
        [error[0], (1 - processed) * error[2] + processed * error[3], error[1]]
    )
    unique, reverse = np.unique(times, return_inverse=True)
    values = np.array([expm(q * t) @ p0 for t in unique]) @ mapping
    if (
        not np.isfinite(values).all()
        or values.min() < -1e-12
        or not np.allclose(values.sum(axis=1), 1, atol=1e-10)
    ):
        raise ArithmeticError("Invalid assay probabilities")
    return values[reverse]


def negative_log_likelihood(theta, times, counts, error):
    values = predict(theta, times, error)
    counts = np.asarray(counts)
    if (
        counts.shape != values.shape
        or np.any(counts < 0)
        or not np.isfinite(counts).all()
    ):
        raise ValueError("Invalid observation counts")
    if np.any((values <= 0) & (counts > 0)):
        return float("inf")
    return float(-xlogy(counts, values).sum())


def fit(times, counts, error, *, precise_fixed=None, starts=8):
    """Bounded multistart MLE; caller supplies training rows only."""
    if isinstance(starts, bool) or not isinstance(starts, int) or not 1 <= starts <= 32:
        raise ValueError("Invalid multistart budget")
    if precise_fixed is not None and (
        not np.isfinite(precise_fixed) or precise_fixed < 0
    ):
        raise ValueError("Invalid fixed precise-repair rate")
    bounds = [(-12, 0), (-12, np.log(8)), (-12, 0), (1e-9, 1 - 1e-9), (1e-9, 0.2)]
    indices = [0, 1, 2, 3, 4] if precise_fixed is None else [0, 2, 3, 4]
    rng = np.random.default_rng(6103)

    def unpack(x):
        full = np.zeros(5)
        full[indices] = x
        full[:3] = np.exp(full[:3])
        if precise_fixed is not None:
            full[1] = precise_fixed
        return full

    def objective(x):
        return negative_log_likelihood(unpack(x), times, counts, error)

    successful, attempts = [], []
    for index in range(starts):
        physical = np.array([0.01, 0.08, 0.035, 0.2, 0.025])
        physical[:3] *= np.exp(rng.uniform(-2, 2, 3)) if index else 1
        full = np.r_[np.log(physical[:3]), physical[3:]]
        result = minimize(
            objective,
            full[indices],
            method="L-BFGS-B",
            bounds=[bounds[i] for i in indices],
            options={"maxiter": 400, "ftol": 1e-12, "gtol": 1e-5},
        )
        record = {
            "start": index,
            "success": bool(result.success),
            "message": str(result.message),
            "nll": float(result.fun),
            "parameters": unpack(result.x).tolist(),
        }
        attempts.append(record)
        if result.success and np.isfinite(result.fun):
            successful.append(record)
    if not successful:
        raise ArithmeticError("No calibration start converged")
    best = min(successful, key=lambda x: x["nll"])
    return {
        "parameters": best["parameters"],
        "training_nll": best["nll"],
        "attempts": attempts,
        "successful_starts": len(successful),
        "fitted_parameter_count": len(indices),
    }


def run_calibration():
    times, counts, error, source = assay_data()
    train, test = times < 72, times == 72
    full = fit(times[train], counts[train], error)
    ablation = fit(times[train], counts[train], error, precise_fixed=0)
    for result in (full, ablation):
        predicted = predict(result["parameters"], times[test], error)
        result["held_out_cross_entropy_per_read"] = negative_log_likelihood(
            result["parameters"], times[test], counts[test], error
        ) / int(counts[test].sum())
        result["held_out_mean_absolute_fraction_error"] = float(
            abs(
                predicted - counts[test] / counts[test].sum(axis=1, keepdims=True)
            ).mean()
        )
    profile = []
    for precise in (0, 0.003, 0.01, 0.03, 0.0791, 0.2, 0.6, 2, 8):
        result = fit(
            times[train], counts[train], error, precise_fixed=precise, starts=2
        )
        profile.append(
            {
                "precise_repair_per_h": precise,
                "training_nll": result["training_nll"],
                "parameters": result["parameters"],
            }
        )
    return {
        "model": "Constant-cutting I <-> DSB -> indel model; effective early induction; initial DSB fraction estimated from training data",
        "status": "example_assay_fit_completed",
        "source": source,
        "training_times_h": sorted(set(times[train].tolist())),
        "held_out_times_h": [72.0],
        "training_rows": int(train.sum()),
        "held_out_rows": int(test.sum()),
        "parameter_order": [
            "cut_per_h",
            "precise_repair_per_h",
            "mutation_repair_per_h",
            "processed_DSB_fraction",
            "initial_DSB_fraction",
        ],
        "full": full,
        "no_precise_repair": ablation,
        "precise_repair_profile": profile,
        "profile_interpretation": "Descriptive conditional read-count likelihood; not confidence intervals. Two replicates per time do not establish population uncertainty. Finite multistart optimization does not certify a global optimum.",
        "times_h": times.tolist(),
        "counts_in_source_order": counts.tolist(),
        "predicted_observed_fractions": predict(
            full["parameters"], times, error
        ).tolist(),
        "biological_transfer": "These are tomato example-assay estimates, not calibrated Arabidopsis or ABCG18/19 parameters. No cell-lineage or promoter data were fitted.",
    }
