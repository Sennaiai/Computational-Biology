"""Single-target, single-cleavage Cas9 kinetics with fixed free-RNP concentration.

Column convention: dp/dt = Q @ p; Q[destination, source] is a transition
rate. Repair is a separate, lumped layer after cleavage. Precisely repaired,
indel and HDR outcomes are terminal in this baseline: re-cutting, enzyme
depletion, chromatin and delivery are not modeled.
"""

from dataclasses import asdict, dataclass
import math
import numbers

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

STATES = ("unbound", "PAM_bound", "R_loop", "cleaved", "precise_repair", "indel", "HDR")


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a finite real number")
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")


@dataclass(frozen=True)
class Rates:
    """Explicit units; defaults are synthetic demonstrations, not fitted rates."""

    association_per_nM_s: float = 0.01
    free_rnp_nM: float = 10.0
    dissociation_per_s: float = 0.03
    r_loop_per_s: float = 0.06
    unwind_per_s: float = 0.015
    cleavage_per_s: float = 0.025
    repair_per_s: float = 0.01
    precise_fraction: float = 0.4
    indel_fraction: float = 0.5
    hdr_fraction: float = 0.1

    def __post_init__(self):
        for name, value in asdict(self).items():
            _real(value, name)
        branches = [self.precise_fraction, self.indel_fraction, self.hdr_fraction]
        if any(v > 1 for v in branches) or not math.isclose(
            sum(branches), 1, rel_tol=0, abs_tol=1e-12
        ):
            raise ValueError("Repair fractions must be probabilities summing to one")
        if not math.isfinite(self.association_per_nM_s * self.free_rnp_nM):
            raise ValueError("Pseudo-first-order association rate overflowed")


def generator(rates=Rates()):
    if not isinstance(rates, Rates):
        raise ValueError("Expected a validated Rates configuration")
    q = np.zeros((7, 7))
    transitions = [
        (0, 1, rates.association_per_nM_s * rates.free_rnp_nM),
        (1, 0, rates.dissociation_per_s),
        (1, 2, rates.r_loop_per_s),
        (2, 1, rates.unwind_per_s),
        (2, 3, rates.cleavage_per_s),
        (3, 4, rates.repair_per_s * rates.precise_fraction),
        (3, 5, rates.repair_per_s * rates.indel_fraction),
        (3, 6, rates.repair_per_s * rates.hdr_fraction),
    ]
    for source, destination, rate in transitions:
        q[destination, source] += rate
    np.fill_diagonal(q, -q.sum(axis=0))
    if not np.isfinite(q).all():
        raise ValueError("Generator rates overflowed")
    return q


def _times(values):
    raw = np.asarray(values)
    if raw.dtype.kind not in "iuf" or np.iscomplexobj(raw):
        raise ValueError("Times must be real numeric values")
    times = np.asarray(values, dtype=float)
    if times.ndim != 1 or not 1 <= len(times) <= 10000 or not np.isfinite(times).all():
        raise ValueError("Expected a nonempty bounded one-dimensional time grid")
    if np.any(times < 0) or np.any(np.diff(times) <= 0):
        raise ValueError(
            "Times must be nonnegative and strictly increasing, in seconds"
        )
    return times


def _initial(values):
    if values is None:
        return np.eye(7)[0]
    raw = np.asarray(values)
    if raw.dtype.kind not in "iuf" or np.iscomplexobj(raw):
        raise ValueError("Initial probabilities must be real numeric values")
    p = np.asarray(values, dtype=float)
    if p.shape != (7,) or not np.isfinite(p).all() or np.any(p < 0) or np.any(p > 1):
        raise ValueError("Expected seven finite initial probabilities in [0,1]")
    if not math.isclose(float(p.sum()), 1, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Initial probability must sum to one")
    return p.copy()


def diagnostics(probabilities):
    p = np.asarray(probabilities)
    if (
        p.dtype.kind not in "iuf"
        or p.ndim != 2
        or p.shape[0] == 0
        or p.shape[1] != 7
        or not np.isfinite(p).all()
    ):
        raise ArithmeticError("Non-finite or malformed probability solution")
    result = {
        "max_mass_error": float(np.max(abs(p.sum(axis=1) - 1))),
        "min_probability": float(p.min()),
        "max_probability": float(p.max()),
    }
    if result["max_mass_error"] > 1e-8 or p.min() < -1e-10 or p.max() > 1 + 1e-10:
        raise ArithmeticError("Numerical solution violates probability invariants")
    return result


def solve(rates, times_s, *, method="expm", initial=None):
    """Return time-by-state probabilities without clipping or renormalization."""
    times, p0, q = _times(times_s), _initial(initial), generator(rates)
    if method == "expm":
        p = np.array([expm(q * time) @ p0 for time in times])
    elif method == "ode":
        if times[-1] == 0:
            p = p0[None, :]
        else:
            fit = solve_ivp(
                lambda _t, values: q @ values,
                (0, times[-1]),
                p0,
                t_eval=times,
                method="DOP853",
                rtol=1e-10,
                atol=1e-12,
            )
            if not fit.success:
                raise ArithmeticError("ODE integration failed: " + fit.message)
            p = fit.y.T
    else:
        raise ValueError("Method must be expm or ode")
    diagnostics(p)
    return p


def readouts(probabilities):
    diagnostics(probabilities)
    p = np.asarray(probabilities)
    return {
        "cleaved_ever": p[:, 3:].sum(axis=1),
        "currently_cleaved_unrepaired": p[:, 3],
        "repaired_total": p[:, 4:].sum(axis=1),
        "sequence_altered_repair_proxy": p[:, 5:].sum(axis=1),
    }


def first_cleavage_time(rates):
    """Mean first-passage time from unbound; None means cleavage is unreachable."""
    q = generator(rates)
    if (
        rates.association_per_nM_s * rates.free_rnp_nM == 0
        or rates.r_loop_per_s == 0
        or rates.cleavage_per_s == 0
    ):
        return None
    backward = -q[:3, :3].T
    condition = float(np.linalg.cond(backward))
    if not np.isfinite(condition) or condition > 1e12:
        raise ArithmeticError(
            "First-passage system is too ill-conditioned for this numerical study"
        )
    means = np.linalg.solve(backward, np.ones(3))
    if not np.isfinite(means).all() or np.any(means <= 0):
        raise ArithmeticError("Invalid first-passage solution")
    return {
        "mean_seconds": float(means[0]),
        "condition_number": condition,
        "linear_system_dimension": 3,
        "generator_is_Hermitian": bool(np.allclose(q, q.T)),
    }


def simulate(rates, times_s, *, paths=4000, seed=7, initial=None, event_budget=100000):
    """Exact Gillespie trajectories for the same CTMC; return aggregate counts only."""
    times, p0, q = _times(times_s), _initial(initial), generator(rates)
    for name, value, low, high in (
        ("paths", paths, 1, 100000),
        ("seed", seed, 0, 2**32 - 1),
        ("event_budget", event_budget, 1, 1000000),
    ):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or not low <= value <= high
        ):
            raise ValueError(f"Invalid integer {name}")
    if paths * len(times) > 20000000:
        raise ValueError(
            "Trajectory/time grid exceeds the experiment observation budget"
        )
    rng = np.random.default_rng(seed)
    counts = np.zeros((len(times), 7), dtype=np.int64)
    for _ in range(paths):
        state, now, index, events = int(rng.choice(7, p=p0)), 0.0, 0, 0
        while index < len(times):
            exit_rate = -q[state, state]
            if exit_rate == 0:
                counts[index:, state] += 1
                break
            event_time = now + float(rng.exponential(1 / exit_rate))
            stop = int(np.searchsorted(times, event_time, side="left"))
            counts[index:stop, state] += 1
            index = stop
            if index == len(times):
                break
            if events >= event_budget:
                raise RuntimeError(
                    "SSA event budget exhausted; no partial result returned"
                )
            hazards = q[:, state].copy()
            hazards[state] = 0
            state = int(rng.choice(7, p=hazards / exit_rate))
            now, events = event_time, events + 1
    if not np.all(counts.sum(axis=1) == paths):
        raise ArithmeticError("SSA trajectory ledger lost probability mass")
    return counts
