"""Conditional antigen-escape model: populations are millions of cells, time days."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .validation import number, time_grid


PHENOTYPES = ("CD19+CD22+", "CD19+CD22-", "CD19-CD22+", "CD19-CD22-")
STATES = (
    "effector",
    "memory",
    *PHENOTYPES,
    "cytokine_proxy",
    "proxy_auc",
    "tumor_growth",
    "tumor_kills",
)


@dataclass(frozen=True)
class Parameters:
    growth: float = 0.15
    capacity: float = 100.0
    proliferation: float = 0.8
    reduction: float = 0.35
    memory_formation: float = 0.15
    memory_death: float = 0.005
    recall: float = 1.0
    suppression: float = 0.03
    killing: float = 2.0
    activation_half: float = 0.2
    loss19: float = 0.002
    loss22: float = 0.002
    cytokine_production: float = 1.0
    cytokine_clearance: float = 2.0

    def __post_init__(self):
        for name, value in vars(self).items():
            number(
                value,
                name,
                positive=name in ("capacity", "activation_half", "cytokine_clearance"),
            )
        if self.memory_formation > self.reduction:
            raise ValueError("memory formation is part of total effector reduction")


def loss_generator(loss19, loss22):
    number(loss19, "loss19")
    number(loss22, "loss22")
    return np.array(
        [
            [-loss19 - loss22, loss22, loss19, 0],
            [0, -loss19, 0, loss19],
            [0, 0, -loss22, loss22],
            [0, 0, 0, 0],
        ]
    )


def validate_signals(signals):
    signals = np.asarray(signals, dtype=float)
    if (
        signals.shape != (4,)
        or not np.isfinite(signals).all()
        or np.any((signals < 0) | (signals > 1))
    ):
        raise ValueError("signals must contain four recognition probabilities in [0,1]")
    if signals[3] != 0:
        raise ValueError(
            "double-negative cells cannot be recognized by this antigen-specific model"
        )
    return signals


def rhs(t, state, parameters, signals):
    p = parameters
    e, m = state[:2]
    tumor = state[2:6]
    total = tumor.sum()
    recognized = signals @ tumor
    activation = recognized / (p.activation_half + recognized)
    recall = p.recall * recognized * m
    growth = p.growth * tumor * (1 - total / p.capacity)
    kills = p.killing * e * signals * tumor
    changes = loss_generator(p.loss19, p.loss22).T @ tumor
    cytokine = state[6]
    return np.concatenate(
        (
            [
                (p.proliferation * activation - p.reduction - p.suppression * total) * e
                + recall,
                p.memory_formation * e - recall - p.memory_death * m,
            ],
            growth - kills + changes,
            [
                p.cytokine_production * e * activation
                - p.cytokine_clearance * cytokine,
                cytokine,
                growth.sum(),
                kills.sum(),
            ],
        )
    )


@dataclass
class Trajectory:
    times: np.ndarray
    values: np.ndarray
    evaluations: int
    ledger_error: float
    minimum: float

    @property
    def tumor(self):
        return self.values[:, 2:6].sum(axis=1)

    def summary(self):
        return {
            "final_tumor_million_cells": float(self.tumor[-1]),
            "minimum_tumor_million_cells": float(self.tumor.min()),
            "final_double_negative_fraction": float(self.values[-1, 5] / self.tumor[-1])
            if self.tumor[-1] > 0
            else None,
            "final_effector_million_cells": float(self.values[-1, 0]),
            "final_memory_million_cells": float(self.values[-1, 1]),
            "peak_cytokine_proxy_arbitrary_units": float(self.values[:, 6].max()),
            "cytokine_proxy_auc_arbitrary_unit_days": float(self.values[-1, 7]),
            "tumor_ledger_max_error_million_cells": self.ledger_error,
            "minimum_integrated_state": self.minimum,
            "rhs_evaluations": self.evaluations,
        }


def simulate(
    initial, times, signals, parameters=None, *, method="DOP853", rtol=1e-8, atol=1e-10
):
    """Initial E, M and four tumor pools; counters and cytokine start at zero."""
    times = time_grid(times)
    signals = validate_signals(signals)
    p = parameters or Parameters()
    initial = np.asarray(initial, dtype=float)
    if initial.shape != (6,) or not np.isfinite(initial).all() or np.any(initial < 0):
        raise ValueError("initial must contain six finite nonnegative populations")
    if initial[2:].sum() > p.capacity:
        raise ValueError("initial tumor must not exceed carrying capacity")
    if method not in ("DOP853", "Radau"):
        raise ValueError("method must be DOP853 or Radau")
    for value, name in ((rtol, "rtol"), (atol, "atol")):
        number(value, name, positive=True)
        if value > 1e-3:
            raise ValueError("numerical tolerances must not exceed 1e-3")
    result = solve_ivp(
        rhs,
        (0, times[-1]),
        np.r_[initial, np.zeros(4)],
        args=(p, signals),
        method=method,
        t_eval=times,
        rtol=rtol,
        atol=atol,
        max_step=0.5,
    )
    if (
        not result.success
        or result.y.shape != (10, times.size)
        or not np.isfinite(result.y).all()
    ):
        raise RuntimeError("population integration failed")
    values = result.y.T
    minimum = float(values.min())
    ledger = (
        values[:, 2:6].sum(axis=1) - initial[2:].sum() - values[:, 8] + values[:, 9]
    )
    ledger_error = float(np.max(np.abs(ledger)))
    if minimum < -1e-8 or ledger_error > 1e-7 * max(1, p.capacity):
        raise RuntimeError(
            "population integration violated positivity or tumor accounting"
        )
    if np.max(values[:, 2:6].sum(axis=1)) > p.capacity + 1e-6:
        raise RuntimeError(
            "population integration violated the carrying-capacity bound"
        )
    return Trajectory(times, values, result.nfev, ledger_error, minimum)
