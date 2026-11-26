"""Independent implementation of Barros et al. (2021), equations 1–3."""

from dataclasses import dataclass
from importlib.resources import files
import json

import numpy as np
from scipy.integrate import solve_ivp

from .validation import number, time_grid


@dataclass(frozen=True)
class ReferenceParameters:
    proliferation: float
    reduction: float
    memory_formation: float
    recall: float
    suppression: float
    memory_death: float
    growth: float
    inverse_capacity: float
    killing: float

    def __post_init__(self):
        for name, value in vars(self).items():
            number(value, name, positive=True)
        if (
            not self.reduction - self.memory_formation
            < self.proliferation
            < self.reduction
        ):
            raise ValueError("reference requires rho-epsilon < phi < rho")
        if self.memory_formation > self.reduction:
            raise ValueError("memory formation cannot exceed total effector reduction")


def published_parameters():
    data = json.loads(files("cart_dynamics").joinpath("data/sources.json").read_text())
    return ReferenceParameters(**data["barros2021"]["hdlm2_parameters"])


def rhs(t, state, parameters):
    """Time in days, E/M/T in cells; coefficients in Table 2's original units."""
    e, m, tumor = state
    p = parameters
    recall = p.recall * tumor * m
    return np.array(
        [
            (p.proliferation - p.reduction - p.suppression * tumor) * e + recall,
            p.memory_formation * e - recall - p.memory_death * m,
            p.growth * tumor * (1 - p.inverse_capacity * tumor) - p.killing * e * tumor,
        ]
    )


def coexistence_equilibria(parameters=None):
    p = parameters or published_parameters()
    coefficients = [
        -p.suppression * p.recall,
        (p.proliferation - p.reduction + p.memory_formation) * p.recall
        - p.suppression * p.memory_death,
        (p.proliferation - p.reduction) * p.memory_death,
    ]
    equilibria = []
    for tumor in np.roots(coefficients):
        if np.isreal(tumor) and 0 < tumor.real < 1 / p.inverse_capacity:
            tumor = float(tumor.real)
            e = p.growth * (1 - p.inverse_capacity * tumor) / p.killing
            m = p.memory_formation * e / (p.recall * tumor + p.memory_death)
            equilibria.append([e, m, tumor])
    return np.asarray(sorted(equilibria, key=lambda row: row[2]))


def simulate(initial, times, parameters=None):
    times = time_grid(times)
    initial = np.asarray(initial, dtype=float)
    if initial.shape != (3,) or not np.isfinite(initial).all() or np.any(initial < 0):
        raise ValueError(
            "reference initial state must contain three nonnegative populations"
        )
    p = parameters or published_parameters()
    # Integrate in millions of cells for balanced absolute error control.
    scale = 1_000_000.0
    result = solve_ivp(
        lambda t, y: rhs(t, y * scale, p) / scale,
        (0, times[-1]),
        initial / scale,
        t_eval=times,
        method="Radau",
        rtol=1e-9,
        atol=1e-11,
    )
    if (
        not result.success
        or not np.isfinite(result.y).all()
        or np.min(result.y) < -1e-8
    ):
        raise RuntimeError(
            "published reference integration failed its numerical checks"
        )
    return result.y.T * scale
