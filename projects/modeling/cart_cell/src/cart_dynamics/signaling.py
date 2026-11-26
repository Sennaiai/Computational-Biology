"""Finite receptor-state model and explicit idealized dual-recognition rules."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from .validation import integer, number


@dataclass(frozen=True)
class Receptor:
    """Effective contact rates in min^-1; density is a dimensionless input."""

    association: float = 2.0
    dissociation: float = 0.2
    proofreading: float = 1.0
    steps: int = 2

    def __post_init__(self):
        number(self.association, "association")
        number(self.dissociation, "dissociation", positive=True)
        number(self.proofreading, "proofreading", positive=True)
        integer(self.steps, "steps", 1, 8)

    def generator(self, density):
        number(density, "density")
        matrix = np.zeros((self.steps + 2, self.steps + 2))
        matrix[0, 1] = self.association * density
        for j in range(1, self.steps + 2):
            matrix[j, 0] = self.dissociation
            if j < self.steps + 1:
                matrix[j, j + 1] = self.proofreading
        matrix -= np.diag(matrix.sum(axis=1))
        if not np.isfinite(matrix).all():
            raise ValueError("receptor rates overflowed")
        return matrix

    def equilibrium_active(self, density):
        matrix = self.generator(density)
        on = matrix[0, 1]
        return (on / (on + self.dissociation)) * (
            self.proofreading / (self.proofreading + self.dissociation)
        ) ** self.steps

    def probabilities(self, density, minutes, initial=None):
        number(minutes, "minutes")
        matrix = self.generator(density)
        if initial is None:
            initial = np.eye(self.steps + 2)[0]
        initial = np.asarray(initial, dtype=float)
        if (
            initial.shape != (self.steps + 2,)
            or not np.isfinite(initial).all()
            or np.any(initial < 0)
            or not np.isclose(initial.sum(), 1.0, rtol=0, atol=1e-12)
        ):
            raise ValueError("initial must be a normalized receptor probability vector")
        return expm(matrix.T * minutes) @ initial

    def sample_counts(self, density, minutes, receptors, seed, max_events=1_000_000):
        """Exact finite-pool Gillespie process, initially all unbound."""
        number(minutes, "minutes")
        integer(receptors, "receptors", 1, 1_000_000)
        integer(seed, "seed", 0, 2**63 - 1)
        integer(max_events, "max_events", 1, 100_000_000)
        matrix = self.generator(density)
        transitions = np.argwhere(matrix > 0)
        counts = np.zeros(matrix.shape[0], dtype=int)
        counts[0] = receptors
        rng = np.random.default_rng(seed)
        t = 0.0
        for _ in range(max_events):
            propensities = np.array([counts[i] * matrix[i, j] for i, j in transitions])
            total = propensities.sum()
            if total == 0:
                return counts
            t += rng.exponential(1 / total)
            if t > minutes:
                return counts
            index = np.searchsorted(
                np.cumsum(propensities), rng.random() * total, side="right"
            )
            i, j = transitions[min(index, len(transitions) - 1)]
            counts[i] -= 1
            counts[j] += 1
        raise RuntimeError("receptor simulation exhausted its event budget")


def recognition(first, second, policy):
    """Independent-opportunity OR/AND approximation, not a CAR construct fit."""
    first, second = np.broadcast_arrays(
        np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    )
    if not (np.isfinite(first).all() and np.isfinite(second).all()):
        raise ValueError("recognition probabilities must be finite")
    if np.any((first < 0) | (first > 1) | (second < 0) | (second > 1)):
        raise ValueError("recognition probabilities must lie in [0, 1]")
    if policy == "cd19":
        return first
    if policy == "cd22":
        return second
    if policy == "or":
        return first + second - first * second
    if policy == "and":
        return first * second
    if policy == "none":
        return np.zeros_like(first)
    raise ValueError("policy must be cd19, cd22, or, and or none")


def phenotype_signals(policy, receptor=Receptor()):
    """Order: double-positive, CD19-only, CD22-only, double-negative."""
    active = receptor.equilibrium_active(1.0)
    return recognition([active, active, 0, 0], [active, 0, active, 0], policy)
