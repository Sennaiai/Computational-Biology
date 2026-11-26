"""Conditional adapter PK / cellular PD model for the separate FR-positive system.

Plasma exposure is specified in nM concentration equivalents, not clinical doses.
Local contact-site binding does not deplete the systemic reservoir. Receptor
totals describe a representative contact, not all receptors in the organism.
"""

from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.integrate import solve_ivp

from .adapter import AdapterBinding
from . import population
from .validation import number, time_grid


@dataclass(frozen=True)
class Kinetics:
    half_life_minutes: float = 60
    car_off_per_minute: float = 0.02
    target_off_per_minute: float = 0.2
    bridge_half_nM: float = 0.00005

    def __post_init__(self):
        for name, value in vars(self).items():
            number(value, name, positive=True)

    @property
    def elimination_per_day(self):
        return math.log(2) * 1440 / self.half_life_minutes


@dataclass(frozen=True)
class Protocol:
    name: str
    boluses: tuple = ()
    infusions: tuple = ()

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("protocol needs a name")
        for t, amount in self.boluses:
            number(t, "bolus time")
            number(amount, "concentration-equivalent bolus")
        for start, stop, rate in self.infusions:
            number(start, "infusion start")
            number(stop, "infusion stop", positive=True)
            number(rate, "infusion rate")
            if stop <= start:
                raise ValueError("infusion stop must exceed start")

    @property
    def total_input(self):
        return sum(amount for _, amount in self.boluses) + sum(
            (end - start) * rate for start, end, rate in self.infusions
        )


def protocols():
    weights = np.array([1, 2, 4, 8, 16, 16, 16], dtype=float)
    ramp = weights * (0.7 / weights.sum())
    return (
        Protocol("no_adapter"),
        Protocol("daily_bolus", tuple((float(day), 0.1) for day in range(7))),
        Protocol(
            "split_bolus",
            tuple(
                (day + fraction / 4, 0.025) for day in range(7) for fraction in range(4)
            ),
        ),
        Protocol(
            "ramped_bolus",
            tuple((float(day), float(amount)) for day, amount in enumerate(ramp)),
        ),
        Protocol("continuous", infusions=((0.0, 7.0, 0.1),)),
    )


def binding_rhs(adapter, occupied, binding, kinetics):
    """RA, FA and bridge; both binding orders obey the same equilibrium Kds."""
    ra, fa, bridge = occupied
    free_r = binding.car_total - ra - bridge
    free_f = binding.target_total - fa - bridge
    off_r = kinetics.car_off_per_minute * 1440
    off_f = kinetics.target_off_per_minute * 1440
    on_r = off_r / binding.car_kd
    on_f = off_f / binding.target_kd
    first_r = on_r * adapter * free_r - off_r * ra
    first_f = on_f * adapter * free_f - off_f * fa
    bridge_from_ra = on_f * ra * free_f - off_f * bridge
    bridge_from_fa = on_r * fa * free_r - off_r * bridge
    return np.array(
        [
            first_r - bridge_from_ra,
            first_f - bridge_from_fa,
            bridge_from_ra + bridge_from_fa,
        ]
    )


def coupled_rhs(t, state, parameters, binding, kinetics, infusion_rate):
    # 0:10 population states; 10 plasma A; 11:14 local RA/FA/B;
    # 14 cumulative systemic elimination; 15 cumulative input.
    bridge = state[13]
    activation = bridge / (kinetics.bridge_half_nM + bridge)
    signals = np.array([activation, 0, 0, 0])
    pop = population.rhs(t, state[:10], parameters, signals)
    elimination = kinetics.elimination_per_day * state[10]
    return np.r_[
        pop,
        infusion_rate - elimination,
        binding_rhs(state[10], state[11:14], binding, kinetics),
        elimination,
        infusion_rate,
    ]


@dataclass
class PKPDTrajectory:
    times: np.ndarray
    values: np.ndarray
    protocol: Protocol
    pk_ledger_error: float
    tumor_ledger_error: float
    minimum_site_pool: float
    evaluations: int
    elimination_per_day: float

    def summary(self):
        tumor = self.values[:, 2:6].sum(axis=1)
        return {
            "final_tumor_million_cells": float(tumor[-1]),
            "minimum_tumor_million_cells": float(tumor.min()),
            "peak_cytokine_proxy_arbitrary_units": float(self.values[:, 6].max()),
            "cytokine_proxy_auc_arbitrary_unit_days": float(self.values[-1, 7]),
            "final_car_t_million_cells": float(self.values[-1, :2].sum()),
            "cumulative_concentration_equivalent_input_nM": float(self.values[-1, 15]),
            "free_adapter_auc_nM_days": float(
                self.values[-1, 14] / self.elimination_per_day
            ),
            "pk_ledger_error_nM": self.pk_ledger_error,
            "tumor_ledger_error_million_cells": self.tumor_ledger_error,
            "minimum_receptor_pool_nM": self.minimum_site_pool,
            "rhs_evaluations": self.evaluations,
        }


def simulate(
    protocol,
    times,
    *,
    kinetics=None,
    binding=None,
    parameters=None,
    initial=None,
    rtol=2e-7,
    atol=1e-11,
):
    times = time_grid(times)
    kinetics = kinetics or Kinetics()
    binding = binding or AdapterBinding.published_affinities()
    parameters = parameters or replace(population.Parameters(), loss19=0, loss22=0)
    # Only the first (FR+) and last (FR-) tumor pools are used in this arm.
    # CD19/CD22 phenotype loss must not be smuggled into the FR architecture.
    if parameters.loss19 != 0 or parameters.loss22 != 0:
        raise ValueError("the separate FR model does not implement CD19/CD22 loss")
    initial = np.asarray(
        initial if initial is not None else [0.2, 0, 0.99, 0, 0, 0.01], dtype=float
    )
    if initial.shape != (6,) or not np.isfinite(initial).all() or np.any(initial < 0):
        raise ValueError("PKPD initial state needs six finite nonnegative populations")
    if initial[3] != 0 or initial[4] != 0 or initial[2:].sum() > parameters.capacity:
        raise ValueError("only FR-positive and FR-negative tumor pools are allowed")
    for value, name in ((rtol, "rtol"), (atol, "atol")):
        number(value, name, positive=True)
        if value > 1e-3:
            raise ValueError("PKPD tolerances must not exceed 1e-3")
    end = times[-1]
    if any(t >= end for t, _ in protocol.boluses) or any(
        stop > end for _, stop, _ in protocol.infusions
    ):
        raise ValueError("all input must finish before or at the simulation horizon")
    boundaries = sorted(
        {
            0.0,
            end,
            *(t for t, _ in protocol.boluses),
            *(t for pair in protocol.infusions for t in pair[:2]),
        }
    )
    state = np.r_[initial, np.zeros(10)]
    values = np.full((len(times), 16), np.nan)
    evaluations = 0
    internal_states = []
    for start, stop in zip(boundaries[:-1], boundaries[1:], strict=True):
        bolus = sum(amount for t, amount in protocol.boluses if t == start)
        state[10] += bolus
        state[15] += bolus
        rate = sum(
            rate
            for begin, finish, rate in protocol.infusions
            if begin <= start < finish
        )
        result = solve_ivp(
            coupled_rhs,
            (start, stop),
            state,
            args=(parameters, binding, kinetics, rate),
            method="Radau",
            rtol=rtol,
            atol=atol,
            dense_output=True,
            max_step=0.1,
        )
        if not result.success or not np.isfinite(result.y).all():
            raise RuntimeError("PKPD integration failed")
        evaluations += result.nfev
        internal_states.append(result.y.T)
        mask = (times >= start) & (times < stop)
        values[mask] = result.sol(times[mask]).T
        state = result.y[:, -1].copy()
    values[-1] = state
    checked = np.vstack([values, *internal_states])
    if not np.isfinite(checked).all() or checked.min() < -1e-8:
        raise RuntimeError("PKPD integration returned invalid states")
    pk_error = float(np.max(np.abs(checked[:, 10] + checked[:, 14] - checked[:, 15])))
    tumor_ledger = (
        checked[:, 2:6].sum(axis=1) - initial[2:].sum() - checked[:, 8] + checked[:, 9]
    )
    tumor_error = float(np.max(np.abs(tumor_ledger)))
    free_r = binding.car_total - checked[:, 11] - checked[:, 13]
    free_f = binding.target_total - checked[:, 12] - checked[:, 13]
    minimum_site_pool = float(min(free_r.min(), free_f.min(), checked[:, 11:14].min()))
    if pk_error > 1e-8 or tumor_error > 1e-7 or minimum_site_pool < -1e-10:
        raise RuntimeError("PKPD accounting or finite-receptor constraints failed")
    if not math.isclose(
        values[-1, 15], protocol.total_input, rel_tol=1e-8, abs_tol=1e-10
    ):
        raise RuntimeError("PKPD schedule did not deliver the declared input")
    trajectory = PKPDTrajectory(
        times,
        values,
        protocol,
        pk_error,
        tumor_error,
        minimum_site_pool,
        evaluations,
        kinetics.elimination_per_day,
    )
    return trajectory


def nondominated(outcomes, risk="peak_cytokine_proxy_arbitrary_units"):
    """Minimize final tumor and a declared activation proxy on a finite grid."""
    names = list(outcomes)
    points = np.array(
        [
            [outcomes[name]["final_tumor_million_cells"], outcomes[name][risk]]
            for name in names
        ]
    )
    if not np.isfinite(points).all() or np.any(points < 0):
        raise ValueError("tradeoff endpoints must be finite and nonnegative")
    return [
        name
        for i, name in enumerate(names)
        if not any(
            np.all(point <= points[i]) and np.any(point < points[i]) for point in points
        )
    ]
