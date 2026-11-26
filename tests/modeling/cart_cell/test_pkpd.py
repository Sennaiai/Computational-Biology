"""Analytical PK, binding detailed balance and input-matched protocol checks."""

from dataclasses import replace
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

from cart_dynamics.adapter import AdapterBinding
from cart_dynamics import pkpd, population


class PKPDTests(unittest.TestCase):
    def test_binding_equilibrium_obeys_dynamic_mass_action(self):
        binding = AdapterBinding.published_affinities()
        for adapter in (0.001, 0.1, 10):
            s = binding.species(adapter)
            occupied = [s["car_adapter"], s["target_adapter"], s["bridge"]]
            assert_allclose(
                pkpd.binding_rhs(adapter, occupied, binding, pkpd.Kinetics()),
                0,
                atol=1e-13,
            )

    def test_same_kd_does_not_fix_residence_time(self):
        binding = AdapterBinding()
        eq = binding.species(0.1)
        start = [eq["car_adapter"], eq["target_adapter"], eq["bridge"]]
        slow = pkpd.Kinetics(car_off_per_minute=0.002, target_off_per_minute=0.02)
        fast = replace(slow, car_off_per_minute=0.2, target_off_per_minute=2)

        def remaining(kinetics):
            result = solve_ivp(
                lambda t, y: pkpd.binding_rhs(0, y, binding, kinetics),
                (0, 0.02),
                start,
                method="Radau",
                rtol=1e-10,
                atol=1e-13,
            )
            self.assertTrue(result.success)
            return result.y[2, -1]

        self.assertGreater(remaining(slow), remaining(fast))

    def test_bolus_pk_is_exact_exponential_and_mass_balances(self):
        times = np.linspace(0, 1, 51)
        k = pkpd.Kinetics()
        result = pkpd.simulate(
            pkpd.Protocol("one", boluses=((0, 0.1),)), times, kinetics=k
        )
        expected = 0.1 * np.exp(-k.elimination_per_day * times)
        assert_allclose(result.values[:, 10], expected, atol=2e-10, rtol=2e-6)
        assert_allclose(result.values[:, 10] + result.values[:, 14], 0.1, atol=1e-12)

    def test_continuous_pk_matches_closed_form_at_start_and_stop(self):
        times = np.linspace(0, 2, 101)
        k = pkpd.Kinetics().elimination_per_day
        p = pkpd.Protocol("infusion", infusions=((0, 1, 0.1),))
        result = pkpd.simulate(p, times)
        expected = np.where(
            times <= 1,
            0.1 / k * (1 - np.exp(-k * times)),
            0.1 / k * (1 - np.exp(-k)) * np.exp(-k * (times - 1)),
        )
        assert_allclose(result.values[:, 10], expected, atol=1e-10, rtol=3e-6)

    def test_interior_bolus_boundary_is_not_skipped(self):
        result = pkpd.simulate(
            pkpd.Protocol("late", boluses=((0.5, 0.1),)), [0, 0.25, 0.5, 1]
        )
        assert_allclose(result.values[:2, 10], 0)
        self.assertAlmostEqual(result.values[2, 10], 0.1)

    def test_no_adapter_reduces_to_logistic_growth(self):
        times = np.linspace(0, 14, 51)
        result = pkpd.simulate(pkpd.Protocol("none"), times)
        p = population.Parameters()
        expected = 1 / (
            1 / p.capacity + (1 - 1 / p.capacity) * np.exp(-p.growth * times)
        )
        assert_allclose(result.values[:, 2:6].sum(axis=1), expected, rtol=1e-7)
        assert_allclose(result.values[:, [6, 7, 9, 10, 11, 12, 13]], 0)

    def test_no_antigen_does_not_generate_cytokine_or_killing(self):
        result = pkpd.simulate(
            pkpd.Protocol("one", boluses=((0, 0.1),)),
            [0, 0.5, 1],
            initial=[0.2, 0, 0, 0, 0, 1],
        )
        assert_allclose(result.values[:, [6, 7, 9]], 0)

    def test_matching_input_and_auc_for_all_active_protocols(self):
        kinetics = pkpd.Kinetics()
        for protocol in pkpd.protocols()[1:]:
            self.assertAlmostEqual(protocol.total_input, 0.7)
            result = pkpd.simulate(protocol, np.linspace(0, 14, 141), kinetics=kinetics)
            summary = result.summary()
            self.assertAlmostEqual(
                summary["free_adapter_auc_nM_days"],
                0.7 / kinetics.elimination_per_day,
                places=9,
            )
            self.assertLess(summary["pk_ledger_error_nM"], 1e-12)
            self.assertLess(summary["tumor_ledger_error_million_cells"], 1e-10)
            self.assertGreaterEqual(summary["minimum_receptor_pool_nM"], -1e-12)

    def test_tighter_solver_tolerances_preserve_pd_endpoints(self):
        protocol = pkpd.Protocol("one", boluses=((0, 0.1),))
        first = pkpd.simulate(protocol, np.linspace(0, 2, 101))
        second = pkpd.simulate(protocol, np.linspace(0, 2, 101), rtol=1e-9, atol=1e-13)
        assert_allclose(first.values, second.values, rtol=3e-5, atol=1e-8)

    def test_protocol_contracts(self):
        for kwargs in (
            {"boluses": ((-1, 0.1),)},
            {"infusions": ((1, 0, 1),)},
            {"boluses": ((0, math.nan),)},
        ):
            with self.assertRaises(ValueError):
                pkpd.Protocol("invalid", **kwargs)
        with self.assertRaises(ValueError):
            pkpd.Kinetics(half_life_minutes=0)
        with self.assertRaises(ValueError):
            pkpd.simulate(pkpd.Protocol("late", boluses=((1, 0.1),)), [0, 1])
        with self.assertRaises(ValueError):
            pkpd.simulate(
                pkpd.Protocol("none"), [0, 1], parameters=population.Parameters()
            )

    def test_pareto_rule_keeps_tradeoffs_and_rejects_dominated_points(self):
        def point(tumor, cytokine):
            return {
                "final_tumor_million_cells": tumor,
                "peak_cytokine_proxy_arbitrary_units": cytokine,
            }

        outcomes = {
            "low": point(5, 1),
            "high": point(1, 5),
            "bad": point(6, 2),
            "same": point(5, 1),
        }
        self.assertEqual(pkpd.nondominated(outcomes), ["low", "high", "same"])


if __name__ == "__main__":
    unittest.main()
