"""Analytical limits and independent scientific contracts for CAR-T components."""

from dataclasses import replace
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm

from cart_dynamics.adapter import AdapterBinding, free_adapter_profile
from cart_dynamics import population, reference
from cart_dynamics.signaling import Receptor, phenotype_signals, recognition


class ReceptorTests(unittest.TestCase):
    def test_generator_and_stationary_closed_form(self):
        for steps in (1, 2, 5):
            receptor = Receptor(steps=steps)
            matrix = receptor.generator(0.7)
            assert_allclose(matrix.sum(axis=1), 0, atol=1e-15)
            self.assertTrue(np.all((matrix - np.diag(matrix.diagonal())) >= 0))
            p = receptor.probabilities(0.7, 200)
            assert_allclose(p.sum(), 1, atol=1e-13)
            self.assertGreaterEqual(p.min(), 0)
            assert_allclose(matrix.T @ p, 0, atol=1e-13)
            self.assertAlmostEqual(p[-1], receptor.equilibrium_active(0.7), places=13)

    def test_binding_occupancy_matches_two_state_solution(self):
        receptor = Receptor()
        density, t = 0.8, 3.0
        on = receptor.association * density
        bound = (
            on
            / (on + receptor.dissociation)
            * (1 - np.exp(-(on + receptor.dissociation) * t))
        )
        self.assertAlmostEqual(
            receptor.probabilities(density, t)[1:].sum(), bound, places=13
        )

    def test_zero_antigen_and_withdrawal(self):
        receptor = Receptor()
        self.assertEqual(receptor.equilibrium_active(0), 0)
        self.assertEqual(receptor.probabilities(0, 100)[-1], 0)
        initial = receptor.probabilities(1, 200)
        final = receptor.probabilities(0, 10, initial)
        self.assertAlmostEqual(
            final[1:].sum(), initial[1:].sum() * np.exp(-2), places=13
        )

    def test_proofreading_penalizes_short_residence(self):
        slow = Receptor(dissociation=0.2)
        fast = replace(slow, dissociation=2)
        self.assertLess(fast.equilibrium_active(1), slow.equilibrium_active(1))

    def test_gillespie_agrees_with_matrix_exponential(self):
        receptor = Receptor()
        n = 4000
        first = receptor.sample_counts(0.6, 2, n, 238)
        second = receptor.sample_counts(0.6, 2, n, 238)
        assert_allclose(first, second)
        self.assertEqual(first.sum(), n)
        self.assertGreaterEqual(first.min(), 0)
        bound = np.sqrt(np.log(2 * len(first) / 1e-6) / (2 * n))
        self.assertLess(
            np.max(np.abs(first / n - receptor.probabilities(0.6, 2))), bound
        )

    def test_event_budget_fails_loudly(self):
        with self.assertRaises(RuntimeError):
            Receptor().sample_counts(1, 10, 100, 8, max_events=1)

    def test_recognition_truth_tables_and_null_target(self):
        a, b = [1, 1, 0, 0], [1, 0, 1, 0]
        assert_allclose(recognition(a, b, "or"), [1, 1, 1, 0])
        assert_allclose(recognition(a, b, "and"), [1, 0, 0, 0])
        for policy in ("cd19", "cd22", "or", "and", "none"):
            self.assertEqual(phenotype_signals(policy)[3], 0)

    def test_invalid_inputs(self):
        for kwargs in (
            {"steps": 1.5},
            {"steps": True},
            {"dissociation": 0},
            {"association": -1},
        ):
            with self.assertRaises(ValueError):
                Receptor(**kwargs)
        for density in (-1, np.inf, np.nan, True):
            with self.assertRaises(ValueError):
                Receptor().generator(density)
        with self.assertRaises(ValueError):
            Receptor().probabilities(1, 2, [0, 0, 0, 0])
        with self.assertRaises(ValueError):
            recognition(0.1, np.nan, "or")
        with self.assertRaises(ValueError):
            recognition(1.1, 0.2, "or")
        with self.assertRaises(ValueError):
            recognition(0.1, 0.2, "unknown")


class PublishedReferenceTests(unittest.TestCase):
    def test_table_values_and_reported_equilibria(self):
        p = reference.published_parameters()
        self.assertEqual(p.growth, 0.05650026)
        self.assertEqual(p.killing, 3.715843e-6)
        equilibria = reference.coexistence_equilibria(p)
        # Table and text differ by 0.202% at the lower tumor root; keep that discrepancy visible.
        expected = [[15205, 197519, 1089], [15205, 263, 1442522]]
        assert_allclose(equilibria, expected, rtol=0.003)
        for state in equilibria:
            assert_allclose(reference.rhs(0, state, p), 0, atol=1e-8)

    def test_tumor_only_is_logistic(self):
        p = reference.published_parameters()
        t = np.linspace(0, 100, 101)
        start = 2e5
        values = reference.simulate([0, 0, start], t, p)
        expected = (
            start
            * np.exp(p.growth * t)
            / (1 + p.inverse_capacity * start * np.expm1(p.growth * t))
        )
        assert_allclose(values[:, 2], expected, rtol=1e-8)
        assert_allclose(values[:, :2], 0)

    def test_no_tumor_effector_memory_closed_form(self):
        p = reference.published_parameters()
        t = np.linspace(0, 50, 31)
        e0, m0 = 2e6, 3e5
        rate = p.proliferation - p.reduction
        e = e0 * np.exp(rate * t)
        m = m0 * np.exp(-p.memory_death * t) + p.memory_formation * e0 * (
            np.exp(rate * t) - np.exp(-p.memory_death * t)
        ) / (rate + p.memory_death)
        values = reference.simulate([e0, m0, 0], t, p)
        assert_allclose(values[:, :2], np.c_[e, m], rtol=2e-8)
        assert_allclose(values[:, 2], 0)


class PopulationTests(unittest.TestCase):
    def setUp(self):
        self.times = np.linspace(0, 30, 121)
        self.initial = [0.2, 0, 0.7, 0.15, 0.14, 0.01]

    def test_no_car_total_tumor_is_logistic_despite_switching(self):
        p = population.Parameters()
        values = population.simulate(
            [0, 0, *self.initial[2:]], self.times, phenotype_signals("or"), p
        )
        expected = 1 / (
            1 / p.capacity + (1 - 1 / p.capacity) * np.exp(-p.growth * self.times)
        )
        assert_allclose(values.tumor, expected, rtol=1e-8)
        assert_allclose(values.values[:, [0, 1, 6, 7, 9]], 0)
        self.assertLess(values.ledger_error, 1e-10)

    def test_phenotype_loss_has_exact_markov_limit(self):
        p = population.Parameters(growth=0, killing=0)
        t = self.times
        start = np.array([0.7, 0.15, 0.14, 0.01])
        values = population.simulate([0, 0, *start], t, phenotype_signals("none"), p)
        q = population.loss_generator(p.loss19, p.loss22)
        expected = np.array([expm(q.T * time) @ start for time in t])
        assert_allclose(values.values[:, 2:6], expected, atol=1e-12)
        assert_allclose(values.tumor, 1, atol=1e-12)

    def test_double_negative_cannot_be_killed_or_activate(self):
        values = population.simulate(
            [0.2, 0, 0, 0, 0, 1], self.times, phenotype_signals("or")
        )
        control = population.simulate(
            [0, 0, 0, 0, 0, 1], self.times, phenotype_signals("none")
        )
        assert_allclose(values.tumor, control.tumor, rtol=1e-8)
        assert_allclose(values.values[:, [6, 7, 9]], 0)

    def test_zero_killing_does_not_fake_cytotoxicity(self):
        p = population.Parameters(killing=0)
        values = population.simulate(
            self.initial, self.times, phenotype_signals("or"), p
        )
        control = population.simulate(
            [0, 0, *self.initial[2:]], self.times, phenotype_signals("none"), p
        )
        assert_allclose(values.tumor, control.tumor, rtol=1e-8)
        assert_allclose(values.values[:, 9], 0)
        self.assertGreater(values.values[:, 6].max(), 0)

    def test_solvers_agree_and_accounting_closes(self):
        signal = phenotype_signals("or")
        first = population.simulate(self.initial, self.times, signal)
        second = population.simulate(
            self.initial, self.times, signal, method="Radau", rtol=1e-10, atol=1e-12
        )
        assert_allclose(first.values, second.values, rtol=3e-6, atol=2e-8)
        self.assertGreaterEqual(first.minimum, -1e-10)
        self.assertLess(first.ledger_error, 1e-10)

    def test_nonnegative_boundary_derivatives(self):
        p = population.Parameters()
        signal = phenotype_signals("or")
        for index in range(10):
            state = np.array([0.2, 0.1, 0.7, 0.15, 0.14, 0.01, 0.2, 1, 1, 1])
            state[index] = 0
            self.assertGreaterEqual(population.rhs(0, state, p, signal)[index], 0)

    def test_production_zero_gives_zero_cytokine(self):
        values = population.simulate(
            self.initial,
            self.times,
            phenotype_signals("or"),
            population.Parameters(cytokine_production=0),
        )
        assert_allclose(values.values[:, 6:8], 0)

    def test_contracts_reject_invalid_models(self):
        for kwargs in (
            {"growth": -1},
            {"capacity": 0},
            {"memory_formation": 0.5},
            {"loss19": np.inf},
        ):
            with self.assertRaises(ValueError):
                population.Parameters(**kwargs)
        for times in ([1, 2], [0, 0], [0, np.nan], [0]):
            with self.assertRaises(ValueError):
                population.simulate(self.initial, times, phenotype_signals("or"))
        for signal in ([1, 1, 1, 1], [np.nan, 0, 0, 0], [1, 0, 0], [1.1, 0, 0, 0]):
            with self.assertRaises(ValueError):
                population.simulate(self.initial, self.times, signal)
        with self.assertRaises(ValueError):
            population.simulate(
                [0, 0, 101, 0, 0, 0], self.times, phenotype_signals("or")
            )
        with self.assertRaises(ValueError):
            population.simulate(
                self.initial, self.times, phenotype_signals("or"), rtol=0.1
            )

    def test_failed_integration_cannot_report_success(self):
        with patch.object(population, "solve_ivp") as solver:
            solver.return_value.success = False
            with self.assertRaises(RuntimeError):
                population.simulate(self.initial, self.times, phenotype_signals("or"))


class AdapterTests(unittest.TestCase):
    def test_published_cell_affinities_and_receptor_accounting(self):
        model = AdapterBinding.published_affinities()
        self.assertEqual(model.car_kd, 0.03)
        self.assertEqual(model.target_kd, 1)
        for a in (0, 1e-5, 0.01, 0.2, 100, 1e6):
            s = model.species(a)
            self.assertAlmostEqual(
                s["free_car"] + s["car_adapter"] + s["bridge"],
                model.car_total,
                places=13,
            )
            self.assertAlmostEqual(
                s["free_target"] + s["target_adapter"] + s["bridge"],
                model.target_total,
                places=13,
            )
            self.assertAlmostEqual(
                s["bridge"],
                s["free_car"] * s["free_target"] * a / (model.car_kd * model.target_kd),
                places=13,
            )
            self.assertTrue(all(np.isfinite(v) and v >= 0 for v in s.values()))
            self.assertLessEqual(model.bridge_fraction(a), 1)

    def test_hook_effect_and_analytic_peak(self):
        model = AdapterBinding()
        peak = model.peak_free_adapter
        self.assertEqual(model.bridge_fraction(0), 0)
        for factor in (0.001, 0.5, 2, 1000):
            self.assertLess(
                model.bridge_fraction(peak * factor), model.bridge_fraction(peak)
            )
        self.assertAlmostEqual(
            model.bridge_fraction(peak / 3), model.bridge_fraction(peak * 3), places=13
        )

    def test_pulse_superposition_and_half_life(self):
        times = np.array([0, 60, 120, 180])
        assert_allclose(
            free_adapter_profile(times, [(0, 1)], 60), [1, 0.5, 0.25, 0.125]
        )
        assert_allclose(
            free_adapter_profile(times, [(0, 1), (120, 2)], 60), [1, 0.5, 2.25, 1.125]
        )

    def test_invalid_inputs(self):
        for a in (-1, np.nan, True):
            with self.assertRaises(ValueError):
                AdapterBinding().species(a)
        with self.assertRaises(ValueError):
            AdapterBinding(car_kd=0)
        with self.assertRaises(ValueError):
            free_adapter_profile([0, 1], [(0, 1)], 0)


if __name__ == "__main__":
    unittest.main()
