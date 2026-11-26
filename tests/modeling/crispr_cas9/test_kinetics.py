import unittest

import numpy as np

from crispr_kinetics.kinetics import (
    Rates,
    generator,
    solve,
    simulate,
    readouts,
    first_cleavage_time,
    diagnostics,
)


class KineticsTests(unittest.TestCase):
    def test_reversible_first_passage_matches_closed_form(self):
        rates = Rates()
        a = rates.association_per_nM_s * rates.free_rnp_nM
        b, c = rates.dissociation_per_s, rates.r_loop_per_s
        d, e = rates.unwind_per_s, rates.cleavage_per_s
        expected = 1 / a + (1 + b / a) / c + (1 + d * (1 + b / a) / c) / e
        self.assertAlmostEqual(first_cleavage_time(rates)["mean_seconds"], expected)

    def test_mixed_initial_distribution_and_nonzero_start_grid(self):
        initial = np.array([0.5, 0, 0, 0.25, 0, 0.25, 0])
        rates = Rates(free_rnp_nM=0, repair_per_s=0)
        times = [2, 17, 99]
        reference = np.tile(initial, (len(times), 1))
        for method in ("expm", "ode"):
            np.testing.assert_allclose(
                solve(rates, times, initial=initial, method=method), reference
            )
        sampled = simulate(rates, times, initial=initial, paths=4000, seed=39) / 4000
        self.assertLess(float(abs(sampled - reference).max()), 0.035)

    def test_zero_time_and_invalid_probability_readouts(self):
        for method in ("expm", "ode"):
            np.testing.assert_array_equal(
                solve(Rates(), [0], method=method), np.eye(7)[:1]
            )
        for values in ([], np.zeros((0, 7)), [["x"] * 7], np.eye(7, dtype=complex)):
            with self.assertRaises(ArithmeticError):
                diagnostics(values)

    def test_generator_conserves_probability_and_has_terminal_repair_states(self):
        q = generator()
        np.testing.assert_allclose(q.sum(axis=0), 0, rtol=0, atol=1e-15)
        off = q.copy()
        np.fill_diagonal(off, 0)
        self.assertTrue(np.all(off >= 0))
        np.testing.assert_array_equal(q[:, 4:], 0)

    def test_association_unit_conversion(self):
        q = generator(Rates(association_per_nM_s=0.02, free_rnp_nM=5))
        self.assertAlmostEqual(q[1, 0], 0.1)

    def test_binding_only_matches_closed_form(self):
        rates = Rates(dissociation_per_s=0, r_loop_per_s=0)
        times = np.array([0, 1, 3, 15, 60.0])
        expected = np.zeros((len(times), 7))
        expected[:, 0] = np.exp(-0.1 * times)
        expected[:, 1] = 1 - expected[:, 0]
        for method in ("ode", "expm"):
            np.testing.assert_allclose(
                solve(rates, times, method=method), expected, rtol=0, atol=1e-9
            )

    def test_equal_rate_cascade_matches_erlang_cdf_and_mean(self):
        r = 0.2
        rates = Rates(
            association_per_nM_s=r,
            free_rnp_nM=1,
            dissociation_per_s=0,
            r_loop_per_s=r,
            unwind_per_s=0,
            cleavage_per_s=r,
            repair_per_s=0,
        )
        times = np.linspace(0, 40, 50)
        z = r * times
        expected = np.zeros((len(times), 7))
        expected[:, 0] = np.exp(-z)
        expected[:, 1] = z * np.exp(-z)
        expected[:, 2] = z**2 * np.exp(-z) / 2
        expected[:, 3] = 1 - expected[:, :3].sum(axis=1)
        for method in ("expm", "ode"):
            np.testing.assert_allclose(
                solve(rates, times, method=method), expected, rtol=0, atol=1e-9
            )
        self.assertAlmostEqual(first_cleavage_time(rates)["mean_seconds"], 3 / r)

    def test_repair_layer_matches_independent_branching_formula(self):
        rates = Rates(repair_per_s=0.04)
        initial = np.eye(7)[3]
        times = np.array([0, 2, 10, 80.0])
        expected = np.zeros((len(times), 7))
        expected[:, 3] = np.exp(-0.04 * times)
        for i, branch in enumerate((0.4, 0.5, 0.1), start=4):
            expected[:, i] = branch * (1 - expected[:, 3])
        p = solve(rates, times, initial=initial)
        np.testing.assert_allclose(p, expected, atol=1e-12)
        np.testing.assert_allclose(readouts(p)["cleaved_ever"], 1)
        self.assertTrue(
            np.all(np.diff(readouts(p)["currently_cleaved_unrepaired"]) < 0)
        )

    def test_zero_concentration_or_blocked_path_cannot_cleave(self):
        for rates in (
            Rates(free_rnp_nM=0),
            Rates(r_loop_per_s=0),
            Rates(cleavage_per_s=0),
        ):
            p = solve(rates, [0, 1, 300])
            np.testing.assert_array_equal(readouts(p)["cleaved_ever"], 0)
            self.assertIsNone(first_cleavage_time(rates))

    def test_ode_and_expm_agree_without_probability_clipping(self):
        times = np.linspace(0, 600, 61)
        for rates in (Rates(), Rates(free_rnp_nM=1), Rates(repair_per_s=0)):
            a = solve(rates, times)
            b = solve(rates, times, method="ode")
            np.testing.assert_allclose(a, b, rtol=0, atol=2e-9)
            self.assertTrue(np.all(np.diff(readouts(a)["cleaved_ever"]) >= -1e-12))

    def test_repair_parameters_do_not_change_first_cleavage_statistics(self):
        a, b = (
            Rates(),
            Rates(repair_per_s=0, precise_fraction=0, indel_fraction=1, hdr_fraction=0),
        )
        self.assertEqual(first_cleavage_time(a), first_cleavage_time(b))
        for t in ([0, 10, 100],):
            np.testing.assert_allclose(
                readouts(solve(a, t))["cleaved_ever"],
                readouts(solve(b, t))["cleaved_ever"],
                atol=1e-13,
            )

    def test_ssa_agrees_with_probability_reference_under_union_bound(self):
        times = np.linspace(0, 180, 19)
        paths = 2500
        counts = simulate(Rates(), times, paths=paths, seed=2026)
        reference = solve(Rates(), times)
        epsilon = np.sqrt(np.log(2 * reference.size / 1e-6) / (2 * paths))
        self.assertLess(float(abs(counts / paths - reference).max()), epsilon)
        np.testing.assert_array_equal(counts.sum(axis=1), paths)

    def test_rng_is_local_and_reproducible(self):
        np.random.seed(42)
        before = np.random.get_state()
        a = simulate(Rates(), [0, 20, 100], paths=20, seed=17)
        b = simulate(Rates(), [0, 20, 100], paths=20, seed=17)
        np.testing.assert_array_equal(a, b)
        after = np.random.get_state()
        np.testing.assert_array_equal(before[1], after[1])
        self.assertEqual(before[2:], after[2:])

    def test_ssa_absorbing_states_and_event_budget(self):
        counts = simulate(Rates(free_rnp_nM=0), [0, 5, 100], paths=12)
        np.testing.assert_array_equal(counts[:, 0], 12)
        with self.assertRaises(RuntimeError):
            simulate(Rates(cleavage_per_s=0), [1e5], paths=2, event_budget=1)

    def test_bad_inputs_are_rejected(self):
        for kwargs in (
            {"free_rnp_nM": True},
            {"association_per_nM_s": -1},
            {"cleavage_per_s": float("nan")},
            {"repair_per_s": 1j},
            {"precise_fraction": 0.2},
            {"hdr_fraction": 2},
        ):
            with self.assertRaises(ValueError):
                Rates(**kwargs)
        for times in (
            [],
            [-1, 2],
            [2, 1],
            [1, 1],
            [0, float("inf")],
            [[1]],
            [True],
            ["1"],
        ):
            with self.assertRaises(ValueError):
                solve(Rates(), times)
        with self.assertRaises(ValueError):
            solve(Rates(), [0, 1], initial=np.ones(7))
        for kwargs in (
            {"paths": 0},
            {"paths": True},
            {"seed": -1},
            {"event_budget": 0},
        ):
            with self.assertRaises(ValueError):
                simulate(Rates(), [0, 1], **kwargs)


if __name__ == "__main__":
    unittest.main()
