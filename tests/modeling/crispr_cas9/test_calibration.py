import unittest

import numpy as np

from crispr_kinetics.calibration import (
    assay_data,
    fit,
    negative_log_likelihood,
    predict,
    run_calibration,
)


class CalibrationTests(unittest.TestCase):
    def test_pinned_data_schema_and_true_column_order(self):
        times, counts, error, source = assay_data()
        self.assertEqual(source["raw_columns"]["y2"], "indels")
        np.testing.assert_array_equal(counts[0], [3095, 0, 62, 9])
        self.assertEqual(int((times == 72).sum()), 2)
        np.testing.assert_allclose(error.sum(axis=1), 1, atol=1e-14)

    def test_assay_model_matches_analytic_irreversible_cascade(self):
        t = np.array([0, 1, 3, 8])
        r = 0.2
        p = predict([r, 0, r, 0.25, 0], t, np.eye(4))
        intact = np.exp(-r * t)
        broken = r * t * intact
        mutant = 1 - intact - broken
        np.testing.assert_allclose(p[:, 0], intact, atol=1e-12)
        np.testing.assert_allclose(p[:, 1], mutant, atol=1e-12)
        np.testing.assert_allclose(p[:, 2], 0.75 * broken, atol=1e-12)
        np.testing.assert_allclose(p[:, 3], 0.25 * broken, atol=1e-12)

    def test_fit_recovers_synthetic_assay_parameters(self):
        t = np.array([0, 1, 3, 7, 15, 30, 60, 100])
        truth = [0.04, 0.08, 0.03, 0.2, 0.025]
        counts = predict(truth, t, np.eye(4)) * 1000000
        result = fit(t, counts, np.eye(4), starts=4)
        np.testing.assert_allclose(result["parameters"], truth, rtol=0.003, atol=1e-5)

    def test_observation_mapping_and_initial_conditions(self):
        initial = predict([0, 0, 0, 0.3, 0.2], [0], np.eye(4))
        np.testing.assert_allclose(initial, [[0.8, 0, 0.14, 0.06]])
        with self.assertRaises(ValueError):
            predict([0, 0, 0, 2, 0], [0], np.eye(4))
        self.assertTrue(
            np.isinf(
                negative_log_likelihood([0, 0, 0, 0, 0], [0], [[0, 1, 0, 0]], np.eye(4))
            )
        )

    def test_whole_timepoint_holdout_and_profile(self):
        report = run_calibration()
        self.assertEqual(report["training_rows"], 12)
        self.assertEqual(report["held_out_rows"], 2)
        self.assertNotIn(72, report["training_times_h"])
        self.assertEqual(report["held_out_times_h"], [72.0])
        self.assertGreaterEqual(report["full"]["successful_starts"], 1)
        self.assertLessEqual(
            report["full"]["training_nll"],
            report["no_precise_repair"]["training_nll"] + 1e-5,
        )
        self.assertTrue(
            all(
                np.isfinite(row["training_nll"])
                for row in report["precise_repair_profile"]
            )
        )


if __name__ == "__main__":
    unittest.main()
