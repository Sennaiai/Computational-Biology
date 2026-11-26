import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from crispr_kinetics.cell_study import bounded_mean, run


class CellularStudyTests(unittest.TestCase):
    def test_full_pipeline_exports_consistent_conditional_results(self):
        with tempfile.TemporaryDirectory() as directory:
            report = run(directory, founders=2)
            saved = json.loads((Path(directory) / "cellular_results.json").read_text())
            self.assertEqual(saved, report)
            self.assertEqual(len(report["scenarios"]), 7)
            for row in report["scenarios"].values():
                self.assertTrue(row["all_molecular_ledgers_exact"])
                for key in (
                    "target1_mutant_allele_fraction",
                    "target2_mutant_allele_fraction",
                    "double_biallelic_disruption_proxy",
                    "unresolved_break_cell_fraction",
                ):
                    expected = np.mean(
                        [
                            r["terminal_observables"][key]
                            for r in row["founder_outcomes"]
                        ]
                    )
                    self.assertEqual(row["final"][key]["mean"], expected)
                    self.assertLessEqual(row["final"][key]["MC_lower"], expected)
                    self.assertGreaterEqual(row["final"][key]["MC_upper"], expected)
                self.assertLessEqual(
                    row["final"]["mosaicism_proxy_lower_event"]["mean"],
                    row["final"]["mosaicism_proxy_upper_event"]["mean"],
                )
            control = report["scenarios"]["no_Cas9_control"]
            self.assertEqual(
                control["final"]["target1_mutant_allele_fraction"]["mean"], 0
            )
            self.assertEqual(
                control["final"]["target2_mutant_allele_fraction"]["mean"], 0
            )
            self.assertEqual(control["terminal_cell_counts"], [8, 8])
            self.assertIn(
                "<svg", (Path(directory) / "cellular_simulation.svg").read_text()
            )

    def test_sampling_bounds_do_not_claim_certainty_after_zero_events(self):
        result = bounded_mean([0] * 64)
        self.assertEqual(result["MC_lower"], 0)
        self.assertGreater(result["MC_upper"], 0)
        for values in ([], [float("nan")], [2], [-1]):
            with self.assertRaises(ValueError):
                bounded_mean(values)


if __name__ == "__main__":
    unittest.main()
