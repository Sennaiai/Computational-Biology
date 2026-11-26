"""Ensure public numerical evidence is current and internally consistent."""

import json
from pathlib import Path
import unittest

from cart_dynamics import pkpd
from cart_dynamics.study import source_digest


PROJECT = Path(__file__).resolve().parents[3] / "projects/modeling/cart_cell"


class EvidenceTests(unittest.TestCase):
    def test_saved_evidence_matches_current_implementation(self):
        for name in ("cart_results.json", "pkpd_results.json"):
            report = json.loads((PROJECT / "results" / name).read_text())
            self.assertEqual(report["source_sha256"], source_digest())
            self.assertGreater(report["end_to_end_seconds_including_figures"], 0)

    def test_every_paired_scenario_retains_its_parameters_and_accounting(self):
        report = json.loads((PROJECT / "results/cart_results.json").read_text())
        rows = report["uncertainty"]["runs"]
        self.assertEqual(len(rows), 64)
        for row in rows:
            self.assertEqual(
                set(row["parameters"]), set(report["configuration"]["parameters"])
            )
            for outcome in row["results"].values():
                self.assertLess(outcome["tumor_ledger_max_error_million_cells"], 1e-7)
                self.assertGreaterEqual(outcome["minimum_integrated_state"], -1e-8)

    def test_pkpd_sensitivity_and_pareto_reports_are_consistent(self):
        report = json.loads((PROJECT / "results/pkpd_results.json").read_text())
        rows = report["sensitivity"]["runs"]
        self.assertEqual(len(rows), 16)
        for row in rows:
            self.assertEqual(
                row["peak_proxy_frontier"], pkpd.nondominated(row["outcomes"])
            )
            for outcome in row["outcomes"].values():
                self.assertAlmostEqual(
                    outcome["cumulative_concentration_equivalent_input_nM"], 0.7
                )
                self.assertLess(outcome["pk_ledger_error_nM"], 1e-8)
                self.assertLess(outcome["tumor_ledger_error_million_cells"], 1e-7)


if __name__ == "__main__":
    unittest.main()
