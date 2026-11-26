import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

import crispr_kinetics
from crispr_kinetics.study import run


class StudyIntegrationTests(unittest.TestCase):
    def test_full_experiment_exports_consistent_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            os.environ.setdefault("MPLCONFIGDIR", str(Path(directory) / "mpl"))
            report = run(directory)
            saved = json.loads((Path(directory) / "kinetics.json").read_text())
            self.assertEqual(saved, report)
            self.assertEqual(report["biological_calibration"], "not_performed")
            self.assertEqual(
                report["parameter_source"]["provenance"], "synthetic_demonstration"
            )
            reference = np.asarray(report["reference_probabilities"])
            for ensemble in report["stochastic_validation"]:
                counts = np.asarray(ensemble["counts"])
                np.testing.assert_array_equal(counts.sum(axis=1), ensemble["paths"])
                self.assertAlmostEqual(
                    float(abs(counts / ensemble["paths"] - reference).max()),
                    ensemble["max_absolute_probability_gap"],
                )
            package = Path(crispr_kinetics.__file__).parent
            for name, digest in report["provenance"]["source_sha256"].items():
                self.assertEqual(
                    hashlib.sha256((package / name).read_bytes()).hexdigest(), digest
                )
            figure = (Path(directory) / "kinetics.svg").read_text()
            self.assertIn("Synthetic kinetic baseline", figure)
            self.assertIn("<svg", figure)


if __name__ == "__main__":
    unittest.main()
