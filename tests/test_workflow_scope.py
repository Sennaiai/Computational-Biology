"""Guard project isolation without executing completed scientific studies."""

from fnmatch import fnmatchcase
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def paths(workflow, event):
    text = (ROOT / ".github/workflows" / workflow).read_text()
    block = re.search(
        rf"^  {event}:\n(.*?)(?=^  [a-z_]+:|^permissions:)", text, re.M | re.S
    )
    if block is None:
        raise AssertionError(f"Missing {event} trigger")
    return re.findall(r"^      - '([^']+)'$", block.group(1), re.M)


class WorkflowScopeTests(unittest.TestCase):
    def test_cart_execution_is_isolated_from_completed_studies(self):
        cart_paths = (
            "projects/modeling/cart_cell/src/cart_dynamics/pkpd.py",
            "tests/modeling/cart_cell/test_models.py",
            ".github/workflows/cart-dynamics.yml",
        )
        for event in ("pull_request", "push"):
            active = paths("cart-dynamics.yml", event)
            for changed in cart_paths:
                self.assertTrue(any(fnmatchcase(changed, p) for p in active))
                for completed in ("quality-control.yml", "crispr-kinetics.yml"):
                    self.assertFalse(
                        any(fnmatchcase(changed, p) for p in paths(completed, event))
                    )
            for changed in (
                "README.md",
                "tests/test_workflow_scope.py",
                "projects/modeling/crispr_cas9/MODEL.md",
                "projects/Genomic_Ancestry/Snakefile",
            ):
                self.assertFalse(any(fnmatchcase(changed, p) for p in active))
        self.assertIn(
            "  workflow_dispatch:",
            (ROOT / ".github/workflows/cart-dynamics.yml").read_text(),
        )

    def test_crispr_changes_do_not_execute_bioinformatics(self):
        for event in ("pull_request", "push"):
            triggers = paths("quality-control.yml", event)
            self.assertTrue(
                any(
                    fnmatchcase("projects/Genomic_Ancestry/Snakefile", p)
                    for p in triggers
                )
            )
            for changed in (
                "projects/modeling/crispr_cas9/src/crispr_kinetics/kinetics.py",
                ".github/workflows/quality-control.yml",
                "tests/test_workflow_scope.py",
                "README.md",
            ):
                self.assertFalse(any(fnmatchcase(changed, p) for p in triggers))

    def test_crispr_execution_is_scoped_and_manual_reruns_remain_available(self):
        for event in ("pull_request", "push"):
            triggers = paths("crispr-kinetics.yml", event)
            self.assertTrue(
                any(
                    fnmatchcase(
                        "projects/modeling/crispr_cas9/src/crispr_kinetics/kinetics.py",
                        p,
                    )
                    for p in triggers
                )
            )
            self.assertFalse(
                any(
                    fnmatchcase("projects/Genomic_Ancestry/Snakefile", p)
                    for p in triggers
                )
            )
        for name in ("quality-control.yml", "crispr-kinetics.yml"):
            text = (ROOT / ".github/workflows" / name).read_text()
            self.assertIn("  workflow_dispatch:", text)
            self.assertNotIn("Sennaiai/Quality-Control-", text)

    def test_shared_checks_only_run_scope_guards(self):
        text = (ROOT / ".github/workflows/repository-checks.yml").read_text()
        self.assertIn(
            "python -m unittest discover -s tests -p test_workflow_scope.py -v", text
        )
        self.assertNotIn("pip install", text)
        self.assertNotIn("execute_notebook.py", text)


if __name__ == "__main__":
    unittest.main()
