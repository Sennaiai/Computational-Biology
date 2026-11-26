"""Execute the study notebook using the current Python interpreter."""

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
import nbformat
from nbclient import NotebookClient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    project = Path(__file__).resolve().parents[1]
    notebook = nbformat.read(project / "CRISPR_Cas9_Kinetics.ipynb", as_version=4)
    # Never reuse saved outputs or a kernel from a different environment.
    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    with tempfile.TemporaryDirectory(prefix="crispr-kernel-") as temp:
        root = Path(temp)
        spec = root / "crispr-verified"
        spec.mkdir()
        (spec / "kernel.json").write_text(
            json.dumps(
                {
                    "argv": [
                        sys.executable,
                        "-m",
                        "ipykernel_launcher",
                        "-f",
                        "{connection_file}",
                    ],
                    "display_name": "CRISPR verified interpreter",
                    "language": "python",
                }
            )
        )
        manager = KernelManager(
            kernel_name="crispr-verified",
            kernel_spec_manager=KernelSpecManager(kernel_dirs=[temp]),
        )
        client = NotebookClient(
            notebook,
            km=manager,
            timeout=120,
            allow_errors=False,
            resources={"metadata": {"path": str(project)}},
        )
        env = dict(os.environ)
        env["CRISPR_OUTPUT_DIR"] = str(output)
        env["MPLCONFIGDIR"] = str(root / "matplotlib")
        client.execute(env=env)
    notebook.metadata.kernelspec = {
        "display_name": "Python 3 (ipykernel)",
        "language": "python",
        "name": "python3",
    }
    nbformat.validate(notebook)
    nbformat.write(notebook, output / "CRISPR_Cas9_Kinetics.ipynb")
    print("Executed notebook and numerical evidence saved successfully.")


if __name__ == "__main__":
    main()
