import importlib.util
from pathlib import Path

import pandas as pd


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "projects/Genomic_Ancestry/scripts/assign_population.py"
)
spec = importlib.util.spec_from_file_location("assign_population", SCRIPT)
assign_population = importlib.util.module_from_spec(spec)
spec.loader.exec_module(assign_population)


def test_assign_populations_supports_fewer_than_five_reference_samples(tmp_path):
    reference = tmp_path / "reference.tsv"
    study = tmp_path / "study.tsv"
    output = tmp_path / "populations.tsv"
    reference.write_text(
        "ID\tPOPULATION\tPC1\tPC2\tPC3\tPC4\n"
        "ref1\tEUR\t0.0\t0.0\t0.0\t0.0\n"
        "ref2\tEUR\t0.1\t0.1\t0.1\t0.1\n"
    )
    study.write_text("SAMPLE\tPC1\tPC2\tPC3\tPC4\nsample\t0.02\t0.02\t0.02\t0.02\n")

    assign_population.assign_populations(str(reference), str(study), str(output))

    result = pd.read_csv(output, sep="\t")
    assert result.loc[0, "POPULATION"] == "EUR"
