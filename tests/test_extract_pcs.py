import importlib.util
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "projects/Genomic_Ancestry/scripts/extract_pcs.py"
)
spec = importlib.util.spec_from_file_location("extract_pcs", SCRIPT)
extract_pcs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extract_pcs)


def test_reference_pcs_with_four_columns(tmp_path):
    values = tmp_path / "reference.V"
    populations = tmp_path / "populations.txt"
    values.write_text("0.1 0.2 0.3 0.4\n0.5 0.6 0.7 0.8\n")
    populations.write_text("sample1 POP1\nsample2 POP2\n")
    result = extract_pcs.load_reference_pcs(str(values), str(populations))
    assert result["ID"].tolist() == ["sample1", "sample2"]
    assert result["POPULATION"].tolist() == ["POP1", "POP2"]
    assert result["PC1"].tolist() == [0.1, 0.5]
    assert result["PC4"].tolist() == [0.4, 0.8]


def test_reference_pcs_reject_too_few_columns(tmp_path):
    values = tmp_path / "reference.V"
    populations = tmp_path / "populations.txt"
    values.write_text("0.1 0.2 0.3\n")
    populations.write_text("sample1 POP1\n")
    with pytest.raises(ValueError, match="Not enough columns"):
        extract_pcs.load_reference_pcs(str(values), str(populations))


def test_reference_pcs_aligns_sample_ids_instead_of_row_order(tmp_path):
    values = tmp_path / "reference.V"
    populations = tmp_path / "populations.txt"
    values.write_text(
        "sample2 0.5 0.6 0.7 0.8\n"
        "sample1 0.1 0.2 0.3 0.4\n"
    )
    populations.write_text("sample1 POP1\nsample2 POP2\n")

    result = extract_pcs.load_reference_pcs(str(values), str(populations))

    assert result["ID"].tolist() == ["sample1", "sample2"]
    assert result["PC1"].tolist() == [0.1, 0.5]


def test_reference_pcs_rejects_mismatched_row_counts(tmp_path):
    values = tmp_path / "reference.V"
    populations = tmp_path / "populations.txt"
    values.write_text("0.1 0.2 0.3 0.4\n0.5 0.6 0.7 0.8\n")
    populations.write_text("sample1 POP1\n")

    with pytest.raises(ValueError, match="positional alignment is unsafe"):
        extract_pcs.load_reference_pcs(str(values), str(populations))


def test_sample_id_prefers_input_and_falls_back_to_filename():
    assert extract_pcs.infer_sample_id_from_row_or_path(
        pd.Series({"SampleID": "actual"}), Path("fallback.Ancestry")
    ) == "actual"
    assert extract_pcs.infer_sample_id_from_row_or_path(
        pd.Series(dtype=object), Path("fallback.Ancestry")
    ) == "fallback"


def test_wide_study_pcs(tmp_path):
    ancestry = tmp_path / "sample.Ancestry"
    ancestry.write_text("SampleID PC1 PC2 PC3 PC4\nsample 0.1 0.2 0.3 0.4\n")
    result = extract_pcs.load_study_pcs([ancestry])
    assert result.iloc[0].to_dict() == {
        "SAMPLE": "sample",
        "PC1": 0.1,
        "PC2": 0.2,
        "PC3": 0.3,
        "PC4": 0.4,
    }


def test_wide_study_pcs_sort_pc_columns_numerically(tmp_path):
    ancestry = tmp_path / "sample.Ancestry"
    ancestry.write_text(
        "SampleID PC1 PC2 PC3 PC10 PC4\n"
        "sample 0.1 0.2 0.3 10.0 0.4\n"
    )

    result = extract_pcs.load_study_pcs([ancestry])

    assert result.loc[0, "PC4"] == 0.4


def test_study_pcs_rejects_missing_columns(tmp_path):
    ancestry = tmp_path / "sample.Ancestry"
    ancestry.write_text("SampleID PC1 PC2\nsample 0.1 0.2\n")
    with pytest.raises(ValueError, match="Could not find 4 PC columns"):
        extract_pcs.load_study_pcs([ancestry])
