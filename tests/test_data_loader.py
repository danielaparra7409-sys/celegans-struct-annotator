from pathlib import Path

import pandas as pd
import pytest

from src.data_loader import load_metadata_table


def test_load_metadata_csv(tmp_path: Path) -> None:
    file_path = tmp_path / "metadata.csv"
    file_path.write_text("id,value\nA,1\n", encoding="utf-8")

    frame = load_metadata_table(file_path)

    assert list(frame.columns) == ["id", "value"]
    assert frame.iloc[0]["id"] == "A"


def test_load_metadata_tsv(tmp_path: Path) -> None:
    file_path = tmp_path / "metadata.tsv"
    file_path.write_text("id\tvalue\nA\t1\n", encoding="utf-8")

    frame = load_metadata_table(file_path)

    assert frame.iloc[0]["value"] == 1


def test_load_metadata_xlsx(tmp_path: Path) -> None:
    file_path = tmp_path / "metadata.xlsx"
    pd.DataFrame([{"id": "A", "value": 1}]).to_excel(file_path, index=False)

    frame = load_metadata_table(file_path)

    assert frame.iloc[0]["id"] == "A"


def test_load_metadata_rejects_unknown_extension(tmp_path: Path) -> None:
    file_path = tmp_path / "metadata.txt"
    file_path.write_text("id,value\nA,1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported metadata file format"):
        load_metadata_table(file_path)
