from pathlib import Path

import pytest

from src.result_parser import (
    FOLDSEEK_COMPACT_COLUMNS,
    FOLDSEEK_DEFAULT_COLUMNS,
    parse_foldseek_results,
)


def test_parse_foldseek_results_uses_default_columns_for_12_col_output(
    tmp_path: Path,
) -> None:
    file_path = tmp_path / "results.tsv"
    file_path.write_text(
        "q1\tt1\t99.0\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n",
        encoding="utf-8",
    )

    frame = parse_foldseek_results(file_path)

    assert list(frame.columns) == FOLDSEEK_DEFAULT_COLUMNS
    assert frame.iloc[0]["query"] == "q1"
    assert frame.iloc[0]["fident"] == pytest.approx(99.0)
    assert frame.iloc[0]["alnlen"] == 10
    assert frame.iloc[0]["evalue"] == pytest.approx(1e-5)
    assert frame.iloc[0]["bits"] == 50


def test_parse_foldseek_results_uses_compact_columns_for_6_col_output(
    tmp_path: Path,
) -> None:
    file_path = tmp_path / "results.m8"
    file_path.write_text(
        "q1\tt1\t0.491\t169\t1.186E-18\t620\n",
        encoding="utf-8",
    )

    frame = parse_foldseek_results(file_path)

    assert list(frame.columns) == FOLDSEEK_COMPACT_COLUMNS
    assert frame.iloc[0]["query"] == "q1"
    assert frame.iloc[0]["fident"] == pytest.approx(0.491)
    assert frame.iloc[0]["alnlen"] == 169
    assert frame.iloc[0]["evalue"] == pytest.approx(1.186e-18)
    assert frame.iloc[0]["bits"] == 620


def test_parse_foldseek_results_rejects_empty_file(tmp_path: Path) -> None:
    file_path = tmp_path / "results.tsv"
    file_path.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="result file is empty"):
        parse_foldseek_results(file_path)
