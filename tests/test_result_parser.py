from pathlib import Path

import pytest

from src.result_parser import parse_foldseek_results


def test_parse_foldseek_results_uses_default_columns(tmp_path: Path) -> None:
    file_path = tmp_path / "results.tsv"
    file_path.write_text(
        "q1\tt1\t99.0\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n",
        encoding="utf-8",
    )

    frame = parse_foldseek_results(file_path)

    assert list(frame.columns[:4]) == ["query", "target", "fident", "alnlen"]
    assert frame.iloc[0]["query"] == "q1"


def test_parse_foldseek_results_rejects_empty_file(tmp_path: Path) -> None:
    file_path = tmp_path / "results.tsv"
    file_path.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="result file is empty"):
        parse_foldseek_results(file_path)
