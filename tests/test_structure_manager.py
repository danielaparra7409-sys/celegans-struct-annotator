from pathlib import Path

import pytest

from src.structure_manager import validate_cif_directory


def test_validate_cif_directory_builds_inventory(tmp_path: Path) -> None:
    query_dir = tmp_path / "queries"
    query_dir.mkdir()
    (query_dir / "valid.cif").write_text("data_example\n_entry.id 1\n", encoding="utf-8")
    (query_dir / "empty.cif").write_text("", encoding="utf-8")
    (query_dir / "bad.cif").write_text("not_a_cif\n", encoding="utf-8")

    frame = validate_cif_directory(query_dir)

    statuses = dict(zip(frame["file_name"], frame["status"]))
    assert statuses["valid.cif"] == "ok"
    assert statuses["empty.cif"] == "empty"
    assert statuses["bad.cif"] == "invalid"


def test_validate_cif_directory_requires_cif_files(tmp_path: Path) -> None:
    query_dir = tmp_path / "queries"
    query_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="No .cif files found"):
        validate_cif_directory(query_dir)
