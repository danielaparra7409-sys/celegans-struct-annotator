from pathlib import Path

import pytest

from src.foldseek_runner import FoldseekRunError, build_foldseek_command, run_foldseek


def test_build_foldseek_command_includes_top_hits(tmp_path: Path) -> None:
    command = build_foldseek_command(
        foldseek_bin=tmp_path / "foldseek",
        query_dir=tmp_path / "queries",
        target_database=tmp_path / "db",
        output_tsv=tmp_path / "results.tsv",
        tmp_dir=tmp_path / "tmp",
        top_hits=7,
    )

    assert command.args[-1] == "7"


def test_run_foldseek_creates_tmp_dir_and_checks_output(tmp_path: Path) -> None:
    foldseek_bin = tmp_path / "foldseek"
    foldseek_bin.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    foldseek_bin.chmod(0o755)

    query_dir = tmp_path / "queries"
    query_dir.mkdir()
    target_database = tmp_path / "db"
    target_database.write_text("db", encoding="utf-8")
    output_tsv = tmp_path / "out" / "results.tsv"
    tmp_dir = tmp_path / "tmp"

    with pytest.raises(FoldseekRunError, match="expected output file"):
        run_foldseek(
            foldseek_bin=foldseek_bin,
            query_dir=query_dir,
            target_database=target_database,
            output_tsv=output_tsv,
            tmp_dir=tmp_dir,
        )

    assert tmp_dir.exists()


def test_run_foldseek_rejects_missing_binary(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="Foldseek binary not found"):
        run_foldseek(
            foldseek_bin=tmp_path / "missing-foldseek",
            query_dir=tmp_path / "queries",
            target_database=tmp_path / "db",
            output_tsv=tmp_path / "results.tsv",
            tmp_dir=tmp_path / "tmp",
        )
