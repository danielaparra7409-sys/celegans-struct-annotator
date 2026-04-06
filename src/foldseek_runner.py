from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess


class FoldseekRunError(RuntimeError):
    """Raised when Foldseek execution or output validation fails."""


@dataclass(frozen=True)
class FoldseekCommand:
    args: list[str]
    output_tsv: Path


def validate_foldseek_inputs(
    foldseek_bin: Path,
    query_dir: Path,
    target_database: Path,
    tmp_dir: Path,
) -> None:
    if not foldseek_bin.exists():
        raise FileNotFoundError(f"Foldseek binary not found: {foldseek_bin}")
    if not foldseek_bin.is_file():
        raise FileNotFoundError(f"Foldseek binary path is not a file: {foldseek_bin}")
    if not query_dir.exists():
        raise FileNotFoundError(f"Foldseek query directory not found: {query_dir}")
    if not query_dir.is_dir():
        raise NotADirectoryError(f"Foldseek query path is not a directory: {query_dir}")
    if not target_database.exists():
        raise FileNotFoundError(f"Foldseek target database not found: {target_database}")
    tmp_dir.mkdir(parents=True, exist_ok=True)


def build_foldseek_command(
    foldseek_bin: Path,
    query_dir: Path,
    target_database: Path,
    output_tsv: Path,
    tmp_dir: Path,
    top_hits: int,
) -> FoldseekCommand:
    args = [
        str(foldseek_bin),
        "easy-search",
        str(query_dir),
        str(target_database),
        str(output_tsv),
        str(tmp_dir),
        "--format-mode",
        "4",
        "--max-seqs",
        str(top_hits),
    ]
    return FoldseekCommand(args=args, output_tsv=output_tsv)


def _validate_output_file(output_tsv: Path) -> None:
    if not output_tsv.exists():
        raise FoldseekRunError(
            f"Foldseek finished without creating the expected output file: {output_tsv}"
        )
    if output_tsv.stat().st_size == 0:
        raise FoldseekRunError(f"Foldseek created an empty output file: {output_tsv}")


def run_foldseek(
    foldseek_bin: Path,
    query_dir: Path,
    target_database: Path,
    output_tsv: Path,
    tmp_dir: Path,
    top_hits: int = 10,
) -> Path:
    foldseek_bin = Path(foldseek_bin).expanduser().resolve()
    query_dir = Path(query_dir).expanduser().resolve()
    target_database = Path(target_database).expanduser().resolve()
    output_tsv = Path(output_tsv).expanduser().resolve()
    tmp_dir = Path(tmp_dir).expanduser().resolve()

    validate_foldseek_inputs(
        foldseek_bin=foldseek_bin,
        query_dir=query_dir,
        target_database=target_database,
        tmp_dir=tmp_dir,
    )

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = build_foldseek_command(
        foldseek_bin=foldseek_bin,
        query_dir=query_dir,
        target_database=target_database,
        output_tsv=output_tsv,
        tmp_dir=tmp_dir,
        top_hits=top_hits,
    )

    completed = subprocess.run(
        command.args,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise FoldseekRunError(
            "Foldseek execution failed with a non-zero exit code.\n"
            f"Command: {' '.join(command.args)}\n"
            f"stdout: {completed.stdout.strip()}\n"
            f"stderr: {completed.stderr.strip()}"
        )

    _validate_output_file(command.output_tsv)
    return command.output_tsv
