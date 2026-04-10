from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import ConfigError, load_config
from .data_loader import load_metadata_table
from .foldseek_runner import run_foldseek
from .result_parser import parse_foldseek_results
from .structure_manager import validate_cif_directory


def _build_parser(command_name: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=command_name)
    parser.add_argument("--config", required=True, type=Path)
    return parser


def validate_local_entrypoint(config_path: Path) -> None:
    config = load_config(config_path)
    if config.paths.metadata_file is None:
        raise ConfigError("Config is missing 'paths.metadata_file' for validate-local.")
    metadata = load_metadata_table(config.paths.metadata_file)
    cif_inventory = validate_cif_directory(config.paths.query_cif_dir)
    print(f"metadata_rows={len(metadata)}")
    print(cif_inventory.to_string(index=False))


def run_foldseek_entrypoint(config_path: Path) -> None:
    config = load_config(config_path)
    if config.paths.target_database is None:
        raise ConfigError(
            "Config is missing 'paths.target_database' or 'paths.test_database_dir' for run-foldseek."
        )
    output_path = run_foldseek(
        foldseek_bin=config.paths.foldseek_bin,
        query_dir=config.paths.query_cif_dir,
        target_database=config.paths.target_database,
        output_tsv=config.paths.foldseek_output_tsv,
        tmp_dir=config.paths.tmp_dir,
        top_hits=config.run.top_hits,
    )
    print(output_path)


def parse_results_entrypoint(config_path: Path) -> None:
    config = load_config(config_path)
    table = parse_foldseek_results(config.paths.foldseek_output_tsv)
    print(table.to_string(index=False))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.main")
    subparsers = parser.add_subparsers(dest="command")

    for command in ("validate-local", "run-foldseek", "parse-results"):
        command_parser = subparsers.add_parser(command)
        command_parser.add_argument("--config", required=True, type=Path)

    args = parser.parse_args(argv)
    if args.command == "validate-local":
        validate_local_entrypoint(args.config)
        return 0
    if args.command == "run-foldseek":
        run_foldseek_entrypoint(args.config)
        return 0
    if args.command == "parse-results":
        parse_results_entrypoint(args.config)
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
