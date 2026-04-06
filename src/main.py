from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config
from .data_loader import load_metadata_table
from .foldseek_runner import run_foldseek
from .result_parser import parse_foldseek_results
from .structure_manager import validate_cif_directory


def _build_parser(command_name: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=command_name)
    parser.add_argument("--config", required=True, type=Path)
    return parser


def validate_local_entrypoint() -> None:
    args = _build_parser("validate-local").parse_args()
    config = load_config(args.config)
    metadata = load_metadata_table(config.paths.metadata_file)
    cif_inventory = validate_cif_directory(config.paths.query_cif_dir)
    print(f"metadata_rows={len(metadata)}")
    print(cif_inventory.to_string(index=False))


def run_foldseek_entrypoint() -> None:
    args = _build_parser("run-foldseek").parse_args()
    config = load_config(args.config)
    output_path = run_foldseek(
        foldseek_bin=config.paths.foldseek_bin,
        query_dir=config.paths.query_cif_dir,
        target_database=config.paths.target_database,
        output_tsv=config.paths.foldseek_output_tsv,
        tmp_dir=config.paths.tmp_dir,
        top_hits=config.run.top_hits,
    )
    print(output_path)


def parse_results_entrypoint() -> None:
    args = _build_parser("parse-results").parse_args()
    config = load_config(args.config)
    table = parse_foldseek_results(config.paths.foldseek_output_tsv)
    print(table.to_string(index=False))
