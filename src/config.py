from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


class ConfigError(ValueError):
    """Raised when configuration content is invalid."""


@dataclass(frozen=True)
class PathsConfig:
    metadata_file: Path
    query_cif_dir: Path
    foldseek_bin: Path
    target_database: Path
    foldseek_output_tsv: Path
    tmp_dir: Path


@dataclass(frozen=True)
class RunConfig:
    top_hits: int = 10


@dataclass(frozen=True)
class AppConfig:
    project_name: str
    paths: PathsConfig
    run: RunConfig


def _resolve_path(raw_value: str, config_dir: Path) -> Path:
    path = Path(raw_value).expanduser()
    if not path.is_absolute():
        path = (config_dir / path).resolve()
    return path


def load_config(config_path: Path) -> AppConfig:
    config_path = Path(config_path).expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}

    try:
        raw_paths = payload["paths"]
    except KeyError as exc:
        raise ConfigError("Missing required 'paths' section in config.") from exc

    config_dir = config_path.parent
    paths = PathsConfig(
        metadata_file=_resolve_path(raw_paths["metadata_file"], config_dir),
        query_cif_dir=_resolve_path(raw_paths["query_cif_dir"], config_dir),
        foldseek_bin=_resolve_path(raw_paths["foldseek_bin"], config_dir),
        target_database=_resolve_path(raw_paths["target_database"], config_dir),
        foldseek_output_tsv=_resolve_path(raw_paths["foldseek_output_tsv"], config_dir),
        tmp_dir=_resolve_path(raw_paths["tmp_dir"], config_dir),
    )

    run_section = payload.get("run") or {}
    run = RunConfig(top_hits=int(run_section.get("top_hits", 10)))
    project_name = str(payload.get("project_name", "celegans-struct-annotator"))
    return AppConfig(project_name=project_name, paths=paths, run=run)
