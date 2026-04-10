from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


class ConfigError(ValueError):
    """Raised when configuration content is invalid."""


@dataclass(frozen=True)
class PathsConfig:
    metadata_file: Path | None
    query_cif_dir: Path
    foldseek_bin: Path
    target_database: Path | None
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


def _get_required_path(
    raw_paths: dict[str, str], config_dir: Path, *keys: str
) -> Path:
    for key in keys:
        value = raw_paths.get(key)
        if value:
            return _resolve_path(value, config_dir)
    joined_keys = ", ".join(keys)
    raise ConfigError(f"Missing required path config key. Expected one of: {joined_keys}")


def _get_optional_path(
    raw_paths: dict[str, str], config_dir: Path, *keys: str
) -> Path | None:
    for key in keys:
        value = raw_paths.get(key)
        if value:
            return _resolve_path(value, config_dir)
    return None


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
        metadata_file=_get_optional_path(raw_paths, config_dir, "metadata_file"),
        query_cif_dir=_get_required_path(raw_paths, config_dir, "query_cif_dir", "local_cif_dir"),
        foldseek_bin=_get_required_path(raw_paths, config_dir, "foldseek_bin"),
        target_database=_get_optional_path(raw_paths, config_dir, "target_database", "test_database_dir"),
        foldseek_output_tsv=_get_required_path(
            raw_paths, config_dir, "foldseek_output_tsv", "local_foldseek_results"
        ),
        tmp_dir=_get_required_path(raw_paths, config_dir, "tmp_dir"),
    )

    run_section = payload.get("run") or {}
    run = RunConfig(top_hits=int(run_section.get("top_hits", 10)))
    project_name = str(payload.get("project_name", "celegans-struct-annotator"))
    return AppConfig(project_name=project_name, paths=paths, run=run)
