from __future__ import annotations

from pathlib import Path

import pandas as pd


def _inspect_cif_file(file_path: Path) -> tuple[str, str]:
    if file_path.stat().st_size == 0:
        return "empty", "File exists but has zero size."

    with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
        first_chunk = handle.read(512)

    if "data_" not in first_chunk:
        return "invalid", "Missing expected CIF header token 'data_'."

    return "ok", ""


def validate_cif_directory(query_dir: Path) -> pd.DataFrame:
    directory = Path(query_dir).expanduser().resolve()
    if not directory.exists():
        raise FileNotFoundError(f"Query CIF directory not found: {directory}")
    if not directory.is_dir():
        raise NotADirectoryError(f"Query CIF path is not a directory: {directory}")

    records: list[dict[str, object]] = []
    for cif_file in sorted(directory.glob("*.cif")):
        status, message = _inspect_cif_file(cif_file)
        records.append(
            {
                "file_name": cif_file.name,
                "file_path": str(cif_file),
                "status": status,
                "message": message,
                "size_bytes": cif_file.stat().st_size,
            }
        )

    if not records:
        raise FileNotFoundError(f"No .cif files found in query directory: {directory}")

    return pd.DataFrame.from_records(records)
