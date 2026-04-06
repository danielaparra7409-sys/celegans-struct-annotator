from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_metadata_table(path: Path) -> pd.DataFrame:
    file_path = Path(path).expanduser().resolve()
    if not file_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {file_path}")

    suffix = file_path.suffix.lower()
    if suffix == ".xlsx":
        return pd.read_excel(file_path)
    if suffix == ".csv":
        return pd.read_csv(file_path)
    if suffix == ".tsv":
        return pd.read_csv(file_path, sep="\t")

    raise ValueError(
        "Unsupported metadata file format. Expected one of: .xlsx, .csv, .tsv"
    )
