from __future__ import annotations

from pathlib import Path

import pandas as pd

FOLDSEEK_DEFAULT_COLUMNS = [
    "query",
    "target",
    "fident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
]


def parse_foldseek_results(result_path: Path) -> pd.DataFrame:
    file_path = Path(result_path).expanduser().resolve()
    if not file_path.exists():
        raise FileNotFoundError(f"Foldseek result file not found: {file_path}")
    if file_path.stat().st_size == 0:
        raise ValueError(f"Foldseek result file is empty: {file_path}")

    frame = pd.read_csv(file_path, sep="\t", header=None)
    extra_columns = frame.shape[1] - len(FOLDSEEK_DEFAULT_COLUMNS)
    column_names = list(FOLDSEEK_DEFAULT_COLUMNS)
    if extra_columns > 0:
        column_names.extend(
            [f"extra_column_{index + 1}" for index in range(extra_columns)]
        )
    frame.columns = column_names[: frame.shape[1]]
    return frame
