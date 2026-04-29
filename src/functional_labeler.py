from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd


@dataclass(frozen=True)
class LabelRule:
    label: str
    patterns: tuple[str, ...]


DEFAULT_RULES: tuple[LabelRule, ...] = (
    LabelRule("Chaplin", ("chaplin",)),
    LabelRule("PE/PPE", ("pe/ppe", "pe35", "ppe", "pe family")),
    LabelRule("PKD", ("pkd", "polycystic kidney disease")),
    LabelRule("ACP", ("acyl carrier protein", "acp")),
    LabelRule("NlpC/P60", ("nlpc/p60", "nlpc", "p60")),
    LabelRule("GH16", ("gh16", "glycosyl hydrolase family 16")),
    LabelRule("PIN/VapC", ("pin", "vapc")),
)


def _normalize_text(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def infer_family_label(
    text: object,
    rules: Iterable[LabelRule] = DEFAULT_RULES,
) -> str:
    normalized = _normalize_text(text)
    if not normalized:
        return "unknown"

    for rule in rules:
        for pattern in rule.patterns:
            if pattern in normalized:
                return rule.label

    return "unknown"


def add_functional_labels(
    frame: pd.DataFrame,
    text_column: str | None = None,
    rules: Iterable[LabelRule] = DEFAULT_RULES,
) -> pd.DataFrame:
    """
    Add a family_label column to a parsed Foldseek dataframe.

    Notes
    -----
    This is intentionally conservative.
    It does not invent annotations from AFDB IDs alone.
    It only labels rows when a known textual pattern is present.
    """

    labeled = frame.copy()

    candidate_columns: list[str] = []
    if text_column is not None:
        candidate_columns.append(text_column)

    candidate_columns.extend(
        [
            "hit_description",
            "target_description",
            "annotation",
            "target",
        ]
    )

    chosen_column = next((col for col in candidate_columns if col in labeled.columns), None)

    if chosen_column is None:
        labeled["family_label"] = "unknown"
        labeled["label_source"] = "no_text_column"
        return labeled

    labeled["family_label"] = labeled[chosen_column].apply(
        lambda value: infer_family_label(value, rules)
    )
    labeled["label_source"] = chosen_column
    return labeled
