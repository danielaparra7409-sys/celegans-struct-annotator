from __future__ import annotations

import pandas as pd

from src.functional_labeler import (
    add_functional_labels,
    infer_family_label,
    merge_annotations_and_label,
)

def test_infer_family_label_acp() -> None:
    assert infer_family_label("Acyl carrier protein") == "ACP"


def test_infer_family_label_chaplin() -> None:
    assert infer_family_label("Chaplin domain protein") == "Chaplin"


def test_infer_family_label_unknown() -> None:
    assert infer_family_label("Completely unrelated description") == "unknown"


def test_add_functional_labels_uses_target_when_no_description_exists() -> None:
    frame = pd.DataFrame(
        [
            {
                "query": "q1",
                "target": "Acyl carrier protein",
                "fident": 0.49,
                "alnlen": 169,
                "evalue": 1.186e-18,
                "bits": 620,
            }
        ]
    )

    labeled = add_functional_labels(frame)

    assert "family_label" in labeled.columns
    assert "label_source" in labeled.columns
    assert labeled.iloc[0]["family_label"] == "ACP"
    assert labeled.iloc[0]["label_source"] == "target"


def test_add_functional_labels_prefers_explicit_text_column() -> None:
    frame = pd.DataFrame(
        [
            {
                "query": "q1",
                "target": "AF-A0A0000000-F1-model_v6",
                "hit_description": "Chaplin domain-containing protein",
            }
        ]
    )

    labeled = add_functional_labels(frame, text_column="hit_description")

    assert labeled.iloc[0]["family_label"] == "Chaplin"
    assert labeled.iloc[0]["label_source"] == "hit_description"


def test_add_functional_labels_returns_unknown_when_no_text_columns_exist() -> None:
    frame = pd.DataFrame([{"query": "q1", "fident": 0.5}])

    labeled = add_functional_labels(frame)

    assert labeled.iloc[0]["family_label"] == "unknown"
    assert labeled.iloc[0]["label_source"] == "no_text_column"
def test_merge_annotations_and_label_uses_lookup_table() -> None:
    frame = pd.DataFrame(
        [
            {
                "query": "q1",
                "target": "AF-TEST-1",
                "fident": 0.49,
                "alnlen": 169,
                "evalue": 1.186e-18,
                "bits": 620,
            }
        ]
    )

    annotations = pd.DataFrame(
        [
            {
                "target": "AF-TEST-1",
                "hit_description": "Acyl carrier protein",
            }
        ]
    )

    labeled = merge_annotations_and_label(frame, annotations)

    assert labeled.iloc[0]["family_label"] == "ACP"
    assert labeled.iloc[0]["label_source"] == "hit_description"
    assert labeled.iloc[0]["hit_description"] == "Acyl carrier protein"
