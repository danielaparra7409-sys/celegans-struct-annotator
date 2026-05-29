from pathlib import Path
import pandas as pd

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
repo = Path.home() / "code" / "celegans-struct-annotator"

unique_table = base / "tables" / "proteome_unique_sequence_representatives.tsv"
cluster_dir = repo / "work" / "mmseqs" / "clusters"

cluster_files = {
    "90": cluster_dir / "proteome_cluster_90_cluster.tsv",
    "70": cluster_dir / "proteome_cluster_70_cluster.tsv",
    "50": cluster_dir / "proteome_cluster_50_cluster.tsv",
}

out = base / "tables" / "proteome_cluster_membership_strategy.tsv"

def clean_mmseqs_id(value):
    if pd.isna(value):
        return None
    v = str(value).strip()
    # Our FASTA headers start with:
    # sequence_hash|representative=...|organism=...
    # We only need the sequence_hash part.
    return v.split("|")[0]

unique = pd.read_csv(unique_table, sep="\t")

membership = unique.copy()

for level, path in cluster_files.items():
    cluster = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[f"cluster_{level}_representative_raw", f"cluster_{level}_member_raw"]
    )

    cluster[f"cluster_{level}_representative"] = cluster[f"cluster_{level}_representative_raw"].map(clean_mmseqs_id)
    cluster["sequence_hash"] = cluster[f"cluster_{level}_member_raw"].map(clean_mmseqs_id)

    cluster = cluster[
        [
            "sequence_hash",
            f"cluster_{level}_representative",
            f"cluster_{level}_representative_raw",
            f"cluster_{level}_member_raw",
        ]
    ]

    membership = membership.merge(
        cluster,
        on="sequence_hash",
        how="left"
    )

    membership[f"is_representative_{level}"] = (
        membership["sequence_hash"] == membership[f"cluster_{level}_representative"]
    )

cols_first = [
    "sequence_hash",
    "representative_canonical_id",
    "representative_record_id",
    "representative_organism",
    "sequence_length",
    "protein_count",
    "organism_count",
    "contains_nls_candidate",
    "nls_candidate_count",
    "cluster_90_representative",
    "is_representative_90",
    "cluster_70_representative",
    "is_representative_70",
    "cluster_50_representative",
    "is_representative_50",
]

cols_existing = [c for c in cols_first if c in membership.columns]
remaining = [c for c in membership.columns if c not in cols_existing]

membership = membership[cols_existing + remaining]

membership.to_csv(out, sep="\t", index=False)

print("Wrote:", out)
print("Rows:", len(membership))
print()

print("Missing cluster representatives:")
for level in ["90", "70", "50"]:
    print(f"{level}% missing:", membership[f"cluster_{level}_representative"].isna().sum())

print()
print("Representative counts:")
for level in ["90", "70", "50"]:
    print(f"{level}% representatives:", membership[f"is_representative_{level}"].sum())

print()
print("NLS-containing unique sequence groups:", membership["contains_nls_candidate"].sum())

print()
print("Preview NLS candidates:")
print(
    membership[membership["contains_nls_candidate"]][
        [
            "representative_canonical_id",
            "representative_organism",
            "sequence_length",
            "cluster_90_representative",
            "cluster_70_representative",
            "cluster_50_representative",
            "is_representative_50",
        ]
    ].head(20).to_string(index=False)
)
