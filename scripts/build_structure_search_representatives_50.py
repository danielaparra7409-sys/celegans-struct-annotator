from pathlib import Path
import pandas as pd
import textwrap

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
repo = Path.home() / "code" / "celegans-struct-annotator"

membership_file = base / "tables" / "proteome_cluster_membership_strategy.tsv"
cluster_file = repo / "work" / "mmseqs" / "clusters" / "proteome_cluster_50_cluster.tsv"

out_table = base / "tables" / "proteome_structure_search_representatives_50.tsv"
out_fasta = base / "tables" / "proteome_structure_search_representatives_50.fasta"

df = pd.read_csv(membership_file, sep="\t")

def clean_mmseqs_id(value):
    if pd.isna(value):
        return None
    return str(value).strip().split("|")[0]

cluster = pd.read_csv(
    cluster_file,
    sep="\t",
    header=None,
    names=["cluster_50_representative_raw", "cluster_50_member_raw"]
)

cluster["cluster_50_representative"] = cluster["cluster_50_representative_raw"].map(clean_mmseqs_id)
cluster["sequence_hash"] = cluster["cluster_50_member_raw"].map(clean_mmseqs_id)

cluster_sizes = (
    cluster.groupby("cluster_50_representative")
    .agg(cluster_50_member_count=("sequence_hash", "count"))
    .reset_index()
)

cluster_nls = (
    df[["sequence_hash", "contains_nls_candidate"]]
    .merge(cluster[["sequence_hash", "cluster_50_representative"]], on="sequence_hash", how="left")
    .groupby("cluster_50_representative")
    .agg(cluster_50_contains_nls_candidate=("contains_nls_candidate", "any"),
         cluster_50_nls_sequence_count=("contains_nls_candidate", "sum"))
    .reset_index()
)

reps = df[df["is_representative_50"]].copy()

reps = reps.merge(cluster_sizes, on="cluster_50_representative", how="left")
reps = reps.merge(cluster_nls, on="cluster_50_representative", how="left")

reps["structure_search_level"] = "cluster_50_representative"
reps["recommended_first_strategy"] = "sequence_based_structure_search"
reps["structure_search_status"] = "pending"
reps["matched_uniprot_id"] = ""
reps["matched_structure_id"] = ""
reps["structure_file_path"] = ""
reps["propagate_to_cluster_members"] = "yes"

cols_first = [
    "sequence_hash",
    "representative_canonical_id",
    "representative_record_id",
    "representative_organism",
    "sequence_length",
    "cluster_50_member_count",
    "cluster_50_contains_nls_candidate",
    "cluster_50_nls_sequence_count",
    "protein_count",
    "organism_count",
    "contains_nls_candidate",
    "nls_candidate_count",
    "cluster_50_representative",
    "is_representative_50",
    "structure_search_level",
    "recommended_first_strategy",
    "structure_search_status",
    "matched_uniprot_id",
    "matched_structure_id",
    "structure_file_path",
    "propagate_to_cluster_members",
    "sequence",
]

cols_existing = [c for c in cols_first if c in reps.columns]
remaining = [c for c in reps.columns if c not in cols_existing]
reps = reps[cols_existing + remaining]

reps.to_csv(out_table, sep="\t", index=False)

with out_fasta.open("w", encoding="utf-8") as f:
    for _, row in reps.iterrows():
        header = (
            f">{row['sequence_hash']}"
            f"|representative={row['representative_canonical_id']}"
            f"|organism={row['representative_organism']}"
            f"|cluster_50_members={row['cluster_50_member_count']}"
            f"|cluster_contains_nls={row['cluster_50_contains_nls_candidate']}"
        )
        seq = str(row["sequence"]).strip()
        f.write(header + "\n")
        f.write("\n".join(textwrap.wrap(seq, width=80)) + "\n")

print("Wrote:", out_table)
print("Wrote:", out_fasta)
print("Representatives 50:", len(reps))
print("Unique sequences represented by 50% clusters:", int(reps["cluster_50_member_count"].sum()))
print("Representatives whose cluster contains NLS candidates:", int(reps["cluster_50_contains_nls_candidate"].sum()))
print("NLS sequence memberships across 50% clusters:", int(reps["cluster_50_nls_sequence_count"].sum()))
print()
print("Top representatives by cluster size:")
print(
    reps[
        [
            "representative_canonical_id",
            "representative_organism",
            "sequence_length",
            "cluster_50_member_count",
            "cluster_50_contains_nls_candidate",
        ]
    ].sort_values("cluster_50_member_count", ascending=False).head(20).to_string(index=False)
)
