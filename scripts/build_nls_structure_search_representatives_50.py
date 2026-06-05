from pathlib import Path
import pandas as pd
import textwrap

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")

reps50_file = base / "tables" / "proteome_structure_search_representatives_50.tsv"
membership_file = base / "tables" / "proteome_cluster_membership_strategy.tsv"
nls_unique_file = base / "tables" / "nls_unique_candidates_no_sequence.tsv"

out_table = base / "tables" / "nls_structure_search_representatives_50.tsv"
out_fasta = base / "tables" / "nls_structure_search_representatives_50.fasta"
out_note = base / "notes" / "nls_structure_search_representatives_50_summary.txt"

reps50 = pd.read_csv(reps50_file, sep="\t")
membership = pd.read_csv(membership_file, sep="\t")
nls = pd.read_csv(nls_unique_file, sep="\t")

nls_ids = set(nls["canonical_id"].astype(str))

# Select 50% cluster representatives whose cluster contains at least one NLS candidate.
nls_reps = reps50[reps50["cluster_50_contains_nls_candidate"] == True].copy()

# Recover the NLS candidate IDs inside each 50% cluster.
nls_members = membership[membership["contains_nls_candidate"] == True].copy()

cluster_to_nls_ids = (
    nls_members.groupby("cluster_50_representative")["representative_canonical_id"]
    .apply(lambda x: ";".join(sorted(set(map(str, x)))))
    .reset_index(name="nls_candidate_ids_in_cluster")
)

nls_reps = nls_reps.merge(
    cluster_to_nls_ids,
    on="cluster_50_representative",
    how="left"
)

nls_reps["structure_search_scope"] = "nls_containing_cluster_50_representative"
nls_reps["structure_search_status"] = "pending"
nls_reps["matched_uniprot_id"] = ""
nls_reps["matched_structure_id"] = ""
nls_reps["structure_file_path"] = ""
nls_reps["notes"] = (
    "This representative belongs to a 50 percent identity cluster containing at least one validated NLS candidate."
)

cols_first = [
    "sequence_hash",
    "representative_canonical_id",
    "representative_record_id",
    "representative_organism",
    "sequence_length",
    "cluster_50_member_count",
    "cluster_50_contains_nls_candidate",
    "cluster_50_nls_sequence_count",
    "nls_candidate_ids_in_cluster",
    "cluster_50_representative",
    "structure_search_scope",
    "structure_search_status",
    "matched_uniprot_id",
    "matched_structure_id",
    "structure_file_path",
    "notes",
    "sequence",
]

cols_existing = [c for c in cols_first if c in nls_reps.columns]
remaining = [c for c in nls_reps.columns if c not in cols_existing]
nls_reps = nls_reps[cols_existing + remaining]

nls_reps.to_csv(out_table, sep="\t", index=False)

with out_fasta.open("w", encoding="utf-8") as f:
    for _, row in nls_reps.iterrows():
        header = (
            f">{row['sequence_hash']}"
            f"|representative={row['representative_canonical_id']}"
            f"|organism={row['representative_organism']}"
            f"|cluster_50_members={row['cluster_50_member_count']}"
            f"|nls_ids={row['nls_candidate_ids_in_cluster']}"
        )
        seq = str(row["sequence"]).strip()
        f.write(header + "\n")
        f.write("\n".join(textwrap.wrap(seq, width=80)) + "\n")

with out_note.open("w", encoding="utf-8") as f:
    f.write("NLS structure search representatives at 50 percent identity\n")
    f.write("===========================================================\n\n")
    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file documents the pilot structure-search input for validated NLS candidates using 50 percent MMseqs2 clusters.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"Validated unique NLS candidates: {len(nls)}\n")
    f.write(f"NLS-containing 50 percent clusters: {len(nls_reps)}\n")
    f.write(f"NLS sequence memberships represented: {int(nls_reps['cluster_50_nls_sequence_count'].sum())}\n")
    f.write(f"Total cluster members represented by these NLS-containing clusters: {int(nls_reps['cluster_50_member_count'].sum())}\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("The 51 validated NLS candidates fall into 49 clusters at the 50 percent identity level.\n")
    f.write("This creates a small pilot FASTA for sequence-based structure search before scaling to the full proteome-wide representative set.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Use nls_structure_search_representatives_50.fasta as the first pilot input for UniProt/SwissProt or structure-homology search.\n")
    f.write("If the pilot works, the same procedure can be scaled to proteome_structure_search_representatives_50.fasta.\n")

print("Wrote:", out_table)
print("Wrote:", out_fasta)
print("Wrote:", out_note)
print("NLS unique candidates:", len(nls))
print("NLS-containing 50% cluster representatives:", len(nls_reps))
print("NLS sequence memberships represented:", int(nls_reps["cluster_50_nls_sequence_count"].sum()))
print("Total cluster members represented:", int(nls_reps["cluster_50_member_count"].sum()))
