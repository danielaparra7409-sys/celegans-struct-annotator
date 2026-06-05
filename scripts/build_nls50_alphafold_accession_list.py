from pathlib import Path
import pandas as pd

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")

annot_file = base / "tables" / "nls50_vs_swissprot_best_hits_annotated.tsv"

out_table = base / "tables" / "nls50_swissprot_accessions_for_alphafold.tsv"
out_note = base / "notes" / "nls50_swissprot_accessions_for_alphafold_summary.txt"

df = pd.read_csv(annot_file, sep="\t")

priority_order = {
    "strong_sequence_hit": 1,
    "moderate_sequence_hit": 2,
    "weak_but_possible_hit": 3,
    "low_confidence_hit": 4,
}

df["alphafold_lookup_priority"] = df["manual_review_confidence"].map(priority_order)

# Keep one row per SwissProt accession, prioritizing stronger hits.
df = df.sort_values(
    ["alphafold_lookup_priority", "evalue", "bits", "pident"],
    ascending=[True, True, False, False]
)

acc = df.drop_duplicates(subset=["swissprot_accession"], keep="first").copy()

acc["alphafold_status"] = "pending"
acc["alphafold_cif_url"] = ""
acc["alphafold_pdb_url"] = ""
acc["alphafold_cif_path"] = ""
acc["alphafold_pdb_path"] = ""
acc["alphafold_lookup_note"] = ""

cols_first = [
    "swissprot_accession",
    "swissprot_entry_from_fasta",
    "swissprot_protein_name",
    "swissprot_organism",
    "swissprot_gene",
    "manual_review_confidence",
    "alphafold_lookup_priority",
    "representative_canonical_id",
    "representative_organism",
    "nls_candidate_ids_in_cluster",
    "pident",
    "alnlen",
    "evalue",
    "bits",
    "alphafold_status",
    "alphafold_cif_url",
    "alphafold_pdb_url",
    "alphafold_cif_path",
    "alphafold_pdb_path",
    "alphafold_lookup_note",
]

cols_existing = [c for c in cols_first if c in acc.columns]
remaining = [c for c in acc.columns if c not in cols_existing]
acc = acc[cols_existing + remaining]

acc.to_csv(out_table, sep="\t", index=False)

with out_note.open("w", encoding="utf-8") as f:
    f.write("NLS 50 percent SwissProt accessions for AlphaFold lookup\n")
    f.write("=======================================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file lists unique SwissProt accessions recovered from the NLS 50 percent representative pilot search and prioritizes them for AlphaFold DB lookup.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"Best SwissProt hits: {len(df)}\n")
    f.write(f"Unique SwissProt accessions selected: {len(acc)}\n\n")

    f.write("Confidence categories\n")
    f.write("---------------------\n")
    f.write(acc["manual_review_confidence"].value_counts().to_string())
    f.write("\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("These accessions are not original CeMbio protein IDs. They are SwissProt homologs recovered by sequence search.\n")
    f.write("They can be used as entry points for AlphaFold DB lookup and homolog-based structural interpretation.\n")
    f.write("Any structure retrieved for these accessions should be interpreted as a homolog structure, not necessarily the exact structure of the original CeMbio protein.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Attempt AlphaFold DB retrieval first for strong and moderate SwissProt hits.\n")
    f.write("Weak or low-confidence hits should be reviewed manually before being used for structural propagation.\n")

print("Wrote:", out_table)
print("Wrote:", out_note)
print()
print("Unique accessions:", len(acc))
print("Confidence counts:")
print(acc["manual_review_confidence"].value_counts())
print()
print("Preview:")
print(
    acc[
        [
            "swissprot_accession",
            "swissprot_entry_from_fasta",
            "swissprot_protein_name",
            "swissprot_organism",
            "manual_review_confidence",
            "representative_canonical_id",
            "representative_organism",
            "pident",
            "evalue",
            "bits",
        ]
    ].to_string(index=False)
)
