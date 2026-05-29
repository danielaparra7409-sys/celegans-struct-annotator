from pathlib import Path
import pandas as pd
import hashlib
import textwrap

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")

master_file = base / "tables" / "cembio_master_table_minimal.tsv"
nls_file = base / "tables" / "nls_unique_candidates_no_sequence.tsv"

out_table = base / "tables" / "proteome_unique_sequence_representatives.tsv"
out_fasta = base / "tables" / "proteome_unique_sequence_representatives.fasta"

master = pd.read_csv(master_file, sep="\t")
nls = pd.read_csv(nls_file, sep="\t")

nls_ids = set(nls["canonical_id"].astype(str))

def canonical_id(value):
    if pd.isna(value):
        return None
    v = str(value).strip().lstrip(">")
    if not v or v.lower() == "nan":
        return None
    v = v.split()[0]
    if v.startswith("ENA|"):
        parts = v.split("|")
        if len(parts) >= 3:
            return parts[2]
    return v

def clean_sequence(seq):
    if pd.isna(seq):
        return ""
    return str(seq).strip().replace("*", "").upper()

def seq_hash(seq):
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()

master["canonical_id"] = master["record_id"].map(canonical_id)
master["sequence_clean"] = master["sequence"].map(clean_sequence)
master["sequence_hash"] = master["sequence_clean"].map(seq_hash)
master["has_sequence"] = master["sequence_clean"].str.len() > 0
master["is_unique_nls_candidate"] = master["canonical_id"].astype(str).isin(nls_ids)

valid = master[master["has_sequence"]].copy()

grouped = (
    valid.groupby("sequence_hash")
    .agg(
        representative_record_id=("record_id", "first"),
        representative_canonical_id=("canonical_id", "first"),
        representative_organism=("organism_guess", "first"),
        sequence_length=("sequence_length", "first"),
        sequence=("sequence_clean", "first"),
        protein_count=("record_id", "count"),
        organism_count=("organism_guess", lambda x: x.nunique()),
        contains_nls_candidate=("is_unique_nls_candidate", "any"),
        nls_candidate_count=("is_unique_nls_candidate", "sum"),
    )
    .reset_index()
)

grouped = grouped.sort_values(
    ["contains_nls_candidate", "protein_count", "sequence_length"],
    ascending=[False, False, False]
)

grouped.to_csv(out_table, sep="\t", index=False)

with out_fasta.open("w", encoding="utf-8") as f:
    for _, row in grouped.iterrows():
        header = (
            f">{row['sequence_hash']}|representative={row['representative_canonical_id']}"
            f"|organism={row['representative_organism']}"
            f"|protein_count={row['protein_count']}"
            f"|contains_nls={row['contains_nls_candidate']}"
        )
        f.write(header + "\n")
        f.write("\n".join(textwrap.wrap(row["sequence"], width=80)) + "\n")

print("Master proteins:", len(master))
print("Proteins with sequence:", len(valid))
print("Unique sequences:", len(grouped))
print("Repeated sequence groups:", (grouped["protein_count"] > 1).sum())
print("Proteins represented in repeated groups:", grouped.loc[grouped["protein_count"] > 1, "protein_count"].sum())
print("Unique sequence groups containing NLS candidates:", grouped["contains_nls_candidate"].sum())
print("Wrote:", out_table)
print("Wrote:", out_fasta)
