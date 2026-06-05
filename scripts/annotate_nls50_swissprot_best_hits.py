from pathlib import Path
import pandas as pd

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
repo = Path.home() / "code" / "celegans-struct-annotator"

best_file = base / "tables" / "nls50_vs_swissprot_best_hits.tsv"
swissprot_fasta = repo / "work" / "databases" / "swissprot" / "uniprot_sprot.fasta"

out_table = base / "tables" / "nls50_vs_swissprot_best_hits_annotated.tsv"
out_note = base / "notes" / "nls50_vs_swissprot_best_hits_annotated_summary.txt"

best = pd.read_csv(best_file, sep="\t")

# Force accession fields to string
best["swissprot_accession"] = best["swissprot_accession"].astype(str)
if "swissprot_entry" in best.columns:
    best["swissprot_entry"] = best["swissprot_entry"].astype(str)

def parse_header(header):
    h = header.strip().lstrip(">")
    first, *rest = h.split(" ", 1)
    desc = rest[0] if rest else ""

    parts = first.split("|")
    db = parts[0] if len(parts) > 0 else ""
    acc = parts[1] if len(parts) > 1 else first
    entry = parts[2] if len(parts) > 2 else ""

    protein_name = desc.split(" OS=")[0].strip()

    organism = ""
    if " OS=" in desc:
        organism = desc.split(" OS=", 1)[1].split(" OX=", 1)[0].strip()

    gene = ""
    if " GN=" in desc:
        gene = desc.split(" GN=", 1)[1].split(" ", 1)[0].strip()

    pe = ""
    if " PE=" in desc:
        pe = desc.split(" PE=", 1)[1].split(" ", 1)[0].strip()

    sv = ""
    if " SV=" in desc:
        sv = desc.split(" SV=", 1)[1].split(" ", 1)[0].strip()

    return {
        "swissprot_db_from_fasta": db,
        "swissprot_accession": str(acc),
        "swissprot_entry_from_fasta": entry,
        "swissprot_protein_name": protein_name,
        "swissprot_organism": organism,
        "swissprot_gene": gene,
        "swissprot_PE": pe,
        "swissprot_SV": sv,
        "swissprot_full_header": h,
    }

records = []
with swissprot_fasta.open("r", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            records.append(parse_header(line))

headers = pd.DataFrame(records)
headers["swissprot_accession"] = headers["swissprot_accession"].astype(str)

# Merge only by accession. This avoids dtype problems with swissprot_entry.
annot = best.merge(
    headers,
    on="swissprot_accession",
    how="left",
)

def confidence(row):
    pident = row["pident"]
    alnlen = row["alnlen"]
    evalue = row["evalue"]

    if pident >= 50 and alnlen >= 100 and evalue <= 1e-20:
        return "strong_sequence_hit"
    if pident >= 30 and alnlen >= 100 and evalue <= 1e-10:
        return "moderate_sequence_hit"
    if evalue <= 1e-5:
        return "weak_but_possible_hit"
    return "low_confidence_hit"

annot["manual_review_confidence"] = annot.apply(confidence, axis=1)

cols_first = [
    "representative_canonical_id",
    "representative_organism",
    "sequence_length",
    "nls_candidate_ids_in_cluster",
    "swissprot_accession",
    "swissprot_entry",
    "swissprot_entry_from_fasta",
    "swissprot_protein_name",
    "swissprot_organism",
    "swissprot_gene",
    "pident",
    "alnlen",
    "evalue",
    "bits",
    "manual_review_confidence",
    "query_hash",
    "target",
    "swissprot_full_header",
]

cols_existing = [c for c in cols_first if c in annot.columns]
remaining = [c for c in annot.columns if c not in cols_existing]
annot = annot[cols_existing + remaining]

annot.to_csv(out_table, sep="\t", index=False)

with out_note.open("w", encoding="utf-8") as f:
    f.write("Annotated NLS 50 percent SwissProt best hits summary\n")
    f.write("===================================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file summarizes the annotated best SwissProt hits recovered for NLS-containing 50 percent cluster representatives.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"Best hits annotated: {len(annot)}\n")
    f.write(f"Unique SwissProt accessions: {annot['swissprot_accession'].nunique()}\n")
    f.write(f"Hits with parsed protein name: {annot['swissprot_protein_name'].notna().sum()}\n\n")

    f.write("Manual confidence categories\n")
    f.write("----------------------------\n")
    f.write(annot["manual_review_confidence"].value_counts().to_string())
    f.write("\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("These SwissProt hits provide curated homolog candidates for manual functional review and possible AlphaFold DB lookup.\n")
    f.write("The confidence labels are heuristic and should not be treated as definitive biological assignments.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Use nls50_vs_swissprot_best_hits_annotated.tsv for manual review before attempting AlphaFold DB retrieval or broader UniProt/TrEMBL search.\n")

print("Wrote:", out_table)
print("Wrote:", out_note)
print()
print("Best hits annotated:", len(annot))
print("Confidence counts:")
print(annot["manual_review_confidence"].value_counts())
print()
print("Preview:")
print(
    annot[
        [
            "representative_canonical_id",
            "representative_organism",
            "swissprot_accession",
            "swissprot_entry_from_fasta",
            "swissprot_protein_name",
            "swissprot_organism",
            "pident",
            "alnlen",
            "evalue",
            "bits",
            "manual_review_confidence",
        ]
    ].head(20).to_string(index=False)
)
