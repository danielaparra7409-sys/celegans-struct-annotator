from pathlib import Path
import pandas as pd
import textwrap

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")

sequence_table = base / "tables" / "nls_unique_candidates.tsv"
category_table = base / "tables" / "nls_unique_candidates_functional_categories.tsv"

out_fasta = base / "tables" / "nls_unique_candidates.fasta"
out_summary = base / "notes" / "nls_unique_candidates_fasta_summary.txt"

seq_df = pd.read_csv(sequence_table, sep="\t")
cat_df = pd.read_csv(category_table, sep="\t")

# Keep only the category column from the categorized table.
cat_small = cat_df[["canonical_id", "functional_category_initial"]].copy()

df = seq_df.merge(cat_small, on="canonical_id", how="left")

records_written = 0
empty_sequences = 0

with out_fasta.open("w", encoding="utf-8") as f:
    for _, row in df.iterrows():
        canonical_id = str(row["canonical_id"]).strip()
        organism = str(row["organism_guess"]).strip()
        category = str(row["functional_category_initial"]).strip()
        seq = str(row["sequence"]).strip().replace("*", "")

        if not seq or seq.lower() == "nan":
            empty_sequences += 1
            continue

        header = f">{canonical_id}|organism={organism}|category={category}"
        f.write(header + "\n")
        f.write("\n".join(textwrap.wrap(seq, width=80)) + "\n")
        records_written += 1

with out_summary.open("w", encoding="utf-8") as f:
    f.write("NLS unique candidates FASTA summary\n")
    f.write("==================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file documents the FASTA export of the validated unique NLS candidate proteins.\n")
    f.write("The FASTA file is intended for sequence-based searches against UniProt, BLAST databases, AlphaFold-related workflows, or other homology tools.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"Input unique candidate rows: {len(df)}\n")
    f.write(f"FASTA records written: {records_written}\n")
    f.write(f"Empty sequences skipped: {empty_sequences}\n\n")

    f.write("FASTA file\n")
    f.write("----------\n")
    f.write(str(out_fasta) + "\n\n")

    f.write("Header format\n")
    f.write("-------------\n")
    f.write(">canonical_id|organism=organism_guess|category=functional_category_initial\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("All exported sequences correspond to validated unique NLS candidate proteins.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Use nls_unique_candidates.fasta as the clean input FASTA for sequence-based structure retrieval or homology search workflows.\n")

print("Wrote:", out_fasta)
print("Wrote:", out_summary)
print("Input rows:", len(df))
print("FASTA records written:", records_written)
print("Empty sequences skipped:", empty_sequences)
