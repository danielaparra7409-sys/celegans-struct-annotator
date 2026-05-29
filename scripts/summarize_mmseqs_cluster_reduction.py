from pathlib import Path
import pandas as pd

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
repo = Path.home() / "code" / "celegans-struct-annotator"

input_fasta = repo / "work" / "mmseqs" / "input" / "proteome_unique_sequence_representatives.fasta"
cluster_dir = repo / "work" / "mmseqs" / "clusters"

out_table = base / "tables" / "proteome_mmseqs_cluster_reduction_summary.tsv"
out_note = base / "notes" / "proteome_mmseqs_cluster_reduction_summary.txt"

def count_fasta_records(path):
    count = 0
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

exact_unique = count_fasta_records(input_fasta)

levels = [
    ("exact_unique", None, input_fasta),
    ("cluster_90", 0.90, cluster_dir / "proteome_cluster_90_rep_seq.fasta"),
    ("cluster_70", 0.70, cluster_dir / "proteome_cluster_70_rep_seq.fasta"),
    ("cluster_50", 0.50, cluster_dir / "proteome_cluster_50_rep_seq.fasta"),
]

rows = []

for label, identity, fasta_path in levels:
    representatives = count_fasta_records(fasta_path)
    searches_saved = exact_unique - representatives
    reduction_percent = round((searches_saved / exact_unique) * 100, 2)

    rows.append({
        "level": label,
        "min_seq_id": identity if identity is not None else "",
        "representative_sequences": representatives,
        "searches_saved_vs_exact_unique": searches_saved,
        "reduction_percent_vs_exact_unique": reduction_percent,
        "fasta_path": str(fasta_path),
    })

summary = pd.DataFrame(rows)
summary.to_csv(out_table, sep="\t", index=False)

with out_note.open("w", encoding="utf-8") as f:
    f.write("Proteome MMseqs cluster reduction summary\n")
    f.write("========================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file documents how much the full CeMbio/Jose proteome structure-search workload can be reduced by clustering protein sequences with MMseqs2.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"Exact unique protein sequences: {exact_unique}\n\n")

    f.write("Cluster reduction table\n")
    f.write("-----------------------\n")
    f.write(summary.to_string(index=False))
    f.write("\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("Exact sequence deduplication reduced the dataset only slightly, so exact deduplication alone is not enough to substantially reduce structure retrieval time.\n")
    f.write("MMseqs2 clustering provides a more useful reduction. The 50 percent identity clustering level gives the strongest reduction and can be used as an exploratory proteome-wide structure retrieval layer.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Use 50 percent identity clusters for an initial broad structure-search strategy across the full proteome.\n")
    f.write("Use 70 percent identity clusters for a more conservative intermediate layer.\n")
    f.write("Use 90 percent identity clusters for fine-grained analyses or high-confidence propagation.\n")
    f.write("The next step is to map each protein to its 90, 70, and 50 percent cluster representative so structure hits can be propagated from representatives to cluster members.\n")

print("Wrote:", out_table)
print("Wrote:", out_note)
print()
print(summary.to_string(index=False))
