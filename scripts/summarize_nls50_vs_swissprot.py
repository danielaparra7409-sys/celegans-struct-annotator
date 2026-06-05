from pathlib import Path
import pandas as pd
import re

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
repo = Path.home() / "code" / "celegans-struct-annotator"

hits_file = repo / "work" / "mmseqs_search" / "nls_vs_swissprot" / "nls50_vs_swissprot.tsv"
query_table = base / "tables" / "nls_structure_search_representatives_50.tsv"

out_all = base / "tables" / "nls50_vs_swissprot_all_hits.tsv"
out_best = base / "tables" / "nls50_vs_swissprot_best_hits.tsv"
out_nohit = base / "tables" / "nls50_vs_swissprot_no_hits.tsv"
out_note = base / "notes" / "nls50_vs_swissprot_summary.txt"

cols = [
    "query",
    "target",
    "pident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
]

hits = pd.read_csv(hits_file, sep="\t", header=None, names=cols)
queries = pd.read_csv(query_table, sep="\t")

def clean_query_id(value):
    return str(value).split("|")[0]

def parse_swissprot_accession(target):
    """
    SwissProt FASTA headers usually look like:
    sp|Q6GZX4|001R_FRG3G ...
    MMseqs may keep the full FASTA ID.
    """
    value = str(target)
    parts = value.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[1]
    return parts[0]

def parse_swissprot_entry(target):
    value = str(target)
    parts = value.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[2].split()[0]
    return ""

hits["query_hash"] = hits["query"].map(clean_query_id)
hits["swissprot_accession"] = hits["target"].map(parse_swissprot_accession)
hits["swissprot_entry"] = hits["target"].map(parse_swissprot_entry)

# Add query metadata.
query_meta_cols = [
    "sequence_hash",
    "representative_canonical_id",
    "representative_organism",
    "sequence_length",
    "cluster_50_member_count",
    "cluster_50_nls_sequence_count",
    "nls_candidate_ids_in_cluster",
]

query_meta_cols = [c for c in query_meta_cols if c in queries.columns]

merged = hits.merge(
    queries[query_meta_cols],
    left_on="query_hash",
    right_on="sequence_hash",
    how="left",
)

# Sort best hits by evalue ascending, then bits descending, then pident descending.
merged = merged.sort_values(
    ["query_hash", "evalue", "bits", "pident"],
    ascending=[True, True, False, False],
)

best = merged.groupby("query_hash", as_index=False).first()

all_query_hashes = set(queries["sequence_hash"].astype(str))
hit_query_hashes = set(best["query_hash"].astype(str))
nohit_hashes = sorted(all_query_hashes - hit_query_hashes)

nohit = queries[queries["sequence_hash"].astype(str).isin(nohit_hashes)].copy()

merged.to_csv(out_all, sep="\t", index=False)
best.to_csv(out_best, sep="\t", index=False)
nohit.to_csv(out_nohit, sep="\t", index=False)

with out_note.open("w", encoding="utf-8") as f:
    f.write("NLS 50 percent representatives versus SwissProt summary\n")
    f.write("======================================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file summarizes the first pilot sequence-based search of NLS-containing 50 percent cluster representatives against SwissProt using MMseqs2.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"NLS 50 percent representatives searched: {len(queries)}\n")
    f.write(f"Total SwissProt hits passing MMseqs thresholds: {len(merged)}\n")
    f.write(f"Representatives with at least one SwissProt hit: {len(best)}\n")
    f.write(f"Representatives without SwissProt hit: {len(nohit)}\n\n")

    if len(best) > 0:
        f.write("Best-hit identity summary\n")
        f.write("-------------------------\n")
        f.write(f"Min percent identity: {best['pident'].min():.2f}\n")
        f.write(f"Max percent identity: {best['pident'].max():.2f}\n")
        f.write(f"Mean percent identity: {best['pident'].mean():.2f}\n")
        f.write(f"Median percent identity: {best['pident'].median():.2f}\n\n")

        f.write("Best-hit bitscore summary\n")
        f.write("-------------------------\n")
        f.write(f"Min bitscore: {best['bits'].min():.2f}\n")
        f.write(f"Max bitscore: {best['bits'].max():.2f}\n")
        f.write(f"Mean bitscore: {best['bits'].mean():.2f}\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("This pilot tests whether the NLS-containing representatives can recover curated UniProtKB/SwissProt homologs.\n")
    f.write("SwissProt hits can provide reviewed accessions, functional clues, and possible entry points for AlphaFold DB retrieval or homolog-based structural interpretation.\n")
    f.write("However, SwissProt hits should be interpreted carefully: a hit does not automatically prove identical function, localization, secretion, or nuclear targeting.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Use nls50_vs_swissprot_best_hits.tsv as the first table for manual review of candidate homologs.\n")
    f.write("Representatives without SwissProt hits may require broader UniProt/TrEMBL search, structure-homology search, or selected structure prediction.\n")

print("Wrote:", out_all)
print("Wrote:", out_best)
print("Wrote:", out_nohit)
print("Wrote:", out_note)
print()
print("Queries:", len(queries))
print("Total hits:", len(merged))
print("Queries with hit:", len(best))
print("Queries without hit:", len(nohit))
print()
print("Best hits preview:")
print(
    best[
        [
            "representative_canonical_id",
            "representative_organism",
            "swissprot_accession",
            "swissprot_entry",
            "pident",
            "alnlen",
            "evalue",
            "bits",
        ]
    ].head(20).to_string(index=False)
)
