"""Summarize proteome-wide 50% representative vs SwissProt MMseqs2 search results."""
from pathlib import Path
import pandas as pd

BASE = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
REPO = Path.home() / "code" / "celegans-struct-annotator"

hits_file   = REPO / "work/mmseqs_search/proteome50_vs_swissprot/proteome50_vs_swissprot.tsv"
query_table = REPO / "docs/tables/proteome_structure_search_representatives_50.tsv"

out_all   = BASE / "tables/proteome50_vs_swissprot_all_hits.tsv"
out_best  = BASE / "tables/proteome50_vs_swissprot_best_hits.tsv"
out_nohit = BASE / "tables/proteome50_vs_swissprot_no_hits.tsv"
out_note  = BASE / "notes/proteome50_vs_swissprot_summary.txt"
out_best_repo = REPO / "docs/tables/proteome50_vs_swissprot_best_hits.tsv"

COLS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]

print("Loading hits...")
hits = pd.read_csv(hits_file, sep="\t", header=None, names=COLS)
print(f"  Total hits loaded: {len(hits):,}")

print("Loading query table...")
queries = pd.read_csv(query_table, sep="\t")
print(f"  Total queries: {len(queries):,}")

# Parse sequence hash (first field before |) and SwissProt accession
hits["query_hash"] = hits["query"].str.split("|").str[0]
hits["swissprot_accession"] = hits["target"].apply(
    lambda v: v.split("|")[1] if len(v.split("|")) >= 3 else v.split("|")[0]
)
hits["swissprot_entry"] = hits["target"].apply(
    lambda v: v.split("|")[2].split()[0] if len(v.split("|")) >= 3 else ""
)

# Confidence categories — heuristic only, consistent with NLS pilot labels
def confidence(row):
    if row["evalue"] < 1e-50 and row["pident"] >= 40:
        return "strong_sequence_hit"
    elif row["evalue"] < 1e-20 and row["pident"] >= 25:
        return "moderate_sequence_hit"
    elif row["evalue"] < 1e-5:
        return "weak_but_possible_hit"
    else:
        return "low_confidence_hit"

print("Assigning confidence labels...")
hits["confidence"] = hits.apply(confidence, axis=1)

# Merge with query metadata
META_COLS = [
    "sequence_hash", "representative_canonical_id", "representative_organism",
    "sequence_length", "cluster_50_member_count",
    "cluster_50_contains_nls_candidate", "cluster_50_nls_sequence_count",
]
META_COLS = [c for c in META_COLS if c in queries.columns]

print("Merging with query metadata...")
merged = hits.merge(
    queries[META_COLS],
    left_on="query_hash",
    right_on="sequence_hash",
    how="left",
)
merged = merged.sort_values(
    ["query_hash", "evalue", "bits", "pident"],
    ascending=[True, True, False, False],
)

# Best hit per representative (lowest evalue, highest bits, highest pident)
best = merged.groupby("query_hash", as_index=False).first()

all_hashes  = set(queries["sequence_hash"].astype(str))
hit_hashes  = set(best["query_hash"].astype(str))
nohit_hashes = all_hashes - hit_hashes
nohit = queries[queries["sequence_hash"].astype(str).isin(nohit_hashes)].copy()

pct_hit   = len(hit_hashes) / len(all_hashes) * 100
pct_nohit = 100 - pct_hit

# Best hits for repo: exclude the `sequence` column (too large for git)
REPO_COLS = [c for c in best.columns if c != "sequence"]
best_repo = best[REPO_COLS]

print("Writing output files...")
merged.to_csv(out_all,       sep="\t", index=False)
best.to_csv(out_best,        sep="\t", index=False)
nohit.to_csv(out_nohit,      sep="\t", index=False)
best_repo.to_csv(out_best_repo, sep="\t", index=False)
print(f"  {out_all}")
print(f"  {out_best}")
print(f"  {out_nohit}")
print(f"  {out_best_repo}")

# Summary text
with out_note.open("w", encoding="utf-8") as f:
    f.write("Proteome 50 percent representatives versus SwissProt — summary\n")
    f.write("==============================================================\n\n")

    f.write("Overview\n")
    f.write("--------\n")
    f.write(f"Representatives searched:         {len(all_hashes):>8,}\n")
    f.write(f"Representatives with >= 1 hit:    {len(hit_hashes):>8,}  ({pct_hit:.1f}%)\n")
    f.write(f"Representatives without hit:      {len(nohit_hashes):>8,}  ({pct_nohit:.1f}%)\n")
    f.write(f"Total hits (all thresholds):      {len(merged):>8,}\n")
    f.write(f"Unique SwissProt accessions:      {hits['swissprot_accession'].nunique():>8,}\n\n")

    f.write("Confidence distribution (best hits per representative)\n")
    f.write("------------------------------------------------------\n")
    for label, n in best["confidence"].value_counts().items():
        f.write(f"  {label:<30} {n:>7,}\n")
    f.write("\n")

    f.write("Identity distribution — best hits (% identity)\n")
    f.write("-----------------------------------------------\n")
    bins   = [0, 20, 30, 40, 50, 60, 70, 80, 90, 101]
    labels = ["<20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">90%"]
    best["pident_bin"] = pd.cut(best["pident"], bins=bins, labels=labels, right=False)
    for b, n in best["pident_bin"].value_counts(sort=False).items():
        f.write(f"  {str(b):<10} {n:>7,}\n")
    f.write("\n")

    f.write("E-value distribution — best hits\n")
    f.write("--------------------------------\n")
    ev_thresholds = [
        (best["evalue"] < 1e-100,              "< 1e-100  (very strong)"),
        ((best["evalue"] >= 1e-100) & (best["evalue"] < 1e-50),  "1e-100 to 1e-50"),
        ((best["evalue"] >= 1e-50)  & (best["evalue"] < 1e-20),  "1e-50  to 1e-20"),
        ((best["evalue"] >= 1e-20)  & (best["evalue"] < 1e-5),   "1e-20  to 1e-5"),
        ((best["evalue"] >= 1e-5)   & (best["evalue"] < 1e-3),   "1e-5   to 1e-3"),
    ]
    for mask, label in ev_thresholds:
        f.write(f"  {label:<30} {mask.sum():>7,}\n")
    f.write("\n")

    f.write("Bitscore distribution — best hits\n")
    f.write("---------------------------------\n")
    f.write(f"  Min:    {best['bits'].min():>8.1f}\n")
    f.write(f"  Median: {best['bits'].median():>8.1f}\n")
    f.write(f"  Mean:   {best['bits'].mean():>8.1f}\n")
    f.write(f"  Max:    {best['bits'].max():>8.1f}\n\n")

    f.write("Output files\n")
    f.write("------------\n")
    f.write(f"  All hits:        {out_all}\n")
    f.write(f"  Best hits:       {out_best}\n")
    f.write(f"  No hits:         {out_nohit}\n")
    f.write(f"  Best hits (repo):{out_best_repo}\n\n")

    f.write("Scientific caution\n")
    f.write("------------------\n")
    f.write("SwissProt hits represent probable functional homology, NOT confirmed function.\n")
    f.write("50% identity clustering does not imply identical function across members.\n")
    f.write("Confidence labels are heuristic thresholds, not experimentally validated.\n")
    f.write("This search is exploratory — intended to prioritize representatives for\n")
    f.write("structural analysis, not to assign definitive functional annotations.\n")

print(f"\nWrote summary: {out_note}")

# Print to stdout for immediate review
print("\n" + "=" * 62)
print("PROTEOME 50% vs SwissProt — RESULTS")
print("=" * 62)
print(f"Representatives searched:         {len(all_hashes):>8,}")
print(f"Representatives with hit:         {len(hit_hashes):>8,}  ({pct_hit:.1f}%)")
print(f"Representatives without hit:      {len(nohit_hashes):>8,}  ({pct_nohit:.1f}%)")
print(f"Total hits (all):                 {len(merged):>8,}")
print(f"Unique SwissProt accessions:      {hits['swissprot_accession'].nunique():>8,}")
print()
print("Confidence (best hits):")
for label, n in best["confidence"].value_counts().items():
    print(f"  {label:<30} {n:>7,}")
print()
print("Identity (best hits):")
for b, n in best["pident_bin"].value_counts(sort=False).items():
    print(f"  {str(b):<10} {n:>7,}")
