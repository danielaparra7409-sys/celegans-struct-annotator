from pathlib import Path
import pandas as pd

base = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")

table = base / "tables" / "proteome_structure_search_representatives_50.tsv"
out = base / "notes" / "proteome_structure_search_representatives_50_summary.txt"

df = pd.read_csv(table, sep="\t")

representatives = len(df)
unique_sequences_represented = int(df["cluster_50_member_count"].sum())
searches_saved = unique_sequences_represented - representatives
reduction_percent = round((searches_saved / unique_sequences_represented) * 100, 2)

clusters_with_nls = int(df["cluster_50_contains_nls_candidate"].sum())
nls_memberships = int(df["cluster_50_nls_sequence_count"].sum())

with out.open("w", encoding="utf-8") as f:
    f.write("Proteome structure search representatives at 50 percent identity\n")
    f.write("================================================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write("This file documents the first practical structure-search input derived from MMseqs2 clustering of the CeMbio/Jose proteome.\n")
    f.write("The goal is to reduce structure retrieval time by searching representative sequences instead of searching every unique protein sequence individually.\n\n")

    f.write("Current status\n")
    f.write("--------------\n")
    f.write(f"Unique protein sequences represented: {unique_sequences_represented}\n")
    f.write(f"Cluster representatives at 50 percent identity: {representatives}\n")
    f.write(f"Searches saved compared with unique sequences: {searches_saved}\n")
    f.write(f"Reduction percent: {reduction_percent}%\n")
    f.write(f"Clusters whose representative group contains at least one NLS candidate: {clusters_with_nls}\n")
    f.write(f"NLS candidate sequence memberships across 50 percent clusters: {nls_memberships}\n\n")

    f.write("Input files\n")
    f.write("-----------\n")
    f.write("proteome_cluster_membership_strategy.tsv: maps each unique sequence to its 90, 70, and 50 percent cluster representatives.\n")
    f.write("proteome_structure_search_representatives_50.tsv: table of 50 percent representatives selected for initial structure search.\n")
    f.write("proteome_structure_search_representatives_50.fasta: FASTA input for sequence-based structure search.\n\n")

    f.write("Interpretation\n")
    f.write("--------------\n")
    f.write("Exact deduplication alone was not enough to reduce the workload substantially, because the proteome still contained 55,174 unique sequences.\n")
    f.write("MMseqs2 clustering at 50 percent identity reduced the first-pass structure search space to 36,488 representative sequences.\n")
    f.write("This represents a reduction of approximately 33.87 percent compared with searching every unique sequence independently.\n")
    f.write("The 51 validated NLS candidates are still represented within this clustering system and fall into 49 clusters at the 50 percent identity level.\n\n")

    f.write("Decision\n")
    f.write("--------\n")
    f.write("Use proteome_structure_search_representatives_50.fasta as the first broad input for sequence-based structure retrieval.\n")
    f.write("Searches performed on these representatives can later be propagated to cluster members using proteome_cluster_membership_strategy.tsv.\n")
    f.write("For more conservative propagation, use the 70 percent or 90 percent clustering layers.\n\n")

    f.write("Next technical route\n")
    f.write("--------------------\n")
    f.write("1. Search representatives against UniProt or SwissProt using BLAST or MMseqs.\n")
    f.write("2. If UniProt IDs are found, attempt AlphaFold DB retrieval.\n")
    f.write("3. If no direct UniProt/AlphaFold match is available, search for homologs with known structures.\n")
    f.write("4. If no structure is recoverable, predict structures only for selected high-priority candidates.\n")

print("Wrote:", out)
print("Representatives:", representatives)
print("Unique sequences represented:", unique_sequences_represented)
print("Searches saved:", searches_saved)
print("Reduction percent:", reduction_percent)
print("Clusters with NLS candidates:", clusters_with_nls)
print("NLS memberships:", nls_memberships)
