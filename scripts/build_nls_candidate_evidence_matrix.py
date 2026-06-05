from pathlib import Path
import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(
    "/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia"
    "/Escritorio/Tareas/Trabajo C.elegans/CeMbio"
)
REPO = Path(__file__).resolve().parents[1]

# ── Output paths ─────────────────────────────────────────────────────────────
out_table_full   = BASE  / "tables"  / "nls_candidate_evidence_matrix.tsv"
out_note         = BASE  / "notes"   / "nls_candidate_evidence_matrix_summary.txt"
out_table_slim   = REPO  / "docs" / "tables" / "nls_candidate_evidence_matrix_summary.tsv"
fig_dir          = REPO  / "docs" / "figures" / "nls_evidence"
fig_dir.mkdir(parents=True, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load inputs
# ══════════════════════════════════════════════════════════════════════════════
func_df = pd.read_csv(
    BASE / "tables" / "nls_functional_categories_review_table.tsv", sep="\t"
)  # 51 rows — one per NLS candidate

reps_df = pd.read_csv(
    BASE / "tables" / "nls_structure_search_representatives_50.tsv", sep="\t"
)  # 49 rows — one per 50%-identity cluster containing NLS candidates

hits_df = pd.read_csv(
    BASE / "tables" / "nls50_vs_swissprot_best_hits_annotated.tsv", sep="\t"
)  # 19 rows — best SwissProt hit per representative (19 reps matched)

af_df = pd.read_csv(
    BASE / "tables" / "nls50_swissprot_accessions_for_alphafold.tsv", sep="\t"
)  # 18 rows — one per unique SwissProt accession, with AlphaFold status

# ── Extract pLDDT from lookup note ──────────────────────────────────────────
def _plddt(note):
    m = re.search(r"pLDDT=([0-9.]+)", str(note))
    return float(m.group(1)) if m else None

af_df["alphafold_plddt"] = af_df["alphafold_lookup_note"].apply(_plddt)


# ══════════════════════════════════════════════════════════════════════════════
# 2. Build NLS-candidate → 50%-representative mapping
#    nls_candidate_ids_in_cluster is semicolon-separated; may or may not include
#    the representative itself depending on whether the rep is an NLS candidate.
# ══════════════════════════════════════════════════════════════════════════════
rep_rows = []
for _, row in reps_df.iterrows():
    rep_id    = row["representative_canonical_id"]
    rep_org   = row["representative_organism"]
    mem_count = row["cluster_50_member_count"]
    for nls_id in str(row["nls_candidate_ids_in_cluster"]).split(";"):
        nls_id = nls_id.strip()
        if nls_id:
            rep_rows.append({
                "canonical_id":             nls_id,
                "cluster_50_representative": rep_id,
                "cluster_50_rep_organism":   rep_org,
                "cluster_50_member_count":   mem_count,
            })

rep_map = pd.DataFrame(rep_rows)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Prepare SwissProt and AlphaFold slices
# ══════════════════════════════════════════════════════════════════════════════
hits_slim = hits_df[[
    "representative_canonical_id", "swissprot_accession", "swissprot_protein_name",
    "swissprot_organism", "manual_review_confidence", "pident", "alnlen", "evalue", "bits",
]].rename(columns={"representative_canonical_id": "cluster_50_representative"})

af_slim = af_df[[
    "swissprot_accession", "alphafold_status", "alphafold_plddt", "alphafold_lookup_note",
]]


# ══════════════════════════════════════════════════════════════════════════════
# 4. Build matrix — left-join everything onto the 51 NLS candidates
# ══════════════════════════════════════════════════════════════════════════════
func_cols = [
    "canonical_id", "organism_guess", "functional_category_initial",
    "combined_functional_annotation", "Puntaje Extracelular", "Score NLS",
    "Substrate", "Possible Secretion System",
]
matrix = func_df[func_cols].copy()

# 4a. Add cluster representative
matrix = matrix.merge(rep_map, on="canonical_id", how="left")

# 4b. Add SwissProt best hit (one per representative, so 1:1 at candidate level)
matrix = matrix.merge(hits_slim, on="cluster_50_representative", how="left")

# 4c. Add AlphaFold status (keyed by SwissProt accession)
matrix = matrix.merge(af_slim, on="swissprot_accession", how="left")


# ══════════════════════════════════════════════════════════════════════════════
# 5. Derived columns
# ══════════════════════════════════════════════════════════════════════════════
def _interp_type(row):
    status  = str(row.get("alphafold_status", ""))
    acc     = str(row.get("swissprot_accession", ""))
    has_acc = acc not in ("", "nan", "None")
    if status == "found":
        return "homolog_SwissProt_AlphaFold_structure"
    if status == "not_found":
        return "SwissProt_homolog_not_in_AlphaFold_DB"
    if has_acc:
        return "SwissProt_homolog_no_structure_yet"
    return "no_SwissProt_hit"

def _priority(row):
    conf    = str(row.get("manual_review_confidence", ""))
    status  = str(row.get("alphafold_status", ""))
    acc     = str(row.get("swissprot_accession", ""))
    has_acc = acc not in ("", "nan", "None")

    if status == "found" and conf in ("strong_sequence_hit", "moderate_sequence_hit"):
        return "high_priority_review", (
            f"AlphaFold homolog structure available ({conf}); "
            "structure is of SwissProt homolog, not exact CeMbio protein"
        )
    if has_acc and status == "not_found":
        return "medium_priority_review", (
            f"SwissProt hit ({conf}) but homolog not in AlphaFold DB; "
            "consider direct prediction from CeMbio representative sequence"
        )
    if has_acc:
        return "medium_priority_review", (
            f"SwissProt hit present ({conf}); AlphaFold lookup pending or not attempted"
        )
    return "needs_broader_search", (
        "No SwissProt hit at 50% identity threshold; "
        "consider TrEMBL search or Foldseek on CeMbio representative"
    )

matrix["structure_interpretation_type"] = matrix.apply(_interp_type, axis=1)
matrix[["preliminary_review_priority", "priority_reason"]] = matrix.apply(
    lambda r: pd.Series(_priority(r)), axis=1
)


# ══════════════════════════════════════════════════════════════════════════════
# 6. Final column order and write outputs
# ══════════════════════════════════════════════════════════════════════════════
full_cols = [
    "canonical_id", "organism_guess",
    "functional_category_initial", "combined_functional_annotation",
    "Puntaje Extracelular", "Score NLS", "Substrate", "Possible Secretion System",
    "cluster_50_representative", "cluster_50_rep_organism", "cluster_50_member_count",
    "swissprot_accession", "swissprot_protein_name", "swissprot_organism",
    "manual_review_confidence", "pident", "alnlen", "evalue", "bits",
    "alphafold_status", "alphafold_plddt", "alphafold_lookup_note",
    "structure_interpretation_type",
    "preliminary_review_priority", "priority_reason",
]
matrix = matrix[[c for c in full_cols if c in matrix.columns]]
matrix.to_csv(out_table_full, sep="\t", index=False)

# Slim version for repo: drop sequence-containing, path-containing, or heavy cols
slim_drop = {"alphafold_lookup_note", "Substrate", "Possible Secretion System"}
slim_cols = [c for c in full_cols if c not in slim_drop and c in matrix.columns]
slim = matrix[slim_cols].copy()
# Trim long text fields
slim["combined_functional_annotation"] = (
    slim["combined_functional_annotation"].astype(str).str.slice(0, 120)
)
slim.to_csv(out_table_slim, sep="\t", index=False)


# ══════════════════════════════════════════════════════════════════════════════
# 7. Summary note
# ══════════════════════════════════════════════════════════════════════════════
pri_counts  = matrix["preliminary_review_priority"].value_counts()
af_counts   = matrix["alphafold_status"].fillna("no_SwissProt_hit").value_counts()
conf_counts = matrix["manual_review_confidence"].fillna("no_hit").value_counts()

with out_note.open("w", encoding="utf-8") as f:
    f.write("NLS candidate evidence matrix summary\n")
    f.write("======================================\n\n")

    f.write("Purpose\n-------\n")
    f.write(
        "Integrated evidence table for the 51 NLS candidate proteins from CeMbio, combining\n"
        "functional annotations, SwissProt homology, AlphaFold homolog structure availability,\n"
        "and cluster membership. No function or structure is confirmed here — all annotations\n"
        "are exploratory and heuristic.\n\n"
    )

    f.write(f"Rows: {len(matrix)}  (one per NLS candidate)\n")
    f.write(f"Columns: {len(matrix.columns)}\n\n")

    f.write("Preliminary review priority\n---------------------------\n")
    f.write(pri_counts.to_string())
    f.write("\n\n")

    f.write("AlphaFold status (homolog SwissProt)\n------------------------------------\n")
    f.write(af_counts.to_string())
    f.write("\n\n")

    f.write("SwissProt confidence categories\n--------------------------------\n")
    f.write(conf_counts.to_string())
    f.write("\n\n")

    f.write("Caution\n-------\n")
    f.write(
        "AlphaFold structures in this table are for SwissProt HOMOLOGS, not the original\n"
        "CeMbio proteins. Sequence identity to SwissProt ranges from ~32% to ~84%.\n"
        "Low-pLDDT regions indicate structural disorder and must be interpreted carefully.\n"
        "Functional categories are keyword-based and exploratory only.\n"
    )

print("Matrix built.")
print(f"  Rows: {len(matrix)}")
print(f"  Cols: {len(matrix.columns)}")
print()
print("Preliminary review priority:")
print(pri_counts.to_string())
print()
print("AlphaFold status:")
print(af_counts.to_string())


# ══════════════════════════════════════════════════════════════════════════════
# 8. Figures
# ══════════════════════════════════════════════════════════════════════════════
COLORS = {
    "found":                "#4C8BB5",
    "not_found":            "#E07B4F",
    "no_SwissProt_hit":     "#AAAAAA",
    "pending":              "#C4C48A",
    "high_priority_review": "#2E7D32",
    "medium_priority_review":"#F9A825",
    "needs_broader_search": "#AAAAAA",
    "strong_sequence_hit":  "#1565C0",
    "moderate_sequence_hit":"#42A5F5",
    "weak_but_possible_hit":"#FFCA28",
    "low_confidence_hit":   "#EF9A9A",
    "no_hit":               "#CCCCCC",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.labelsize":    11,
    "axes.titlesize":    12,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
})


# ── Fig 1: functional category × AlphaFold status ──────────────────────────
fig1_col = matrix["alphafold_status"].fillna("no_SwissProt_hit")
ct = pd.crosstab(matrix["functional_category_initial"], fig1_col)
status_order = ["found", "not_found", "pending", "no_SwissProt_hit"]
ct = ct.reindex(columns=[s for s in status_order if s in ct.columns], fill_value=0)

fig, ax = plt.subplots(figsize=(10, 5))
bottom = pd.Series(0, index=ct.index)
for status in ct.columns:
    color = COLORS.get(status, "#BBBBBB")
    ax.bar(ct.index, ct[status], bottom=bottom, label=status,
           color=color, edgecolor="white", linewidth=0.5)
    bottom += ct[status]

ax.set_xlabel("Functional category (initial)")
ax.set_ylabel("Number of NLS candidates")
ax.set_title("Functional category by AlphaFold homolog retrieval status")
ax.legend(title="AlphaFold status", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
plt.xticks(rotation=35, ha="right")
plt.tight_layout()
fig.savefig(fig_dir / "functional_category_by_alphafold_status.png", dpi=150)
plt.close()

# ── Fig 2: review priority counts ───────────────────────────────────────────
pri_order = ["high_priority_review", "medium_priority_review", "needs_broader_search"]
pri_vals  = [pri_counts.get(p, 0) for p in pri_order]
pri_colors = [COLORS.get(p, "#BBBBBB") for p in pri_order]
pri_labels = ["High priority\n(AlphaFold found)", "Medium priority\n(SwissProt hit)", "Needs broader\nsearch"]

fig, ax = plt.subplots(figsize=(7, 4))
bars = ax.barh(pri_labels, pri_vals, color=pri_colors, edgecolor="white", height=0.5)
for bar, val in zip(bars, pri_vals):
    ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
            str(val), va="center", fontsize=10)
ax.set_xlabel("Number of NLS candidates")
ax.set_title("Preliminary review priority — NLS pilot (n=51)")
ax.set_xlim(0, max(pri_vals) + 5)
plt.tight_layout()
fig.savefig(fig_dir / "review_priority_counts.png", dpi=150)
plt.close()

# ── Fig 3: SwissProt confidence counts ──────────────────────────────────────
conf_order = ["strong_sequence_hit", "moderate_sequence_hit",
              "weak_but_possible_hit", "low_confidence_hit", "no_hit"]
conf_vals   = [conf_counts.get(c, 0) for c in conf_order]
conf_colors = [COLORS.get(c, "#BBBBBB") for c in conf_order]
conf_labels = ["Strong\nhit", "Moderate\nhit", "Weak/\npossible", "Low\nconfidence", "No\nSwissProt hit"]

fig, ax = plt.subplots(figsize=(7, 4))
bars = ax.bar(conf_labels, conf_vals, color=conf_colors, edgecolor="white")
for bar, val in zip(bars, conf_vals):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
            str(val), ha="center", fontsize=10)
ax.set_ylabel("Number of NLS candidates")
ax.set_title("SwissProt sequence hit confidence — NLS 50% pilot (n=51)")
ax.set_ylim(0, max(conf_vals) + 4)
plt.tight_layout()
fig.savefig(fig_dir / "swissprot_confidence_counts.png", dpi=150)
plt.close()

# ── Fig 4: Extracellular score vs NLS score, colored by priority ─────────────
pri_palette = {
    "high_priority_review":  ("#2E7D32", "High priority", 80),
    "medium_priority_review": ("#F9A825", "Medium priority", 60),
    "needs_broader_search":  ("#AAAAAA", "Needs broader search", 40),
}

fig, ax = plt.subplots(figsize=(7, 6))
for pri, (color, label, zorder) in pri_palette.items():
    sub = matrix[matrix["preliminary_review_priority"] == pri]
    ax.scatter(
        sub["Puntaje Extracelular"], sub["Score NLS"],
        c=color, label=f"{label} (n={len(sub)})",
        s=60, alpha=0.85, edgecolors="white", linewidths=0.4, zorder=zorder,
    )

ax.set_xlabel("Extracellular score (Puntaje Extracelular)")
ax.set_ylabel("NLS score (Score NLS)")
ax.set_title("Extracellular vs NLS score by review priority\n(NLS candidates, n=51)")
ax.legend(fontsize=9, loc="lower right")
plt.tight_layout()
fig.savefig(fig_dir / "extracellular_vs_nls_by_priority.png", dpi=150)
plt.close()

print()
print("Figures saved:")
for f in sorted(fig_dir.glob("*.png")):
    print(f"  {f}")
print()
print("Output files:")
print(f"  {out_table_full}")
print(f"  {out_note}")
print(f"  {out_table_slim}")
