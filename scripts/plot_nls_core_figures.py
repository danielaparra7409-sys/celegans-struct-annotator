from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio")
MASTER = BASE / "tables" / "cembio_master_table_minimal_no_sequence.tsv"
NLS = BASE / "tables" / "nls_functional_categories_review_table.tsv"
OUTDIR = BASE / "results" / "figures" / "nls_summary"

OUTDIR.mkdir(parents=True, exist_ok=True)

master_df = pd.read_csv(MASTER, sep="\t")
nls_df = pd.read_csv(NLS, sep="\t")

print("Master rows:", len(master_df))
print("NLS rows:", len(nls_df))
print("Output dir:", OUTDIR)

def save_bar(series, title, xlabel, ylabel, outfile, rotation=45):
    fig, ax = plt.subplots(figsize=(10, 6))
    series.plot(kind="bar", ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=rotation)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

def save_hist(series, title, xlabel, ylabel, outfile, bins=10):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(series.dropna(), bins=bins)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

def save_scatter(df, xcol, ycol, title, xlabel, ylabel, outfile):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df[xcol], df[ycol])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

def save_heatmap(matrix, title, xlabel, ylabel, outfile):
    fig, ax = plt.subplots(figsize=(12, 7))
    im = ax.imshow(matrix.values, aspect="auto")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Count")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix.iloc[i, j]), ha="center", va="center", fontsize=8)

    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

# 1. Proteínas totales por organismo
master_counts = (
    master_df["organism_guess"]
    .value_counts()
    .sort_values(ascending=False)
)
save_bar(
    master_counts,
    "Total proteins per organism in CeMbio master dataset",
    "Organism",
    "Protein count",
    OUTDIR / "01_master_proteins_per_organism.png"
)

# 2. Candidatas NLS únicas por organismo
nls_counts = (
    nls_df["organism_guess"]
    .value_counts()
    .sort_values(ascending=False)
)
save_bar(
    nls_counts,
    "Validated unique NLS candidates per organism",
    "Organism",
    "NLS candidate count",
    OUTDIR / "02_nls_candidates_per_organism.png"
)

# 3. Categorías funcionales iniciales
category_counts = (
    nls_df["functional_category_initial"]
    .value_counts()
    .sort_values(ascending=False)
)
save_bar(
    category_counts,
    "Initial functional categories of validated NLS candidates",
    "Functional category",
    "Count",
    OUTDIR / "03_functional_category_counts.png"
)

# 4. Heatmap organismo x categoría
heatmap_df = pd.crosstab(
    nls_df["organism_guess"],
    nls_df["functional_category_initial"]
)
heatmap_df = heatmap_df.loc[
    sorted(heatmap_df.index),
    sorted(heatmap_df.columns)
]
save_heatmap(
    heatmap_df,
    "Validated NLS candidates by organism and functional category",
    "Functional category",
    "Organism",
    OUTDIR / "04_heatmap_organism_by_category.png"
)

# 5. Distribución puntaje extracelular
save_hist(
    nls_df["Puntaje Extracelular"],
    "Distribution of extracellular scores",
    "Extracellular score",
    "Frequency",
    OUTDIR / "05_extracellular_score_distribution.png",
    bins=10
)

# 6. Distribución score NLS
save_hist(
    nls_df["Score NLS"],
    "Distribution of NLS scores",
    "NLS score",
    "Frequency",
    OUTDIR / "06_nls_score_distribution.png",
    bins=10
)

# 7. Dispersión score extracelular vs score NLS
save_scatter(
    nls_df,
    "Puntaje Extracelular",
    "Score NLS",
    "Extracellular score vs NLS score",
    "Extracellular score",
    "NLS score",
    OUTDIR / "07_extracellular_vs_nls_scatter.png"
)

summary_lines = []
summary_lines.append("NLS core figures summary")
summary_lines.append("========================")
summary_lines.append("")
summary_lines.append(f"Master rows: {len(master_df)}")
summary_lines.append(f"NLS rows: {len(nls_df)}")
summary_lines.append("")
summary_lines.append("Figures generated:")
for p in sorted(OUTDIR.glob("*.png")):
    summary_lines.append(str(p.name))

summary_path = OUTDIR / "figures_summary.txt"
summary_path.write_text("\n".join(summary_lines), encoding="utf-8")

print()
print("Generated files:")
for p in sorted(OUTDIR.glob("*")):
    print(p.name)
