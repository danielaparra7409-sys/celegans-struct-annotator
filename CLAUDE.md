# CLAUDE.md — celegans-struct-annotator

## Project goal
Build a reproducible bioinformatics workflow to identify, annotate, and retrieve structures for bacterial NLS-candidate proteins associated with *C. elegans* gut microbiome (CeMbio). End goal: scientific article by **2026-10-30**.

---

## Environment setup (always run first)
```bash
cd ~/code/celegans-struct-annotator
source .venv/bin/activate
BASE="/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/CeMbio"
```

---

## Key paths
| Purpose | Path |
|---|---|
| Repo | `~/code/celegans-struct-annotator` |
| Data root | `$BASE` |
| Scripts | `scripts/` |
| Checkpoints | `docs/cembio_checkpoints/` |
| Figures | `docs/figures/nls_summary/` |
| Tables (repo) | `docs/tables/` |
| Tables (local) | `$BASE/tables/` |
| Notes (local) | `$BASE/notes/` |
| MMseqs work | `work/mmseqs/`, `work/mmseqs_search/` |
| SwissProt DB | `work/databases/swissprot/` |

**Never commit:** `work/`, `*.fasta`, `*.fa`, `*.faa`, `*.xlsx`, large derived tables, SwissProt files.

---

## Current state (as of last checkpoint)

### Numbers
| Item | Count |
|---|---|
| Master proteins | 55,534 |
| Organisms | 12 |
| Exact unique sequences | 55,174 |
| 90% cluster reps | 52,399 |
| 70% cluster reps | 46,324 |
| 50% cluster reps | 36,488 |
| Validated NLS candidates | 51 |
| NLS-containing 50% clusters | 49 |
| Members in NLS clusters | 58 |
| NLS reps vs SwissProt hits | 19/49 |
| Unique SwissProt accessions | 18 |
| Priority accessions (strong+moderate) | **9** |

### Completed steps
1. Master proteome table (55,534 proteins, 12 organisms, no duplicates)
2. NLS candidate validation (51 unique, all matched master table)
3. Integration with José's Excel annotations (100% match)
4. Keyword-based functional classification (9 categories)
5. Exploratory figures (7 PNGs in `docs/figures/nls_summary/`)
6. MMseqs2 clustering at 90/70/50% identity
7. Cluster membership strategy table
8. 50% representative FASTA + table (proteome-wide and NLS-specific)
9. SwissProt DB setup (574,627 records)
10. NLS 50% reps vs SwissProt search → 585 hits, 19 reps matched
11. Annotated best hits + confidence categories
12. AlphaFold accession list (18 accessions, 9 priority)

### Next step
**AlphaFold DB retrieval** for 9 priority SwissProt accessions (4 strong + 5 moderate confidence).

Route:
1. Retrieve `.cif`/`.pdb` for priority accessions → save locally
2. Update AlphaFold accession table with status + file paths
3. Manual review for weak/low-confidence hits
4. For 30 no-hit NLS reps: broader TrEMBL or Foldseek search
5. Structure prediction only for selected high-priority candidates
6. Scale to 36,488 proteome-wide reps after pilot validates

---

## Key tables and their purpose
| File | Purpose |
|---|---|
| `tables/proteome_cluster_membership_strategy.tsv` | Maps each unique seq → 90/70/50% rep |
| `tables/proteome_structure_search_representatives_50.tsv` | 36,488 reps for full search |
| `tables/nls_structure_search_representatives_50.tsv` | 49 NLS reps for pilot |
| `tables/nls50_vs_swissprot_best_hits_annotated.tsv` | 19 annotated SwissProt hits |
| `tables/nls50_swissprot_accessions_for_alphafold.tsv` | 18 accessions + confidence labels |

---

## Confidence labels (heuristic, not final biology)
| Label | Meaning |
|---|---|
| `strong_sequence_hit` | Reliable homolog candidate |
| `moderate_sequence_hit` | Useful, needs review |
| `weak_but_possible_hit` | Distant/partial, treat cautiously |
| `low_confidence_hit` | Caution, verify before use |

---

## Scientific cautions (always apply)
- NLS prediction ≠ proven nuclear localization
- Extracellular score ≠ proven secretion
- SwissProt homolog ≠ identical function
- 50% identity cluster ≠ identical function
- AlphaFold structure of SwissProt hit = **homolog structure**, not exact CeMbio protein
- Functional categories are keyword-based/exploratory only

---

## Tools
| Tool | Version | Use |
|---|---|---|
| Python | — | All data processing and figures |
| MMseqs2 | 15-6f452+ds-2 | Clustering + sequence search |
| Git/GitHub | — | Version control (scripts, notes, figures, small tables) |

---

## Git workflow
```bash
git status --short
git add <file>
git commit -m "concise message"
git push origin main
```

---

## Quick inspection commands
```bash
# Priority AlphaFold accessions
cat "$BASE/notes/nls50_swissprot_accessions_for_alphafold_summary.txt"

# SwissProt search results
cat "$BASE/notes/nls50_vs_swissprot_best_hits_annotated_summary.txt"

# Cluster reduction summary
cat "$BASE/notes/proteome_mmseqs_cluster_reduction_summary.txt"

# Open figures
explorer.exe "$(wslpath -w "$BASE/results/figures/nls_summary")"
```

---

## Article deadline: 2026-10-30
