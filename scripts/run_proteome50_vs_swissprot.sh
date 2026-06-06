#!/usr/bin/env bash
# run_proteome50_vs_swissprot.sh
# Search all 36,488 proteome-wide 50% cluster representatives against SwissProt using MMseqs2.
# Output: work/mmseqs_search/proteome50_vs_swissprot/proteome50_vs_swissprot.tsv
set -euo pipefail

REPO="$HOME/code/celegans-struct-annotator"
QUERY_FASTA="$REPO/work/mmseqs/clusters/proteome_cluster_50_rep_seq.fasta"
SWISSPROT_DB="$REPO/work/databases/swissprot/mmseqs/swissprot_db"
OUTDIR="$REPO/work/mmseqs_search/proteome50_vs_swissprot"
QUERY_DB="$OUTDIR/proteome50_query_db"
RESULT_DB="$OUTDIR/proteome50_vs_swissprot_result"
OUT_TSV="$OUTDIR/proteome50_vs_swissprot.tsv"
TMP="$OUTDIR/tmp"

mkdir -p "$OUTDIR" "$TMP"

START=$(date +%s)
echo "[$(date '+%H:%M:%S')] Starting proteome50 vs SwissProt search"
echo "  Query FASTA : $QUERY_FASTA"
echo "  Target DB   : $SWISSPROT_DB"
echo "  Output TSV  : $OUT_TSV"
echo ""

echo "[1/3] Building query DB..."
mmseqs createdb "$QUERY_FASTA" "$QUERY_DB"

echo "[2/3] Searching against SwissProt (this is the slow step)..."
mmseqs search "$QUERY_DB" "$SWISSPROT_DB" "$RESULT_DB" "$TMP" \
    -e 1e-3 \
    --threads 8 \
    -v 2

echo "[3/3] Converting results to TSV..."
mmseqs convertalis "$QUERY_DB" "$SWISSPROT_DB" "$RESULT_DB" "$OUT_TSV" \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

END=$(date +%s)
ELAPSED=$(( END - START ))
echo ""
echo "[$(date '+%H:%M:%S')] Done. Elapsed: ${ELAPSED}s ($(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s)"
echo "Output: $OUT_TSV"
echo "Rows:   $(wc -l < "$OUT_TSV")"
