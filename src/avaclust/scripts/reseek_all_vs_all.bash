#!/bin/bash
# All-vs-all protein structure alignment with reseek (https://github.com/rcedgar/reseek).
#
# Standalone for now; intended to be wired into the avaclust CLI later as an
# alternative to usalign_all_vs_all.bash.
#
# Scores are reseek p-values, NOT TM-scores: lower = more similar (~1e-12 for
# near-identical, up to 1 for unrelated). Anything consuming the output must use
# a p-value-aware distance instead of 1 - TM.
set -euo pipefail

DIR="${1:?Usage: $0 <struct_dir> [output_file] [stats_model] [threads]}"
DIR="${DIR%/}"
MODEL="${3:-superfamily}"     # -stats model for p-values: family, superfamily or fold
THREADS="${4:-$(nproc)}"

case "$MODEL" in
  family|fam|superfamily|sf|fold) ;;
  *) echo "ERROR: stats model must be family, superfamily or fold (got '$MODEL')" >&2; exit 1 ;;
esac

command -v reseek >/dev/null 2>&1 || { echo "ERROR: reseek not found on PATH" >&2; exit 1; }
[[ -d "$DIR" ]] || { echo "ERROR: $DIR is not a directory" >&2; exit 1; }

TMPDIR_RS=$(mktemp -d /tmp/reseek_ava.XXXXXX)
trap 'rm -rf "$TMPDIR_RS"' EXIT
DB="$TMPDIR_RS/db.bcb"
HITS="$TMPDIR_RS/hits.tsv"
LOG="$TMPDIR_RS/reseek.log"

# reseek is chatty on stderr; keep its log and only show it if a step fails.
run_reseek() {
  if ! reseek "$@" >>"$LOG" 2>&1; then
    echo "ERROR: reseek $1 failed:" >&2
    tail -n 20 "$LOG" >&2
    exit 1
  fi
}

# Input structures, named the way reseek labels them (basename minus extension).
# reseek scans the directory recursively for these formats itself.
find "$DIR" -type f ! -name '._*' \
     \( -name '*.pdb' -o -name '*.ent' -o -name '*.cif' -o -name '*.mmcif' \
        -o -name '*.pdb.gz' -o -name '*.ent.gz' -o -name '*.cif.gz' -o -name '*.mmcif.gz' \) \
  | sed -E 's#.*/##; s/\.gz$//; s/\.[^.]+$//' | sort -u > "$TMPDIR_RS/inputs"
N_IN=$(wc -l < "$TMPDIR_RS/inputs")
(( N_IN > 0 )) || { echo "ERROR: no .pdb/.ent/.cif structures found in $DIR" >&2; exit 1; }

# 1. Build the database. Every protein chain becomes an entry labelled <name>_<chain>;
#    non-protein chains, unreadable files and chains under reseek's minimum length
#    are silently skipped. Multithreaded -convert intermittently segfaults in
#    reseek v3.0, so run it single-threaded (it is fast relative to the search).
echo "Converting $N_IN structures → reseek DB" >&2
run_reseek -convert "$DIR" -bcb "$DB" -threads 1
run_reseek -convert "$DB" -fasta "$TMPDIR_RS/db.fa" -threads 1
grep '^>' "$TMPDIR_RS/db.fa" | sed -E 's/^>//; s/_[^_]*$//' | sort -u > "$TMPDIR_RS/names"
N=$(wc -l < "$TMPDIR_RS/names")

dropped=$(comm -23 "$TMPDIR_RS/inputs" "$TMPDIR_RS/names")
if [[ -n "$dropped" ]]; then
  echo "WARNING: $(wc -l <<<"$dropped") structure(s) have no usable protein chain" \
       "(unreadable, non-protein, or too short) and are excluded:" >&2
  sed 's/^/  /' <<<"$dropped" >&2
fi
(( N >= 2 )) || { echo "ERROR: need at least 2 usable structures, got $N" >&2; exit 1; }

TOTAL=$(( N * (N - 1) / 2 ))
echo "Searching $N structures → $TOTAL pairs (stats: $MODEL, threads: $THREADS)" >&2

# 2. Search the DB against itself. -pvalue 1 disables the reporting cutoff so we
#    keep weak hits too; reseek still omits pairs with no detectable similarity.
run_reseek -search "$DB" -db "$DB" -output "$HITS" -sensitive -stats "$MODEL" \
           -pvalue 1 -threads "$THREADS" \
           -columns query+target+pvalue+pctid+ql+tl+qlo+qhi+tlo+thi

# 3. Collapse to one row per unordered structure pair, keeping the best (lowest)
#    p-value over both search directions and all chain combinations. Pairs reseek
#    did not report are written with PVALUE=1 and NA for the alignment fields, so
#    every structure and every pair is present in the output.
collapse_pairs() {
  printf '#PDB1\tPDB2\tPVALUE\tCH1\tCH2\tPCTID\tL1\tL2\tQLO\tQHI\tTLO\tTHI\n'
  awk -F'\t' 'BEGIN{OFS="\t"}
  NR == FNR { names[++n] = $1; next }
  {
    q = $1; t = $2
    qc = q; tc = t
    sub(/_[^_]*$/, "", q); sub(/.*_/, "", qc)
    sub(/_[^_]*$/, "", t); sub(/.*_/, "", tc)
    if (q == t) next
    p = $3 + 0
    # Store in canonical (sorted) orientation so both directions share a key
    if (q < t) { key = q SUBSEP t; rest = qc OFS tc OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 }
    else       { key = t SUBSEP q; rest = tc OFS qc OFS $4 OFS $6 OFS $5 OFS $9 OFS $10 OFS $7 OFS $8 }
    if (!(key in best) || p < best[key]) { best[key] = p; line[key] = rest }
  }
  END {
    for (i = 1; i < n; i++)
      for (j = i + 1; j <= n; j++) {
        key = names[i] SUBSEP names[j]
        if (key in best) { print names[i], names[j], best[key], line[key]; hit++ }
        else print names[i], names[j], 1, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
      }
    printf "%d / %d pairs with a reseek hit; the rest set to PVALUE=1\n", hit, n * (n - 1) / 2 > "/dev/stderr"
  }' "$TMPDIR_RS/names" "$HITS"
}

if [[ -n "${2:-}" ]]; then
  collapse_pairs > "$2"
  echo "Done → $2" >&2
else
  collapse_pairs
fi
