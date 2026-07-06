#!/bin/bash
set -euo pipefail

DIR="${1:?Usage: $0 <pdb_dir> [output_file] [chains] [threads]}"
DIR="${DIR%/}"
SUFFIX=".pdb"
CHAINS="${3:-}"        # empty = let USalign use the first chain of each structure
THREADS="${4:-100%}"   # GNU parallel syntax: 100% = one job per CPU core

# With explicit chains, align those chains as a complex (-mm 1); otherwise leave
# USalign on its single-chain default, which aligns the first chain of each file.
if [[ -n "$CHAINS" ]]; then
  USALIGN_OPTS="-chain1 ${CHAINS} -chain2 ${CHAINS} -mm 1 -ter 1"
else
  USALIGN_OPTS=""
fi

command -v USalign  >/dev/null 2>&1 || { echo "ERROR: USalign not found on PATH"        >&2; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel not found on PATH"    >&2; exit 1; }

TMPLIST=$(mktemp /tmp/usalign_list.XXXXXX)
trap 'rm -f "$TMPLIST"' EXIT

# List basenames without suffix
for f in "$DIR"/*"$SUFFIX"; do
  basename "$f" "$SUFFIX"
done > "$TMPLIST"

N=$(wc -l < "$TMPLIST")
TOTAL=$(( N * (N - 1) / 2 ))
echo "Found $N structures → $TOTAL pairs (threads: $THREADS, chains: ${CHAINS:-first})" >&2

header="#PDB1\tPDB2\tTM\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali"

# Run USalign on every unordered pair in parallel.
#
# Each job filters its own output with `awk 'NF>=11 && $1 !~ /^#/'` so that only
# well-formed -outfmt 2 rows (11 tab-separated columns) survive: USalign prints
# warnings/errors (e.g.
# "ERROR! 0 chain in complex 2") to *stdout*, and those have no tab separators,
# so they are dropped instead of leaking into the results as garbage rows.
# This also keeps every job's exit status at 0, so a single unparseable pair
# cannot trip `set -e`/`pipefail` and abort the whole run.
run_alignments() {
  local counter=0
  awk '{a[NR]=$1} END{for(i=1;i<NR;i++) for(j=i+1;j<=NR;j++) print a[i],a[j]}' \
      "$TMPLIST" \
    | parallel --will-cite --colsep ' ' -j "$THREADS" \
      "USalign '${DIR}'/{1}${SUFFIX} '${DIR}'/{2}${SUFFIX} ${USALIGN_OPTS} -outfmt 2 2>/dev/null | awk -F'\t' 'NF>=11 && \$1 !~ /^#/'" \
    | while IFS= read -r line; do
        counter=$((counter + 1))
        printf '\r%d / %d pairs aligned' "$counter" "$TOTAL" >&2
        printf '%s\n' "$line"
      done
  printf '\n' >&2
}

# Strip chain suffixes from PDB names and add TM=max(TM1,TM2) column.
# USalign -outfmt 2 produces entries like path/file.pdb:A:B — we keep only path/file.pdb.
# When multiple chain combinations exist for the same PDB pair, keep only the best TM.
add_tm_col() {
  echo -e "$header"
  awk -F'\t' 'BEGIN{OFS="\t"}
  NF >= 11 {
    # Strip path, extension, and chain suffixes → bare basename
    pdb1=$1; pdb2=$2
    sub(/.*\//, "", pdb1); sub(/\.pdb.*/, "", pdb1)
    sub(/.*\//, "", pdb2); sub(/\.pdb.*/, "", pdb2)

    # Compute TM = max(TM1, TM2)
    tm=($3>$4)?$3:$4

    # Build canonical pair key (sorted)
    if (pdb1 < pdb2) key = pdb1 SUBSEP pdb2
    else             key = pdb2 SUBSEP pdb1

    # Keep only the best TM per pair
    if (!(key in best) || tm > best[key]) {
      best[key] = tm
      line[key] = pdb1 OFS pdb2 OFS tm OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
    }
  }
  END {
    for (k in line) print line[k]
  }'
}

if [[ -n "${2:-}" ]]; then
  run_alignments | add_tm_col > "$2"
  aligned=$(grep -vc '^#' "$2")
  echo "Done → $2 ($aligned/$TOTAL pairs)" >&2
else
  out=$(run_alignments | add_tm_col)
  printf '%s\n' "$out"
  aligned=$(printf '%s\n' "$out" | grep -vc '^#')
fi

if (( aligned < TOTAL )); then
  echo "WARNING: $((TOTAL - aligned)) of $TOTAL pairs produced no alignment — " \
       "unparseable structures or no matching '$CHAINS' chains" >&2
fi
