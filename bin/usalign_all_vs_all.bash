#!/bin/bash
set -euo pipefail

DIR="${1:?Usage: $0 <pdb_dir> [output_file] [chains] [threads]}"
DIR="${DIR%/}"
SUFFIX=".pdb"
CHAINS=${3:-"A,B"}
THREADS="${4:-$(nproc)}"

TMPLIST=$(mktemp /tmp/usalign_list.XXXXXX)
trap 'rm -f "$TMPLIST"' EXIT

# List basenames without suffix
for f in "$DIR"/*"$SUFFIX"; do
  basename "$f" "$SUFFIX"
done > "$TMPLIST"

N=$(wc -l < "$TMPLIST")
echo "Found $N structures → $(( N*(N-1)/2 )) pairs" >&2

header="#PDB1\tPDB2\tTM\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali"

run_alignments() {
  local total=$(( N*(N-1)/2 ))
  local counter=0
  awk '{a[NR]=$1} END{for(i=1;i<NR;i++) for(j=i+1;j<=NR;j++) print a[i],a[j]}' \
      "$TMPLIST" \
    | parallel --will-cite --colsep ' ' -j "$THREADS" \
      "USalign '${DIR}'/{1}${SUFFIX} '${DIR}'/{2}${SUFFIX} -chain1 ${CHAINS} -chain2 ${CHAINS} -mm 1 -ter 1 -outfmt 2  2>/dev/null | grep -v '^#'" \
    | while IFS= read -r line; do
        counter=$((counter + 1))
        printf '\r%d / %d completed' "$counter" "$total" >&2
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
  {
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
  echo "Done → $2" >&2
else
  run_alignments | add_tm_col
fi


