#!/bin/bash
set -euo pipefail

# Align PDB structures within each cluster.
#
# Input format (wide):
#   Each line = one cluster, members space-separated
#   Each member = bare basename (e.g. run004_final)
#
# Usage:
#   align_clusters.bash <clusters.dat> <output_dir> <pdb_dir> [min_members] [chains]
#
# USalign superimposes each member onto the cluster representative (first entry).

CLUSTERS="${1:?Usage: $0 <clusters.dat> <output_dir> <pdb_dir> [min_members] [chains]}"
OUTDIR="${2:?Provide output directory}"
PDBDIR="${3:?Provide PDB directory}"
PDBDIR="${PDBDIR%/}"
MIN_MEMBERS="${4:-3}"
CHAINS="${5:-}"   # empty = align the first chain of each structure (USalign default)

if [[ -n "$CHAINS" ]]; then
  USALIGN_OPTS="-chain1 ${CHAINS} -chain2 ${CHAINS} -mm 1 -ter 1"
else
  USALIGN_OPTS=""
fi

if [ -d "$OUTDIR" ]; then
  rm -rf "$OUTDIR"
fi

mkdir -p "$OUTDIR"

# Structures in clusters below MIN_MEMBERS go here instead of being aligned.
UNCLUST="$OUTDIR/unclust"
mkdir -p "$UNCLUST"

# Resolve bare basename to PDB file path
resolve_pdb() {
  local name="$1"
  local pdb="$PDBDIR/${name}.pdb"
  if [[ -f "$pdb" ]]; then
    echo "$pdb"
  else
    echo ""
  fi
}

clust_idx=0
unclust_n=0
while IFS= read -r line; do
  [[ -z "$line" ]] && continue

  read -r -a members <<< "$line"

  # Too small to align: copy every member into unclust/.
  if [[ ${#members[@]} -lt $MIN_MEMBERS ]]; then
    for m in "${members[@]}"; do
      pdb="$(resolve_pdb "$m")"
      if [[ -n "$pdb" ]]; then
        cp "$pdb" "$UNCLUST/${m}.pdb"
        unclust_n=$((unclust_n + 1))
      fi
    done
    continue
  fi

  clust_idx=$((clust_idx + 1))
  cdir="$OUTDIR/clust$(printf '%04d' $clust_idx)"
  mkdir -p "$cdir"

  ref_pdb="$(resolve_pdb "${members[0]}")"
  if [[ -z "$ref_pdb" ]]; then
    echo "WARNING: reference PDB not found: ${members[0]} — skipping cluster $clust_idx" >&2
    continue
  fi

  cp "$ref_pdb" "$cdir/${members[0]}_ref.pdb"

  for ((i=1; i<${#members[@]}; i++)); do
    target_pdb="$(resolve_pdb "${members[$i]}")"
    if [[ -z "$target_pdb" ]]; then
      echo "WARNING: target PDB not found: ${members[$i]} — skipping" >&2
      continue
    fi

    out_prefix="$cdir/${members[$i]}_sup"

    USalign "$target_pdb" "$ref_pdb" ${USALIGN_OPTS} \
            -o "$out_prefix" >/dev/null || \
      echo "WARNING: USalign failed for ${members[$i]}" >&2
  done

  rm -f "$cdir"/*.pml

  echo "clust$(printf '%04d' $clust_idx): ${#members[@]} members, rep=${members[0]}" >&2
done < "$CLUSTERS"

# Drop the unclust/ dir if nothing landed there.
rmdir "$UNCLUST" 2>/dev/null || true

echo "Done: $clust_idx clusters written to $OUTDIR ($unclust_n unclustered → ${UNCLUST})" >&2
