#!/usr/bin/env bash
# Fetch chr20 smoke data from SiYangming/Ribo-seq (not stored in this repo).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-$HERE/resources}"
mkdir -p "$DEST"
URL="https://github.com/SiYangming/Ribo-seq.git"
TMP="$(mktemp -d)"
echo "[INFO] sparse-checkout test/ from $URL"
git clone --depth 1 --filter=blob:none --sparse "$URL" "$TMP/Ribo-seq"
git -C "$TMP/Ribo-seq" sparse-checkout set test
rsync -a "$TMP/Ribo-seq/test/" "$DEST/"
rm -rf "$TMP"
echo "[INFO] testdata → $DEST"
echo "[INFO] point sample_sheet / fasta_dir at this tree (see test/README in Ribo-seq)."
