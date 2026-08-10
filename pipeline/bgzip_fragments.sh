#!/usr/bin/env bash
# ============================================================================
# bgzip_fragments.sh
# ----------------------------------------------------------------------------
# Compress the per-cell-type fragment BEDs in place. They are ~125 GB raw and
# compress ~4x, which matters twice over: disk here, and the transfer to a GPU
# host, which is the real friction in getting ChromBPNet trained.
#
# bgzip rather than gzip: BGZF is block-compressed, so the output stays
# tabix-indexable and seekable. ChromBPNet reads bgzipped fragments for -ifrag
# (its own bias step passes a .tsv.gz), and MACS2 reads .gz natively, so
# nothing downstream needs the plain text back.
#
# pipeline/call_celltype_peaks.sh and pipeline/run_chrombpnet.sh both resolve
# <name>.bed.gz in preference to <name>.bed and fall back to plain, so the
# pipeline works before or after this runs. Do not compress these BEDs without
# those two scripts in place: their existence guards would skip every cell type
# with only a log line, and `sort` would silently sort compressed bytes.
#
# SAFETY. An original is deleted only after all three of:
#   1. bgzip exits 0
#   2. bgzip -t passes (BGZF structure + CRC, so truncation is caught)
#   3. the decompressed line count equals the original's
# Any failure keeps the original and removes the partial .gz. Nothing here is
# irreplaceable in any case -- `Rscript pipeline/R/prepare_ml_input.R` rebuilds
# all of it in ~30 min -- but a half-written file that looks complete would be
# worse than an obvious failure.
#
# Idempotent: a file whose .gz already exists and passes bgzip -t is skipped,
# and any leftover plain copy removed. Safe to re-run after an interruption.
#
# Run from the repo root:
#   bash pipeline/bgzip_fragments.sh
#   BGZIP_THREADS=4 bash pipeline/bgzip_fragments.sh
#
# This takes ~50 min over 125 GB. This Mac idle-sleeps and has killed a long
# run before, so launch it detached under caffeinate:
#   caffeinate -i -m -s bash pipeline/bgzip_fragments.sh > logs/bgzip.log 2>&1 &
# ============================================================================
# Deliberately NOT `set -e`: a single unreadable file should be logged and
# skipped, not abort a 50-minute run over the other 40.
set -uo pipefail

THREADS=${BGZIP_THREADS:-10}
DIRS=${DIRS:-"ML/celltype_fragments ML/celltype_fragments_injury"}

log(){ echo "$(date '+%H:%M:%S') | $*"; }

command -v bgzip >/dev/null || { echo "bgzip not found (htslib). Install it or run in the container." >&2; exit 1; }

# stat is not portable: BSD/macOS takes -f%z, GNU/Linux -c%s. This repo runs
# both on this Mac and in an amd64 Linux container, so pick once up front.
if stat -f%z . >/dev/null 2>&1; then fsize(){ stat -f%z "$1"; }
else                                 fsize(){ stat -c%s "$1"; }; fi

gb(){ echo "scale=1; $1/1073741824" | bc; }

tot_before=0; tot_after=0; failed=0; done_n=0; skipped=0

for d in $DIRS; do
  [ -d "$d" ] || { log "no such directory: $d -- skipping"; continue; }
  for f in "$d"/*.bed; do
    [ -e "$f" ] || continue          # unmatched glob when a dir is fully compressed
    gz="$f.gz"

    if [ -s "$gz" ] && bgzip -t "$gz" 2>/dev/null; then
      log "SKIP $(basename "$f") (.gz already present and valid)"
      [ -e "$f" ] && { log "  removing leftover plain $(basename "$f")"; rm -f "$f"; }
      skipped=$((skipped+1)); continue
    fi

    b=$(fsize "$f")
    log "COMPRESS $(basename "$f") ($(gb "$b") GB)"
    n_orig=$(wc -l < "$f" | tr -d ' ')

    if ! bgzip -@ "$THREADS" -c "$f" > "$gz"; then
      log "  !! bgzip FAILED, keeping original"; rm -f "$gz"; failed=$((failed+1)); continue
    fi
    if ! bgzip -t "$gz" 2>/dev/null; then
      log "  !! integrity check FAILED, keeping original"; rm -f "$gz"; failed=$((failed+1)); continue
    fi
    n_gz=$(bgzip -dc "$gz" | wc -l | tr -d ' ')
    if [ "$n_orig" != "$n_gz" ]; then
      log "  !! LINE MISMATCH $n_orig vs $n_gz, keeping original"
      rm -f "$gz"; failed=$((failed+1)); continue
    fi

    a=$(fsize "$gz")
    rm -f "$f"
    tot_before=$((tot_before+b)); tot_after=$((tot_after+a)); done_n=$((done_n+1))
    log "  ok: $n_gz lines, $(echo "scale=2;$b/$a"|bc)x -> $(gb "$a") GB"
  done
done

log "================================================"
log "compressed $done_n, skipped $skipped, failed $failed"
[ "$tot_before" -gt 0 ] && log "total $(gb "$tot_before") GB -> $(gb "$tot_after") GB"
[ "$failed" -gt 0 ] && { log "DONE WITH FAILURES"; exit 1; }
log "DONE"
