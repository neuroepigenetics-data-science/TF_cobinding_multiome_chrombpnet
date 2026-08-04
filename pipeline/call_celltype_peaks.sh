#!/usr/bin/env bash
# ============================================================================
# call_celltype_peaks.sh
# ----------------------------------------------------------------------------
# Track A, step 1 of run_chrombpnet.sh: per-cell-type MACS2 peaks from the
# fragment BEDs, blacklist-filtered, plus the DAR x peak intersect that
# run_chrombpnet.sh expects but neither published script produces.
#
# Config-driven version of code/machine_learning/run_chrombpnet.sh lines 40-62.
#
# Per cell type:
#   1. drop non-canonical contigs (GL*/JH* scaffolds)  -> <ct>_chr.bed
#   2. macs2 callpeak                                  -> <ct>_fragments_peaks.narrowPeak
#   3. subtract the slopped blacklist                  -> <ct>_peaks_no_blacklist[_chr].bed
#   4. intersect the injury DARs with those peaks      -> <ct>_intersect_dars.bed
#
# Deliberate departures from run_chrombpnet.sh:
#   1. `-B` (bedGraph pileup) is OFF by default. It does not affect the peak
#      calls at all -- it only writes extra genome-wide pileup tracks, several
#      GB per cell type. Nothing downstream in ChromBPNet reads them. Set
#      MAKE_BEDGRAPH=1 to restore the published behaviour.
#   2. `<ct>_chr.bed` is written to ML/celltype_fragments/, not
#      celltype_fragments_multiome_only/. The published script writes the former
#      path but its training step reads `ML/celltype_fragments/${i}_chr.bed`
#      (line 67) -- the two never agree. We use the path training reads.
#   3. Step 4 does not exist in either published script. run_chrombpnet.sh reads
#      dars_intersect_peaks/<ct>_intersect_dars.bed for contribution scoring;
#      this produces it. `bedtools intersect -u` keeps the WHOLE DAR region
#      (uniform 500 bp) when it overlaps a called peak, rather than the clipped
#      sub-interval, because ChromBPNet expects fixed-width inputs.
#   4. Idempotent for the expensive steps: existing <ct>_chr.bed and MACS2
#      narrowPeak outputs are skipped, so an interrupted run resumes instead of
#      redoing MACS2. Steps 3 and 4 are cheap and always re-run.
#   5. Step 4 emits 10-column narrowPeak, not 3-column BED. ChromBPNet's
#      contribs_bw/pred_bw read -r with NARROWPEAK_SCHEMA and centre each
#      sequence at start(col2)+summit(col10); a 3-column file cannot be parsed.
#      Verified against chrombpnet/evaluation/interpret/interpret.py.
#
# Usage (from the repo root):
#   bash pipeline/call_celltype_peaks.sh                 # all 10 cell types
#   bash pipeline/call_celltype_peaks.sh Endothelial     # named subset
#   MAKE_BEDGRAPH=1 bash pipeline/call_celltype_peaks.sh # with pileup tracks
# ============================================================================
set -euo pipefail

# Tool locations, in precedence order:
#   1. MACS2 / BEDTOOLS environment variables
#   2. pipeline/tools.local.env, if present (host-specific, gitignored) --
#      needed on this Mac because macs2 and bedtools live in their own conda
#      envs and are therefore not on PATH
#   3. whatever is on PATH (the container case: both are installed into it)
# Hardcoded absolute paths were removed: they made this script unrunnable
# anywhere but the machine it was written on.
if [ -f pipeline/tools.local.env ]; then
  # shellcheck disable=SC1091
  . pipeline/tools.local.env
fi
# Resolve to something actually executable. Checking only for "non-empty" is not
# enough: a stale tools.local.env (or an inherited variable) can name a path that
# does not exist on this machine -- exactly what happens if the host file is
# baked into a container -- and we would then fail deep inside the run instead
# of here. So validate, and fall back to PATH when the named path is unusable.
resolve_tool() {
  local given="$1" name="$2"
  if [ -n "$given" ] && [ -x "$given" ]; then echo "$given"; return 0; fi
  if [ -n "$given" ] && [ ! -x "$given" ]; then
    echo "note: $name=$given is not executable here; falling back to PATH" >&2
  fi
  command -v "$name" 2>/dev/null || true
}
MACS2=$(resolve_tool "${MACS2:-}" macs2)
BEDTOOLS=$(resolve_tool "${BEDTOOLS:-}" bedtools)
if [ -z "$MACS2" ] || [ -z "$BEDTOOLS" ]; then
  echo "macs2 and/or bedtools not found." >&2
  echo "  Set MACS2=/path/to/macs2 BEDTOOLS=/path/to/bedtools," >&2
  echo "  or create pipeline/tools.local.env exporting them," >&2
  echo "  or run inside the container where both are on PATH." >&2
  exit 1
fi
echo "using macs2:    $MACS2"
echo "using bedtools: $BEDTOOLS"

FRAG_DIR=ML/celltype_fragments
DAR_DIR=ML/celltype_dars
PEAK_DIR=ML/macs2_celltype_peaks
FILT_DIR=ML/filtered_peaks
ISEC_DIR=ML/dars_intersect_peaks

CHROM_SIZES=meta/mm10.chr.sizes
BLACKLIST=meta/mm10-blacklist.v2.bed
GENOME_SIZE=1.87e9
SLOP=1057
MAKE_BEDGRAPH=${MAKE_BEDGRAPH:-0}

ALL_CELLTYPES=(Astrocytes Ependymal Microglia Oligodendrocytes Macrophages
               Neurons_V Perivascular Endothelial OPCs Neurons_D)
CELLTYPES=("${@:-${ALL_CELLTYPES[@]}}")
if [ "$#" -gt 0 ]; then CELLTYPES=("$@"); fi

log() { echo "$(date '+%H:%M:%S') | $*"; }

for f in "$MACS2" "$BEDTOOLS" "$CHROM_SIZES" "$BLACKLIST"; do
  [ -e "$f" ] || { echo "missing required file: $f" >&2; exit 1; }
done
mkdir -p "$PEAK_DIR" "$FILT_DIR" "$ISEC_DIR"

# --- blacklist, slopped once ------------------------------------------------
SLOPPED="$FILT_DIR/blacklist_slop${SLOP}.bed"
if [ ! -s "$SLOPPED" ]; then
  log "slopping blacklist by ${SLOP}bp -> $SLOPPED"
  "$BEDTOOLS" slop -i "$BLACKLIST" -g "$CHROM_SIZES" -b "$SLOP" > "$SLOPPED"
fi
log "blacklist regions: $(wc -l < "$SLOPPED")"

for CT in "${CELLTYPES[@]}"; do
  log "================ $CT ================"
  FRAG="$FRAG_DIR/${CT}.bed"
  CHRBED="$FRAG_DIR/${CT}_chr.bed"
  [ -s "$FRAG" ] || { log "  MISSING $FRAG -- skipping"; continue; }

  # 1. canonical contigs only
  if [ -s "$CHRBED" ]; then
    log "  1. $CHRBED exists, skipping filter"
  else
    log "  1. filtering to ^chr ..."
    grep '^chr' "$FRAG" > "$CHRBED"
    log "     $(wc -l < "$FRAG") -> $(wc -l < "$CHRBED") fragments"
  fi

  # 2. MACS2
  NARROW="$PEAK_DIR/${CT}_fragments_peaks.narrowPeak"
  if [ -s "$NARROW" ]; then
    log "  2. $NARROW exists, skipping macs2"
  else
    log "  2. macs2 callpeak ..."
    # plain string, not an array: macOS bash 3.2 errors on empty array
    # expansion under `set -u`. Unquoted and empty, this expands to no word.
    BDG=""
    [ "$MAKE_BEDGRAPH" = "1" ] && BDG="-B"
    # shellcheck disable=SC2086
    "$MACS2" callpeak -f AUTO -t "$CHRBED" \
      -g "$GENOME_SIZE" $BDG -p 0.01 --nomodel --extsize 200 \
      --outdir "$PEAK_DIR" -n "${CT}_fragments"
    log "     peaks called: $(wc -l < "$NARROW")"
  fi

  # 3. blacklist subtraction
  NOBL="$FILT_DIR/${CT}_peaks_no_blacklist.bed"
  NOBL_CHR="$FILT_DIR/${CT}_peaks_no_blacklist_chr.bed"
  log "  3. subtracting blacklist ..."
  "$BEDTOOLS" intersect -v -a "$NARROW" -b "$SLOPPED" > "$NOBL"
  grep '^chr' "$NOBL" > "$NOBL_CHR" || true
  log "     $(wc -l < "$NARROW") -> $(wc -l < "$NOBL_CHR") peaks after blacklist"

  # 4. DAR x peak intersect (the gap; see header notes 3 and 5)
  DARS="$DAR_DIR/${CT}_I_dars.bed"
  ISEC="$ISEC_DIR/${CT}_intersect_dars.bed"
  if [ -s "$DARS" ]; then
    # Emitted as 10-column narrowPeak, NOT 3-column BED: chrombpnet's
    # contribs_bw/pred_bw parse -r with
    #   NARROWPEAK_SCHEMA = [chr,start,end,1,2,3,4,5,6,summit]
    # and extract each sequence centred at start(col2)+summit(col10). A 3-column
    # file has no col10 and the parse fails. Columns 4-9 are never read, so they
    # carry placeholders (-1 = "no data", the narrowPeak convention).
    # summit is derived from the region width rather than hardcoded to 250, so
    # this stays correct if the DAR width ever changes; for these uniform 500 bp
    # regions it yields 250, i.e. the region centre -- the same point verified
    # to sit >1057 bp from every chromosome edge.
    "$BEDTOOLS" intersect -u -a "$DARS" -b "$NOBL_CHR" \
      | awk -v ct="$CT" 'BEGIN{OFS="\t"}
                         {print $1,$2,$3,ct"_dar_"NR,0,".",-1,-1,-1,int(($3-$2)/2)}' \
      > "$ISEC"
    log "  4. DARs $(wc -l < "$DARS") -> $(wc -l < "$ISEC") overlapping a called peak (10-col narrowPeak)"
  else
    log "  4. no DAR file ($DARS) -- skipping intersect"
  fi
done

log "done"
