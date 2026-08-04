#!/usr/bin/env bash
# ============================================================================
# run_chrombpnet.sh  (corrected, config-driven)
# ----------------------------------------------------------------------------
# Track A, steps 2-6: bias model -> per-cell-type ChromBPNet models ->
# contribution scores over the injury DARs -> MoDISco motifs.
#
# Corrected rewrite of code/machine_learning/run_chrombpnet.sh. That script is
# not runnable as published; every departure below is a fix for a specific
# defect in it, verified against the chrombpnet source/wiki (see FIXES).
#
# Steps 1 (per-cell-type MACS2 peaks + DAR intersect) is pipeline/call_celltype_peaks.sh.
# Run that first. This script consumes its outputs.
#
# ---------------------------------------------------------------------------
# FIXES vs the published script (published line numbers in brackets)
# ---------------------------------------------------------------------------
#  1. [92,97,99,100] Literal `{i}` instead of `${i}` in every interpret/MoDISco
#     output path. All 10 cell types wrote into one directory literally named
#     "{i}", each silently overwriting the last -- you would finish a full GPU
#     run holding a single cell type's results and no error anywhere.
#  2. [88,89] Model paths read `models/${i}_fl0/models/...` but training [75]
#     writes to `-o models/${i}`. Per the wiki the layout is
#     <outdir>/models/{chrombpnet.h5,chrombpnet_nobias.h5,bias_model_scaled.h5};
#     the `_fl0` suffix does not exist.
#  3. [94] Contribution scores were computed with `chrombpnet.h5`. The wiki is
#     explicit: "Please use the chrombpnet_nobias.h5 model for this."
#     chrombpnet.h5 still carries the Tn5 bias component, so MoDISco run on its
#     contributions rediscovers the Tn5 insertion motif rather than TF motifs.
#     This is a silent scientific error, not a crash.
#  4. [74] `-b small_frag_bias.h5` names a file the script never produces. The
#     bias step [27-37] writes <bias_outdir>/models/bias.h5 (-fp does not rename
#     the model file).
#  5. [13,17] macs2 writes macs2_from_frag/bias_fragments_peaks.narrowPeak, but
#     the blacklist step read `bias_fragments.narrowPeak` -- wrong directory and
#     missing the `_peaks` infix.
#  6. [16,21,24,31,96] Referenced `mm10.chrom.sizes`; only `mm10.chr.sizes`
#     exists. Mixed with [54,61,70,91] which use the correct name.
#  7. [72] Stray `A` in `-n negatives/A${i}_nonpeaks_negatives.bed`.
#  8. [44 vs 67] `<ct>_chr.bed` written to celltype_fragments_multiome_only/ but
#     read from ML/celltype_fragments/. Handled in call_celltype_peaks.sh.
#  9. [61 vs 71] nonpeaks built from the `_chr` peak file, training from the
#     non-`_chr` one. Both are used consistently here. (They are identical for
#     this dataset -- 0 non-canonical contigs survive MACS2 -- but relying on
#     that is fragile.)
# 10. [86 vs 90,95] Created `<ct>_intersect_dars_chr.bed` then scored the
#     non-`_chr` file. Dead code; dropped.
#
# ---------------------------------------------------------------------------
# NOTES ON THE DATA (carry these into interpretation)
# ---------------------------------------------------------------------------
#  * MoDISco is SKIPPED for cell types with fewer than MODISCO_MIN_REGIONS
#    regions. Macrophages has 273 (its uninjured group is 22 cells) -- far too
#    few for seqlet clustering. Its models still train; only motif calling is
#    skipped, and the skip is logged rather than silently succeeding on noise.
#  * Endothelial keeps only 28% of its DARs after the peak intersect (thin
#    contrasts: 9-cell 3dpi, 30-cell 7dpi). Treat its motifs cautiously.
#  * The fold JSONs cover chr1-19 only -- chrX is in neither train, valid nor
#    test. ~1-2% of the DAR regions are on chrX and are scored by a model that
#    never saw that chromosome. Kept (published behaviour), but logged per cell
#    type so the number is visible. Set DROP_UNSPLIT_CHROMS=1 to exclude them.
#  * DARs use the published raw-p threshold (p_val<0.05, no multiple-testing
#    correction). See RUNNING_LOCALLY.md "DAR significance threshold" -- under
#    Bonferroni several cell types lose >95% of their regions.
#
# ---------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------
#   bash pipeline/run_chrombpnet.sh <stage> [celltype ...]
#
#   stages:  bias-peaks | splits | nonpeaks | bias-train | train | interpret
#            prep   = bias-peaks + splits + nonpeaks   (CPU only)
#            all    = everything
#
#   GPU is required for bias-train / train / interpret. `prep` is CPU-only and
#   can be done on the laptop (emulated amd64 container is fine).
#
#   env overrides: CHROMBPNET, MACS2, BEDTOOLS, THREADS, MODISCO_MIN_REGIONS,
#                  DROP_UNSPLIT_CHROMS, FOLD, DRY_RUN
# ============================================================================
set -euo pipefail

# --- tools ------------------------------------------------------------------
if [ -f pipeline/tools.local.env ]; then
  # shellcheck disable=SC1091
  . pipeline/tools.local.env
fi
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
CHROMBPNET=$(resolve_tool "${CHROMBPNET:-}" chrombpnet)
MODISCO=$(resolve_tool "${MODISCO:-}" modisco)

# --- paths ------------------------------------------------------------------
FRAG_DIR=ML/celltype_fragments
FILT_DIR=ML/filtered_peaks
ISEC_DIR=ML/dars_intersect_peaks
BIAS_PEAK_DIR=ML/macs2_bias_peaks
NEG_DIR=ML/negatives
MODEL_DIR=ML/models
BIAS_DIR=ML/bias
INTERP_DIR=ML/interpret
LOG_DIR=logs/chrombpnet

GENOME=meta/mm10.fa
CHROM_SIZES=meta/mm10.chr.sizes
BLACKLIST=meta/mm10-blacklist.v2.bed
MOTIFS=meta/motifs.meme.txt
SPLIT_DIR=meta/splits
FOLD=${FOLD:-0}
FOLD_JSON="$SPLIT_DIR/fold_${FOLD}.json"

BIAS_FRAG_RAW="$FRAG_DIR/1_input_bias.bed"
BIAS_FRAG_GZ="$FRAG_DIR/1_input_bias.sorted.tsv.gz"

GENOME_SIZE=1.87e9
SLOP=1057
BIAS_THRESHOLD=0.5              # 0.5 for ATAC, 0.8 for DNase
DATA_TYPE=ATAC
THREADS=${THREADS:-8}
MODISCO_MIN_REGIONS=${MODISCO_MIN_REGIONS:-1000}
MODISCO_MAX_SEQLETS=${MODISCO_MAX_SEQLETS:-1000000}
DROP_UNSPLIT_CHROMS=${DROP_UNSPLIT_CHROMS:-0}
DRY_RUN=${DRY_RUN:-0}
SORT_TMP=${SORT_TMP:-$(pwd)/.sorttmp}

ALL_CELLTYPES=(Astrocytes Ependymal Microglia Oligodendrocytes Macrophages
               Neurons_V Perivascular Endothelial OPCs Neurons_D)

log()  { echo "$(date '+%H:%M:%S') | $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
# Safe line count: a bare `wc -l < missing` is a shell redirection error that
# `|| echo 0` cannot catch, which breaks logging under DRY_RUN.
nlines() { if [ -s "$1" ]; then wc -l < "$1" | tr -d ' '; else echo 0; fi; }

# --- bgzipped fragment support ----------------------------------------------
# Fragment BEDs are stored bgzipped (~28 GB vs ~125 GB raw). ChromBPNet accepts
# bgzipped fragments for -ifrag, so training needs no change beyond resolving
# the path. `frag_path` prefers .gz and falls back to plain, so both layouts
# work; without it, the `[ -s "$frag" ]` guard in stage_train would silently
# skip every cell type once the BEDs were compressed.
frag_path() { if [ -s "$1.gz" ]; then echo "$1.gz"; elif [ -s "$1" ]; then echo "$1"; else echo ""; fi; }
frag_cat()  { case "$1" in *.gz) bgzip -dc "$1";; *) cat "$1";; esac; }
# Same decision as frag_cat, but emitted as a command string for `bash -c`
# pipelines (needed so DRY_RUN can print the whole pipeline verbatim).
frag_cat_cmd() { case "$1" in *.gz) echo "bgzip -dc '$1'";; *) echo "cat '$1'";; esac; }
run()  {
  if [ "$DRY_RUN" = "1" ]; then echo "  [dry-run] $*"; else "$@"; fi
}

need_tool() { [ -n "$1" ] || die "$2 not found. Set $3=/path/to/$2 or put it on PATH."; }

# --- args -------------------------------------------------------------------
STAGE=${1:-}
[ -n "$STAGE" ] || die "usage: bash pipeline/run_chrombpnet.sh <stage> [celltype ...]
  stages: bias-peaks | splits | nonpeaks | prep | bias-train | train | interpret | all"
shift || true
if [ "$#" -gt 0 ]; then CELLTYPES=("$@"); else CELLTYPES=("${ALL_CELLTYPES[@]}"); fi

mkdir -p "$LOG_DIR"

# Regions actually scored for a cell type: optionally restricted to chromosomes
# present in the fold, since a model never trained on a chromosome is
# extrapolating there.
regions_for() {
  # NB: separate `local` statements. Assignments on a single `local` line are
  # not guaranteed to see each other, which trips `set -u`.
  local ct="$1"
  local src="$ISEC_DIR/${ct}_intersect_dars.bed"
  [ -s "$src" ] || { echo ""; return 0; }
  if [ "$DROP_UNSPLIT_CHROMS" != "1" ]; then echo "$src"; return 0; fi
  local out="$ISEC_DIR/${ct}_intersect_dars_infold.bed"
  if [ ! -s "$out" ]; then
    python3 - "$FOLD_JSON" "$src" "$out" <<'PY'
import json, sys
fold, src, out = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(fold))
keep = set(d.get("train", []) + d.get("valid", []) + d.get("test", []))
n = k = 0
with open(src) as fh, open(out, "w") as oh:
    for line in fh:
        n += 1
        if line.split("\t", 1)[0] in keep:
            oh.write(line); k += 1
print(f"     in-fold regions: {k}/{n}", file=sys.stderr)
PY
  fi
  echo "$out"
}

# ============================================================================
# stage: bias-peaks -- sort/bgzip/tabix the bias fragments, call peaks on them
# ============================================================================
stage_bias_peaks() {
  need_tool "$MACS2" macs2 MACS2
  need_tool "$BEDTOOLS" bedtools BEDTOOLS
  local raw; raw=$(frag_path "$BIAS_FRAG_RAW")
  [ -n "$raw" ] || die "missing $BIAS_FRAG_RAW[.gz] (run pipeline/R/prepare_ml_input.R)"
  mkdir -p "$BIAS_PEAK_DIR" "$SORT_TMP"

  if [ -s "$BIAS_FRAG_GZ" ] && [ -s "${BIAS_FRAG_GZ}.tbi" ]; then
    log "bias fragments already sorted+indexed: $BIAS_FRAG_GZ"
  else
    # 6.0 GB / ~133M fragments. -T must point somewhere roomy; /tmp is not.
    # Streamed through frag_cat so a bgzipped source works: `sort` cannot read
    # .gz itself, and handing it one silently sorts the compressed bytes.
    log "sorting bias fragments from $raw (tmp=$SORT_TMP) ..."
    run bash -c "$(frag_cat_cmd "$raw") | LC_ALL=C sort -k1,1V -k2,2n -T '$SORT_TMP' \
                   --parallel='$THREADS' | bgzip -@ '$THREADS' > '$BIAS_FRAG_GZ'"
    log "tabix ..."
    run tabix -f -p bed "$BIAS_FRAG_GZ"
  fi

  local narrow="$BIAS_PEAK_DIR/bias_fragments_peaks.narrowPeak"
  if [ -s "$narrow" ]; then
    log "bias peaks exist, skipping macs2"
  else
    log "macs2 on bias fragments ..."
    run "$MACS2" callpeak -f AUTO -t "$BIAS_FRAG_GZ" \
      -g "$GENOME_SIZE" -p 0.01 --nomodel --extsize 200 \
      --outdir "$BIAS_PEAK_DIR" -n bias_fragments
  fi

  local slopped="$FILT_DIR/blacklist_slop${SLOP}.bed"
  mkdir -p "$FILT_DIR"
  [ -s "$slopped" ] || run "$BEDTOOLS" slop -i "$BLACKLIST" -g "$CHROM_SIZES" -b "$SLOP" > "$slopped"

  local biaspk="$FILT_DIR/bias_peaks_no_blacklist.bed"
  if [ -s "$biaspk" ]; then
    log "bias peaks already blacklist-filtered ($(nlines "$biaspk") peaks), skipping"
  else
    log "subtracting blacklist from bias peaks ..."
    run bash -c "'$BEDTOOLS' intersect -v -a '$narrow' -b '$slopped' \
                   | grep '^chr' > '$biaspk'"
    log "bias peaks after blacklist: $(nlines "$biaspk")"
  fi
}

# ============================================================================
# stage: splits
# ============================================================================
stage_splits() {
  need_tool "$CHROMBPNET" chrombpnet CHROMBPNET
  mkdir -p "$SPLIT_DIR"
  if [ -s "$FOLD_JSON" ]; then
    log "$FOLD_JSON exists, skipping prep splits"
  else
    log "chrombpnet prep splits -> $FOLD_JSON"
    run "$CHROMBPNET" prep splits -c "$CHROM_SIZES" \
      -tcr chr1 chr3 -vcr chr8 chr19 -op "$SPLIT_DIR/fold_${FOLD}"
  fi
}

# ============================================================================
# stage: nonpeaks -- GC-matched negatives for the bias model and each cell type
# ============================================================================
stage_nonpeaks() {
  need_tool "$CHROMBPNET" chrombpnet CHROMBPNET
  mkdir -p "$NEG_DIR"
  [ -s "$FOLD_JSON" ] || die "missing $FOLD_JSON (run stage: splits)"

  local biaspk="$FILT_DIR/bias_peaks_no_blacklist.bed"
  if [ -s "$biaspk" ] && [ ! -s "$NEG_DIR/bias_nonpeaks_negatives.bed" ]; then
    log "nonpeaks: bias"
    run "$CHROMBPNET" prep nonpeaks -g "$GENOME" -p "$biaspk" -c "$CHROM_SIZES" \
      -fl "$FOLD_JSON" -br "$BLACKLIST" -o "$NEG_DIR/bias_nonpeaks"
  fi

  for CT in "${CELLTYPES[@]}"; do
    local pk="$FILT_DIR/${CT}_peaks_no_blacklist_chr.bed"
    local out="$NEG_DIR/${CT}_nonpeaks"
    [ -s "$pk" ] || { log "  $CT: missing $pk -- skipping"; continue; }
    if [ -s "${out}_negatives.bed" ]; then
      log "  $CT: negatives exist, skipping"; continue
    fi
    log "nonpeaks: $CT"
    run "$CHROMBPNET" prep nonpeaks -g "$GENOME" -p "$pk" -c "$CHROM_SIZES" \
      -fl "$FOLD_JSON" -br "$BLACKLIST" -o "$out"
  done
}

# ============================================================================
# stage: bias-train  (GPU)
# ============================================================================
stage_bias_train() {
  need_tool "$CHROMBPNET" chrombpnet CHROMBPNET
  local pk="$FILT_DIR/bias_peaks_no_blacklist.bed"
  local neg="$NEG_DIR/bias_nonpeaks_negatives.bed"
  [ -s "$pk" ]  || die "missing $pk (stage: bias-peaks)"
  [ -s "$neg" ] || die "missing $neg (stage: nonpeaks)"
  [ -s "$BIAS_FRAG_GZ" ] || die "missing $BIAS_FRAG_GZ (stage: bias-peaks)"

  if [ -s "$BIAS_DIR/models/bias.h5" ]; then
    log "bias model exists ($BIAS_DIR/models/bias.h5), skipping"
    return 0
  fi
  mkdir -p "$BIAS_DIR"
  log "training bias model -> $BIAS_DIR/models/bias.h5"
  run "$CHROMBPNET" bias pipeline \
    -ifrag "$BIAS_FRAG_GZ" \
    -d "$DATA_TYPE" \
    -g "$GENOME" \
    -c "$CHROM_SIZES" \
    -p "$pk" \
    -n "$neg" \
    -fl "$FOLD_JSON" \
    -b "$BIAS_THRESHOLD" \
    -o "$BIAS_DIR" \
    -fp bias_model 2>&1 | tee "$LOG_DIR/bias_train.log"
}

# ============================================================================
# stage: train  (GPU)  -- one ChromBPNet model per cell type
# ============================================================================
stage_train() {
  need_tool "$CHROMBPNET" chrombpnet CHROMBPNET
  local bias_model="$BIAS_DIR/models/bias.h5"
  [ -s "$bias_model" ] || die "missing $bias_model (stage: bias-train)"

  for CT in "${CELLTYPES[@]}"; do
    # frag_path so a bgzipped fragment file resolves; ChromBPNet reads .gz for
    # -ifrag. Without this the guard below would skip every cell type silently.
    local frag; frag=$(frag_path "$FRAG_DIR/${CT}_chr.bed")
    local pk="$FILT_DIR/${CT}_peaks_no_blacklist_chr.bed"
    local neg="$NEG_DIR/${CT}_nonpeaks_negatives.bed"
    local out="$MODEL_DIR/${CT}"
    [ -n "$frag" ] || { log "  $CT: missing $FRAG_DIR/${CT}_chr.bed[.gz] -- skipping"; continue; }
    for f in "$pk" "$neg"; do
      [ -s "$f" ] || { log "  $CT: missing $f -- skipping"; continue 2; }
    done
    if [ -s "$out/models/chrombpnet_nobias.h5" ]; then
      log "  $CT: model exists, skipping"; continue
    fi
    mkdir -p "$out"
    log "================ training $CT ================"
    run "$CHROMBPNET" pipeline \
      -ifrag "$frag" \
      -d "$DATA_TYPE" \
      -g "$GENOME" \
      -c "$CHROM_SIZES" \
      -p "$pk" \
      -n "$neg" \
      -fl "$FOLD_JSON" \
      -b "$bias_model" \
      -o "$out" 2>&1 | tee "$LOG_DIR/${CT}_train.log"
  done
}

# ============================================================================
# stage: interpret  (GPU) -- predictions, contributions, MoDISco
# ============================================================================
stage_interpret() {
  need_tool "$CHROMBPNET" chrombpnet CHROMBPNET

  for CT in "${CELLTYPES[@]}"; do
    local mdir="$MODEL_DIR/${CT}/models"
    # FIX 3: contributions must come from the bias-corrected model.
    local nobias="$mdir/chrombpnet_nobias.h5"
    local full="$mdir/chrombpnet.h5"
    local scaled="$mdir/bias_model_scaled.h5"
    local regions; regions=$(regions_for "$CT")

    [ -n "$regions" ] || { log "  $CT: no DAR regions -- skipping"; continue; }
    [ -s "$nobias" ]  || { log "  $CT: missing $nobias (train first) -- skipping"; continue; }

    local odir="$INTERP_DIR/${CT}"      # FIX 1: ${CT}, not the literal {i}
    mkdir -p "$odir"
    local n; n=$(nlines "$regions")
    local nx; nx=$(awk '$1=="chrX"' "$regions" | wc -l | tr -d ' ')
    log "================ $CT ($n regions, $nx on chrX) ================"

    if [ -s "$odir/${CT}_dars_bigwig.bw" ]; then
      log "  predictions exist, skipping pred_bw"
    else
      log "  pred_bw ..."
      run "$CHROMBPNET" pred_bw \
        -bm "$scaled" \
        -cm "$full" \
        -r "$regions" \
        -c "$CHROM_SIZES" -g "$GENOME" \
        -op "$odir/${CT}_dars_bigwig" 2>&1 | tee "$LOG_DIR/${CT}_predbw.log"
    fi

    if [ -s "$odir/${CT}_dars_contrib.counts_scores.h5" ]; then
      log "  contributions exist, skipping contribs_bw"
    else
      log "  contribs_bw (chrombpnet_nobias.h5) ..."
      run "$CHROMBPNET" contribs_bw \
        -m "$nobias" \
        -r "$regions" \
        -c "$CHROM_SIZES" -g "$GENOME" \
        -op "$odir/${CT}_dars_contrib" 2>&1 | tee "$LOG_DIR/${CT}_contribs.log"
    fi

    # MoDISco: refuse to cluster seqlets from too few regions.
    if [ "$n" -lt "$MODISCO_MIN_REGIONS" ]; then
      log "  SKIPPING MoDISco: $n regions < MODISCO_MIN_REGIONS=$MODISCO_MIN_REGIONS"
      log "  (expected for Macrophages -- its uninjured group is 22 cells)"
      continue
    fi
    need_tool "$MODISCO" modisco MODISCO
    if [ -s "$odir/${CT}_modisco.h5" ]; then
      log "  modisco output exists, skipping"
    else
      log "  modisco motifs ..."
      run "$MODISCO" motifs \
        -i "$odir/${CT}_dars_contrib.counts_scores.h5" \
        -n "$MODISCO_MAX_SEQLETS" \
        -o "$odir/${CT}_modisco.h5" -v 2>&1 | tee "$LOG_DIR/${CT}_modisco.log"
    fi
    log "  modisco report ..."
    mkdir -p "$odir/modisco_report"
    run "$MODISCO" report \
      -i "$odir/${CT}_modisco.h5" \
      -o "$odir/modisco_report/" \
      -s "$odir/modisco_report/" \
      -m "$MOTIFS"
  done
}

# ============================================================================
# preflight + dispatch
# ============================================================================
preflight() {
  local missing=0
  for f in "$GENOME" "$CHROM_SIZES" "$BLACKLIST"; do
    [ -e "$f" ] || { echo "missing reference: $f" >&2; missing=1; }
  done
  [ "$missing" = "0" ] || die "fix the missing references above"
  log "genome=$GENOME  chrom_sizes=$CHROM_SIZES  fold=$FOLD_JSON"
  log "cell types: ${CELLTYPES[*]}"
  [ "$DRY_RUN" = "1" ] && log "DRY RUN -- no commands will be executed"
  return 0
}

preflight
case "$STAGE" in
  bias-peaks) stage_bias_peaks ;;
  splits)     stage_splits ;;
  nonpeaks)   stage_nonpeaks ;;
  prep)       stage_bias_peaks; stage_splits; stage_nonpeaks ;;
  bias-train) stage_bias_train ;;
  train)      stage_train ;;
  interpret)  stage_interpret ;;
  all)        stage_bias_peaks; stage_splits; stage_nonpeaks
              stage_bias_train; stage_train; stage_interpret ;;
  *)          die "unknown stage: $STAGE
  stages: bias-peaks | splits | nonpeaks | prep | bias-train | train | interpret | all" ;;
esac
log "done: $STAGE"
