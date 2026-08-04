#!/usr/bin/env bash
# ============================================================================
# entrypoint.sh — one-click Track A: annotated object + fragments -> ChromBPNet inputs
# ----------------------------------------------------------------------------
#   dge      step 2  injury DARs (peaks + RNA)      -> markers/
#   split    step 3  fragment BEDs by cell type     -> ML/celltype_fragments*/
#   dars     step 3b DAR BEDs                       -> ML/celltype_dars/
#   peaks    step 1  MACS2 peaks + DAR intersect    -> ML/{macs2_celltype_peaks,filtered_peaks,dars_intersect_peaks}/
#   all              every step above, in order
#   check            validate inputs and exit (no compute)
#   shell            drop into bash
#
# Steps are ordered by dependency, not by the numbering in run_chrombpnet.sh:
# `dars` needs `dge`, and `peaks` needs both `split` and `dars`.
#
# Every step is idempotent, so a re-run resumes rather than redoing work --
# except `split`, which refuses to write into a non-empty output dir because
# SplitFragments appends and would otherwise silently double every BED.
# ============================================================================
set -euo pipefail

CONFIG=${CONFIG:-pipeline/config.yaml}
RUN_R="mamba run --no-capture-output -n zamboni-r Rscript"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

check_inputs() {
  log "checking inputs (config: $CONFIG)"
  [ -f "$CONFIG" ] || die "config not found: $CONFIG"

  local obj
  obj=$($RUN_R -e "cat(yaml::read_yaml('$CONFIG')\$ml\$merged_obj)" 2>/dev/null | tail -1)
  [ -n "$obj" ] || die "could not read ml.merged_obj from $CONFIG"
  [ -f "$obj" ] || die "merged object not mounted: $obj
  Mount the directory holding it, e.g.:
    -v /host/path/TF_cobinding_multiome_chrombpnet:/work/host  (then point config at it)"
  log "  merged object: $obj ($(du -h "$obj" | cut -f1))"

  local n
  n=$(ls -1 data/*/outs/atac_fragments.tsv.gz 2>/dev/null | wc -l | tr -d ' ')
  [ "$n" -gt 0 ] || die "no fragment files under data/*/outs/ -- mount the data dir"
  log "  fragment files: $n"
  local missing_idx=0
  for f in data/*/outs/atac_fragments.tsv.gz; do
    [ -f "${f}.tbi" ] || { echo "  missing index: ${f}.tbi" >&2; missing_idx=1; }
  done
  [ "$missing_idx" = "0" ] || die "fragment .tbi indexes missing -- run: tabix -p bed <file>"

  for f in meta/mm10.chr.sizes meta/mm10-blacklist.v2.bed; do
    [ -f "$f" ] || die "missing reference: $f"
  done
  log "  references OK"
  log "  tools: macs2 $(macs2 --version 2>&1 | awk '{print $2}'), bedtools $(bedtools --version | awk '{print $2}')"
  log "inputs OK"
}

step_dge()   { log "== step 2: injury DARs =="
               $RUN_R pipeline/R/differential_analysis_injury.R "$CONFIG" --only peaks
               $RUN_R pipeline/R/differential_analysis_injury.R "$CONFIG" --only rna; }
step_split() { log "== step 3: fragment BEDs =="
               $RUN_R pipeline/R/prepare_ml_input.R "$CONFIG"; }
step_dars()  { log "== step 3b: DAR BEDs =="
               $RUN_R pipeline/R/prepare_ml_dars.R "$CONFIG"; }
step_peaks() { log "== step 1: MACS2 peaks + DAR intersect =="
               bash pipeline/call_celltype_peaks.sh; }

case "${1:-all}" in
  check) check_inputs ;;
  dge)   check_inputs; step_dge ;;
  split) check_inputs; step_split ;;
  dars)  step_dars ;;
  peaks) step_peaks ;;
  all)   check_inputs; step_dge; step_split; step_dars; step_peaks
         log "== Track A complete ==" ;;
  shell) exec /bin/bash ;;
  *)     die "unknown command: $1 (expected: all|check|dge|split|dars|peaks|shell)" ;;
esac
