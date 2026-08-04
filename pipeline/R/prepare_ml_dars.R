# ============================================================================
# prepare_ml_dars.R  (parameterised)
# ----------------------------------------------------------------------------
# Track A. Turns the per-cell-type injury DARs into BED files, the regions
# ChromBPNet computes contribution scores over (`contribs_bw`) before MoDISco
# clusters them into motifs.
#
# Config-driven version of the tail of code/machine_learning/prepare_ml_input.R
# (its lines ~48-63). Split into its own script deliberately: the fragment
# splits in prepare_ml_input.R take ~30 min over 21 GB, and there is no reason
# to redo them every time the DAR thresholds change.
#
# Reads:   config ml.dars_markers  (from differential_analysis_injury.R)
# Writes:  config ml.dars_outdir   <celltype>_I_dars.bed
#
# Deliberate departures from the published code:
#   1. The peak coordinate comes from the `peak` column rather than being
#      re-derived with strsplit(rownames(x), "[.]")[[2]]. After rbind() the
#      rownames are "<timepoint>.<peak>", but R also de-duplicates repeated
#      rownames by appending suffixes, which corrupts that parse. The `peak`
#      column is written before the rbind and is unambiguous.
#   2. Peaks recurring at several timepoints are collapsed (ml.dars_dedupe).
#      The original emits one BED line per DAR *per timepoint*, so a region
#      opening at both 1dpi and 7dpi appears twice; contribution scores would
#      then be computed twice over it and MoDISco would double-count its
#      seqlets. The count dropped is logged. Set dars_dedupe: false for the
#      original behaviour.
#   3. BEDs are written to ml.dars_outdir, not the working directory.
#
# NOTE: run_chrombpnet.sh reads dars_intersect_peaks/<celltype>_intersect_dars.bed
# -- i.e. these DARs intersected with that cell type's MACS2 peak set. That
# intersect step exists in neither published script and needs the per-cell-type
# peaks (run_chrombpnet.sh's earlier macs2 loop), so it is NOT done here.
# Contribution scores should only be computed over regions the model saw.
#
# Run from the repo root:
#   Rscript pipeline/R/prepare_ml_dars.R [config.yaml]
# ============================================================================

source("pipeline/R/pipeline_utils.R")

args        <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "pipeline/config.yaml"
cfg <- load_config(config_path)
m   <- cfg$ml

suppressPackageStartupMessages(library(dplyr))
say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

if (!file.exists(m$dars_markers)) {
  stop("markers file not found: ", m$dars_markers,
       "\nRun pipeline/R/differential_analysis_injury.R first.")
}

dars <- readRDS(m$dars_markers)
say("loaded ", m$dars_markers, ": ", length(dars), " cell types")

dir.create(m$dars_outdir, recursive = TRUE, showWarnings = FALSE)

total <- 0L
for (i in names(dars)) {
  df <- dars[[i]]
  if (is.null(df) || nrow(df) == 0) {
    say("  ", i, ": no DARs, skipped")
    next
  }
  if (!"peak" %in% colnames(df)) {
    stop("cell type ", i, " has no `peak` column; was it produced by ",
         "differential_analysis_injury.R?")
  }

  # positive only: regions opening after injury
  keep <- df %>%
    dplyr::filter(p_val < m$dars_max_p & avg_log2FC > m$dars_min_log2fc)
  n_before <- nrow(keep)

  peaks <- keep$peak
  if (isTRUE(m$dars_dedupe)) peaks <- unique(peaks)
  n_dropped <- n_before - length(peaks)

  # "chr1-3000-3500" -> c(chr1, 3000, 3500)
  parts <- strsplit(peaks, "-", fixed = TRUE)
  bad   <- lengths(parts) != 3
  if (any(bad)) {
    stop("unparseable peak name(s) for ", i, ": ",
         paste(head(peaks[bad], 3), collapse = ", "))
  }
  bed <- data.frame(
    chr   = vapply(parts, `[`, character(1), 1),
    start = vapply(parts, `[`, character(1), 2),
    end   = vapply(parts, `[`, character(1), 3),
    stringsAsFactors = FALSE
  )
  # genome order, so downstream bedtools/macs2 behave predictably
  bed <- bed[order(bed$chr, as.numeric(bed$start)), ]

  out <- file.path(m$dars_outdir, paste0(i, "_I_dars.bed"))
  write.table(bed, out, row.names = FALSE, col.names = FALSE,
              sep = "\t", quote = FALSE)
  total <- total + nrow(bed)
  say("  ", i, ": ", nrow(bed), " regions",
      if (n_dropped > 0) paste0(" (", n_dropped, " duplicate rows collapsed)") else "",
      " -> ", out)
}

say("done: ", total, " DAR regions across ", length(dars), " cell types")
