# ============================================================================
# differential_analysis_injury.R  (parameterised)
# ----------------------------------------------------------------------------
# Track A, step 2. Per-cell-type differential analysis of each injury timepoint
# against uninjured, on the ATAC peak assay and on RNA. The peaks output is the
# file the ChromBPNet DAR / MoDISco-interpretation branch consumes.
#
# Config-driven version of the "differential analysis injury" section of
# code/explore/differential_analysis.R (its lines ~85-143).
#
# Reads:   config dge.merged_obj (the annotated merged object)
# Writes:  config dge.peaks_out  markers_Uvsothers_peaks_all.rds
#          config dge.rna_out    markers_Uvsothers_rna_all.rds
#          -- each a named list (one element per cell type) of rbind-ed
#             FindMarkers tables, filtered as in the original.
#
# Deliberately does NOT source code/source.R (unconditional library(ArchR)).
#
# Deliberate departures from code/explore/differential_analysis.R:
#   1. Idents are set explicitly from config dge.celltype_col. The original
#      relies on whatever the object's active ident happens to be
#      (`clusters <- Idents(multiome) %>% unique()`), which silently changes the
#      entire analysis if the object was saved with a different active ident.
#   2. Group sizes are logged, and a loud warning is emitted for any comparison
#      resting on fewer than dge.warn_below_cells cells. On this dataset
#      Macrophages-vs-U has 22 cells and Endothelial 3dpi has 9, so some DAR
#      sets are statistically thin. Nothing is skipped -- the published
#      thresholds are preserved -- but it is recorded rather than invisible.
#   3. `subset(multiome, ident = i)` -> `idents = i`. The original relies on R
#      partial-argument matching to reach Seurat's `idents`; spelled out here.
#   4. Which halves run is selectable (--only peaks|rna|both). The peaks half is
#      what ChromBPNet needs; RNA is independent of it.
#   5. The object is slimmed with DietSeurat() to just the assay under test
#      before looping. Without it FindMarkers exhausts memory on the largest
#      cell type on a 24 GB machine. A memory measure only -- the retained
#      assay's data is untouched, so results are unchanged.
#
# Run from the repo root:
#   Rscript pipeline/R/differential_analysis_injury.R [config.yaml] [--only peaks]
# ============================================================================

source("pipeline/R/pipeline_utils.R")

args <- commandArgs(trailingOnly = TRUE)
only <- "both"
if ("--only" %in% args) {
  only <- args[[which(args == "--only") + 1]]
  if (!only %in% c("peaks", "rna", "both")) {
    stop("--only must be one of: peaks, rna, both")
  }
}
cfg_args    <- args[!grepl("^--", args) & !(args %in% c("peaks", "rna", "both"))]
config_path <- if (length(cfg_args) >= 1) cfg_args[[1]] else "pipeline/config.yaml"
cfg <- load_config(config_path)
d   <- cfg$dge

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
})

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

# --- load ------------------------------------------------------------------
say("loading ", d$merged_obj)
multiome <- readRDS(d$merged_obj)
say("loaded: ", ncol(multiome), " cells")

for (col in c(d$celltype_col, d$condition_col)) {
  if (!col %in% colnames(multiome@meta.data)) {
    stop("metadata column not found: ", col)
  }
}

# --- slim to fit in RAM (departure 5) --------------------------------------
# The object carries 4 assays (peaks, RNA, RNA_raw, chromvar), 5 reductions and
# the WNN graphs. On a 24 GB machine that leaves no room to subset the largest
# cell type: FindMarkers died with "vector memory exhausted" on Oligodendrocytes
# (20,816 cells). FindMarkers only ever touches the assay under test, so drop
# everything else up front. Purely a memory measure -- the retained assay's
# counts/data are untouched, so results are identical.
assays_needed <- switch(only,
  peaks = d$peaks_assay,
  rna   = d$rna_assay,
  both  = c(d$peaks_assay, d$rna_assay))

say("assays present: ", paste(Assays(multiome), collapse = ", "))
g <- Graphs(multiome)
say("graphs present: ", if (length(g)) paste(g, collapse = ", ") else "none")
say("size before slim: ", format(object.size(multiome), units = "auto"))

# DietSeurat() refuses to drop the object's DefaultAssay ("peaks" as saved), so
# point it at an assay we are keeping first -- otherwise --only rna fails with
# "The default assay is slated to be removed".
DefaultAssay(multiome) <- assays_needed[[1]]
multiome <- DietSeurat(multiome, assays = assays_needed,
                       dimreducs = NULL, graphs = NULL)
invisible(gc())
say("size after slim:  ", format(object.size(multiome), units = "auto"),
    "  (kept: ", paste(assays_needed, collapse = ", "), ")")
if (only == "both") {
  say("NOTE: --only both keeps two assays; if this OOMs, run ",
      "--only peaks and --only rna as separate processes")
}

# Departure 1: set idents explicitly rather than trusting the saved active ident.
Idents(multiome) <- d$celltype_col
clusters <- levels(Idents(multiome))
if (is.null(clusters)) clusters <- unique(as.character(Idents(multiome)))
say("cell types: ", paste(clusters, collapse = ", "))

multiome@meta.data[[d$condition_col]] <-
  factor(multiome@meta.data[[d$condition_col]], levels = d$condition_levels)

# --- group sizes (departure 2) ---------------------------------------------
tab <- table(multiome@meta.data[[d$celltype_col]],
             multiome@meta.data[[d$condition_col]])
say("celltype x condition group sizes:")
print(tab)

thin <- which(tab < d$warn_below_cells, arr.ind = TRUE)
if (nrow(thin)) {
  msg <- apply(thin, 1, function(r)
    paste0(rownames(tab)[r["row"]], " x ", colnames(tab)[r["col"]],
           " = ", tab[r["row"], r["col"]]))
  warning("comparisons resting on < ", d$warn_below_cells,
          " cells (results are statistically thin):\n  ",
          paste(msg, collapse = "\n  "))
  say("THIN GROUPS: ", paste(msg, collapse = "; "))
}

dir.create(d$outdir, recursive = TRUE, showWarnings = FALSE)

# --- the per-cell-type / per-timepoint loop --------------------------------
run_markers <- function(assay, extra_args = list()) {
  out <- list()
  for (i in clusters) {
    celltype <- subset(multiome, idents = i)
    Idents(celltype) <- d$condition_col
    DefaultAssay(celltype) <- assay

    per_tp <- list()
    for (j in d$timepoints) {
      n1 <- sum(celltype@meta.data[[d$condition_col]] == j, na.rm = TRUE)
      n2 <- sum(celltype@meta.data[[d$condition_col]] == d$baseline, na.rm = TRUE)
      if (n1 < 3 || n2 < 3) {
        say("  SKIP ", i, " ", j, " vs ", d$baseline,
            " (n=", n1, " vs ", n2, "; FindMarkers needs >= 3 per group)")
        next
      }
      fm_args <- c(list(object = celltype, ident.1 = j, ident.2 = d$baseline,
                        only.pos = FALSE, min.pct = d$min_pct),
                   extra_args)
      res <- do.call(FindMarkers, fm_args)
      res$comparison <- paste0(i, "_", j)
      res$peak_or_gene <- rownames(res)
      # keep the original's column names so downstream code is unchanged
      if (assay == d$peaks_assay) res$peak <- rownames(res) else res$gene <- rownames(res)
      per_tp[[j]] <- res
      say("  ", i, " ", j, " vs ", d$baseline, ": n=", n1, " vs ", n2,
          " -> ", nrow(res), " rows")
    }
    out[[i]] <- if (length(per_tp)) do.call("rbind", per_tp) else NULL
    rm(celltype); gc(verbose = FALSE)
  }
  out
}

apply_filter <- function(lst) {
  for (i in names(lst)) {
    if (is.null(lst[[i]])) next
    lst[[i]] <- lst[[i]] %>%
      dplyr::filter(p_val < d$filter_p & abs(avg_log2FC) > d$filter_abs_log2fc)
  }
  lst
}

if (only %in% c("peaks", "both")) {
  say("=== PEAKS (", d$peaks_assay, ") ===")
  peaks_all <- apply_filter(run_markers(d$peaks_assay))
  saveRDS(peaks_all, d$peaks_out)
  say("wrote ", d$peaks_out)
  for (i in names(peaks_all)) {
    say("  ", i, ": ", if (is.null(peaks_all[[i]])) 0 else nrow(peaks_all[[i]]),
        " DARs after filter")
  }
}

if (only %in% c("rna", "both")) {
  say("=== RNA (", d$rna_assay, ") ===")
  rna_all <- apply_filter(
    run_markers(d$rna_assay,
                extra_args = list(logfc.threshold = d$rna_logfc_threshold))
  )
  saveRDS(rna_all, d$rna_out)
  say("wrote ", d$rna_out)
}

say("done")
