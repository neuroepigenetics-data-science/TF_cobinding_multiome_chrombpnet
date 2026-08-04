# ============================================================================
# filtering.R  (parameterised)
# ----------------------------------------------------------------------------
# Per-sample QC on the unmerged Seurat objects from create_object.R:
#   1. drop nuclei with incomplete metadata
#   2. flag per-cluster outliers (IQR rule) on RNA/ATAC QC metrics, then remove
#   3. manually drop specific low-quality clusters (per sample)
#   4. recompute WNN neighbours / UMAPs / clustering on the cleaned cells
#   5. RNA doublets (DoubletFinder) + ATAC doublets (scDblFinder)
#   6. keep only relevant clusters for a few samples (drop cells from another
#      experiment)
#
# Config-driven version of code/preprocessing/filtering.R. Behavioural notes
# where this departs from the published script (both intentional):
#   * The original referenced peak assays "macs2_peaks" / "cellranger_peaks"
#     that create_object.R never builds. Both are the single ATAC assay named
#     in config `assays.peaks`.
#   * The original doublet loop did `best_pK <- best_pK[i]`, which overwrites the
#     vector with a scalar on the first iteration, so every later sample got
#     pK = NA. Per-sample pK is read here from samples.csv (`best_pK`) instead.
#   * DoubletFinder is called via the 2.0.6 API (doubletFinder(), reuse.pANN=NULL)
#     rather than the paper's doubletFinder_v3(); 2.0.6 is the Seurat-v5-compatible
#     release. Same algorithm/params, only the entry-point name/default changed.
#
# Reads:   config paths.unmerged_obj          (named list, one Seurat obj/sample)
# Writes:  config paths.unmerged_obj_cb_clean  (checkpoint, after step 4)
#          config paths.unmerged_obj_clean     (final, after step 6)
#
# Run from the repo root:
#   Rscript pipeline/R/filtering.R [path/to/config.yaml]
# ============================================================================

# --- load helpers + config -------------------------------------------------
source("pipeline/R/pipeline_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "pipeline/config.yaml"
cfg <- load_config(config_path)

# shared preamble (libraries incl. DoubletFinder/scDblFinder, seed, etc.)
source(cfg$project$source_r)

samples     <- load_samples(cfg)
p           <- cfg$filtering
peaks_assay <- cfg$assays$peaks
rna_assay   <- cfg$assays$rna
cluster_var <- p$cluster_var

object <- readRDS(cfg$paths$unmerged_obj)
stopifnot(length(object) == nrow(samples))

# --- 1-2. per-cluster outlier detection & removal --------------------------
# A nucleus is an outlier if any metric falls outside `outlier_iqr_mult` IQRs of
# its cluster, in the tail direction configured for that metric.
is_outlier <- function(feature, tail = c("both", "upper", "lower"), iqr_mult = 3) {
  tail  <- match.arg(tail)
  upper <- feature > quantile(feature, .75, na.rm = TRUE) + iqr_mult * IQR(feature, na.rm = TRUE)
  lower <- feature < quantile(feature, .25, na.rm = TRUE) - iqr_mult * IQR(feature, na.rm = TRUE)
  switch(tail, upper = upper, lower = lower, both = upper | lower)
}

metrics      <- p$outlier_metrics          # named list: metadata column -> tail
iqr_mult     <- p$outlier_iqr_mult
outlier_cols <- paste0("is_outlier_", names(metrics))

for (i in seq_along(object)) {
  # remove cells with missing values in metadata
  na_cells <- rownames(object[[i]]@meta.data)[!complete.cases(object[[i]]@meta.data)]
  if (length(na_cells) > 0) {
    object[[i]] <- subset(object[[i]], cells = na_cells, invert = TRUE)
  }

  # groupwise outlier flags, one column per metric
  md <- object[[i]]@meta.data %>% group_by(.data[[cluster_var]])
  for (m in names(metrics)) {
    md <- md %>%
      mutate(!!paste0("is_outlier_", m) := is_outlier(.data[[m]], metrics[[m]], iqr_mult))
  }
  md <- md %>% ungroup() %>% as.data.frame()
  rownames(md) <- colnames(object[[i]])

  md$is_outlier <- rowSums(md[outlier_cols], na.rm = TRUE) > 0
  object[[i]]@meta.data <- md
}

# removes ~10% of cells
object <- lapply(object, function(x) subset(x, subset = is_outlier == FALSE))

# --- 3. manual removal of specific low-quality clusters --------------------
for (sid in names(p$manual_cluster_drop)) {
  if (is.null(object[[sid]])) {
    warning("manual_cluster_drop: sample '", sid, "' not in object; skipping")
    next
  }
  drop <- as.character(p$manual_cluster_drop[[sid]])
  Idents(object[[sid]]) <- cluster_var
  object[[sid]] <- subset(object[[sid]], idents = drop, invert = TRUE)
}

# --- 4. recompute WNN neighbours / UMAPs / clustering ----------------------
# Reuses the pca / lsi reductions carried over from create_object.R.
dp <- p$wnn_dims_pca
dl <- p$wnn_dims_lsi
object <- lapply(object, function(x) {
  x <- x %>%
    FindMultiModalNeighbors(
      reduction.list = list("pca", "lsi"),
      dims.list = list(dp[1]:dp[2], dl[1]:dl[2]),
      modality.weight.name = "RNA.weight") %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, resolution = p$wnn_resolutions[[1]])

  x <- RunUMAP(x, reduction = "pca", dims = dp[1]:dp[2], assay = rna_assay,
               reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  x <- RunUMAP(x, reduction = "lsi", dims = dl[1]:dl[2], assay = peaks_assay,
               reduction.name = "lsi.umap", reduction.key = "lsiUMAP_")

  x <- FindClusters(x, graph.name = "wsnn", algorithm = 3, resolution = p$wnn_resolutions[[2]])
  x
})

saveRDS(object, cfg$paths$unmerged_obj_cb_clean)
message("Wrote checkpoint ", cfg$paths$unmerged_obj_cb_clean)

# --- 5a. RNA doublets (DoubletFinder) --------------------------------------
sig_pc <- p$sig_pc
for (i in seq_along(object)) {
  x   <- object[[i]]
  sid <- names(object)[i]
  best_pK <- as.numeric(samples$best_pK[match(sid, samples$sample_id)])  # per-sample, from the sheet

  homotypic.prop <- modelHomotypic(x@meta.data[[cluster_var]])
  doublet_rate   <- (p$doublet_rate_slope * nrow(x@meta.data) - p$doublet_rate_intercept) / 100
  nExp_poi       <- round(doublet_rate * nrow(x@meta.data))
  nExp_poi.adj   <- round(nExp_poi * (1 - homotypic.prop))

  # DoubletFinder 2.0.6 API (Seurat-v5 compatible): doubletFinder() not
  # doubletFinder_v3(), and reuse.pANN defaults to NULL (passing FALSE breaks).
  x <- doubletFinder(x, PCs = 1:sig_pc, pN = p$doublet_pN, pK = best_pK,
                     nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
  x <- doubletFinder(x, PCs = 1:sig_pc, pN = p$doublet_pN, pK = best_pK,
                     nExp = nExp_poi.adj,
                     reuse.pANN = colnames(x@meta.data)[grep("pANN_", colnames(x@meta.data))],
                     sct = FALSE)

  df_class <- grep("DF.classifications_", colnames(x@meta.data))
  x$DF_doublet <- ifelse(x[[colnames(x@meta.data)[df_class[1]]]] == "Doublet" &
                         x[[colnames(x@meta.data)[df_class[2]]]] == "Doublet",
                         "Doublet", "Singlet")
  object[[i]] <- x
}

# --- 5b. ATAC doublets (scDblFinder) ---------------------------------------
for (i in seq_along(object)) {
  x <- object[[i]]
  doublet_rate       <- (p$doublet_rate_slope * nrow(x@meta.data) - p$doublet_rate_intercept) / 100
  artificialDoublets <- round(doublet_rate * nrow(x@meta.data))

  scdbl_obj <- scDblFinder(
    GetAssayData(x, assay = peaks_assay, slot = "counts"),
    clusters           = x@meta.data[[cluster_var]],
    nfeatures          = p$scdbl_nfeatures,   # aggregated features
    dims               = 1:sig_pc,
    artificialDoublets = artificialDoublets,
    aggregateFeatures  = TRUE,
    verbose            = TRUE,
    processing         = "normFeatures")

  scdbl_out <- as.data.frame(scdbl_obj@colData@listData)
  rownames(scdbl_out) <- scdbl_obj@colData@rownames
  x <- AddMetaData(x, scdbl_out)
  object[[i]] <- x
}

# --- 6. keep only relevant clusters (drop cells from another experiment) ----
for (sid in names(p$manual_cluster_keep)) {
  if (is.null(object[[sid]])) {
    warning("manual_cluster_keep: sample '", sid, "' not in object; skipping")
    next
  }
  keep <- as.character(p$manual_cluster_keep[[sid]])
  Idents(object[[sid]]) <- cluster_var
  object[[sid]] <- subset(object[[sid]], idents = keep)
}

saveRDS(object, cfg$paths$unmerged_obj_clean)
message("Wrote ", cfg$paths$unmerged_obj_clean)
