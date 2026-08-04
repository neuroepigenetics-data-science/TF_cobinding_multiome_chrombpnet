# ============================================================================
# create_object.R  (parameterised)
# ----------------------------------------------------------------------------
# Build per-sample Seurat objects (RNA + ATAC peaks) from Cell Ranger ARC
# output and the ArchR feature matrix. Config-driven version of
# code/preprocessing/create_object.R.
#
# Reads:   config paths.feature_matrix  (peak x cell counts from archr step)
#          per sample: filtered_feature_bc_matrix, cellbender h5, atac fragments
# Writes:  config paths.unmerged_obj     (named list of per-sample Seurat objects)
#
# DELIBERATE DEPARTURE — degenerate-dimension guard (added 2026-07-28)
# ---------------------------------------------------------------------------
# The published script picks WNN dimensionality adaptively:
#     sig_pc_pca <- length(which(pca@stdev > 2));  dims 1:sig_pc_pca
#     sig_pc_lsi <- length(which(lsi@stdev > 2));  dims 2:sig_pc_lsi
# Applying one cutoff (2) to both is fragile, because PCA and LSI stdev are on
# very different scales. Signac's RunSVD normalises the singular values by
# sqrt(n_features - 1) -- verified identical in EVERY Signac release from 0.2.4
# to 1.17.0, so this is NOT version drift -- which puts LSI stdev around
# 6.6 / 2.8 / 2.0 / 1.4 ... Only the first two or three clear a cutoff of 2.
#
# Because LSI discards component 1 (sequencing depth) and indexes 2:sig_pc_lsi,
# sig_pc_lsi == 2 yields the single dimension 2:2. The embedding then collapses
# from a matrix to a vector and FindMultiModalNeighbors fails inside L2Norm with
# "dim(X) must have a positive length". Observed on both smoke samples, where
# LSI component 3 landed at 1.97 and 1.63 -- i.e. the published code runs or
# crashes depending on which side of 2.0 component 3 happens to fall.
#
# The guard below only raises sig_pc_* when it is too small to index at all, so
# it is a NO-OP on any dataset where the published code produced a result. It
# cannot silently shift results: either it never fires (output is faithful) or
# it fires in a case the original could not have produced output for.
# The per-sample sig_pc values are logged so you can see whether it fired.
#
# ** TO IMPROVE: this heuristic should be replaced with a scale-appropriate
#    criterion for LSI (e.g. depth-correlation cutoff, or the fixed dims the
#    paper's own filtering.R / merge_sample_objects.R use: 2:15 / 2:30) --
#    deferred so the reproduction is not confounded by a methodology change. **
#
# Run from the repo root:
#   Rscript pipeline/R/create_object.R [path/to/config.yaml]
# ============================================================================

# --- load helpers + config -------------------------------------------------
source("pipeline/R/pipeline_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "pipeline/config.yaml"
cfg <- load_config(config_path)

# shared preamble (libraries, mm10 `annotation`, JASPAR, palette, seed)
source(cfg$project$source_r)

samples <- load_samples(cfg)
p <- cfg$create_object

# --- inputs ----------------------------------------------------------------
# feature matrices (list of peak x cell count matrices, one per sample, in the
# same order as the ArchR step produced them = sample-sheet order)
mtx_cluster <- readRDS(cfg$paths$feature_matrix)
stopifnot(length(mtx_cluster) == nrow(samples))

# --- build one Seurat object per sample ------------------------------------
unmerged_obj <- vector("list", nrow(samples))
names(unmerged_obj) <- samples$sample_id

for (i in seq_len(nrow(samples))) {
  sid      <- samples$sample_id[i]
  spath    <- samples$path[i]
  frag_tsv <- outs_path(cfg, spath, "atac_fragments")

  # ---- ATAC peaks assay ----
  assay <- CreateChromatinAssay(
    counts     = mtx_cluster[[i]],
    fragments  = frag_tsv,
    annotation = annotation           # from source.R
  )
  x <- CreateSeuratObject(assay, assay = "peaks")
  x$sample <- sid

  message(sid, " - all cells: ", ncol(x))
  x <- subset(x, subset = nCount_peaks > p$min_ncount_peaks)
  message(sid, " - after nCount_peaks > ", p$min_ncount_peaks, ": ", ncol(x))

  # ---- raw RNA assay ----
  counts_gex <- Read10X(outs_path(cfg, spath, "filtered_matrix"))$`Gene Expression`
  rownames(counts_gex) <- fix_gene_names(rownames(counts_gex), cfg)
  counts_gex <- counts_gex[, colnames(x)]
  x[["RNA_raw"]] <- CreateAssayObject(counts_gex)
  x <- subset(x, subset = nCount_RNA_raw > p$min_ncount_rna_raw)

  # ---- cellbender (ambient-corrected) RNA assay ----
  counts_cb <- Read10X_h5(outs_path(cfg, spath, "cellbender_h5"))
  rownames(counts_cb) <- fix_gene_names(rownames(counts_cb), cfg)
  counts_cb <- counts_cb[, colnames(x)]
  x[["RNA"]] <- CreateAssayObject(counts_cb)

  unmerged_obj[[i]] <- x
}

# --- QC metrics ------------------------------------------------------------
unmerged_obj <- lapply(unmerged_obj, function(x) {
  DefaultAssay(x) <- "peaks"
  x <- NucleosomeSignal(x)
  x <- TSSEnrichment(x)
  DefaultAssay(x) <- "RNA"
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = cfg$genome$mt_pattern)
  x
})

# --- process RNA + ATAC ----------------------------------------------------
unmerged_obj <- lapply(unmerged_obj, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- x %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

  DefaultAssay(x) <- "peaks"
  x <- x %>% FindTopFeatures(min.cutoff = p$top_features_cutoff) %>% RunTFIDF() %>% RunSVD()
  x
})

# --- joint (WNN) neighbours, UMAPs, clustering -----------------------------
unmerged_obj <- lapply(unmerged_obj, function(x) {
  sid <- as.character(x$sample[1])

  sig_pc_pca <- length(which(x@reductions[["pca"]]@stdev > p$pca_stdev_cutoff))
  sig_pc_lsi <- length(which(x@reductions[["lsi"]]@stdev > p$pca_stdev_cutoff))

  # Degenerate-dimension guard (see header). No-op unless the adaptive counts
  # are too small to index a >=2-dimension embedding.
  guarded_pca <- max(sig_pc_pca, p$min_dims_pca)
  guarded_lsi <- max(sig_pc_lsi, p$min_dims_lsi)
  if (guarded_pca != sig_pc_pca || guarded_lsi != sig_pc_lsi) {
    warning(sid, ": degenerate-dimension GUARD FIRED — ",
            "sig_pc_pca ", sig_pc_pca, "->", guarded_pca, ", ",
            "sig_pc_lsi ", sig_pc_lsi, "->", guarded_lsi,
            ". Published dims would have collapsed; see create_object.R header.",
            call. = FALSE, immediate. = TRUE)
  }
  message(sid, " - dims: pca 1:", guarded_pca, "  lsi 2:", guarded_lsi,
          "  (adaptive counts: pca ", sig_pc_pca, ", lsi ", sig_pc_lsi, ")")
  sig_pc_pca <- guarded_pca
  sig_pc_lsi <- guarded_lsi

  x <- x %>%
    FindMultiModalNeighbors(
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:sig_pc_pca, 2:sig_pc_lsi),
      modality.weight.name = "RNA.weight") %>%
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") %>%
    FindClusters(graph.name = "wsnn", algorithm = 3, resolution = p$wnn_resolutions[[1]])

  x <- RunUMAP(x, reduction = "pca", dims = 1:sig_pc_pca, assay = "RNA",
               reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
  x <- RunUMAP(x, reduction = "lsi", dims = 2:sig_pc_lsi, assay = "peaks",
               reduction.name = "lsi.umap", reduction.key = "lsiUMAP_")

  x <- FindClusters(x, graph.name = "wsnn", algorithm = 3, resolution = p$wnn_resolutions[[2]])
  x
})

# --- save ------------------------------------------------------------------
saveRDS(unmerged_obj, cfg$paths$unmerged_obj)
message("Wrote ", cfg$paths$unmerged_obj)
