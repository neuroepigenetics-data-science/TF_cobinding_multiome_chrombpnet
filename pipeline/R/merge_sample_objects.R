# ============================================================================
# merge_sample_objects.R  (parameterised)
# ----------------------------------------------------------------------------
# Merge the cleaned per-sample objects from filtering.R into one dataset, drop
# RNA doublets, reprocess RNA + ATAC, build WNN neighbours / UMAPs / clustering,
# annotate cell types, and add per-modality (RNA-only / ATAC-only) UMAPs.
#
# Config-driven version of code/preprocessing/merge_sample_objects.R. This is the
# step that produces the critical handoff object everything downstream depends on
# (config `paths.merged_obj`, i.e. merged_multiome.rds).
#
# Behavioural notes where this departs from the published script:
#   * Per-cell experiment/condition/batch are read from samples.csv (keyed by the
#     cell's `sample`) instead of strsplit(sample, "_") — robust to sample names
#     that don't follow the EXPERIMENT_CONDITION_BATCH convention.
#   * The published script annotated off `wsnn_res.0.1` but only ran FindClusters
#     at resolution 0.2, so that column never existed. Here the resolutions in
#     `merge.cluster_resolutions` are all computed, and the annotation reads the
#     configured `annotation.cluster_col`.
#   * The published per-modality RunUMAP calls omitted `reduction=`, so both
#     defaulted to "pca" (the ATAC-only UMAP was built from pca, not lsi). Here
#     the RNA-only UMAP uses pca and the ATAC-only UMAP uses lsi.
#
# ** The cluster -> cell-type map in config is DATA-SPECIFIC. Re-derive it for
#    your own data before trusting `cluster_ids`. **
#
# Reads:   config paths.unmerged_obj_clean  (named list from filtering.R)
# Writes:  config paths.merged_obj          (single annotated Seurat object)
#
# Run from the repo root:
#   Rscript pipeline/R/merge_sample_objects.R [path/to/config.yaml]
# ============================================================================

# --- load helpers + config -------------------------------------------------
source("pipeline/R/pipeline_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "pipeline/config.yaml"
cfg <- load_config(config_path)

# shared preamble (libraries, seed, mm10 annotation, etc.)
source(cfg$project$source_r)

samples     <- load_samples(cfg)
m           <- cfg$merge
rna_assay   <- cfg$assays$rna
peaks_assay <- cfg$assays$peaks

# --- merge -----------------------------------------------------------------
object <- readRDS(cfg$paths$unmerged_obj_clean)
multiome <- merge(x = object[[1]], y = object[-1])   # any number of samples
rm(object); gc()

# --- per-cell metadata from the sample sheet -------------------------------
idx <- match(multiome$sample, samples$sample_id)
if (anyNA(idx)) {
  warning("cells with `sample` not in samples.csv: ",
          paste(unique(multiome$sample[is.na(idx)]), collapse = ", "))
}
for (col in m$sample_metadata_cols) {
  multiome[[col]] <- samples[[col]][idx]
}

# --- drop RNA doublets (flagged in filtering.R) ----------------------------
keep <- multiome@meta.data[[m$doublet_col]] %in% m$doublet_keep
multiome <- subset(multiome, cells = colnames(multiome)[keep])

# --- reprocess RNA + ATAC on the merged object -----------------------------
DefaultAssay(multiome) <- rna_assay
multiome <- multiome %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(multiome) <- peaks_assay
multiome <- multiome %>%
  FindTopFeatures(min.cutoff = m$top_features_cutoff) %>% RunTFIDF() %>% RunSVD()

# --- joint (WNN) neighbours, UMAP, clustering ------------------------------
dp <- m$wnn_dims_pca
dl <- m$wnn_dims_lsi
multiome <- multiome %>%
  FindMultiModalNeighbors(
    reduction.list = list("pca", "lsi"),
    dims.list = list(dp[1]:dp[2], dl[1]:dl[2]),
    modality.weight.name = "RNA.weight") %>%
  RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

for (res in m$cluster_resolutions) {
  multiome <- FindClusters(multiome, graph.name = "wsnn", algorithm = 3, resolution = res)
}

multiome <- RunUMAP(multiome, reduction = "pca", dims = dp[1]:dp[2], assay = rna_assay,
                    reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
multiome <- RunUMAP(multiome, reduction = "lsi", dims = dl[1]:dl[2], assay = peaks_assay,
                    reduction.name = "lsi.umap", reduction.key = "lsiUMAP_")

# --- annotate clusters -> cell types (DATA-SPECIFIC map) -------------------
a <- m$annotation
Idents(multiome) <- a$cluster_col
multiome <- do.call(RenameIdents, c(list(object = multiome), a$map))
Idents(multiome) <- factor(Idents(multiome), levels = unlist(a$levels))
multiome[[a$ident_col]] <- Idents(multiome)

# --- per-modality clustering (RNA-only / ATAC-only) — Fig 1 ----------------
mr <- m$modality_resolution
multiome <- multiome %>%
  FindNeighbors(reduction = "pca", dims = dp[1]:dp[2], graph.name = "RNA.weight") %>%
  RunUMAP(reduction = "pca", dims = dp[1]:dp[2],
          reduction.name = "rna.umap", reduction.key = "rnaUMAP_") %>%
  FindClusters(graph.name = "RNA.weight", algorithm = 3, resolution = mr)

multiome <- multiome %>%
  FindNeighbors(reduction = "lsi", dims = dl[1]:dl[2], graph.name = "ATAC.weight") %>%
  RunUMAP(reduction = "lsi", dims = dl[1]:dl[2],
          reduction.name = "lsi.umap", reduction.key = "lsiUMAP_") %>%
  FindClusters(graph.name = "ATAC.weight", algorithm = 3, resolution = mr)

# --- save ------------------------------------------------------------------
saveRDS(multiome, cfg$paths$merged_obj)
message("Wrote ", cfg$paths$merged_obj)
