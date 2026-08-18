# ============================================================================
# archr_peak_calling.R  (parameterised)
# ----------------------------------------------------------------------------
# ArchR-based peak calling — the FIRST step of the pipeline. From the per-sample
# ATAC fragments it builds arrow files, iterative-LSI dimensionality reduction,
# Harmony batch correction, clustering, a UMAP, then calls a reproducible MACS2
# peak set per cluster and quantifies per-sample counts over that peak set.
#
# Config-driven version of code/preprocessing/archr_peak_calling.R.
#
# Reads:   per sample: atac_fragments.tsv.gz + a barcode-metrics csv (valid cells)
# Writes:  config paths.archr_peaks     (GRanges peak set)
#          config paths.feature_matrix  (list of peak x cell matrices, one per
#                                        sample, in sample-sheet order — this is
#                                        the input to create_object.R)
#          plus the saved ArchRProject at archr.proj_out_dir
#
# Deliberate departure from code/preprocessing/archr_peak_calling.R:
#   createArrowFiles() is called with force = config archr.force_arrows
#   (default FALSE) instead of a hard-coded TRUE, so an interrupted run resumes
#   rather than rebuilding all 21 arrows from scratch — the arrow step is the
#   multi-hour part and was previously all-or-nothing. This cannot change
#   results: ArchR stamps Metadata/Completed only on an arrow it moved out of
#   tmp/ intact, reuses exactly those, and deletes+rebuilds any file missing the
#   stamp, so a partial arrow is never adopted. Arrow content is fully
#   determined by the fragments, valid barcodes, and min_tss/min_frags — all
#   unchanged across a restart. Set archr.force_arrows: true to force rebuilds.
#
# Requires ArchR + MACS2. ArchR is NOT on the ChromBPNet critical path -- the
# peak set it builds is already baked into the author's annotated object, and
# run_chrombpnet.sh re-calls peaks with MACS2 regardless. Install it only for
# Track B (reproduction, or processing our own data).
#
# Run from the repo root:
#   Rscript pipeline/R/archr_peak_calling.R [path/to/config.yaml]
# ============================================================================

# --- load helpers + config -------------------------------------------------
source("pipeline/R/pipeline_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "pipeline/config.yaml"
cfg <- load_config(config_path)

# shared preamble (libraries incl. ArchR/Signac, seed, mm10 annotation, etc.)
source(cfg$project$source_r)

samples <- load_samples(cfg)
a       <- cfg$archr

addArchRThreads(threads = a$threads)
addArchRGenome(a$genome)

# --- inputs ----------------------------------------------------------------
inputFiles <- vapply(samples$path, function(sp) outs_path(cfg, sp, "atac_fragments"),
                     character(1))
names(inputFiles) <- samples$sample_id

# valid (called-cell) barcodes per sample — layout handled by the config
barcodeList <- lapply(samples$path, function(sp) valid_barcodes(cfg, sp))
names(barcodeList) <- samples$sample_id

# --- arrow files -----------------------------------------------------------
# Resume support (see header): reuse arrows a previous interrupted run finished.
# Absent from an older config -> FALSE, the resumable behaviour.
force_arrows <- isTRUE(a$force_arrows)
message("createArrowFiles: force = ", force_arrows,
        if (!force_arrows) "  (completed arrow files will be reused)" else "")

createArrowFiles(
  inputFiles       = inputFiles,
  sampleNames      = names(inputFiles),
  validBarcodes    = barcodeList,
  geneAnnotation   = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  QCDir            = a$qc_dir,
  logFile          = createLogFile("createArrows"),
  minTSS           = a$min_tss,     # don't set too high; can filter up later
  minFrags         = a$min_frags,
  addTileMat       = a$add_tile_mat,
  addGeneScoreMat  = a$add_gene_score_mat,
  force            = force_arrows
)

# arrow files are written as "<sampleName>.arrow" in the working directory
ArrowFiles <- paste0(samples$sample_id, ".arrow")
archrproj <- ArchRProject(
  ArrowFiles      = ArrowFiles,
  outputDirectory = a$proj_tmp_dir,
  copyArrows      = a$copy_arrows   # keep an original copy of the arrows
)

# --- dimensionality reduction / clustering / UMAP --------------------------
lsi_dims <- a$lsi_dims
archrproj <- addIterativeLSI(
  ArchRProj     = archrproj,
  useMatrix     = "TileMatrix",
  name          = "IterativeLSI",
  iterations    = a$lsi_iterations,
  clusterParams = list(                 # see Seurat::FindClusters
    resolution  = a$lsi_resolution,
    sampleCells = a$lsi_sample_cells,
    n.start     = a$lsi_n_start
  ),
  varFeatures   = a$lsi_var_features,
  dimsToUse     = lsi_dims[1]:lsi_dims[2]
)

archrproj <- addHarmony(
  ArchRProj   = archrproj,
  reducedDims = "IterativeLSI",
  name        = "Harmony",
  groupBy     = a$harmony_group_by
)

archrproj <- addClusters(
  input       = archrproj,
  reducedDims = "IterativeLSI",
  method      = "Seurat",
  name        = "Clusters",
  resolution  = a$cluster_resolution
)

archrproj <- addUMAP(
  ArchRProj   = archrproj,
  reducedDims = "IterativeLSI",
  name        = "UMAP",
  nNeighbors  = a$umap_n_neighbors,
  minDist     = a$umap_min_dist,
  metric      = a$umap_metric
)

# --- MACS2 reproducible peak set -------------------------------------------
archrproj <- addGroupCoverages(ArchRProj = archrproj, groupBy = a$group_by)

pathToMacs2 <- if (!is.null(a$macs2_path) && nzchar(a$macs2_path)) a$macs2_path else findMacs2()
archrproj <- addReproduciblePeakSet(
  ArchRProj   = archrproj,
  groupBy     = a$group_by,
  pathToMacs2 = pathToMacs2
)

archr_peaks <- getPeakSet(archrproj)

# --- per-sample feature matrices over the peak set -------------------------
# One peak x cell count matrix per sample, in sample-sheet order (the order
# create_object.R indexes into).
mtx_cluster <- vector("list", nrow(samples))
names(mtx_cluster) <- samples$sample_id
for (i in seq_len(nrow(samples))) {
  frag <- CreateFragmentObject(inputFiles[i], cells = barcodeList[[i]])
  mtx_cluster[[i]] <- FeatureMatrix(
    fragments = frag,
    features  = archr_peaks,
    process_n = a$feature_matrix_process_n
  )
}

# --- save ------------------------------------------------------------------
saveArchRProject(ArchRProj = archrproj, outputDirectory = a$proj_out_dir, load = FALSE)
saveRDS(archr_peaks, cfg$paths$archr_peaks)
saveRDS(mtx_cluster, cfg$paths$feature_matrix)
message("Wrote ", cfg$paths$archr_peaks, " and ", cfg$paths$feature_matrix)
