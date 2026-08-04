# ============================================================================
# pipeline_utils.R
# ----------------------------------------------------------------------------
# Small helpers shared by the parameterised pipeline scripts. Loading config +
# sample sheet lives here so every step reads the same two files instead of
# carrying hard-coded paths/filenames.
# ============================================================================

suppressPackageStartupMessages(library(yaml))

# Load the YAML config. Default location is pipeline/config.yaml relative to the
# repo root (the expected working directory).
load_config <- function(config_path = "pipeline/config.yaml") {
  if (!file.exists(config_path)) {
    stop("Config not found: ", normalizePath(config_path, mustWork = FALSE),
         "\nRun scripts from the repo root, or pass an explicit path.")
  }
  yaml::read_yaml(config_path)
}

# Load the sample sheet as a data.frame (one row per sample).
load_samples <- function(cfg) {
  samples_path <- cfg$project$samples
  if (!file.exists(samples_path)) {
    stop("Sample sheet not found: ", samples_path)
  }
  df <- read.csv(samples_path, stringsAsFactors = FALSE)
  required <- c("sample_id", "path")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Sample sheet is missing required column(s): ", paste(missing, collapse = ", "))
  }
  df
}

# Build the full path to a Cell Ranger `outs/` artefact for one sample.
# `which` is a key under cellranger_outs in the config
# (e.g. "atac_fragments", "filtered_matrix", "cellbender_h5").
outs_path <- function(cfg, sample_path, which) {
  suffix <- cfg$cellranger_outs[[which]]
  if (is.null(suffix)) {
    stop("Unknown cellranger_outs key: '", which,
         "'. Available: ", paste(names(cfg$cellranger_outs), collapse = ", "))
  }
  file.path(sample_path, suffix)
}

# Return the valid (called-cell) barcodes for one sample. Cell Ranger ARC and
# older Cell Ranger name this file differently and use a different cell-flag
# column, so cfg$barcode_metrics is a list of {file, cell_col} candidates tried
# in order; the first file that exists wins. Rows with cell_col == 1 are kept.
valid_barcodes <- function(cfg, sample_path) {
  for (bm in cfg$barcode_metrics) {
    f <- file.path(sample_path, bm$file)
    if (file.exists(f)) {
      meta <- utils::read.csv(f, row.names = 1)
      return(rownames(meta)[meta[[bm$cell_col]] == 1])
    }
  }
  stop("No barcode-metrics file found for sample at ", sample_path,
       " (tried: ", paste(vapply(cfg$barcode_metrics, `[[`, character(1), "file"),
                          collapse = ", "), ")")
}

# Apply the configured transgene/gene-name fixes to a vector of feature names.
# gene_name_fixes in the config is a named list: "raw name" -> "safe name". might create a problem!!!! reiterating the same names???
fix_gene_names <- function(gene_names, cfg) {
  fixes <- cfg$gene_name_fixes
  if (is.null(fixes)) return(gene_names)
  for (raw in names(fixes)) {
    gene_names <- gsub(pattern = raw, replacement = fixes[[raw]],
                       x = gene_names, fixed = TRUE)
  }
  gene_names
}
