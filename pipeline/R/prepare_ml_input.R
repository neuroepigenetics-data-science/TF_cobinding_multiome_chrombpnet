# ============================================================================
# prepare_ml_input.R  (parameterised)
# ----------------------------------------------------------------------------
# Track A — the ChromBPNet critical path. Splits the 21 ATAC fragment files by
# cell type into the per-cell-type BEDs that run_chrombpnet.sh trains on.
#
# Config-driven version of code/machine_learning/prepare_ml_input.R.
#
# Reads:   config ml.merged_obj (the annotated merged object) + the local
#          per-sample data/<SAMPLE>/outs/atac_fragments.tsv.gz
# Writes:  config ml.celltype_outdir  <celltype>.bed          (per cell type)
#          config ml.celltype_outdir  1_input_bias.bed        (bias model input)
#          config ml.injury_outdir    <celltype>_<U|I>.bed    (per type x injury)
#
# Deliberately does NOT source code/source.R: that preamble does
# library(ArchR) unconditionally, and nothing here needs ArchR. The peak set is
# already baked into the object's `peaks` assay, and run_chrombpnet.sh re-calls
# peaks with MACS2 from these BEDs anyway.
#
# Deliberate departures from code/machine_learning/prepare_ml_input.R:
#   1. Fragment paths are repointed before splitting. The object still records
#      the author's Dropbox paths, so the original script cannot run anywhere
#      but her machine. Each path's <SAMPLE> segment is mapped to our local
#      data/<SAMPLE>/outs/ and the fragment file is md5-validated against the
#      hash stored in the Fragment object. Note this does NOT use
#      Signac::UpdatePath(): that validates the .tbi hash too, and our indexes
#      were regenerated locally (see the block below).
#   2. The injury label is NOT recomputed as ifelse(condition == "Uninj", ...).
#      This object's `condition` levels are 1dpi/28dpi/3dpi/7dpi/U -- "Uninj"
#      never matches, so that line would silently label every uninjured cell as
#      injured. We use the object's precomputed injury column and verify it
#      against `condition`. On a mismatch we WARN and skip only the injury
#      split: the per-cell-type BEDs are the actual ChromBPNet training input
#      and do not depend on the injury labels.
#   3. Output dirs must be empty. SplitFragments(append = TRUE) is required for
#      correctness (it appends across the 21 fragment files), which also means a
#      second run silently doubles every BED. Guarded rather than discovered.
#
# Run from the repo root:
#   Rscript pipeline/R/prepare_ml_input.R [config.yaml] [--check-only]
# ============================================================================

source("pipeline/R/pipeline_utils.R")

args        <- commandArgs(trailingOnly = TRUE)
check_only  <- "--check-only" %in% args
cfg_args    <- setdiff(args, "--check-only")
config_path <- if (length(cfg_args) >= 1) cfg_args[[1]] else "pipeline/config.yaml"
cfg         <- load_config(config_path)
m           <- cfg$ml

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

# --- load ------------------------------------------------------------------
say("loading ", m$merged_obj)
multiome <- readRDS(m$merged_obj)
say("loaded: ", ncol(multiome), " cells, assay '", m$assay, "' with ",
    nrow(multiome[[m$assay]]), " peaks")

stopifnot(m$assay %in% Assays(multiome))
for (col in c(m$celltype_col, m$bias_group_col, "condition")) {
  if (!col %in% colnames(multiome@meta.data)) {
    stop("metadata column not found on the object: ", col)
  }
}

# --- repoint fragment files onto our local copies --------------------------
# The stored paths are the author's; map each one's <SAMPLE> segment onto
# data/<SAMPLE>/outs/ and validate before touching anything.
frags <- Fragments(multiome[[m$assay]])
say("repointing ", length(frags), " fragment objects")

new_paths <- character(length(frags))
for (i in seq_along(frags)) {
  old <- GetFragmentData(frags[[i]], slot = "path")
  # .../<SAMPLE>/outs/atac_fragments.tsv.gz -> <SAMPLE>
  sid  <- basename(dirname(dirname(old)))
  new_paths[i] <- outs_path(cfg, file.path("data", sid), "atac_fragments")
}

missing <- new_paths[!file.exists(new_paths)]
if (length(missing)) {
  stop("local fragment files not found:\n  ", paste(missing, collapse = "\n  "))
}

# Signac::UpdatePath() cannot be used here. It calls ValidateHash(), which
# requires md5sum(fragment file) AND md5sum(.tbi) to match the stored pair --
# but the author's .tbi files were never distributed, so ours were regenerated
# locally with `tabix -p bed` and necessarily hash differently. Repoint the path
# directly and validate the fragment file itself: that is the hash which
# actually proves our data is the data the annotations were derived from. The
# index is a derived artifact, reproducible from the fragment file.
say("validating fragment-file hashes (index hashes are expected to differ)")
idx_missing <- character(0)
bad_hash    <- character(0)
for (i in seq_along(frags)) {
  idx <- paste0(new_paths[i], ".tbi")
  if (!file.exists(idx)) idx_missing <- c(idx_missing, idx)

  stored <- GetFragmentData(frags[[i]], slot = "hash")[1]
  got    <- unname(tools::md5sum(new_paths[i]))
  if (!identical(got, stored)) bad_hash <- c(bad_hash, new_paths[i])

  slot(frags[[i]], name = "path") <- normalizePath(new_paths[i], mustWork = TRUE)
  # Re-stamp the hash pair to our local files. Required: Fragments<- ->
  # AddFragments() calls ValidateHash() again, so leaving the author's index
  # hash in place makes the object unassignable. Safe because the fragment file
  # was just checked against her stored hash above; only the derived index
  # differs, and we record ours so later Signac validation is meaningful.
  slot(frags[[i]], name = "hash") <- unname(tools::md5sum(c(new_paths[i], idx)))
}
if (length(idx_missing)) {
  stop("fragment index (.tbi) missing -- run `tabix -p bed` on:\n  ",
       paste(idx_missing, collapse = "\n  "))
}
if (length(bad_hash)) {
  stop("fragment file md5 does not match the hash stored in the object:\n  ",
       paste(bad_hash, collapse = "\n  "),
       "\nThese are not the files the annotations were derived from.")
}
say("all ", length(frags), " fragment files matched the author's md5")

Fragments(multiome[[m$assay]]) <- NULL
Fragments(multiome[[m$assay]]) <- frags

# --- cell-type + injury labels ---------------------------------------------
say("cell types (", m$celltype_col, "):")
print(table(multiome@meta.data[[m$celltype_col]]))

# Verify the precomputed injury column instead of recomputing it (see header).
# NOTE: this only gates split 3. The per-cell-type BEDs are the actual
# ChromBPNet training input and do not depend on the injury labels at all, so a
# bad injury column must not prevent them from being written.
injury_col <- m$injury_col
injury_ok  <- FALSE
if (!injury_col %in% colnames(multiome@meta.data)) {
  warning("injury column '", injury_col, "' not on the object; skipping the injury split")
} else {
  expected_injury <- ifelse(multiome$condition == "U", "U", "I")
  expected_combo  <- paste0(multiome@meta.data[[m$celltype_col]], "_", expected_injury)
  actual_combo    <- as.character(multiome@meta.data[[injury_col]])
  n_mismatch      <- sum(actual_combo != expected_combo)

  say("injury check: ", n_mismatch, " / ", ncol(multiome),
      " cells disagree with paste0(", m$celltype_col, ", '_', condition=='U' ? U : I)")
  if (n_mismatch > 0) {
    print(head(data.frame(condition = multiome$condition,
                          on_object = actual_combo,
                          expected  = expected_combo)[actual_combo != expected_combo, ], 10))
    warning("'", injury_col, "' does not match the condition-derived labels; ",
            "SKIPPING the injury split. Inspect that column before using it.")
  } else {
    say("injury labels verified against `condition`")
    injury_ok <- TRUE
  }
}

if (check_only) {
  say("--check-only: stopping before SplitFragments (injury_ok = ", injury_ok, ")")
  quit(save = "no", status = 0)
}

# --- output dirs (must be empty; see header note 3) ------------------------
ensure_empty_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(d, pattern = "\\.bed$")
  if (length(existing)) {
    stop("output dir is not empty: ", d, "\n  ", length(existing),
         " .bed file(s) present. SplitFragments appends, so re-running would ",
         "double their contents. Move or delete them first.")
  }
}
ensure_empty_dir(m$celltype_outdir)
if (injury_ok) ensure_empty_dir(m$injury_outdir)

# --- 1. bias-model input: batch 1 (unsorted, all timepoints) ---------------
say("SplitFragments -> bias input (", m$bias_group_col, " == ", m$bias_ident, ")")
SplitFragments(
  multiome,
  assay         = m$assay,
  group.by      = m$bias_group_col,
  idents        = m$bias_ident,
  outdir        = m$celltype_outdir,
  append        = TRUE,
  buffer_length = as.integer(m$buffer_length),
  verbose       = TRUE,
  file.suffix   = m$bias_suffix
)

# --- 2. per-cell-type fragments (the ChromBPNet training input) ------------
say("SplitFragments -> per cell type (", m$celltype_col, ")")
SplitFragments(
  multiome,
  assay         = m$assay,
  group.by      = m$celltype_col,
  outdir        = m$celltype_outdir,
  append        = TRUE,
  buffer_length = as.integer(m$buffer_length),
  verbose       = TRUE
)

# --- 3. per cell type x injury (only if the labels checked out) ------------
if (injury_ok) {
  say("SplitFragments -> per cell type x injury (", injury_col, ")")
  SplitFragments(
    multiome,
    assay         = m$assay,
    group.by      = injury_col,
    outdir        = m$injury_outdir,
    append        = TRUE,
    buffer_length = as.integer(m$buffer_length),
    verbose       = TRUE
  )
} else {
  say("SKIPPED the injury split -- '", injury_col, "' failed verification")
}

say("done. wrote:")
for (d in c(m$celltype_outdir, m$injury_outdir)) {
  for (f in list.files(d, pattern = "\\.bed$", full.names = TRUE)) {
    say("  ", f, "  (", format(file.size(f), big.mark = ","), " bytes)")
  }
}
