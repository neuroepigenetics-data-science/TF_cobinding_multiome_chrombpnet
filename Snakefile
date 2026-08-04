# ============================================================================
# Snakefile — orchestration for the R preprocessing chain
# ----------------------------------------------------------------------------
#   archr_peak_calling.R -> create_object.R -> filtering.R -> merge_sample_objects.R
#
# The four scripts are already config-driven (pipeline/config.yaml +
# pipeline/samples.csv); this layer only sequences them, declares their real
# inputs/outputs so a failed run resumes from the step that broke instead of
# redoing ArchR, and keeps per-step logs.
#
# Run from the REPO ROOT:
#
#   snakemake -j1 -n                                          # dry run (full, 21 samples)
#   snakemake -j1                                             # full run
#   snakemake -j1 --config cfg=pipeline/config.smoke.yaml     # 2-sample smoke run
#   snakemake -j1 --until create_object                       # stop after a step
#   snakemake -j1 --rerun-incomplete                          # after an interrupted run
#
# Overrides (all optional, via --config key=value):
#   cfg=<path>       pipeline config to use          (default pipeline/config.yaml)
#   rscript=<cmd>    how to invoke Rscript           (default: conda run in zamboni-r)
#
# NOTE: -j1 is deliberate in the examples. The steps are strictly sequential, and
# ArchR already parallelises internally via config archr.threads.
# ============================================================================

import csv
import os
import yaml

# ---- config ---------------------------------------------------------------
# Read the pipeline config ourselves rather than via Snakemake's `configfile:`,
# so we hold the config *path* (the R scripts take it as argv[1]) and its
# *contents* (output paths, sample sheet) from one unambiguous source.
CONFIG_PATH = config.get("cfg", "pipeline/config.yaml")
if not os.path.exists(CONFIG_PATH):
    raise WorkflowError(f"pipeline config not found: {CONFIG_PATH}")

with open(CONFIG_PATH) as fh:
    CFG = yaml.safe_load(fh)

# `conda run --no-capture-output` streams the R output live; without it conda
# buffers everything until the process exits, which is painful on multi-hour steps.
RSCRIPT = config.get("rscript", "conda run --no-capture-output -n zamboni-r Rscript")

PATHS = CFG["paths"]
OUTS = CFG["cellranger_outs"]

# Logs land next to the run's outputs, so a smoke run doesn't overwrite full-run logs.
LOGDIR = os.path.join(os.path.dirname(PATHS["merged_obj"]) or ".", "logs")


# ---- sample sheet ---------------------------------------------------------
with open(CFG["project"]["samples"], newline="") as fh:
    SAMPLES = list(csv.DictReader(fh))

if not SAMPLES:
    raise WorkflowError(f"sample sheet is empty: {CFG['project']['samples']}")

SAMPLE_IDS = [s["sample_id"] for s in SAMPLES]
SAMPLE_PATH = {s["sample_id"]: s["path"] for s in SAMPLES}


def outs_path(sample_id, which):
    """Mirror of outs_path() in pipeline/R/pipeline_utils.R."""
    return os.path.join(SAMPLE_PATH[sample_id], OUTS[which])


def barcode_metrics(sample_id):
    """First existing barcode-metrics candidate, mirroring valid_barcodes() in
    pipeline_utils.R (Cell Ranger ARC per_barcode_metrics.csv vs older singlecell.csv)."""
    for bm in CFG["barcode_metrics"]:
        candidate = os.path.join(SAMPLE_PATH[sample_id], bm["file"])
        if os.path.exists(candidate):
            return candidate
    raise WorkflowError(
        f"no barcode-metrics file for {sample_id}; tried: "
        + ", ".join(bm["file"] for bm in CFG["barcode_metrics"])
    )


def filtered_matrix_files(sample_id):
    """Read10X() needs all three files; list them individually so Snakemake
    notices a partial/interrupted download rather than just seeing a directory."""
    d = outs_path(sample_id, "filtered_matrix")
    return [os.path.join(d, f) for f in
            ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")]


# Files every step genuinely depends on: the sample sheet and the shared R code.
# The config file itself is deliberately NOT listed here. If it were, editing an
# unrelated section (say `merge:`) would invalidate archr_peak_calling and force
# a multi-hour re-run of ArchR. Instead each rule declares the config sections it
# actually reads under `params:`, so Snakemake's `params` rerun-trigger
# invalidates precisely the steps a given edit affects.
COMMON = [CFG["project"]["samples"],
          "pipeline/R/pipeline_utils.R", CFG["project"]["source_r"]]


# ---- targets --------------------------------------------------------------
rule all:
    input:
        PATHS["merged_obj"]


# ---- 1. ArchR peak calling ------------------------------------------------
# Heaviest step: arrow files from fragments, iterative LSI, Harmony, clustering,
# MACS2 reproducible peak set, then per-sample counts over that peak set.
# Side effects not declared as outputs (ArchR manages them itself, and declaring
# them would make Snakemake delete them between attempts): <sample>.arrow in the
# working directory, ArchRLogs/, and the ArchRProject dirs under config archr.*
rule archr_peak_calling:
    input:
        COMMON,
        script="pipeline/R/archr_peak_calling.R",
        fragments=[outs_path(s, "atac_fragments") for s in SAMPLE_IDS],
        barcodes=[barcode_metrics(s) for s in SAMPLE_IDS],
    output:
        peaks=PATHS["archr_peaks"],
        feature_matrix=PATHS["feature_matrix"],
    params:
        archr=CFG["archr"],
        barcode_metrics=CFG["barcode_metrics"],
        cellranger_outs=OUTS,
    threads: CFG["archr"]["threads"]
    log:
        os.path.join(LOGDIR, "01_archr.log")
    shell:
        "{RSCRIPT} {input.script} {CONFIG_PATH} > {log} 2>&1"


# ---- 2. per-sample Seurat objects -----------------------------------------
rule create_object:
    input:
        COMMON,
        script="pipeline/R/create_object.R",
        feature_matrix=PATHS["feature_matrix"],
        fragments=[outs_path(s, "atac_fragments") for s in SAMPLE_IDS],
        cellbender=[outs_path(s, "cellbender_h5") for s in SAMPLE_IDS],
        filtered=[f for s in SAMPLE_IDS for f in filtered_matrix_files(s)],
    output:
        unmerged=PATHS["unmerged_obj"],
    params:
        create_object=CFG["create_object"],
        assays=CFG["assays"],
        genome=CFG["genome"],
        gene_name_fixes=CFG["gene_name_fixes"],
        cellranger_outs=OUTS,
    log:
        os.path.join(LOGDIR, "02_create_object.log")
    shell:
        "{RSCRIPT} {input.script} {CONFIG_PATH} > {log} 2>&1"


# ---- 3. QC filtering + doublet removal ------------------------------------
rule filtering:
    input:
        COMMON,
        script="pipeline/R/filtering.R",
        unmerged=PATHS["unmerged_obj"],
    output:
        cb_clean=PATHS["unmerged_obj_cb_clean"],
        clean=PATHS["unmerged_obj_clean"],
    params:
        filtering=CFG["filtering"],
        assays=CFG["assays"],
    log:
        os.path.join(LOGDIR, "03_filtering.log")
    shell:
        "{RSCRIPT} {input.script} {CONFIG_PATH} > {log} 2>&1"


# ---- 4. merge + annotate --------------------------------------------------
# Produces the critical handoff object everything downstream depends on.
rule merge_sample_objects:
    input:
        COMMON,
        script="pipeline/R/merge_sample_objects.R",
        clean=PATHS["unmerged_obj_clean"],
    output:
        merged=PATHS["merged_obj"],
    params:
        merge=CFG["merge"],
        assays=CFG["assays"],
    log:
        os.path.join(LOGDIR, "04_merge.log")
    shell:
        "{RSCRIPT} {input.script} {CONFIG_PATH} > {log} 2>&1"
