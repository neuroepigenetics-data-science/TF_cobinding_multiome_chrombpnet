pi# Running the Zamboni et al. pipeline locally

A practical, step-by-step guide to what this repository does, the software each
step needs, and its inputs/outputs — written so you can (a) reproduce the paper's
analysis, (b) **skip preprocessing and start from published processed objects**,
or (c) adapt the pipeline to your **own multiome data**.

Paper: *Decoding injury-responsive enhancers in the CNS for cell state targeting.*
Data: single-nucleus **multiome** (paired snRNA-seq + snATAC-seq) of mouse spinal
cord across injury timepoints (Uninjured, 1/3/7/28 dpi), plus CUT&Tag and an
enhancer-AAV reporter assay. Genome: **mm10 (mouse)**.

> This is a collection of analysis scripts, **not a turnkey pipeline**. Paths are
> hard-coded to the authors' HPC cluster (e.g. `/cfs/klemming/...`, `/proj/...`).
> You will edit paths in every script. Run scripts interactively/section by
> section rather than expecting `Rscript script.R` to work end-to-end.

---

## Project status — updated 2026-07-21

Living log of **our** local reproduction/adaptation effort. (The reference guide in
Sections 0–6 below describes the original repo generically; this section records
where *we* actually are. Where the two disagree, this section is newer.)

### ✅ Done
- **Got the merged, annotated handoff object from the author** (Margherita Zamboni,
  2026-07-17): `multiome_obj_clean_240326.rds` (7.2 GB) at the repo root. This *is*
  the ★ `merged_multiome.rds` of Section 0 (same object referenced elsewhere as
  `multiome_clean_240925.rds`) — 67,072 cells, 21 samples, `cluster_ids`
  cell-type annotations already present.
- **Got all 21 per-sample `atac_fragments.tsv.gz` from the author's Dropbox**
  (2026-07-21) → `data/<SAMPLE>/outs/atac_fragments.tsv.gz` (~21 GB). This is the
  file GEO does *not* provide and that Section 4 flagged as the ChromBPNet blocker.
  **That blocker is now cleared.** The shared folder held 31 samples (author's
  extras); we pulled only the 21 in `pipeline/samples.csv`.
- **Regenerated the `.tbi` index for all 21** fragment files (`tabix -p bed`).
- **Built the R environment (Tier 1)** — isolated conda env **`zamboni-r`**: R 4.3.3,
  Seurat 5.3.0, Signac 1.14.0, + `EnsDb.Mmusculus.v79`, `BSgenome.Mmusculus.UCSC.mm10`,
  `JASPAR2020`, `TFBSTools`, `DropletUtils`, `scDblFinder`, dplyr/ggplot2/patchwork/Matrix.
  Whole `source.R` stack loads **except ArchR (intentionally deferred)**. System R 4.2.1
  left untouched. Run scripts with `conda run -n zamboni-r Rscript ...`.
- **Started the de-hardcoding refactor** in `pipeline/` — see `pipeline/README.md`.
  `create_object.R` is fully parameterised; config + sample sheet + utils are in place.

### 🔎 Checked / verified
- Object inspected **without Seurat installed** (via `attr()` slot access): Seurat v5,
  67,072 cells, `cluster_ids` with 0 NAs across 10 cell types, 21 samples, assays
  `peaks`/`RNA_raw`/`RNA`/`chromvar`, reductions pca/lsi/wnn.umap/rna.umap/lsi.umap.
- **All 21 fragment files md5-verified byte-identical to the author's originals** —
  each file's md5 matches the hash stored in that sample's Signac `Fragment` slot
  inside the object (21/21 OK). `gzip -t` passes on all 21; BGZF confirmed; a
  `tabix` random-access query returns real fragment records.
- `pipeline/` validated earlier without data: config loads, sample sheet parses,
  `best_pK` matches the original, scripts parse.

### 🚧 To do  *(revised 2026-07-30 — the Track-A items below are all now DONE)*
- ~~Set up the R environment~~ ✅ **DONE.** conda env `zamboni-r`; **ArchR also
  installed** (1.0.3), though Track A never needs it.
- ~~Wire fragments into `SplitFragments`~~ ✅ **DONE.** Note the gotcha:
  `UpdatePath()` does **not** work here — Signac validates the `.tbi` hash too and
  ours are locally regenerated. See `pipeline/R/prepare_ml_input.R` header.
- ~~Finish the refactor + Snakemake~~ ✅ **DONE.** All four Track-B scripts are
  config-driven, `Snakefile` orchestrates them, and Track A added four more
  scripts in `pipeline/`.
- ~~ChromBPNet prep~~ ✅ **DONE**, including the injury DARs, the DAR BEDs, the
  per-cell-type MACS2 peaks, and the DAR×peak intersect that no published script
  produces.
- ~~mm10 `genome.fa`~~ ✅ **DONE** — `meta/mm10.fa` + `.fai`, md5-verified
  (`db005b65828db31735f384e4c5787be5`).
- **Still needed for ML: a GPU.** That is the sole blocker. Also still needed for
  **SCENIC+** (a separate branch, untouched): the cisTarget `.feather` motif
  databases.
- **Before any transfer to a GPU host:** bgzip the fragment BEDs (~43 GB raw,
  Microglia alone 12.8 GB) and sort/bgzip/tabix the 6.0 GB bias file.
- **Not downloaded (only needed to reproduce preprocessing from scratch):** per-sample
  cellbender `.h5` and `per_barcode_metrics.csv` — actually these *were* later
  pulled and are present under `data/<SAMPLE>/outs/`. Track B is runnable, it is
  simply **not on the critical path**: we already hold the merged object, and
  `run_chrombpnet.sh` re-calls peaks with MACS2 regardless. Track B matters only
  to reproduce the paper or to process **our own** data (where we cannot inherit
  Zamboni's peak set or annotations).

### ▶️ How to proceed

> ## ✅ TRACK A IS COMPLETE (2026-07-30) — GPU access is the only blocker.
>
> Every non-GPU step is done and containerised. **Read `2026-07-30.txt` first**;
> it has the full session log, the measured results, and a "RESUME HERE" list.
> Do not rebuild these — they are on disk (`ML/` is 125 GB):
>
> | output | contents |
> |---|---|
> | `markers/` | injury DARs (peaks) + RNA (81,161 DEGs) |
> | `ML/celltype_fragments/` + `_injury/` | fragment BEDs, bias input, 21 + 20 files |
> | `ML/celltype_dars/` | 226,614 unique regions |
> | `ML/dars_intersect_peaks/` | 190,742 regions — ChromBPNet's contribs input |
> | `ML/macs2_celltype_peaks/`, `ML/filtered_peaks/` | per-cell-type peaks |
> | `meta/mm10.fa` + `.fai` | md5-verified vs UCSC |
>
> Reproduce any of it with `./run_tracka.sh all` (image `zamboni-tracka:1.0`,
> verified to regenerate host outputs byte-for-byte).
>
> **Next:** GPU. Local training is impossible — Docker on macOS cannot pass
> through the Metal GPU and Apptainer has no macOS build. Use HPC via
> `apptainer build chrombpnet.sif docker://kundajelab/chrombpnet:latest` +
> `apptainer exec --nv`, a rented cloud GPU, or Colab for one cell type.
> bgzip the ~43 GB of fragment BEDs before any transfer.

Two independent tracks; pick per goal.

- **Track A — reach the science fastest. ✅ DONE through the pre-GPU steps.** We
  already had the ★ object + fragments, so all of Section 2 was skipped. Order was:
  (1) set up the R env; (2) injury DARs → `markers_Uvsothers_peaks_all.rds`;
  (3) `prepare_ml_input.R` → per-cell-type fragment BEDs; (3b) DAR BEDs;
  (1 of run_chrombpnet.sh) per-cell-type MACS2 peaks + the DAR×peak intersect;
  (4) `run_chrombpnet.sh` on a GPU box ← **only this remains**.
  The refactored, config-driven versions live in `pipeline/`
  (`differential_analysis_injury.R`, `prepare_ml_input.R`, `prepare_ml_dars.R`,
  `call_celltype_peaks.sh`) — not the originals in `code/`, which have bugs
  documented in each script's header. SCENIC+ can still start in parallel from
  the object.
- **Track B — reproducibility / adapt to our own data.** Continue the `pipeline/`
  refactor (`filtering` → `merge` → `archr`), then Snakemake. Verify each step's
  output against the author's object. Cell Ranger ARC itself is **not runnable on
  this 24 GB machine** — that stays a cluster job.

**Local compute reality:** 24 GB RAM, no GPU, ~660 GB free disk. Fine for inspecting
the object, `tabix`/fragment ops, and light R once the env exists. Cell Ranger ARC
(~64 GB) and ChromBPNet training (GPU) belong on a cluster.

**Artifacts from this work:** downloads at `data/<SAMPLE>/outs/`, download log at
`data/_download.log`; helper scripts (`download_fragments.sh`, `extract_hashes.R`,
`finalize_fragments.sh`) and the expected-md5 list are in the session scratchpad.

### R environment — why each package, and the plan

The three heavy R packages are **not** interchangeable, and only two are on our
critical path:

- **`Seurat` — required.** The handoff object *is* a Seurat v5 object; without
  `Seurat`/`SeuratObject` its S4 classes aren't defined (hence the `attr()`-only
  peeking we've done so far). Needed to actually operate on it — `subset`,
  `Idents()`, `FindMarkers`, normalize/PCA, plotting. `differential_analysis.R`
  (which makes `markers_Uvsothers_peaks_all.rds` for the ML branch) is pure Seurat.
- **`Signac` — required; unblocks the immediate step.** ATAC extension of Seurat.
  Defines the `ChromatinAssay`/`Fragment` classes the `peaks` assay is built from,
  and provides **`SplitFragments`** — the `prepare_ml_input.R` function that splits
  our 21 fragment files by `cluster_ids` into the per-cell-type BEDs ChromBPNet
  trains on. This is the whole point of having chased down the fragments.
- **`ArchR` — NOT on our critical path.** Used in exactly one script,
  `archr_peak_calling.R` (consensus peak calling via MACS2). Its output is already
  baked into the object we received, so ArchR is only needed to *re-call peaks*
  (Track B / own data), never for downstream or ML.

**Two gotchas to remember:**
1. **`code/source.R` does `library(ArchR)` unconditionally**, and every script
   sources it — so a missing ArchR blocks even scripts that never use it. Fix: a
   trimmed preamble (e.g. `source_min.R`) without the ArchR line for Track-A scripts.
   This lets us **defer ArchR** (the single riskiest install) for now.
2. **Fragment paths need repointing before `SplitFragments`.** The object's
   `Fragment` objects still record the author's Dropbox paths
   (`/Users/margherita.zamboni/...`). Use `Signac::UpdatePath()` to repoint them at
   our local `data/<SAMPLE>/outs/...`. Because our files are md5-identical to hers,
   the stored hashes/barcode maps validate cleanly against the repointed paths.

**Install plan — two tiers (chosen approach: isolated conda env, `r-base ≈ 4.3`,
prebuilt binaries; leaves system R 4.2.1 untouched):**

- **Tier 1 (do now — unblocks Track A):** `Seurat` (v5) + `SeuratObject`, `Signac`,
  `EnsDb.Mmusculus.v79`, `BSgenome.Mmusculus.UCSC.mm10`, `JASPAR2020`, `TFBSTools`,
  `DropletUtils`, `scDblFinder`, `DoubletFinder` (GitHub), + `dplyr`/`ggplot2`/
  `patchwork`/`Matrix`. Gets us to `SplitFragments` and `differential_analysis.R`.
- **Tier 2 (defer):** `ArchR` (+ `MACS2`, `cmake`) — GitHub-only, version-finicky,
  weak support on recent R, and unnecessary while we start from the merged object.

Why not just upgrade system R: R 4.2.1 is below Seurat v5's ≥ 4.3 requirement, and
source-compiling Seurat's large dep tree on macOS is slow/failure-prone — conda
binaries sidestep most of it.

---

## 0. The big picture (data flow)

```
RAW SEQUENCING (fastq)
   │  mgi_demultiplex.sh → Cell Ranger ARC (not in repo) → cellbender.sh
   ▼
PER-SAMPLE 10x matrices + ATAC fragments
   │  archr_peak_calling.R → create_object.R → filtering.R → merge_sample_objects.R
   ▼
★ merged_multiome.rds  ◄─── THE KEY HANDOFF FILE (one merged, annotated Seurat object)
   │
   ├─► Explore (differential_analysis.R, etc.)         → marker tables (.rds/.csv), figures
   ├─► Machine learning (prepare_ml_input.R → ChromBPNet) → per-cell-type deep-learning models + motifs
   └─► SCENIC+ (pycistopic → scenicplus notebooks)      → eRegulons / gene-regulatory networks
```

If you only want the downstream analyses, you need **one file: the merged Seurat
object** (`merged_multiome.rds`). See [Section 4 — Skipping preprocessing](#4-skipping-preprocessing).

---

## 1. Environments to set up

You need **three** separate environments. Install them once.

### A. R environment (Seurat/Signac/ArchR world)
Used by: everything in `code/preprocessing/*.R`, `code/machine_learning/prepare_ml_input.R`, `code/explore/*.R`.

- R ≥ 4.2, Bioconductor.
- All packages are `library()`-loaded at the top of `code/source.R` (every R
  script does `source("../source.R")`). Key ones:
  `Seurat`, `Signac`, `ArchR`, `GenomicRanges`, `EnsDb.Mmusculus.v79`,
  `BSgenome.Mmusculus.UCSC.mm10`, `DropletUtils`, `JASPAR2020`, `TFBSTools`,
  `DoubletFinder`, `scDblFinder`, `TCseq`, `gprofiler2`, `UpSetR`, `dplyr`,
  `ggplot2`, `Matrix`, `patchwork`.
- ArchR and DoubletFinder install from GitHub (commented install lines are in
  `source.R`). MACS2 must be on your PATH (ArchR calls `findMacs2()`).

### B. `chrombpnet` conda env (deep learning)
Used by: `code/machine_learning/run_chrombpnet.sh`. GPU strongly recommended.
The env recipe is at the top of that script:
```bash
conda create -n chrombpnet python=3.8
conda activate chrombpnet
conda install -y -c conda-forge -c bioconda samtools bedtools ucsc-bedgraphtobigwig pybigwig meme
pip install chrombpnet
```

### C. `scenicplus` conda env (gene regulatory networks)
Used by: `pycistopic_multiome.ipynb`, `multiome_scenicplus.ipynb`, `config.yaml`.
Install [SCENIC+](https://scenicplus.readthedocs.io/) (which brings `pycisTopic`,
`pycistarget`, `scanpy`, `mudata`). A fourth mini-env with `modisco-lite` +
`h5py` is used by `peak_interpretation.ipynb` to read ChromBPNet motif outputs.

### Reference/meta files
- **Provided in this repo** (`meta/`): `motifs.meme.txt`, `mm10.chr.sizes`,
  `mm10-blacklist.v2.bed`, `splits/fold_0.json … fold_4.json` (chromosome
  train/val/test folds for ChromBPNet).
- **You must download separately:**
  - `genome.fa` — the mm10 genome FASTA (UCSC). Needed by ChromBPNet and CUT&Tag.
  - Cell Ranger ARC `refdata` for mm10.
  - SCENIC+ motif databases (cisTarget rankings/scores `.feather` + motif
    annotation `.tbl`) — the `ctx_db_fname` / `dem_db_fname` / `path_to_motif_annotations`
    entries in `config.yaml`.

---

## 2. Preprocessing — raw reads → merged Seurat object

Run these in order. All `.R` scripts assume a working dir containing a `data/`
folder with one subfolder per sample, each holding Cell Ranger ARC `outs/`.

| # | Script | Software | Input | Output |
|---|--------|----------|-------|--------|
| 1 | `mgi_demultiplex.sh` | `pigz`, `seqtk`, `deML` | Raw MGI/DNB fastq | Demultiplexed per-sample fastq (Illumina-style names). **Skip if your data is already demultiplexed Illumina fastq.** |
| — | *(Cell Ranger ARC — not in repo)* | `cellranger-arc count` | Demuxed fastq + mm10 ref | Per-sample `outs/`: `filtered_feature_bc_matrix`, `gex_raw_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`, `per_barcode_metrics.csv`. **You run this yourself.** |
| 2 | `cellbender.sh` | `cellbender` (GPU) | `gex_raw_feature_bc_matrix.h5` | `cellbended_gex_matrix_fpr001_seurat.h5` (ambient-RNA-corrected counts) |
| 3 | `archr_peak_calling.R` | R: ArchR + MACS2 | `atac_fragments.tsv.gz` + `per_barcode_metrics.csv` per sample | `archr_peaks_multiome.rds` (consensus peak set), `featurematrix_archr_clusterpeaks_multiome.rds` (peak×cell counts), saved ArchR project |
| 4 | `create_object.R` | R: Seurat + Signac | feature matrix (step 3) + `filtered_feature_bc_matrix` + cellbender `.h5` | `unmerged_obj_multiome.rds` (list of per-sample Seurat objects, RNA+ATAC, QC'd, per-sample UMAP/clusters) |
| 5 | `filtering.R` | R: DoubletFinder, scDblFinder | `unmerged_obj_multiome.rds` | `unmerged_obj_multiome_clean.rds` (outliers removed, doublets flagged) |
| 6 | `merge_sample_objects.R` | R: Seurat + Signac | `unmerged_obj_multiome_clean.rds` | **`merged_multiome.rds`** ★ — one merged object, doublets removed, WNN UMAP, **cell-type annotations** in `cluster_ids` |

> **The cluster→cell-type renaming in `merge_sample_objects.R` (lines ~57–75) is
> specific to this dataset.** For your own data the numeric clusters will differ —
> you must re-derive marker-based annotations, not reuse these labels.

### CUT&Tag (independent branch, Fig. S11)
`cutandtag_preprocess.sh` (bowtie2, picard, samtools, bedtools, SEACR, MACS2,
deeptools): fastq → aligned/dedup BAM → fragment beds → peaks (SEACR + MACS2) →
bigwigs. Analyzed downstream by `code/explore/cutandtag.R`. Independent of the
multiome object — only needed if you have CUT&Tag data.

---

## 3. Downstream analyses (all start from the merged object)

### 3a. Explore — differential expression/accessibility, motifs, reporter assay
R scripts in `code/explore/`, all `readRDS("merged_multiome.rds")` (or a
dated copy of it) and write marker tables / figures:

| Script | Reads | Produces |
|--------|-------|----------|
| `differential_analysis.R` | `merged_multiome.rds` | `markers_rna_multiome_celltypes.rds`, `markers_peaks_multiome_celltypes.rds`, `*_linkpeaks.rds`, **`markers_Uvsothers_peaks_all.rds`** (needed by ML), etc. |
| `explore_injury_differential_analysis.R` | merged obj + `markers_Uvsothers_*` | injury time-course clustering (TCseq), GO tables `top_go_*.csv` |
| `ml_interpretations.R` / `ml_injury_interpretation.R` | ChromBPNet/MoDISco CSV reports + marker rds | per-region motif tables, `markers/markers_injury_peaks.rds` |
| `explore_scenicplus.R` | merged obj + `eRegulons_direct.tsv` | SCENIC+ eRegulon exploration |
| `enhancer_aav_reporter_assay.R` | separate AAV-reporter Cell Ranger `outs/` + `Enhancer_coordinates_BCs.csv` | `astros_aav_obj_240625.rds`, enhancer beds |
| `cutandtag.R` | CUT&Tag QC/peak tables from `cutandtag_preprocess.sh` | CUT&Tag figures |

### 3b. Machine learning — ChromBPNet per-cell-type models
1. **`prepare_ml_input.R`** (R). Reads the merged object
   (`objects/multiome_clean_240925.rds`) + `markers/markers_Uvsothers_peaks_all.rds`.
   Uses `Signac::SplitFragments` to write **per-cell-type fragment BED files**
   (`ML/celltype_fragments/*.bed`) and per-cell-type/injury DAR BED files
   (`*_I_dars.bed`).
2. **`run_chrombpnet.sh`** (chrombpnet env, GPU). Inputs: the fragment beds
   above, `genome.fa`, `meta/mm10.chr.sizes`, `meta/mm10-blacklist.v2.bed`,
   `meta/splits/fold_0.json`, `meta/motifs.meme.txt`. Steps: MACS2 peak calling →
   train a **bias model** → train one **ChromBPNet model per cell type** →
   `pred_bw`/`contribs_bw` (predictions + base-pair contribution scores) →
   `modisco motifs`/`report` (discover motifs, link to TFs). Outputs:
   `models/<celltype>/`, bigwigs, `*_modisco.h5`, MoDISco HTML/CSV reports.
3. **`peak_interpretation.ipynb`** (modisco-lite env). Reads
   `*_contrib.counts_scores.h5` + `*_modisco.h5` to visualize/interpret motifs
   at specific enhancers and count motifs per region.

> Use **`pipeline/run_chrombpnet.sh`**, not `code/machine_learning/run_chrombpnet.sh`.
> The published script is not runnable as written; the corrected version fixes
> 10 defects, documented in its header. Two are silent rather than fatal:
> all cell types wrote into one directory literally named `{i}`, and
> contribution scores were computed from `chrombpnet.h5` instead of
> `chrombpnet_nobias.h5` (which leaves the Tn5 bias in, so MoDISco recovers the
> insertion motif rather than TF motifs).

---

### ⚠️ DAR significance threshold — unresolved, decide before interpreting motifs

**The DAR region sets fed to ChromBPNet are filtered on the raw `p_val`, with no
multiple-testing correction.** This reproduces the paper exactly
(`code/machine_learning/prepare_ml_input.R:53` and
`code/explore/differential_analysis.R:134,137` both use
`p_val < 0.05 & avg_log2FC > 0.5`), so it is *faithful*, not a porting mistake.
We keep it as the default (`dge.filter_p`, `ml.dars_max_p` in
`pipeline/config.yaml`).

Two things make it worth a decision rather than a shrug:

1. **The same repository is inconsistent.** Its *cell-type marker* path filters
   on `p_val_adj` (`differential_analysis.R:46,49`); only the *injury* path uses
   raw `p_val`. ~466,834 peaks are tested per contrast.
2. **The effect is not uniform** — it is catastrophic exactly where the
   contrasts are thin. Under Seurat's `p_val_adj` (Bonferroni, n = 466,834,
   confirmed against the stored values), unique regions surviving
   `p_val_adj < 0.05 & avg_log2FC > 0.5`:

   | cell type | raw p | Bonferroni | kept |
   |---|---|---|---|
   | Microglia | 42,940 | 35,112 | 82% |
   | Oligodendrocytes | 7,037 | 4,928 | 70% |
   | Astrocytes | 19,325 | 11,200 | 58% |
   | Ependymal | 86,999 | 42,301 | 49% |
   | Neurons_V | 6,161 | 269 | 4.4% |
   | Endothelial | 30,740 | 559 | 1.8% |
   | Neurons_D | 3,155 | 43 | 1.4% |
   | **OPCs** | **25,501** | **200** | **0.8%** |
   | Perivascular | 4,481 | 0 | 0% |
   | Macrophages | 275 | 0 | 0% |
   | **total** | **226,614** | **94,612** | |

**The important finding is OPCs.** Our existing QC metric — what fraction of
DARs land inside independently-called per-cell-type MACS2 peaks — gave OPCs a
clean bill of health (89% retention, comparable to Astrocytes). Under
multiple-testing correction it loses 99.2% of its regions. **Peak-retention does
not detect multiple-testing inflation**; the two checks are independent and both
are needed. Ependymal, Microglia, Astrocytes and Oligodendrocytes pass both and
are the defensible cell types for motif claims.

**Blocker if we ever want a middle ground:** `differential_analysis_injury.R`
applies the filter *before* `saveRDS` (line 174 vs 182), so `markers/*.rds`
contains only already-significant rows. A BH-FDR threshold — which would sit
between raw p and Bonferroni — **cannot be recomputed from the saved files**;
computing `p.adjust` on them trivially returns everything. Getting BH requires
re-running the DGE with the filter disabled and saving unfiltered output. That
is a ~1 h re-run, not a re-analysis of what we have.

*Decision still open. Options: (a) keep raw p, report the caveat alongside every
motif result; (b) restrict motif claims to the four cell types that survive
both checks; (c) re-run DGE unfiltered and adopt BH-FDR.*

### 3c. SCENIC+ — enhancer-driven gene regulatory networks
1. **`pycistopic_multiome.ipynb`** — builds the cisTopic object + RNA AnnData +
   candidate region sets from exported ATAC/RNA matrices and fragments.
   Outputs (see `config.yaml` `input_data`): `cistopic_obj.pkl`, `rna_adata.h5ad`,
   `region_sets/`.
2. **SCENIC+ Snakemake workflow** driven by **`config.yaml`** — motif enrichment
   (cisTarget/DEM), TF→gene and region→gene inference, eGRN/eRegulon inference,
   AUCell. Outputs (see `config.yaml` `output_data`): `scplusmdata.h5mu`,
   `eRegulons_direct.tsv`, `AUCell_*.h5mu`, etc.
3. **`multiome_scenicplus.ipynb`** — reads `scplusmdata.h5mu` for regulon
   specificity scores (RSS), dotplots, downstream plots.

---

## 4. Skipping preprocessing

If you don't want to run Section 2, you only need the **outputs** of preprocessing
as inputs to the downstream stages. The public data lives in **GEO under three
sub-series** (raw HPC paths in the scripts, `/cfs/klemming/...`, are the authors'
private copies, not public URLs).

### What GEO actually provides (and what it does *not*)

| GEO accession | Assay | Files deposited | Lets you skip… |
|---------------|-------|-----------------|----------------|
| **GSE304399** | multiome (main) | RNA count matrix (`rna_matrix.mtx.gz` + barcodes + features) and ATAC **peak** count matrix (`atac_matrix.mtx.gz` + barcodes + `atac_peaks.tsv.gz`). RNA and ATAC barcode files are identical → a **matched, already-filtered** cell set. | Raw-fastq/demux, Cell Ranger, cellbender, **QC/doublet filtering, and peak calling** for the multiome — the deposited cells are the clean, matched set (~the cells in `merged_multiome.rds`). These `.mtx`/`.tsv` are exactly what `pycistopic_multiome.ipynb` ingests → **SCENIC+ can start from GEO** (once cells are annotated). |
| **GSE304196** | CUT&Tag | `cutandtag_merged_fragments_sorted.tsv.gz` | **All of `cutandtag_preprocess.sh`** — this file *is* its final output. Go straight to `cutandtag.R`. |
| **GSE304349** | enhancer-AAV | injured + uninjured `matrix.mtx.gz` / `barcodes` / `features` | Cell Ranger for the AAV assay (the script reads a raw `.h5`; GEO gives filtered `.mtx` — minor adaptation). |

> **Superseded for this project (2026-07-21):** we obtained both the merged
> annotated object *and* all 21 `atac_fragments.tsv.gz` directly from the author.
> The GEO limitations below still describe what a fresh public-data start faces,
> but we are no longer blocked by them — see **Project status** at the top.

**Important — GEO does NOT let you fully skip preprocessing for the core analyses:**

1. **No cell-type annotations / no merged object.** There is no
   `merged_multiome.rds`, `.h5ad`, or metadata file on GEO — the barcode files are
   plain barcode lists with **no `cluster_ids`/`merged_id` column**. Every
   `explore/` and ML script, and SCENIC+'s `merged_id` grouping, keys off those
   labels. So you must build a Signac object from the `.mtx`, **re-cluster, and
   manually re-annotate cell types** (the hard-coded cluster→celltype map in
   `merge_sample_objects.R` will not match a freshly clustered object). Also check
   whether sample/timepoint is even recoverable — likely only from a prefix in
   `atac_barcodes.tsv.gz`, since there is no metadata column.
2. **No fragment files — ML branch is a hard blocker.** ChromBPNet needs
   per-cell-type ATAC **fragment** files (`prepare_ml_input.R` `SplitFragments`).
   GEO ships only a peak×cell **count matrix** (`atac_matrix.mtx.gz`), and a count
   matrix **cannot be converted back into fragments**. No `atac_fragments.tsv.gz`
   is deposited at the series level (and none per-GSM in this series). Without the
   authors' fragment files, the **ML branch cannot be reproduced from public data
   at all.**

**Net:** GEO gives a clean, filtered, matched multiome **count matrix** — enough
to rebuild a Signac object and skip fastq → Cell Ranger → cellbender → QC →
doublet → peak-calling. But you must **re-cluster + re-annotate cell types**
before differential analysis or SCENIC+, and the **ChromBPNet branch is not
reproducible without the authors' fragment files**. CUT&Tag (GSE304196) and the
AAV assay (GSE304349) are fully covered by their deposits.

### Can you *reproduce* `merged_multiome.rds` and the fragment files?

**Raw reads are also public** — GSE304399 links to SRA under BioProject
**PRJNA1301245**. So yes, both are reproducible, with an important distinction:

- **`merged_multiome.rds` — reproducible (functionally, not byte-identical).**
  - *From the GEO matrices (fast):* build a Signac object directly from the
    `.mtx` files (they're the clean, matched cell set), then re-run
    normalize → PCA/LSI → WNN → clustering → **manual cell-type annotation**.
    Equivalent for all count-based analyses (DE, DA, marker peaks, SCENIC+).
  - *From raw fastq (complete):* run all of Section 2.
  - Caveats: (1) not byte-identical — LSI subsampling, UMAP, clustering,
    DoubletFinder, Harmony are stochastic (scripts set `set.seed(1234)`, but
    tool versions still shift results); (2) cell-type labels must be
    **re-derived manually** — the hard-coded cluster→celltype map is tied to
    the authors' cluster numbering. Canonical markers make the labels
    reproducible in practice, but it is a manual step.

- **Per-cell-type fragment files — reproducible ONLY from raw fastq.**
  `SplitFragments` partitions a raw `atac_fragments.tsv.gz` by cell type. You
  **cannot** derive fragments from the GEO peak×cell count matrix (lossy
  aggregation — per-read coordinates are gone). The only route:
  1. Download raw fastq from SRA (**PRJNA1301245**).
  2. Run **Cell Ranger ARC** per sample → `atac_fragments.tsv.gz`.
  3. Rebuild + annotate the merged object.
  4. `prepare_ml_input.R` → `SplitFragments` → per-cell-type BEDs → ChromBPNet.
  This reproduces everything including the ML branch — it is just the full,
  compute-heavy Section 2 (Cell Ranger ARC ×21 samples, cellbender on GPU).

> **Recommended:** if you want to *reproduce their result* rather than *re-run
> the whole pipeline*, email the corresponding author for the annotated
> `merged_multiome.rds` and the per-cell-type `atac_fragments.tsv.gz` (or
> `ML/celltype_fragments/*.bed`) — far cheaper than re-running Cell Ranger ARC on
> all 21 samples. Otherwise, everything is reproducible from raw fastq
> (PRJNA1301245).

If the authors later deposit the processed `.rds`/`.h5mu` objects (Zenodo/Figshare
or a GSM supplementary), point the `readRDS()` calls at your copy and skip all of
Section 2.

**Minimum file to start each downstream stage:**

| To run… | You need (produced by preprocessing) | Notes |
|---------|----------------------------------------|-------|
| Explore (`differential_analysis.R` and everything depending on it) | **`merged_multiome.rds`** — the merged, annotated Seurat object | This single object unlocks almost all `explore/` scripts. In ML/explore scripts it appears under dated names (`multiome_clean_240925.rds`, `multiome_obj_clean_240326.rds`) — it's the same object; just point the `readRDS()` at your copy. |
| Machine learning (`prepare_ml_input.R` → ChromBPNet) | `merged_multiome.rds` **+** `markers/markers_Uvsothers_peaks_all.rds` | The markers file is itself produced by `differential_analysis.R`; either grab it from the deposit or regenerate it from the merged object. |
| ChromBPNet training only (skip R entirely) | Per-cell-type fragment BEDs (`ML/celltype_fragments/*.bed`) + DAR BEDs | These are the outputs of `prepare_ml_input.R`. With them + `genome.fa` + the `meta/` files you can run `run_chrombpnet.sh` directly. |
| SCENIC+ notebooks/plots | `cistopic_obj.pkl`, `rna_adata.h5ad`, and/or the final `scplusmdata.h5mu` + `eRegulons_direct.tsv` | Depends how far along you start; `multiome_scenicplus.ipynb` needs only `scplusmdata.h5mu`. |
| Enhancer-AAV analysis | Its own Cell Ranger `outs/` + `Enhancer_coordinates_BCs.csv` | Independent of the multiome object. |
| CUT&Tag analysis | QC/peak tables from `cutandtag_preprocess.sh` | Independent branch. |

**Bottom line:** get `merged_multiome.rds` and you can skip *all* of Section 2 and
go straight to differential analysis, ML input prep, and SCENIC+ exploration.

---

## 5. Adapting to your own data

1. **Start from raw fastq → Cell Ranger ARC** to get per-sample `outs/`. You
   cannot skip preprocessing for your own samples — the merged object is
   dataset-specific.
2. Lay out `data/<sample>/outs/...` exactly as the R scripts expect
   (`archr_peak_calling.R` globs `data/*/outs/atac_fragments.tsv.gz` and
   `per_barcode_metrics.csv`).
3. Run Section 2 steps 2→6. **Re-do cell-type annotation** in
   `merge_sample_objects.R` for your clusters — don't reuse the hard-coded
   `RenameIdents` labels.
4. Adjust anything mm10/mouse-specific if your data isn't mouse: `source.R`
   annotations (`EnsDb.Mmusculus.v79`, `BSgenome.Mmusculus.UCSC.mm10`, JASPAR
   `tax_group='vertebrates'`), the `meta/` genome/blacklist/chrom-size files,
   `genome.fa`, and the `species`/`biomart_host` fields in `config.yaml`.
5. Fix hard-coded paths and dated filenames in each script before running.

---

## 6. Notes for future sessions (e.g. Opus 4.8)

- `code/source.R` is the shared R preamble; every `.R` script sources it. Look
  there for the loaded packages, mm10 gene annotation, JASPAR motif set, plotting
  palette, and the human→mouse TF-name dictionary used for MoDISco reports.
- Scripts are **not parameterized** — inputs/outputs are literal filenames inside
  each script. When wiring the pipeline, grep each script for
  `readRDS` / `saveRDS` / `Read10X` / `read.csv` / `write.table` to see its exact
  I/O contract (that's how this guide's tables were built).
- The same merged object is referenced under several dated names across scripts;
  treat them as one file and repoint `readRDS()`.
- Three compute stacks: R (Seurat/Signac/ArchR), `chrombpnet` conda env (GPU),
  `scenicplus` conda env. Keep them separate.
- `meta/` already ships the ChromBPNet folds, mm10 chrom sizes, blacklist, and
  MEME motifs; the mm10 `genome.fa` and SCENIC+ `.feather` motif databases are
  the big external downloads.
