# `pipeline/` — parameterised, config-driven version of the scripts

This directory is the **de-hardcoding layer**. The original analysis code in
`code/` is untouched (keep it as the reference implementation). Here, paths,
sample lists, and thresholds live in two editable files instead of being baked
into every script:

- **`config.yaml`** — all paths, references, and per-step parameters.
- **`samples.csv`** — one row per sample (`sample_id, path, experiment,
  condition, batch, best_pK`). Replaces the `list.dirs("data/")` globbing and the
  hard-coded 21-element vectors (cellbender dir list, per-sample `best_pK`, the
  cluster→celltype map).

Scripts read those two files via **`R/pipeline_utils.R`** and are meant to be run
**from the repo root**:

```bash
Rscript pipeline/R/create_object.R                 # uses pipeline/config.yaml
Rscript pipeline/R/create_object.R my_config.yaml  # or an explicit config
```

## Status

**All four preprocessing scripts are refactored (config-driven).** Run order —
`archr_peak_calling.R` runs **first** (it produces `paths.archr_peaks` +
`paths.feature_matrix`, the input to `create_object.R`):

```
archr_peak_calling.R  →  create_object.R  →  filtering.R  →  merge_sample_objects.R
```

| Step | Parameterised here | Run-ready? |
|------|--------------------|-----------|
| archr peak calling | ✅ `R/archr_peak_calling.R` | ArchR install PENDING (see below) |
| create_object | ✅ `R/create_object.R` (worked template) | ✅ inputs + hdf5r ready |
| filtering | ✅ `R/filtering.R` (uses DoubletFinder 2.0.6 API) | ✅ DoubletFinder installed |
| merge | ✅ `R/merge_sample_objects.R` — **cluster→celltype map in `config.yaml:merge.annotation` must be re-derived for new data** | ✅ |

`create_object.R` is the pattern the rest follow: `source()` the utils,
`load_config()` + `load_samples()`, then use `cfg$...` and the sample-sheet
columns in place of literals.

### Data — all local & verified (2026-07-24)

Every per-sample input for the whole chain is under `data/<SAMPLE>/outs/`:
`atac_fragments.tsv.gz` (+`.tbi`), `per_barcode_metrics.csv`,
`cellbended_gex_matrix_fpr001_seurat.h5`, and
`filtered_feature_bc_matrix/{matrix,barcodes,features}`. The merged annotated
object `multiome_obj_clean_240326.rds` (Zamboni's) is at the repo root.

### Environment (conda env `zamboni-r`; run `conda run -n zamboni-r Rscript ...`)

| Package | State |
|---------|-------|
| R 4.3.3, Seurat 5.3.0, Signac 1.14.0, scDblFinder | ✅ |
| hdf5r 1.3.12 | ✅ (Read10X_h5 works) |
| DoubletFinder 2.0.6 | ✅ — **new API**; `filtering.R` updated to `doubletFinder()` / `reuse.pANN=NULL` |
| macs2 2.2.9.1 | ✅ in **separate** env `macs2`; path in `config.yaml:archr.macs2_path` |
| devtools + ArchR 1.0.3 | ✅ (`library(ArchR)` loads; first attempt needed `r-devtools` installed first) |

**All data and all packages are now in place** — `code/source.R` (which does
`library(ArchR)` unconditionally) loads end-to-end, so every Track-A script can be
run. Nothing has been executed as a pipeline yet.

### ✅ Smoke test PASSED end-to-end (2026-07-28)

The chain has now actually run. 2-sample subset via `samples.smoke.csv` +
`config.smoke.yaml` (generated *from* `config.yaml` by a yaml round-trip so it
can't drift; all outputs redirected under `smoke/`). Samples: `ALL_7dpi_2` +
`ALL_1dpi_3` — smallest fragments, and chosen because neither has
`manual_cluster_drop`/`manual_cluster_keep` entries (those cluster IDs are only
meaningful at 21-sample cluster numbering). ≈50 min wall clock.

```
archr  → 110,625 peaks, 10 clusters, median width 501
create → 2,478 cells
filter → 1,433 + 861 cells
merge  → 2,277-cell smoke/merged_multiome.rds, cluster_ids populated, 0 NAs
```

Two refactor fixes moved from *statically validated* to *empirically proven*:
the per-sample `best_pK` fix (visible in DoubletFinder's column names —
`pANN_0.25_0.01_16` vs `pANN_0.25_0.18_6`, where the original's aliasing bug
would have given the second sample `pK=NA`), and `merge`'s sample-sheet-driven
metadata + the `wsnn_res.0.1` column the original never actually computed.

Benign in a 2-sample run: `Cannot find identity 8..15` warnings, because
`merge.annotation.map` is data-specific to the full 21-sample clustering.

### Snakemake orchestration (added 2026-07-28)

`Snakefile` at the repo root sequences the four scripts and makes a failed run
resumable from the step that broke, instead of redoing ArchR. Run from the
**repo root**:

```bash
snakemake -j1 -n                                       # dry run (full, 21 samples)
snakemake -j1                                          # full run
snakemake -j1 --config cfg=pipeline/config.smoke.yaml  # 2-sample smoke run
snakemake -j1 -n -r                                    # dry run, explain WHY each step reruns
snakemake -j1 --until create_object                    # stop after a given step
snakemake -j1 --rerun-incomplete                       # resume after an interrupted run
snakemake --dag | dot -Tpng > dag.png                  # picture of the DAG
```

`-j1` is deliberate: the steps are strictly sequential and ArchR already
parallelises internally via `config archr.threads`.

Per-step logs land beside that run's outputs (`logs/` for the full run,
`smoke/logs/` for the smoke run), so a smoke run can't clobber full-run logs.

**Invalidation is scoped per step.** The config file is deliberately *not* a
declared input of every rule — if it were, editing `merge:` would invalidate
`archr_peak_calling` and force a multi-hour ArchR re-run. Instead each rule
declares the config sections it actually reads under `params:`. Verified
behaviour:

| edit | reruns |
|---|---|
| `merge:` | `merge_sample_objects` only |
| `filtering:` | `filtering` → `merge_sample_objects` |
| `archr:` | all four (correct: everything derives from the peak set) |

⚠️ **Gotcha:** the `params` trigger hashes the config *structure*, so
**reformatting `config.yaml` invalidates everything** even with identical
values — e.g. a round-trip through Python's `yaml.safe_dump`, which sorts keys
alphabetically. Edit values in place. `config.smoke.yaml` is generated from
`config.yaml` by an R `yaml` round-trip; regenerating it is a reformat, so
follow it with `snakemake --touch --config cfg=pipeline/config.smoke.yaml` if
you want to keep existing smoke outputs.

Because the smoke outputs were produced by hand (before the Snakefile existed),
they were adopted into Snakemake's provenance with `--touch`.

### ▶ Where to resume

> ⚠️ **Everything in this file is Track B (reproducibility / our own data).**
> **None of it is on the critical path to ChromBPNet.** If the goal is the ML
> branch, go to `R/prepare_ml_input.R` (below) and skip this chain entirely —
> see `RUNNING_LOCALLY.md` § "How to proceed". On 2026-07-29 a full 21-sample
> ArchR run was started off this section by mistake and abandoned an hour in;
> it advances reproduction, not the science.

**Track A — ChromBPNet (the critical path).** `R/prepare_ml_input.R` reads the
author's already-annotated `multiome_obj_clean_240326.rds` plus the 21 local
fragment files and writes the per-cell-type BEDs `run_chrombpnet.sh` trains on.
It needs **no ArchR** and none of the four steps below: the peak set is already
baked into the object's `peaks` assay, and `run_chrombpnet.sh` re-calls peaks
with MACS2 from the BEDs anyway.

**✅ ALL OF TRACK A IS DONE (2026-07-30).** Everything below has been run; the
outputs are on disk (`ML/` is 125 GB). See `2026-07-30.txt` for the full session
log. Easiest way to re-run any of it is the container — `./run_tracka.sh all`.

```bash
# step 2 — per-cell-type injury DARs (peaks half is what ChromBPNet needs)
#   --only, never `both`: the object is 18.8 GB loaded and DietSeurat trims it
#   to one assay (11.5 GB peaks / 3.8 GB RNA) to fit in 24 GB.
Rscript pipeline/R/differential_analysis_injury.R pipeline/config.yaml --only peaks
Rscript pipeline/R/differential_analysis_injury.R pipeline/config.yaml --only rna

# step 3 — fragment BEDs (~30 min, reads all 21 GB)
Rscript pipeline/R/prepare_ml_input.R pipeline/config.yaml --check-only  # validate only
Rscript pipeline/R/prepare_ml_input.R pipeline/config.yaml               # + write BEDs

# step 3b — DAR BEDs for contribution scoring / MoDISco (seconds)
Rscript pipeline/R/prepare_ml_dars.R pipeline/config.yaml

# step 1 — per-cell-type MACS2 peaks, blacklist filter, and the DAR x peak
#   intersect that run_chrombpnet.sh needs but no published script produces.
#   ~24 min for all 10. Idempotent: existing outputs are skipped.
bash pipeline/call_celltype_peaks.sh              # all 10
bash pipeline/call_celltype_peaks.sh Endothelial  # a subset
```

| output | contents |
|---|---|
| `markers/` | injury DARs (peaks, ~448,600 rows) + RNA (81,161 DEGs) |
| `ML/celltype_fragments/` | 10 cell types + bias input + `_chr` filtered |
| `ML/celltype_fragments_injury/` | 20 (cell type × U/I) |
| `ML/celltype_dars/` | 226,614 unique regions |
| `ML/macs2_celltype_peaks/`, `ML/filtered_peaks/` | 10 each |
| `ML/dars_intersect_peaks/` | **190,742 regions — ChromBPNet's contribs input** |
| `meta/mm10.fa` + `.fai` | md5-verified against UCSC |

Verified along the way: the injury fragment BEDs partition each cell type's total
**byte-exactly**; DAR BEDs are 3-column, 500 bp, duplicate-free; and 8/10 cell
types retain 89–99% of DARs when intersected with independently-called
per-cell-type peaks (two methods agreeing).

⚠️ **Interpretation caveats — read before trusting motif output.**
- **Endothelial: only 28% DAR retention** (vs 89–99% elsewhere). Second-most DARs
  (30,740) from the fewest cells (439), including 9-cell `3dpi` and 30-cell
  `7dpi` contrasts. Perivascular is the control — nearly as few cells (673) but
  81% retention — so it's the thin *contrasts*, not cell count. Treat its motifs
  cautiously.
- **Macrophages: 273 regions.** Uninjured group is 22 cells. Regions are
  trustworthy (99% in peaks) but far too few for MoDISco. Exclude it rather than
  reporting a negative.
- **Ependymal is the best candidate** — richest (86,999) *and* most validated
  (96%), and it's the paper's reactive population. RNA agrees independently.
- 48% of Ependymal/Microglia DAR rows were duplicates across timepoints and are
  now collapsed. Leaving them (as the published script does) makes MoDISco
  double-count seqlets in persistently-open regions. `ml.dars_dedupe: false`
  restores the original behaviour.

### Reproducibility — `docker/`

`zamboni-tracka:1.0` (5.4 GB, amd64) runs all of the above on any machine.
**Proven, not asserted:** the container regenerated all 10 DAR BEDs to a combined
md5 identical to the host's, byte-for-byte. See `docker/README.md` for mounts and
the ≥16 GB Docker Desktop memory requirement.

```bash
./run_tracka.sh check   # validate inputs, no compute
./run_tracka.sh all     # full pipeline
```

⚠️ **Two things to know before using these.**
- **`--only peaks` / `--only rna`, not `both`.** The object is 18.8 GB loaded
  (4 assays + 5 reductions + WNN graphs); `DietSeurat()` to one assay brings it
  to 11.5 GB, which is what makes a 20,816-cell subset fit in 24 GB. Without the
  slim, `FindMarkers` dies with `vector memory exhausted` on Oligodendrocytes.
- **Thin comparisons.** `Macrophages x U = 22` cells and `Endothelial x 3dpi = 9`.
  Macrophages yield only 275 DARs as a result — likely too few for meaningful
  motif discovery. Biologically expected (macrophages infiltrate after injury),
  but don't read those two like the others.

**Done — the DAR/peak intersect.** `run_chrombpnet.sh` reads
`dars_intersect_peaks/<celltype>_intersect_dars.bed`, i.e. these DARs intersected
with that cell type's MACS2 peak set. That step is in neither published script
and needs the per-cell-type peaks, so it lives at the end of
`call_celltype_peaks.sh` (step 4). Contribution scores should only cover regions
the model was trained on. 190,742 regions across the 10 cell types.

It is emitted as **10-column narrowPeak, not 3-column BED**. ChromBPNet parses
`-r` with `NARROWPEAK_SCHEMA = [chr,start,end,1..6,summit]` and extracts each
sequence centred at `start(col2) + summit(col10)`; a 3-column file cannot be
parsed. Verified against `chrombpnet/evaluation/interpret/interpret.py`. Our
regions are uniformly 500 bp so `summit = 250`, and every one of the 190,742
centres was checked to sit >1057 bp from its chromosome edge (the `inputlen`
2114 half-window), so none will overflow.

**GPU access is now the only blocker.** `meta/` has `motifs.meme.txt`,
`mm10.chr.sizes`, the blacklist, `splits/fold_0.json`, and the md5-verified
`mm10.fa` + `.fai`. What remains needs an NVIDIA GPU:

```bash
apptainer build chrombpnet.sif docker://kundajelab/chrombpnet:latest
apptainer exec --nv chrombpnet.sif chrombpnet pipeline ...   # --nv passes the GPU
```

Local training is **impossible, not just slow**: Docker on macOS cannot pass
through the Metal GPU, and Apptainer has no macOS build at all (it needs Linux
kernel namespaces), so a Lima VM would not help either. Use HPC (Apptainer —
clusters disallow Docker's root daemon), a rented cloud GPU, or Colab to validate
one cell type.

**A corrected `pipeline/run_chrombpnet.sh` now exists** — use it, not
`code/machine_learning/run_chrombpnet.sh`, which is not runnable as published.
It is staged so the CPU-only work can be done here and only training needs a GPU:

```bash
bash pipeline/run_chrombpnet.sh prep        # bias-peaks + splits + nonpeaks (CPU)
bash pipeline/run_chrombpnet.sh bias-train  # GPU
bash pipeline/run_chrombpnet.sh train       # GPU
bash pipeline/run_chrombpnet.sh interpret   # GPU
DRY_RUN=1 bash pipeline/run_chrombpnet.sh all   # print commands, run nothing
```

It fixes 10 defects, all listed in its header. The two that matter most are
*silent* — they produce wrong results rather than an error:

- **`{i}` instead of `${i}`** in every interpret/MoDISco output path (published
  lines 92, 97, 99, 100). All 10 cell types wrote into one directory literally
  named `{i}`, each overwriting the last. A full GPU run would finish holding a
  single cell type's results, with nothing in the logs to say so.
- **contribution scores from `chrombpnet.h5`** (line 94). The wiki is explicit:
  *"Please use the `chrombpnet_nobias.h5` model for this."* `chrombpnet.h5` still
  carries the Tn5 bias component, so MoDISco run on its contributions recovers
  the Tn5 insertion motif instead of TF motifs.

The rest: model paths read a nonexistent `_fl0` suffix; `-b small_frag_bias.h5`
names a file nothing produces (it is `<outdir>/models/bias.h5`); the bias
narrowPeak is read from the wrong directory without its `_peaks` infix;
`mm10.chrom.sizes` vs `mm10.chr.sizes`; the stray `A` in
`-n negatives/A${i}_...`; `_chr` peaks used for negatives but non-`_chr` for
training; and dead `_intersect_dars_chr.bed` code.

Two remaining local tasks, both cheap and required anyway:
1. `sort | bgzip | tabix` the 6.0 GB bias fragment file — now `run_chrombpnet.sh
   bias-peaks`, which points `sort -T` at `$SORT_TMP` (defaults to a repo-local
   dir, not `/tmp`, which is too small).
2. bgzip the per-cell-type fragment BEDs. They total ~43 GB (Microglia alone is
   12.8 GB), which is the real friction in any off-machine plan; ChromBPNet
   accepts bgzipped fragments (its own bias step passes `.tsv.gz`).

**Open question, not a bug — the DAR significance threshold.** The regions are
filtered on raw `p_val` with no multiple-testing correction (faithful to the
paper). Under Bonferroni, OPCs keeps 0.8% of its regions, Endothelial 1.8%,
Perivascular and Macrophages 0% — while Ependymal/Microglia/Astrocytes/
Oligodendrocytes keep 49–82%. Note OPCs passes the peak-retention check (89%)
and fails this one, so the two checks are independent. Full table and the
options in `RUNNING_LOCALLY.md` § "DAR significance threshold".

**Track B — reproduction, below.**

1. **Full 21-sample run** with `pipeline/config.yaml`, then compare against
   Zamboni's `multiome_obj_clean_240326.rds` (repo root; 67,072 cells, 10 cell
   types) as ground truth. Watch the `dims: pca 1:N  lsi 2:M` log lines to see
   whether the dimension guard fires (see Known weakness above). The
   cluster→cell-type map will likely need re-deriving.
   *Budget:* ~21 GB of fragments vs 0.5 GB in the smoke run, and merge handles
   ~67k cells on a 24 GB machine — expect many hours and possible RAM pressure
   at the merge step. Launch it detached under `caffeinate` (this machine
   idle-sleeps and has already killed one run), and use miniconda's
   `snakemake` — `conda run -n zamboni-r snakemake` resolves to a broken
   system-Python install.
2. Deferred: the LSI dimensionality improvement (see Known weakness above) —
   only after the faithful full run has been validated.

`archr.force_arrows` (default `false`) makes the arrow step resumable: an
interrupted run reuses the `<sample>.arrow` files it already finished instead of
rebuilding all 21. ArchR stamps `Metadata/Completed` only on an intact arrow and
rebuilds any file missing the stamp, so a partial arrow is never adopted.

### Environment fixes required to make the published code run (2026-07-28)

None of these are refactor bugs; the published scripts hit them identically.

| Break | Fix |
|---|---|
| `TCseq`, `gprofiler2`, `UpSetR`, `rvest` missing — `source.R` `library()`s them unconditionally, so every script died at the preamble | conda-forge/bioconda install |
| `biovizBase` missing — soft dep of `Signac::GetGRangesFromEnsDb()` | `bioconductor-biovizbase` |
| `seqlevelsStyle(x) <- "UCSC"` fails: GenomeInfoDb 1.38.1 hard-codes the retired host `hgdownload.cse.ucsc.edu`, whose cert no longer matches | `setHook(packageEvent("GenomeInfoDb","onLoad"), …)` repointing to `hgdownload.soe.ucsc.edu`, in the env's `Rprofile.site` (setting the option up front doesn't work — `.onLoad` overwrites it) |
| `ArchR::addHarmony()` → `object 'HarmonyMatrix' not found` — harmony ≥1.0 removed it; CRAN's 2.0.5 had been pulled in as an ArchR dep | pin **harmony 0.1.1** from the CRAN archive (conda-forge's oldest, 1.2.0, is already too new) |

(Sanity re-check any time: `conda run -n zamboni-r Rscript -e 'requireNamespace("ArchR")'`.)

## For your own data

1. Edit `samples.csv` — list *your* samples and the folder each one's Cell Ranger
   ARC output lives in.
2. Edit `config.yaml` — set `genome.species`/`mt_pattern`, add any transgenes to
   `gene_name_fixes`, adjust thresholds if needed.
3. Run the steps from the repo root.

## ⚠ Known weakness in the published method — TO IMPROVE (do not fix silently)

**`create_object.R`'s adaptive WNN dimensionality selection is fragile.** It applies
a single stdev cutoff of `2` to *both* PCA and LSI:

```r
sig_pc_pca <- length(which(pca@stdev > 2));  dims 1:sig_pc_pca
sig_pc_lsi <- length(which(lsi@stdev > 2));  dims 2:sig_pc_lsi
```

The two are on very different scales. Signac's `RunSVD` normalises singular values
by `sqrt(n_features - 1)` — **verified byte-identical in every Signac release from
0.2.4 to 1.17.0, so this is not version drift** — leaving LSI stdev around
`6.6 / 2.8 / 2.0 / 1.4 / …`. Typically only components 1–3 clear a cutoff of 2,
while PCA clears it 42–50 times.

Since LSI drops component 1 (sequencing depth) and indexes `2:sig_pc_lsi`,
**`sig_pc_lsi == 2` produces the single dimension `2:2`**, the embedding collapses
from a matrix to a vector, and `FindMultiModalNeighbors` fails inside `L2Norm`
with `dim(X) must have a positive length`. On the 2-sample smoke set, LSI
component 3 landed at **1.97** and **1.63** — i.e. whether the published code runs
at all turns on which side of `2.0` component 3 happens to fall.

**What we did:** added only a *structural floor* (`create_object.min_dims_pca` /
`min_dims_lsi`) that raises the counts when they are too small to index a
≥2-dimension embedding. It is a **no-op on any dataset where the published code
produced a result**, so it cannot silently shift results — either it never fires
(output is faithful) or it fires where the original would have crashed. Per-sample
`sig_pc` values are logged, and a loud warning is emitted whenever it fires.

**What still needs doing:** replace the heuristic with a scale-appropriate
criterion for LSI — e.g. a depth-correlation cutoff (`DepthCor`), or simply the
fixed dims the paper's *own* other scripts use (`filtering.R` → `1:15, 2:15`;
`merge_sample_objects.R` → `1:30, 2:30`). **Deliberately deferred:** changing the
methodology now would confound the reproduction, since any deviation from
`multiome_obj_clean_240326.rds` could no longer be attributed to a single cause.
Revisit only *after* the full 21-sample run has been compared against the
author's object.

## Notes / intentional changes vs the original

- `create_object.R` builds a **named list** (keyed by `sample_id`) and uses a
  single consistent loop index. The original mixed integer and character
  indexing into an unnamed list; the named-list version is equivalent in result
  but avoids that ambiguity. No statistical behaviour changed.
- Nothing here alters analysis logic, algorithms, or seeds — it only moves
  hard-coded values into `config.yaml` / `samples.csv`.
- **`filtering.R`** has two deliberate departures from the published script, both
  documented in the file header:
  - The original referenced peak assays `macs2_peaks` / `cellranger_peaks` that
    `create_object.R` never builds. Both now resolve to the single ATAC assay in
    `config.yaml:assays.peaks`.
  - The original doublet loop did `best_pK <- best_pK[i]`, which overwrites the
    21-element vector with a scalar on the first iteration, giving every later
    sample `pK = NA`. Per-sample `best_pK` is now read from `samples.csv`, which
    both de-hardcodes and fixes that aliasing bug. Per-cluster IQR outlier logic
    was verified identical to the original on synthetic data.
- **`merge_sample_objects.R`** has three deliberate departures, all in the file
  header:
  - `experiment`/`condition`/`batch` are read from `samples.csv` (keyed by each
    cell's `sample`) instead of `strsplit(sample, "_")`, so it survives sample
    names that don't follow the `EXPERIMENT_CONDITION_BATCH` convention.
  - The original annotated off `wsnn_res.0.1` but only ran `FindClusters` at
    resolution 0.2, so that column never existed. The refactor computes every
    resolution in `merge.cluster_resolutions` and annotates off
    `merge.annotation.cluster_col`.
  - The original per-modality `RunUMAP` calls omitted `reduction=`, so both
    defaulted to `pca` (the ATAC-only UMAP came from pca, not lsi). The refactor
    uses `pca` for the RNA-only UMAP and `lsi` for the ATAC-only UMAP.
  - **The cluster→cell-type map (`merge.annotation.map`) is data-specific and
    must be re-derived for your own data.**
- **`archr_peak_calling.R`**: the sample loop, thresholds, LSI/UMAP/clustering
  params, and MACS2 peak-calling settings all move to `config.yaml:archr`. The
  called-cell barcode source is now `config.yaml:barcode_metrics` (a list of
  `{file, cell_col}` candidates tried in order), read via the new
  `valid_barcodes()` helper — supporting both Cell Ranger ARC
  (`per_barcode_metrics.csv` / `is_cell`) and older Cell Ranger
  (`singlecell.csv` / `is__cell_barcode`) without editing the script.
