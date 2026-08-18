# TF_cobinding_multiome_chrombpnet

Regulatory genomics pipeline for discovering co-binding transcription factors in
cell-type-specific enhancers, built on **ChromBPNet + TF-MoDISco**.

Starting point is the single-nucleus multiome (paired snRNA-seq + snATAC-seq)
dataset of mouse spinal cord across injury timepoints (uninjured, 1/3/7/28 dpi)
from *Decoding injury-responsive enhancers in the CNS for cell state targeting*
(Zamboni et al.). 21 samples, 67,072 cells, 10 cell types, mm10.
GEO: `GSE304196`, `GSE304349`, `GSE304399`.

---

## Provenance — read this first

This repository contains two different things, and the distinction matters:

| directory | what it is |
|---|---|
| **`code/`** | **The authors' published scripts, unmodified.** Not our work. Kept verbatim as the reference implementation, under the original `LICENSE`. |
| **`pipeline/`** | **Our config-driven re-implementation.** Same analysis, but parameterised, logged, runnable off the authors' HPC, and with a number of defects in the published code corrected. |

The published scripts have paths hard-coded to the authors' cluster
(`/cfs/klemming/...`, `/proj/...`) and are written to be run interactively,
section by section — not as `Rscript script.R`. That is the main reason
`pipeline/` exists.

**Every departure from the published behaviour is documented in the header of
the corresponding `pipeline/` script**, with the reason and the published line
numbers. Where a change alters results rather than just mechanics, it is
switchable via `pipeline/config.yaml` (e.g. `ml.dars_dedupe: false` restores the
original behaviour).

---

## What it produces

Per cell type, from the annotated multiome object plus raw ATAC fragments:

```
fragments ──> MACS2 peaks ──> blacklist filter ──┐
                                                 ├──> ChromBPNet model ──> contribution
injury DARs ─────────────────> DAR ∩ peak ───────┘                          scores ──> MoDISco motifs
```

Current state: **all non-GPU inputs are built and verified.** 190,742
contribution-scoring regions across 10 cell types — 10-column narrowPeak,
uniform 500 bp, zero duplicates, every window verified edge-safe for ChromBPNet's
2114 bp input length. Model training is the remaining step and needs a GPU.

---

## Layout

```
pipeline/                     our corrected, config-driven implementation
  config.yaml                 single source of truth for paths and thresholds
  R/differential_analysis_injury.R   injury DARs (peaks + RNA)
  R/prepare_ml_input.R               per-cell-type fragment BEDs (SplitFragments)
  R/prepare_ml_dars.R                DAR BEDs
  call_celltype_peaks.sh             MACS2 peaks, blacklist, DAR ∩ peak intersect
  run_chrombpnet.sh                  bias model -> models -> contribs -> MoDISco
  README.md                          detailed per-step notes
code/                         the authors' published scripts (reference, unmodified)
docker/                       reproducibility image (zamboni-tracka:1.0, amd64)
meta/                         chrom sizes, blacklist, motif DB, chromosome folds
```

Not in the repo (see `.gitignore`): the merged Seurat object (7.2 GB), raw
fragments, `mm10.fa`, and all `ML/` outputs (~133 GB). All are regenerable —
each ignore rule names the command that rebuilds it.

---

## Running it

Track A (CPU, local) builds every ChromBPNet input:

```bash
Rscript pipeline/R/differential_analysis_injury.R pipeline/config.yaml --only peaks
Rscript pipeline/R/prepare_ml_input.R  pipeline/config.yaml
Rscript pipeline/R/prepare_ml_dars.R   pipeline/config.yaml
bash    pipeline/call_celltype_peaks.sh
```

Then ChromBPNet, split so only training needs a GPU:

```bash
bash pipeline/run_chrombpnet.sh prep        # bias peaks, splits, negatives (CPU)
bash pipeline/run_chrombpnet.sh bias-train  # GPU
bash pipeline/run_chrombpnet.sh train       # GPU
bash pipeline/run_chrombpnet.sh interpret   # GPU
DRY_RUN=1 bash pipeline/run_chrombpnet.sh all   # print every command, run nothing
```

Or reproduce Track A in the container (`docker/README.md`):

```bash
docker build -f docker/Dockerfile -t zamboni-tracka:1.0 .
./run_tracka.sh all
```

> ⚠️ **Use `pipeline/run_chrombpnet.sh`, not `code/machine_learning/run_chrombpnet.sh`.**
> The published script is not runnable as written; the corrected version fixes 10
> defects, listed in its header. Two are *silent* — they produce wrong results
> rather than an error:
> - `{i}` instead of `${i}` in every interpret/MoDISco output path, so all 10
>   cell types wrote into one directory literally named `{i}`, each overwriting
>   the last.
> - Contribution scores computed from `chrombpnet.h5` instead of
>   `chrombpnet_nobias.h5`. The former still carries the Tn5 bias component, so
>   MoDISco recovers the insertion motif rather than TF motifs.

Requirements: R 4.x with Seurat 5.3 / Signac 1.14, MACS2, bedtools, and the
`chrombpnet` conda environment. Docker Desktop needs ≥16 GB for the container —
the slimmed object alone is 11.5 GB. See `docker/README.md`.

---

## Data caveats

These are properties of the dataset, not bugs, and they become invisible once you
are looking only at motif output. **Two independent checks matter, and passing
one says nothing about the other.**

*Peak retention* — do the DARs (called on the pooled consensus peak set) land
inside peaks called independently per cell type from raw fragments? *Multiple
testing* — do they survive correction across the ~466,834 peaks tested?

| cell type | regions | peak retention | survives Bonferroni |
|---|---|---|---|
| Ependymal | 83,906 | 96% | 49% |
| Microglia | 38,309 | 89% | 82% |
| Astrocytes | 18,027 | 93% | 58% |
| Oligodendrocytes | 6,716 | 95% | 70% |
| OPCs | 22,671 | 89% | **0.8%** |
| Neurons_V | 5,646 | 92% | 4.4% |
| Neurons_D | 3,000 | 95% | 1.4% |
| Perivascular | 3,648 | 81% | 0% |
| Endothelial | 8,546 | **28%** | 1.8% |
| Macrophages | 273 | 99% | 0% |

- **Ependymal, Microglia, Astrocytes and Oligodendrocytes pass both** and are the
  defensible cell types for motif claims.
- **OPCs is the cautionary case**: it looks healthy on peak retention (89%) and
  loses 99.2% of its regions under correction.
- **Endothelial** has the second-most DARs from the fewest cells (439), including
  9-cell and 30-cell contrasts. Perivascular is the control — nearly as few cells
  (673) but 81% retention — so the cause is thin *contrasts*, not cell count.
- **Macrophages** is excluded from MoDISco automatically (273 regions is far too
  few for seqlet clustering). Biologically expected: macrophages infiltrate after
  injury, so its uninjured group is only 22 cells.

**Open question — the DAR significance threshold.** The region set fed to
ChromBPNet is filtered on raw `p_val` with no multiple-testing correction
(`code/machine_learning/prepare_ml_input.R:53`, on a table already filtered the
same way at `code/explore/differential_analysis.R:134,137`). ~466,834 peaks are
tested per contrast. This is faithful to the paper and is the default here.

Note that Seurat computes a Bonferroni `p_val_adj` inside every `FindMarkers`
call, so the corrected values are present in all of the authors' tables — the
question is only what each downstream filter *uses*, and the published repository
is not consistent about it. Adjusted p is used for the cell-type marker analysis
(`differential_analysis.R:46,49`) and for the shared-glial-gene figure of the
injury analysis (`explore_injury_differential_analysis.R:192`); raw p is used for
the DAR set that feeds the deep-learning model. The `avg_log2FC > 0.5` gate that
accompanies it is an effect-size filter, not error control, and it is weakest
exactly where it is needed most — in a thin contrast a large fold-change is *more*
likely by chance, because the estimate is unstable.

A BH-FDR middle ground cannot be recovered from the saved files: the DGE step
applies its filter before `saveRDS`, and BH is rank-based over the *complete* set
of p-values, so `p.adjust` on a pre-filtered table passes essentially everything.
Recovering it needs a ~1 h re-run of the DGE with filtering disabled. Bonferroni,
by contrast, stays valid after subsetting, which is why it could be used to
measure the problem in the table above.

The chromosome folds cover chr1–19 only, so ~1–2% of regions (on chrX) are scored
by a model that never trained on that chromosome. Kept by default and logged per
cell type; `DROP_UNSPLIT_CHROMS=1` excludes them.

---

## License

`code/` is the authors' published work and is covered by the original `LICENSE`
carried over from their repository. Attribution for `pipeline/` and the corrections
documented in it is not yet settled — treat this repository as internal for now.
