# Containerised Track A — reproducible ChromBPNet input preparation

One image that turns **the annotated merged object + the 21 fragment files** into
everything ChromBPNet needs, on any machine. No GPU is involved at any stage
here; this covers the whole pre-training half of the pipeline.

## What it produces

| output | contents |
|---|---|
| `markers/markers_Uvsothers_peaks_all.rds` | per-cell-type injury DARs (40 comparisons) |
| `markers/markers_Uvsothers_rna_all.rds` | the same contrasts on RNA |
| `ML/celltype_fragments/` | fragment BEDs per cell type + bias-model input |
| `ML/celltype_fragments_injury/` | per cell type × injury |
| `ML/celltype_dars/` | DAR BEDs |
| `ML/macs2_celltype_peaks/`, `ML/filtered_peaks/` | MACS2 peaks, blacklist-filtered |
| `ML/dars_intersect_peaks/` | DARs ∩ called peaks — ChromBPNet's contribution-scoring regions |

## Build

From the repo root (`TF_cobinding_multiome_chrombpnet/`), **not** from `docker/`:

```bash
docker build -f docker/Dockerfile -t zamboni-tracka:1.0 .
```

Add ArchR (only needed to re-call peaks for *new* data — Track B):

```bash
docker build -f docker/Dockerfile --build-arg WITH_ARCHR=true -t zamboni-tracka:1.0-archr .
```

**amd64 only, on purpose.** bioconda's linux-aarch64 coverage is too thin to pin
these versions (`bedtools`: 3 ARM builds vs 72 for amd64), and a single
architecture means identical bits on a laptop and on an amd64 cluster. On Apple
Silicon it runs emulated — fine, since every step here is I/O-bound.

## Run

Data is **mounted, never baked in** (the inputs are ~30 GB). The container
expects the repo layout under `/work`:

```bash
docker run --rm \
  -v "$PWD/data:/work/data:ro" \
  -v "$PWD/multiome_obj_clean_240326.rds:/work/multiome_obj_clean_240326.rds:ro" \
  -v "$PWD/ML:/work/ML" \
  -v "$PWD/markers:/work/markers" \
  -v "$PWD/logs:/work/logs" \
  zamboni-tracka:1.0 all
```

Or use the wrapper, which wires those mounts for you:

```bash
./run_tracka.sh check     # validate inputs, no compute
./run_tracka.sh all       # everything
```

### Commands

| command | does |
|---|---|
| `check` | verify object, fragments, `.tbi` indexes, references, tools — then exit |
| `dge` | step 2 — injury DARs (peaks + RNA) |
| `split` | step 3 — fragment BEDs (~30 min, reads all 21 GB) |
| `dars` | step 3b — DAR BEDs (seconds) |
| `peaks` | step 1 — MACS2 peaks + DAR intersect (~25 min) |
| `all` | all of the above, in dependency order |
| `shell` | interactive bash |

**Always run `check` first.** It catches the failure modes that otherwise
surface an hour in: an unmounted object, missing `.tbi` indexes, absent
references.

Order is by *dependency*, not the numbering in `run_chrombpnet.sh`: `dars`
needs `dge`, and `peaks` needs both `split` and `dars`.

## Resource requirements

- **RAM: 16 GB minimum, 24 GB recommended.** The merged object is 18.8 GB
  loaded; `DietSeurat()` trims it to 11.5 GB (peaks) / 3.8 GB (RNA), which is
  what makes the largest cell-type subset fit. Below ~16 GB expect
  `vector memory exhausted` on Oligodendrocytes. Raise Docker Desktop's memory
  limit — its default (often 8 GB) is **not** enough.
- **Disk: ~100 GB** for `ML/` (fragment BEDs are ~43 GB, plus `_chr.bed` copies).
- **Time: ~1 hour** for `all` natively; longer emulated on ARM.

## Idempotency

Every step resumes except one. `split` **refuses to write into a non-empty
output directory**, because `SplitFragments` appends across the 21 fragment
files and a second run would silently double every BED. Clear
`ML/celltype_fragments*/` to redo it.

`peaks` skips cell types whose `.narrowPeak` already exists, so an interrupted
run does not redo MACS2.

## What is pinned, and why

The **image is the lockfile**. `environment.yml` only needs to solve once; these
are the versions that mattered:

| pin | reason |
|---|---|
| `r-base=4.3` | fixes Bioconductor 3.18 → GenomeInfoDb 1.38.x |
| `r-seurat=5.3.0` | the object is a Seurat v5 object |
| `r-signac=1.14.0` | `SplitFragments`, `ChromatinAssay` |
| `harmony==0.1.1` | `ArchR::addHarmony()` calls `HarmonyMatrix()`, removed in ≥1.0; conda's oldest (1.2.0) is already too new |
| `DoubletFinder ≥2.0.6` | `filtering.R` targets the 2.0.6 API (`doubletFinder()`, `reuse.pANN=NULL`) |
| `macs2=2.2.9.1`, `bedtools=2.31.1` | what the results were produced with |

`Rprofile.site` carries a required shim: GenomeInfoDb 1.38.1 hard-codes the
retired UCSC host `hgdownload.cse.ucsc.edu`, whose certificate no longer
matches. Setting the option up front does not work — `.onLoad` overwrites it —
so the fix is a load hook repointing to `hgdownload.soe.ucsc.edu`. Without it
any call reaching UCSC fails.

`macs2` lives in a **separate conda env** from R: it needs its own Python, which
is what forced the same split on the host.

## Known data caveats

Carried by the outputs, not the container — worth knowing before interpreting:

- **Macrophages: 273 DAR regions.** Its uninjured group has only 22 cells. The
  regions are trustworthy (99% land in called peaks) but there are far too few
  for MoDISco seqlet clustering.
- **Endothelial: only 28% of DARs survive the peak intersect**, versus 89–99%
  elsewhere. It has the second-most DARs from the fewest cells (439), including
  9-cell and 30-cell contrasts. Those underpowered comparisons produce DARs that
  do not sit in real peaks. Treat its motifs cautiously.
- Duplicate DARs across timepoints are collapsed (`ml.dars_dedupe`). 48% of
  Ependymal and Microglia rows were duplicates; leaving them would make MoDISco
  double-count seqlets in persistently-open regions. Set `dars_dedupe: false`
  for the published behaviour.

## Not included

ChromBPNet itself. Training needs a GPU and the authors' own image:

```bash
apptainer build chrombpnet.sif docker://kundajelab/chrombpnet:latest
apptainer exec --nv chrombpnet.sif chrombpnet pipeline ...
```

`mm10.fa` is also excluded (2.8 GB) — fetch and verify with:

```bash
curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
# md5 db005b65828db31735f384e4c5787be5
gunzip -c mm10.fa.gz > meta/mm10.fa && samtools faidx meta/mm10.fa
```
