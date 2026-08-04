# ============================================================================
# Zamboni multiome — Track A (ChromBPNet input preparation), containerised.
# ----------------------------------------------------------------------------
# Build (from the repo root, i.e. TF_cobinding_multiome_chrombpnet/):
#   docker build -f docker/Dockerfile -t zamboni-tracka:1.0 .
#   docker build -f docker/Dockerfile --build-arg WITH_ARCHR=true \
#         -t zamboni-tracka:1.0-archr .        # also installs ArchR (Track B)
#
# amd64 ONLY, deliberately. bioconda's linux-aarch64 coverage is far thinner
# (e.g. bedtools has 3 aarch64 builds vs 72 for amd64), so exact pins do not
# solve on ARM. One architecture also means identical bits on a laptop and on
# an amd64 cluster, which is the point. On Apple Silicon this runs emulated:
# fine for these I/O-bound steps, and no GPU is involved at any stage here.
#
# Code and small references are baked in; DATA IS MOUNTED (see docker/README.md).
# ============================================================================
FROM --platform=linux/amd64 condaforge/miniforge3:24.9.2-0

LABEL org.opencontainers.image.title="Zamboni multiome Track A" \
      org.opencontainers.image.description="Per-cell-type ChromBPNet input prep: injury DARs, fragment BEDs, MACS2 peaks" \
      org.opencontainers.image.licenses="MIT"

# procps: Snakemake/R subprocess handling. build-essential + git: needed to
# compile the GitHub R packages in install_extras.R.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential git procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "-lc"]

# --- R environment ----------------------------------------------------------
COPY docker/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

# GenomeInfoDb 1.38.1 hard-codes a retired UCSC host; this shim repoints it.
# Must land in the env's R etc/ dir to apply to every R invocation.
COPY docker/Rprofile.site /opt/conda/envs/zamboni-r/lib/R/etc/Rprofile.site

# --- non-conda R packages (pinned; see the script for why each) -------------
ARG WITH_ARCHR=false
COPY docker/install_extras.R /tmp/install_extras.R
RUN mamba run -n zamboni-r Rscript /tmp/install_extras.R "${WITH_ARCHR}"

# --- command-line tools (separate env: macs2 needs its own Python) ----------
COPY docker/environment-tools.yml /tmp/environment-tools.yml
RUN mamba env create -f /tmp/environment-tools.yml && mamba clean -afy

# Put the tools on PATH so call_celltype_peaks.sh finds them via `command -v`
# instead of the host-specific absolute paths it used to hardcode.
ENV PATH=/opt/conda/envs/tools/bin:$PATH

# --- project code + small references ---------------------------------------
# Baked in on purpose: image = exact code + exact environment. Large inputs
# (data/, the merged .rds, mm10.fa) are excluded via .dockerignore and mounted.
WORKDIR /work
COPY pipeline/ /work/pipeline/
COPY code/     /work/code/
COPY meta/     /work/meta/

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Fail fast and loudly if the env is broken, rather than midway through a
# multi-hour step.
RUN mamba run -n zamboni-r Rscript -e \
      'suppressPackageStartupMessages({library(Seurat);library(Signac);library(dplyr);library(yaml)}); \
       cat("R stack OK\n")' \
 && macs2 --version && bedtools --version && bgzip --version | head -1

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["all"]
