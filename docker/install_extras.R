# ============================================================================
# install_extras.R — the R packages conda cannot provide at the needed version.
# ----------------------------------------------------------------------------
# Each of these was a real breakage during local setup; the comments record why
# the version is what it is, so nobody "helpfully" upgrades them later.
#
#   Rscript install_extras.R [with_archr]
# ============================================================================
args       <- commandArgs(trailingOnly = TRUE)
with_archr <- length(args) >= 1 && tolower(args[[1]]) %in% c("true", "1", "yes")

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1800)

need <- function(pkg) requireNamespace(pkg, quietly = TRUE)
say  <- function(...) cat("[install_extras]", ..., "\n")

# --- harmony 0.1.1 ----------------------------------------------------------
# ArchR::addHarmony() calls HarmonyMatrix(), which harmony >= 1.0 removed
# ("object 'HarmonyMatrix' not found"). conda-forge's oldest build is 1.2.0 --
# already too new -- so this comes from the CRAN archive. Only ArchR (Track B)
# needs it, but it is cheap and keeps the two tracks' envs identical.
if (!need("harmony") || as.character(packageVersion("harmony")) != "0.1.1") {
  say("installing harmony 0.1.1 from the CRAN archive")
  remotes::install_version("harmony", version = "0.1.1",
                           upgrade = "never", quiet = FALSE)
}
stopifnot(as.character(packageVersion("harmony")) == "0.1.1")
say("harmony", as.character(packageVersion("harmony")))

# --- DoubletFinder 2.0.6 ----------------------------------------------------
# GitHub-only. 2.0.6 changed the API to doubletFinder() / reuse.pANN = NULL,
# which pipeline/R/filtering.R targets. Asserted rather than tagged: the repo
# does not carry a reliable v2.0.6 git tag, so we install and then verify.
if (!need("DoubletFinder")) {
  say("installing DoubletFinder from GitHub")
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder",
                          upgrade = "never", quiet = FALSE)
}
df_v <- as.character(packageVersion("DoubletFinder"))
if (utils::compareVersion(df_v, "2.0.6") < 0) {
  stop("DoubletFinder ", df_v, " is older than 2.0.6; filtering.R expects the ",
       "2.0.6 API (doubletFinder(), reuse.pANN=NULL)")
}
say("DoubletFinder", df_v)

# --- ArchR 1.0.3 (optional; Track B only) -----------------------------------
# Not needed for Track A: the peak set is already baked into the merged object,
# and run_chrombpnet.sh re-calls peaks with MACS2 anyway. Build with
# --build-arg WITH_ARCHR=true when you need to re-call peaks for new data.
if (with_archr) {
  if (!need("ArchR")) {
    say("installing ArchR 1.0.3 from GitHub")
    remotes::install_github("GreenleafLab/ArchR", ref = "v1.0.3",
                            upgrade = "never", quiet = FALSE)
  }
  stopifnot(need("ArchR"))
  say("ArchR", as.character(packageVersion("ArchR")))
} else {
  say("skipping ArchR (WITH_ARCHR not set) -- Track A does not need it")
}

# --- verify the pipeline's own imports load --------------------------------
for (p in c("Seurat", "Signac", "dplyr", "yaml")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  say(p, as.character(packageVersion(p)))
}
say("OK")
