library(tidyverse)
library(data.table)

source(file.path("..", "functions", "run-analysis.R"), chdir = T)

RESULTS_DIR <- file.path("..", "results")

args <- commandArgs(trailingOnly = TRUE)
pheno_name <- as.character(args[1]) # name of phenotype
if (length(args) >= 2) { # number of snps to use as features
  nsnps <- as.numeric(args[2])
} else {
  nsnps <- 1000
}
if (length(args) >= 3) { # whether or not to include demographic information
  include_dems <- as.numeric(args[3]) == 1
} else {
  include_dems <- TRUE
}

# fit models using top 1000 gwas snps
pheno_binary <- FALSE
n <- 15000
n_test <- 5000
predx <- NULL
include_npcs <- 5
keep_genes <- NULL
n_trees <- 500
n_trials <- 1
if (include_dems) {
  models <- c(
    "lasso_baseline", "ridge_baseline", "rf_baseline",
    "lasso", "lasso_std", "ridge", "ridge_std", "rf", "irf"
  )
  save_path <- file.path(
    RESULTS_DIR,
    paste0(
      pheno_name, "_gwas_filtered_nsnps", nsnps,
      "_dems"
    )
  )
} else {
  models <- c("lasso", "lasso_std", "ridge", "ridge_std", "rf", "irf")
  save_path <- file.path(
    RESULTS_DIR,
    paste0(pheno_name, "_gwas_filtered_nsnps", nsnps)
  )
}

# get top snps from gwas
bolt_dir <- file.path(RESULTS_DIR, "gwas_bolt_lmm")
plink_dir <- file.path(RESULTS_DIR, "gwas_plink")
keep_snps <- map(
  c(bolt_dir, plink_dir),
  function(fdir) {
    fpheno <- str_remove(pheno_name, "_binary_thr.*$")
    fpath <- file.path(
      fdir,
      paste0(fpheno, "_norm.stats_annot")
    )
    snps <- fread(fpath) %>%
      dplyr::slice(1:nsnps) %>%
      pull(SNP)
  }
) %>%
  purrr::reduce(c) %>%
  unique()

runAnalysis(
  pheno_name = pheno_name, pheno_binary = pheno_binary,
  n = n, n_test = n_test, models = models, predx = predx,
  include_dems = include_dems, include_npcs = include_npcs,
  keep_genes = keep_genes, keep_snps = keep_snps,
  n_trees = n_trees, n_trials = n_trials, save_path = save_path,
  n_cores = 1
)
