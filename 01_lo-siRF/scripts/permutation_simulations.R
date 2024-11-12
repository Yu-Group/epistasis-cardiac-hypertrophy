rm(list = ls())
library(magrittr)
library(simChef)
library(future)

source(file.path("..", "functions", "local_stability.R"), chdir = TRUE)

RESULTS_DIR <- file.path("..", "results")
SAVE <- TRUE
USE_CACHED <- TRUE
N_REPS <- 200

options(simChef.plot_theme = "vthemes")

set.seed(331)
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
print(sprintf("Num. Cores = %s", n_cores))
plan(multicore, workers = n_cores)

#### Create experiment parts ####

gaussian_X_dgp_fun <- function(n, p, s = NULL, beta = NULL) {
  if (is.null(beta)) {
    beta <- rep(0, p)
  }
  if (!is.null(s) && (length(beta) == 1)) {
    beta <- c(rep(beta, s), rep(0, p - s))
  } else if (length(beta) == 1) {
    beta <- rep(beta, p)
  }
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  ylin <- X %*% beta
  prob <- 1 / (1 + exp(-ylin))
  y <- rbinom(n, size = 1, prob = prob)
  return(list(X = X, y = y))
}

generate_interaction <- function(X, int_order = 0, num_int = 0, beta = 0.4) {
  if ((int_order > 0) && (num_int > 0)) {
    prob <- purrr::map(
      0:(num_int - 1),
      ~ apply(
        X[, (.x * int_order + 1):(.x * int_order + int_order), drop = FALSE] > 0, 1, all
      ) * beta
    ) |>
      setNames(paste0("X", 1:num_int)) |>
      dplyr::bind_cols() |>
      rowSums()
  } else {
    prob <- beta
  }
  y <- rbinom(nrow(X), size = 1, prob = prob)
  return(y)
}

snp_X_dgp_fun <- function(n, p, int_order = 0, num_int = 0, beta = 0.4) {
  X <- matrix(sample(0:2, size = n * p, replace = TRUE), nrow = n, ncol = p)
  y <- generate_interaction(
    X,
    int_order = int_order,
    num_int = num_int,
    beta = beta
  )
  return(list(X = X, y = y))
}

snp_X_rwd_dgp_fun <- function(n, p, int_order = 0, num_int = 0, beta = 0.4) {
  
  # # uncomment and run once to generate data
  # pheno_name <- "iLVM_norm"
  # pheno_binary <- FALSE
  # n <- 15000
  # n_test <- 5000
  # predx <- NULL
  # include_dems <- TRUE
  # include_npcs <- 0
  # keep_genes <- NULL
  # 
  # # get top snps from gwas
  # bolt_dir <- file.path(RESULTS_DIR, "gwas_bolt_lmm")
  # plink_dir <- file.path(RESULTS_DIR, "gwas_plink")
  # keep_snps <- map(
  #   c(bolt_dir, plink_dir),
  #   function(fdir) {
  #     fpheno <- str_remove(pheno_name, "_binary_thr.*$")
  #     fpath <- file.path(fdir, paste0(fpheno, "_norm.stats_annot"))
  #     snps <- fread(fpath) |>
  #       slice(1:nsnps) |>
  #       pull(SNP)
  #   }
  # ) |>
  #   purrr::reduce(c) |>
  #   unique()
  # 
  # # load data
  # out <- loadData(
  #   pheno_name = pheno_name, pheno_binary = pheno_binary,
  #   n = n, n_test = n_test, predx = predx,
  #   include_dems = include_dems, include_npcs = include_npcs,
  #   keep_genes = keep_genes, keep_snps = keep_snps
  # )
  # saveRDS(out, file.path(DATA_DIR, "iLVM_norm_data.rds"))
  # 
  # # only keep genes with > 10 snps
  # snps_df <- loadSNPInfo(colnames(out$geno_train))
  # keep_genes <- snps_df |>
  #   dplyr::group_by(Gene) |>
  #   dplyr::summarise(
  #     n_snps = dplyr::n()
  #   ) |>
  #   dplyr::filter(
  #     n_snps >= 10
  #   ) |>
  #   dplyr::pull(Gene)
  # keep_snps_df <- snps_df |>
  #   dplyr::filter(Gene %in% keep_genes)
  # saveRDS(keep_snps_df, file.path(DATA_DIR, "snps_df_abridged.rds"))
  
  DATA_DIR <- file.path("..", "data")
  snps_df <- readRDS(file.path(DATA_DIR, "snps_df_abridged.rds"))
  data_out <- readRDS(file.path(DATA_DIR, "iLVM_norm_data.rds"))
  keep_samples <- sample(1:nrow(data_out$geno_train), size = n, replace = FALSE)
  keep_snps_df <- snps_df |>
    dplyr::slice_sample(n = p) |>
    dplyr::group_by(Gene) |>
    dplyr::mutate(
      order = 1:dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(order) |>
    dplyr::mutate(
      Gene = paste0("Gene", as.numeric(forcats::fct_inorder(Gene)))
    )
  X <- data_out$geno_train[keep_samples, keep_snps_df$Name]
  y <- generate_interaction(
    X,
    int_order = int_order,
    num_int = num_int,
    beta = beta
  )
  return(
    list(X = X, y = y, varnames_grp = keep_snps_df$Gene)
  )
}

losirf_method_fun <- function(X, y,
                              varnames_grp = NULL,
                              max_ints = 50,
                              cpe_thr = 0, 
                              fsd_thr = 0, 
                              mip_thr = 0,
                              stability_thr = 0.5,
                              stability_cpe_thr = 0.5,
                              stability_fsd_thr = 0.5,
                              stability_mip_thr = 0.5, ...) {
  
  result <- NULL
  
  # split data
  n <- nrow(X)
  train_idx <- sample(1:n, size = round(n * 2/3), replace = FALSE)
  X_train <- data.frame(X[train_idx, , drop = FALSE])
  y_train <- y[train_idx]
  X_test <- data.frame(X[-train_idx, , drop = FALSE])
  y_test <- y[-train_idx]
  
  irf_gene_out <- iRF::iRF(
    x = as.matrix(X_train),
    y = as.factor(y_train),
    varnames.grp = varnames_grp,
    n.iter = 3,
    iter.return = 1:3,
    int.return = 3,
    select.iter = FALSE,
    type = 'ranger',
    n.bootstrap = 50,
    n.core = 1,
    ...
  )
  
  keep_ints <- irf_gene_out$interaction |>
    dplyr::filter(
      cpe > !!cpe_thr,
      fsd > !!fsd_thr,
      mip > !!mip_thr,
      stability > !!stability_thr,
      sta.cpe > !!stability_cpe_thr,
      sta.fsd > !!stability_fsd_thr,
      sta.mip > !!stability_mip_thr
    ) |>
    dplyr::arrange(desc(prevalence))
  
  if (nrow(keep_ints) > 0) {
    if (!is.null(max_ints)) {
      keep_ints <- keep_ints |>
        dplyr::slice(1:min(dplyr::n(), max_ints))
    }
    
    if (!is.null(varnames_grp)) {
      feature_groups <- data.frame(
        feature = colnames(X_train), 
        group = varnames_grp
      )
    } else {
      feature_groups <- NULL
    }
    
    # compute local interaction stability
    lsi_test <- localIntStabilityRF(
      rf_fit = irf_gene_out$rf.list[[3]],
      X = X_test,
      ints = keep_ints$int,
      feature_groups = feature_groups
    )
    
    pvals <- purrr::map_dbl(
      colnames(lsi_test),
      ~ runPermutationTest(
        stab_df = lsi_test,
        y = y_test,
        feature = .x, 
        nperm = 1e4
      )$pval
    )
    
    result <- tibble::tibble(int = colnames(lsi_test), pval = pvals)
  }
  return(
    list(
      result = result,
      rf_fit = irf_gene_out$rf.list[[3]],
      int_df = irf_gene_out$interaction,
      X_test = X_test,
      y_test = y_test,
      varnames_grp = varnames_grp,
      n_int_total = nrow(irf_gene_out$interaction),
      n_int_tested = nrow(keep_ints)
    )
  )
}

#### Experiment Parts ####
gaussian_X_dgp <- create_dgp(
  .dgp_fun = gaussian_X_dgp_fun,
  .name = "Gaussian X",
  n = 1000, p = 100, s = 2, beta = 0
)

snp_X_dgp <- create_dgp(
  .dgp_fun = snp_X_dgp_fun,
  .name = "SNP X",
  n = 1000, p = 100, int_order = 0, num_int = 0, beta = 0.4
)
snp_X_dgp1 <- create_dgp(
  .dgp_fun = snp_X_dgp_fun,
  .name = "SNP X (order 1)",
  n = 1000, p = 100, int_order = 1, num_int = 2, beta = 0.4
)
snp_X_dgp2 <- create_dgp(
  .dgp_fun = snp_X_dgp_fun,
  .name = "SNP X (order 2)",
  n = 1000, p = 100, int_order = 2, num_int = 2, beta = 0.4
)

snp_X_rwd_dgp <- create_dgp(
  .dgp_fun = snp_X_rwd_dgp_fun,
  .name = "Real SNP X",
  n = 2500, p = 500, int_order = 0, num_int = 0, beta = 0.4
)

losirf_filter_method <- create_method(
  .method_fun = losirf_method_fun,
  .name = "lo-siRF"
)

losirf_method <- create_method(
  .method_fun = losirf_method_fun,
  .name = "lo-siRF (no filtering)",
  stability_cpe_thr = 0,
  stability_fsd_thr = 0,
  stability_mip_thr = 0
)

#### Run Experiments ####

exp_name <- "Permutation Validity"

## Null - Gaussian X DGP
experiment <- create_experiment(
  name = exp_name, save_dir = file.path(RESULTS_DIR, exp_name)
) |>
  add_dgp(gaussian_X_dgp) |>
  add_method(losirf_method) |>
  add_vary_across(
    .dgp = gaussian_X_dgp$name,
    n = c(1000, 1500, 2000)
  )
out <- run_experiment(
  experiment, n_reps = N_REPS, use_cached = USE_CACHED, save = SAVE
)

## Null - SNP X DGP
experiment <- create_experiment(
  name = exp_name, save_dir = file.path(RESULTS_DIR, exp_name)
) |>
  add_dgp(snp_X_dgp) |>
  add_method(losirf_method) |>
  add_vary_across(
    .dgp = snp_X_dgp$name,
    n = c(1000, 1500, 2000)
  )
out <- run_experiment(
  experiment, n_reps = N_REPS, use_cached = USE_CACHED, save = SAVE
)

## Null - Real SNP X DGP
experiment <- create_experiment(
  name = exp_name, save_dir = file.path(RESULTS_DIR, exp_name)
) |>
  add_dgp(snp_X_rwd_dgp) |>
  add_method(losirf_method) |>
  add_vary_across(
    .dgp = snp_X_rwd_dgp$name,
    n = c(1000, 2000, 3000)
  )
out <- run_experiment(
  experiment, n_reps = N_REPS, use_cached = USE_CACHED, save = SAVE
)

## Marginal - SNP X DGP
experiment <- create_experiment(
  name = exp_name, save_dir = file.path(RESULTS_DIR, exp_name)
) |>
  add_dgp(snp_X_dgp1) |>
  add_method(losirf_filter_method) |>
  add_vary_across(
    .dgp = snp_X_dgp1$name,
    n = c(1000, 2000, 3000)
  )
out <- run_experiment(
  experiment, n_reps = N_REPS, use_cached = USE_CACHED, save = SAVE
)

## Interaction - SNP X DGP
experiment <- create_experiment(
  name = exp_name, save_dir = file.path(RESULTS_DIR, exp_name)
) |>
  add_dgp(snp_X_dgp2) |>
  add_method(losirf_filter_method) |>
  add_vary_across(
    .dgp = snp_X_dgp2$name,
    n = c(1000, 2000, 3000)
  )
out <- run_experiment(
  experiment, n_reps = N_REPS, use_cached = USE_CACHED, save = SAVE
)
