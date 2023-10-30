rm(list = ls())
library(tidyverse)
library(data.table)
library(iRF)

source(file.path("..", "functions", "load-functions.R"), chdir = TRUE)
source(file.path("..", "functions", "eval-functions.R"), chdir = TRUE)
source(file.path("..", "functions", "local_stability.R"), chdir = TRUE)
source(file.path("..", "functions", "utils.R"), chdir = TRUE)

RESULTS_DIR <- file.path("..", "results")

args <- commandArgs(trailingOnly = TRUE)
pheno_name <- as.character(args[1]) # name of phenotype
if (length(args) >= 2) { # number of snps to use as features
  nsnps <- as.numeric(args[2])
} else {
  nsnps <- 1000
}

print(paste0("Phenotype: ", pheno_name, "; # GWAS SNPs: ", nsnps))

save_path <- file.path(
  RESULTS_DIR, 
  paste0(pheno_name, "_gwas_filtered_nsnps", nsnps)
)
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

n_cores <- 5
test_ngenes <- 25
test_nints <- 25
n_iter <- 3

###############################################################################
## Run iRF
###############################################################################

if (!file.exists(file.path(save_path, "irf_interaction_fit.Rdata"))) {
  print("Runing iRF...")
  pheno_binary <- FALSE
  n <- 15000
  n_test <- 5000
  predx <- NULL
  include_dems <- TRUE
  include_npcs <- 0
  keep_genes <- NULL

  # get top snps from gwas
  bolt_dir <- file.path(RESULTS_DIR, "gwas_bolt_lmm")
  plink_dir <- file.path(RESULTS_DIR, "gwas_plink")
  keep_snps <- map(
    c(bolt_dir, plink_dir),
    function(fdir) {
      fpheno <- str_remove(pheno_name, "_binary_thr.*$")
      fpath <- file.path(fdir, paste0(fpheno, "_norm.stats_annot"))
      snps <- fread(fpath) %>%
        slice(1:nsnps) %>%
        pull(SNP)
    }
  ) %>%
    purrr::reduce(c) %>%
    unique()

  # load data
  out <- loadData(
    pheno_name = pheno_name, pheno_binary = pheno_binary,
    n = n, n_test = n_test, predx = predx,
    include_dems = include_dems, include_npcs = include_npcs,
    keep_genes = keep_genes, keep_snps = keep_snps
  )

  # omit additional covariates except for gender
  if (str_detect(pheno_name, "male")) {
    # no demographic covariates
    omit_vars <- c("gender", "height", "weight", "age", paste0("PC", 1:5))
  } else {
    # use gender only
    omit_vars <- c("height", "weight", "age", paste0("PC", 1:5))
  }
  out$geno_train <- out$geno_train[, !(colnames(out$geno_train) %in% omit_vars)]
  out$geno_test <- out$geno_test[, !(colnames(out$geno_test) %in% omit_vars)]

  # split training data into two training sets
  train1_idx <- 1:(round(nrow(out$geno_train) * 2 / 3))
  out$geno_train1 <- out$geno_train[train1_idx, ]
  out$pheno_train1 <- out$pheno_train[train1_idx]
  out$geno_train2 <- out$geno_train[-train1_idx, ]
  out$pheno_train2 <- out$pheno_train[-train1_idx]
  out$geno_train <- NULL
  out$pheno_train <- NULL
  str(out)

  # load in info about SNPs
  snps.df <- loadSNPInfo(colnames(out$geno_train1))

  # flip snp encoding
  flip_snp_names <- flipSNPData(
    geno_train = out$geno_train1,
    snps_df = snps.df
  )$flip_snps
  out$geno_train1[, flip_snp_names] <- 2 - out$geno_train1[, flip_snp_names]
  out$geno_train2[, flip_snp_names] <- 2 - out$geno_train2[, flip_snp_names]
  out$geno_test[, flip_snp_names] <- 2 - out$geno_test[, flip_snp_names]

  # run iRF at gene level
  x <- out$geno_train1
  y <- as.factor(out$pheno_train1)
  irf_gene_out <- iRF(
    x = x,
    y = y,
    varnames.grp = snps.df$Gene,
    n.iter = n_iter,
    iter.return = 1:n_iter,
    int.return = 1:n_iter,
    select.iter = FALSE,
    type = "ranger",
    respect.unordered.factors = "partition",
    rit.param = list(
      ntree = 500, depth = 3, nchild = 5,
      class.id = 1, class.cut = NULL,
      min.nd = 1
    ),
    ntree = 500,
    n.bootstrap = 50,
    n.core = n_cores
  )

  # make predictions
  yhat_tr2 <- predict(irf_gene_out$rf.list[[n_iter]],
    data = out$geno_train2,
    predict.all = T, num.threads = 1
  )$predictions %>%
    rowMeans()
  yhat_va <- predict(irf_gene_out$rf.list[[n_iter]],
    data = out$geno_test,
    predict.all = T, num.threads = 1
  )$predictions %>%
    rowMeans()

  # get interactions at SNV level for demonstration of instability
  ints.eval <- iRF::gRIT(irf_gene_out$rf.list[[n_iter]],
    x = x, y = y,
    rit.param = list(
      ntree = 500, depth = 3, nchild = 5,
      class.id = 1, class.cut = NULL,
      min.nd = 1
    ),
    n.core = n_cores
  )
  if (length(ints.eval) > 0) {
    bs.sample <- iRF:::lreplicate(50, iRF:::bsSample(y))
    irf_snv_int_out <- iRF::stabilityScore(
      x = x, y = y,
      ntree = 500,
      mtry.select.prob = irf_gene_out$rf.list[[2]][["variable.importance"]],
      ints.idx.eval = ints.eval$int.idx,
      rit.param = list(
        ntree = 500, depth = 3, nchild = 5,
        class.id = 1, class.cut = NULL,
        min.nd = 1
      ),
      bs.sample = bs.sample,
      type = "ranger",
      n.core = n_cores,
      respect.unordered.factors = "partition"
    )
  } else {
    irf_snv_int_out <- NULL
  }

  save(out, snps.df, irf_gene_out, yhat_tr2, yhat_va, irf_snv_int_out,
    file = file.path(save_path, "irf_interaction_fit.Rdata")
  )
  print("Completed.")
}

load(file.path(save_path, "irf_interaction_fit.Rdata"))

print(sprintf(
  "OOB Classification Error: %s",
  irf_gene_out$rf.list[[n_iter]]$prediction.error
))
print("OOB Confusion Matrix:")
print(irf_gene_out$rf.list[[n_iter]]$confusion.matrix)

print("Training2 Error")
bind_rows(
  evalPreds(
    y = out$pheno_train2, yhat = yhat_tr2,
    metric = c("AUC", "PR"),
    group = as.factor(out$geno_train2[, "gender"])
  ),
  evalPreds(
    y = out$pheno_train2, yhat = yhat_tr2 >= 0.5,
    metric = c("Class"),
    group = as.factor(out$geno_train2[, "gender"])
  )
) %>%
  print()

print("Validation Error")
bind_rows(
  evalPreds(
    y = out$pheno_test, yhat = yhat_va,
    metric = c("AUC", "PR"),
    group = as.factor(out$geno_test[, "gender"])
  ),
  evalPreds(
    y = out$pheno_test, yhat = yhat_va >= 0.5,
    metric = c("Class"),
    group = as.factor(out$geno_test[, "gender"])
  )
) %>%
  print()

###############################################################################
##  Evaluate local rf feature stability
###############################################################################

print("Evaluting local RF feature stability...")
load(file.path(save_path, "irf_interaction_fit.Rdata"))

feature_groups <- snps.df %>%
  select(feature = Name, group = Gene)

# compute local stability
stab_tr2_ls <- map(
  irf_gene_out$rf.list,
  ~ localFeatureStabilityRF(
    rf_fit = .x,
    X = out$geno_train2,
    feature_groups = feature_groups
  )
)

stab_va_ls <- map(
  irf_gene_out$rf.list,
  ~ localFeatureStabilityRF(
    rf_fit = .x,
    X = out$geno_test,
    feature_groups = feature_groups
  )
)

# evaluate via permutation test
test_genes_ls <- map(
  stab_tr2_ls,
  ~ data.frame(
    gene = colnames(.),
    mean_stability = apply(., 2, mean)
  ) %>%
    arrange(desc(mean_stability)) %>%
    dplyr::slice(1:test_ngenes) %>%
    pull(gene) %>%
    as.character() %>%
    setNames(., .)
)

perm_out_ls <- map2(
  test_genes_ls, stab_va_ls,
  function(test_genes, stab_va) {
    map(
      test_genes,
      ~ runPermutationTest(
        stab_df = stab_va,
        y = out$pheno_test,
        feature = .x, nperm = 1e4
      )
    )
  }
)

# save results
save(stab_tr2_ls, stab_va_ls, test_genes_ls, perm_out_ls,
  file = file.path(save_path, "local_feature_stability_results.Rdata")
)
print("Completed.")

###############################################################################
##  Evaluate local rf interaction stability
###############################################################################

print("Evaluating local RF interaction stability...")
load(file.path(save_path, "irf_interaction_fit.Rdata"))

feature_groups <- snps.df %>%
  select(feature = Name, group = Gene)

# select interactions to tests
irf_ints_ls <- map(irf_gene_out$interaction, ~ annotInts(.x, snps.df, F))
irf_keep_ints_ls <- map(
  irf_ints_ls,
  ~ .x %>%
    filter(
      stability >= 0.5,
      sta.fsd > 0,
      sta.mip > 0
    ) %>%
    arrange(desc(prevalence)) %>%
    slice(1:min(n(), test_nints))
)

# compute local interaction stability
int_stab_tr2_ls <- map2(
  irf_gene_out$rf.list, irf_keep_ints_ls,
  ~ localIntStabilityRF(
    rf_fit = .x,
    X = out$geno_train2,
    ints = .y %>% pull(int),
    feature_groups = feature_groups
  )
)

int_stab_va_ls <- map2(
  irf_gene_out$rf.list, irf_keep_ints_ls,
  ~ localIntStabilityRF(
    rf_fit = .x,
    X = out$geno_test,
    ints = .y %>% pull(int),
    feature_groups = feature_groups
  )
)

# evaluate via permutation test
test_ints_ls <- map(
  int_stab_tr2_ls,
  ~ data.frame(
    int = colnames(.),
    mean_stability = apply(., 2, mean)
  ) %>%
    arrange(desc(mean_stability)) %>%
    dplyr::slice(1:min(n(), test_nints)) %>%
    pull(int) %>%
    as.character() %>%
    setNames(., .)
)

int_perm_out_ls <- map2(
  test_ints_ls, int_stab_va_ls,
  function(test_ints, stab_va) {
    map(
      test_ints,
      ~ runPermutationTest(
        stab_df = stab_va,
        y = out$pheno_test,
        feature = .x, nperm = 1e4
      )
    )
  }
)

# save results
save(int_stab_tr2_ls, int_stab_va_ls, test_ints_ls, int_perm_out_ls, irf_ints_ls,
  file = file.path(save_path, "local_interaction_stability_results.Rdata")
)
print("Completed.")

###############################################################################
##  Evaluate local rf feature stability using first split only
###############################################################################

print("Evaluting local RF feature stability using the first split only...")
load(file.path(save_path, "irf_interaction_fit.Rdata"))

feature_groups <- snps.df %>%
  select(feature = Name, group = Gene)

# compute local stability
stab_tr2_ls <- map(
  irf_gene_out$rf.list,
  ~ localFeatureStabilityRF(
    rf_fit = .x,
    X = out$geno_train2,
    feature_groups = feature_groups,
    first_only = TRUE
  )
)

stab_va_ls <- map(
  irf_gene_out$rf.list,
  ~ localFeatureStabilityRF(
    rf_fit = .x,
    X = out$geno_test,
    feature_groups = feature_groups,
    first_only = TRUE
  )
)

# evaluate via permutation test
test_genes_ls <- map(
  stab_tr2_ls,
  ~ data.frame(
    gene = colnames(.),
    mean_stability = apply(., 2, mean)
  ) %>%
    arrange(desc(mean_stability)) %>%
    dplyr::slice(1:test_ngenes) %>%
    pull(gene) %>%
    as.character() %>%
    setNames(., .)
)

perm_out_ls <- map2(
  test_genes_ls, stab_va_ls,
  function(test_genes, stab_va) {
    map(
      test_genes,
      ~ runPermutationTest(
        stab_df = stab_va,
        y = out$pheno_test,
        feature = .x, nperm = 1e4
      )
    )
  }
)

# save results
save(stab_tr2_ls, stab_va_ls, test_genes_ls, perm_out_ls,
  file = file.path(
    save_path,
    "local_feature_stability_first_only_results.Rdata"
  )
)
print("Completed.")

###############################################################################
##  Evaluate local rf interaction stability using first split only
###############################################################################

print("Evaluating local RF interaction stability using first split only...")
load(file.path(save_path, "irf_interaction_fit.Rdata"))

feature_groups <- snps.df %>%
  select(feature = Name, group = Gene)

# select interactions to tests
irf_ints_ls <- map(irf_gene_out$interaction, ~ annotInts(.x, snps.df, F))
irf_keep_ints_ls <- map(
  irf_ints_ls,
  ~ .x %>%
    filter(
      stability >= 0.5,
      sta.fsd > 0,
      sta.mip > 0
    ) %>%
    arrange(desc(prevalence)) %>%
    slice(1:min(n(), test_nints))
)

# compute local interaction stability
int_stab_tr2_ls <- map2(
  irf_gene_out$rf.list, irf_keep_ints_ls,
  ~ localIntStabilityRF(
    rf_fit = .x,
    X = out$geno_train2,
    ints = .y %>% pull(int),
    feature_groups = feature_groups,
    first_only = TRUE
  )
)

int_stab_va_ls <- map2(
  irf_gene_out$rf.list, irf_keep_ints_ls,
  ~ localIntStabilityRF(
    rf_fit = .x,
    X = out$geno_test,
    ints = .y %>% pull(int),
    feature_groups = feature_groups,
    first_only = TRUE
  )
)

# evaluate via permutation test
test_ints_ls <- map(
  int_stab_tr2_ls,
  ~ data.frame(
    int = colnames(.),
    mean_stability = apply(., 2, mean)
  ) %>%
    arrange(desc(mean_stability)) %>%
    dplyr::slice(1:test_nints) %>%
    pull(int) %>%
    as.character() %>%
    setNames(., .)
)

int_perm_out_ls <- map2(
  test_ints_ls, int_stab_va_ls,
  function(test_ints, stab_va) {
    map(
      test_ints,
      ~ runPermutationTest(
        stab_df = stab_va,
        y = out$pheno_test,
        feature = .x, nperm = 1e4
      )
    )
  }
)

# save results
save(int_stab_tr2_ls, int_stab_va_ls, test_ints_ls, int_perm_out_ls,
  file = file.path(
    save_path,
    "local_interaction_stability_first_only_results.Rdata"
  )
)
print("Completed.")
