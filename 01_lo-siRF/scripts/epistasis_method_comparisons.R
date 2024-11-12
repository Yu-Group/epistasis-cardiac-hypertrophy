rm(list = ls())

source(file.path("..", "functions", "load-functions.R"), chdir = TRUE)

DATA_DIR <- file.path("..", "data")
RESULTS_DIR <- file.path("..", "results", "epistasis_comparisons")

maf_df <- data.table::fread(file.path(DATA_DIR, "maf.frq")) |>
  tibble::as_tibble()

#### Regression-based SNP x SNP epistasis scan ####
pheno_names <- c(
  "iLVM",
  "iLVM_norm",
  "iLVM_binary_thr0.15",
  "iLVM_binary_thr0.2",
  "iLVM_binary_thr0.25"
)
nsnps <- 1000

for (pheno_name in pheno_names) {
  print(pheno_name)
  if (stringr::str_detect(pheno_name, "binary")) {
    outcome_type <- "binary"
  } else {
    outcome_type <- "continuous"
  }
  
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
      snps <- fread(fpath) |> 
        slice(1:nsnps) |> 
        pull(SNP)
    }
  ) |> 
    purrr::reduce(c) |> 
    unique()
  
  # load data
  out <- loadData(
    pheno_name = pheno_name, pheno_binary = pheno_binary,
    n = n, n_test = n_test, predx = predx,
    include_dems = include_dems, include_npcs = include_npcs,
    keep_genes = keep_genes, keep_snps = keep_snps
  )
  
  X <- out$geno_train |> 
    dplyr::select(weight, height, gender, age)
  Z <- out$geno_train |>
    dplyr::select(-tidyselect::all_of(colnames(X)))
  y <- out$pheno_train

  # remove snps with MAF < 0.05
  snps_df <- loadSNPInfo(colnames(Z))
  rm_snps <- maf_df |>
    dplyr::filter(MAF < 0.05) |>
    dplyr::left_join(snps_df, by = c("SNP" = "rsID")) |>
    dplyr::pull(Name)
  Z <- Z |>
    dplyr::select(-tidyselect::all_of(rm_snps))
  
  snp_pairs_mat <- combn(colnames(Z), 2)
  pvals_df <- purrr::map(
    1:ncol(snp_pairs_mat),
    function(j) {
      fit_df <- dplyr::bind_cols(
        .y = y,
        X,
        snp1 = Z[[snp_pairs_mat[1, j]]],
        snp2 = Z[[snp_pairs_mat[2, j]]]
      )
      if (outcome_type == "binary") {
        lm_fit <- glm(
          .y ~ weight + height + gender + age + snp1 + snp2 + snp1 * snp2, 
          data = fit_df, family = "binomial"
        )
      } else if (outcome_type == "continuous") {
        lm_fit <- lm(
          .y ~ weight + height + gender + age + snp1 + snp2 + snp1 * snp2, 
          data = fit_df
        )
      }
      lm_out <- broom::tidy(lm_fit) |>
        dplyr::select(term, p.value) |>
        tidyr::pivot_wider(names_from = term, values_from = p.value) |>
        dplyr::select(snp1, snp2, `snp1:snp2`) |>
        dplyr::mutate(
          snp1_name = snp_pairs_mat[1, j],
          snp2_name = snp_pairs_mat[2, j]
        )
    }
  ) |>
    purrr::list_rbind() |>
    dplyr::relocate(snp1_name, snp2_name, .before = snp1) |>
    dplyr::left_join(
      snps_df |> 
        dplyr::select(snp1_name = Name, `SNP1 Chr` = Chr, `SNP1 Gene` = Gene), 
      by = "snp1_name"
    ) |>
    dplyr::left_join(
      snps_df |> 
        dplyr::select(snp2_name = Name, `SNP2 Chr` = Chr, `SNP2 Gene` = Gene), 
      by = "snp2_name"
    )
  
  saveRDS(
    pvals_df, 
    file = file.path(
      RESULTS_DIR, "pairwise_interaction_scan", sprintf("%s.rds", pheno_name)
    )
  )
}

#### MapIt ####
# install.packages("mvMAPIT")

pheno_name <- "iLVM_norm"
pheno_binary <- FALSE
n <- 15000
n_test <- 5000
predx <- NULL
include_dems <- TRUE
include_npcs <- 0
keep_genes <- NULL
nsnps <- 1000

# get top snps from gwas
bolt_dir <- file.path(RESULTS_DIR, "gwas_bolt_lmm")
plink_dir <- file.path(RESULTS_DIR, "gwas_plink")
keep_snps <- map(
  c(bolt_dir, plink_dir),
  function(fdir) {
    fpheno <- str_remove(pheno_name, "_binary_thr.*$")
    fpath <- file.path(fdir, paste0(fpheno, "_norm.stats_annot"))
    snps <- fread(fpath) |> 
      slice(1:nsnps) |> 
      pull(SNP)
  }
) |> 
  purrr::reduce(c) |> 
  unique()

# load data
out <- loadData(
  pheno_name = pheno_name, pheno_binary = pheno_binary,
  n = n, n_test = n_test, predx = predx,
  include_dems = include_dems, include_npcs = include_npcs,
  keep_genes = keep_genes, keep_snps = keep_snps
)

X <- out$geno_train |> 
  dplyr::select(weight, height, gender, age)
Z <- out$geno_train |>
  dplyr::select(-tidyselect::all_of(colnames(X)))
y <- out$pheno_train

# remove snps with MAF < 0.05
snps_df <- loadSNPInfo(colnames(Z))
rm_snps <- maf_df |>
  dplyr::filter(MAF < 0.05) |>
  dplyr::left_join(snps_df, by = c("SNP" = "rsID")) |>
  dplyr::pull(Name)
Z <- Z |>
  dplyr::select(-tidyselect::all_of(rm_snps))

# run mapit
start_time <- Sys.time()
mvmapit_out <- mvMAPIT::mvmapit(
  X = t(as.matrix(Z)),
  Y = t(as.matrix(y)),
  Z = t(as.matrix(X)),
  test = "normal"
)
end_time <- Sys.time()
print(end_time - start_time)

saveRDS(mvmapit_out, file.path(RESULTS_DIR, "mvmapit_normal_nsnps1000.rds"))

#### MAPIT + GSEA ####
pheno_name <- "iLVM_norm"
pheno_binary <- FALSE
n <- 15000
n_test <- 5000
predx <- NULL
include_dems <- TRUE
include_npcs <- 0
keep_genes <- NULL
nsnps <- 10000

# get top snps from gwas
bolt_dir <- file.path(RESULTS_DIR, "gwas_bolt_lmm")
plink_dir <- file.path(RESULTS_DIR, "gwas_plink")
keep_snps <- map(
  c(bolt_dir, plink_dir),
  function(fdir) {
    fpheno <- str_remove(pheno_name, "_binary_thr.*$")
    fpath <- file.path(fdir, paste0(fpheno, "_norm.stats_annot"))
    snps <- fread(fpath) |> 
      slice(1:nsnps) |> 
      pull(SNP)
  }
) |> 
  purrr::reduce(c) |> 
  unique()

# load data
out <- loadData(
  pheno_name = pheno_name, pheno_binary = pheno_binary,
  n = n, n_test = n_test, predx = predx,
  include_dems = include_dems, include_npcs = include_npcs,
  keep_genes = keep_genes, keep_snps = keep_snps
)

X <- out$geno_train |> 
  dplyr::select(weight, height, gender, age)
Z <- out$geno_train |>
  dplyr::select(-tidyselect::all_of(colnames(X)))
y <- out$pheno_train

# run mapit
start_time <- Sys.time()
mvmapit_out <- mvMAPIT::mvmapit(
  X = t(as.matrix(Z)),
  Y = t(as.matrix(y)),
  Z = t(as.matrix(X)),
  test = "normal"
)
end_time <- Sys.time()
print(end_time - start_time)

# # run mapit in parallel
# library(Rmpi)
# start_time <- Sys.time()
# mvmapit_out <- mvMAPIT::mvmapit(
#   X = t(as.matrix(Z)),
#   Y = t(as.matrix(y)),
#   Z = t(as.matrix(X)),
#   cores = 10,
#   test = "normal"
# )
# end_time <- Sys.time()
# print(end_time - start_time)

saveRDS(mvmapit_out, file.path(RESULTS_DIR, "mvmapit_normal_nsnps10000.rds"))

# Run GSEA on mvMAPIT results