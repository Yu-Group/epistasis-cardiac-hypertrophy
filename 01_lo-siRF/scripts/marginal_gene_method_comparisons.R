rm(list = ls())
library(SKAT)

source(file.path("..", "functions", "load-functions.R"), chdir = TRUE)

DATA_DIR <- file.path("..", "data")
RESULTS_DIR <- file.path("..", "results", "marginal_comparisons")

#### SKAT-O ####
maf_df <- data.table::fread(file.path(DATA_DIR, "maf.frq")) |>
  tibble::as_tibble()

# continuous outcome
out_type <- "C"

# load data
pheno_name <- "iLVM_norm"
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

# get snp sets by gene
snps_df <- loadSNPInfo(colnames(Z))
snp_sets <- snps_df |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(snp_set = list(Name))

# run skat
null_train_df <- dplyr::bind_cols(.y = y, X)
train_df <- dplyr::bind_cols(.y = y, X, Z)
null_mod <- SKAT::SKAT_Null_Model(
  .y ~ ., data = null_train_df, out_type = out_type
)
skato_out <- purrr::map(
  snp_sets |> dplyr::pull(snp_set),
  function(snp_set) {
    skat_out <- SKAT::SKAT(
      Z = as.matrix(Z[, snp_set]), obj = null_mod, method = "SKATO"
    )
    skat_out$p.value
  }
) |>
  setNames(snp_sets |> dplyr::pull(Gene))
skat_out <- purrr::map(
  snp_sets |> dplyr::pull(snp_set),
  function(snp_set) {
    skat_out <- SKAT::SKAT(
      Z = as.matrix(Z[, snp_set]), obj = null_mod
    )
    skat_out$p.value
  }
) |>
  setNames(snp_sets |> dplyr::pull(Gene))

gene_pos_df <- snps_df |>
  dplyr::group_by(Gene) |>
  dplyr::summarise(
    Chr = Chr[1],
    min_pos = min(`hg19 Pos`),
    max_pos = max(`hg19 Pos`),
    mean_pos = mean(`hg19 Pos`),
    median_pos = median(`hg19 Pos`)
  )

skat_df <- dplyr::bind_rows(
  tibble::tibble(
    Gene = names(skato_out),
    pval = unlist(skato_out),
    mode = "SKATO"
  ),
  tibble::tibble(
    Gene = names(skat_out),
    pval = unlist(skat_out),
    mode = "SKAT"
  )
) |>
  dplyr::left_join(gene_pos_df, by = "Gene")

save(
  skato_out, skat_out, skat_df,
  file = file.path(RESULTS_DIR, "skat_results.Rdata")
)

#### MAGMA ####
# Run using FUMA