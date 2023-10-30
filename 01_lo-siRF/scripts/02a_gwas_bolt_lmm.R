rm(list = ls())
library(tidyverse)
library(data.table)

source(file.path("..", "functions", "load-functions.R"), chdir = TRUE)
source(file.path("..", "functions", "utils-gwas.R"), chdir = TRUE)

# set file paths
OUT_DIR <- file.path("..", "results", "gwas_bolt_lmm")
DATA_DIR <- file.path("..", "data")
BIM_DATA_PATH <- file.path(DATA_DIR, "ukbb_wbr_imp_morgan_chr")
FAM_DATA_PATH <- file.path(DATA_DIR, "ukbb_wbr_imp_morgan_chr1.fam")
BOLT_LMM_PATH <- file.path("..", "software", "BOLT-LMM_v2.3.4")

args <- commandArgs(trailingOnly = TRUE)
pheno_name <- as.character(args[1]) # name of phenotype

# create output directory
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

########################### Format Data for BOLT-LMM ###########################

print("Formatting input data for BOLT-LMM...")
if (!file.exists(file.path(OUT_DIR, "nontraining_ids.remove"))) {
  pheno_binary <- FALSE
  n <- 15000
  n_test <- 5000
  keep_genes <- c("TTN")
  include_dems <- TRUE
  include_npcs <- 5

  # read in training data to get sample ids
  snp_data <- loadData(
    pheno_name = pheno_name, pheno_binary = pheno_binary,
    n = n, n_test = n_test,
    keep_genes = keep_genes,
    include_dems = include_dems, include_npcs = include_npcs
  )
  # str(snp_data)

  # remove non-training ids
  geno_ids <- fread(FAM_DATA_PATH)
  rm_ids <- geno_ids[!(geno_ids$V1 %in% rownames(snp_data$geno_train)), 1:2]
  write.table(rm_ids, file.path(OUT_DIR, "nontraining_ids.remove"),
    quote = F, row.names = F, col.names = F
  )

  # get pheno file
  covar_df <- fread(FAM_DATA_PATH) %>%
    select(FID = V1, IID = V2) %>%
    cbind(
      ., loadDemographics(.$IID), loadLVDemographics(.$IID),
      loadPCs(.$IID, npcs = 5)
    ) %>%
    mutate(
      LVM = loadPheno(geno_id = IID, pheno_name = "LVM"),
      iLVM = loadPheno(geno_id = IID, pheno_name = "iLVM"),
      LVM_norm = loadPheno(geno_id = IID, pheno_name = "LVM_norm"),
      iLVM_norm = loadPheno(geno_id = IID, pheno_name = "iLVM_norm")
    )
  write.table(covar_df, file.path(OUT_DIR, "pheno.covars"),
    quote = F, row.names = F, col.names = T
  )
}

print("Completed")

################################ Run BOLT-LMM ##################################
bolt_lmm <- file.path(BOLT_LMM_PATH, "bolt")
remove_ids <- file.path(OUT_DIR, "nontraining_ids.remove")
pheno_file <- file.path(OUT_DIR, "pheno.covars")
covar_file <- file.path(OUT_DIR, "pheno.covars")
ld_scores <- file.path(BOLT_LMM_PATH, "tables", "LDSCORE.1000G_EUR.tab.gz")
out_file <- file.path(OUT_DIR, pheno_name)

print("Running BOLT-LMM...")
system(paste0(
  bolt_lmm,
  " --bim=", BIM_DATA_PATH, "{1:22}.bim",
  " --bed=", BIM_DATA_PATH, "{1:22}.bed",
  " --fam=", BIM_DATA_PATH, "1.fam",
  " --remove=", remove_ids,
  " --phenoFile=", pheno_file,
  " --phenoCol=", pheno_name,
  " --covarFile=", covar_file,
  " --covarCol=gender",
  " --qCovarCol=age",
  " --qCovarCol=weight",
  " --qCovarCol=height",
  " --qCovarCol=PC{1:5}",
  " --lmm",
  " --LDscoresFile=", ld_scores,
  " --numThreads=10",
  " --statsFile=", out_file, ".stats",
  " --maxModelSnps=15300000",
  " 2>&1 | tee ", out_file, ".log"
))
print("Completed")

######################## Format BOLT-LMM GWAS Results ##########################

print("Formatting GWAS results...")
gwas <- fread(paste0(out_file, ".stats"))
str(gwas)
gwas_annot <- annotGWAS(gwas,
  gwas_type = "bolt-lmm", save = TRUE,
  save_path = paste0(out_file, ".stats_annot")
)
print("Completed.")
