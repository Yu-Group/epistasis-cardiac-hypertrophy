rm(list = ls())
library(tidyverse)
library(data.table)

source(file.path("..", "functions", "load-functions.R"), chdir = TRUE)
source(file.path("..", "functions", "utils-gwas.R"), chdir = TRUE)

# set file paths
OUT_DIR <- file.path("..", "results", "gwas_plink")
DATA_DIR <- file.path("..", "data")
BIM_DATA_PATH <- file.path(DATA_DIR, "ukbb_wbr_imp_morgan_chr")
FAM_DATA_PATH <- file.path(DATA_DIR, "ukbb_wbr_imp_morgan_chr1.fam")

args <- commandArgs(trailingOnly = TRUE)
pheno_name <- as.character(args[1]) # name of phenotype
chr <- as.numeric(args[2]) # chromosome

# create output directory
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

########################## Format Data for Plink GWAS ##########################

print("Formatting input data for plink...")
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

  # get pheno files
  pheno_names <- c("LVM_norm", "iLVM_norm", "hiLVM_norm", "lmiLVM_norm")
  for (pname in pheno_names) {
    pheno_df <- fread(FAM_DATA_PATH) %>%
      select(FID = V1, IID = V2) %>%
      mutate(pheno = loadPheno(geno_id = IID, pheno_name = pname)) %>%
      filter(!is.na(pheno)) %>%
      setNames(c("FID", "IID", pname))
    write.table(pheno_df,
      file.path(OUT_DIR, paste0(pname, ".phe")),
      quote = F, row.names = F, col.names = T
    )
  }

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
print("Completed.")

############################### Run Plink GWAS #################################
remove_ids <- file.path(OUT_DIR, "nontraining_ids.remove")
pheno_file <- file.path(OUT_DIR, paste0(pheno_name, ".phe"))
covar_file <- file.path(OUT_DIR, "pheno.covars")
out_file <- file.path(OUT_DIR, paste0("chr", chr))

print("Running GWAS via plink...")
system(paste0(
  "plink2",
  " --bfile ", BIM_DATA_PATH, chr,
  " --glm firth-fallback hide-covar",
  " --maf 0.01",
  " --remove ", remove_ids,
  " --pheno ", pheno_file,
  " --out ", out_file,
  " --covar ", covar_file,
  " --covar-name gender age weight height PC1-PC5",
  " --covar-variance-standardize"
))
print("Completed.")

########################## Format Plink GWAS Results ###########################

out_files <- file.path(
  OUT_DIR,
  paste0("chr", 1:22, ".", pheno_name, ".glm.linear")
)
if (all(file.exists(out_files))) { # wait until all chr. have been processed
  print("Formatting GWAS results...")
  gwas <- map_dfr(1:22, ~ fread(out_files[.x]))
  str(gwas)
  gwas_annot <- annotGWAS(gwas,
    gwas_type = "plink", save = TRUE,
    save_path = file.path(
      OUT_DIR,
      paste0(pheno_name, ".stats_annot")
    )
  )
  print("Completed.")
}
