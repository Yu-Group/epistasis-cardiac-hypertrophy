source("fit-functions.R", chdir = TRUE)
source("predict-functions.R", chdir = TRUE)

#' Fits and predicts from various models (LASSO, ridge, kernel ridge, RF,
#' iRF, etc.) using training Predixcan or SNP data for a given phenotype.
#'
#' @param pheno_name Character string; phenotype name.
#' @param pheno_binary Logical; whether or not phenotype is binary.
#' @param n Number of training samples; if phenotype is binary, this can be a
#'   vector of length 2 with the # of positive and negative cases respectively.
#' @param n_test Number of test samples; if phenotype is binary, this can be a
#'   vector of length 2 with the # of positive and negative cases respectively
#' @param models Vector of models to fit; if NULL, fit all models that have been
#'   implemented.
#' @param predx NULL or a character string (e.g., "ao" or "skin"); if NULL,
#'   predixcan is not used; otherwise, specifies whether to use imputed gene
#'   expression from the aorta, skin, or some other tissue.
#' @param include_dems Logical; whether to include age and gender in training
#'   data.
#' @param include_npcs Number of pcs to include in training data.
#' @param keep_genes (Optional) character vector with gene names to keep;
#'   used if predx = NULL.
#' @param keep_snps (Optional) character vector with snp names to keep; used
#'   if predx = NULL.
#' @param n_trees Number of trees for ranger and iRF prediction models.
#' @param n_trials Number of bootstrap trials to run with different training 
#'   control samples if pheno_binary = TRUE.
#' @param save_path Patth to folder, where results are to be saved.
#' @param n_cores Integer; number of cores.
runAnalysis <- function(pheno_name, pheno_binary, n, n_test,
                        models = NULL, predx = NULL,
                        include_dems = T, include_npcs = 0,
                        keep_genes = NULL, keep_snps = NULL,
                        n_trees = 500, n_trials = 1,
                        save_path = "results", n_cores = 1) {
  fitModels(
    pheno_name = pheno_name, pheno_binary = pheno_binary,
    n = n, models = models, predx = predx,
    include_dems = include_dems, include_npcs = include_npcs,
    keep_genes = keep_genes, keep_snps = keep_snps,
    n_trees = n_trees, n_trials = n_trials, save_path = save_path,
    n_cores = 1
  )

  predictModels(
    pheno_name = pheno_name, pheno_binary = pheno_binary,
    n = n, n_test = n_test, models = models, predx = predx,
    include_dems = include_dems, include_npcs = include_npcs,
    keep_genes = keep_genes, keep_snps = keep_snps,
    n_trials = n_trials, save_path = save_path
  )
}
