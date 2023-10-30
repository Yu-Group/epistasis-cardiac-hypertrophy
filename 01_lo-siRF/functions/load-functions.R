source("utils.R", chdir = TRUE)

#' Main function to load in training/test X and y data
#'
#' @inheritParams runAnalysis
loadData <- function(pheno_name, pheno_binary = T, n, n_test,
                     predx = NULL, include_dems = T, include_npcs = 0,
                     keep_genes = NULL, keep_snps = NULL) {
  options(stringsAsFactors = FALSE)

  if (is.null(predx)) { # read in SNP data
    snp_data <- loadSNPData(
      pheno_name = pheno_name,
      pheno_binary = pheno_binary,
      n = n, n_test = n_test,
      keep_snps = keep_snps,
      keep_genes = keep_genes
    )
    geno_train <- snp_data$geno
    pheno_train <- snp_data$pheno
    geno_test <- snp_data$geno_test
    pheno_test <- snp_data$pheno_test
    rm(snp_data)
  } else { # read in predixcan data
    # get predixcan paths
    PREDX_DATA_DIR <- file.path("..", "data", "ukbb_predixcan")
    path_excl <- file.path(PREDX_DATA_DIR, "geneExclusionList.txt")
    if (predx == "ao") {
      path_predx <- file.path(PREDX_DATA_DIR, "aorta_predicted_expression.txt")
    } else if (predx == "skin") {
      path_predx <- file.path(PREDX_DATA_DIR, "skin_predicted_expression.txt")
    } else if (predx == "lv") {
      path_predx <- file.path(PREDX_DATA_DIR, "lv_predicted_expression.txt")
    } else if (predx == "aa") {
      path_predx <- file.path(PREDX_DATA_DIR, "aa_predicted_expression.txt")
    } else if (predx == "lv_utmost") {
      path_predx <- file.path(
        PREDX_DATA_DIR, "lv_predicted_expression_utmost.txt"
      )
    }
    train_data <- loadTrainData(
      pheno_name = pheno_name,
      pheno_binary = pheno_binary,
      n = n,
      predx = predx,
      path_predx = path_predx,
      path_excl = path_excl,
      keep_genes = keep_genes
    )
    geno_train <- train_data$geno
    pheno_train <- train_data$pheno
    rm(train_data)
    if (n_test[1] > 0) {
      test_data <- loadTestData(
        pheno_name = pheno_name,
        pheno_binary = pheno_binary,
        n = n, n_test = n_test,
        predx = predx,
        path_predx = path_predx,
        path_excl = path_excl,
        keep_genes = keep_genes
      )
      geno_test <- test_data$geno_test
      pheno_test <- test_data$pheno_test
      rm(test_data)
    } else {
      geno_test <- NULL
      pheno_test <- NULL
    }
  }

  if (stringr::str_detect(pheno_name, "binary_thr")) { # binarize phenotype
    thr <- stringr::str_extract(pheno_name, "(?<=binary_thr)(\\d*\\.?\\d*)") %>%
      as.numeric()
    gender_train <- loadDemographics(rownames(geno_train))[, "gender"]
    if (is.null(geno_test)) {
      gender_test <- NULL
    } else {
      gender_test <- loadDemographics(rownames(geno_test))[, "gender"]
    }
    binned_data <- binarizePhenotype(
      pheno_train = pheno_train,
      gender_train = gender_train,
      pheno_test = pheno_test,
      gender_test = gender_test,
      low_q = thr, hi_q = thr
    )
    print("Low Thresholds:")
    print(binned_data$low_thr)
    print("High Thresholds:")
    print(binned_data$hi_thr)
    tmp <- removeNASamples(
      geno_train, binned_data$pheno_train,
      geno_test, binned_data$pheno_test
    )
    geno_train <- tmp$geno_train
    geno_test <- tmp$geno_test
    pheno_train <- tmp$pheno_train
    pheno_test <- tmp$pheno_test
    rm(tmp)
  }

  if (include_dems) {
    geno_train <- cbind(loadDemographics(rownames(geno_train)), geno_train)
    # remove NAs
    na_idx <- is.na(geno_train[, "gender"])
    geno_train <- geno_train[!na_idx, ]
    pheno_train <- pheno_train[!na_idx]
    if (stringr::str_detect(pheno_name, "LV")) {
      geno_train <- cbind(loadLVDemographics(rownames(geno_train)), geno_train)
      # remove NAs
      na_idx <- is.na(geno_train[, "weight"]) | is.na(geno_train[, "height"])
      geno_train <- geno_train[!na_idx, ]
      pheno_train <- pheno_train[!na_idx]
    }
  }
  if (include_npcs > 0) {
    geno_train <- cbind(loadPCs(rownames(geno_train), include_npcs), geno_train)
  }
  geno_train <- as.matrix(geno_train)

  if (!is.null(geno_test)) {
    if (include_dems) {
      geno_test <- cbind(loadDemographics(rownames(geno_test)), geno_test)
      # remove NAs
      na_idx <- is.na(geno_test[, "gender"])
      geno_test <- geno_test[!na_idx, ]
      pheno_test <- pheno_test[!na_idx]
      if (stringr::str_detect(pheno_name, "LV")) {
        geno_test <- cbind(loadLVDemographics(rownames(geno_test)), geno_test)
        # remove NAs
        na_idx <- is.na(geno_test[, "weight"]) | is.na(geno_test[, "height"])
        geno_test <- geno_test[!na_idx, ]
        pheno_test <- pheno_test[!na_idx]
      }
    }
    if (include_npcs > 0) {
      geno_test <- cbind(loadPCs(rownames(geno_test), include_npcs), geno_test)
    }
    geno_test <- as.matrix(geno_test)
  }

  # remove constant columns
  const_cols <- apply(geno_train, 2, var) == 0
  print(paste0("Number of constant columns: ", sum(const_cols)))
  geno_train <- geno_train[, !const_cols]

  # remove duplicated columns
  geno_train_t <- t(geno_train)
  print(paste0("Number of duplicated columns: ", sum(duplicated(geno_train_t))))
  geno_train_t <- geno_train_t[!duplicated(geno_train_t), ]
  geno_train <- t(geno_train_t)

  names(pheno_train) <- rownames(geno_train)
  geno_train <- as.matrix(geno_train)

  if (!is.null(geno_test)) {
    geno_test <- geno_test[, colnames(geno_train)]
    names(pheno_test) <- rownames(geno_test)
    geno_test <- as.matrix(geno_test)
  }

  return(list(
    geno_train = geno_train,
    pheno_train = pheno_train,
    geno_test = geno_test,
    pheno_test = pheno_test
  ))
}

#' Load SNPs and phenotype training data
#'
#' @inheritParams runAnalysis
#' @param verbose Logical; whether or not to print summary stats of loaded data.
#'
#' @returns A list of 4:
#' \describe{
#' \item{geno}{SNP training data matrix.}
#' \item{pheno}{Phenotype training data vector.}
#' \item{geno_test}{SNP test data matrix.}
#' \item{pheno_test}{Phenotype test data vector.}
#' }
loadSNPData <- function(pheno_name, pheno_binary = FALSE, n, n_test = 0,
                        keep_snps = NULL, keep_genes = NULL, verbose = TRUE) {
  print("Loading data")

  # get data file paths
  BIM_DATA_PATH <- file.path(
    "..", "data", paste0("ukbb_wbr_imp_morgan_chr", 1:22)
  )
  FAM_DATA_PATH <- file.path("..", "data", "ukbb_wbr_imp_morgan_chr1.fam")
  SNPINFO_DATA_PATH <- file.path("..", "results", "annovar", "snp2gene_df.csv")

  # get genotype data
  if (is.numeric(keep_genes)) { # get data for selected chromosomes
    chrs <- keep_genes
    geno_ls <- list()
    for (chr in chrs) {
      snps <- read.table(paste0(BIM_DATA_PATH[chr], ".bim"))[, 2]
      is_dup <- duplicated(snps)
      geno <- snpStats::read.plink(
        BIM_DATA_PATH[chr],
        select.snps = which(!is_dup)
      )$genotypes
      colnames(geno) <- purrr::map_chr(colnames(geno), ~ renameSNP(.x, chr))
    }
  } else { # get data for selected genes or snps
    if (!is.null(keep_genes)) {
      # get SNP information data frame (snp -> loci mapping)
      snps_df <- data.table::fread(SNPINFO_DATA_PATH)
      # get list of SNPs to keep
      keep_snps_vec <- snps_df %>%
        dplyr::filter(Gene %in% keep_genes) %>%
        dplyr::pull(rsID)
    } else if (!is.null(keep_snps)) {
      # get list of SNPs to keep
      keep_snps_vec <- keep_snps
    } else {
      stop("Must provide either keep_genes or keep_snps")
    }
    print(paste0("Number of SNPs in specified genes: ", length(keep_snps_vec)))

    # generate name of scratch file to write temporary snplist
    flag <- TRUE
    while (flag) {
      rand_id <- sample(1:10000, 1)
      path_snplist <- file.path("scratch", paste0("keep.snplist_", rand_id))
      if (!file.exists(path_snplist)) {
        flag <- FALSE
      }
    }
    write.table(keep_snps_vec,
      file = path_snplist,
      row.names = F, col.names = F, quote = F
    )
    snp_files <- file.path("scratch", paste0("snps_chr", 1:22, "_", rand_id))

    # get data for selected snps
    geno_ls <- list()
    for (chr in 1:22) {
      # subset to specified snps using plink
      system(paste0(
        "plink --bfile ", BIM_DATA_PATH[chr],
        " --indiv-sort f ", FAM_DATA_PATH,
        " --extract ", path_snplist,
        " --make-bed --out ", snp_files[chr]
      ))
      # load in specified snps
      if (file.exists(paste0(snp_files[chr], ".bim"))) {
        snps <- read.table(paste0(snp_files[chr], ".bim"))[, 2]
        is_dup <- duplicated(snps)
        geno <- snpStats::read.plink(
          snp_files[chr],
          select.snps = which(!is_dup)
        )$genotypes
        colnames(geno) <- purrr::map_chr(colnames(geno), ~ renameSNP(.x, chr))
        rnames <- rownames(geno)
        geno_ls[[chr]] <- geno
      }
      system(paste0("rm ", snp_files[chr], "*"))
    }
    file.remove(path_snplist)
  }

  if (verbose) {
    str(geno_ls)
  }

  # get phenotype data
  pheno <- loadPheno(geno_id = rnames, pheno_name = pheno_name)

  # remove samples with NA phenotype
  if (any(is.na(pheno))) {
    # print(sum(is.na(pheno)))
    geno_ls <- lapply(geno_ls, function(geno) geno[!is.na(pheno), ])
    pheno <- pheno[!is.na(pheno)]
  }

  # get desired number of samples
  if (pheno_binary) {
    n_pos <- n[1]
    n_pos_test <- n_test[1]
    if (length(n) > 1) {
      n_neg <- n[2]
    } else {
      n_neg <- n_pos
    }
    if (length(n_test) > 1) {
      n_neg_test <- n_test[2]
    } else {
      n_neg_test <- n_pos_test
    }
    num_pos <- sum(pheno == 1, na.rm = T)
    num_neg <- sum(pheno == 0, na.rm = T)
    if (num_pos < (n_pos + n_pos_test)) {
      stop(sprintf(
        "Number of positive cases did not reach %s",
        n_pos + n_pos_test
      ))
    } else if (num_neg < (n_neg + n_neg_test)) {
      stop(sprintf(
        "Number of negative cases did not reach %s",
        n_neg + n_neg_test
      ))
    }
    train_id <- c(which(pheno == 0)[1:n_neg], which(pheno == 1)[1:n_pos])
    test_id <- c(
      which(pheno == 0)[(n_neg + 1):(n_neg + n_neg_test)],
      which(pheno == 1)[(n_pos + 1):(n_pos + n_pos_test)]
    )
  } else {
    train_id <- 1:n
    test_id <- (n + 1):min(n + n_test, length(pheno))
  }

  # get test data
  if (n_test[1] > 0) {
    geno_test_ls <- lapply(geno_ls, function(geno) geno[test_id, ])
    pheno_test <- pheno[test_id]

    # clean data
    geno_test <- purrr::reduce(geno_test_ls, cbind)
    gnames <- rownames(geno_test)
    geno_test <- apply(geno_test, 2, as.numeric)
    rownames(geno_test) <- gnames

    # median imputation for missing values
    geno_test[geno_test == 0] <- NA
    geno_test <- 3 - geno_test # 0 = no mutation, 2 = two mutations
    geno_test <- apply(geno_test, 2,
      FUN = function(x) {
        if (any(is.na(x))) {
          x[is.na(x)] <- round(median(x, na.rm = T))
        }
        return(x)
      }
    )
  } else {
    geno_test <- NULL
    pheno_test <- NULL
  }

  # get train data
  geno_ls <- lapply(geno_ls, function(geno) geno[train_id, ])
  pheno <- pheno[train_id]

  # clean data
  geno <- purrr::reduce(geno_ls, cbind)
  gnames <- rownames(geno)
  geno <- apply(geno, 2, as.numeric)
  rownames(geno) <- gnames

  # median imputation for missing values
  geno[geno == 0] <- NA
  geno <- 3 - geno # 0 = no mutation, 2 = two mutations
  geno <- apply(geno, 2,
    FUN = function(x) {
      if (any(is.na(x))) {
        x[is.na(x)] <- round(median(x, na.rm = T))
      }
      return(x)
    }
  )

  # remove duplicate columns
  geno_t <- t(geno)
  geno_t <- geno_t[!duplicated(geno_t), ]
  geno <- t(geno_t)

  # remove constant and na columns
  const_cols <- apply(geno, 2, var) == 0
  na_cols <- apply(geno, 2, function(x) all(is.na(x)))
  geno <- geno[, !const_cols & !na_cols]

  if (n_test > 0) {
    geno_test <- geno_test[, colnames(geno)]
  }

  if (verbose) {
    print(paste("memory used:"))
    print(pryr::mem_used())
    print(sprintf(
      "dim of geno: (%s, %s)",
      nrow(geno), ncol(geno)
    ))
    print(sprintf("length of pheno: %s", length(pheno)))
    if (n_test > 0) {
      print(sprintf(
        "dim of geno_test: (%s, %s)",
        nrow(geno_test), ncol(geno_test)
      ))
      print(sprintf("length of pheno_test: %s", length(pheno_test)))
    }
  }

  return(list(
    geno = geno, pheno = pheno,
    geno_test = geno_test, pheno_test = pheno_test
  ))
}

#' Load Predixcan imputed gene expression and phenotype training data
#'
#' @inheritParams runAnalysis
#'
#' @returns A list of 4:
#' \describe{
#' \item{geno}{Imputed gene expression training data matrix.}
#' \item{pheno}{Phenotype training data vector.}
#' }
loadTrainData <- function(pheno_name, pheno_binary = T, n,
                          predx, path_predx, path_excl, keep_genes = NULL,
                          verbose = T) {
  print("Load training data")

  n_pos <- n[1]
  if (pheno_binary) {
    if (length(n) > 1) {
      n_neg <- n[2]
    } else {
      n_neg <- n_pos
    }
  }

  ## load batch 1
  print("Load batch 1")
  load_id <- 1:100000
  geno <- loadPredixcan(
    load_id = load_id, path_predx = path_predx,
    path_excl = path_excl
  )
  pheno <- loadPheno(geno_id = rownames(geno), pheno_name = pheno_name)
  if (any(is.na(pheno))) {
    geno <- geno[!is.na(pheno), ]
    pheno <- pheno[!is.na(pheno)]
  }

  if (verbose) {
    print(paste("memory after loading training data batch 1:"))
    print(pryr::mem_used())
    print(paste0("dim of geno: (", nrow(geno), ", ", ncol(geno), ")"))
    print(paste("length of pheno:", length(pheno)))
    print(paste("mean of pheno:", mean(pheno)))
  }

  if (pheno_binary) {
    # get balanced or specified number of +/- training set if binary phenotype
    num_pos <- sum(pheno == 1)
    num_neg <- sum(pheno == 0)
    if (verbose) {
      print(paste("number of positive cases in batch 1", num_pos))
      print(paste("number of negative cases in batch 1", num_neg))
    }

    # get specified number of +/- samples if possible
    load_id <- c(
      which(pheno == 0)[1:min(n_neg, num_neg)],
      which(pheno == 1)[1:min(n_pos, num_pos)]
    )
    geno <- geno[load_id, ]
    pheno <- pheno[load_id]
    if (verbose) {
      print(paste0(
        "after balance, dim of geno: (",
        nrow(geno), ", ", ncol(geno), ")"
      ))
      print(paste("length of pheno", length(pheno)))
      print(paste("mean of pheno", mean(pheno)))
      print(paste(
        "after 1 batch we have", num_pos, "positive cases and",
        num_neg, "negative cases"
      ))
    }

    # load in more data if haven't reached number of desired positive cases
    batch <- 2
    while ((num_pos < n_pos | num_neg < n_neg) & batch <= 4) {
      print(paste("load batch", batch))
      geno_old <- geno
      pheno_old <- pheno
      num_pos_old <- num_pos
      num_neg_old <- num_neg

      # load in 100000 more samples (max number of samples: 337535)
      load_id <- ((batch - 1) * 100000 + 1):(min(batch * 100000, 337535))
      geno <- loadPredixcan(
        load_id = load_id, path_predx = path_predx,
        path_excl = path_excl
      )
      pheno <- loadPheno(geno_id = rownames(geno), pheno_name = pheno_name)
      if (any(is.na(pheno))) {
        geno <- geno[!is.na(pheno), ]
        pheno <- pheno[!is.na(pheno)]
      }

      # get balanced classes if possible
      num_pos <- sum(pheno == 1)
      num_neg <- sum(pheno == 0)
      load_id <- NULL
      if (num_pos_old < n_pos) {
        load_id <- c(
          load_id,
          which(pheno == 1)[1:min(n_pos - num_pos_old, num_pos)]
        )
      }
      if (num_neg_old < n_neg) {
        load_id <- c(
          load_id,
          which(pheno == 0)[1:min(n_neg - num_neg_old, num_neg)]
        )
      }
      geno <- rbind(geno_old, geno[load_id, ])
      pheno <- c(pheno_old, pheno[load_id])
      num_pos <- sum(pheno == 1)
      num_neg <- sum(pheno == 0)

      if (verbose) {
        print(paste0("memory after loading training data batch ", batch, ": "))
        print(pryr::mem_used())
        print(paste(
          "after batch", batch, "we have", num_pos, "positive cases",
          "and", num_neg, "negative cases"
        ))
        print(paste0("dim of geno: (", nrow(geno), ", ", ncol(geno), ")"))
        print(paste("length of pheno:", length(pheno)))
        print(paste("mean of pheno:", mean(pheno)))
      }

      batch <- batch + 1
    }

    if (num_pos < n_pos) {
      stop(paste0("Number of positive cases did not reach ", n_pos))
    } else if (num_neg < n_neg) {
      stop(paste0("Number of negative cases did not reach ", n_neg))
    }
  } else {
    num_pos <- length(pheno)
    if (num_pos >= n_pos) {
      geno <- geno[1:n_pos, ]
      pheno <- pheno[1:n_pos]
    }

    # load in additional batches of data in order to get up to n_pos samples
    batch <- 2
    while ((num_pos < n_pos) & batch <= 4) {
      print(paste("load batch", batch))
      geno_old <- geno
      pheno_old <- pheno
      num_pos_old <- num_pos

      # load in 100000 more samples (max number of samples: 337535)
      load_id <- ((batch - 1) * 100000 + 1):(min(batch * 100000, 337535))
      geno <- loadPredixcan(
        load_id = load_id, path_predx = path_predx,
        path_excl = path_excl
      )
      pheno <- loadPheno(geno_id = rownames(geno), pheno_name = pheno_name)
      if (any(is.na(pheno))) {
        geno <- geno[!is.na(pheno), ]
        pheno <- pheno[!is.na(pheno)]
      }
      num_pos <- length(pheno)
      geno <- rbind(geno_old, geno[1:min(n_pos - num_pos_old, num_pos), ])
      pheno <- c(pheno_old, pheno[1:min(n_pos - num_pos_old, num_pos)])
      num_pos <- length(pheno)

      if (verbose) {
        print(paste("memory after loading training data batch ", batch, ": "))
        print(pryr::mem_used())
        print(paste0("dim of geno: (", nrow(geno), ", ", ncol(geno), ")"))
        print(paste("length of pheno:", length(pheno)))
        print(paste("mean of pheno:", mean(pheno)))
      }

      batch <- batch + 1
    }

    if (num_pos < n_pos) {
      stop(paste0("Number of cases did not reach ", n_pos))
    }
  }

  if (!is.null(keep_genes)) {
    if (stringr::str_detect(predx, "utmost")) {
      geno <- geno[, intersect(keep_genes, colnames(geno))]
    } else {
      genenames <- ENSGtoGene(colnames(geno), predx = predx)
      geno <- geno[, which(genenames %in% keep_genes)]
    }
  }

  return(list(geno = geno, pheno = pheno))
}

#' Load Predixcan imputed gene expression and phenotype test data
#'
#' @inheritParams runAnalysis
#'
#' @returns A list of 4:
#' \describe{
#' \item{geno_test}{Imputed gene expression test data matrix.}
#' \item{pheno_test}{Phenotype test data vector.}
#' }
loadTestData <- function(pheno_name, pheno_binary = T, n, n_test,
                         predx, path_predx, path_excl, keep_genes = NULL,
                         verbose = T) {
  print("Load test data")

  n_pos <- n[1]
  n_pos_test <- n_test[1]
  if (pheno_binary) {
    if (length(n) > 1) {
      n_neg <- n[2]
    } else {
      n_neg <- n_pos
    }
    if (length(n_test) > 1) {
      n_neg_test <- n_test[2]
    } else {
      n_neg_test <- n_pos_test
    }
  }

  ## load batch 1
  print("Load batch 1")
  load_id <- 1:100000
  geno <- loadPredixcan(
    load_id = load_id, path_predx = path_predx,
    path_excl = path_excl
  )
  pheno <- loadPheno(geno_id = rownames(geno), pheno_name = pheno_name)
  if (any(is.na(pheno))) {
    geno <- geno[!is.na(pheno), ]
    pheno <- pheno[!is.na(pheno)]
  }

  if (verbose) {
    print(paste("memory after loading training data batch 1:"))
    print(pryr::mem_used())
    print(paste0("dim of geno: (", nrow(geno), ", ", ncol(geno), ")"))
    print(paste("length of pheno:", length(pheno)))
    print(paste("mean of pheno:", mean(pheno)))
  }

  # get balanced training set if phenotype is binary
  if (pheno_binary) {
    n_neg_all <- n_neg + n_neg_test
    n_pos_all <- n_pos + n_pos_test

    num_pos <- sum(pheno == 1)
    num_neg <- sum(pheno == 0)
    if (verbose) {
      print(paste("total number of positive cases in batch 1", num_pos))
      print(paste("total number of negative cases in batch 1", num_neg))
    }

    # get balanced classes if possible
    load_id <- NULL
    if (num_pos > n_pos) {
      load_id <- c(
        load_id,
        which(pheno == 1)[(n_pos + 1):min(n_pos_all, num_pos)]
      )
    }
    if (num_neg > n_neg) {
      load_id <- c(
        load_id,
        which(pheno == 0)[(n_neg + 1):min(n_neg_all, num_neg)]
      )
    }

    if (is.null(load_id)) {
      geno <- NULL
      pheno <- NULL
    } else {
      geno <- geno[load_id, ]
      pheno <- pheno[load_id]
    }

    if (verbose) {
      print(paste0(
        "after balance, dim of geno: (",
        nrow(geno), ", ", ncol(geno), ")"
      ))
      print(paste("length of pheno", length(pheno)))
      print(paste("mean of pheno", mean(pheno)))
      print(paste(
        "after 1 batch we have", sum(pheno == 1),
        "positive test cases and", sum(pheno == 0),
        "negative test cases"
      ))
    }

    # load in more data if haven't reached number of desired positive cases
    batch <- 2
    while ((num_pos < n_pos_all | num_neg < n_neg_all) & batch <= 4) {
      print(paste("load batch", batch))
      geno_old <- geno
      pheno_old <- pheno
      num_pos_old <- num_pos
      num_neg_old <- num_neg

      # load in 100000 more samples (max number of samples: 337535)
      load_id <- ((batch - 1) * 100000 + 1):(min(batch * 100000, 337535))
      geno <- loadPredixcan(
        load_id = load_id, path_predx = path_predx,
        path_excl = path_excl
      )
      pheno <- loadPheno(geno_id = rownames(geno), pheno_name = pheno_name)
      if (any(is.na(pheno))) {
        geno <- geno[!is.na(pheno), ]
        pheno <- pheno[!is.na(pheno)]
      }

      # get balanced classes if possible
      num_pos <- sum(pheno == 1)
      num_neg <- sum(pheno == 0)
      load_id <- NULL
      if ((num_pos_old + num_pos > n_pos) & (num_pos_old < n_pos_all)) {
        load_id <- c(
          load_id,
          which(pheno == 1)[max(n_pos + 1 - num_pos_old, 1):
          min(n_pos_all - num_pos_old, num_pos)]
        )
      }
      if ((num_neg_old + num_neg > n_neg) & (num_neg_old < n_neg_all)) {
        load_id <- c(
          load_id,
          which(pheno == 0)[max(n_neg + 1 - num_neg_old, 1):
          min(n_neg_all - num_neg_old, num_neg)]
        )
      }

      if (is.null(load_id)) {
        geno <- geno_old
        pheno <- pheno_old
      } else {
        geno <- rbind(geno_old, geno[load_id, ])
        pheno <- c(pheno_old, pheno[load_id])
      }

      num_pos <- num_pos + num_pos_old
      num_neg <- num_neg + num_neg_old

      if (verbose) {
        print(paste0("memory after loading training data batch ", batch, ": "))
        print(pryr::mem_used())
        print(paste(
          "after batch", batch, "we have", sum(pheno == 1),
          "positive test cases and", sum(pheno == 0),
          "negative test cases"
        ))
        print(paste0("dim of geno: (", nrow(geno), ", ", ncol(geno), ")"))
        print(paste("length of pheno:", length(pheno)))
        print(paste("mean of pheno:", mean(pheno)))
      }

      batch <- batch + 1
    }

    if (num_pos < n_pos_all) {
      stop(paste0("Number of positive test cases did not reach ", n_pos_test))
    } else if (num_neg < n_neg_all) {
      stop(paste0("Number of negative test cases did not reach ", n_neg_test))
    }
  } else {
    n_pos_all <- n_pos + n_pos_test

    num_pos <- length(pheno)
    if (num_pos > n_pos) {
      geno <- geno[(n_pos + 1):min(n_pos_all, num_pos), ]
      pheno <- pheno[(n_pos + 1):min(n_pos_all, num_pos)]
    } else {
      geno <- NULL
      pheno <- NULL
    }

    # load in additional batches of data in order to get up to n_pos samples
    batch <- 2
    while ((num_pos < n_pos_all) & batch <= 4) {
      print(paste("load batch", batch))
      geno_old <- geno
      pheno_old <- pheno
      num_pos_old <- num_pos

      # load in 100000 more samples (max number of samples: 337535)
      load_id <- ((batch - 1) * 100000 + 1):(min(batch * 100000, 337535))
      geno <- loadPredixcan(
        load_id = load_id, path_predx = path_predx,
        path_excl = path_excl
      )
      pheno <- loadPheno(geno_id = rownames(geno), pheno_name = pheno_name)
      if (any(is.na(pheno))) {
        geno <- geno[!is.na(pheno), ]
        pheno <- pheno[!is.na(pheno)]
      }

      num_pos <- length(pheno)
      if (num_pos_old + num_pos > n_pos) {
        geno <- geno[max(n_pos + 1 - num_pos_old, 1):
        min(n_pos_all - num_pos_old, num_pos), ]
        pheno <- pheno[max(n_pos + 1 - num_pos_old, 1):
        min(n_pos_all - num_pos_old, num_pos)]
        geno <- rbind(geno_old, geno)
        pheno <- c(pheno_old, pheno)
      } else {
        geno <- geno_old
        pheno <- pheno_old
      }

      num_pos <- num_pos + num_pos_old

      if (verbose) {
        print(paste("memory after loading training data batch ", batch, ": "))
        print(pryr::mem_used())
        print(paste0("dim of geno: (", nrow(geno), ", ", ncol(geno), ")"))
        print(paste("length of pheno:", length(pheno)))
        print(paste("mean of pheno:", mean(pheno)))
      }

      batch <- batch + 1
    }

    if (num_pos < n_pos_all) {
      stop(paste0("Number of cases did not reach ", n_pos_test))
    }
  }

  if (!is.null(keep_genes)) {
    if (stringr::str_detect(predx, "utmost")) {
      geno <- geno[, intersect(keep_genes, colnames(geno))]
    } else {
      genenames <- ENSGtoGene(colnames(geno), predx = predx)
      geno <- geno[, which(genenames %in% keep_genes)]
    }
  }

  return(list(geno_test = geno, pheno_test = pheno))
}

#' Load in Predixcan genotype data
#'
#' @param load_id Vector of sample indices.
#' @param path_predx Path to PrediXcan data.
#' @param path_excl Path to file with genes to exclude
loadPredixcan <- function(load_id, path_predx, path_excl) {
  options(stringsAsFactors = FALSE)

  if (missing(load_id)) stop('Must specify indices to load as: "load_id"')
  if (missing(path_predx)) stop('Must specify PrediXcan path as: "path_predx"')
  if (all(diff(load_id) != 1)) stop("Only consecutive indices supported")


  # Load gene identifiers
  gnames <- data.table::fread(path_predx[1], nrow = 1, header = F)
  gnames <- as.character(gnames)[-c(1, 2)]

  # Load genes to be excluded
  if (!missing(path_excl)) {
    print("load gene-exclusion-list")
    gene_ex <- data.table::fread(path_excl, skip = 1, header = F)
    gene_ex <- gene_ex[[1]]

    ind_ex <- fastmatch::fmatch(gene_ex, gnames)
    ind_ex <- ind_ex[!is.na(ind_ex)]
  }

  skip <- min(load_id) - 1 + 4
  N <- max(load_id) - (skip - 4)

  geno <- 0
  active <- numeric(0)
  gchrom <- list()

  for (i in 1:length(path_predx)) {
    print(paste("load PrediXcan file", i))
    xin <- data.table::fread(path_predx[i], skip = skip, nrow = N)
    if (i == 1) {
      print("load subject names")
      subnames <- as.character(unlist(xin[, 2]))
    }
    if (!identical(subnames, as.character(unlist(xin[, 2])))) {
      warning("Subjects in files not identical!")
    }
    xin <- Matrix::Matrix(as.matrix(xin[, -(1:2)]), sparse = TRUE)
    # Update genotype matrix with gene expression from current chromosome
    geno <- geno + xin
    gc(verbose = FALSE)
  }

  # Remove genes to be excluded
  if (!missing(path_excl)) {
    print(paste("Removing", length(ind_ex), "out of", ncol(geno), "genes"))
    if (length(ind_ex) > 0) {
      print(paste("Removed genes:", gnames[ind_ex]))
      geno <- geno[, -ind_ex]
      gnames <- gnames[-ind_ex]
    }
  }

  colnames(geno) <- gnames
  rownames(geno) <- subnames
  return(geno)
}

#' Load in phenotype
#'
#' @param geno_id Vector of patient IDs from genotype files, e.g., 
#'   rownames(geno).
#' @param pheno_name Character string; phenotype name.
#'
#' @returns A vector of phenotypes in the same order as geno_id.
loadPheno <- function(geno_id, pheno_name) {
  options(stringsAsFactors = FALSE)

  if (missing(geno_id)) stop("Must provide geno_id")
  if (missing(pheno_name)) stop("Must provide pheno_name")

  if (stringr::str_detect(pheno_name, "HCM")) {
    #### HCM ####
    PHENO_PATH <- file.path("..", "data", "ukbb_icd_diagnoses_data.csv")

    # load in ukb phenotype data
    pheno_orig <- data.table::fread(PHENO_PATH)

    # icd codes of interest
    icd10_codes_keep <- c("I421", "I422")
    icd9_codes_keep <- c("4251")
    icd10_codes_omit <- NULL
    icd9_codes_omit <- NULL

    # ids of positive phenotype
    icd_out <- ICD2Pheno(pheno_orig,
      icd10_codes_keep = icd10_codes_keep,
      icd9_codes_keep = icd9_codes_keep,
      icd10_codes_omit = icd10_codes_omit,
      icd9_codes_omit = icd9_codes_omit,
      controls_only = T
    )
    positive_ids <- icd_out$positive_ids
    omit_ids <- icd_out$omit_ids

    # merge phenotypes with geno_id
    pheno_df <- data.frame(id = geno_id, pheno = 0)
    pheno_df[pheno_df$id %in% omit_ids, "pheno"] <- NA
    pheno_df[pheno_df$id %in% positive_ids, "pheno"] <- 1
    pheno <- pheno_df$pheno
  } else if (stringr::str_detect(pheno_name, "LV")) {
    #### LV ####
    PHENO_PATH <- file.path(
      "..", "data", "table_ventricular_volume_with_indexing.csv"
    )

    # load in cardiac mri-derived phenotype data
    pheno_orig <- data.table::fread(PHENO_PATH)
    pheno_df <- dplyr::left_join(
      x = data.frame(id = as.numeric(geno_id)), y = pheno_orig, by = "id"
    )
    if (stringr::str_detect(pheno_name, "hiLVM")) {
      pheno_df$pheno <- pheno_df$`hiLVM (g/cm)`
    } else if (stringr::str_detect(pheno_name, "iLVM")) {
      pheno_df$pheno <- pheno_df$`iLVM (g/m2)`
    } else if (stringr::str_detect(pheno_name, "LVM")) {
      pheno_df$pheno <- pheno_df$`LVM (g)`
    } else if (stringr::str_detect(pheno_name, "iLVEDV")) {
      pheno_df$pheno <- pheno_df$`LVEDV (mL)` /
        pheno_df$`body surface area estimate (m2)`
    } else if (stringr::str_detect(pheno_name, "LVEDV")) {
      pheno_df$pheno <- pheno_df$`LVEDV (mL)`
    } else if (stringr::str_detect(pheno_name, "iLVESV")) {
      pheno_df$pheno <- pheno_df$`LVESV (mL)` /
        pheno_df$`body surface area estimate (m2)`
    } else if (stringr::str_detect(pheno_name, "LVESV")) {
      pheno_df$pheno <- pheno_df$`LVESV (mL)`
    } else if (stringr::str_detect(pheno_name, "LVSV")) {
      pheno_df$pheno <- pheno_df$`LVSV (mL)`
    } else if (stringr::str_detect(pheno_name, "LVEF")) {
      pheno_df$pheno <- pheno_df$`LVEF (%)`
    } else if (stringr::str_detect(pheno_name, "LVCO")) {
      pheno_df$pheno <- pheno_df$`LVCO (L/min)`
    }
    
    if (stringr::str_detect(pheno_name, "norm")) {
      # rank-based inverse normal transform (from RNOmni package)
      rankNorm <- function(u, k = 0.375) {
        # Input checks. 
        if (!is.vector(u)) {
          stop("A numeric vector is expected for u.")
        }
        if ((k < 0) || (k > 0.5)) {
          stop("Select the offset within the interval (0,0.5).")
        }
        
        na_idx <- is.na(u)
        
        # Observations.
        n <- sum(!na_idx)
        
        # Ranks.
        r <- rank(u[!na_idx])
        
        # Apply transformation.
        rt <- qnorm((r - k) / (n - 2 * k + 1))
        
        # Format output.
        out <- rep(NA, length(u))
        out[!na_idx] <- rt
        
        return(out)
      }
      pheno_df$pheno <- rankNorm(pheno_df$pheno)
    }
    
    pheno <- pheno_df$pheno
  }

  return(pheno)
}

#' Helper function for loadPheno(). Convets ICD9/10 codes to phenotype
#' indicators (TRUE, FALSE, NA).
#'
#' @param pheno_orig Phenotype dataframe with icd codes (original ukb data).
#' @param icd10_codes_keep ICD10 codes to keep.
#' @param icd10_codes_u60 ICD10 codes to keep if age < 60.
#' @param icd10_codes_omit ICD10 codes to omit.
#' @param icd9_codes_keep ICD9 codes to keep.
#' @param icd9_codes_u60 ICD9 codes to keep if age < 60.
#' @param icd9_codes_omit ICD9 codes to omit.
#' @param controls_only Logical; whether to omit codes for controls only.
#'
#' @returns If \code{controls_only == FALSE}, returns a vector of patient IDs
#'   in the positive group. Otherwise, returns a list of two:
#'   \describe{
#'   \item{positive_ids}{Patient IDs in the positive group.}
#'   \item{omit_ids}{Patient IDs that should be omitted.}
#'   }
ICD2Pheno <- function(pheno_orig,
                      icd10_codes_keep = NULL, icd10_codes_u60 = NULL,
                      icd10_codes_omit = NULL, icd9_codes_keep = NULL,
                      icd9_codes_u60 = NULL, icd9_codes_omit = NULL,
                      controls_only = FALSE) {

  # get id, year of birth, icd10 codes
  icd10 <- pheno_orig %>%
    dplyr::select(eid,
      gender = "31-0.0", birth_year = "34-0.0",
      tidyselect::starts_with("41202"),
      tidyselect::starts_with("41204")
    )
  # get id, year of birth, icd9 codes
  icd9 <- pheno_orig %>%
    dplyr::select(eid,
      gender = "31-0.0", birth_year = "34-0.0",
      tidyselect::starts_with("41203"),
      tidyselect::starts_with("41205")
    )

  ## determine if icd codes of interest appear at least once in icd9/10 datasets

  # icd10 codes to keep
  if (!is.null(icd10_codes_keep)) {
    icd10_keep <- apply(
      icd10[, -(1:3)], 1,
      function(x) {
        return(any(x %in% icd10_codes_keep))
      }
    )
  } else {
    icd10_keep <- F
  }

  # icd10 codes to keep if under 60
  if (!is.null(icd10_codes_u60)) {
    icd10_u60 <- apply(
      icd10[, -(1:3)], 1,
      function(x) {
        return(any(x %in% icd10_codes_u60))
      }
    )
  } else {
    icd10_u60 <- F
  }

  # icd10 codes to omit
  if (!is.null(icd10_codes_omit)) {
    icd10_omit <- apply(
      icd10[, -(1:3)], 1,
      function(x) {
        return(any(x %in% icd10_codes_omit))
      }
    )
  } else {
    icd10_omit <- F
  }

  # icd9 codes to keep
  if (!is.null(icd9_codes_keep)) {
    icd9_keep <- apply(
      icd9[, -(1:3)], 1,
      function(x) {
        return(any(x %in% icd9_codes_keep))
      }
    )
  } else {
    icd9_keep <- F
  }

  # icd9 codes to keep if under 60
  if (!is.null(icd9_codes_u60)) {
    icd9_u60 <- apply(
      icd9[, -(1:3)], 1,
      function(x) {
        return(any(x %in% icd9_codes_u60))
      }
    )
  } else {
    icd9_u60 <- F
  }

  # icd9 codes to omit
  if (!is.null(icd9_codes_omit)) {
    icd9_omit <- apply(
      icd9[, -(1:3)], 1,
      function(x) {
        return(any(x %in% icd9_codes_omit))
      }
    )
  } else {
    icd9_omit <- F
  }

  # get positive ids
  if (!controls_only) {
    positive_ids <- icd10 %>%
      dplyr::select(eid, gender, birth_year) %>%
      dplyr::mutate(age = 2020 - birth_year) %>%
      dplyr::filter((icd10_keep + icd10_u60 * (age < 60) +
        icd9_keep + icd9_u60 * (age < 60)) *
        (1 - icd10_omit) * (1 - icd9_omit) > 0) %>%
      dplyr::pull(eid)
    return(positive_ids)
  } else {
    positive_ids <- icd10 %>%
      dplyr::select(eid, gender, birth_year) %>%
      dplyr::mutate(age = 2020 - birth_year) %>%
      dplyr::filter((icd10_keep + icd10_u60 * (age < 60) +
        icd9_keep + icd9_u60 * (age < 60)) > 0) %>%
      dplyr::pull(eid)
    omit_ids <- icd10 %>%
      dplyr::select(eid, gender, birth_year) %>%
      dplyr::mutate(age = 2020 - birth_year) %>%
      dplyr::filter(icd10_omit | icd9_omit) %>%
      dplyr::pull(eid)
    return(list(positive_ids = positive_ids, omit_ids = omit_ids))
  }
}

#' Load in gender and age demographics
#'
#' @inheritParams loadPheno
#' @returns A matrix with columns "gender" and "age"
loadDemographics <- function(geno_id) {
  DEMO_DATA_PATH <- file.path("..", "data", "ukbb_demographics_data.csv")
  demo_df <- data.table::fread(DEMO_DATA_PATH)
  pheno_df <- data.frame(id = as.numeric(as.character(geno_id))) %>%
    dplyr::left_join(
      y = demo_df %>%
        dplyr::select(id = eid, gender = "31-0.0", birth_year = "34-0.0"),
      by = "id"
    ) %>%
    dplyr::mutate(age = 2020 - birth_year) %>%
    dplyr::select(-birth_year, -id)
  return(as.matrix(pheno_df))
}

#' Load in height and weight from cMRI visit
#'
#' @inheritParams loadPheno
#' @returns A matrix with columns "weight" and "height"
loadLVDemographics <- function(geno_id) {
  DEMO_DATA_PATH <- file.path(
    "..", "data", "table_ventricular_volume_with_indexing.csv"
  )
  demo_df <- data.table::fread(DEMO_DATA_PATH)
  pheno_df <- data.frame(id = as.numeric(as.character(geno_id))) %>%
    dplyr::left_join(
      y = demo_df %>%
        dplyr::select(id, weight = `weight (kg)`, height = `height (cm)`),
      by = "id"
    ) %>%
    dplyr::select(-id)
  return(as.matrix(pheno_df))
}

#' Load in PCs to account for population stratification
#'
#' @inheritParams loadPheno
#' @param npcs Number of PCs
#'
#' @returns A matrix with columns PC1, ..., PC{npcs}.
loadPCs <- function(geno_id, npcs) {
  PCS_DATA_PATH <- file.path("..", "data", "ukbb_pc_data.csv")
  pcs_df <- data.table::fread(PCS_DATA_PATH)
  keep_pcs <- paste0("PC", 1:npcs)
  pcs_out <- data.frame(id = as.numeric(as.character(geno_id))) %>%
    dplyr::left_join(y = pc_df, by = "id") %>%
    dplyr::select(tidyselect::all_of(keep_pcs))
  return(as.matrix(pcs_out))
}

#' Flip SNP encoding (0->2 and 2->0).
#'
#' @description SNP encodings are set so that within a gene, SNPs are (mostly)
#'   positively correlated.
#'
#' @param geno_train Training SNP data.
#' @param snps_df SNP to gene mapping data frame.
#' @param cor_thr Correlation threshold. If the correlation is below this
#'   threshold, the SNP encoding is flipped.
#'
#' @returns A list of 2:
#' \describe{
#' \item{geno_train}{Cleaned training SNP data (same size as input).}
#' \item{flip_snps}{Vector of SNPs that were flipped.}
#' }
flipSNPData <- function(geno_train, snps_df, cor_thr = -0.15) {
  genes <- unique(snps_df$Gene)
  names(genes) <- genes
  flip_snps_all <- purrr::map(
    genes,
    function(gene) {
      x <- geno_train[, snps_df$Gene == gene, drop = FALSE]
      if (ncol(x) == 1) {
        return(NULL)
      } else {
        cor_mat <- cor(x)
        cor_mat[is.na(cor_mat)] <- 1
        flip_snps <- cor_mat[, 1] < cor_thr
        if (!any(flip_snps)) {
          return(NULL)
        } else {
          if (sum(!flip_snps) < sum(flip_snps)) {
            flip_snps <- !flip_snps
          }
          return(colnames(x)[flip_snps])
        }
      }
    }
  ) %>%
    purrr::reduce(c)

  geno_train_flipped <- geno_train
  geno_train_flipped[, flip_snps_all] <- 2 - geno_train[, flip_snps_all]

  return(list(geno_train = geno_train_flipped, flip_snps = flip_snps_all))
}

#' Binarize phenotypes (stratified by gender if provided)
#'
#' @param pheno_train Training phenotype vector.
#' @param gender_train (Optional) Training gender vector if discretization is
#'   meant to be gender-specific.
#' @param pheno_test Test phenotype vector.
#' @param gender_test (Optional) Test gender vector if discretization is meant
#'   to be gender-specific.
#' @param low_q Lower quantile threshold.
#' @param hi_q Upper quantile threshold.
#'
#' @returns A list of 4:
#' \describe{
#' \item{pheno_train}{Discretized training phenotype vector where 1 = high
#'   phenotype value (>hi_q), 0 = low phenotype value (<low_q), and NA o/w.}
#' \item{pheno_test}{Discretized test phenotype vector, similar ot pheno_train.}
#' \item{low_thr}{Lower threshold value.}
#' \item{hi_thr}{Upper threshold value.}
#' }
binarizePhenotype <- function(pheno_train, gender_train = NULL,
                              pheno_test = NULL, gender_test = NULL,
                              low_q = 0.2, hi_q = 0.2) {
  if (is.null(gender_train)) {
    low_thr <- quantile(pheno_train, low_q)
    hi_thr <- quantile(pheno_train, 1 - hi_q)
    pheno_bin <- dplyr::case_when(
      pheno_train <= low_thr ~ 0,
      pheno_train >= hi_thr ~ 1,
      TRUE ~ NA_real_
    )

    if (!is.null(pheno_test)) {
      pheno_va_bin <- dplyr::case_when(
        pheno_test <= low_thr ~ 0,
        pheno_test >= hi_thr ~ 1,
        TRUE ~ NA_real_
      )
    } else {
      pheno_va_bin <- NULL
    }
  } else {
    low_thr0 <- quantile(pheno_train[gender_train == 0], low_q)
    low_thr1 <- quantile(pheno_train[gender_train == 1], low_q)
    low_thr <- c("0" = low_thr0, "1" = low_thr1)

    hi_thr0 <- quantile(pheno_train[gender_train == 0], 1 - hi_q)
    hi_thr1 <- quantile(pheno_train[gender_train == 1], 1 - hi_q)
    hi_thr <- c("0" = hi_thr0, "1" = hi_thr1)

    pheno_bin <- dplyr::case_when(
      (gender_train == 0) & (pheno_train <= low_thr0) ~ 0,
      (gender_train == 0) & (pheno_train >= hi_thr0) ~ 1,
      (gender_train == 1) & (pheno_train <= low_thr1) ~ 0,
      (gender_train == 1) & (pheno_train >= hi_thr1) ~ 1,
      TRUE ~ NA_real_
    )

    if (!is.null(pheno_test) & !is.null(gender_test)) {
      pheno_va_bin <- dplyr::case_when(
        (gender_test == 0) & (pheno_test <= low_thr0) ~ 0,
        (gender_test == 0) & (pheno_test >= hi_thr0) ~ 1,
        (gender_test == 1) & (pheno_test <= low_thr1) ~ 0,
        (gender_test == 1) & (pheno_test >= hi_thr1) ~ 1,
        TRUE ~ NA_real_
      )
    } else {
      pheno_va_bin <- NULL
    }
  }

  return(list(
    pheno_train = pheno_bin,
    pheno_test = pheno_va_bin,
    low_thr = low_thr,
    hi_thr = hi_thr
  ))
}

#' Remove samples with NA phenotype
#'
#' @param geno_train Training genotype data.
#' @param pheno_train Training phenotype vector.
#' @param geno_test Test genotype data.
#' @param pheno_test Test phenotype vector.
#'
#' @returns A list of 4:
#' \describe{
#' \item{geno_train}{Training genotype data without NA phenotype samples.}
#' \item{pheno_train}{Training phenotype data without NA phenotype samples.}
#' \item{geno_test}{Test genotype data without NA phenotype samples.}
#' \item{pheno_test}{Test phenotype data without NA phenotype samples.}
#' }
removeNASamples <- function(geno_train, pheno_train, geno_test, pheno_test) {
  tr_keep_idx <- !is.na(pheno_train)
  geno_train <- geno_train[tr_keep_idx, ]
  pheno_train <- pheno_train[tr_keep_idx]

  ts_keep_idx <- !is.na(pheno_test)
  geno_test <- geno_test[ts_keep_idx, ]
  pheno_test <- pheno_test[ts_keep_idx]

  return(list(
    geno_train = geno_train, pheno_train = pheno_train,
    geno_test = geno_test, pheno_test = pheno_test
  ))
}

#' Load in SNP information
#'
#' @param snps Vector of SNPs that have been named via renameSNPs().
#' @param by Name of column used to merge with annotated SNP data frame.
#' @returns A data.frame with the SNP information and other gene-related
#'   information.
loadSNPInfo <- function(snps, by = "Name",
                        fpath = file.path(
                          "..", "results", "annovar", "snp2gene_df.csv"
                        )) {
  snps_orig <- data.table::fread(fpath)
  snps_df <- data.frame(Name = snps) %>%
    setNames(by) %>%
    dplyr::left_join(y = snps_orig, by = by) %>%
    dplyr::mutate(Gene = stringr::str_replace_all(Gene, "[\\+\\-]", "."))
  snps_df$Gene[is.na(snps_df$Gene)] <- snps_df[[by]][is.na(snps_df$Gene)]
  return(snps_df)
}
