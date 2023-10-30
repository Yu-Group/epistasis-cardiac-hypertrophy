library(iRF)
library(ranger)
library(glmnet)
library(caret)
library(KRLS)

source("load-functions.R")

#' Fits various models (lasso, ridge, kernel ridge, rf, irf, etc) using
#' training Predixcan or SNP data for a given phenotype
#'
#' @inheritParams runAnalysis
fitModels <- function(pheno_name, pheno_binary = T, n, models = NULL,
                      predx = NULL, include_dems = T, include_npcs = 0,
                      keep_genes = NULL, keep_snps = NULL,
                      n_trees = 500, n_trials = 1,
                      save_path = "results", n_cores = 1) {
  options(stringsAsFactors = FALSE)
  out_file <- "model_fits.Rdata"
  dems <- c("age", "gender", "height", "weight")

  if (!file.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  if (is.null(models)) {
    models <- c(
      "lasso", # lasso
      "lasso_std", # lasso (standardized X)
      "ridge", # ridge
      "ridge_std", # ridge (standardized X)
      "svm", # svm
      "rf", # rf
      "irf"
    ) # irf
  } else if (identical(models, "")) {
    return() # don't run any models
  }

  ##############################################################################
  ############################  Load training  data ############################
  ##############################################################################

  train_data <- loadData(
    pheno_name = pheno_name,
    pheno_binary = pheno_binary,
    n = n, n_test = 0,
    predx = predx,
    include_dems = include_dems,
    include_npcs = include_npcs,
    keep_genes = keep_genes,
    keep_snps = keep_snps
  )
  pheno <- train_data$pheno_train
  geno <- train_data$geno_train
  rm(train_data)

  if (stringr::str_detect(pheno_name, "binary_thr")) {
    # do binarized version
    pheno_binary <- TRUE
  }

  if (pheno_binary) {
    n_pos <- sum(pheno == 1)
    prop_pos <- sum(pheno == 1) / length(pheno)
    n_neg <- sum(pheno == 0)
    prop_neg <- sum(pheno == 0) / length(pheno)
  } else {
    n_pos <- length(pheno)
  }

  if (!stringr::str_detect(save_path, "bootstrap") & (n_trials != 1)) {
    if (n_neg < n_trials * n_pos) {
      stop("(# of neg. samples) < n_trials * (# of pos. samples)")
    }
  }

  ##############################################################################
  #################################### Fit models ##############################
  ##############################################################################

  for (model in models) {
    path_model <- file.path(save_path, paste0(model, "_", out_file))
    if (file.exists(path_model)) { # load in cached model fit
      load(path_model)
    } else { # fit model
      cat(paste0("Running ", model, "...\n"))

      if (stringr::str_detect(model, "baseline")) {
        geno_orig <- geno
        geno <- geno[, colnames(geno) %in% dems]
      }

      if ((stringr::str_detect(model, "lasso") |
        stringr::str_detect(model, "ridge")) &
        !stringr::str_detect(model, "kernel")) {
        #### Fit Lasso/Ridge ####
        if (stringr::str_detect(model, "lasso")) { # fit lasso
          alpha <- 1
        } else { # fit ridge
          alpha <- 0
        }

        if (stringr::str_detect(model, "std")) { # standardize data
          standardize <- TRUE
        } else {
          standardize <- FALSE
        }

        if (pheno_binary) {
          glm_fam <- "binomial"
          if (stringr::str_detect(model, "auc")) { # metric for picking lambda in cv
            type_measure <- "auc"
          } else if (stringr::str_detect(model, "class")) {
            type_measure <- "class"
          } else {
            type_measure <- "deviance"
          }
        } else {
          glm_fam <- "gaussian"
          type_measure <- "deviance"
        }

        # do one hot encoding if categorical variables
        geno_flag <- F
        if (is.data.frame(geno)) {
          geno_orig <- geno
          geno <- model.matrix(~., geno)[, -1]
          geno_flag <- T
        }

        if (n_trials == 1) {
          if (pheno_binary) {
            samp_idx <- c(
              which(pheno == 0)[1:n_neg],
              which(pheno == 1)[1:n_pos]
            )
          } else {
            samp_idx <- 1:n_pos
          }

          cvfit <- glmnet::cv.glmnet(
            x = geno[samp_idx, ], y = pheno[samp_idx],
            alpha = alpha, family = glm_fam,
            standardize = standardize, type.measure = type_measure
          )
          fit <- glmnet::glmnet(
            geno[samp_idx, ], pheno[samp_idx],
            alpha = alpha, family = glm_fam,
            lambda = cvfit$lambda.min, standardize = standardize
          )
          if (alpha == 1) {
            print(paste0("Number of features selected by Lasso: ", fit$df))
          }
          print(paste0(
            "CV Accuracy: ",
            cvfit$cvm[cvfit$lambda == cvfit$lambda.min]
          ))
        } else {
          fit <- list()
          cvfit <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              if (pheno_binary) {
                samp_idx <- c(
                  sample(which(pheno == 0), n_pos, replace = F),
                  which(pheno == 1)
                )
              } else {
                samp_idx <- sample(1:length(pheno), n_pos, replace = T)
              }
            } else {
              if (pheno_binary) {
                samp_idx <- c(
                  which(pheno == 0)[((trial - 1) * n_pos + 1):
                  (n_pos * trial)],
                  which(pheno == 1)[1:n_pos]
                )
              } else {
                samp_idx <- ((trial - 1) * n_pos + 1):(n_pos * trial)
              }
            }
            cvfit[[trial]] <- glmnet::cv.glmnet(
              x = geno[samp_idx, ], y = pheno[samp_idx],
              alpha = alpha, family = glm_fam,
              standardize = standardize, type.measure = type_measure
            )
            fit[[trial]] <- glmnet::glmnet(
              geno[samp_idx, ], pheno[samp_idx],
              alpha = alpha, family = glm_fam,
              lambda = cvfit[[trial]]$lambda.min, standardize = standardize
            )
          }
        }

        # save fits
        if (alpha == 1) {
          if (standardize) {
            lasso_std_cvfit <- cvfit
            lasso_std <- fit
            save(file = path_model, lasso_std_cvfit, lasso_std, compress = T)
          } else {
            lasso_cvfit <- cvfit
            lasso <- fit
            save(file = path_model, lasso_cvfit, lasso, compress = T)
          }
        } else {
          if (standardize) {
            ridge_std_cvfit <- cvfit
            ridge_std <- fit
            save(file = path_model, ridge_std_cvfit, ridge_std, compress = T)
          } else {
            ridge_cvfit <- cvfit
            ridge <- fit
            save(file = path_model, ridge_cvfit, ridge, compress = T)
          }
        }

        if (geno_flag) {
          geno <- geno_orig
        }
      } else if (stringr::str_detect(model, "svm")) {
        #### Fit SVM ####
        if (!pheno_binary) {
          stop("SVM has not been implemented for non-binary responses.")
        }
        svm_cols <- which(apply(as.matrix(geno), 2, sd) != 0)
        svm_trcontrol <- caret::trainControl(
          method = "cv",
          number = 5,
          allowParallel = F,
          verboseIter = T
        )
        # svm_tunegrid <- expand.grid(
        #   C = c(.0001, .001, .01, .1, 1, 10),
        #   sigma = c(1, 10, 100, 250, 500, 750,1000, 5000)
        # )
        svm_tunegrid <- expand.grid(
          C = c(1e-3, 1e-2, 1e-1, 1, 1e2, 1e3),
          sigma = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2)
        )

        # do one hot encoding if categorical variables
        geno_flag <- F
        if (is.data.frame(geno)) {
          geno_orig <- geno
          geno <- model.matrix(~., geno)[, -1]
          geno_flag <- T
        }

        if (n_trials == 1) {
          svm_idx <- c(which(pheno == 0)[1:n_neg], which(pheno == 1)[1:n_pos])
          svmfit <- caret::train(
            as.matrix(geno)[svm_idx, svm_cols], as.factor(pheno)[svm_idx],
            method = "svmRadial",
            trControl = svm_trcontrol,
            tuneGrid = svm_tunegrid,
            type = "C-svc",
            prob.model = T
          )$finalModel
        } else {
          svmfit <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              svm_idx <- c(
                sample(which(pheno == 0), n_pos, replace = F),
                which(pheno == 1)
              )
            } else {
              svm_idx <- c(
                which(pheno == 0)[((trial - 1) * n_pos + 1):
                (n_pos * trial)],
                which(pheno == 1)[1:n_pos]
              )
            }
            svmfit[[trial]] <- caret::train(
              as.matrix(geno)[svm_idx, svm_cols], as.factor(pheno)[svm_idx],
              method = "svmRadial",
              trControl = svm_trcontrol,
              tuneGrid = svm_tunegrid,
              type = "C-svc",
              prob.model = T
            )$finalModel
          }
        }
        save(file = path_model, svmfit, compress = T)

        if (geno_flag) {
          geno <- geno_orig
        }
      } else if (stringr::str_detect(model, "xgb")) {
        #### Fit XGB ####
        if (!pheno_binary) {
          stop("XGB has not been implemented for not binary responses.")
        }
        xgb_trcontrol <- caret::trainControl(
          method = "cv",
          number = 5,
          allowParallel = FALSE,
          verboseIter = FALSE
        )

        xgbGrid <- expand.grid(
          nrounds = c(100, 500),
          max_depth = c(6, 10),
          colsample_bytree = seq(.33, .5),
          eta = c(0.1, 0.3),
          gamma = 0,
          min_child_weight = 1,
          subsample = .6,
          objective = "binary:logistic"
        )

        # do one hot encoding if categorical variables
        geno_flag <- F
        if (is.data.frame(geno)) {
          geno_orig <- geno
          geno <- model.matrix(~., geno)[, -1]
          geno_flag <- T
        }

        if (n_trials == 1) {
          xgb_idx <- c(which(pheno == 0)[1:n_neg], which(pheno == 1)[1:n_pos])
          xgb <- caret::train(
            as.matrix(geno)[xgb_idx, ], as.factor(pheno)[xgb_idx],
            trControl = xgb_trcontrol,
            tuneGrid = xgbGrid,
            method = "xgbTree",
            print_every_n = 50, verbose = 1,
            weight = c(prop_neg * sum(pheno == 1) / sum(pheno == 0), prop_pos)
          )
          print("Best xgb hyperparameters:")
          print(xgb_train$bestTune)
        } else {
          xgb <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              xgb_idx <- c(
                sample(which(pheno == 0), n_pos, replace = F),
                which(pheno == 1)
              )
            } else {
              xgb_idx <- c(
                which(pheno == 0)[((trial - 1) * n_pos + 1):
                (n_pos * trial)],
                which(pheno == 1)[1:n_pos]
              )
            }
            xgb[[trial]] <- caret::train(
              as.matrix(geno)[xgb_idx, ], as.factor(pheno)[xgb_idx],
              trControl = xgb_trcontrol,
              tuneGrid = xgbGrid,
              method = "xgbTree",
              print_every_n = 50, verbose = 1,
              weight = c(prop_neg * sum(pheno == 1) / sum(pheno == 0), prop_pos)
            )
          }
        }
        save(file = path_model, xgb, compress = T)

        if (geno_flag) {
          geno <- geno_orig
        }
      } else if (stringr::str_detect(model, "rf") &
        !stringr::str_detect(model, "irf")) {
        #### Fit Ranger ####
        if (pheno_binary) {
          rang_pheno <- as.factor(pheno)
          rang_importance <- "impurity_corrected"
        } else {
          rang_pheno <- pheno
          rang_importance <- "impurity"
        }

        if (is.data.frame(geno)) {
          rang_df <- data.frame(geno, pheno = rang_pheno)
        } else {
          rang_df <- data.frame(as.matrix(geno), pheno = rang_pheno)
        }

        if (n_trials == 1) {
          # Run full data RF to filter out unimportant features
          if (pheno_binary) {
            rang_idx <- c(
              which(pheno == 0)[1:n_neg],
              which(pheno == 1)[1:n_pos]
            )
          } else {
            rang_idx <- 1:n_pos
          }

          rang <- ranger::ranger(
            data = rang_df[rang_idx, ],
            dependent.variable.name = "pheno",
            classification = pheno_binary,
            importance = rang_importance,
            respect.unordered.factors = "partition",
            # need to manually edit ranger source code to run this
            # uncommented (comment out lines 507-509 in R/ranger.R)
            # sample.fraction = c(prop_neg * sum(pheno == 1) /
            #                       sum(pheno == 0),
            #                     prop_pos),  # to achieve balance
            keep.inbag = T, num.trees = n_trees, num.threads = n_cores
          )
          print(paste0("OOB Error from ranger: ", rang$prediction.error))
          if (pheno_binary) {
            print("Confusion matrix from ranger:")
            print(rang$confusion.matrix)
          }
        } else {
          rang <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              if (pheno_binary) {
                rang_idx <- c(
                  sample(which(pheno == 0), n_pos, replace = F),
                  which(pheno == 1)
                )
              } else {
                rang_idx <- sample(1:length(pheno), n_pos, replace = T)
              }
            } else {
              if (pheno_binary) {
                rang_idx <- c(
                  which(pheno == 0)[((trial - 1) * n_pos + 1):
                  (n_pos * trial)],
                  which(pheno == 1)[1:n_pos]
                )
              } else {
                rang_idx <- ((trial - 1) * n_pos + 1):(n_pos * trial)
              }
            }
            rang[[trial]] <- ranger::ranger(
              data = rang_df[rang_idx, ],
              dependent.variable.name = "pheno",
              classification = pheno_binary,
              respect.unordered.factors = "partition",
              importance = rang_importance,
              num.trees = n_trees,
              num.threads = n_cores
            )
          }
        }
        save(file = path_model, rang, compress = T)
      } else if (stringr::str_detect(model, "irf")) {
        #### Fit iRF on top K features ####
        gnames <- colnames(geno)
        rit_param <- list(
          ntree = 5000, depth = 3, nchild = 5,
          class.id = 1, class.cut = NULL, min.nd = 1
        )
        if (pheno_binary) {
          irf_pheno <- as.factor(pheno)
          irf_importance <- "impurity_corrected"
        } else {
          irf_pheno <- pheno
          irf_importance <- "impurity"
        }

        if (ncol(geno) > 2500) {
          load(file.path(save_path, paste0("rf_", out_file))) # > rang
          print("Select k such that it minimizes oob error")
          KTop <- c(10, 50, 100, 500, 1000, 2500)
          oobError <- numeric(length(KTop))

          if (is.data.frame(geno)) {
            X_irf <- data.frame(geno)
          } else {
            X_irf <- data.frame(as.matrix(geno))
          }

          for (i in 1:length(KTop)) {
            # Select top kTop features from ranger
            kTop <- KTop[i]

            if (n_trials == 1) {
              idkp <- order(rang$variable.importance, decreasing = TRUE)[1:kTop]

              # Run iRF to find specific interactions
              print(paste("Run iRF with kTop =", kTop))
              fit <- iRF::iRF(
                x = X_irf[, idkp],
                y = irf_pheno,
                varnames.grp = gnames[idkp],
                n.iter = 3,
                iter.return = 3,
                select.iter = FALSE,
                type = "ranger",
                respect.unordered.factors = "partition",
                ntree = n_trees,
                # importance = importance,
                # sample_fraction = c(prop_neg * sum(pheno == 1) /
                #                       sum(pheno == 0),
                #                     prop_pos),  # to achieve balance
                n.core = n_cores
              )
              oobError[i] <- fit$rf.list$prediction.error
              print(paste("oob error with kTop =", kTop, "is", oobError[i]))
            } else {
              fit <- list()
              for (trial in 1:n_trials) {
                idkp <- order(rang[[trial]]$variable.importance,
                  decreasing = TRUE
                )[1:kTop]

                if (stringr::str_detect(save_path, "bootstrap")) {
                  if (pheno_binary) {
                    irf_idx <- c(
                      sample(which(pheno == 0), n_pos, replace = F),
                      which(pheno == 1)
                    )
                  } else {
                    irf_idx <- sample(1:length(pheno), n_pos, replace = T)
                  }
                } else {
                  if (pheno_binary) {
                    irf_idx <- c(
                      which(pheno == 0)[((trial - 1) * n_pos + 1):
                      (n_pos * trial)],
                      which(pheno == 1)[1:n_pos]
                    )
                  } else {
                    irf_idx <- ((trial - 1) * n_pos + 1):(n_pos * trial)
                  }
                }

                # Run iRF to find specific interactions
                print(paste("Run iRF with kTop =", kTop))
                fit[[trial]] <- iRF::iRF(
                  x = X_irf[irf_idx, idkp],
                  y = irf_pheno[irf_idx],
                  varnames.grp = gnames[idkp],
                  n.iter = 3,
                  iter.return = 3,
                  select.iter = FALSE,
                  type = "ranger",
                  respect.unordered.factors = "partition",
                  # importance = importance,
                  ntree = n_trees,
                  n.core = n_cores
                )
              }
              oobError[i] <- mean(mapply(
                X = fit,
                function(X) {
                  return(X$rf.list$prediction.error)
                }
              ))
              print(paste("oob error with kTop =", kTop, "is", oobError[i]))
            }
          }

          kTop <- KTop[which.min(oobError)]
          print(paste("Minimal oob error obtained for kTop = ", kTop))
        }

        print("Run iRF with interactions with best kTop")

        if (is.data.frame(geno)) {
          X_irf <- data.frame(geno)
        } else {
          X_irf <- data.frame(as.matrix(geno))
        }

        if (n_trials == 1) {
          if (ncol(geno) > 2500) {
            idkp <- order(rang$variable.importance, decreasing = TRUE)[1:kTop]
            irf_geno <- X_irf[, idkp]
            irf_gnames <- gnames[idkp]
          } else {
            irf_geno <- X_irf
            irf_gnames <- gnames
          }
          print(dim(irf_geno))

          # Run iRF to find specific interactions
          fit <- iRF::iRF(
            x = irf_geno,
            y = irf_pheno,
            varnames.grp = irf_gnames,
            n.iter = 3,
            iter.return = 1:3,
            # int.return = 3,
            select.iter = FALSE,
            n.bootstrap = 50,
            rit.param = rit_param,
            type = "ranger",
            respect.unordered.factors = "partition",
            # importance = importance,
            ntree = n_trees,
            # sample_fraction = c(prop_neg * sum(pheno == 1) /
            #                       sum(pheno == 0),
            #                     prop_pos),  # to achieve balance
            n.core = n_cores
          )
          print(paste0(
            "OOB Error from iRF with best kTop: ",
            fit$rf.list[[3]]$prediction.error
          ))
          if (pheno_binary) {
            print("Confusion matrix from iRF with best kTop:")
            print(fit$rf.list[[3]]$confusion.matrix)
          }

          # Read forest
          print("Read Forest")
          rdForest <- iRF::readForest(
            fit$rf.list[[3]],
            data.frame(as.matrix(irf_geno))
          )
        } else {
          fit <- list()
          rdForest <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              if (pheno_binary) {
                irf_idx <- c(
                  sample(which(pheno == 0), n_pos, replace = F),
                  which(pheno == 1)
                )
              } else {
                irf_idx <- sample(1:length(pheno), n_pos, replace = T)
              }
            } else {
              if (pheno_binary) {
                irf_idx <- c(
                  which(pheno == 0)[((trial - 1) * n_pos + 1):
                  (n_pos * trial)],
                  which(pheno == 1)[1:n_pos]
                )
              } else {
                irf_idx <- ((trial - 1) * n_pos + 1):(n_pos * trial)
              }
            }

            if (ncol(geno) > 2500) {
              idkp <- order(rang[[trial]]$variable.importance,
                decreasing = TRUE
              )[1:kTop]
              irf_geno <- X_irf[irf_idx, idkp]
              irf_gnames <- gnames[idkp]
            } else {
              irf_geno <- X_irf[irf_idx, ]
              irf_gnames <- gnames
            }

            # Run iRF to find specific interactions
            fit[[trial]] <- iRF::iRF(
              x = irf_geno,
              y = irf_pheno[irf_idx],
              varnames.grp = irf_gnames,
              n.iter = 3,
              iter.return = 1:3,
              # int.return = 3,
              select.iter = FALSE,
              n.bootstrap = 50,
              rit.param = rit_param,
              type = "ranger",
              respect.unordered.factors = "partition",
              # importance = importance,
              ntree = n_trees,
              n.core = n_cores
            )

            # Read forest
            print("Read Forest")
            rdForest[[trial]] <- iRF::readForest(
              fit[[trial]]$rf.list[[3]],
              irf_geno
            )
          }
        }

        irf <- fit
        save(file = path_model, irf, rdForest, compress = T)
      } else if (stringr::str_detect(model, "kernel_ridge")) {
        #### Fit Kernel Ridge ####
        if (pheno_binary) {
          stop("Kernel ridge has not been implemented for binary responses.")
        }

        # do one hot encoding if categorical variables
        geno_flag <- F
        if (is.data.frame(geno)) {
          geno_orig <- geno
          geno <- model.matrix(~., geno)[, -1]
          geno_flag <- T
        }

        if (n_trials == 1) {
          samp_idx <- 1:n_pos
          kernel_fit <- KRLS::krls(
            X = geno[samp_idx, ], y = pheno[samp_idx],
            whichkernel = "gaussian", print.level = 0
          )
          print(paste0("Lambda: ", kernel_fit$lambda))
          print(paste0("Sigma: ", kernel_fit$sigma))
          print(paste0("LOOE: ", kernel_fit$Looe))
        } else {
          kernel_fit <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              samp_idx <- sample(1:length(pheno), n_pos, replace = T)
            } else {
              samp_idx <- ((trial - 1) * n_pos + 1):(n_pos * trial)
            }
            kernel_fit[[trial]] <- KRLS::krls(
              X = geno[samp_idx, ],
              y = pheno[samp_idx],
              whichkernel = "gaussian",
              print.level = 0
            )
          }
        }

        # save fits
        save(file = path_model, kernel_fit, compress = T)

        if (geno_flag) {
          geno <- geno_orig
        }
      } else if (stringr::str_detect(model, "logistic")) {
        #### Fit Logistic Regression ####
        if (!pheno_binary) {
          stop("Logistic regression is not applicable for continuous responses.")
        }

        # do one hot encoding if categorical variables
        geno_flag <- F
        if (is.data.frame(geno)) {
          geno_orig <- geno
          geno <- model.matrix(~., geno)[, -1]
          geno_flag <- T
        }

        tr_data <- cbind(as.data.frame(as.matrix(geno)), y = pheno)
        if (n_trials == 1) {
          samp_idx <- c(which(pheno == 0)[1:n_neg], which(pheno == 1)[1:n_pos])
          fit <- glm(y ~ ., data = tr_data[samp_idx, ], family = "binomial")
        } else {
          fit <- list()
          for (trial in 1:n_trials) {
            if (stringr::str_detect(save_path, "bootstrap")) {
              samp_idx <- c(
                sample(which(pheno == 0), n_pos, replace = F),
                which(pheno == 1)
              )
            } else {
              samp_idx <- c(
                which(pheno == 0)[((trial - 1) * n_pos + 1):
                (n_pos * trial)],
                which(pheno == 1)[1:n_pos]
              )
            }
            fit[[trial]] <- glm(y ~ .,
              data = tr_data[samp_idx, ],
              family = "binomial"
            )
          }
        }

        # save fits
        log_fit <- fit
        save(file = path_model, log_fit, compress = T)

        if (geno_flag) {
          geno <- geno_orig
        }
      }

      if (stringr::str_detect(model, "baseline")) {
        geno <- geno_orig
      }
    }
  }

  cat("Completed model fitting.\n")
}
