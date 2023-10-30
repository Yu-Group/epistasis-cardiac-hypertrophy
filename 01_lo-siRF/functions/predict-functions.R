source("load-functions.R", chdir = TRUE)

#' Make predictions from various models
#'
#' @inheritParams runAnalysis
predictModels <- function(pheno_name, pheno_binary = T, n, n_test,
                          models = NULL, predx = NULL,
                          include_dems = T, include_npcs = 0,
                          keep_genes = NULL, keep_snps = NULL,
                          n_trials = 1, save_path = "results") {
  options(stringsAsFactors = FALSE)
  out_file <- "model_fits.Rdata"
  dems <- c("age", "gender", "height", "weight")

  # load in genotype and phenotype data
  data <- loadData(
    pheno_name = pheno_name,
    pheno_binary = pheno_binary,
    n = n, n_test = n_test,
    predx = predx,
    include_dems = include_dems,
    include_npcs = include_npcs,
    keep_genes = keep_genes,
    keep_snps = keep_snps
  )

  # save data for python models
  write.csv(data$geno_train, file.path(save_path, "Xtrain.csv"), row.names = F)
  write.csv(data$geno_test, file.path(save_path, "Xtest.csv"), row.names = F)
  write.csv(data$pheno_train, file.path(save_path, "ytrain.csv"), row.names = F)
  write.csv(data$pheno_test, file.path(save_path, "ytest.csv"), row.names = F)

  pheno_train <- data$pheno_train
  geno_train <- data$geno_train
  data$geno_train <- NULL # for memory reasons
  pheno <- data$pheno_test
  geno <- data$geno_test
  rm(data)

  if (stringr::str_detect(pheno_name, "binary_thr")) {
    # do binarized version
    pheno_binary <- TRUE
  }

  print(paste("memory after loading all training and test data:", mem_used()))

  save(file = file.path(save_path, "validation_pheno.Rdata"), pheno)

  ##############################################################################
  ##############################  Test Model  ##################################
  ##############################################################################

  ypred_ls <- list()
  for (model in models) {
    if (model == "") {
      next
    }
    cat(paste0("Predicting and evaluating ", model, "...\n"))
    path_model <- file.path(save_path, paste0(model, "_", out_file))
    load(path_model)

    if (stringr::str_detect(model, "baseline")) {
      geno_train_orig <- geno_train
      geno_orig <- geno
      geno_train <- geno_train[, colnames(geno_train) %in% dems]
      geno <- geno[, colnames(geno) %in% dems]
    }

    if ((stringr::str_detect(model, "lasso") |
      stringr::str_detect(model, "ridge")) &
      !stringr::str_detect(model, "kernel")) {
      #### Predict Lasso/Ridge ####
      if (stringr::str_detect(model, "lasso")) {
        if (stringr::str_detect(model, "std")) {
          fit <- lasso_std
        } else {
          fit <- lasso
        }
      } else {
        if (stringr::str_detect(model, "std")) {
          fit <- ridge_std
        } else {
          fit <- ridge
        }
      }

      if (n_trials == 1) {
        ypred <- predict(fit, geno, type = "response")
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = fit,
            function(X) {
              return(predict(X, geno, type = "response"))
            }
          )
        } else {
          ypred <- mapply(
            X = fit, trial = 1:length(fit),
            function(X, trial) {
              if (pheno_binary) {
                samp_idx <- c(
                  which(pheno == 0)[((trial - 1) * n_pos_test + 1):
                  (n_pos_test * trial)],
                  which(pheno == 1)[1:n_pos_test]
                )
              } else {
                samp_idx <- ((trial - 1) * n_pos_test + 1):(n_pos_test * trial)
              }
              return(predict(X, geno[samp_idx, ], type = "response"))
            }
          )
        }
        # ypred <- rowMeans(ypred)
      }
    } else if (stringr::str_detect(model, "svm")) {
      #### Predict SVM ####
      if (!pheno_binary) {
        stop("SVM has not been implemented for non-binary responses.")
      }
      svm_cols <- which(apply(as.matrix(geno_train), 2, sd) != 0)

      if (n_trials == 1) {
        ypred <- predict(svmfit, as.matrix(geno)[, svm_cols],
          type = "probabilities"
        )[, 2]
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = svmfit,
            function(X) {
              return(predict(X, as.matrix(geno)[, svm_cols],
                type = "probabilities"
              )[, 2])
            }
          )
        } else {
          ypred <- mapply(
            X = svmfit, trial = 1:length(svmfit),
            function(X, trial) {
              samp_idx <- c(
                which(pheno == 0)[((trial - 1) * n_pos_test + 1):
                (n_pos_test * trial)],
                which(pheno == 1)[1:n_pos_test]
              )
              return(predict(X, as.matrix(geno)[samp_idx, svm_cols],
                type = "probabilities"
              )[, 2])
            }
          )
        }
        # ypred <- rowMeans(ypred)
      }
    } else if (stringr::str_detect(model, "xgb")) {
      #### Predict XGB ####
      if (!pheno_binary) {
        stop("XGB has not been implemented for non-binary responses.")
      }
      if (n_trials == 1) {
        ypred <- predict(xgb, as.matrix(geno), type = "prob")$`1`
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = xgb,
            function(X) {
              return(predict(X, as.matrix(geno), type = "prob")$`1`)
            }
          )
        } else {
          ypred <- mapply(
            X = xgb, trial = 1:length(xgb),
            function(X, trial) {
              samp_idx <- c(
                which(pheno == 0)[((trial - 1) * n_pos_test + 1):
                (n_pos_test * trial)],
                which(pheno == 1)[1:n_pos_test]
              )
              return(predict(X, as.matrix(geno)[samp_idx, ], type = "prob")$`1`)
            }
          )
        }
        # ypred <- rowMeans(ypred)
      }
    } else if (stringr::str_detect(model, "rf") & !stringr::str_detect(model, "irf")) {
      #### Predict Ranger ####
      if (stringr::str_detect(save_path, "discrete")) {
        rang_df <- geno_df
      } else {
        rang_df <- data.frame(as.matrix(geno))
      }
      if (n_trials == 1) {
        if (pheno_binary) {
          ypred <- predict(rang,
            data = rang_df,
            predict.all = TRUE, num.threads = 1
          )
          ypred <- rowMeans(ypred$predictions - 1) # turn (1, 2)'s to (0, 1)'s
        } else {
          ypred <- predict(rang,
            data = rang_df, predict.all = FALSE,
            num.threads = 1
          )$predictions
        }
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = rang,
            function(X) {
              if (pheno_binary) {
                pred <- predict(X,
                  data = rang_df,
                  predict.all = T, num.threads = 1
                )
                return(rowMeans(pred$predictions - 1))
              } else {
                pred <- predict(X,
                  data = rang_df, predict.all = FALSE,
                  num.threads = 1
                )$predictions
                return(pred)
              }
            }
          )
        } else {
          ypred <- mapply(
            X = rang, trial = 1:length(rang),
            function(X, trial) {
              if (pheno_binary) {
                samp_idx <- c(
                  which(pheno == 0)[((trial - 1) * n_pos_test + 1):
                  (n_pos_test * trial)],
                  which(pheno == 1)[1:n_pos_test]
                )
                pred <- predict(X,
                  data = rang_df[samp_idx, ],
                  predict.all = T, num.threads = 1
                )
                return(rowMeans(pred$predictions - 1))
              } else {
                samp_idx <- ((trial - 1) * n_pos_test + 1):(n_pos_test * trial)
                pred <- predict(X,
                  data = rang_df[samp_idx, ], predict.all = F,
                  num.threads = 1
                )$predictions
                return(pred)
              }
            }
          )
        }
        # ypred <- rowMeans(ypred)
      }
    } else if (stringr::str_detect(model, "irf")) {
      #### Predict iRF ####
      geno_orig <- geno
      geno_train_orig <- geno_train

      if (stringr::str_detect(save_path, "discrete")) {
        geno_all <- geno_df
        geno_all_train <- geno_train_df
      } else {
        geno_all <- geno
        geno_all_train <- geno_train
      }

      cnames <- colnames(data.frame(as.matrix(geno_all)))

      if (n_trials == 1) {
        if (length(irf$rf.list[[3]]$variable.importance) != ncol(geno_all)) {
          ind <- sapply(
            names(irf$rf.list[[3]]$variable.importance),
            function(x) which(cnames == x)
          )

          geno <- geno_all[, ind]
          geno_train <- geno_all_train[, ind]
        }

        print(paste("dim of geno", dim(geno)))
        print(paste("length of pheno", length(pheno)))
        print(paste("mean of pheno", mean(pheno)))

        if (stringr::str_detect(save_path, "discrete")) {
          X_irf <- geno
        } else {
          X_irf <- data.frame(as.matrix(geno))
        }

        if (pheno_binary) {
          ypred <- predict(irf$rf.list[[3]],
            data = X_irf,
            predict.all = T, num.threads = 1
          )
          ypred <- rowMeans(ypred$predictions)
        } else {
          ypred <- predict(irf$rf.list[[3]],
            data = X_irf, predict.all = F,
            num.threads = 1
          )$predictions
        }
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = irf,
            function(X) {
              imp <- X$rf.list[[3]]$variable.importance
              if (length(imp) != ncol(geno_all)) {
                ind <- sapply(names(imp), function(x) {
                  which(cnames == x)
                })
                geno <- geno_all[, ind]
                geno_train <- geno_all_train[, ind]
              }
              if (pheno_binary) {
                pred <- predict(X$rf.list[[3]],
                  data = data.frame(as.matrix(geno)),
                  predict.all = T, num.threads = 1
                )
                return(rowMeans(pred$predictions))
              } else {
                pred <- predict(X$rf.list[[3]],
                  data = data.frame(as.matrix(geno)),
                  predict.all = F, num.threads = 1
                )$predictions
                return(pred)
              }
            }
          )
        } else {
          ypred <- mapply(
            X = irf, trial = 1:length(irf),
            function(X, trial) {
              imp <- X$rf.list[[3]]$variable.importance
              if (length(imp) != ncol(geno_all)) {
                ind <- sapply(names(imp), function(x) {
                  which(cnames == x)
                })
                geno <- geno_all[, ind]
                geno_train <- geno_all_train[, ind]
              }

              if (stringr::str_detect(save_path, "discrete")) {
                X_irf <- geno
              } else {
                X_irf <- data.frame(as.matrix(geno))
              }

              if (pheno_binary) {
                samp_idx <- c(
                  which(pheno == 0)[((trial - 1) * n_pos_test + 1):
                  (n_pos_test * trial)],
                  which(pheno == 1)[1:n_pos_test]
                )
                pred <- predict(X$rf.list[[3]],
                  data = X_irf[samp_idx, ],
                  predict.all = T, num.threads = 1
                )
                return(rowMeans(pred$predictions))
              } else {
                samp_idx <- ((trial - 1) * n_pos_test + 1):(n_pos_test * trial)
                pred <- predict(X$rf.list[[3]],
                  data = X_irf[samp_idx, ],
                  predict.all = F, num.threads = 1
                )$predictions
                return(pred)
              }
            }
          )
        }
        # ypred <- rowMeans(ypred)
      }

      geno <- geno_orig
      geno_train <- geno_train_orig
    } else if (stringr::str_detect(model, "kernel_ridge")) {
      #### Predict Kernel Ridge ####
      if (pheno_binary) {
        stop("Kernel Ridge has not been implemented for binary responses.")
      }

      if (n_trials == 1) {
        ypred <- predict(kernel_fit, newdata = geno)$fit
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = kernel_fit,
            function(X) {
              return(predict(X, newdata = geno)$fit)
            }
          )
        } else {
          ypred <- mapply(
            X = kernel_fit, trial = 1:length(kernel_fit),
            function(X, trial) {
              samp_idx <- ((trial - 1) * n_pos_test + 1):(n_pos_test * trial)
              return(predict(X, geno[samp_idx, ])$fit)
            }
          )
        }
      }
    } else if (stringr::str_detect(model, "logistic")) {
      #### Predict Logistic ####
      if (!pheno_binary) {
        stop("Logistic regression is not applicable for continuous responses.")
      }
      if (n_trials == 1) {
        ypred <- predict(log_fit, as.data.frame(geno), type = "response")
      } else {
        if (stringr::str_detect(save_path, "bootstrap")) {
          ypred <- mapply(
            X = log_fit,
            function(X) {
              return(predict(X, as.data.frame(geno), type = "response"))
            }
          )
        } else {
          ypred <- mapply(
            X = log_fit, trial = 1:length(log_fit),
            function(X, trial) {
              samp_idx <- c(
                which(pheno == 0)[((trial - 1) * n_pos_test + 1):
                (n_pos_test * trial)],
                which(pheno == 1)[1:n_pos_test]
              )
              return(predict(X, as.data.frame(geno[samp_idx, ]),
                type = "response"
              ))
            }
          )
        }
        # ypred <- rowMeans(ypred)
      }
    }

    ypred_ls[[model]] <- ypred

    if (stringr::str_detect(model, "baseline")) {
      geno_train <- geno_train_orig
      geno <- geno_orig
    }
  }

  if (length(ypred_ls) > 0) {
    saveRDS(ypred_ls, file.path(save_path, "ypred.rds"))
  }
}
