#' Calculate the prediction error between the predicted and observed responses
#' across various metrics
#'
#' @param y Vector, matrix, or data.frame of the true response values
#' @param yhat Vector, matrix, or data.frame of the estimated response values
#' @param metric Character vector of prediction error metrics to compute;
#'   elements should be one of "RMSE", "MSE", "R2", "MAE", "Correlation",
#'   "Class", "BalancedClass", "AUC", "PR".
#' @param group (Optional) vector of factors to group prediction errors by
#' @param na_rm = logical; whether or not to remove NAs
#'
#' @returns A data.frame with the following columns:
#' \describe{
#' \item{Group}{(if \code{group} is not NULL); Name of group for grouped
#'   prediction errors.}
#' \item{Meric}{Prediciton error metric.}
#' \item{Value}{Predicttion error value.}
#' }
evalPreds <- function(y, yhat, metric, group = NULL, na_rm = F) {

  # error checking
  isvec <- is.null(dim(y))
  if ((isvec & (length(y) != length(yhat))) |
    (!isvec & any(dim(y) != dim(yhat)))) {
    stop("y and yhat must be the same size.")
  }
  if (!all(metric %in% c(
    "RMSE", "MSE", "R2", "MAE",
    "Correlation", "Class", "BalancedClass",
    "AUC", "PR"
  ))) {
    stop("metric has not been implemented.")
  }
  if (("AUC" %in% metric) | ("PR" %in% metric)) {
    if (length(unique(y)) != 2) {
      stop("y must be binary to evaluate AUC and PR metrics.")
    }
    if ((min(yhat) < 0) | (max(yhat) > 1)) {
      stop("yhat must give the class proportions.")
    }
    ylevels <- levels(as.factor(y))
    Y0 <- ylevels[1]
    Y1 <- ylevels[2]
  }

  # create (long) grouped prediction data frame with groups, y, and yhat
  if (isvec) {
    pred_df <- data.frame(Group = "all", y = y, yhat = yhat)
    if (!is.null(group)) {
      pred_df <- rbind(pred_df, data.frame(Group = group, y = y, yhat = yhat))
    }
    pred_df <- pred_df %>%
      dplyr::group_by(Group)
  } else {
    y_long <- data.frame(Group = "all", y) %>%
      tidyr::gather(key = "column", value = "y", -Group)
    yhat_long <- data.frame(Group = "all", yhat) %>%
      tidyr::gather(key = "column", value = "yhat", -Group)
    if (!is.null(group)) {
      y_long <- rbind(
        y_long,
        data.frame(Group = group, y) %>%
          tidyr::gather(key = "column", value = "y", -Group)
      )
      yhat_long <- rbind(
        yhat_long,
        data.frame(Group = group, yhat) %>%
          tidyr::gather(key = "column", value = "yhat", -Group)
      )
    }
    pred_df <- dplyr::left_join(y_long, yhat_long,
      by = c("Group", "column")
    ) %>%
      dplyr::group_by(Group, column)
  }

  # compute error metrics between y and yhat
  err_out <- NULL
  for (m in metric) {
    if (m == "RMSE") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = sqrt(mean((y - yhat)^2, na.rm = na_rm))
        )
    } else if (m == "MSE") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = mean((y - yhat)^2, na.rm = na_rm)
        )
    } else if (m == "R2") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = 1 - mean((y - yhat)^2, na.rm = na_rm) /
            mean((y - mean(y))^2, na.rm = na_rm)
        )
    } else if (m == "MAE") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = mean(abs(y - yhat), na.rm = na_rm)
        )
    } else if (m == "Correlation") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = cor(y, yhat, use = "pairwise.complete.obs")
        )
    } else if (m == "Class") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = mean(y == yhat, na.rm = na_rm)
        )
    } else if (m == "BalancedClass") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = mean(sapply(
            unique(y),
            function(y0) {
              mean(y[y == y0] == yhat[y == y0],
                na.rm = na_rm
              )
            }
          ), na.rm = na_rm)
        )
    } else if (m == "AUC") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = PRROC::roc.curve(yhat[(y == Y1) & !(is.na(y))],
            yhat[(y == Y0) & !(is.na(y))],
            curve = F
          )$auc
        )
    } else if (m == "PR") {
      err <- pred_df %>%
        dplyr::summarise(
          Metric = m,
          Value = PRROC::pr.curve(yhat[(y == Y1) & !(is.na(y))],
            yhat[(y == Y0) & !(is.na(y))],
            curve = F
          )$auc.integral
        )
    } else {
      stop("Metric has not been implemented.")
    }
    err_out <- rbind(err_out, err)
  }

  # clean up output formatting
  if (!isvec) {
    err_out <- err_out %>%
      tidyr::spread(key = "column", value = "Value") %>%
      dplyr::select(Group, Metric, tidyselect::all_of(colnames(data.frame(y))))
  }
  if (is.null(group)) {
    err_out <- err_out %>%
      dplyr::ungroup() %>%
      dplyr::select(-Group)
  }

  return(err_out)
}

#' Evalute confusion matrix for binary classification problem.
#'
#' @param y Observed binary response vector.
#' @param yhat Predicted binary response vector.
evalConfusion <- function(y, yhat) {
  conf_tab <- table(round(yhat), y)
  if (nrow(conf_tab) != 2) {
    if (!("0" %in% rownames(conf_tab))) {
      conf_tab <- rbind(c(0, 0), conf_tab)
    } else if (!("1" %in% rownames(conf_tab))) {
      conf_tab <- rbind(conf_tab, c(0, 0))
    }
  }
  rownames(conf_tab) <- c("Predicted 0", "Predicted 1")
  colnames(conf_tab) <- c("Observed 0", "Observed 1")
  return(conf_tab)
}

#' Evaluate AUC (for ROC or PR) between observed and predicted responses.
#'
#' @inheritParams evalConfusion
#' @param metric One of "roc" or "pr"
evalAUC <- function(y, yhat, metric = "roc") {
  if (all(yhat == yhat[1])) {
    warning("Predictions are all the same.")
    out <- NULL
  } else {
    if (metric == "roc") {
      out <- PRROC::roc.curve(yhat[y == 1], yhat[y == 0], curve = T)
    } else if (metric == "pr") {
      out <- PRROC::pr.curve(yhat[y == 1], yhat[y == 0], curve = T)
    } else {
      stop("metric is unknown. metric must be one of 'roc' or 'pr'.")
    }
  }
  return(out)
}

#' Evaluate variable importance scores for a given model.
#'
#' @param res_dir Path to results directory.
#' @param method Name of method
#' @param snps_df SNP to gene mapping data frame.
evalVimp <- function(res_dir, method, snps_df) {

  # load in model fit
  load(file.path(res_dir, paste0(method, "_model_fits.Rdata")))

  if (stringr::str_detect(method, "lasso")) {
    if (stringr::str_detect(method, "std")) {
      fit <- lasso_std
    } else {
      fit <- lasso
    }
  } else if (stringr::str_detect(method, "ridge") &
    !stringr::str_detect(method, "kernel")) {
    if (stringr::str_detect(method, "std")) {
      fit <- ridge_std
    } else {
      fit <- ridge
    }
  } else if (stringr::str_detect(method, "svm")) {
    fit <- svmfit
  } else if (stringr::str_detect(method, "xgb")) {
    fit <- xgb
  } else if (stringr::str_detect(method, "irf")) {
    fit <- irf
  } else if (stringr::str_detect(method, "rf")) {
    fit <- rang
  } else if (stringr::str_detect(method, "kernel_ridge")) {
    fit <- kernel_fit
  } else if (stringr::str_detect(method, "logistic")) {
    fit <- log_fit
  }

  # compute variable importance
  if ((stringr::str_detect(method, "ridge")) |
    (stringr::str_detect(method, "lasso"))) {
    imp_df <- as.data.frame(as.matrix(fit$beta)) %>%
      setNames("Importance") %>%
      tibble::rownames_to_column("var") %>%
      dplyr::arrange(dplyr::desc(abs(Importance)))
  } else if (stringr::str_detect(method, "rf")) {
    if (stringr::str_detect(method, "irf")) {
      fit <- fit$rf.list[[length(fit$rf.list)]]
    }
    imp_df <- as.data.frame(fit$variable.importance) %>%
      setNames("Importance") %>%
      tibble::rownames_to_column("var") %>%
      dplyr::mutate(
        var = stringr::str_remove(var, "^X")
      ) %>%
      dplyr::arrange(dplyr::desc(Importance))
  } else {
    imp_df <- NULL
  }

  # annotate SNPs
  if (!is.null(imp_df)) {
    imp_df <- dplyr::left_join(
      x = imp_df, y = snps_df,
      by = c("var" = "Name")
    ) %>%
      dplyr::mutate(Gene = ifelse(is.na(Gene), var, Gene)) %>%
      dplyr::select(-var) %>%
      dplyr::relocate(Importance, .after = last_col())
  }

  return(imp_df)
}
