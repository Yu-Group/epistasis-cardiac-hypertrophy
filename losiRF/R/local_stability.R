#' Compute the local stability importance scores for features and interactions
#'   in a random forest
#'
#' @param rf_fit Fitted random forest object. Must be of class "ranger" (i.e.,
#'   the output of [ranger::ranger()]).
#' @param X Data frame of sample design matrix, used to evaluate local feature
#'   stability.
#' @param y Vector of observed responses in the same order as the rows of
#'   \code{X}. Must be binary vector with 0s and 1s.
#' @param features Vector of (signed) features, for which to evaluate stability.
#'   Default (\code{NULL}) is to evaluate the local feature stability for all
#'   features in \code{X}.
#' @param ints Vector of (signed) interactions to evaluate local stability.
#' @param feature_groups Data frame of feature-to-group/superfeature mapping.
#'   Should have the two columns: "feature" and "group". Default (\code{NULL})
#'   performs no grouping so that each feature is its own group.
#' @param first_only Logical; whether or not to include only the first
#'   appearance of the feature in a tree when computing local stability scores.
#'   Default (\code{FALSE}) includes all occurrences of a feature in the
#'   local stability score computation.
#' @param nperm Number of permutations used in the permutation test.
#'
#' @returns A list with the following components:
#' \describe{
#'   \item{feature_stability_df}{Data frame of size (\code{nrow(X)}) x
#'     (\code{length(features)}) containing the local feature stability scores
#'     from the fitted random forest.}
#'   \item{feature_stability_pvals}{Tibble of p-values from the permutation test
#'     for each feature in \code{feature_stability_df}.}
#'   \item{int_stability_df}{Data frame of size (\code{nrow(X)}) x
#'     (\code{length(ints)}) containing the local interaction stability scores
#'     from the fitted random forest.}
#'   \item{int_stability_pvals}{Tibble of p-values from the permutation test
#'     for each interaction in \code{int_stability_df}.}
#' }
#'
#' @export
local_rf_stability_importance <- function(rf_fit, X, y,
                                          features = NULL,
                                          ints = NULL,
                                          feature_groups = NULL,
                                          first_only = FALSE,
                                          nperm = 1e4) {
  if (!inherits(rf_fit, "ranger")) {
    stop("rf_fit must be of class 'ranger'.")
  }
  if (!is.null(feature_groups)) {
    if (!isTRUE(all(c("feature", "group") %in% colnames(feature_groups)))) {
      stop("feature_groups must have columns 'feature' and 'group'.")
    }
  }
  if (!is.null(features) && 
      !isTRUE(all(remove_signs(features) %in% colnames(X)))) {
    stop("(unsigned) features must be a subset of colnames(X).")
  }

  # compute local stability importance for features
  feature_stability_df <- local_rf_feature_stability(
    rf_fit = rf_fit, X = X, features = features,
    feature_groups = feature_groups, first_only = first_only
  )
  feature_stability_pvals <- purrr::map(
    colnames(feature_stability_df),
    function(feature) {
      perm_out <- twosample_permutation_test(
        local_stability_df = feature_stability_df,
        y = y, feature = feature, nperm = nperm
      )
      perm_out$perm_dist <- list(perm_out$perm_dist)
      return(tibble::as_tibble_row(perm_out))
    }
  ) |>
    setNames(colnames(feature_stability_df)) |>
    dplyr::bind_rows(.id = "feature")

  # compute local stability importance for interactions
  if (!is.null(ints)) {
    int_stability_df <- local_rf_int_stability(
      rf_fit = rf_fit, X = X, ints = ints,
      feature_groups = feature_groups, first_only = first_only
    )
    int_stability_pvals <- purrr::map(
      colnames(int_stability_df),
      function(int) {
        perm_out <- twosample_permutation_test(
          local_stability_df = int_stability_df,
          y = y, feature = int, nperm = nperm
        )
        perm_out$perm_dist <- list(perm_out$perm_dist)
        return(tibble::as_tibble_row(perm_out))
      }
    ) |>
      setNames(colnames(int_stability_df)) |>
      dplyr::bind_rows(.id = "int")
  } else {
    int_stability_df <- NULL
    int_stability_pvals <- NULL
  }

  out <- list(
    feature_stability_df = feature_stability_df,
    feature_stability_pvals = feature_stability_pvals,
    int_stability_df = int_stability_df,
    int_stability_pvals = int_stability_pvals
  )
  return(out)
}


#' Get most stable features according to average local stability importance
#' 
#' @inheritParams local_rf_stability_importance
#' @param k Number of top features to return. If both \code{k} and 
#'   \code{min_stability} are \code{NULL}, a data frame with all features along 
#'   with their average local stability importance is returned.
#' @param min_stability Minimum threshold for the average local stability 
#'   importance. Any feature with a lower average local stability importance
#'   will be filtered out. Only used if \code{k} is \code{NULL}. If both 
#'   \code{k} and \code{min_stability} are \code{NULL}, a data frame with all 
#'   features along with their average local stability importance is returned.
#'   
#' @returns If \code{k} and \code{min_stability} are \code{NULL}, returns data 
#'   frame with all features along with their average local stability 
#'   importance, which can be used to rank the features manually. Otherwise,
#'   returns vector of the most stable features.
#'   
#' @export
get_top_stable_features <- function(rf_fit, X, y, 
                                    k = NULL, min_stability = NULL,
                                    features = NULL,
                                    feature_groups = NULL,
                                    first_only = FALSE) {
  if (!inherits(rf_fit, "ranger")) {
    stop("rf_fit must be of class 'ranger'.")
  }
  if (!is.null(feature_groups)) {
    if (!isTRUE(all(c("feature", "group") %in% colnames(feature_groups)))) {
      stop("feature_groups must have columns 'feature' and 'group'.")
    }
  }
  if (!is.null(features) && 
      !isTRUE(all(remove_signs(features) %in% colnames(X)))) {
    stop("(unsigned) features must be a subset of colnames(X).")
  }
  
  # compute local stability importance for features
  feature_stability_df <- local_rf_feature_stability(
    rf_fit = rf_fit, X = X, features = features,
    feature_groups = feature_groups, first_only = first_only
  ) |>
    dplyr::summarise(
      dplyr::across(
        tidyselect::everything(),
        ~ mean(.x)
      )
    ) |>
    tidyr::pivot_longer(
      cols = tidyselect::everything(), 
      names_to = "feature", 
      values_to = "stability"
    ) |>
    dplyr::arrange(
      dplyr::desc(stability)
    )
  
  if (is.null(k) && is.null(min_stability)) {
    return(feature_stability_df)
  } else if (!is.null(k)) {
    return(feature_stability_df$feature[1:k])
  } else {
    return(
      feature_stability_df |> 
        dplyr::filter(stability > min_stability) |>
        dplyr::pull(feature)
    )
  }
}


#' Compute local feature stability scores in a random forest.
#'
#' @inheritParams local_rf_stability_importance
#'
#' @returns A data frame of the same size as \code{X} containing the local
#'   feature stability scores from the fitted random forest.
#'
#' @export
local_rf_feature_stability <- function(rf_fit, X,
                                       features = NULL,
                                       feature_groups = NULL,
                                       first_only = FALSE) {

  if (!inherits(rf_fit, "ranger")) {
    stop("rf_fit must be of class 'ranger'.")
  }
  if (!is.null(feature_groups)) {
    if (!isTRUE(all(c("feature", "group") %in% colnames(feature_groups)))) {
      stop("feature_groups must have columns 'feature' and 'group'.")
    }
  }
  if (!is.null(features) && 
      !isTRUE(all(remove_signs(features) %in% colnames(X)))) {
    stop("(unsigned) features must be a subset of colnames(X).")
  }

  ntrees <- rf_fit$num.trees
  if (is.null(feature_groups)) {
    tree_infos <- purrr::map(1:ntrees, ~ranger::treeInfo(rf_fit, .x))
  } else {
    tree_infos <- purrr::map(
      1:ntrees,
      ~ranger::treeInfo(rf_fit, .x) |>
        dplyr::left_join(
          y = feature_groups,
          by = c("splitvarName" = "feature")
        ) |>
        dplyr::mutate(splitvarName = group) |>
        dplyr::select(-group)
    )
  }
  forest_paths <- get_forest_paths(tree_infos)
  terminal_node_ids <- predict(rf_fit, X, type = "terminalNodes")$predictions

  if (is.null(features)) {
    features <- intersect(
      colnames(X),
      rf_fit$forest$independent.variable.names
    )
    if (!is.null(feature_groups)) {
      features <- feature_groups |>
        dplyr::filter(
          feature %in% tidyselect::all_of(features)
        ) |>
        dplyr::pull(group) |>
        unique()
    }
    features <- c(paste0(features, c("+")), paste0(features, "-"))
  }

  out <- apply(
    terminal_node_ids, 1,
    function(x_terminal_nodes) {
      purrr::map2(
        forest_paths, x_terminal_nodes,
        function(x, y) {
          unique_path <- rev(unique(x[[as.character(y)]]))
          if (first_only) {
            rm_ids <- duplicated(remove_signs(unique_path))
            unique_path <- unique_path[!rm_ids]
          }
          return(unique_path)
        }
      ) |>
        purrr::reduce(c) |>
        factor(levels = features) |>
        table()
    }
  )
  if (is.vector(out)) {
    out <- data.frame(out / ntrees) |> setNames(features)
  } else {
    out <- as.data.frame(t(out) / ntrees)
  }
  return(out)
}


#' Compute local interacttion stability scores in a random forest.
#'
#' @inheritParams local_rf_stability_importance
#'
#' @returns A data frame of the same size as \code{X} containing the local
#'   interaction stability scores from the fitted random forest.
#'
#' @export
local_rf_int_stability <- function(rf_fit, X, ints,
                                   feature_groups = NULL,
                                   first_only = FALSE) {

  if (!inherits(rf_fit, "ranger")) {
    stop("rf_fit must be of class 'ranger'.")
  }
  if (!is.null(feature_groups)) {
    if (!isTRUE(all(c("feature", "group") %in% colnames(feature_groups)))) {
      stop("feature_groups must have columns 'feature' and 'group'.")
    }
  }

  ntrees <- rf_fit$num.trees
  if (is.null(feature_groups)) {
    tree_infos <- purrr::map(1:ntrees, ~ranger::treeInfo(rf_fit, .x))
  } else {
    tree_infos <- purrr::map(
      1:ntrees,
      ~ranger::treeInfo(rf_fit, .x) |>
        dplyr::left_join(
          y = feature_groups,
          by = c("splitvarName" = "feature")
        ) |>
        dplyr::mutate(splitvarName = group) |>
        dplyr::select(-group)
    )
  }
  forest_paths <- get_forest_paths(tree_infos)
  terminal_node_ids <- predict(rf_fit, X, type = "terminalNodes")$predictions

  if (is.null(ints)) {
    stop("Must provide ints argument.")
  } else {
    ints_name <- ints
    ints <- stringr::str_split(ints, "_")
    names(ints) <- ints_name
  }

  out <- apply(
    terminal_node_ids, 1,
    function(x_terminal_nodes) {
      purrr::map2_dfr(
        forest_paths, x_terminal_nodes,
        function(x, y) {
          unique_path <- rev(unique(x[[as.character(y)]]))
          if (first_only) {
            rm_ids <- duplicated(remove_signs(unique_path))
            unique_path <- unique_path[!rm_ids]
          }
          purrr::map_lgl(ints, ~all(.x %in% unique_path))
        }
      ) |>
        colMeans()
    }
  )
  if (is.vector(out)) {
    out <- data.frame(out) |> setNames(names(ints))
  } else {
    out <- as.data.frame(t(out))
  }
  return(out)
}


#' Run a two-sample permutation test with the difference-in-means test statistic.
#'
#' @inheritParams local_rf_stability_importance
#' @param local_stability_df Data frame with the local stability scores.
#'   Generally the output of \code{local_rf_feature_stability()} or
#'   \code{local_rf_int_stability()}.
#' @param y Vector of observed responses in the same order as the rows of
#'   \code{local_stability_df}. Must be binary vector with 0s and 1s.
#' @param feature Character string indicating the feature/interaction to test.
#'
#' @returns A list of 3:
#' \describe{
#' \item{pval}{Two-sided permutation p-value.}
#' \item{T_obs}{Observed difference-in-means test statistic.}
#' \item{perm_dist}{Null permutation distribution.}
#' }
#'
#' @export
twosample_permutation_test <- function(local_stability_df, y, feature,
                                       nperm = 1e4) {

  if (nrow(local_stability_df) != length(y)) {
    stop("local_stability_df and y must have the same number of rows.")
  }
  if (!isTRUE(all(y %in% c(0, 1)))) {
    stop("y must be a binary vector.")
  }
  if (!feature %in% colnames(local_stability_df)) {
    stop("feature must be a column in local_stability_df.")
  }
  if (nperm < 1) {
    stop("nperm must be a positive integer.")
  }

  x <- local_stability_df[, feature]
  T_obs <- mean(x[y == 1]) - mean(x[y == 0])
  perm_out <- replicate(
    n = nperm,
    expr = {
      y_perm <- sample(y, size = length(y), replace = FALSE)
      mean(x[y_perm == 1]) - mean(x[y_perm == 0])
    }
  )

  out <- list(
    pval = mean(abs(perm_out) > abs(T_obs)),
    T_obs = T_obs,
    perm_dist = perm_out
  )
  return(out)
}
