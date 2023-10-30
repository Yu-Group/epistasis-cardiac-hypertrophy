#' Compute local feature stability scores in a RF.
#' 
#' @param rf_fit Random forest fit. Should be the output of [ranger::ranger()] 
#'   or of class "ranger"
#' @param X Data frame of sample design matrix to evaluate local feature 
#'   stability.
#' @param features Vector of features to evaluate stability. Default (NULL) is
#'   to evaluate the local feature stability for all features.
#' @param feature_groups Data frame of feature to group/superfeature mapping.
#'   Should have the two columns: "feature" and "group".
#' @param first_only Logical; whether or not to only include the first
#'   appearance of the feature in a tree when computing local stability scores.
#'   
#' @returns A data frame of the same size as \code{X} containing the local 
#'   feature stability scores from the fitted random forest.
#' 
localFeatureStabilityRF <- function(rf_fit, X, features = NULL, 
                                    feature_groups = NULL, first_only = FALSE) {
  
  ntrees <- rf_fit$num.trees
  if (is.null(feature_groups)) {
    tree_infos <- purrr::map(1:ntrees, ~ranger::treeInfo(rf_fit, .x))
  } else {
    tree_infos <- purrr::map(
      1:ntrees, 
      ~ranger::treeInfo(rf_fit, .x) %>%
        dplyr::left_join(y = feature_groups, 
                         by = c("splitvarName" = "feature")) %>%
        dplyr::mutate(splitvarName = group) %>%
        dplyr::select(-group)
    )
  }
  forest_paths <- getForestPaths(tree_infos)
  terminal_node_ids <- predict(rf_fit, X, type = "terminalNodes")$predictions
  
  if (is.null(features)) {
    features <- intersect(colnames(X), 
                          rf_fit$forest$independent.variable.names)
    if (!is.null(feature_groups)) {
      features <- feature_groups %>%
        dplyr::filter(feature %in% tidyselect::all_of(features)) %>%
        dplyr::pull(group) %>%
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
            rm_ids <- duplicated(cleanInt(unique_path))
            unique_path <- unique_path[!rm_ids]
          }
          return(unique_path)
        }
      ) %>%
        purrr::reduce(c) %>%
        factor(levels = features) %>%
        table()
    }
  )
  if (is.vector(out)) {
    out <- data.frame(out / ntrees) %>% setNames(features)
  } else {
    out <- as.data.frame(t(out) / ntrees)
  }
  return(out)
}

#' Compute local interacttion stability scores in a RF.
#' 
#' @inheritParams localFeatureStabilityRF
#' @param ints Vector of interactitons to evaluate local stability.
#' 
#' @returns A data frame of the same size as \code{X} containing the local
#'   interaction stability scores from the fitted random forest.
localIntStabilityRF <- function(rf_fit, X, ints, 
                                feature_groups = NULL, first_only = FALSE) {

  ntrees <- rf_fit$num.trees
  if (is.null(feature_groups)) {
    tree_infos <- purrr::map(1:ntrees, ~ranger::treeInfo(rf_fit, .x))
  } else {
    tree_infos <- purrr::map(
      1:ntrees, 
      ~ranger::treeInfo(rf_fit, .x) %>%
        dplyr::left_join(y = feature_groups,
                         by = c("splitvarName" = "feature")) %>%
        dplyr::mutate(splitvarName = group) %>%
        dplyr::select(-group)
    )
  }
  forest_paths <- getForestPaths(tree_infos)
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
            rm_ids <- duplicated(cleanInt(unique_path))
            unique_path <- unique_path[!rm_ids]
          }
          purrr::map_lgl(ints, ~all(.x %in% unique_path))
        }
      ) %>%
        colMeans()
    }
  ) 
  if (is.vector(out)) {
    out <- data.frame(out) %>% setNames(names(ints))
  } else {
    out <- as.data.frame(t(out))
  }
  return(out)
}

#' Run permutation test for the local stability scores using the difference in
#' means test statistic.
#' 
#' @param stab_df Data frame with the local stability scores. Typically the
#'   output of \code{localFeatureStabilityRF()} or \code{localIntStabilityRF()}.
#' @param y Vector of observed responses in the same order as the rows of 
#'   \code{stab_df}.
#' @param feature Character string indicating the feature/interaction to test.
#' @param nperm Number of permutations.
#' 
#' @returns A list of 3:
#' \describe{
#' \item{pval}{Permutation p-value.}
#' \item{T_obs}{Observed difference-in-means test statistic.}
#' \item{perm_dist}{Null permutation distribution.}
#' }
runPermutationTest <- function(stab_df, y, feature, nperm = 1e4) {

  x <- stab_df[, feature]
  
  T_obs <- mean(x[y == 1]) - mean(x[y == 0])
  
  perm_out <- replicate(
    n = nperm,
    expr = {
      y_perm <- sample(y, size = length(y), replace = F)
      mean(x[y_perm == 1]) - mean(x[y_perm == 0])
    }
  )
  
  return(list(pval = mean(abs(perm_out) > abs(T_obs)),
              T_obs = T_obs,
              perm_dist = perm_out))
}

#' Extract all root-to-leaf paths in a forest.
#' 
#' @param tree_infos List of size num.trees with each entry being the output of
#'   [ranger::treeInfo()].
#'   
#' @returns A list of size num.trees with each entry being a list of decision
#'   paths in each tree.
getForestPaths <- function(tree_infos) {
  purrr::map(tree_infos, ~getTreePaths(.x))
}

#' Extract all root-to-leaf paths in a tree.
#' 
#' @param tree_info Output of [ranger::treeInfo()] for a single tree.
#' 
#' @returns A list of decision paths in a single tree.
getTreePaths <- function(tree_info) {

  terminal_node_ids <- tree_info$nodeID[tree_info$terminal]
  inner_tree_info <- tree_info %>%
    dplyr::filter(!terminal)
  tree_paths <- list()
  for (terminal_node_id in terminal_node_ids) {
    node_id <- terminal_node_id
    tree_path <- c()
    while (node_id != 0) {
      if (node_id %in% inner_tree_info$leftChild) {
        idx <- inner_tree_info$leftChild == node_id
        node_id <- inner_tree_info$nodeID[idx]
        tree_path <- c(tree_path, 
                       paste0(inner_tree_info$splitvarName[idx], "-"))
      } else {
        idx <- inner_tree_info$rightChild == node_id
        node_id <- inner_tree_info$nodeID[idx]
        tree_path <- c(tree_path,
                       paste0(inner_tree_info$splitvarName[idx], "+"))
      }
    }
    tree_paths[[as.character(terminal_node_id)]] <- tree_path
  }
  return(tree_paths)
}
