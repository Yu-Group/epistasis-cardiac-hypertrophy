#' Removes sign information from interactions
#'
#' @param x A character vector of interactions
#'
#' @return A character vector of interactions without sign information
#'
#' @keywords internal
clean_ints <- function(x) {
  return(stringr::str_remove_all(x, "[-\\+]"))
}


#' Extract all root-to-leaf paths in a forest.
#'
#' @param tree_infos List of size num.trees with each entry being the output of
#'   [ranger::treeInfo()].
#'
#' @returns A list of size num.trees with each entry being a list of decision
#'   paths in each tree.
#'
#' @keywords internal
get_forest_paths <- function(tree_infos) {
  purrr::map(tree_infos, ~get_tree_paths(.x))
}


#' Extract all root-to-leaf paths in a tree.
#'
#' @param tree_info Output of [ranger::treeInfo()] for a single tree.
#'
#' @returns A list of decision paths in a single tree.
#'
#' @keywords internal
get_tree_paths <- function(tree_info) {

  terminal_node_ids <- tree_info$nodeID[tree_info$terminal]
  inner_tree_info <- tree_info |>
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
