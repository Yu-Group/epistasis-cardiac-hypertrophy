source("local_stability.R", chdir = T)
source("utils.R", chdir = T)

#' Plot distribution of individual SNP frequencies in the fitted RF
#'
#' @param rf_fit Random forest fit; object of class ranger.
#' @param snps_df SNP to gene mapping data fame.
#' @param genes Vector of genes under consideration.
#' @param prop Logical; whether or not to return occurrence proportion.
#' @param return_df Logical; whether or not to return a data frame with the
#'   SNP frequencies.
#'
#' @returns Returns the plot if \code{return_df = FALSE} and returns a list of
#'   two (the plot and the data frame of SNP frequencies) if
#'   \code{return_df = TRUE}.
plotRFSnpDist <- function(rf_fit, snps_df, genes,
                          prop = FALSE, return_df = FALSE) {
  genes <- stringr::str_remove(genes, "[[+/-]]")
  tree_infos_snps <- purrr::map(
    1:rf_fit$num.trees,
    ~ ranger::treeInfo(rf_fit, .x) %>%
      dplyr::left_join(
        y = snps_df %>% dplyr::select(Name, Gene),
        by = c("splitvarName" = "Name")
      )
  )
  tree_infos_gene <- purrr::map(
    tree_infos_snps,
    ~ .x %>%
      dplyr::mutate(splitvarName = Gene) %>%
      dplyr::select(-Gene)
  )
  forest_paths_snps <- getForestPaths(tree_infos_snps)
  forest_paths_gene <- getForestPaths(tree_infos_gene)

  forest_paths <- purrr::map2(
    forest_paths_gene, forest_paths_snps,
    function(fpaths_gene, fpaths_snps) {
      purrr::map2(
        fpaths_gene, fpaths_snps,
        function(fpath_gene, fpath_snp) {
          x <- substr(fpath_snp, start = 1, stop = nchar(fpath_snp) - 1)
          names(x) <- substr(fpath_gene,
            start = 1, stop = nchar(fpath_gene) - 1
          )
          return(x)
        }
      )
    }
  )

  getSNPDist <- function(forest_paths, gene) {
    forest_paths_filt <- purrr::map(
      forest_paths,
      function(fpaths) {
        purrr::map(
          fpaths,
          function(fpath) {
            keep_idx <- names(fpath) %in% gene
            if (any(keep_idx)) {
              return(fpath[keep_idx])
            } else {
              return(NULL)
            }
          }
        )
      }
    )

    gene_counts <- purrr::map_dbl(
      forest_paths_filt,
      ~ mean(!sapply(.x, is.null))
    )
    snp_dist_by_tree <- purrr::map(forest_paths_filt, ~ purrr::reduce(.x, c))
    snp_dist <- purrr::reduce(snp_dist_by_tree, c)

    if (is.null(snp_dist)) {
      snp_freq_by_tree <- NULL
      snp_freq <- NULL
    } else {
      snp_freq_by_tree <- purrr::map2(
        snp_dist_by_tree, 1:length(snp_dist_by_tree),
        function(snps, i) {
          if (is.null(snps)) {
            return(data.frame(
              SNP = snp_dist[1],
              Gene = names(snp_dist)[1],
              n = NA
            ))
          }
          out <- data.frame(SNP = snps, Gene = names(snps)) %>%
            dplyr::group_by(SNP, Gene) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
            dplyr::arrange(desc(n)) %>%
            dplyr::ungroup()
          colnames(out)[colnames(out) == "n"] <- paste0("Tree", i)
          return(out)
        }
      ) %>%
        purrr::reduce(dplyr::full_join, by = c("SNP", "Gene"))
      snp_freq <- cbind(
        snp_freq_by_tree %>% dplyr::select(SNP, Gene),
        n = rowSums(snp_freq_by_tree %>% dplyr::select(-SNP, -Gene), na.rm = T)
      )
    }

    if (return_df) {
      return(list(
        gene_counts = gene_counts,
        snp_freq_by_tree = snp_freq_by_tree,
        snp_freq = snp_freq
      ))
    } else {
      return(snp_freq)
    }
  }

  out <- purrr::map(genes, ~ getSNPDist(forest_paths = forest_paths, gene = .x))

  if (return_df) {
    plt_df <- purrr::map_dfr(out, "snp_freq")
  } else {
    plt_df <- purrr::map_dfr(out, ~.x)
  }

  plt_df <- plt_df %>%
    dplyr::left_join(y = snps_df, by = c("SNP" = "Name", "Gene")) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(Pos) %>%
    dplyr::mutate(
      Name = paste0("SNP = ", rsID, "\nPosition = ", Pos),
      Pos = 1:dplyr::n()
    ) %>%
    # mutate(Pos = Pos - min(Pos)) %>%
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
    dplyr::rename("# Occurrences in RF Paths" = "n") %>%
    dplyr::ungroup()
  if (prop) {
    n_decision_paths <- purrr::map_int(
      1:rf_fit$num.trees,
      ~ ranger::treeInfo(rf_fit, .x) %>%
        dplyr::pull(terminal) %>%
        sum()
    ) %>%
      sum()
    plt_df <- plt_df %>%
      dplyr::mutate(
        `# Occurrences in RF Paths` = `# Occurrences in RF Paths` / n_decision_paths
      )
  }

  plt <- vdocs::plot_bar(
    plt_df,
    x_str = "Pos", y_str = "# Occurrences in RF Paths",
    stat = "identity", size = 0.1
  ) +
    ggplot2::aes(text = Name) +
    ggplot2::facet_wrap(~Gene, scales = "free", ncol = 3) +
    ggplot2::labs(x = "SNV") +
    vthemes::theme_vmodern() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14))

  if (prop) {
    plt <- plt +
      ggplot2::labs(y = "Occurrence Proportion across RF Paths")
  }

  if (return_df) {
    return(list(plot = plt, df = out))
  } else {
    return(plt)
  }
}

#' Plot distribution of individual SNP interaction frequencies in the fitted RF
#'
#' @inheritParams plotRFSnpDist
#' @param ints Vector of interactions under consideration.
#'
#' @returns Returns the plot if \code{return_df = FALSE} and returns a list of
#'   two (the plot and the data frame of SNP interaction frequencies) if
#'   \code{return_df = TRUE}.
plotRFSnpIntDist <- function(rf_fit, snps_df, ints,
                             prop = FALSE, return_df = FALSE) {
  tree_infos_snps <- purrr::map(
    1:rf_fit$num.trees,
    ~ ranger::treeInfo(rf_fit, .x) %>%
      dplyr::left_join(
        y = snps_df %>% dplyr::select(Name, Gene),
        by = c("splitvarName" = "Name")
      )
  )
  tree_infos_gene <- purrr::map(
    tree_infos_snps,
    ~ .x %>%
      dplyr::mutate(splitvarName = Gene) %>%
      dplyr::select(-Gene)
  )
  forest_paths_snps <- getForestPaths(tree_infos_snps)
  forest_paths_gene <- getForestPaths(tree_infos_gene)

  forest_paths <- purrr::map2(
    forest_paths_gene, forest_paths_snps,
    function(fpaths_gene, fpaths_snps) {
      purrr::map2(
        fpaths_gene, fpaths_snps,
        function(fpath_gene, fpath_snp) {
          x <- substr(fpath_snp, start = 1, stop = nchar(fpath_snp) - 1)
          names(x) <- substr(fpath_gene,
            start = 1, stop = nchar(fpath_gene) - 1
          )
          return(x)
        }
      )
    }
  )

  getSNPIntDist <- function(forest_paths, int) {
    int_genes <- stringr::str_remove(strsplit(int, "_")[[1]], "[[+/-]]")
    forest_paths_filt <- purrr::map(
      forest_paths,
      function(fpaths) {
        purrr::map(
          fpaths,
          function(fpath) {
            if (all(int_genes %in% names(fpath))) {
              keep_idx <- names(fpath) %in% int_genes
              return(fpath[keep_idx])
            } else {
              return(NULL)
            }
          }
        )
      }
    )

    int_counts <- purrr::map_dbl(forest_paths_filt, ~ mean(!sapply(.x, is.null)))
    snp_dist_by_tree <- purrr::map(forest_paths_filt, ~ purrr::reduce(.x, c))
    snp_dist <- purrr::reduce(snp_dist_by_tree, c)

    if (is.null(snp_dist)) {
      snp_freq_by_tree <- NULL
      snp_freq <- NULL
      int_dist_by_tree <- NULL
      int_dist <- NULL
      int_freq_by_tree <- NULL
      int_freq <- NULL
    } else {
      snp_freq_by_tree <- purrr::map2(
        snp_dist_by_tree, 1:length(snp_dist_by_tree),
        function(snps, i) {
          if (is.null(snps)) {
            return(data.frame(
              SNP = snp_dist[1],
              Gene = names(snp_dist)[1],
              n = NA
            ))
          }
          out <- data.frame(
            SNP = snps,
            Gene = names(snps)
          ) %>%
            dplyr::group_by(SNP, Gene) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
            dplyr::arrange(dplyr::desc(n)) %>%
            dplyr::ungroup()
          colnames(out)[colnames(out) == "n"] <- paste0("Tree", i)
          return(out)
        }
      ) %>%
        purrr::reduce(dplyr::full_join, by = c("SNP", "Gene"))
      snp_freq <- cbind(
        snp_freq_by_tree %>% dplyr::select(SNP, Gene),
        n = rowSums(snp_freq_by_tree %>% dplyr::select(-SNP, -Gene), na.rm = T)
      )

      int_dist_by_tree <- purrr::map(
        forest_paths_filt,
        function(fpaths) {
          purrr::map(
            fpaths,
            function(fpath) {
              if (is.null(fpath)) {
                return(NULL)
              }
              tapply(fpath, names(fpath), FUN = function(x) x, simplify = F) %>%
                expand.grid() %>%
                tidyr::unite("int") %>%
                dplyr::pull(int)
            }
          ) %>%
            purrr::reduce(c)
        }
      )
      int_dist <- purrr::reduce(int_dist_by_tree, c)

      int_freq_by_tree <- purrr::map2(
        int_dist_by_tree, 1:length(int_dist_by_tree),
        function(ints, i) {
          if (is.null(ints)) {
            return(data.frame(Int = int_dist[1], n = NA))
          }
          out <- data.frame(Int = ints) %>%
            dplyr::group_by(Int) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
            dplyr::arrange(dplyr::desc(n)) %>%
            dplyr::ungroup()
          colnames(out)[colnames(out) == "n"] <- paste0("Tree", i)
          return(out)
        }
      ) %>%
        purrr::reduce(dplyr::full_join, by = c("Int"))
      int_freq <- cbind(
        int_freq_by_tree %>% dplyr::select(Int),
        n = rowSums(int_freq_by_tree %>% dplyr::select(-Int), na.rm = T)
      )
    }

    if (return_df) {
      return(list(
        int_counts = int_counts,
        int_freq_by_tree = int_freq_by_tree,
        int_freq = int_freq
      ))
    } else {
      return(int_freq)
    }
  }

  out <- purrr::map(
    ints, ~ getSNPIntDist(forest_paths = forest_paths, int = .x)
  ) %>%
    setNames(ints)

  if (return_df) {
    plt_df <- purrr::map_dfr(out, "int_freq", .id = "Genes")
  } else {
    plt_df <- dplyr::bind_rows(out, .id = "Genes")
  }

  plt_df <- plt_df %>%
    dplyr::group_by(Genes) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::mutate(
      Name = paste0("Name = ", Int %>%
        stringr::str_remove("[0-9]+ch") %>%
        stringr::str_replace("_[0-9]+ch", "-")),
      Rank = 1:dplyr::n()
    ) %>%
    dplyr::filter(Rank <= 50) %>%
    dplyr::rename("# Occurrences in RF Paths" = "n") %>%
    dplyr::ungroup()
  if (prop) {
    n_decision_paths <- purrr::map_int(
      1:rf_fit$num.trees,
      ~ ranger::treeInfo(rf_fit, .x) %>%
        dplyr::pull(terminal) %>%
        sum()
    ) %>%
      sum()
    plt_df <- plt_df %>%
      dplyr::mutate(
        `# Occurrences in RF Paths` = `# Occurrences in RF Paths` / n_decision_paths
      )
  }

  plt <- vdocs::plot_bar(
    plt_df,
    x_str = "Rank", y_str = "# Occurrences in RF Paths",
    stat = "identity", size = 0.1
  ) +
    ggplot2::aes(text = Name) +
    ggplot2::facet_wrap(~Genes, scales = "free", ncol = 3) +
    ggplot2::labs(x = "Interaction (Ranked by Frequency)") +
    vthemes::theme_vmodern() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14))

  if (prop) {
    plt <- plt +
      ggplot2::labs(y = "Occurrence Proportion across RF Paths")
  }

  if (return_df) {
    return(list(plot = plt, df = out))
  } else {
    return(plt)
  }
}
