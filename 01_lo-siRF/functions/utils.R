#' Rename SNPs by replacing ':' with 'chr' for chromosome, '_' with 'b' for
#' base pair, and add chromosome number to name
renameSNP <- function(snpName, chr) {
  snpName <- stringr::str_replace_all(snpName, ":", "ch")
  snpName <- stringr::str_replace_all(snpName, "_", "b")
  rs_idx <- startsWith(snpName, "rs")
  snpName[rs_idx] <- paste0(chr[rs_idx], "ch", snpName[rs_idx])
  return(snpName)
}

#' Removes sign information from interactions
cleanInt <- function(x) {
  return(stringr::str_remove_all(x, "[-\\+]"))
}

#' Determine index of features in interaction
getIntIndex <- function(int, varnames, signed) {
  if (!signed) {
    int <- stringr::str_remove_all(int, "[\\+\\-]")
  }
  lapply(stringr::str_split(int, "_"),
    FUN = function(xint) sapply(xint, function(x) which(varnames == x))
  )
}

#' Count number of snps/genes in interaction (i.e., interaction order)
intOrder <- function(int,
                     omit_vars = c(
                       "gender", "age", "height", "weight",
                       paste0("PC", 1:5)
                     )) {
  cleaned_ints <- lapply(
    cleanInt(int),
    function(x) unlist(stringr::str_split(x, "_"))
  )
  int_order <- sapply(cleaned_ints, function(x) sum(!(x %in% omit_vars)))
  return(int_order)
}

#' Annotate interactions
annotInts <- function(int_df, snps_df, snp = TRUE) {
  chr_merge_df <- snps_df %>%
    dplyr::distinct(Chr, Gene)
  int_annot <- int_df %>%
    mutate(
      order = intOrder(int),
      genes = lapply(int, function(x) {
        tmp <- setdiff(
          unlist(stringr::str_split(cleanInt(x), "_")),
          c(
            "gender", "age", "height", "weight",
            paste0("PC", 1:5)
          )
        ) %>%
          sort()
        if (snp) {
          res <- data.frame(Name = tmp[!is.na(tmp)]) %>%
            dplyr::left_join(y = snps_df, by = "Name") %>%
            dplyr::pull(Gene) %>%
            sort()
        } else {
          res <- tmp[!is.na(tmp)]
        }
        return(res)
      }),
      chrs = lapply(int, function(x) {
        chr_merge_df %>%
          dplyr::filter(
            Gene %in% unlist(stringr::str_split(cleanInt(x), "_"))
          ) %>%
          dplyr::pull(Chr) %>%
          unique()
      }),
      chrs_order = lapply(int, function(x) {
        tmp <- chr_merge_df %>%
          dplyr::filter(
            Gene %in% unlist(stringr::str_split(cleanInt(x), "_"))
          ) %>%
          dplyr::pull(Chr) %>%
          unique()
        return(sum(!is.na(tmp)))
      })
    )
  return(int_annot)
}
