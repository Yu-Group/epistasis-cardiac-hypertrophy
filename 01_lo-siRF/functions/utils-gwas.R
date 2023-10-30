#' Annotate GWAS results
#'
#' @param gwas Data frame of GWAS results.
#' @param gwas_type One of "bolt-lmm" or "plink", indicating the GWAS method.
#' @param av_path File path to snp-to-gene mapping.
#' @param save Logical indicating whether or not tot save the annotation result.
#' @param save_path Where to save the output.
annotGWAS <- function(gwas, gwas_type = c("bolt-lmm", "plink"),
                      av_path = file.path(
                        "..", "results", "annovar", "snp2gene_df.csv"
                      ),
                      save = FALSE, save_path) {
  gwas_type <- match.arg(gwas_type)
  av_out <- data.table::fread(av_path)
  if (identical(gwas_type, "plink")) {
    gwas_annot <- gwas %>%
      dplyr::rename(
        "CHR" = "#CHROM", "BP" = "POS", "SNP" = "ID", "P_GLM" = "P"
      ) %>%
      dplyr::left_join(
        av_out %>% dplyr::select(-Chr, -Pos, -Alt, -Ref),
        by = c("SNP" = "rsID")
      ) %>%
      dplyr::arrange(P_GLM)
  } else {
    gwas_annot <- gwas %>%
      dplyr::left_join(
        av_out %>% dplyr::select(-Chr, -Pos, -Alt, -Ref),
        by = c("SNP" = "rsID")
      ) %>%
      dplyr::arrange(P_BOLT_LMM_INF)
  }

  if (save) {
    write.csv(gwas_annot, save_path)
  }

  return(gwas_annot)
}
