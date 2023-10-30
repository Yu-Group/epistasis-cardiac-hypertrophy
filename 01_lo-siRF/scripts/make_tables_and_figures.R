library(magrittr)

DATA_DIR <- file.path("..", "data")
FUNCTIONS_DIR <- file.path("..", "functions")
RESULTS_DIR <- file.path("..", "results")
FIGS_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
TABLES_SUPP_DIR <- file.path(RESULTS_DIR, "tables_supplementary")

source(file.path(FUNCTIONS_DIR, "load-functions.R"), chdir = TRUE)
source(file.path(FUNCTIONS_DIR, "eval-functions.R"), chdir = TRUE)
source(file.path(FUNCTIONS_DIR, "plot-functions.R"), chdir = TRUE)

save_fig <- function(plot, filename, width, height, devices = c("pdf", "png")) {
  for (device in devices) {
    ggplot2::ggsave(
      plot = plot,
      filename = file.path(FIGS_DIR, sprintf("%s.%s", filename, device)),
      device = device, width = width, height = height
    )
  }
}

#### Get and save summary of lo-siRF results ####
phenotype <- "iLVM"
iter <- 3
PVAL_THR <- 0.1

gene_pvals_df <- NULL
int_pvals_df <- NULL
for (thr in c(0.15, 0.2, 0.25)) {
  fpath <- file.path(
    RESULTS_DIR,
    paste0(phenotype, "_binary_thr", thr, "_gwas_filtered_nsnps1000")
  )
  # > out, yhat_tr2, yhat_va, irf_gene_out, snps.df
  load(file.path(fpath, "irf_interaction_fit.Rdata"))

  # local interaction stability results
  # > int_stab_tr2_ls, int_stab_va_ls, test_ints_ls, int_perm_out_ls
  load(file.path(fpath, "local_interaction_stability_results.Rdata"))
  int_out <- list(
    stab_tr2_ls = int_stab_tr2_ls,
    stab_va_ls = int_stab_va_ls,
    test_ints_ls = test_ints_ls,
    perm_out_ls = int_perm_out_ls
  )

  pval_int <- purrr::map(
    1:length(int_out$perm_out_ls),
    ~ purrr::map_dfr(
      int_out$perm_out_ls[[.x]],
      function(x) {
        data.frame(pval = x$pval) %>%
          setNames(paste0("P-value (Iter", .x, ")"))
      },
      .id = "int"
    )
  ) %>%
    purrr::reduce(dplyr::full_join, by = "int")

  # local feature stability results
  # > stab_tr2_ls, stab_va_ls, test_genes_ls, perm_out_ls
  load(file.path(fpath, "local_feature_stability_results.Rdata"))
  stab_out <- list(
    stab_tr2_ls = stab_tr2_ls,
    stab_va_ls = stab_va_ls,
    test_genes_ls = test_genes_ls,
    perm_out_ls = perm_out_ls
  )

  ft_pvals <- purrr::map_dfr(
    stab_out$perm_out_ls[[iter]],
    ~ data.frame(pval = .x$pval),
    .id = "feature"
  ) %>%
    dplyr::arrange(pval)

  gene_pvals_df <- dplyr::bind_rows(
    gene_pvals_df,
    data.frame(
      thr = thr, phenotype = phenotype,
      gene = ft_pvals$feature,
      pval = ft_pvals$pval
    )
  )

  int_pvals <- purrr::map_dfr(
    int_out$perm_out_ls[[iter]],
    ~ data.frame(pval = .x$pval),
    .id = "int"
  ) %>%
    dplyr::arrange(pval)

  int_pvals_df <- dplyr::bind_rows(
    int_pvals_df,
    data.frame(
      thr = thr, phenotype = phenotype,
      int = int_pvals$int,
      pval = int_pvals$pval
    )
  )
}

# get genes/interactions that appear in all three binarization runs and with p-val < 0.1
keep_ints <- int_pvals_df %>%
  dplyr::filter(!is.na(pval), pval <= !!PVAL_THR) %>%
  # dplyr::filter(!is.na(pval)) %>%
  dplyr::count(int) %>%
  dplyr::filter(n >= 3) %>%
  dplyr::pull(int)
keep_genes <- gene_pvals_df %>%
  dplyr::filter(!is.na(pval), pval <= !!PVAL_THR) %>%
  # dplyr::filter(!is.na(pval)) %>%
  dplyr::count(gene) %>%
  dplyr::filter(n >= 3) %>%
  dplyr::pull(gene)
c(keep_ints, keep_genes)

results <- int_pvals_df %>%
  dplyr::rename(gene = int) %>%
  dplyr::bind_rows(gene_pvals_df) %>%
  dplyr::filter(gene %in% c(keep_ints, keep_genes)) %>%
  dplyr::mutate(thr = as.factor(thr))

# save results
saveRDS(results, file.path(TABLES_DIR, "losiRF_results.rds"))

#### Main lo-siRF results table ####
results <- readRDS(file.path(TABLES_DIR, "losiRF_results.rds"))

mean_results <- results %>%
  dplyr::group_by(phenotype, gene) %>%
  dplyr::summarise(pval = round(mean(pval), 4)) %>%
  dplyr::mutate(thr = "mean")

table_wide <- dplyr::bind_rows(results, mean_results) %>%
  tidyr::pivot_wider(names_from = thr, values_from = pval, id_cols = -phenotype) %>%
  dplyr::arrange(mean) %>%
  dplyr::rename(
    `Loci / Interaction` = gene,
    `Mean p-value` = mean,
    `0.20` = "0.2"
  ) %>%
  # remove signs and take min p-value for each signed feature group
  dplyr::mutate(
    `Loci / Interaction` = stringr::str_remove_all(`Loci / Interaction`, "[[+/-]]")
  ) %>%
  dplyr::distinct(`Loci / Interaction`, .keep_all = TRUE)

## add snp information
snps <- c(
  "rs7591091", "rs62024491", "rs569550", "rs621679",
  "rs9401921", "rs66733621", "rs12465459",
  "rs6999852", "rs71537873", "rs28657002", "rs55686521", "rs11633294"
)
snp_names <- c(
  "2chrs7591091", "15chrs62024491", "11chrs569550", "11chrs621679",
  "6chrs9401921", "2chrs66733621", "2chrs12465459",
  "8chrs6999852", "8chrs71537873", "15chrs28657002", "15chrs55686521", "15chrs11633294"
)

snps_df_orig <- data.table::fread(
  file.path(RESULTS_DIR, "annovar", "snp2gene_df.csv")
)

snps_df <- snps_df_orig %>%
  tibble::as_tibble() %>%
  dplyr::filter(rsID %in% snps) %>%
  dplyr::rename(
    "SNP" = "rsID",
    "Position" = "Pos",
    "Function" = "Gene Function"
  ) %>%
  dplyr::select(-Name, -`Exonic Function`, -Ref, -Alt)

gwas_plink <- data.table::fread(
  file.path(RESULTS_DIR, "gwas_plink", "iLVM_norm.stats_annot")
)

top_snps <- c(
  CCDC141 = "rs7591091",
  IGF1R = "rs62024491",
  LSP1 = "rs569550",
  `MIR588;RSPO3` = "rs9401921",
  TTN = "rs66733621"
)
top_snps_df <- tibble::tibble(
  `Loci / Interaction` = names(top_snps),
  SNP = top_snps
) %>%
  dplyr::left_join(
    snps_df,
    by = "SNP"
  ) %>%
  dplyr::select(`Loci / Interaction`, SNP, Chr, Position, Function) %>%
  dplyr::left_join(
    gwas_plink %>%
      dplyr::select(
        SNP,
        `GWAS p-value` = P_GLM, `GWAS Beta` = BETA, `GWAS SE` = SE
      ),
    by = "SNP"
  ) %>%
  dplyr::left_join(
    table_wide %>%
      dplyr::select(`Loci / Interaction`, `Mean p-value`),
    by = "Loci / Interaction"
  ) %>%
  dplyr::arrange(`Mean p-value`) %>%
  dplyr::mutate(
    `Mean p-value` = ifelse(
      `Mean p-value` < 1e-4, "$< 10<sup>-4</sup>$",
      ifelse(`Mean p-value` < 1e-3, "< 10<sup>-3</sup>",
        sprintf("%s", formatC(`Mean p-value`, digits = 3, format = "f"))
      )
    ),
    # convert to hg38 positions
    Position = c(
      "98733068",
      "126925592",
      "178799323",
      "1865838",
      "178889467"
    )
  ) %>%
  dplyr::rename(
    "Loci" = "Loci / Interaction",
    "Top SNV" = "SNP",
    "Mean lo-siRF p-value" = "Mean p-value"
  )

int_snps <- c(
  CCDC141 = "rs7591091",
  TTN = "rs66733621",
  CCDC141 = "rs7591091",
  IGF1R = "rs62024491",
  CCDC141 = "rs7591091",
  `LOC157273;TNKS` = "rs6999852"
)
top_int_snps_df <- tibble::tibble(
  Gene = names(int_snps),
  SNP = int_snps
) %>%
  dplyr::left_join(
    snps_df,
    by = c("SNP", "Gene")
  ) %>%
  dplyr::mutate(
    Interaction = c(
      "CCDC141_TTN", "CCDC141_TTN",
      "CCDC141_IGF1R", "CCDC141_IGF1R",
      "CCDC141_LOC157273;TNKS", "CCDC141_LOC157273;TNKS"
    )
  ) %>%
  dplyr::left_join(
    table_wide %>%
      dplyr::select(`Loci / Interaction`, `Mean p-value`),
    by = c("Interaction" = "Loci / Interaction")
  ) %>%
  dplyr::arrange(`Mean p-value`) %>%
  dplyr::mutate(
    Interaction = stringr::str_replace_all(Interaction, "_", "&#8211;"),
    `Mean p-value` = ifelse(
      `Mean p-value` < 1e-4, "$< 10<sup>-4</sup>$",
      ifelse(`Mean p-value` < 1e-3, "< 10<sup>-3</sup>",
        sprintf("%s", formatC(`Mean p-value`, digits = 3, format = "f"))
      )
    )
  ) %>%
  dplyr::group_by(Interaction) %>%
  dplyr::summarise(
    `Top SNV` = paste(sprintf("%s", `SNP`), collapse = "<br>"),
    `Mean lo-siRF p-value` = unique(`Mean p-value`),
    dplyr::across(Chr:Function, ~ paste(.x, collapse = "<br>"))
  ) %>%
  dplyr::mutate(
    # convert to hg38 positions
    Position = c(
      "178889467<br>98727212",
      "178889467<br>178799323",
      "178889467<br>9478458"
    )
  )

write.csv(
  top_snps_df, file.path(TABLES_SUPP_DIR, "losiRF_loci_results_table.csv")
)
write.csv(
  top_int_snps_df, file.path(TABLES_SUPP_DIR, "losiRF_int_results_table.csv")
)

top_snps_df %>%
  vthemes::pretty_DT(
    rownames = FALSE,
    options = list(
      dom = "t", ordering = FALSE,
      columnDefs = list(
        list(width = "220px", targets = 0),
        list(className = "dt-center", targets = "_all")
      ),
      initComplete = DT::JS(
        "function(settings, json) {",
        "$('body').css({'font-family': 'sans-serif'});",
        "}"
      )
    ),
    escape = FALSE, na_disp = "--"
  )

top_int_snps_df %>%
  vthemes::pretty_DT(
    rownames = FALSE,
    options = list(
      dom = "t", ordering = FALSE,
      columnDefs = list(
        list(width = "220px", targets = 0),
        list(className = "dt-center", targets = "_all")
      ),
      initComplete = DT::JS(
        "function(settings, json) {",
        "$('body').css({'font-family': 'sans-serif'});",
        "}"
      )
    ),
    escape = FALSE, na_disp = "--"
  )

#### Distribution of LVM/LVMi ####
# get phenotype data for white british unrelated individuals only
geno_ids <- readRDS(file.path(DATA_DIR, "geno_ids.rds"))
pheno_data_orig <- data.table::fread(
  file.path(DATA_DIR, "table_ventricular_volume_with_indexing.csv")
) %>%
  dplyr::filter(id %in% !!geno_ids)

lvm_plot <- vdocs::plot_histogram(pheno_data, x_str = "LVM (g)", bins = 30)
lvmi_plot <- vdocs::plot_histogram(
  pheno_data,
  x_str = "iLVM (g/m2)", bins = 30
) +
  ggplot2::labs(x = expression(bold(paste(LVMi ~ (g / m^2)))))
lvm_v_lvmi_plot <- vdocs::plot_point(
  pheno_data,
  x_str = "LVM (g)", y_str = "iLVM (g/m2)", size = 0.1, alpha = 0.2
) +
  ggplot2::labs(y = expression(bold(paste(LVMi ~ (g / m^2)))))

save_fig(
  lvm_plot,
  filename = "lvm_distribution", width = 3.5, height = 3.25
)
save_fig(
  lvmi_plot,
  filename = "lvmi_distribution", width = 3.5, height = 3.25
)
save_fig(
  lvm_v_lvmi_plot,
  filename = "lvm_v_lvmi_scatter", width = 3.5, height = 3.25
)

#### Distribution of Local Stability Importance Scores for Interactions ####
keep_ints <- c(
  "CCDC141-_IGF1R-", "CCDC141-_TTN-", "CCDC141-_LOC157273;TNKS-"
)
keep_genes <- c(
  "IGF1R-", "MIR588;RSPO3+", "TTN-", "TTN+", "MIR588;RSPO3-", "LSP1-", "CCDC141-"
)

phenotype <- "iLVM"
iter <- 3
type <- "Using all splits"
int_lstab_df <- NULL
int_lstab_ls <- list()
ft_lstab_df <- NULL
ft_lstab_ls <- list()
for (thr in c(0.15, 0.2, 0.25)) {
  fpath <- file.path(
    RESULTS_DIR,
    paste0(phenotype, "_binary_thr", thr, "_gwas_filtered_nsnps1000")
  )
  # > out, yhat_tr2, yhat_va, irf_gene_out, snps.df
  load(file.path(fpath, "irf_interaction_fit.Rdata"))

  # local interaction stability results
  # > int_stab_tr2_ls, int_stab_va_ls, test_ints_ls, int_perm_out_ls
  load(file.path(fpath, "local_interaction_stability_results.Rdata"))
  int_out <- list(
    stab_tr2_ls = int_stab_tr2_ls,
    stab_va_ls = int_stab_va_ls,
    test_ints_ls = test_ints_ls,
    perm_out_ls = int_perm_out_ls
  )

  int_all_out <- list(`Using all splits` = int_out)

  pval_int <- purrr::map(
    1:length(int_out$perm_out_ls),
    ~ purrr::map_dfr(
      int_out$perm_out_ls[[.x]],
      function(x) {
        data.frame(pval = x$pval) %>%
          setNames(paste0("P-value (Iter", .x, ")"))
      },
      .id = "int"
    )
  ) %>%
    purrr::reduce(dplyr::full_join, by = "int")

  if (any(pval_int$int %in% keep_ints)) {
    # distribution plots
    int_pvals <- purrr::map_dfr(
      int_all_out[[type]]$perm_out_ls[[iter]],
      ~ data.frame(pval = .x$pval),
      .id = "int"
    ) %>%
      dplyr::filter(int %in% keep_ints) %>%
      dplyr::arrange(pval) %>%
      dplyr::mutate(int_pval = paste0(int, " (p-val = ", round(pval, 3), ")"))
    plt_df <- int_all_out[[type]]$stab_va_ls[[iter]] %>%
      dplyr::bind_cols(y = as.factor(out$pheno_test)) %>%
      tidyr::gather(key = "int", value = "Stability", -y) %>%
      dplyr::filter(int %in% keep_ints) %>%
      dplyr::left_join(y = int_pvals, by = "int")
    int_lstab_df <- dplyr::bind_rows(
      int_lstab_df,
      plt_df %>% dplyr::mutate(thr = !!thr)
    )
    int_lstab_ls[[as.character(thr)]] <- int_all_out[[type]]$stab_va_ls[[iter]] %>%
      dplyr::select(tidyselect::all_of(keep_ints)) %>%
      dplyr::bind_cols(y = as.factor(out$pheno_test))
  }

  # local feature stability results
  # > stab_tr2_ls, stab_va_ls, test_genes_ls, perm_out_ls
  load(file.path(fpath, "local_feature_stability_results.Rdata"))
  stab_out <- list(
    stab_tr2_ls = stab_tr2_ls,
    stab_va_ls = stab_va_ls,
    test_genes_ls = test_genes_ls,
    perm_out_ls = perm_out_ls
  )

  stab_all_out <- list(`Using all splits` = stab_out)

  stab_df <- stab_all_out[[type]]$stab_va_ls[[iter]] %>%
    dplyr::select(all_of(stab_all_out[[type]]$test_genes_ls[[iter]]))
  ft_pvals <- purrr::map_dfr(
    stab_all_out[[type]]$perm_out_ls[[iter]],
    ~ data.frame(pval = .x$pval),
    .id = "feature"
  ) %>%
    dplyr::filter(feature %in% keep_genes) %>%
    dplyr::arrange(pval) %>%
    dplyr::mutate(ft_pval = paste0(feature, " (p-val = ", round(pval, 3), ")"))
  plt_df <- stab_df %>%
    dplyr::bind_cols(y = as.factor(out$pheno_test)) %>%
    tidyr::gather(key = "feature", value = "Stability", -y) %>%
    dplyr::filter(feature %in% keep_genes) %>%
    dplyr::left_join(y = ft_pvals, by = "feature")
  ft_lstab_df <- dplyr::bind_rows(
    ft_lstab_df,
    plt_df %>% dplyr::mutate(thr = !!thr)
  )
  ft_lstab_ls[[as.character(thr)]] <- stab_df %>%
    dplyr::select(tidyselect::all_of(keep_genes)) %>%
    dplyr::bind_cols(y = as.factor(out$pheno_test))
}

# feature stability kernel densities
plt_df <- ft_lstab_df %>%
  dplyr::mutate(
    feature = factor(feature, levels = keep_genes),
    thr = paste0("Threshold = ", thr),
    y = factor(forcats::fct_recode(y, "Low" = "0", "High" = "1"),
      levels = c("High", "Low")
    ),
    pval_str = dplyr::case_when(
      pval < 1e-4 ~ "p < 10^-4",
      pval < 1e-3 ~ "p < 10^-3",
      TRUE ~ sprintf("p~`=`~%s", formatC(pval, digits = 3, format = "f"))
    )
  )

plt <- ggplot2::ggplot(plt_df) +
  ggplot2::aes(x = Stability, fill = y) +
  ggplot2::geom_density(alpha = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(x = 0.65, y = 7.25, label = pval_str),
    data = plt_df %>% dplyr::distinct(feature, thr, .keep_all = TRUE),
    size = 5, parse = TRUE
  ) +
  ggplot2::facet_grid(thr ~ feature) +
  ggplot2::labs(
    x = "Local Stability Importance Score", y = "Density", fill = "LVMi"
  ) +
  ggplot2::scale_fill_manual(values = c("#2B4550", "#F8AC47")) +
  vthemes::theme_vmodern(size_preset = "large")
save_fig(
  plt,
  filename = "local_stability_score_dist_loci", width = 16, height = 8
)

# interaction stability kernel densities
plt_df <- int_lstab_df %>%
  dplyr::mutate(
    int = factor(int, levels = keep_ints),
    thr = paste0("Threshold = ", thr),
    y = factor(forcats::fct_recode(y, "Low" = "0", "High" = "1"),
      levels = c("High", "Low")
    ),
    pval_str = dplyr::case_when(
      pval < 1e-4 ~ "p < 10^-4",
      pval < 1e-3 ~ "p < 10^-3",
      TRUE ~ sprintf("p~`=`~%s", formatC(pval, digits = 3, format = "f"))
    )
  )
plt <- ggplot2::ggplot(plt_df) +
  ggplot2::aes(x = Stability, fill = y) +
  ggplot2::geom_density(alpha = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(x = 0.6, y = 8, label = pval_str),
    data = plt_df %>% dplyr::distinct(int, thr, .keep_all = TRUE),
    size = 5, parse = TRUE
  ) +
  ggplot2::facet_grid(thr ~ int) +
  ggplot2::labs(
    x = "Local Stability Importance Score", y = "Density", fill = "LVMi"
  ) +
  ggplot2::scale_fill_manual(values = c("#2B4550", "#F8AC47")) +
  vthemes::theme_vmodern(size_preset = "large")
save_fig(
  plt,
  filename = "local_stability_score_dist_int", width = 12.2, height = 8
)

#### SNP Frequencies ####
keep_ints <- c(
  "CCDC141-_IGF1R-", "CCDC141-_TTN-", "CCDC141-_LOC157273;TNKS-"
)
keep_genes <- c(
  "IGF1R", "MIR588;RSPO3", "TTN", "LSP1", "CCDC141", "LOC157273;TNKS"
)

phenotype <- "iLVM"

snp_dist_ls <- list()
int_snp_dist_ls <- list()
for (thr in c(0.15, 0.2, 0.25)) {
  print(thr)
  fpath <- file.path(
    RESULTS_DIR,
    paste0(phenotype, "_binary_thr", thr, "_gwas_filtered_nsnps1000")
  )
  # > out, yhat_tr2, yhat_va, irf_gene_out, snps.df
  load(file.path(fpath, "irf_interaction_fit.Rdata"))

  snp_dist_ls[[as.character(thr)]] <- plotRFSnpDist(
    irf_gene_out$rf.list[[3]], snps.df, keep_genes,
    return_df = TRUE
  )
  int_snp_dist_ls[[as.character(thr)]] <- plotRFSnpIntDist(
    irf_gene_out$rf.list[[3]], snps.df, keep_ints,
    return_df = TRUE
  )
}

# get total number of decision paths per tree
n_decision_paths <- list()
for (thr in c(0.15, 0.2, 0.25)) {
  fpath <- file.path(
    RESULTS_DIR,
    paste0(phenotype, "_binary_thr", thr, "_gwas_filtered_nsnps1000")
  )
  # > out, yhat_tr2, yhat_va, irf_gene_out, snps.df
  load(file.path(fpath, "irf_interaction_fit.Rdata"))

  rf_fit <- irf_gene_out$rf.list[[3]]
  n_decision_paths[[as.character(thr)]] <- purrr::map_int(
    1:rf_fit$num.trees,
    ~ ranger::treeInfo(rf_fit, .x) %>%
      dplyr::pull(terminal) %>%
      sum()
  ) %>%
    sum()
}

snps_df <- data.table::fread(
  file.path(RESULTS_DIR, "annovar", "snp2gene_df.csv")
)

snp_dist_df <- purrr::map2_dfr(
  snp_dist_ls, n_decision_paths,
  function(x, total_paths) {
    purrr::map_dfr(x$df, "snp_freq") %>%
      dplyr::left_join(y = snps_df, by = c("SNP" = "Name", "Gene")) %>%
      dplyr::group_by(Gene) %>%
      dplyr::arrange(Pos) %>%
      dplyr::mutate(
        Name = paste0("Name = ", rsID, "\nPos = ", Pos),
        Order = 1:dplyr::n(),
        Pos = forcats::fct_inorder(as.character(Pos)),
        Prop = n / total_paths
      ) %>%
      dplyr::mutate(
        dplyr::across(tidyselect:::where(is.character), as.factor)
      ) %>%
      dplyr::ungroup()
  },
  .id = "Threshold"
)

int_snp_dist_df <- purrr::map2_dfr(
  int_snp_dist_ls, n_decision_paths,
  function(x, total_paths) {
    purrr::map_dfr(x$df, "int_freq", .id = "Gene") %>%
      dplyr::group_by(Gene) %>%
      dplyr::arrange(-n) %>%
      dplyr::mutate(
        Name = paste0("Name = ", Int),
        Order = 1:dplyr::n(),
        Int = forcats::fct_inorder(as.character(Int)),
        Prop = n / total_paths
      ) %>%
      dplyr::mutate(
        dplyr::across(tidyselect:::where(is.character), as.factor)
      ) %>%
      dplyr::ungroup()
  },
  .id = "Threshold"
)

snp_dist_df_agg <- snp_dist_df %>%
  tidyr::pivot_wider(
    id_cols = c(Gene, SNP, rsID, Name, Pos),
    names_from = "Threshold",
    values_from = "Prop"
  ) %>%
  dplyr::mutate(
    Prop = (`0.15` + `0.2` + `0.25`) / 3
  ) %>%
  dplyr::filter(!is.na(Prop)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::arrange(Pos) %>%
  dplyr::mutate(
    Order = 1:dplyr::n(),
    Pos = forcats::fct_inorder(as.character(Pos))
  ) %>%
  dplyr::ungroup()

int_snp_dist_df_agg <- int_snp_dist_df %>%
  tidyr::pivot_wider(
    id_cols = c(Gene, Int),
    names_from = "Threshold",
    values_from = "Prop"
  ) %>%
  dplyr::mutate(
    Prop = (`0.15` + `0.2` + `0.25`) / 3
  ) %>%
  dplyr::filter(!is.na(Prop)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::arrange(-Prop) %>%
  dplyr::mutate(
    Order = 1:dplyr::n(),
    Int = forcats::fct_inorder(as.character(Int))
  ) %>%
  dplyr::filter(Order <= 50)

snp_plt <- snp_dist_df_agg %>%
  vdocs::plot_bar(x_str = "Order", y_str = "Prop", stat = "identity") +
  ggplot2::aes(text = Name) +
  ggplot2::facet_wrap(~Gene, scales = "free") +
  ggplot2::labs(
    x = "SNV (Ordered by Position)",
    y = "Occurrence Proportion\nacross RF Paths"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "white"),
    panel.background = ggplot2::element_rect(fill = "white")
  )

int_plt <- int_snp_dist_df_agg %>%
  vdocs::plot_bar(x_str = "Int", y_str = "Prop", stat = "identity") +
  ggplot2::aes(text = Int) +
  ggplot2::facet_wrap(~Gene, scales = "free") +
  ggplot2::labs(x = "SNV-SNV Pair", y = "Occurrence Proportion\nacross RF Paths") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "white"),
    panel.background = ggplot2::element_rect(fill = "white")
  )

save_fig(
  snp_plt,
  filename = "snv_frequency_of_occurrence", width = 12, height = 5.5
)
save_fig(
  int_plt,
  filename = "snv_int_frequency_of_occurrence", width = 12, height = 3
)

write.csv(
  snp_dist_df_agg,
  file = file.path(TABLES_DIR, "snv_frequency_of_occurrence.csv"),
  row.names = FALSE
)
write.csv(
  int_snp_dist_df_agg,
  file = file.path(TABLES_DIR, "snv_int_frequency_of_occurrence.csv"),
  row.names = FALSE
)

#### Patient Characteristics Table ####

# see get_patient_characteristics.R file

#### Binarization Thresholds ####
phenotype <- "iLVM"

# pulled numbers from .out files when running 03_prediction_models.R for binarized LVMi phenotype
thresholds_df <- data.frame(
  threshold = rep(c("0.15", "0.20", "0.25"), each = 4),
  gender = rep(c("Male", "Female"), each = 2, times = 3),
  direction = rep(c("Low", "High"), times = 6),
  value = c(
    43.80485, 58.48090, 35.14884, 46.10459,
    45.07532, 56.78818, 36.02134, 44.87951,
    46.02158, 55.41779, 36.75897, 43.83382
  )
)
thresholds_df_wide <- thresholds_df %>%
  tidyr::pivot_wider(
    id_cols = threshold,
    names_from = c(gender, direction), values_from = value
  ) %>%
  dplyr::rename(
    "Binarization Threshold" = "threshold",
    "Male - Low LVMi Threshold" = "Male_Low",
    "Male - High LVMi Threshold" = "Male_High",
    "Female - Low LVMi Threshold" = "Female_Low",
    "Female - High LVMi Threshold" = "Female_High",
  )
write.csv(
  thresholds_df_wide,
  file.path(TABLES_SUPP_DIR, "binarization_thresholds.csv"),
  row.names = FALSE
)

vthemes::pretty_kable(thresholds_df_wide,
  format = "latex",
  col.names = c(
    "Binarization\\nThreshold", "Low LVMi\\nThreshold", "High LVMi\\nThreshold",
    "Low LVMi\\nThreshold", "High LVMi\\nThreshold"
  )
) %>%
  kableExtra::add_header_above(c(" " = 1, "Male" = 2, "Female" = 2)) %>%
  kableExtra::kable_styling(full_width = F)

vthemes::pretty_kable(thresholds_df_wide,
  format = "html",
  col.names = c(
    "Binarization<br>Threshold", "Low LVMi<br>Threshold", "High LVMi<br>Threshold",
    "Low LVMi<br>Threshold", "High LVMi<br>Threshold"
  )
) %>%
  kableExtra::add_header_above(c(" " = 1, "Male" = 2, "Female" = 2)) %>%
  kableExtra::kable_styling(full_width = F)

#### Prediction Accuracies ####
phenotype <- "iLVM"
phenotype_name <- "LVMi"
methods <- c(
  "Lasso" = "lasso", "Ridge" = "ridge", "RF" = "rf", "siRF" = "irf", "SVM" = "svm"
)

eval_ls <- list()
for (thr in c(0.15, 0.2, 0.25)) {
  fpath <- file.path(RESULTS_DIR, paste0(
    phenotype, "_binary_thr", thr,
    "_gwas_filtered_nsnps1000"
  ))
  # read in prediction results
  load(file.path(fpath, "validation_pheno.Rdata")) # > pheno
  preds_data <- list(
    yhat = c(
      readRDS(file.path(fpath, "ypred.rds")),
      data.table::fread(file.path(fpath, "ypred.csv"))
    ),
    y = pheno
  )

  preds_data$yhat <- preds_data$yhat[methods]
  names(preds_data$yhat) <- names(methods)

  # prediction accuracy table
  metrics <- c("Class", "AUC", "PR")
  eval_out <- purrr::map_dfr(
    preds_data$yhat,
    ~ dplyr::bind_rows(
      evalPreds(
        y = preds_data$y, yhat = c(.x),
        metric = c("AUC", "PR")
      ),
      evalPreds(
        y = preds_data$y,
        yhat = round(c(.x)),
        metric = c("BalancedClass", "Class")
      )
    ),
    .id = "Method"
  ) %>%
    tidyr::spread(key = "Metric", value = "Value")
  eval_ls[[as.character(thr)]] <- eval_out
}

eval_df <- dplyr::bind_rows(eval_ls, .id = "Threshold") %>%
  dplyr::rename("Accuracy" = "Class", "AUROC" = "AUC", "AUPRC" = "PR") %>%
  dplyr::mutate(
    Threshold = paste0("Binarization Threshold = ", Threshold),
    Method = factor(Method, levels = c("siRF", "RF", "Lasso", "Ridge", "SVM"))
  ) %>%
  dplyr::select(Threshold, Method, Accuracy, AUROC, AUPRC) %>%
  dplyr::group_by(Threshold)

tab_csv <- eval_df %>%
  dplyr::mutate(
    dplyr::across(
      c(Accuracy, AUROC, AUPRC),
      ~ formatC(.x, digits = 3, format = "g", flag = "#")
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(Accuracy, AUROC, AUPRC), names_to = "Metric", values_to = "Value"
  ) %>%
  tidyr::pivot_wider(
    id_cols = Method, names_from = c(Threshold, Metric), values_from = Value
  ) %>%
  dplyr::arrange(Method)
write.csv(tab_csv, file.path(TABLES_SUPP_DIR, "prediction_accuracy.csv"))

tab_dt <- eval_df %>%
  dplyr::mutate(dplyr::across(
    c(Accuracy, AUROC, AUPRC),
    ~ ifelse(
      .x == max(.x),
      kableExtra::cell_spec(
        formatC(.x, digits = 3, format = "g", flag = "#"),
        bold = TRUE
      ),
      formatC(.x, digits = 3, format = "g", flag = "#")
    )
  )) %>%
  tidyr::pivot_longer(
    cols = c(Accuracy, AUROC, AUPRC), names_to = "Metric", values_to = "Value"
  ) %>%
  tidyr::pivot_wider(
    id_cols = Method, names_from = c(Threshold, Metric), values_from = Value
  ) %>%
  dplyr::arrange(Method) %>%
  vthemes::pretty_kable(
    format = "html",
    col.names = c("Method", rep(c("Accuracy", "AUROC", "AUPRC"), times = 3))
  ) %>%
  kableExtra::add_header_above(c(
    " " = 1,
    "Binarization Threshold = 0.15" = 3,
    "Binarization Threshold = 0.20" = 3,
    "Binarization Threshold = 0.25" = 3
  ))




#### Summary of siRF Metrics ####
phenotype <- "iLVM"
iter <- 3

irf_out_ls <- list()
for (thr in c(0.15, 0.2, 0.25)) {
  fpath <- file.path(
    RESULTS_DIR,
    paste0(phenotype, "_binary_thr", thr, "_gwas_filtered_nsnps1000")
  )
  # > out, yhat_tr2, yhat_va, irf_gene_out, snps.df
  load(file.path(fpath, "irf_interaction_fit.Rdata"))
  irf_out_ls[[as.character(thr)]] <- irf_gene_out
}

keep_ints <- c("CCDC141-_TTN-", "CCDC141-_IGF1R-", "CCDC141-_LOC157273;TNKS-")

irf_out_df <- purrr::map_dfr(
  irf_out_ls,
  ~ .x$interaction[[3]] %>%
    dplyr::filter(int %in% keep_ints),
  .id = "thr"
) %>%
  dplyr::arrange(int)

tab <- irf_out_df %>%
  dplyr::select(
    Threshold = thr, Interaction = int, Prevalence = prevalence,
    Precision = precision, `Class Difference in Prevalence` = cpe,
    `Stability of Class Difference in Prevalence` = sta.cpe,
    `Independence of Feature Selection` = fsd,
    `Stability of Independence of Feature Selection` = sta.fsd,
    `Increase in Precision` = mip,
    `Stability of Increase in Precision` = sta.mip,
    Stability = stability
  )
write.csv(tab, file.path(TABLES_SUPP_DIR, "irf_summary_metrics.csv"))

tab_dt <- tab %>%
  dplyr::mutate(
    Threshold = paste0("Binarization Threshold = ", Threshold),
    Interaction = stringr::str_replace_all(Interaction, "_", "&#8211;") %>%
      stringr::str_replace_all("-", "<sup>&#8211;</sup>")
  ) %>%
  dplyr::group_by(Threshold) %>%
  vthemes::pretty_DT(
    digits = 2, sigfig = TRUE, rownames = FALSE, na_disp = "--",
    extensions = "RowGroup",
    options = list(
      pageLength = nrow(irf_out_df),
      dom = "t",
      ordering = FALSE,
      rowGroup = list(dataSrc = 1),
      columnDefs = list(
        list(className = "dt-center", targets = "_all"),
        list(visible = FALSE, targets = 1),
        list(width = "160px", targets = 0)
      )
    )
  )
tab_dt

#### lo-siRF p-values Table ####
results <- readRDS(file.path(TABLES_DIR, "losiRF_results.rds"))

mean_results <- results %>%
  dplyr::group_by(phenotype, gene) %>%
  dplyr::summarise(pval = round(mean(pval), 4)) %>%
  dplyr::mutate(thr = "mean")

table_wide_signed <- dplyr::bind_rows(results, mean_results) %>%
  tidyr::pivot_wider(
    names_from = thr, values_from = pval, id_cols = -phenotype
  ) %>%
  dplyr::arrange(mean) %>%
  dplyr::rename(
    `Loci / Interaction` = "gene",
    `Mean p-value` = mean,
    `0.20` = "0.2"
  )

tab <- table_wide_signed %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::all_of(c("0.15", "0.20", "0.25", "Mean p-value")),
      ~ ifelse(.x < 1e-4, "< 10^-4",
        ifelse(.x < 1e-3, "< 10^-3",
          sprintf("%s", formatC(.x, digits = 3, format = "f"))
        )
      )
    )
  )
write.csv(
  tab,
  file.path(TABLES_SUPP_DIR, "losiRF_signed_interactions_table.csv")
)

pval_format <- formattable::formatter(
  .tag = "span",
  style = function(x) {
    if (is.character(x)) {
      x <- stringr::str_remove_all(x, "^<") %>%
        stringr::str_remove_all("$") %>%
        stringr::str_trim()
      x <- stringr::str_remove_all(x, "^<|\\$|\\{|\\}") %>%
        stringr::str_remove_all("</sup>") %>%
        stringr::str_trim() %>%
        stringr::str_replace("0<sup>", "e") %>%
        as.numeric()
    }
    formattable::style(
      display = "block",
      padding = "0 4px",
      `border-radius` = "4px",
      `background-color` = formattable::csscolor(
        formattable::gradient(c(
          -log10(c(0, 0.1) + 0.001),
          as.numeric(-log10(x + 0.001))
        ),
        min.color = "#DEF7E9",
        max.color = "#44B876"
        )
      )[-c(1, 2)]
    )
  }
)

gene_format <- formattable::formatter(
  .tag = "span",
  style = x ~ formattable::style(
    `font-style` = "italic",
    color = ifelse(stringr::str_detect(x, "_"), "#236CE1", "black")
  )
)

tab_kable <- table_wide_signed %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::all_of(c("0.15", "0.20", "0.25", "Mean p-value")),
      ~ ifelse(.x < 1e-4, "< 10<sup>-4</sup>",
        ifelse(.x < 1e-3, "< 10<sup>-3</sup>",
          sprintf("%s", formatC(.x, digits = 3, format = "f"))
        )
      )
    )
  ) %>%
  dplyr::mutate(
    `Loci / Interaction` = gene_format(`Loci / Interaction`),
    dplyr::across(
      tidyselect::all_of(c("0.15", "0.20", "0.25", "Mean p-value")),
      ~ pval_format(.x)
    )
  ) %>%
  dplyr::mutate(
    `Loci / Interaction` = `Loci / Interaction` %>%
      stringr::str_replace_all("-_", "<sup>&#8211;</sup>&#8211;") %>%
      stringr::str_replace_all("- ", "<sup>&#8211;</sup>") %>%
      stringr::str_replace_all("\\+", "<sup>+</sup>")
  )

knitr::kable(
  tab_kable,
  format = "html", align = "rcccc", escape = FALSE
) %>%
  kableExtra::kable_styling(full_width = FALSE) %>%
  kableExtra::add_header_above(
    c(" " = 1, "Binarization Threshold" = 3, " " = 1)
  ) %>%
  kableExtra::column_spec(
    5,
    border_left = "1px solid #ddd", bold = TRUE, background = "#fafafa"
  )

#### List of GWAS-filtered SNPs ####
phenotype <- "iLVM"
nsnps <- 1000
bolt_dir <- file.path(RESULTS_DIR, "gwas_bolt_lmm")
plink_dir <- file.path(RESULTS_DIR, "gwas_plink")
keep_snps <- purrr::map(
  c(bolt_dir, plink_dir),
  function(fdir) {
    fpath <- file.path(fdir, paste0(phenotype, "_norm.stats_annot"))
    snps <- data.table::fread(fpath) %>%
      dplyr::slice(1:nsnps) %>%
      dplyr::pull(SNP)
  }
) %>%
  purrr::reduce(c) %>%
  unique()

snps_tab <- data.table::fread(
  file.path(RESULTS_DIR, "annovar", "snp2gene_df.csv")
) %>%
  dplyr::left_join(x = data.frame(rsID = keep_snps), y = ., by = "rsID") %>%
  dplyr::mutate(
    Alleles = purrr::map2(Ref, Alt, ~ sort(c(.x, .y))),
    Allele1 = purrr::map_chr(Alleles, ~ .x[[1]]),
    Allele2 = purrr::map_chr(Alleles, ~ .x[[2]]),
    uniqID = paste(Chr, Pos, Allele1, Allele2, sep = ":")
  ) %>%
  dplyr::select(
    rsID, Chr, Pos, uniqID,
    Loci = Gene,
    Function = `Gene Function`, `Exonic Function`
  )
write.csv(
  snps_tab,
  file.path(TABLES_SUPP_DIR, "gwas_filtered_snps.csv"),
  row.names = FALSE
)
