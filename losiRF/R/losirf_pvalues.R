


# print("Evaluting local RF feature stability...")
# load(file.path(save_path, "irf_interaction_fit.Rdata"))
#
# feature_groups <- snps.df %>%
#   select(feature = Name, group = Gene)
#
# # compute local stability
# stab_tr2_ls <- map(
#   irf_gene_out$rf.list,
#   ~ localFeatureStabilityRF(
#     rf_fit = .x,
#     X = out$geno_train2,
#     feature_groups = feature_groups
#   )
# )
#
# stab_va_ls <- map(
#   irf_gene_out$rf.list,
#   ~ localFeatureStabilityRF(
#     rf_fit = .x,
#     X = out$geno_test,
#     feature_groups = feature_groups
#   )
# )
#
# # evaluate via permutation test
# test_genes_ls <- map(
#   stab_tr2_ls,
#   ~ data.frame(
#     gene = colnames(.),
#     mean_stability = apply(., 2, mean)
#   ) %>%
#     arrange(desc(mean_stability)) %>%
#     dplyr::slice(1:test_ngenes) %>%
#     pull(gene) %>%
#     as.character() %>%
#     setNames(., .)
# )
#
# perm_out_ls <- map2(
#   test_genes_ls, stab_va_ls,
#   function(test_genes, stab_va) {
#     map(
#       test_genes,
#       ~ runPermutationTest(
#         stab_df = stab_va,
#         y = out$pheno_test,
#         feature = .x, nperm = 1e4
#       )
#     )
#   }
# )
#
# print("Evaluating local RF interaction stability...")
# load(file.path(save_path, "irf_interaction_fit.Rdata"))
#
# feature_groups <- snps.df %>%
#   select(feature = Name, group = Gene)
#
# # select interactions to tests
# irf_ints_ls <- map(irf_gene_out$interaction, ~ annotInts(.x, snps.df, F))
# irf_keep_ints_ls <- map(
#   irf_ints_ls,
#   ~ .x %>%
#     filter(
#       stability >= 0.5,
#       sta.fsd > 0,
#       sta.mip > 0
#     ) %>%
#     arrange(desc(prevalence)) %>%
#     slice(1:min(n(), test_nints))
# )
#
# # compute local interaction stability
# int_stab_tr2_ls <- map2(
#   irf_gene_out$rf.list, irf_keep_ints_ls,
#   ~ localIntStabilityRF(
#     rf_fit = .x,
#     X = out$geno_train2,
#     ints = .y %>% pull(int),
#     feature_groups = feature_groups
#   )
# )
#
# int_stab_va_ls <- map2(
#   irf_gene_out$rf.list, irf_keep_ints_ls,
#   ~ localIntStabilityRF(
#     rf_fit = .x,
#     X = out$geno_test,
#     ints = .y %>% pull(int),
#     feature_groups = feature_groups
#   )
# )
#
# # evaluate via permutation test
# test_ints_ls <- map(
#   int_stab_tr2_ls,
#   ~ data.frame(
#     int = colnames(.),
#     mean_stability = apply(., 2, mean)
#   ) %>%
#     arrange(desc(mean_stability)) %>%
#     dplyr::slice(1:min(n(), test_nints)) %>%
#     pull(int) %>%
#     as.character() %>%
#     setNames(., .)
# )
#
# int_perm_out_ls <- map2(
#   test_ints_ls, int_stab_va_ls,
#   function(test_ints, stab_va) {
#     map(
#       test_ints,
#       ~ runPermutationTest(
#         stab_df = stab_va,
#         y = out$pheno_test,
#         feature = .x, nperm = 1e4
#       )
#     )
#   }
# )
#
# # save results
# save(int_stab_tr2_ls, int_stab_va_ls, test_ints_ls, int_perm_out_ls, irf_ints_ls,
#      file = file.path(save_path, "local_interaction_stability_results.Rdata")
# )
# print("Completed.")
#
