# ==============================================================================
# 12_edgeR.R
# run modeling / pairwise contrasts, & export all results to Excel
# ==============================================================================




# primary analysis
# ==============================================================================

# design matrix
# ------------------------------------------------------------------------------
designmat_RUV <-
  model.matrix(~ Group + Age + Sex + Obese + Batch + ruv_correction$W,
    data = metadata_filtered[filter_subjects_to_keep, ]
  )
edgeR_obj_RUV <- estimateDisp(edgeR_obj, designmat_RUV)
edgeR_obj_RUV <- calcNormFactors(edgeR_obj_RUV)
edgeR_fit_RUV <- glmQLFit(edgeR_obj_RUV, designmat_RUV)


# function for pairwise comparisons
# ------------------------------------------------------------------------------
extractPairwiseContrasts <-
  function(input, edgeR_output, method = "coef", fdr_thresh = 0.05) {
    if (method == "coef") {
      edgeR_output$design %>%
        colnames() %>%
        .[input] %>%
        print()
      results <-
        glmQLFTest(edgeR_output, coef = input)
    }

    if (method == "contrast") {
      edgeR_output$design %>%
        colnames() %>%
        .[input != 0] %>%
        print()
      results <-
        glmQLFTest(edgeR_output, contrast = input)
    }


    results$table %>%
      data.frame(ID = row.names(.), .) %>%
      bind_cols(., fdr = p.adjust(.[, "PValue"], method = "fdr")) %>%
      as_tibble() %>%
      filter(fdr < fdr_thresh) %>%
      arrange(fdr) %>%
      left_join(.,
        annotations_by_gene %>% dplyr::select(-length),
        by = c("ID" = "ensembl_gene_id")
      )
  }




# run pairwise comparisons -- FDR < 0.05
# ------------------------------------------------------------------------------

# logFC > 0 means increase in marijuana group relative to tobacco
edgeR_test_results_MvT <-
  extractPairwiseContrasts(3, edgeR_fit_RUV) %>%
  mutate(
    logFC = logFC * -1,
    logCPM = logCPM * -1
  )

# logFC > 0 means increase in marijuana group relative to controls
edgeR_test_results_MvC <-
  extractPairwiseContrasts(2, edgeR_fit_RUV) %>%
  mutate(
    logFC = logFC * -1,
    logCPM = logCPM * -1
  )

# logFC > 0 makes increase in tobacco group relative to controls
edgeR_test_results_TvC <-
  extractPairwiseContrasts(c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
    edgeR_fit_RUV,
    method = "contrast"
  )








# all test results for various plotting purposes
# ------------------------------------------------------------------------------

# logFC > 0 means increase in marijuana group
edgeR_test_results_MvT_all <-
  extractPairwiseContrasts(3, edgeR_fit_RUV, fdr_thresh = 1.1) %>%
  mutate(
    logFC = logFC * -1,
    logCPM = logCPM * -1
  )

# logFC > 0 means increase in marijuana group
edgeR_test_results_MvC_all <-
  extractPairwiseContrasts(2, edgeR_fit_RUV, fdr_thresh = 1.1) %>%
  mutate(
    logFC = logFC * -1,
    logCPM = logCPM * -1
  )

# logFC > 0 makes increase in tobacco group
edgeR_test_results_TvC_all <-
  extractPairwiseContrasts(c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
    edgeR_fit_RUV,
    method = "contrast", fdr_thresh = 1.1
  )





# RUV-Seq associations [scaled to one standard deviation in factor values]
# ==============================================================================
RUV_associations_1 <-
  extractPairwiseContrasts(c(0, 0, 0, 0, 0, 0, 0, 1 / sd(ruv_correction$W[, 1]), 0, 0),
    edgeR_fit_RUV,
    method = "contrast", fdr_thresh = 0.05
  )
RUV_associations_2 <-
  extractPairwiseContrasts(c(0, 0, 0, 0, 0, 0, 0, 0, 1 / sd(ruv_correction$W[, 2]), 0),
    edgeR_fit_RUV,
    method = "contrast", fdr_thresh = 0.05
  )
RUV_associations_3 <-
  extractPairwiseContrasts(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1 / sd(ruv_correction$W[, 3])),
    edgeR_fit_RUV,
    method = "contrast", fdr_thresh = 0.05
  )

# plot associations
plot_corr_with_RUV(edgeR_obj$counts[RUV_associations_1$ID[1], ])




# export differential results
# ==============================================================================

prepare_for_export <-
  function(table) {
    table %>%
      dplyr::select(1, 7, 8, 2, 4, 5, 6) %>%
      set_names(c("Ensembl_ID", "Symbol", "Description", "Log2FC", "F", "PValue", "FDR"))
  }



summary_alltestedgenes <-
  annotations_by_gene %>%
  filter(ensembl_gene_id %in% rownames(edgeR_fit_RUV)) %>%
  dplyr::select(1:3) %>%
  arrange(ensembl_gene_id) %>%
  set_names(., c("ID", "Symbol", "Description")) %>%
  left_join(
    .,
    edgeR_test_results_MvC_all %>%
      transmute(ID, MvC_logFC = logFC, MvC_FDR = fdr,
                MvC_Signif = if_else(fdr < 0.05, "*", ""))
  ) %>%
  left_join(
    .,
    edgeR_test_results_MvT_all %>%
      transmute(ID, MvT_logFC = logFC, MvT_FDR = fdr,
                MvT_Signif = if_else(fdr < 0.05, "*", ""))
  ) %>%
  left_join(
    .,
    edgeR_test_results_TvC_all %>%
      transmute(ID, TvC_logFC = logFC, TvC_FDR = fdr,
                TvC_Signif = if_else(fdr < 0.05, "*", ""))
  ) %>%
  dplyr::rename(Ensembl_ID = ID)











# no covariate-adjustment 
# ==============================================================================


# design matrix
# ------------------------------------------------------------------------------
designmat_RUV <-
  model.matrix(~ Group + Age + Sex + Obese + Batch + ruv_correction$W,
               data = metadata_filtered[filter_subjects_to_keep, ]
  )
edgeR_obj_RUV <- estimateDisp(edgeR_obj, designmat_RUV)
edgeR_obj_RUV <- calcNormFactors(edgeR_obj_RUV)
edgeR_fit_RUV <- glmQLFit(edgeR_obj_RUV, designmat_RUV)




# final export (ResultsTable1 - all DEG results)
write_xlsx(
  list(
    MvC = edgeR_test_results_MvC %>% prepare_for_export(),
    MvT = edgeR_test_results_MvT %>% prepare_for_export(),
    TvC = edgeR_test_results_TvC %>% prepare_for_export(),
    RUV1 = RUV_associations_1 %>% prepare_for_export(),
    RUV2 = RUV_associations_2 %>% prepare_for_export(),
    RUV3 = RUV_associations_3 %>% prepare_for_export(),
    AllTestedGenes = summary_alltestedgenes
  ),
  path = "./FiguresTables/ResultsTable1.xlsx",
  format_headers = F
)

system('open "./FiguresTables/ResultsTable1.xlsx"')
