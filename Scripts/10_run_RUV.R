# ==============================================================================
# 10_run_RUV.R
# afer diagnostics in script 09, run with k = 3 RUVr components
# ==============================================================================



# save primary output of above (commented out)
# k = 3 and "resid_firstpass_batch" is from time-intensive script 09
# ------------------------------------------------------------------------------

# ruv_correction <-
#   RUVr(
#     edgeR_obj$counts,
#     T,
#     k = 3,
#     resid_firstpass_batch)
# 
# save(ruv_correction, file = "Data/Output/RUV.Rdata")


load("Data/Output/RUV.Rdata")



# check correlations/associations between RUVr components
# and known covariates (e.g., cell type, smoking frequency)
# ------------------------------------------------------------------------------
check_corr_with_RUV <- function(variable_to_check) {
  
  ruv_correction$W %>%
    apply(., 2, cor, y = as.numeric(variable_to_check),
          use = "complete", method = "pearson") %>%
    c(., N = sum(!is.na(as.numeric(variable_to_check))))
}

# Table B in Supplemental Material
RUV_Correlations <-
  rbind(
    check_corr_with_RUV(metadata_filtered$CellType_mono_mac),
    check_corr_with_RUV(metadata_filtered$CellType_lymph),
    check_corr_with_RUV(metadata_filtered$MJSmoke_JointYears),
    check_corr_with_RUV(metadata_filtered$Urine_THC_COOH),
    check_corr_with_RUV(metadata_filtered$Urine_THC_gluc),
    check_corr_with_RUV(metadata_filtered$BAL_cell_THC),
    check_corr_with_RUV(metadata_filtered$Tobacco_PackYears),
    check_corr_with_RUV(metadata_filtered$MJSmoke_LastMonth_TimesPerDay_edit),
    check_corr_with_RUV(metadata_filtered$MJVape_LastMonth_TimesPerDay_edit),
    check_corr_with_RUV(metadata_filtered$Tobacco_PacksDay),
    check_corr_with_RUV(metadata_filtered$BMI),
    check_corr_with_RUV(metadata_filtered$FEV1_pred),
    check_corr_with_RUV(metadata_filtered$FVC_pred),
    check_corr_with_RUV(metadata_filtered$FEVFVC)
  ) %>%
  as_tibble() %>%
  bind_cols(., 
            Variable = c("Macrophage/Monocyte (BAL Cell %)", "Lymphocytes (BAL Cell %)", "Marijuana Joint-Years", "THC Concentration: BAL Cell",
                         "THC Concentration: Urine THC-COOH", "THC Concentration: THC-Gluc",
                         "Tobacco Pack-Years", "Marijuana Smoking Frequency (Times/Day)", "Marijuana Vaping Frequency (Times/Day)",
                         "Tobacco Smoking Frequency (Packs/Day)", "BMI", "FEV1 (%p)", "FVC (%p)", "FEV1/FVC")
  ) %>%
  dplyr::select(5, 1:4) %>%
  set_names(c("Variable", "RUV1", "RUV2", "RUV3", "N")) %>%
  mutate_at(2:4, function(x) {
    large <- if_else(abs(x) >= 0.3, "*", "")
    formatC(x, format = "f", digits = 2) %>% paste0(., large) }) %>%
  dplyr::arrange(Variable)  %>%
  dplyr::select(1, 5, 2:4)


# visual inspection of associations above
plot_corr_with_RUV <- function(variable_to_check) {
  par(mfrow = c(3, 1))
  plot(variable_to_check, ruv_correction$W[, 1])
  plot(variable_to_check, ruv_correction$W[, 2])
  plot(variable_to_check, ruv_correction$W[, 3])
}

plot_corr_with_RUV(metadata_filtered$Tobacco_PackYears)
plot_corr_with_RUV(metadata_filtered$Urine_THC_COOH)
