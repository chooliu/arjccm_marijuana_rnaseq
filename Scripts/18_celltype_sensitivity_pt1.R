# ==============================================================================
# 18_celltype_sensitivity_pt1.R
# re-do pre-filtering, RUVSeq k on subjects with cell differentials available
# ==============================================================================



# prep edgeR object & metadata
# ==============================================================================

filter_samples_celltype <- !is.na(metadata_filtered$CellType_lymph)
sum(!filter_samples_celltype)

# do filtering
edgeR_obj_celltype <- edgeR_obj[ , filter_samples_celltype]

# log10-Lymph due to extreme right skewness near zero
# hopefully reduces influence of >10% lymph points
metadata_filtered %<>%
  mutate(Log10Lymph = log10(CellType_lymph + 1))



# determine number of RUV components (RUVr based on first-pass)
# ==============================================================================

dm_firstpass_celltype <-
  model.matrix(~ Group + Age + Sex + Obese + Batch + CellType_lymph,
               data = metadata_filtered[filter_samples_celltype, ])
edgeR_obj_firstpass_celltype <- estimateDisp(edgeR_obj_celltype, dm_firstpass_celltype)
edgeR_obj_firstpass_celltype <- calcNormFactors(edgeR_obj_firstpass_celltype)
edgeR_fit_celltype <- glmQLFit(edgeR_obj_firstpass_celltype, dm_firstpass_celltype)
resid_firstpass_celltype <- residuals(edgeR_fit_celltype, type = "deviance")


evaluate_btwnparticipant_by_k_celltype <-
  function (k, variance_function = iqr,
            resid_in = resid_firstpass, dm_in = dm_firstpass) {
    
    ruv_results <-
      RUVr(
        edgeR_obj_celltype$counts,
        T,
        k = k,
        resid_in)
    
    ruv_results$normalizedCounts %>% # 17602 x 41
      apply(., 1, variance_function) %>% # within-subject variability
      (function(x) { max(x) - min(x) }) # between-subjects
    
  }

plot(1:15,
     sapply(1:15, evaluate_btwnparticipant_by_k_celltype,
            variance_function = iqr,
            resid_in = resid_firstpass_celltype,
            dm_in = dm_firstpass_celltype))

cv_filter_celltype <-
  apply(edgeR_obj_celltype$counts, 1,
        function(x) { sd(x) / mean(x) }) %>%
  order(decreasing = T) %>%
  .[1:2500]

dist_by_sample_celltype <-
  log(edgeR_obj$counts[cv_filter_celltype, filter_samples_celltype] + 1) %>% t %>% dist

evaluate_ruv_k_adonis_celltype <- function(k_to_try) {
  ruv_results <-
    RUVr(
      edgeR_obj_celltype$counts,
      T,
      k = k_to_try,
      resid_firstpass_celltype
    )
  
  adonis(dist_by_sample_celltype ~ ruv_results$W,
         data = NULL, permutations = 1000) %>%
    .$aov.tab %>% .$R2 %>% .[1]
  
}


plot(1:15, sapply(1:15, evaluate_ruv_k_adonis_celltype))









# finally, k = 3 RUVSeq correction
# save counts to new edgeR object (edgeR_obj_RUV)
# ==============================================================================
ruv_correction_celltype <-
  RUVr(
    edgeR_obj_celltype$counts,
    T,
    k = 3,
    resid_firstpass_celltype)

