# ==============================================================================
# 11_PCA_by_batch.R
# to check if batch correction adequate,
# plot PCA with: - resid_firstpass (group + covariates)
#                - resid_examinebatch (group + covariates + RUV + batch)
# ==============================================================================



# batch + RUVr
# ------------------------------------------------------------------------------
designmat_examinebatch <-
  model.matrix(~ Group + Age + Sex + Obese + Batch + ruv_correction$W,
               data = metadata_filtered[filter_subjects_to_keep, ])

edgeR_obj_examinebatch <- estimateDisp(edgeR_obj, designmat_examinebatch)
edgeR_obj_examinebatch <- calcNormFactors(edgeR_obj_examinebatch)
edgeR_fit_examinebatch <- glmQLFit(edgeR_obj_examinebatch, designmat_examinebatch)

resid_examinebatch <- residuals(edgeR_fit_examinebatch, type = "deviance")




# plotting fxn, for top 2,500 genes by coefficient of variation
# ------------------------------------------------------------------------------
cv_pca <-
  apply(edgeR_obj$counts, 1, function(x) { sd(x) / mean(x) }) %>%
  order(decreasing = T) %>%
  .[1:2500]


# exprs_in = genes x sample expression matrix of interest
plot_by_batch <- function(exprs_in)
{
  
  pca_unnorm <-
    prcomp(
      exprs_in[cv_pca, ] %>% t
    )
  
  pca_unnorm_expl <- pca_unnorm$sdev^2 %>% `/`(., sum(.)) %>% `*`(100)
  pca_unnorm_expl <- formatC(pca_unnorm_expl, format = "f", digits = 2) %>% paste0(., "%")
  
  ggplot(data = NULL,
         aes(x = pca_unnorm$x[, 1], y = pca_unnorm$x[, 2])) +
    geom_convexhull(aes(fill = metadata_filtered$Batch %>% as.factor()), alpha = 0.2) +
    geom_text(aes(color = metadata_filtered$Batch %>% as.factor(),
                  label = metadata_filtered$Batch %>% as.factor()), size = 3.5,
              show.legend = F) +
    theme_few(base_size = 10) +
    xlab(paste0("PC1 (", pca_unnorm_expl[1], ")")) +
    ylab(paste0("PC2 (", pca_unnorm_expl[2], ")")) +
    scale_color_manual(name = "Batch", values = c("#ef8a62", "#67a9cf")) +
    scale_fill_manual(name = "Batch", values = c("#ef8a62", "#67a9cf")) +
    theme(legend.position = "right")
}



# make plots
# ------------------------------------------------------------------------------

# no batch adjustment, needs script 10 to run
plot_by_batch(resid_firstpass)
ggsave("./FiguresTables/FigS1B.png", width = 4.5, height = 2.5, dpi = 300)

# post batch and RUV adjustment
plot_by_batch(resid_examinebatch)
ggsave("./FiguresTables/FigS2B.png", width = 4.5, height = 2.5, dpi = 300)









