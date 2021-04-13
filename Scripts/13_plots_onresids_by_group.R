# ==============================================================================
# 13_plots_onresids_by_group.R
# - extract deviance residuals from edgeR for use in downstream analyses
# - plot by PCA by group (residualized, Fig 1)
# ==============================================================================




# no group edgeR model
# ------------------------------------------------------------------------------
designmat_nogroup <-
  model.matrix(~ Age + Sex + BMI + Batch + ruv_correction$W,
               data = metadata_filtered[filter_subjects_to_keep, ])

edgeR_obj_nogroup <- estimateDisp(edgeR_obj, designmat_nogroup)
edgeR_obj_nogroup <- calcNormFactors(edgeR_obj_nogroup)
edgeR_fit_nogroup <- glmFit(edgeR_obj_nogroup, designmat_nogroup)

resid_nogroup <- residuals(edgeR_fit_nogroup, type = "deviance")











# plot by PCA
# ------------------------------------------------------------------------------

# Dark2 colorbrewer palette
palette_color_bygroup <-
  c(M = "#1B9E77", C = "#D95F02", `T` = "#7570B3",
    `Marijuana (M)` = "#1B9E77", `Control (C)` = "#D95F02", `Tobacco (T)` = "#7570B3")

plot_PCA_by_group <- function(exprs_in)
{
  
  pca_unnorm <-
    prcomp(
      exprs_in %>% t
    )
  
  pca_unnorm_expl <- pca_unnorm$sdev^2 %>% `/`(., sum(.)) %>% `*`(100)
  pca_unnorm_expl <- formatC(pca_unnorm_expl, format = "f", digits = 2) %>% paste0(., "%")
  
  ggplot(data = NULL,
         aes(x = pca_unnorm$x[, 1], y = pca_unnorm$x[, 2])) +
    geom_convexhull(aes(fill = metadata_filtered$Group_Detailed), alpha = 0.2) +
    geom_text(aes(color = metadata_filtered$Group_Detailed, label = metadata_filtered$Group),
              size = 3.5, show.legend = F) +
    theme_few(base_size = 9) +
    xlab(paste0("PC1 (", pca_unnorm_expl[1], ")")) +
    ylab(paste0("PC2 (", pca_unnorm_expl[2], ")")) +
    scale_color_manual(name = "Group", values = palette_color_bygroup) +
    scale_fill_manual(name = "Group", values = palette_color_bygroup) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
    guides(size = FALSE)
  
}





# plot gene counts/resiudals by group
# ------------------------------------------------------------------------------
edger_normfact <- edgeR_fit_RUV$samples

plot_gene_by_group <-
  function(edgeR_results_df = edgeR_test_results_MvC,
           index = NULL, ensg = NULL, residuals = T, custom_description = NULL,
           beecex_custom = 0.4, ylabel_custom = NULL, groups_to_include = c("M", "T", "C")) {
    
    if (residuals == T) {
      counts <- resid_nogroup
      ylabel <- if_else(is.null(ylabel_custom),
                        "Residualized Gene Expression", ylabel_custom)
    } else {
      counts <- edgeR_obj$counts
      ylabel <- if_else(is.null(ylabel_custom),
                        "Normalized Gene Counts", ylabel_custom)
    }
    
    if ( is.null(ensg)) {
      ensg <- edgeR_results_df[index, ]$ID
    }
    
    df <-
      bind_cols(
        Group = metadata_filtered$Group,
        y = counts[ensg, ] * edger_normfact$norm.factors
      ) %>%
      filter(Group %in% groups_to_include)
    
    annot <- annotations_by_gene %>% filter(ensembl_gene_id == ensg)
    
    ggboxplotWrapper(df, "Group", "y", font_size = 10, beesize = 2, beecex = beecex_custom) +
      labs(
        title =
          paste0(
            annot$external_gene_name, " (",
            ifelse(is.null(custom_description),
                   annot$description %>% str_trim(),
                   custom_description
            ), ")"
          )) +
      theme_few(base_size = 9) +
      theme(legend.position = "none") +
      theme(legend.title = element_text(size = 9)) +
      scale_color_manual(name = "Smoking Status",
                         values = c("#1b9e77", "#d95f02", "#7570b3")) +
      xlab("Smoking Group") +
      scale_x_discrete(name = "") +
      scale_y_continuous(
        name = ylabel,
        limits = ylimits,
        breaks =
          if (residuals) {
            ylimits <- c(
              floor(min(df$y) / 0.5) * 0.5,
              ceiling(max(df$y) / 0.5) * 0.5
            )
            seq(ylimits[1], ylimits[2], length.out = 5) }
        else { 
          ylimits <- c(
            floor(min(df$y) / 100) * 100,
            ceiling(max(df$y) / 100) * 100
          )
          seq(ylimits[1], ylimits[2], length.out = 5) %>% round(digits = 0) })
    
  }






# make example plots
# ------------------------------------------------------------------------------

# AJRCCM letter Figure 2: PCA by group **
plot_PCA_by_group(resid_nogroup)
ggsave("./FiguresTables/AJRCC_Figure2.png",
       width = 4, height = 3.5)

# examples of checking top DE/selected genes of interest below
genes_to_check <-
  edgeR_test_results_MvC %>%
  filter(external_gene_name %in% c("HIF1A", "CASP1", "IFNG", "CYP1A1")) %>%
  arrange(external_gene_name)

# top MvT gene
plot_gene_by_group(edgeR_test_results_MvT,
                   index = 1)

# plot CASP1 with normalized counts, compare to residualized
plot_gene_by_group(edgeR_test_results_MvC,
                   ensg = genes_to_check$ID[1],
                   custom_description = formatC(genes_to_check$logFC[1]^2, digits = 2) %>%
                     paste0("MvC FC = ", .),
                   residuals = F)
plot_gene_by_group(edgeR_test_results_MvC,
                   ensg = genes_to_check$ID[1],
                   custom_description = formatC(genes_to_check$logFC[1]^2, digits = 2) %>%
                     paste0("MvC FC = ", .))

