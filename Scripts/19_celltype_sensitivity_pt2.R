# ==============================================================================
# 19_celltype_sensitivity_pt2.R
# run the edgeR models with cell-type
# ==============================================================================




# add RUV from script 18 and re-run edgeR for cell-type analysis
# ------------------------------------------------------------------------------

designmat_celltype <-
  model.matrix(~ Group + Age + Sex + Obese + Batch + Log10Lymph + ruv_correction_celltype$W,
               data = metadata_filtered[filter_samples_celltype, ])
edgeR_obj_celltype <- estimateDisp(edgeR_obj_celltype, designmat_celltype)
edgeR_obj_celltype <- calcNormFactors(edgeR_obj_celltype)
edgeR_fit_celltype <- glmQLFit(edgeR_obj_celltype, designmat_celltype)


# associations with log10(lymph)
# ------------------------------------------------------------------------------
results_lymph <-
  glmQLFTest(edgeR_fit_celltype, coef = "Log10Lymph") %>%
  .$table %>%
  data.frame(ID = row.names(.), .) %>%
  bind_cols(., fdr = p.adjust(.[ , "PValue"], method = "fdr")) %>%
  as_tibble() %>%
  filter(fdr < 0.05) %>%
  arrange(fdr) %>%
  left_join(.,
            annotations_by_gene %>% dplyr::select(-length),
            by = c("ID" = "ensembl_gene_id"))

dim(results_lymph)




# between-group comparisons
# ------------------------------------------------------------------------------

# logFC > 0 means increase in marijuana group
sensitivity_results_MvT <-
  extractPairwiseContrasts(3, edgeR_fit_celltype) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1)

# logFC > 0 means increase in marijuana group
sensitivity_results_MvC <-
  extractPairwiseContrasts(2, edgeR_fit_celltype) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1)

# logFC > 0 makes increase in tobacco group
sensitivity_results_TvC <-
  extractPairwiseContrasts(c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                           edgeR_fit_celltype,
                           method = "contrast")






# logFC > 0 means increase in marijuana group
sensitivity_results_MvT_all <-
  extractPairwiseContrasts(3, edgeR_fit_celltype, fdr_thresh = 1.1) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1) %>%
  full_join(edgeR_test_results_MvT_all, by = "ID") %>%
  mutate(NegLogP.x = -log10(PValue.x),
         NegLogP.y = -log10(PValue.y))

# logFC > 0 means increase in marijuana group
sensitivity_results_MvC_all <-
  extractPairwiseContrasts(2, edgeR_fit_celltype, fdr_thresh = 1.1) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1) %>%
  full_join(edgeR_test_results_MvC_all, by = "ID") %>%
  mutate(NegLogP.x = -log10(PValue.x),
         NegLogP.y = -log10(PValue.y))

# logFC > 0 makes increase in tobacco group
sensitivity_results_TvC_all <-
  extractPairwiseContrasts(c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                           edgeR_fit_celltype,
                           fdr_thresh = 1.1,
                           method = "contrast") %>%
  full_join(edgeR_test_results_TvC_all, by = "ID") %>%
  mutate(NegLogP.x = -log10(PValue.x),
         NegLogP.y = -log10(PValue.y))





examine_signif <- function(joined_table) {
  tmp <-
    joined_table %>%
  transmute(Significance = 
              case_when(fdr.x > 0.05 & fdr.y <= 0.05 ~ "Signif Primary",
                        fdr.x <= 0.05 & fdr.y > 0.05 ~ "Signif Celltype",
                        fdr.x <= 0.05 & fdr.y <= 0.05 ~ "Signif Both",
                        fdr.x > 0.05 & fdr.y > 0.05 ~ "Signif Neither")) %>%
              .$Significance %>%
              table
  tmp["Signif Both"] / sum(tmp[c("Signif Both", "Signif Primary")])
}




# get simple stats & figures for suppl file
# ------------------------------------------------------------------------------

# number significant in both
list(sensitivity_results_MvC_all,
     sensitivity_results_MvT_all,
     sensitivity_results_TvC_all) %>%
  sapply(examine_signif)


# plot p-values
plot_grid(
ggplot(data = sensitivity_results_TvC_all,
         aes(x = NegLogP.x, y = NegLogP.y)) +
  geom_point(alpha = 0.5, color = "grey") + geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(expression("TvC"~-log[10]~P~", With Cell Type"), limits = c(0, 14), expand = c(0, 0)) +
  scale_y_continuous(expression("TvC"~-log[10]~P~", Primary Analysis"), limits = c(0, 14), expand = c(0, 0)) +
  annotate("text", x = 4, y = 13.5,
           label = paste0("Pearson r = ", cor(sensitivity_results_TvC_all$NegLogP.x, sensitivity_results_TvC_all$NegLogP.y, use = "complete") %>% format(digits = 2)))
,

ggplot(data = sensitivity_results_MvC_all,
       aes(x = NegLogP.x, y = NegLogP.y)) +
  geom_point(alpha = 0.5, color = "grey") + geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(expression("MvC"~-log[10]~P~", With Cell Type"), limits = c(0, 14), expand = c(0, 0)) +
  scale_y_continuous(expression("MvC"~-log[10]~P~", Primary Analysis"), limits = c(0, 14), expand = c(0, 0)) +
  annotate("text", x = 4, y = 13.5,
           label = paste0("Pearson r = ", cor(sensitivity_results_MvC_all$NegLogP.x, sensitivity_results_MvC_all$NegLogP.y, use = "complete") %>% format(digits = 2)))

,
ggplot(data = sensitivity_results_MvT_all,
       aes(x = NegLogP.x, y = NegLogP.y)) +
  geom_point(alpha = 0.5, color = "grey") + geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(expression("MvT"~-log[10]~P~", With Cell Type"), limits = c(0, 14), expand = c(0, 0)) +
  scale_y_continuous(expression("MvT"~-log[10]~P~", Primary Analysis"), limits = c(0, 14), expand = c(0, 0)) +
  annotate("text", x = 4, y = 13.5,
           label = paste0("Pearson r = ", cor(sensitivity_results_MvT_all$NegLogP.x, sensitivity_results_MvT_all$NegLogP.y, use = "complete") %>% format(digits = 2))),

nrow = 1)

ggsave("./FiguresTables/FigS3A.png", width = 9, height = 3)
system('open "./FiguresTables/FigS3A.png"')

plot_grid(
  ggplot(data = sensitivity_results_TvC_all,
         aes(x = logFC.x, y = logFC.y)) +
    geom_point(alpha = 0.3, aes(color = fdr.y < 0.05, shape = fdr.y < 0.05)) + geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, lty = 3) + geom_hline(yintercept = 0, lty = 3) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = c("grey", "#67001f")) + scale_shape_manual(values = c(1, 19)) +
    scale_x_continuous(expression("TvC"~logFC~", With Cell Type"), limits = c(-10, 10), expand = c(0, 0)) +
    scale_y_continuous(expression("TvC"~logFC~", Primary Analysis"), limits = c(-10, 10), expand = c(0, 0)) +
    annotate("text", x = -3.5, y = 9.5,
             label = paste0("Pearson r = ", cor(sensitivity_results_TvC_all$logFC.x, sensitivity_results_TvC_all$logFC.y, use = "complete") %>% format(digits = 2)))
  ,
  
  ggplot(data = sensitivity_results_MvC_all,
         aes(x = logFC.x, y = logFC.y)) +
    geom_point(alpha = 0.3, aes(color = fdr.y < 0.05, shape = fdr.y < 0.05)) + geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, lty = 3) + geom_hline(yintercept = 0, lty = 3) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = c("grey", "#67001f")) + scale_shape_manual(values = c(1, 19)) +
    scale_x_continuous(expression("MvC"~logFC~", With Cell Type"), limits = c(-10, 10), expand = c(0, 0)) +
    scale_y_continuous(expression("MvC"~logFC~", Primary Analysis"), limits = c(-10, 10), expand = c(0, 0)) +
    annotate("text", x = -3.5, y = 9.5,
             label = paste0("Pearson r = ", cor(sensitivity_results_MvC_all$logFC.x, sensitivity_results_MvC_all$logFC.y, use = "complete") %>% format(digits = 2)))
  
  ,
  ggplot(data = sensitivity_results_MvT_all,
         aes(x = logFC.x, y = logFC.y)) +
    geom_point(alpha = 0.3, aes(color = fdr.y < 0.05, shape = fdr.y < 0.05)) + geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, lty = 3) + geom_hline(yintercept = 0, lty = 3) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = c("grey", "#67001f")) + scale_shape_manual(values = c(1, 19)) +
    scale_x_continuous(expression("MvT"~logFC~", With Cell Type"), limits = c(-10, 10), expand = c(0, 0)) +
    scale_y_continuous(expression("MvT"~logFC~", Primary Analysis"), limits = c(-10, 10), expand = c(0, 0)) +
    annotate("text", x = -3.5, y = 9.5,
             label = paste0("Pearson r = ", cor(sensitivity_results_MvT_all$logFC.x, sensitivity_results_MvT_all$logFC.y, use = "complete") %>% format(digits = 2))),
  
  nrow = 1)

ggsave("./FiguresTables/FigS3B.png", width = 9, height = 3)
system('open "./FiguresTables/FigS3B.png"')



# any changes in direction?
# ------------------------------------------------------------------------------
sensitivity_results_MvT_all %>%
  filter(fdr.y < 0.05) %>%
  filter(logFC.x*logFC.y < 0) %>%
  dplyr::select(ID, external_gene_name.y)


