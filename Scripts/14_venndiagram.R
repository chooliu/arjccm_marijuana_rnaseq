# ==============================================================================
# 14_venndiagram.R
# examine overlap between pair-wise comparisons 
# ==============================================================================




# generating venn diagram (using eulerr package)
# ------------------------------------------------------------------------------

# ultimately created venn externally for better graphical cntl
# (20210217-VennDiagram.png = AJRCCM letter Figure 1)
# but numbers in the figure were generated here + checked manually
euler(
  list(`Marijuana vs Non-Smoker (MvC)` = edgeR_test_results_MvC$ID,
       `Marijuana vs Tobacco (MvT)` = edgeR_test_results_MvT$ID,
       `Tobacco vs Control (TvC)` = edgeR_test_results_TvC$ID),
  shape = "ellipse") %>%
  plot(labels = T, quantities = T)




# exploring degree of set overlap
# (i.e., is there a shared vs specific tobacco signature)
# ------------------------------------------------------------------------------


# number of significant results
c(MvC = nrow(edgeR_test_results_MvC),
  MvT = nrow(edgeR_test_results_MvT),
  TvC = nrow(edgeR_test_results_TvC))

# percent of TvC genes also differentially expressed in MvC
intersect(edgeR_test_results_MvC$external_gene_name,
          edgeR_test_results_TvC$external_gene_name) %>%
  length %>%
  `/`(., edgeR_test_results_TvC$external_gene_name %>% length)

# plot effect estimates by pairwise comparison
MvC_TvC_specificity <-
  left_join(edgeR_test_results_MvC_all,
            edgeR_test_results_TvC_all,
            by = c("ID", "external_gene_name"), suffix = c(".MvC", ".TvC")) %>%
  mutate(Significance = case_when(
    fdr.MvC < 0.05 & fdr.TvC < 0.05 ~ "MvC & TvC",
    fdr.MvC < 0.05 & fdr.TvC > 0.05 ~ "MvC",
    fdr.MvC > 0.05 & fdr.TvC < 0.05 ~ "TvC",
    fdr.MvC > 0.05 & fdr.TvC > 0.05 ~ "neither") %>%
      fct_relevel(., c("MvC", "TvC", "MvC & TvC", "neither"))) %>%
  mutate(dirsame = logFC.TvC * logFC.MvC > 0) %>%
  mutate(text_label =
           if_else((Significance != "neither" & dirsame & (abs(logFC.MvC) > 4 | abs(logFC.TvC) > 4)) |
                     (Significance != "neither" & !dirsame & (abs(logFC.MvC) > 2 | abs(logFC.TvC) > 2)), external_gene_name, "")) %>%
  mutate(text_label = if_else(duplicated(text_label), "", text_label))


filter_signif_in_either <- MvC_TvC_specificity$Significance != "neither"

# correlation between effects?
cor(MvC_TvC_specificity$logFC.MvC,
    MvC_TvC_specificity$logFC.TvC,
    method = "pearson", use = "complete")

# in effects
cor(MvC_TvC_specificity$logFC.MvC[filter_signif_in_either],
    MvC_TvC_specificity$logFC.TvC[filter_signif_in_either],
    method = "pearson", use = "complete")

# how many effects have same direction?
MvC_TvC_specificity$dirsame %>% table

# same direction, in genes significant in either comparison only
MvC_TvC_specificity$dirsame[filter_signif_in_either] %>% table



# compare "substance specific" effects
# ------------------------------------------------------------------------------
ggplot(data = MvC_TvC_specificity,
       aes(x = logFC.MvC, y = logFC.TvC,
           color = Significance, shape = Significance)) +
  geom_point(aes(alpha = Significance == "neither"), size = 0.75) +
  geom_text_repel(aes(label = text_label), color = "black", size = 2.5, segment.size = 0.25,
                  box.padding = unit(0.1, "lines")) +
  theme_few(base_size = 9) +
  geom_abline(slope = 1, intercept = 0, color = "#555555") +
  geom_hline(yintercept = 0, size = 0.25, color = "#555555") +
  geom_vline(xintercept = 0, size = 0.25, color = "#555555") +
  scale_color_manual(values = c(neither = "black", `MvC` = "#1b9e77", `TvC` = "#7570b3", `MvC & TvC` = "#980043")) +
  scale_shape_manual(values = c(neither = 1, `MvC` = 19, `TvC` = 19, `MvC & TvC` = 19)) +
  scale_alpha_manual(values = c(`TRUE` = 0.2, `FALSE` = 0.4), guide = F) +
  xlab(expression(log[2]~"FC (MvC)")) + ylab(expression(log[2]~"FC (TvC)")) +
  annotate(x = 5, y = -7, geom = "text", label = "Pearson r = 0.56", size = 4) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.spacing.x = unit(1, "pt"))

