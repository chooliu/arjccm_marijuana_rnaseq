# ==============================================================================
# 17_rtPCR_validation.R
# load data (two plates, each diff probes/gene targets)
# analyze RT-PCR by group and compare to RNA-seq
# ==============================================================================




# load RT-PCR data
# ==============================================================================

# load plate 1 --> convert to wide format for analysis


rtPCR_plate1 <-
  read_excel("./data/rtPCR/Burnham MJ RNA Seq plate 1.xls", sheet = 4, skip = 34) %>%
  filter(!is.na(`Sample Name`) & !(`Sample Name` == "NTC") & !(CT == "Undetermined")) %>%
  dplyr::select(`Sample Name`, `Target Name`, CT) %>%
  mutate(CT = as.numeric(CT)) %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarize(CT = mean(CT))

rtPCR_curvecheck1 <-
  read_excel("./data/rtPCR/Burnham MJ RNA Seq plate 1.xls", sheet = 2, skip = 34)

rtPCR_plate1 <-
  rtPCR_plate1 %>%
  spread(`Target Name`, CT) %>%
  dplyr::rename(HIF1A = HIF1a) %>%
  mutate(dARNT = 2^(TBP - ARNT),
         dHIF1A = 2^(TBP - HIF1A),
         dHMOX1 = 2^(TBP - HMOX1))

rtPCR_plate1$`Sample Name`[rtPCR_plate1$`Sample Name` == "MJ107"] <- "MF107"
rtPCR_plate1 %<>%
  arrange(`Sample Name`)

rtPCR_plate1 %<>%
  inner_join(., metadata_filtered, by = c("Sample Name" = "ID"))



# load plate 2 --> convert to wide format for analysis
# ------------------------------------------------------------------------------
rtPCR_plate2 <-
  read_excel("./data/rtPCR/Burnham RNA Seq Plate 2 050418.xls", sheet = 4, skip = 34)

rtPCR_plate2 %<>%
  filter(!is.na(`Sample Name`) & !(`Sample Name` == "NTC") & !(CT == "Undetermined")) %>%
  dplyr::select(`Sample Name`, `Target Name`, CT, Well) %>% 
  mutate(CT = as.numeric(CT)) %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarize(CT = median(CT))

rtPCR_plate2 <-
  rtPCR_plate2 %>%
  spread(`Target Name`, CT) %>%
  mutate(dADRB2 = 2^(`TBP-VIC` - ADRB2),
         dCASP1 = 2^(`TBP-VIC` - CASP1),
         dSIRT1 = 2^(`TBP-VIC` - SIRT1),
         dTLR6 = 2^(`TBP-VIC` - TLR6))

# typo "MJ" / "MF"
rtPCR_plate2$`Sample Name`[rtPCR_plate2$`Sample Name` == "MJ107"] <- "MF107"
rtPCR_plate2 %<>%
  arrange(`Sample Name`) %>%
  inner_join(., metadata_filtered, by = c("Sample Name" = "ID"))




# finally join plate 1 and 2, in wide format
# extract edgeR normalization factors for samples to plot on count-scale
# ------------------------------------------------------------------------------
rtPCR_data <- inner_join(rtPCR_plate1, rtPCR_plate2)

filter_samples_rtpcr <- metadata_filtered$ID %in% rtPCR_plate1$`Sample Name`
edger_normfact <- edgeR_fit_RUV$samples[filter_samples_rtpcr , 3]









# plot each RT-PCR target (Figure S4)
# ------------------------------------------------------------------------------
dir.create("./FiguresTables/rtPCR")
make_plots_and_tables_rtPCR <- function(gene_symbol, make_images = T) {
  ensg_id <- annotations_by_gene %>% filter(external_gene_name == gene_symbol) %>% .$ensembl_gene_id
  rtpcr_target <-
    tibble(ID = metadata_filtered %>% filter(filter_samples_rtpcr) %>% .$ID,
           rtPCR_val =  rtPCR_data %>% .[ , paste0("d", gene_symbol), drop = T],
           RNAseq_val = edgeR_fit_RUV$counts[ensg_id, filter_samples_rtpcr] * edger_normfact
    ) %>%
    left_join(metadata_filtered, by = "ID")
  
  target_full <- lm(rtPCR_val ~ Age + Sex + Group + Obese, data = rtpcr_target)
  target_partial <- lm(rtPCR_val ~ Age + Sex + Obese, data = rtpcr_target)
  anova(target_full, target_partial)
  summary(target_full)
  
  rtpcr_target %<>% cbind(rtPCR_resid = resid(target_partial),
                          RNAseq_resid = resid_nogroup[ensg_id, filter_samples_rtpcr])
  
  if (make_images) {
  
  fig_target_1 <-
    ggboxplotWrapper(rtpcr_target, "Group", "rtPCR_val") + ylab(expression("RT-PCR"~2^{Delta*Ct})) +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() + theme(legend.position = 'none', legend.title = element_text(size = 10))
  fig_target_2 <-
    ggboxplotWrapper(rtpcr_target, "Group", "RNAseq_val", ylabel = "RNA-seq Counts") +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() + theme(legend.position = 'none', legend.title = element_text(size = 10))
  fig_target_3 <-
    ggboxplotWrapper(rtpcr_target, "Group", "rtPCR_resid", ylabel = "RT-PCR Residuals") +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() + theme(legend.position = 'none', legend.title = element_text(size = 10))
  fig_target_4 <-
    ggboxplotWrapper(rtpcr_target, "Group", "RNAseq_resid", ylabel = "RNA-seq Residuals") +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() + theme(legend.position = 'none', legend.title = element_text(size = 10))
  fig_target_5 <-
    ggplot(data = rtpcr_target, aes(x = rtPCR_val, y = RNAseq_val)) +
    geom_smooth(method = 'lm', alpha = 0.5, color = "#d73027", se = F, size = 1) +
    annotate(geom = "text", label = cor(rtpcr_target$RNAseq_val, rtpcr_target$rtPCR_val, use = "complete")) +
    geom_point(alpha = 0.6) + xlab(expression(2^{Delta*Ct})) + ylab("RNA-seq Counts") +
    theme_classic()
  
    compiled_plots <-
      cowplot::plot_grid(
        fig_target_2, fig_target_1, fig_target_5,
        nrow = 1)
    
    title <- ggdraw() + 
      draw_label(
        gene_symbol,
        x = 0,
        hjust = 0
      ) +
      theme(
        plot.margin = margin(0, 0, 0, 7)
      )
    plot_grid(
      title, compiled_plots,
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
  
  ggsave(filename = paste0("./FiguresTables/rtPCR/", gene_symbol, ".png"), dpi = 300, width = 7.5, height = 2)
}
  
  
  rtpcr_model <-
    paste0("d", gene_symbol, " ~ Group + Age + Sex + Obese") %>%
    as.formula %>%
    lm(., rtPCR_data) %>%
    summary() %>%
    coefficients
  
  rtpcr_model_log <-
    paste0("log2(d", gene_symbol, ") ~ Group + Age + Sex + Obese") %>%
    as.formula %>%
    lm(., rtPCR_data) %>%
    summary() %>%
    coefficients
  
  rtpcr_model_2 <-
    paste0("d", gene_symbol, " ~ I(Group == 'M') + I(Group == 'C') + Age + Sex + Obese") %>%
    as.formula %>%
    lm(., rtPCR_data) %>%
    summary() %>%
    coefficients
  
  rtpcr_model_2_log <-
    paste0("log2(d", gene_symbol, ") ~ I(Group == 'M') + I(Group == 'C') + Age + Sex + Obese") %>%
    as.formula %>%
    lm(., rtPCR_data) %>%
    summary() %>%
    coefficients
  
  c(rtpcr_MvC_b = rtpcr_model_log["GroupC", "Estimate"] %>% `*`(-1),
    rtpcr_MvT_b = rtpcr_model_log["GroupT", "Estimate"] %>% `*`(-1),
    rtpcr_TvC_b = rtpcr_model_2_log['I(Group == "C")TRUE', "Estimate"] %>% `*`(-1),
    rtpcr_MvC_p = rtpcr_model["GroupC", "Pr(>|t|)"],
    rtpcr_MvT_p = rtpcr_model["GroupT", "Pr(>|t|)"],
    rtpcr_TvC_p = rtpcr_model_2['I(Group == "C")TRUE', "Pr(>|t|)"],
    
    rnaseq_MvC_b = edgeR_test_results_MvC_all %>% filter(external_gene_name == gene_symbol) %>% .$logFC,
    rnaseq_MvT_b = edgeR_test_results_MvT_all %>% filter(external_gene_name == gene_symbol) %>% .$logFC,
    rnaseq_TvC_b = edgeR_test_results_TvC_all %>% filter(external_gene_name == gene_symbol) %>% .$logFC,
    rnaseq_MvC_p = edgeR_test_results_MvC_all %>% filter(external_gene_name == gene_symbol) %>% .$fdr,
    rnaseq_MvT_p = edgeR_test_results_MvT_all %>% filter(external_gene_name == gene_symbol) %>% .$fdr,
    rnaseq_TvC_p = edgeR_test_results_TvC_all %>% filter(external_gene_name == gene_symbol) %>% .$fdr,
    pearson = cor(rtpcr_target[ , "rtPCR_val", drop = T], rtpcr_target$RNAseq_val) 
  )
  
}





# summarize p-values / effects in supplementary table
# ------------------------------------------------------------------------------
rtpcr_validation_targets <- c("ADRB2", "ARNT", "CASP1", "HIF1A", "SIRT1", "TLR6")
rtPCR_data %<>% mutate(dTBP = TBP)
make_plots_and_tables_rtPCR("TBP", make_images = T)


TableB <-
  rtpcr_validation_targets %>%
  sapply(make_plots_and_tables_rtPCR, make_images = T) %>%
  as_tibble(rownames = "Feature") %>%
  gather(Gene, Value, -1) %>%
  spread(Feature, Value) %>%
  mutate_at(names(.) %>% .[grepl("_b", .)], function(x) { formatC( x, format = "f", digits = 2) }) %>%
  mutate_at(names(.) %>% .[grepl("_p", .)], function(x) { if_else( x < 0.05, " *", "") }) %>%
  transmute(Gene = Gene,
            RNASeq_MvC = paste0(rnaseq_MvC_b, rnaseq_MvC_p),
            RNASeq_MvT = paste0(rnaseq_MvT_b, rnaseq_MvT_p),
            RNASeq_TvC = paste0(rnaseq_TvC_b, rnaseq_TvC_p),
            RTPCR_MvC = paste0(rtpcr_MvC_b, rtpcr_MvC_p),
            RTPCR_MvT = paste0(rtpcr_MvT_b, rtpcr_MvT_p),
            RTPCR_TvC = paste0(rtpcr_TvC_b, rtpcr_TvC_p),
            Correlation = formatC(pearson, digits = 2))





# check house-keeping gene
# ------------------------------------------------------------------------------

# extract TBP
ensg_tbp <-
  annotations_by_gene %>%
  filter(external_gene_name == "TBP") %>%
  .$ensembl_gene_id
map(list(edgeR_test_results_MvC_all, edgeR_test_results_MvT_all, edgeR_test_results_TvC_all),
    ~ filter(.x, ID == ensg_tbp))

rtPCR_data %<>%
  ungroup() %>%
  mutate(TBP_RNAseq = edgeR_fit_RUV$counts[ensg_tbp, filter_samples_rtpcr] * edger_normfact,
         TBP_RNAseq_resid = resid_nogroup[ensg_tbp, filter_samples_rtpcr])

# TBP reference value correlated between plates (Pearson r = 0.8)
ggplot(data = rtPCR_data, aes(x = TBP, y = `TBP-VIC`)) + geom_point() +
  xlab("rtPCR CT, Plate 1") + ylab("rtPCR CT, Plate 2") + theme_classic()
cor(rtPCR_data$TBP, rtPCR_data$`TBP-VIC`)

# TBP RT-PCR value not highly correlated to RNA-seq counts of TBP
ggplot(data = rtPCR_data, aes(x = TBP, y = TBP_RNAseq)) + geom_point()

# plot by group
ggplot(data = rtPCR_data, aes(x = Group, y = TBP_RNAseq)) + geom_point()


plot_grid(
  ggboxplotWrapper(rtPCR_data, "Group", "TBP_RNAseq") +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() +
    theme(legend.position = 'none', legend.title = element_text(size = 10)) +
    ggtitle("TBP, RNA-Seq") + ylab("RNA-Seq Counts"),
ggboxplotWrapper(rtPCR_data, "Group", "TBP") +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() +
  theme(legend.position = 'none', legend.title = element_text(size = 10)) +
  ggtitle("TBP, RT-PCR Plate 1") + ylab(expression(C[T]~"(Plate 1)")),

ggplot(data = rtPCR_data, aes(x = TBP, y = TBP_RNAseq)) +
  geom_smooth(method = 'lm', alpha = 0.5, color = "#d73027", se = F, size = 1) +
  annotate(geom = "text", label = cor(rtPCR_data$TBP, rtPCR_data$TBP_RNAseq, use = "complete")) +
  geom_point(alpha = 0.6) + xlab(expression(C[T])) + ylab("RNA-seq Counts") +
  theme_classic(),
nrow = 1)
ggsave("./FiguresTables/FigS5.png", width = 7.5, height = 2)

ggboxplotWrapper(rtPCR_data, "Group", "`TBP-VIC`") +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) + theme_classic() + theme(legend.position = 'none', legend.title = element_text(size = 10)) +
  ggtitle("TBP (Housekeeping Gene), rtPCR Plate 2") + ylab("Ct (Plate 2)")


# does TBP vary by group, w/ and w/o covariate adjustment?
anova(lm(TBP ~ Age + Sex + Group + Obese, data = rtPCR_data),
      lm(TBP ~ Age + Sex + Obese, data = rtPCR_data))

anova(lm(TBP ~ Group, data = rtPCR_data),
      lm(TBP ~ 1, data = rtPCR_data))






