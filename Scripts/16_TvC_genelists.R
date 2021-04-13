# ==============================================================================
# 16_TvC_genelists.R
# examine overlap between TvC alveolar macrophage diff ex studies
# note: this was used to explore external validity to other tobacco populations
#       & to assess comparability of our BAL to purified alveolar macs 
# ==============================================================================





# ==============================================================================
# shaykhiev, et al (2008)
# (processing raw GEO data because diff expr results not available)
# ==============================================================================

# Shaykhiev, R. et al. Smoking-Dependent Reprogramming of Alveolar Macrophage
# Polarization: Implication for Pathogenesis of Chronic Obstructive Pulmonary
# Disease. J. Immunol. 183, 2867–2883 (2009).

# # download & extract Affy .cels
# # ------------------------------------------------------------------------------
# library(GEOquery)
# library(oligo)
# library(limma)
#
# dir.create("./Data/GEO_Public_Data")
# dir.create("./Data/GEO_Public_Data/GSE13896")
# gse <- getGEO("GSE13896", GSEMatrix = TRUE, getGPL = T, destdir = "./Data/GEO_Public_Data/GSE13896")
# getGEOSuppFiles("GSE13896", baseDir = "./Data/GEO_Public_Data/GSE13896")
# 
# system("tar xopf ./Data/GEO_Public_Data/GSE13896/GSE13896/GSE13896_RAW.tar")
# system("mv *CEL.gz ./Data/GEO_Public_Data/GSE13896/")
# 
# 
# # oligo pipeline
# # ------------------------------------------------------------------------------# 
# rawData <-
#   tibble(path = paste0("./Data/GEO_Public_Data/GSE13896/")) %>%
#   pmap(.l = ., .f = list.celfiles, listGzipped = T, full.names = T) %>%
#   map(.x = ., .f = read.celfiles)
# 
# normalizedData <-
#   map(.x = rawData, .f = rma)
# 
# hist(rawData[[1]])
# hist(normalizedData[[1]])
# 
# boxplot(rawData[[1]])
# boxplot(normalizedData[[1]])
# 
# save(normalizedData, file = "./Data/GEO_Public_Data/GSE13896/microarray_normalized.Rdata")



# following pre-processing commented out steps above, format for limma
# ------------------------------------------------------------------------------

load("./Data/GEO_Public_Data/GSE13896/microarray_normalized.Rdata")

shaykhiev_exprs <- exprs(normalizedData[[1]])
shaykhiev_pheno <- pData(gse[[1]])
shaykhiev_annotations <- gse$GSE13896_series_matrix.txt.gz@featureData %>% .@data

# use one probe per gene only, top 75%
filter_shaykhiev_genes_noncontrol <-
  shaykhiev_annotations$`Sequence Type` != "Control sequence" &
  shaykhiev_annotations$`Gene Symbol` != ""

filter_shaykhiev_genes_topCV <-
  apply(shaykhiev_exprs, 1, (function(x) { sd(x) / mean(x) } ) ) %>%
  tibble(Probe = shaykhiev_annotations$ID, Gene = shaykhiev_annotations$`Gene Symbol`, CV = .) %>%
  arrange(Gene, -CV) %>%
  filter(!duplicated(Gene)) %>%
  .$Probe %>%
  `%in%`(rownames(shaykhiev_exprs), .)

filter_shaykhiev_genes_topMedian <-
  apply(shaykhiev_exprs, 1, median) %>%
  order(decreasing = T) %>%
  `%in%`(., 1:round(length(.)*0.75))

filter_shaykhiev_genes <-
  filter_shaykhiev_genes_noncontrol &
  filter_shaykhiev_genes_topCV &
  filter_shaykhiev_genes_topMedian


shaykhiev_pheno_plussmokers <-
  pData(gse[[1]]) %>%
  transmute(Age = `Age:ch1` %>% as.numeric,
            Ancestry = `Ancestry:ch1` %>% as.factor,
            Sex = `Sex:ch1`,
            Control = !is.na(shaykhiev_pheno$`Smoking status:ch1`) & (shaykhiev_pheno$`Smoking status:ch1` == "non-smoker"),
            COPD = grepl("COPD", `smoking status:ch1`))
filter_shaykhiev_subjects <-
  !is.na(shaykhiev_pheno_plussmokers$Age) &
  !(shaykhiev_pheno_plussmokers$COPD)
  
# RUVr components (k = 5),
# based on same method of obtaining k for our RNA-seq dat
dm_plussmokers <-
  model.matrix(~ .,
               shaykhiev_pheno_plussmokers %>%
                 filter(filter_shaykhiev_subjects) %>%
                 dplyr::select(Age, Sex, Ancestry, Control))


limma_fit <-
  lmFit(shaykhiev_exprs[filter_shaykhiev_genes, filter_shaykhiev_subjects],
        design = dm_plussmokers)
resid_fit <- residuals.MArrayLM(
  limma_fit, shaykhiev_exprs[filter_shaykhiev_genes, filter_shaykhiev_subjects])

shaykhiev_distmat <-
  shaykhiev_exprs[filter_shaykhiev_genes, filter_shaykhiev_subjects] %>%
  t %>%
  dist()

evaluate_ruv_k <- function(k_to_try) {
  ruv_results <-
    RUVr(
      shaykhiev_exprs[filter_shaykhiev_genes, filter_shaykhiev_subjects],
      T,
      k = k_to_try,
      resid_fit,
      isLog = T
    )
  dm_sva <-
    cbind(dm_plussmokers, ruv_results$W)
  adonis(shaykhiev_distmat ~ dm_sva, data = NULL) %>% .$aov.tab %>% .$R2 %>% .[1]
  
}

ruv_percent_exp <- sapply(1:15, evaluate_ruv_k)
plot(1:15, 1 - ruv_percent_exp)

shaykhiev_ruv <-
  RUVr(
    shaykhiev_exprs[filter_shaykhiev_genes, filter_shaykhiev_subjects],
    T,
    k = 5,
    resid_fit,
    isLog = T
  )


# finally, run limma models & extract results on shaykhiev dat
shaykhiev_pheno_plus_RUV <-
  cbind(dm_plussmokers, RUV = shaykhiev_ruv$W)

limma_fit_plus_RUV <-
  lmFit(shaykhiev_exprs[filter_shaykhiev_genes, filter_shaykhiev_subjects],
        design = shaykhiev_pheno_plus_RUV)
limma_fit_plus_RUV <- eBayes(limma_fit_plus_RUV, trend = T)

TvC_signif_shaykhiev_fdr <- 
  limma_fit_plus_RUV %>%
  topTable(., coef = "ControlTRUE", adjust.method = "fdr",
           number = 20000, p.value = 0.05) %>%
  as_tibble(rownames = "ID") %>%
  left_join(shaykhiev_annotations, by = "ID")

dim(TvC_signif_shaykhiev_fdr)

TvC_signif_shaykhiev_all <- 
  limma_fit_plus_RUV %>%
  topTable(., coef = "ControlTRUE", adjust.method = "fdr",
           number = 30000, p.value = 1) %>%
  as_tibble(rownames = "ID") %>%
  left_join(shaykhiev_annotations, by = "ID")

intersect(edgeR_test_results_TvC$external_gene_name,
          TvC_signif_shaykhiev_fdr$`Gene Symbol`) %>% length

# compare to our TvC results
compare_shaykhiev <-
  left_join(TvC_signif_shaykhiev_all,
            edgeR_test_results_TvC_all, by = c("Gene Symbol" = "external_gene_name")) %>%
  mutate(NeglogP.x = -log10(P.Value),
         NeglogP.y = -log10(PValue)) %>%
  mutate(tested = !is.na(NeglogP.y),
         signif = fdr < 0.05,
         same_dir = logFC.x * logFC.y < 0)

compare_shaykhiev %>%
  filter(!duplicated(`Gene Symbol`)) %>%
  filter(tested) %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(signif, same_dir) %>%
  tally

ggplot(compare_shaykhiev, aes(logFC.x, logFC.y)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1)





# ==============================================================================
# woodruff, et al (2005)
# ==============================================================================

# Woodruff, P. G. et al. A distinctive alveolar macrophage activation state
# induced by cigarette smoking. Am. J. Respir. Crit. Care Med. 172, 1383–1392 (2005).

# download supplemental table E1

woodruff_TvC <-
  read_excel("./Data/Smoker_vs_Unsmokers/Woodruff_2005_PMC2718436.xls", sheet = 1) %>%
  left_join(., edgeR_test_results_TvC_all, c("Symbol" = "external_gene_name")) %>%
  mutate(tested = !is.na(logFC),
         signif = fdr < 0.05,
         same_dir = as.factor(`Fold difference`*logFC > 0))

# note that some gene symbols are duplicated in multiple probes,
# but the direction of fold-change is the same 
woodruff_TvC %>%
  group_by(Symbol) %>%
  summarize(FC = unique(`Fold difference`) %>% paste(collapse = " // ")) %>%
  filter(grepl("//", FC))

woodruff_TvC %>%
  filter(!duplicated(Symbol)) %>%
  filter(tested) %>%
  group_by(signif, same_dir) %>%
  summarize(n())





# ==============================================================================
# morrow, et al (2019)
# ==============================================================================

# Morrow, J. D. et al. RNA-sequencing across three matched tissues reveals shared
# and tissue-specific gene expression and pathway signatures of COPD.
# Respir. Res. 20, 1–12 (2019).

# supplementary file Table S5
morrow_TvC <-
  read_excel("./Data/Emphysema/Table_S5.xls", sheet = "smoking") %>%
  filter(padj < 0.05) %>%
  left_join(., edgeR_test_results_TvC_all, by = c("...1" = "ID")) %>%
  mutate(tested = !is.na(logFC),
         signif = fdr < 0.05,
         same_dir = log2FoldChange * logFC > 0)

morrow_TvC %>%
  filter(tested) %>%
  group_by(signif, same_dir) %>%
  summarize(n())





# ==============================================================================
# 20 TvC genes across all four studies
# (or 94, omitting woodruff due to conservative bonferroni correction)
# ==============================================================================

intersect(morrow_TvC$symbol, TvC_signif_shaykhiev_fdr$`Gene Symbol`) %>%
   intersect(woodruff_TvC$Symbol) %>%
  .[`%in%`(., edgeR_test_results_TvC$external_gene_name)] 
