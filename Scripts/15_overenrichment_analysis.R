# ==============================================================================
# 15_overenrichment_analysis.R
# run WebGestalt on each pairwise gene set
# ==============================================================================




# # wrapper function for webgestalt
# # ==============================================================================
# web_gestalt_wrapper <-
#   function(df, project,
#            db = c("geneontology_Biological_Process"),
#            error_check_on = T) {
# 
#   if (error_check_on) {
#     if (dir.exists(paste0("Output/Project_", project))) {
#       return(paste0("Output/Project_", project, " -- already queried"))
#     }
#   }
# 
#   output <-
#     WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
#                 enrichDatabase = db, enrichDatabaseFile = NULL,
#                 enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL,
#                 interestGeneFile = NULL, interestGene = df,
#                 interestGeneType = "genesymbol", collapseMethod = "mean",
#                 referenceGeneFile = NULL, referenceGene = all_genes,
#                 referenceGeneType = "genesymbol", referenceSet = NULL, minNum = 5,
#                 maxNum = 500, sigMethod = "fdr", fdrMethod = "BH", fdrThr = 1.1,
#                 topThr = 10, reportNum = 100, perNum = 5000, isOutput = TRUE,
#                 outputDirectory = getwd(), projectName = project,
#                 dagColor = "continuous", setCoverNum = 10,
#                 networkConstructionMethod = NULL, neighborNum = 10,
#                 highlightType = "Seeds", highlightSeedNum = 10, nThreads = 1,
#                 cache = NULL, hostName = "http://www.webgestalt.org/")
#   paste0("mv -f Output/Project_", project, " Output") %>% system
# 
#   output
# }
# 
# 
# 
# # run webgestalt enrichment analysis on each pairwise
# # ------------------------------------------------------------------------------
# MvT_pathways <-
#   edgeR_test_results_MvT$external_gene_name %>%
#   web_gestalt_wrapper(df = ., project = "MvT")
# 
# MvC_pathways <-
#   edgeR_test_results_MvC$external_gene_name %>%
#   web_gestalt_wrapper(df = ., project = "MvC")
# 
# TvC_pathways <-
#   edgeR_test_results_MvC$external_gene_name %>%
#   web_gestalt_wrapper(df = ., project = "TvC")






# uncomment out above to run webgestalt,
# if already run, use below to load and format results table #2
# ==============================================================================
GO_MvC <- read_tsv("./Output/Project_MvC/enrichment_results_MvC.txt")
GO_MvT <- read_tsv("./Output/Project_MvT/enrichment_results_MvT.txt")
GO_TvC <- read_tsv("./Output/Project_TvC/enrichment_results_TvC.txt")

format_WebGest <- function(table) {
  table %>%
    transmute(`Gene Set` = geneSet,
              Description = description,
              Size = size,
              Overlap = overlap,
              Expected = expect,
              `Enrichment Ratio` = enrichmentRatio,
              `P-Value` = pValue,
              FDR = FDR,
              `Genes Overlapping` = userId)
}

write_xlsx(
  list(
    MvC = format_WebGest(GO_MvC),
    MvT = format_WebGest(GO_MvT),
    TvC = format_WebGest(GO_TvC)
    ),
  path = "./FiguresTables/ResultsTable2.xlsx",
  format_headers = F
)

system('open "./FiguresTables/ResultsTable2.xlsx"')


# comparing overlap between significant GO sets
# (supplemental file figure, ultimately created in external software)
# ------------------------------------------------------------------------------

euler(
  list(
    MvC = GO_MvC$geneSet,
    MvT = GO_MvT$geneSet,
    TvC = GO_TvC$geneSet),
  shape = "ellipse") %>%
  plot(labels = T, quantities = T)




# generating (former) Table 2 for main text
# genes were selected using this method for RT-PCR -- searching GO terms
# ------------------------------------------------------------------------------
WebGestalt_results <-
  inner_join(GO_MvC, GO_MvT, by = c("geneSet", "description"), suffix = c(".MvC", ".MvT")) %>%
  inner_join(., GO_TvC, by =  c("geneSet", "description"), suffix = c("", ".TvC"))


search_GO_results <- function(query_string, fdr_cutoff = 0.1) {
  GO_MvC %>% filter(grepl(query_string, description, ignore.case = T) & FDR < fdr_cutoff) %>%
    dplyr::select(description, pValue, FDR, userId) %>% print
  GO_MvT %>% filter(grepl(query_string, description, ignore.case = T) & FDR < fdr_cutoff) %>% 
    dplyr::select(description, pValue, FDR, userId) %>% print
  GO_TvC %>% filter(grepl(query_string, description, ignore.case = T) & FDR < fdr_cutoff) %>%
    dplyr::select(description, pValue, FDR, userId) %>% print
}



print_direction_gene <- function(query_string) {
  logFC <-
    c(edgeR_test_results_MvC_all %>%
        filter(external_gene_name == query_string) %>% .$logFC,
    edgeR_test_results_MvT_all %>%
      filter(external_gene_name == query_string) %>% .$logFC,
    edgeR_test_results_TvC_all %>%
      filter(external_gene_name == query_string) %>% .$logFC)
  
  signif <-
    c(edgeR_test_results_MvC_all %>%
        filter(external_gene_name == query_string) %>% .$fdr,
      edgeR_test_results_MvT_all %>%
        filter(external_gene_name == query_string) %>% .$fdr,
      edgeR_test_results_TvC_all %>%
        filter(external_gene_name == query_string) %>% .$fdr) %>%
    `<=`(., 0.05)
  
  output <-
    tibble(logFC = logFC,
         signif = signif) %>%
    mutate(
      direction =
      case_when(logFC < 0 & signif ~ "(↓)",
                logFC > 0 & signif ~ "(↑)",
                T ~ ""))
  
  output %>%
    transmute(logFC %>% formatC(format = "f", digits = 2) %>% paste(direction)) %>%
    unlist %>%
    set_names(., c("MvC", "MvT", "TvC"))
    
}


# search apotosis GO terms, across all three comparisons
search_GO_results("apop|hypox")
search_GO_results("lipopolysaccharide")

Table2 <-
  bind_rows(
  c("CASP1", "IRF5", "TLR6", "TREM2", "IL10RA") %>% sapply(., print_direction_gene) %>% t %>% as_tibble(rownames = "Gene") %>% mutate(Theme = "2_Response to Bacteria"),
  c("ARNT", "HIF1A", "VEGFB", "SIRT1", "CEBPB") %>% sapply(., print_direction_gene) %>% t %>%  as_tibble(rownames = "Gene") %>% mutate(Theme = "1_Apoptosis/Hypoxia"),
  c("TGFA", "TGFB2", "SMAD4", "FGF10", "MMP2") %>% sapply(., print_direction_gene) %>% t %>%  as_tibble(rownames = "Gene") %>% mutate(Theme = "1_Tissue Remodeling"),
  c("CYP1B1", "FOXO3", "ADRB2", "MYB", "MMP9") %>% sapply(., print_direction_gene) %>% t %>%  as_tibble(rownames = "Gene") %>% mutate(Theme = "1_Oxidative Stress"),
  c("ARID5B", "PPARD", "IL4R", "WNT5A", "CELSR1") %>% sapply(., print_direction_gene) %>% t %>% as_tibble(rownames = "Gene") %>% mutate(Theme = "2_Immune Cell Differentiation"),
  c("IFNG", "IRAK3", "SEMA7A", "STAT6", "HSPB1") %>% sapply(., print_direction_gene) %>% t %>% as_tibble(rownames = "Gene") %>% mutate(Theme = "2_Cytokine Production/Response")#,
) %>%
  dplyr::select(5, 1:4) %>%
  arrange(Theme, Gene) %>%
  mutate(Theme = gsub("1_|2_", "", Theme))
Table2$Theme[duplicated(Table2$Theme)] <- ""
View(Table2)


# check some genes manually to make sure print_direction_gene() working properly
edgeR_test_results_MvC %>% filter(grepl("TLR|CASP", external_gene_name))
edgeR_test_results_MvT %>% filter(grepl("TLR|CASP", external_gene_name))
edgeR_test_results_TvC %>% filter(grepl("TLR|CASP", external_gene_name))









