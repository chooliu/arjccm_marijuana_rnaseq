# ==============================================================================
# 06_kallisto_tximport.R
# transcripts --> gene via tximport  --> save as "tximport_counts.Rdata"
# ==============================================================================




# prep for data import
# ------------------------------------------------------------------------------
kal_dirs <- file.path("./Data/kallisto_output", kallisto_ids)



# # load as gene-wise counts
# # ------------------------------------------------------------------------------
# 
# txi <-
#   tximport(
#     file.path(kal_dirs, "abundance.h5") %>%
#       .[gsub("MF107", "MJ107", .) %>% order(.)],
#     type = "kallisto",
#     tx2gene = annotations_by_transcript[ , c(1, 3)],
#     txIdCol = "ensembl_transcript_id",
#     geneIdCol = "ensembl_gene_id",
#     txOut = F,
#     ignoreTxVersion = T
#   )
# 
# # 28775
# txi$counts %>% `==`(., 0) %>% apply(., 1, all) %>% `!` %>% sum
# 
# 
# # double check sample order
# kallisto_ids %>% gsub("MF107", "MJ107", .) %>% .[order(.)] %>%
#   gsub("A|B", "", .) %>%
#   identical(., metadata_filtered$ID)
# 
# # 40320 x 41
# txi$counts %>% dim()
# 
# 
# 
# 
# 
# # count number of non-zero transcripts
# # ------------------------------------------------------------------------------
# txi_transcript <-
#   tximport(
#     file.path(kal_dirs, "abundance.h5") %>%
#       .[gsub("MF107", "MJ107", .) %>% order(.)],
#     type = "kallisto",
#     tx2gene = annotations_by_transcript[ , c(1, 3)],
#     txIdCol = "ensembl_transcript_id",
#     geneIdCol = "ensembl_gene_id",
#     txOut = T,
#     ignoreTxVersion = T
#   )
# 
# # 155188
# txi_transcript$counts %>% `==`(., 0) %>% apply(., 1, all) %>% `!` %>% sum
# 
# 
# 
# 
# 
# # export gene matrix in edgeR (instructions via package vignette)
# # ------------------------------------------------------------------------------
# 
# 
# cts <- txi$counts
# normMat <- txi$length
# normMat <- normMat/exp(rowMeans(log(normMat)))
# 
# o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
# edgeR_obj <- DGEList(cts)
# edgeR_obj <- scaleOffset(edgeR_obj, t(t(log(normMat)) + o))
# 
# 
# save(edgeR_obj, file = "./Data/Output/tximport_counts.Rdata")




# save primary output of above (commented out)
# ------------------------------------------------------------------------------
load("./Data/Output/tximport_counts.Rdata")

