# ==============================================================================
# 05_download_biomart.R
# retreive Ensembl/biomaRt annotations --> save as "biomart.Rdata"
# ==============================================================================



# # Ensembl GRCh38 v96, retrieved April 7th 2018
# # (for consistency w/ kallisto quants)
# # ------------------------------------------------------------------------------
# 
# biomart_obj <-
#   useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#           dataset = "hsapiens_gene_ensembl",
#           host = "apr2019.archive.ensembl.org")
# annotations_by_transcript <- getBM(
#   attributes = c("ensembl_transcript_id", "transcript_version",
#                  "ensembl_gene_id", "external_gene_name", "description",
#                  "transcript_biotype",
#                  "start_position", "end_position"),
#   mart = biomart_obj) %>%
#   as_tibble() %>%
#   mutate(description = description %>% stri_extract(regex = "[^\\[]*")) %>%
#   mutate(length = end_position - start_position)
# 
# 
# annotations_by_gene <-
#   annotations_by_transcript %>%
#   filter(!duplicated(ensembl_gene_id)) %>%
#   dplyr::select(
#     ensembl_gene_id, external_gene_name, description, length)
# 
# save(biomart_obj, annotations_by_transcript, annotations_by_gene,
#      file = "./Data/Annotations/biomart.Rdata")





# save primary output of above (commented out)
# ------------------------------------------------------------------------------
load("./Data/Annotations/biomart.Rdata")
