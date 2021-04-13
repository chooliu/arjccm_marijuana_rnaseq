# ==============================================================================
# 08_prefiltering.R
# filter genes & subjects
# ==============================================================================

# keep gene if it is present at least 1/3 of subjects in each batch
# sum(genes_to_keep) = yields 17,602 genes
genes_to_keep <-
  apply(edgeR_obj$counts, 1,
        function(x) {
          ( (metadata_filtered$Batch == 1) %>% x[.] %>% `<=`(., 0) %>% sum(.) %>% `<=`(10) ) &
          ( (metadata_filtered$Batch == 2) %>% x[.] %>% `<=`(., 0) %>% sum(.) %>% `<=`(3) )
        }
        )

# keep all subjects;
# alternative example, omit those w/o cell differentials
# subjects_to_keep <- !is.na(metadata_filtered$CellType_lymph)
filter_subjects_to_keep <- 1:41

# do filtering
edgeR_obj <- edgeR_obj[genes_to_keep, filter_subjects_to_keep]

