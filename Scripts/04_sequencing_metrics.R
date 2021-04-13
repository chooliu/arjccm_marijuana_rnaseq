# ==============================================================================
# 04_sequencing_metrics.R
# examine sequencing core QC metrics; reads by batch & smoking group
# [note to external reviewers: unfortunately, cannot share all source data files]
# ==============================================================================

# first sequencing run slightly higher quality & read count
# but no obvious difference in seq quality by smoking group


# run 1
# ------------------------------------------------------------------------------
qcseq1 <-
  read_csv("./data/seq_core/161013_K00262_0049_BHGC52BBXX_L1234_Burnham_demux.csv",
        skip = 6) %>%
  dplyr::slice(-1) %>%
  set_names(., make.names(names(.))) %>%
  mutate(Group = as.factor(substring(Sample, 0, 1)),
         Clusters = Clusters %>% gsub(",", "", .) %>% as.numeric,
         Yield_Mb = `Yield..Mbases.` %>% gsub(",", "", .) %>% as.numeric) %>%
  bind_rows(., mutate(., Group = "All")) %>%
  group_by(Group) %>%
  summarise(`Million Clusters` = summaryMedianRange(as.numeric(Clusters/1e6)),
            `Yield (Mb)` = summaryMedianIQR(Yield_Mb),
            `Number Reads (M)` = summaryMedianIQR(Yield_Mb/150),
            `Mean QScore` = summaryMedianRange(as.numeric(Mean.Quality)),
            `% Bases >= Q30` = summaryMedianRange(as.numeric(X.....Q30)),
            `% Barcode Perfect` = summaryMedianRange(as.numeric(X..Perfect)),
            `% Barcode 1-Mismatch` = summaryMedianRange(as.numeric(X..One.mismatch)),
            `% of Lane` = summaryMedianRange(as.numeric(X..of.the))) %>%
  gather(key, value, -Group) %>%
  spread(Group, value)


# run 2
# ------------------------------------------------------------------------------
qcseq2 <-
  read_csv("./data/seq_core/171214_7001413_0377_BCA8MRANXX_L12_Burnham_demux.csv",
        skip = 4) %>%
  dplyr::slice(-c(12, 24)) %>%
  mutate(Group = as.factor(substring(`Sample ID`, 0, 1))) %>%
  bind_rows(., mutate(., Group = "All")) %>%
  group_by(Group) %>%
  mutate(`Yield (Mbases)` = `Yield (Mbases)` %>% gsub(",", "", .) %>% as.numeric()) %>%
  summarise(`Million Clusters` = summaryMedianRange(`% of raw clusters per lane`),
            `Yield (Mb)` = summaryMedianIQR(`Yield (Mbases)`),
            `Number Reads (M)` = summaryMedianIQR(`Yield (Mbases)`/150),
            `Mean QScore` = summaryMedianRange(`Mean Quality Score (PF)`),
            `% Bases >= Q30` = summaryMedianRange(as.numeric(`% of >= Q30 Bases (PF)`)),
            `% Barcode Perfect` = summaryMedianRange(as.numeric(`% Perfect Index Reads`)),
            `% Barcode 1-Mismatch` = summaryMedianRange(as.numeric(`% One Mismatch Reads (Index)`))) %>%
  gather(key, value, -Group) %>%
  spread(Group, value)


