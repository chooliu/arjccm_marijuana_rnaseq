# ==============================================================================
# 20_export_for_GEO.R
# share .fastqs, count table (note: pre-RUV correction), participant metadata
# NCBI GEO accession GSE155213
# ==============================================================================


# fastq checksums
# ------------------------------------------------------------------------------
# # check md5s on a Windows 10 system (where our raw .fastq.gz files were stored)
# Get-FileHash -Algorithm MD5 -Path (Get-ChildItem "*.fastq.gz" -Recurse) | Out-File -FilePath checksum_Run1.txt
# Get-FileHash -Algorithm MD5 -Path (Get-ChildItem "*.fastq.gz" -Recurse) | Out-File -FilePath checksum_Run2.txt

checksums <-
  read_excel("./data/Checksums/checksums.xlsx", col_names = c("fastq", "hash")) %>%
  separate(fastq, into = c("Run", "fastq"), sep = "\\\\") %>%
  transmute(`file name` = fastq,
            `file type` = "fastq",
            `file checksum` = hash,
            `instrument model` = if_else(Run == "Run1", "Illumina HiSeq4000", "Illumina HiSeq2500"),
            `single or paired-end` = "single") %>%
  arrange(`file name`)

checksums_to_join <-
  checksums %>%
  mutate(`Sample name` = str_split(`file name`, "_") %>% map_chr(~ .[1]) %>% gsub("A|B", "", .) %>% gsub("MF", "MJ", .))  %>%
  pivot_wider(id_cols = `Sample name`, values_from = `file name`, names_prefix = "fastqs", values_fn = list ) %>%
  transmute(`Sample name`,
            `raw file 1` = fastqs %>% map_chr( ~ .[1]),
            `raw file 2` = fastqs %>% map_chr( ~ .[2]))


# counts
# ------------------------------------------------------------------------------
edgeR_obj$counts %>%
  set_colnames(metadata_filtered$ID) %>% 
  as_tibble(rownames = "Features") %>%
  write_tsv("./FiguresTables/count_table.txt")

system('md5 FiguresTables/count_table.txt')



# participant-level metadata
# (copy and pasted into GEO submission spreadsheets - Excel)
# ------------------------------------------------------------------------------
# edgeR model:
# model.matrix(~ Group + Age + Sex + Obese + Batch + ruv_correction$W)
geo_metadata <-
  metadata_filtered %>%
  transmute(`Sample name` = ID,
            title = ID,
            `source name` = "bronchoalveolar lavage, alveolar macrophage",
            organism = "Homo sapiens",
            
            # demographics
            `characteristics: group` = Group,
            `characteristics: age` = Age,
            `characteristics: sex` = Sex,
            `characteristics: bmi` = BMI,
            `characteristics: obese` = Obese %>% if_else(., "yes", "no"),

            # batch
            `characteristics: technical_batch` = Batch,            
            
            # cell-type
            `characteristics: celltype_monocyte_macrophage` = CellType_mono_mac,
            `characteristics: celltype_lymphocytes` = CellType_lymph,
            `characteristics: celltype_PMNs` = CellType_pmns,
            `characteristics: celltype_airwaycells` = CellType_airway,
            
            # spirometry
            `characteristics: spirom_FEV1_percentpred` = FEV1_pred,
            `characteristics: spirom_FVC_percentpreed` = FVC_pred,
            `characteristics: spirom_FEVFVCratio` = FEVFVC,
            
            # smoking info
            `characteristics: mjconsumption_smoking` = if_else(Group == "M" & HowUseMJ_Smoking == "Checked", "yes", "no"),
            `characteristics: mjconsumption_edibles` = if_else(Group == "M" & HowUseMJ_Edibles == "Checked", "yes", "no"),
            `characteristics: mjconsumption_vapor`  = if_else(Group == "M" & HowUseMJ_Vapor == "Checked", "yes", "no"),
            `characteristics: mjconsumption_tincture` = if_else(Group == "M" & HowUseMJ_Tincture == "Checked", "yes", "no"),
            
            `characteristics: mjsmoke_agestart` = MJSmoke_AgeStart,
            `characteristics: mjsmoke_yearssmoked` = MJSmoke_Years,
            `characteristics: mjsmoke_jointyears` = MJSmoke_JointYears,
            `characteristics: mjsmoke_lastmonthnumdayssmoked` = MJSmoke_LastMonth_NumDays,
            `characteristics: mjsmoke_lastmonthgramsperday` = MJSmoke_LastMonth_GramsPerDay_edit,
            `characteristics: mjsmoke_lastmonthtimessmokedperday` = MJSmoke_LastMonth_TimesPerDay_edit,
            
            `characteristics: tobaccosmoke_agestart` = Tobacco_AgeStartSmoke,
            `characteristics: tobaccosmoke_yearssmoked` = Tobacco_YearsSmoked,
            `characteristics: tobaccosmoke_currentpacksperday` = Tobacco_PacksDay,
            `characteristics: tobaccosmoke_packyears` = Tobacco_PackYears,

            `characteristics: ruv1` = ruv_correction$W[ , 1],
            `characteristics: ruv2` = ruv_correction$W[ , 2],
            `characteristics: ruv3` = ruv_correction$W[ , 3],
            
            # geo requests
            molecule = if_else(Batch == "1", "polyA RNA", "total RNA"),
            `processed data file` = ""
            ) %>%
  left_join(checksums_to_join, by = "Sample name")


# metadata for Github, remove "characteristics: " label for GEO
write_tsv(geo_metadata %>%
            set_names(., gsub("characteristics: ", "", names(.))),
            "./FiguresTables/subject_metadata.txt")
