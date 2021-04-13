# ==============================================================================
# 01_index.R
# primary script to source R component of analysis
# ==============================================================================

# file structure
sapply(c("Output", "data", "data/Output", "FiguresTables", "Scripts"), dir.create)

# metadata, kallisto pseudocounts loading
source("./Scripts/02_load_libraries_and_customfxns.R")
source("./Scripts/03_load_metadata.R")
source("./Scripts/04_sequencing_metrics.R")
source("./Scripts/05_download_biomart.R")
source("./Scripts/06_kallisto_tximport.R")
source("./Scripts/07_table1_and_exposures.R")

# RUV + DESeq2
source("./Scripts/08_prefiltering.R")
source("./Scripts/09_examine_batcheffect_and_RUV.R") 
source("./Scripts/10_run_RUV.R")
source("./Scripts/11_PCA_by_batch.R")
source("./Scripts/12_edgeR.R")
source("./Scripts/13_plots_onresids_by_group.R")
source("./Scripts/14_venndiagram.R")
source("./Scripts/15_overenrichment_analysis.R")

# validation / sensitivity analyses
source("./Scripts/16_TvC_genelists.R")
source("./Scripts/17_rtPCR_validation.R")
source("./Scripts/18_celltype_sensitivity_pt1.R")
source("./Scripts/19_celltype_sensitivity_pt2.R")

# metadata export for GEO
source("./Scripts/20_export_for_GEO.R")
