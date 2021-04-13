# ==============================================================================
# 03_load_metadata.R
# extract patient metadata into tidy format:
# - loads two separate REDCap files for MJ smokers & controls (tobacco/NS),
# - separate functional/assay data files, etc
# - some manual editing due to database entry issues
# [note to external reviewers: unfortunately, cannot share all source data files]
# ==============================================================================




# REDCap
# ==============================================================================

# MJ dataset (dat1)
# ------------------------------------------------------------------------------

dat1 <-
  read_excel(
    path = "./data/redcap/RedCap MJ and Controls 080919.xlsx") %>%
  dplyr::slice(-1) %>%  # empty row
  set_names(., make.names(names(.)))

# note: two warning msgs due to subjects' responses reported as ranges
# (rather than numerical values)-- e.g., "2-3 times/day" -- to fix after import

# Warning messages:
#   1: In read_fun(path = enc2native(normalizePath(path)), sheet_i = sheet,  :
#      Expecting numeric in DR18 / R18C122: got a date
#   2: In read_fun(path = enc2native(normalizePath(path)), sheet_i = sheet,  :
#      Expecting numeric in EJ29 / R29C140: got a date





# non-smokers & tobacco smokers (dat2)
# ------------------------------------------------------------------------------

dat2 <-
  read_excel(
    path = "./data/redcap/RedCap AUD and Controls 080919.xlsx") %>%
  dplyr::slice(-1) %>%  # empty row
  set_names(., make.names(names(.)))

# Warning messages:
#   1: In read_fun(path = enc2native(normalizePath(path)), sheet_i = sheet,  :
#      Expecting numeric in DR18 / R18C122: got a date
#   2: In read_fun(path = enc2native(normalizePath(path)), sheet_i = sheet,  :
#      Expecting numeric in EJ29 / R29C140: got a date






# extract demographics & smoking frequency
# (normalize variable names in dat1 and dat2, then join)
# ------------------------------------------------------------------------------

dat1_clean <-
  dat1 %>%
  transmute(
    
    # demographics
    ID = Record.ID,
    Age = What.your.age,
    Sex = What.is.your.gender.,
    Height_in = Height..inches.,
    Weight_lbs = Weight..lbs.,
    RaceEth_AmerIndAlasNative = What.is.your.race...choice.American.Indian.Native.Alaskan.,
    RaceEth_Asian = What.is.your.race...choice.Asian.,
    RaceEth_Black = What.is.your.race...choice.Black.African.American.,
    RaceEth_NativeHawaiPacIsl = What.is.your.race...choice.Native.Hawaiian.Other.Pacific.Islander.,
    RaceEth_DontKnow = What.is.your.race...choice.Dont.Know.,
    RaceEth_Refused = What.is.your.race...choice.Refuse.to.Answer.,
    RaceEth_White = What.is.your.race...choice.White.,
    RaceEth_Hispanic = Are.you.Hispanic.or.Latino.,
    
    # marijuana consumption method
    HowUseMJ_Smoking = How.do.you.use.marijuana.cannabis...choice.Smoking.,
    HowUseMJ_Edibles = How.do.you.use.marijuana.cannabis...choice.Edibles.,
    HowUseMJ_Vapor = How.do.you.use.marijuana.cannabis...choice.Vapor.,
    HowUseMJ_Tincture = How.do.you.use.marijuana.cannabis...choice.Tincture.,
    HowUseMJ_NA = How.do.you.use.marijuana.cannabis...choice.Not.Applicable.,
    
    # for smoking, vaping, and edible
    # main readouts on long-term frequency & last month
    MJSmoke_AgeStart = How.old.were.you.when.you.FIRST.started.SMOKING.marijuana.cannabis.,
    MJSmoke_Years = How.many.YEARS.have.you.been.SMOKING.marijuana.,
    MJSmoke_LastMonth_NumDays = Think.specifically.about.the.PAST.30.DAYS..In.these.past.30.days..how.many.days.did.you.SMOKE.marijuana.,
    MJSmoke_LastMonth_TimesPerDay = On.these.days.that.you.smoked.marijuana..how.many.times.did.you.smoke.marijuana.PER.DAY..frequency..,
    MJSmoke_LastMonth_GramsPerDay = On.these.days..how.much.marijuana..in.grams..did.you.smoke.PER.DAY.,
    
    MJVape_AgeStart = How.old.were.you.when.you.first.started.using.vaporized..vaping..marijuana.,
    MJVape_Years = How.many.years.have.you.been.vaporzing..vaping..marijuana.,
    MJVape_LastMonth_NumDays = Now..think.specifically.about.the.PAST.30.DAYS..In.the.past.30.days..on.how.many.days.did.you.vape.marijuana.,
    MJVape_LastMonth_TimesPerDay = On.these.days.that.you.vaped.marijuana..how.many.times.did.you.vape.marijuana.PER.DAY.,
    MJVape_LastMonth_GramsPerDay = On.these.days.that.you.vaped.marijuana..how.much.marijuana.did.you.vape..in.grams..,
    
    MJEdibles_AgeStart = How.old.were.you.when.you.FIRST.started.EATING.CONSUMING.marijuana.products.edibles.,
    MJEdibles_Years = How.many.years.have.you.been.EATING.or.consuming.marijuana.and.marijuana.EDIBLES..in.years..,
    MJEdibles_LastMonth_NumDays = Think.specifically.about.the.PAST.30.DAYS..In.the.past.30.days..how.many.days.did.you.eat.or.consume.marijuana.or.marijuana.edibles.,
    MJEdibles_LastMonth_TimesPerDay = On.these.days.that.you.consumed.edible.marijuana..how.many.times.did.you.eat.or.consume.marijuana.PER.DAY.,
    MJEdibles_LastMonth_GramsPerDay = On.these.days.that.you.consumed.edible.marijuana..how.much.marijuana..in.grams..did.you.eat.PER.DAY.)



dat2_clean <-
  dat2 %>%
  transmute(
    
    # demographics
    ID = Subject.ID,
    Age = X1..Age,
    Sex = X2..Gender.at.Birth,
    Height_in = X3a...Height.in.inches.,
    Weight_lbs = X3b..Weight..in.lbs.,
    RaceEth_AmerIndAlasNative = X2..Race..select.one.or.more...choice.American.Indian.or.Alaska.Native.,
    RaceEth_Asian = X2..Race..select.one.or.more...choice.Asian.,
    RaceEth_Black = X2..Race..select.one.or.more...choice.Black.or.African.American.,
    RaceEth_NativeHawaiPacIsl = X2..Race..select.one.or.more...choice.Native.Hawaiian.or.Other.Pacific.Islander.,
    RaceEth_DontKnow = X2..Race..select.one.or.more...choice.Dont.Know.,
    RaceEth_Refused = X2..Race..select.one.or.more...choice.Refused.,
    RaceEth_White = X2..Race..select.one.or.more...choice.White.,
    RaceEth_Hispanic = X1..Hispanic.Latino,
    
    # tobacco
    Tobacco_YearsSmoked = X8b..How.many.years.have.you.smoked.for.,
    Tobacco_PacksDay = X8a..How.many.packs.per.day.do.you.smoke.,
    Tobacco_AgeStartSmoke = Age - Tobacco_YearsSmoked,
    Tobacco_CurrentSmoke = X8..Do.you.presently.smoke.cigarettes.

  )







# create "metadata" tibble,
# standarized some variables across datasets,
# add some new variables of interest (e.g., obesity)
# ==============================================================================

# merge mj && control/smokers
metadata <-
  bind_rows(dat1_clean, dat2_clean) %>%
  arrange(ID)

metadata %<>%
  mutate(
    BMI = Weight_lbs*0.454 / (Height_in * 0.0254)^2,
    Obese = BMI >= 30,
    ID = toupper(ID),
    Group = substr(ID, 1, 1) %>% as.factor %>% fct_recode(., C = "N", M = "M", T = "S"),
    Group_Detailed = fct_recode(Group, `Marijuana (M)` = "M", `Tobacco (T)` = "T", `Control (C)` = "C"))

# filter to just M, C, T
metadata %<>%
  filter(Group %in% c("M", "C", "T")) %>%
  mutate(Group = fct_drop(Group),
         Group_Detailed = fct_drop(Group_Detailed))

# fixing a few entry issues in consumption frequency (see warning messages on line 17) &
# note: conversion of qualitative values like "8+" --> numbers limits use of means (use medians)
metadata %<>%
  mutate(
    MJSmoke_LastMonth_TimesPerDay_edit =
      if_else(MJSmoke_LastMonth_TimesPerDay == "8 or more", "8", MJSmoke_LastMonth_TimesPerDay) %>% as.numeric,
    MJVape_LastMonth_TimesPerDay_edit =
      if_else(MJVape_LastMonth_TimesPerDay == "8 or more", "8", MJVape_LastMonth_TimesPerDay) %>% as.numeric,
    MJVape_LastMonth_TimesPerDay_edit = if_else(ID == "MJ126", 1.5, MJVape_LastMonth_TimesPerDay_edit),
    MJSmoke_LastMonth_GramsPerDay_edit = gsub(" grams| or greater", "", MJSmoke_LastMonth_GramsPerDay) %>% as.numeric,
    MJVape_LastMonth_GramsPerDay_edit = case_when(MJVape_LastMonth_GramsPerDay == "More than 80 mg" ~ 80,
                                                  MJVape_LastMonth_GramsPerDay == "21 mg to 30 mg" ~ 25.5,
                                                  MJVape_LastMonth_GramsPerDay == "11 mg to 20 mg" ~ 15.5,
                                                  MJVape_LastMonth_GramsPerDay == "0" ~ 0),
    MJSmoke_AgeStart =
           if_else(is.na(MJSmoke_AgeStart), Age - MJSmoke_Years, MJSmoke_AgeStart))


# manually change subjects who do not vape from vape rates of NA/blank --> 0
metadata %<>%
  mutate(
    MJVape_LastMonth_TimesPerDay = ifelse(HowUseMJ_Vapor == "Unchecked" & Group == "M", 0, MJVape_LastMonth_TimesPerDay),
    MJVape_LastMonth_NumDays = ifelse(HowUseMJ_Vapor == "Unchecked" & Group == "M", 0, MJVape_LastMonth_NumDays),
    MJVape_LastMonth_GramsPerDay = ifelse(HowUseMJ_Vapor == "Unchecked" & Group == "M", 0, MJVape_LastMonth_GramsPerDay))




# manual fixes for miscellaneous data entry errors
# ==============================================================================

# S116 duplicated in Redcap, filled out both questionnaires but height, weight different
# e-mail from J.G. 1/23/2018 that first S161 record is correct
metadata %<>% filter(!duplicated(ID))

# data entry mistakes
# e-mail from J.G., 10/6/2017
metadata %<>% mutate(Tobacco_YearsSmoked = if_else(ID == "S161", 10, Tobacco_YearsSmoked))

# MJ120 with 364 years smoked, J.G. correction in Word doc 2/12/2018
metadata %<>% mutate(MJSmoke_Years = if_else(ID == "MJ120", 37, MJSmoke_Years))




# load additional functional data not yet in REDCap
# ==============================================================================

# THC metabolite concentration and marijuana joint-years
# ------------------------------------------------------------------------------
metadata_jmp <-
  read_csv("./Data/redcap/Cannabis_quantification_data_121219.txt") %>%
  set_names(., c("ID", "MJSmoke_JointYears", "Urine_THC_COOH", "Urine_THC_gluc", "BAL_cell_THC"))
metadata <-
  left_join(metadata, metadata_jmp, by = "ID")

metadata %<>%
  mutate(Tobacco_PackYears = Tobacco_YearsSmoked * Tobacco_PacksDay)



# cytokine assay data, spirometry, and cell differentials (in three parts)
# ------------------------------------------------------------------------------

assay_data_A <-
  read_excel(path = "./Data/Subject_Metadata/RedCap Cytokines Apop and Phago Data 9-19-19.xlsx", skip = 4, n_max = 30) %>%
  set_names(., names(.) %>% gsub("-|%|,", "", .) %>% gsub("\\/| ", "_", .) %>%  # remove punctuation
              gsub("Inex", "Index", .) %>% gsub("delta", "Delta", .)) %>% # standardize column names
  mutate(Patient_ID = toupper(Patient_ID)) %>% 
  dplyr::rename(ID = Patient_ID)

assay_data_B <-
  read_excel(path = "./Data/Subject_Metadata/RedCap Cytokines Apop and Phago Data 9-19-19.xlsx", skip = 34, n_max = 30) %>%
  set_names(., names(.) %>% gsub("-|%|,", "", .) %>% gsub("\\/| ", "_", .) %>%  # remove punctuation
              gsub("Inex", "Index", .) %>% gsub("delta", "Delta", .)) %>% # standardize column names
  mutate(Patient_ID = toupper(Patient_ID)) %>% 
  dplyr::rename(ID = Patient_ID)

assay_data_C <-
  read_excel(path = "./Data/Subject_Metadata/RedCap Cytokines Apop and Phago Data 9-19-19.xlsx", skip = 66, n_max = 30) %>%
  set_names(., names(.) %>% gsub("-|%|,", "", .) %>% gsub("\\/| ", "_", .) %>%  # remove punctuation
              gsub("Inex", "Index", .) %>% gsub("delta", "Delta", .)) %>% # standardize column names
  mutate(Patient_ID = toupper(Patient_ID)) %>% 
  dplyr::rename(ID = Patient_ID)


assay_data <-
  bind_rows(assay_data_A,
            assay_data_B,
            assay_data_C %>% mutate(Race = as.numeric(Race))
            )

rm(assay_data_A, assay_data_B, assay_data_C)

assay_data_celldiff <-
  assay_data %>%
  dplyr::select(ID, matches("^_")) %>%
  dplyr::select(-`_Recovery`) %>%
  set_names(., names(.) %>% paste0("CellType", .))



# manual cell diff additions, via
# J.G. email 1 Jun 2020
# ------------------------------------------------------------------------------
assay_data_celldiff <-
  bind_rows(
    assay_data_celldiff,
    tribble(~CellTypeID, ~CellType_mono_mac, ~CellType_lymph,
            ~CellType_pmns, ~CellType_eos, ~CellType_airway,
            "N122", 94.81, 4.72, 0.47, 0.00, NA,
            "S121", 93.66, 5.37, 0.49, 0.49, NA,
            "S124", 90.14, 7.04, 2.35, 0.47, NA,
            "S125", 86.45, 11.21, 2.34, 0.00, NA)
  )

assay_data_fev <-
  assay_data %>%
  dplyr::select(ID, contains("FEV"), contains("FVC")) %>%
  mutate(FEVFVC = if_else( (FEV1_FVC != FEV1_FVC_082019) & !is.na(FEV1_FVC_082019),
                           FEV1_FVC_082019, FEV1_FVC) ) %>%
  dplyr::select(ID, FEV1_pred, FVC_pred, FEVFVC)

metadata %<>%
  left_join(assay_data_celldiff, by = c("ID" = "CellTypeID"))
metadata %<>%
  left_join(assay_data_fev)








# filter metadata to include subjects sequenced only
# uses old metadata files (July 2017) to identify subject IDs in batch 1 (vs 2)
# ==============================================================================

kallisto_ids <-
  dir(file.path("./Data/kallisto_output/"))

metadata_filtered <- 
  metadata %>%
  filter(ID %in% gsub("A|B", "", c(kallisto_ids, "MJ107")))

batch1_names <-
  c(
    read_excel(
      path = "./data/redcap/RedCap Data for RNA-Seq N and S Patients with spirometry 7-11-17.xlsx",
      skip = 1) %>%
      .$`Subject ID`,
    read_excel(
      path = "./data/redcap/RedCap Data for RNA-Seq MJ N and S Patients with spirometry 7-11-17.xlsx",
      skip = 1) %>%
      .$`Record ID`
  )

metadata_filtered %<>%
  mutate(Batch = if_else(ID %in% batch1_names, 1, 2))
metadata_filtered$Batch


# 41 x 59
dim(metadata_filtered)

