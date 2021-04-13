# ==============================================================================
# 07_table1_and_exposures.R
# Table 1: subject features, by group
# also: view tobacco/marijuana smoke exposure
# ==============================================================================



# "Table 1" (actually in online suppl. material Table A)
# clinical/demographic info, stratified by smoke status
# ------------------------------------------------------------------------------

Table1 <-
  metadata_filtered %>%
  bind_rows(., metadata_filtered %>% mutate(Group = "All")) %>%
  group_by(Group) %>%
  summarize(
    
    `00* N` = n(),
    
    `11* Demographic Information` = summaryBlankLine(),
    `12* Sex (Female)` = summaryCountPercent(Sex, "Female", digits = 1),
    `13* Age (years)` = summaryMedianIQR(Age),
    `14* BMI` = summaryMedianIQR(BMI),
    `15* Race/Ethnicity (White)` = summaryCountPercent(RaceEth_White, "Checked", digits = 1),
    
    # note: other cell types all <1%
    `20* % Differential Cell Types` = summaryBlankLine(),
    `21* ... Monocyte/Macrophage` = summaryMedianIQR(CellType_mono_mac, na.rm = T),
    `22* ... Lymphocyte` = summaryMedianIQR(CellType_lymph, na.rm = T),
    `23* ... PMNs` = summaryMedianIQR(CellType_pmns, na.rm = T),
    `24* ... Airway` = summaryMedianIQR(CellType_airway, na.rm = T),
    `25* ... Eosinophils` = summaryMedianIQR(CellType_eos, na.rm = T),

    ) %>%
  gather(`Smoking Status`, values, 2:ncol(.)) %>%
  spread(Group, values) %>%
  mutate(`Smoking Status` = gsub("[0-9][0-9]\\*", "", `Smoking Status`)) %>%
  set_colnames(c("Smoking Status", "All",
                 "Marijuana (M)", "Non-Smoker (C)", "Tobacco (T)")) %>%
  dplyr::select(1, 3:5, 2)


# statistical tests for table
# ------------------------------------------------------------------------------
lm(Age ~ Group, metadata_filtered) %>% anova
fisher.test(metadata_filtered$Sex, metadata_filtered$Group)
lm(BMI ~ Group, metadata_filtered) %>% anova
fisher.test(metadata_filtered$RaceEth_White, metadata_filtered$Group)

lm(log10(CellType_mono_mac + 1e-4) ~ Group, metadata_filtered) %>% anova
lm(log10(CellType_lymph + 1e-4) ~ Group, metadata_filtered) %>% anova
lm(log10(CellType_pmns + 1e-4) ~ Group, metadata_filtered) %>% anova
lm(log10(CellType_airway + 1e-4) ~ Group, metadata_filtered) %>% anova
lm(log10(CellType_eos + 1e-4) ~ Group, metadata_filtered) %>% anova








# degree of MJ & Tobacco smoke exposure
# ------------------------------------------------------------------------------


# tobacco exposure
exposures_T <-
  metadata_filtered %>%
  filter(Group == "T") %>%
  summarize(
    `10* Tobacco Consumption` = summaryBlankLine(),
    `11* ... "How many years have you smoked for?"` = summaryMedianIQR(Tobacco_YearsSmoked),
    `12* ... "How many packs per day do you smoke?"` = summaryMedianIQR(Tobacco_PacksDay, na.rm = T),
    `13* ... Pack-Years*` = summaryMedianIQR(Tobacco_PacksDay*Tobacco_YearsSmoked, na.rm = T)
  ) %>%
  gather(` `, `   `) %>%
  mutate(` ` = gsub("[0-9][0-9]\\*", "", ` `))

# marijuana exposure
exposures_M <-
  metadata_filtered %>%
  filter(Group == "M") %>%
  summarize(
    
    `11* "How do you use marijuana/cannabis?"` = summaryBlankLine(),
    `11* ... ☐ Smoking` = summaryCountPercent(HowUseMJ_Smoking, "Checked"),
    `12* ... ☐ Edibles` = summaryCountPercent(HowUseMJ_Edibles, "Checked"),
    `13* ... ☐ Vapor` = summaryCountPercent(HowUseMJ_Vapor, "Checked"),
    `14* ... ☐ Tincture` = summaryCountPercent(HowUseMJ_Tincture, "Checked"),
    
    `20* Duration Consuming Marijuana` = summaryBlankLine(),
    `21* ... "How old were you when you FIRST started SMOKING marijuana/cannabis?"` = summaryMedianIQR(MJSmoke_AgeStart, na.rm = T),
    `22* ... "How many YEARS have you been SMOKING marijuana?"` = summaryMedianIQR(MJSmoke_Years, na.rm = T),
    
    `30* Recent Use (Last Thirty Days)` = summaryBlankLine(),
    `31* ... "In these past 30 days, how many days did you SMOKE marijuana?"` = summaryMedianIQR(MJSmoke_LastMonth_NumDays, na.rm = T),
    `32* ... "On these days, how many times did you smoke marijuana PER DAY?"` = summaryMedianIQR(MJSmoke_LastMonth_TimesPerDay_edit, na.rm = T),
    `33* ... "On these days, how many grams of marijuana did you smoke PER DAY?"` = summaryMedianIQR(MJSmoke_LastMonth_GramsPerDay_edit, na.rm = T)
    
  ) %>%
  gather(` `, `   `) %>%
  mutate(` ` = gsub("[0-9][0-9]\\*", "", ` `))


sum(metadata_filtered$HowUseMJ_Edibles == "Checked" | metadata_filtered$HowUseMJ_Vapor == "Checked", na.rm = T)

exposure_summary <- rbind(exposures_T, exposures_M)
