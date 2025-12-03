# Import library
library(tidyverse)
library(zipcodeR)

###################################################################### I ### Load data----
path_raw <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "MDSNHSDisparities")

marrow1 <- 
  read_csv(paste0(path_raw,
                  "/RawData/Data_MDSNHS_Aug2025Freeze",
                  "/MDSNHS_marrow1_data_Gillis.csv")) %>% 
  janitor::clean_names()

mutation_data <- 
  readxl::read_xlsx(paste0(path_raw,
                           "/ProcessedData",
                           "/MDSdisparities_MutationsAnnotated_20251114.xlsx")) %>% 
  janitor::clean_names()

blood_count <- 
  read.csv(paste0(path_raw,
                  "/RawData/Data_MDSNHS_Aug2025Freeze",
                  "/MDSNHS_blood_chem_data_Gillis.csv"), encoding = "latin1") %>% 
  janitor::clean_names()

parent_dir_path <- dirname(here::here())
mds_data <- read_rds(paste0(parent_dir_path, "/mds_disparities",
                            "/data/processed data",
                            "/MDSNHS_AnalysisDataset_20251125.rds"))

path_project <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "MDSSerialAnalysis")


###################################################################### II ### Data cleaning----
# Make wide series mutation----
mutation_data1 <- mutation_data %>% 
  # include all mutations in bone marrow samples if available for a patient and 
  # add blood samples if no BM samples tested
  mutate(have_blood = case_when(
    tissue_type == "Blood"                            ~ "Yes"
  ), .after = tissue_type) %>% 
  mutate(have_bm = case_when(
    tissue_type == "Bone Marrow"                      ~ "Yes"
  ), .after = tissue_type) %>% 
  group_by(mdsepid, visit) %>%
  fill(have_blood, have_bm, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(tissue_type_filter = case_when(
    have_blood == "Yes" &
      have_bm == "Yes"                                ~ "need filter"
  ), .after = tissue_type) %>% 
  filter((tissue_type_filter == "need filter" & tissue_type == "Bone Marrow") |
           is.na(tissue_type_filter)) %>% 
  select(-c(have_blood, have_bm, tissue_type_filter)) %>% 
  # filter somatic level 1
  # code mutation level
  mutate(mutation_level = 
           str_match(variant_confidence_level, "(Level [:digit:]) (.*)$")[,2], 
         .after = variant_confidence_level) %>% 
  filter(somatic_germline == "somatic" & 
           mutation_level == "Level 1") %>% 
  # code sequencing visit time
  mutate(sequencing_visit_time_range_days = case_when(
    collection_dt <= (0 + 30) &
      collection_dt >= (0 - 30)                       ~ 0,
    collection_dt <= (180 + 30) &
      collection_dt >= (180 - 30)                     ~ 180,
    collection_dt <= (360 + 30) &
      collection_dt >= (360 - 30)                     ~ 360,
    collection_dt <= (540 + 30) &
      collection_dt >= (540 - 30)                     ~ 540,
    collection_dt <= (720 + 30) &
      collection_dt >= (720 - 30)                     ~ 720,
    collection_dt <= (900 + 30) &
      collection_dt >= (900 - 30)                     ~ 900,
    collection_dt <= (1080 + 30) &
      collection_dt >= (1080 - 30)                    ~ 1080,
    collection_dt <= (1260 + 30) &
      collection_dt >= (1260 - 30)                    ~ 1260,
    collection_dt <= (1440 + 30) &
      collection_dt >= (1440 - 30)                    ~ 1440,
    collection_dt <= (1800 + 30) &
      collection_dt >= (1800 - 30)                    ~ 1800,
    collection_dt <= (1980 + 30) &
      collection_dt >= (1980 - 30)                    ~ 1980,
    collection_dt <= (2160 + 30) &
      collection_dt >= (2160 - 30)                    ~ 2160,
    collection_dt <= (2340 + 30) &
      collection_dt >= (2340 - 30)                    ~ 2340,
    collection_dt <= (2520 + 30) &
      collection_dt >= (2520 - 30)                    ~ 2520
  )) %>% 
  # focus on these
  filter(!is.na(sequencing_visit_time_range_days)) %>% 
  # Create var for mutation in splicing gene
  mutate(slicing_mutation_for_serie = case_when(
    gene %in% c("SF3B1", "SRSF2", 
                "U2AF1", "ZRSR2")                     ~ "splicing"
  )) %>% 
  group_by(mdsepid, sequencing_visit_time_range_days) %>% 
  fill(slicing_mutation_for_serie, .direction = "updown") %>% 
  # Create var for number of mutations
  mutate(number_of_mutation_for_serie = n()) %>% 
  ungroup() %>% 
  # Make data 1 sample time per row
  select(mdsepid, collection_dt_for_serie = collection_dt, 
         sequencing_visit_time_range_days, 
         slicing_mutation_for_serie, number_of_mutation_for_serie) %>% 
  distinct(mdsepid, sequencing_visit_time_range_days, .keep_all = TRUE) %>% 
  arrange(sequencing_visit_time_range_days, mdsepid)

mutation_data2 <- mutation_data1 %>% 
  # Make wide data
  pivot_wider(id_cols = mdsepid, 
              names_from = c(sequencing_visit_time_range_days), 
              values_from = c(collection_dt_for_serie, slicing_mutation_for_serie, number_of_mutation_for_serie), 
              names_vary = "slowest") 


# Mutation other way ----
mutation_data1 <- mutation_data %>% 
  # include all mutations in bone marrow samples if available for a patient and 
  # add blood samples if no BM samples tested
  mutate(have_blood = case_when(
    tissue_type == "Blood"                            ~ "Yes"
  ), .after = tissue_type) %>% 
  mutate(have_bm = case_when(
    tissue_type == "Bone Marrow"                      ~ "Yes"
  ), .after = tissue_type) %>% 
  group_by(mdsepid, visit) %>%
  fill(have_blood, have_bm, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(tissue_type_filter = case_when(
    have_blood == "Yes" &
      have_bm == "Yes"                                ~ "need filter"
  ), .after = tissue_type) %>% 
  filter((tissue_type_filter == "need filter" & tissue_type == "Bone Marrow") |
           is.na(tissue_type_filter)) %>% 
  select(-c(have_blood, have_bm, tissue_type_filter)) %>% 
  # filter somatic level 1
  # code mutation level
  mutate(mutation_level = 
           str_match(variant_confidence_level, "(Level [:digit:]) (.*)$")[,2], 
         .after = variant_confidence_level) %>% 
  filter(somatic_germline == "somatic" & 
           mutation_level == "Level 1") %>% 
  # Create var for mutation in splicing gene
  mutate(slicing_mutation_for_serie = case_when(
    gene %in% c("SF3B1", "SRSF2", 
                "U2AF1", "ZRSR2")                     ~ "splicing"
  )) %>% 
  group_by(mdsepid, collection_dt) %>% 
  fill(slicing_mutation_for_serie, .direction = "updown") %>% 
  # Create var for number of mutations
  mutate(number_of_mutation_for_serie = n()) %>% 
  ungroup() %>% 
  # Make data 1 sample time per row
  select(mdsepid, collection_dt_for_serie = collection_dt, 
         slicing_mutation_for_serie, number_of_mutation_for_serie) %>% 
  distinct(mdsepid, collection_dt_for_serie, .keep_all = TRUE) %>% 
  arrange(mdsepid, collection_dt_for_serie) %>% 
  mutate(has_mutation = "Yes")


# Make wide series blood count - only for Hgb, Plt, ANC ----
# Luckily we have only 1 value within the range of each visit
long_blood_count <- blood_count %>% 
  `colnames<-`(str_remove(colnames(.), "_xx")) %>% 
  rename(series_lab_result = lab_rs,
         series_lab_unit = lab_unit) %>% 
  filter(lab_name %in% c("Absolute Neutrophil Count (ANC), Blood",
                         "Hemoglobin, Blood", "Platelets, Blood")) %>% 
  mutate(lab_name = case_when(
    str_detect(lab_name, "Absolute Neutrophil")       ~ "anc",
    str_detect(lab_name, "Hemoglobin")                ~ "hgb",
    str_detect(lab_name, "Platelets")                 ~ "plt"
  )) %>% 
  # code blood visit time - using the 14 potential visit in sequencing
  mutate(lab_visit_time_range_days = case_when(
    lab_dt <= (0 + 30) &
      lab_dt >= (0 - 30)                              ~ 0,
    lab_dt <= (180 + 30) &
      lab_dt >= (180 - 30)                            ~ 180,
    lab_dt <= (360 + 30) &
      lab_dt >= (360 - 30)                            ~ 360,
    lab_dt <= (540 + 30) &
      lab_dt >= (540 - 30)                            ~ 540,
    lab_dt <= (720 + 30) &
      lab_dt >= (720 - 30)                            ~ 720,
    lab_dt <= (900 + 30) &
      lab_dt >= (900 - 30)                            ~ 900,
    lab_dt <= (1080 + 30) &
      lab_dt >= (1080 - 30)                           ~ 1080,
    lab_dt <= (1260 + 30) &
      lab_dt >= (1260 - 30)                           ~ 1260,
    lab_dt <= (1440 + 30) &
      lab_dt >= (1440 - 30)                           ~ 1440,
    lab_dt <= (1800 + 30) &
      lab_dt >= (1800 - 30)                           ~ 1800,
    lab_dt <= (1980 + 30) &
      lab_dt >= (1980 - 30)                           ~ 1980,
    lab_dt <= (2160 + 30) &
      lab_dt >= (2160 - 30)                           ~ 2160,
    lab_dt <= (2340 + 30) &
      lab_dt >= (2340 - 30)                           ~ 2340,
    lab_dt <= (2520 + 30) &
      lab_dt >= (2520 - 30)                           ~ 2520
  )) %>% 
  filter(!is.na(lab_visit_time_range_days)) %>% 
  select(-lab_dt)

wide_blood_count <- 
  long_blood_count %>% 
  pivot_wider(id_cols = mdsepid, 
              names_from = c(lab_name, lab_visit_time_range_days), 
              values_from = c(series_lab_result, series_lab_unit), 
              names_vary = "slowest") 


# Blood the other way ----
long_blood_count <- blood_count %>% 
  `colnames<-`(str_remove(colnames(.), "_xx")) %>% 
  rename(series_lab_result = lab_rs,
         series_lab_unit = lab_unit) %>% 
  filter(lab_name %in% c("Absolute Neutrophil Count (ANC), Blood",
                         "Hemoglobin, Blood", "Platelets, Blood")) %>% 
  mutate(lab_name = case_when(
    str_detect(lab_name, "Absolute Neutrophil")       ~ "anc",
    str_detect(lab_name, "Hemoglobin")                ~ "hgb",
    str_detect(lab_name, "Platelets")                 ~ "plt"
  ))
wide_blood_data <- 
  long_blood_count %>% 
  pivot_wider(id_cols = c(mdsepid, lab_dt),
              names_from = c(lab_name), 
              values_from = c(series_lab_result, series_lab_unit), 
              names_vary = "slowest") 


# Marrow----
marrow_1 <-
  marrow1 %>% 
  select(mdsepid, bm_bx_dt, bm_blast_pct = blast_pct) %>% 
  # code sequencing visit time
  mutate(blast_visit_time_range_days = case_when(
    bm_bx_dt <= (0 + 30) &
      bm_bx_dt >= (0 - 30)                            ~ 0,
    bm_bx_dt <= (180 + 30) &
      bm_bx_dt >= (180 - 30)                          ~ 180,
    bm_bx_dt <= (360 + 30) &
      bm_bx_dt >= (360 - 30)                          ~ 360,
    bm_bx_dt <= (540 + 30) &
      bm_bx_dt >= (540 - 30)                          ~ 540,
    bm_bx_dt <= (720 + 30) &
      bm_bx_dt >= (720 - 30)                          ~ 720,
    bm_bx_dt <= (900 + 30) &
      bm_bx_dt >= (900 - 30)                          ~ 900,
    bm_bx_dt <= (1080 + 30) &
      bm_bx_dt >= (1080 - 30)                         ~ 1080,
    bm_bx_dt <= (1260 + 30) &
      bm_bx_dt >= (1260 - 30)                         ~ 1260,
    bm_bx_dt <= (1440 + 30) &
      bm_bx_dt >= (1440 - 30)                         ~ 1440,
    bm_bx_dt <= (1800 + 30) &
      bm_bx_dt >= (1800 - 30)                         ~ 1800,
    bm_bx_dt <= (1980 + 30) &
      bm_bx_dt >= (1980 - 30)                         ~ 1980,
    bm_bx_dt <= (2160 + 30) &
      bm_bx_dt >= (2160 - 30)                         ~ 2160,
    bm_bx_dt <= (2340 + 30) &
      bm_bx_dt >= (2340 - 30)                         ~ 2340,
    bm_bx_dt <= (2520 + 30) &
      bm_bx_dt >= (2520 - 30)                         ~ 2520
  )) %>% 
  # focus on these
  filter(!is.na(blast_visit_time_range_days)) %>% 
  # There is multiple values for 1 range, take the closest
  group_by(mdsepid, blast_visit_time_range_days) %>% 
  mutate(int = abs(blast_visit_time_range_days - bm_bx_dt)) %>% 
  # keep the closest to enrolment and earliest if multiple
  arrange(mdsepid, blast_visit_time_range_days, int) %>% 
  distinct(mdsepid, .keep_all = TRUE) %>% 
  select(-int) %>% 
  arrange(blast_visit_time_range_days, mdsepid)

# marrow_1 <- marrow1 %>% 
#   # Make wide data
#   pivot_wider(id_cols = mdsepid, 
#               names_from = c(blast_visit_time_range_days), 
#               values_from = c(bm_bx_dt, bm_blast_pct), 
#               names_vary = "slowest") 

# Marrow the other way----
marrow1_ <-
  marrow1 %>% 
  select(mdsepid, bm_blast_pct_dt = bm_bx_dt, bm_blast_pct = blast_pct)


###################################################################### III ### Join----
# mds_data1 <- mds_data %>% 
#   left_join(., wide_blood_count, by = "mdsepid") %>% 
#   left_join(., marrow_1, by = "mdsepid") %>% 
#   left_join(., mutation_data2, by = "mdsepid")

dat <- wide_blood_data %>% 
  # We have the most patient data in blood so i am using this one as a 
  # "reference" to find other data with close dates
  full_join(., marrow1_, by = "mdsepid") %>% 
  mutate(time_lab_blast = abs(lab_dt - bm_blast_pct_dt)) %>% 
  # keep the closest to blood
  arrange(mdsepid, time_lab_blast) %>% 
  filter(time_lab_blast <= 30 | is.na(time_lab_blast)) %>% 
  # Add mutation data to find close date from blood
  full_join(., mutation_data1, by = "mdsepid") %>% 
  # mutate(has_mutation = case_when(
  #   has_mutation == "Yes"                ~ "Yes",
  #   is.na(has_mutation)                  ~ "No"
  # )) %>% 
  mutate(time_lab_mutation = abs(lab_dt - collection_dt_for_serie)) %>% 
  # keep the closest to blood
  filter(time_lab_mutation <= 30 | is.na(time_lab_mutation)) %>% 
  arrange(mdsepid, time_lab_mutation) %>% 
  left_join(mds_data, ., by = "mdsepid") %>% 
  
  # calculate score
  mutate(cytogenetics_group_final = case_when(
    from_cytogenetics_discordant_file_ipssr_total_new.x == 99            ~ NA_character_,
    !is.na(from_cytogenetics_discordant_file_cyto_group_reviewed.x)      ~ str_to_sentence(from_cytogenetics_discordant_file_cyto_group_reviewed.x),
    TRUE                                                                 ~ mds_cyto_r_ips_prog_sub_g_recoded
  )) %>% 
  mutate(cyto_group_ipssr_point = case_when(
    cytogenetics_group_final == "Very good"     ~ 0,
    cytogenetics_group_final == "Good"          ~ 1,
    cytogenetics_group_final == "Intermediate"  ~ 2,
    cytogenetics_group_final == "Poor"          ~ 3,
    cytogenetics_group_final == "Very poor"     ~ 4
  )) %>% 
  mutate(across(c("bm_blast_pct"), ~ case_when(
    . <= 2         ~ 0,
    . > 2 &
      . < 5        ~ 1,
    . >= 5 &
      . <= 10      ~ 2,
    . > 10         ~ 3
  ), .names = "{.col}_ipssr_point")) %>% 
  mutate(across(c("series_lab_result_hgb"), ~ case_when(
    . < 8         ~ 1.5,
    . >= 8 &
      . < 10        ~ 1,
    . >= 10         ~ 0
  ), .names = "{.col}_ipssr_point")) %>% 
  mutate(across(c("series_lab_result_plt"), ~ case_when(
    . < 50000         ~ 1,
    . >= 50000 &
      . < 100000        ~ 0.5,
    . >= 100000         ~ 0
  ), .names = "{.col}_ipssr_point")) %>% 
  mutate(across(c("series_lab_result_anc"), ~ case_when(
    . < 800         ~ 0.5,
    . >= 800         ~ 0
  ), .names = "{.col}_ipssr_point")) %>% 
  
  mutate(ipssr_score_for_series = rowSums(select(.,c(ends_with("_ipssr_point"))), 
                                              na.rm = FALSE)) %>% 
  # select(mdsepid, ends_with("_ipssr_point"), ipssr_score_for_series)
  # select(mdsepid, collection_dt_for_serie : time_lab_mutation, series_lab_result_plt) %>% 
  
  # Fix for no mut or no splicing
  mutate(across(c(starts_with("number_of_mutation_for_serie")), ~ case_when(
    is.na(.)         ~ 0,
    . < 2            ~ 0,
    . >= 2           ~ 3
  ), .names = "{.col}_ccrs_point")) %>%
  mutate(across(c(starts_with("slicing_mutation_for_serie")), ~ case_when(
    is.na(.)         ~ 0,
    . == "splicing"  ~ 2
  ), .names = "{.col}_ccrs_point")) %>%
  mutate(across(c(starts_with("series_lab_result_plt")), ~ case_when(
    . < 100000           ~ 2.5,
    . >= 100000          ~ 0
  ), .names = "{.col}_ccrs_point")) %>% 
  mutate(ccrs_score_for_series = rowSums(select(.,c(ends_with("_ccrs_point"))), 
                                          na.rm = FALSE))
  # select(mdsepid, ends_with("_ipssr_point"), ipssr_score_for_series,
  #        ends_with("_ccrs_point"), ccrs_score_for_series)
  



# # Add point---- dos not work
# mds_data2 <- mds_data1 %>% 
#   mutate(cyto_group_ipssr_point = case_when(
#     mds_cyto_r_ips_prog_sub_g_recoded == "Very good"     ~ 0,
#     mds_cyto_r_ips_prog_sub_g_recoded == "Good"          ~ 1,
#     mds_cyto_r_ips_prog_sub_g_recoded == "Intermediate"  ~ 2,
#     mds_cyto_r_ips_prog_sub_g_recoded == "Poor"          ~ 3,
#     mds_cyto_r_ips_prog_sub_g_recoded == "Very poor"     ~ 4
#   )) %>% 
#   mutate(across(c(starts_with("bm_blast_pct_")), ~ case_when(
#                   . <= 2         ~ 0,
#                   . > 2 &
#                     . < 5        ~ 1,
#                   . >= 5 &
#                     . <= 10      ~ 2,
#                   . > 10         ~ 3
#                 ), .names = "{.col}_ipssr_point")) %>% 
#   mutate(across(c(starts_with("series_lab_result_hgb_")), ~ case_when(
#     . < 8         ~ 1.5,
#     . >= 8 &
#       . < 10        ~ 1,
#     . >= 10         ~ 0
#   ), .names = "{.col}_ipssr_point")) %>% 
#   mutate(across(c(starts_with("series_lab_result_plt_")), ~ case_when(
#     . < 50         ~ 1,
#     . >= 50 &
#       . < 100        ~ 0.5,
#     . >= 100         ~ 0
#   ), .names = "{.col}_ipssr_point")) %>% 
#   mutate(across(c(starts_with("series_lab_result_anc_")), ~ case_when(
#     . < 0.8         ~ 0.5,
#     . >= 0.8         ~ 0
#   ), .names = "{.col}_ipssr_point")) %>% 
#   
#   mutate(ipssr_score_for_series_180 = rowSums(select(.,c(cyto_group_ipssr_point, 
#                                                          ends_with("180_ipssr_point"))), 
#                                               na.rm = FALSE)) %>% 
#   
#   # mutate(ipssr_score_for_series_180 = rowSums(select(.,c(cyto_group_ipssr_point, ends_with("180_ipssr_point"))), na.rm = FALSE)) %>% 
#   select(mdsepid, cyto_group_ipssr_point, ends_with("180_ipssr_point"), ipssr_score_for_series_180)





new_score_data <- dat %>% 
  mutate(ipssr_2cat_for_series = case_when(
    final_class == "MDS" &
      ipssr_score_for_series <= 3.5                   ~ "Lower risk disease",
    final_class == "MDS" &
      ipssr_score_for_series > 3.5                    ~ "Higher risk disease"
  ), ipssr_2cat_for_series = factor(ipssr_2cat_for_series, levels = c("Lower risk disease", 
                                                                      "Higher risk disease"))
  ) %>% 
  mutate(ccrs_3cat_for_series = case_when(
    final_class %in% c("CCUS") &
      ccrs_score_for_series < 2.5                   ~ "Lower risk disease",
    final_class %in% c("CCUS") &
      ccrs_score_for_series >= 2.5 &
      ccrs_score_for_series < 5                     ~ "Intermediate risk disease",
    final_class %in% c("CCUS") &
      ccrs_score_for_series >= 5                    ~ "Higher risk disease"
  ), ccrs_3cat_for_series = factor(ccrs_3cat_for_series, levels = c("Lower risk disease", 
                                                                    "Intermediate risk disease",
                                                                    "Higher risk disease"))
  )
  
  

# Save
write_rds(new_score_data, 
          paste0(here::here(), 
                 "/data/processed_data",
                 "/MDSSerialAnalysis_NewScore_",
                 str_remove_all(today(), "-"), ".rds"))
write_csv(new_score_data,
          paste0(here::here(), 
                 "/data/processed_data",
                 "/MDSSerialAnalysis_NewScore_",
                 str_remove_all(today(), "-"), ".csv"))

write_csv(new_score_data, 
          paste0(path_project, 
                 "/ProcessedData",
                 "/MDSSerialAnalysis_NewScore_",
                 str_remove_all(today(), "-"), ".csv"))


# End data cleaning

