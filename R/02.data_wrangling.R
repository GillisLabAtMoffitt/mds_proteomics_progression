## Import library
library(tidyverse)
library(readxl)
library(janitor)


###################################################################### I ### CMML Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "MDSProteomics")

manifest <- 
  readxl::read_xlsx(paste0(path, "/ProcessedData",
                           "/MDSProteomics_SampleManifest_v2_20251013.xlsx"))
clinical_data <-
  read_csv(paste0(path, 
                  "/RawData/MDSProteomics_CMMLClinicalData",
                  "/MDSProteomics_CMML_Merged_Long_20251027.csv"))
enrollment_data <-
  read_csv(paste0(path, 
                  "/RawData/MDSProteomics_CMMLClinicalData",
                  "/MDSProteomics_CMML_ClinicalData_20251021.csv"))

# plate1_read <-
#   read.csv(paste0(path, 
#                   "/RawData/MDSProteomics_OlinkData_20251017",
#                   "/Plate1_data_intensity_nmlzd_Extended_NPX_20251017.csv"),
#            sep=";")
# plate2_read <-
#   read.csv(paste0(path, 
#                   "/RawData/MDSProteomics_OlinkData_20251017",
#                   "/Plate2_data_intensity_nmlzd_Extended_NPX_20251017.csv"),
#            sep=";")
plate_read <-
  read.csv(paste0(path, 
                  "/RawData/MDSProteomics_OlinkData_20251021",
                  "/Gillis_3plates_Extended_NPX_2025-10-21.csv"))

###################################################################### II ### CMML data wrangling
manifest <- manifest %>% 
  filter(disease_type == "CMML") %>% 
  select(SampleID = sample_name, WellID = plate_well_id, 
         plate, sample_id, study_id, collection_date) %>% 
  mutate(collection_date = as.Date(collection_date)) %>% 
  distinct(study_id, collection_date, .keep_all = TRUE)

enrollment_data <- enrollment_data %>% 
  select(SampleID, mrn, enrolldate) %>% 
  mutate(mrn = as.character(mrn)) %>% 
  distinct()

clinical_data <- clinical_data %>% 
  filter(str_detect(SampleID, "CMML_")) %>% 
  select(SampleID, MRN, collection_date,
         PD_overall : AIE) %>% 
  distinct() %>% 
  mutate(collection_date = case_when(
    !is.na(collection_date)             ~ collection_date,
    is.na(collection_date)              ~ as.Date("2024-07-10")
  )) %>% 
  mutate(MRN = as.character(MRN))

plate_read <- plate_read %>% 
  # bind_rows(plate1_read, plate2_read, .id = "PlateID") %>%
  filter(str_detect(SampleID, "CMML_")) %>% 
  filter(AssayQC == "PASS") %>% 
  select(SampleID, PlateID, WellID, PCNormalizedNPX, Assay)
# a <- plate_read %>% distinct(SampleID, Assay)
# rm(plate1_read, plate2_read)

wide_reads <- manifest %>% 
  # Use inner join to filter out the duplicate samples
  inner_join(., clinical_data, 
             by = c("SampleID", "collection_date")) %>% 
  inner_join(., enrollment_data, 
             by = c("SampleID", "MRN" = "mrn")) %>% 
  # distinct(SampleID, WellID, PlateID)
  arrange(MRN, collection_date) %>% 
  select(MRN, collection_date, everything()) %>% 
  # distinct(MRN, SampleID, collection_date) %>% 
  group_by(MRN) %>% 
  mutate(sample_seq_num = row_number(), .after = MRN) %>% 
  mutate(sample_seq_num = paste0("date", sample_seq_num)) %>% 
  ungroup() %>% 
  mutate(PlateID = paste0("GillisPadron_plate", plate)) %>% 
  select(-plate) %>% 
  inner_join(., plate_read,
             by = c("SampleID", "WellID", "PlateID")) %>% 
  select(-c(SampleID, WellID, PlateID, sample_id))



wide_reads |>
  dplyr::summarise(n = dplyr::n(), .by = c(MRN, collection_date, Assay)) |>
  dplyr::filter(n > 1L) 

colnames(wide_reads)

wide_cmml_data <- wide_reads %>% 
  pivot_wider(#id_cols = c(mrn, collection_date), 
              names_from = c(Assay), 
              values_from = PCNormalizedNPX
  ) %>% 
  rename(collection = collection_date) %>% 
  pivot_wider(#id_cols = mrn, 
              names_from = sample_seq_num, 
              values_from = -c(MRN, sample_seq_num, study_id, PD_overall, PD_date,
                               enrolldate,
                               os, AIE, # AIE = acuteevent, 
                               # death, 
                               death_status),
              names_vary = "slowest")

write_csv(wide_cmml_data,
          paste0("data/processed_data",
                 "/MDSProteomics_wide_cmml_rawdata_",
                 str_remove_all(today(), "-"), ".csv"))

write_csv(wide_cmml_data,
          paste0(path, "/ProcessedData",
                 "/MDSProteomics_wide_cmml_rawdata_",
                 str_remove_all(today(), "-"), ".csv"))


###################################################################### III ### CMML Add and rename variables
wide_cmml_data <- wide_cmml_data %>% 
  rename(vital_status = death_status,
         os_date = os) %>% 
  mutate(os_event = case_when(
    vital_status == "DEAD"                     ~ 1,
    vital_status == "ALIVE"                    ~ 0
  ), .after = vital_status) %>% 
  mutate(os_time_from_enroll_months = interval(start = enrolldate,
                                           end = os_date)/
           duration(n = 1, unit = "months"), .after = os_date) %>%
  mutate(pfs_event = case_when(
    PD_overall == 1                            ~ 1,
    vital_status == "DEAD"                     ~ 1,
    vital_status == "ALIVE"                    ~ 0
  ), .after = PD_date) %>% 
  mutate(pfs_date = coalesce(PD_date, os_date), .after = pfs_event) %>% 
  mutate(pfs_time_from_enroll_months = interval(start = enrolldate,
                                                  end = pfs_date)/
           duration(n = 1, unit = "months"), .after = pfs_date)

write_csv(wide_cmml_data, 
          paste0("data/processed_data",
                 "/MDSProteomics_CMML_dataWide_", 
                 str_remove_all(today(), "-"), ".csv"))

write_csv(wide_cmml_data, 
          paste0(path, "/ProcessedData",
                 "/MDSProteomics_CMML_dataWide_", 
                 str_remove_all(today(), "-"), ".csv"))


# END CMML----

# Start MDS----

## Import library
library(tidyverse)
library(readxl)
library(janitor)


###################################################################### I ### MDS Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "MDSProteomics")

manifest <- 
  readxl::read_xlsx(paste0(path, "/ProcessedData",
                           "/MDSProteomics_SampleManifest_v2_20251013.xlsx"))
clinical_data <-
  read_csv(paste0(path, 
                  "/RawData/MDSProteomics_NHSClinicalData_20251002",
                  "/demo_disp_data_20250919.csv")) %>% 
  clean_names()
progression_data <-
  read_csv(paste0(path, 
                  "/RawData/MDSProteomics_NHSClinicalData_20251002",
                  "/progression_data_20251017.csv")) %>% 
  clean_names()
plate_read <-
  read.csv(paste0(path, 
                  "/RawData/MDSProteomics_OlinkData_20251021",
                  "/Gillis_3plates_Extended_NPX_2025-10-21.csv"))

###################################################################### II ### MDS data wrangling
manifest <- manifest %>% 
  filter(disease_type == "MDS") %>% 
  select(sample_name, WellID = plate_well_id, 
         plate, SampleID = sample_id, mdsepid = study_id, collectiondt_days,
         original_sampleid, sample_id_used_for_pilot) %>% 
  distinct(mdsepid, collectiondt_days, .keep_all = TRUE)

progression_data <- progression_data %>% 
  group_by(mdsepid) %>% 
  # mutate(n = row_number()) %>% # same as progression_num
  mutate(first_progression_time = case_when(
    progression_num == 1             ~ progression_date
  )) %>% 
  mutate(first_progression_type = case_when(
    progression_num == 1             ~ progression_event
  )) %>% 
  fill(first_progression_time, first_progression_type, .direction = "updown") %>% 
  ungroup() %>% 
  select(-c(progression_category, progression_num))








