## Import library
library(tidyverse)
library(readxl)
library(janitor)


###################################################################### I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "MDSProteomics")

manifest <- 
  readxl::read_xlsx(paste0(path, "/ProcessedData",
                           "/MDSProteomics_SampleManifest_v2_20251013.xlsx"))
clinical_data <-
  read_csv(paste0(path, 
                  "/RawData/MDSProteomics_CMMLClinicalData",
                  "/MDSProteomics_CMML_ClinicalData_20251021.csv"))

plate1_read <-
  read.csv(paste0(path, 
                  "/RawData/MDSProteomics_OlinkData_20251017",
                  "/Plate1_data_intensity_nmlzd_Extended_NPX_20251017.csv"),
           sep=";")
plate2_read <-
  read.csv(paste0(path, 
                  "/RawData/MDSProteomics_OlinkData_20251017",
                  "/Plate2_data_intensity_nmlzd_Extended_NPX_20251017.csv"),
           sep=";")

###################################################################### II ### Data wrangling
manifest <- manifest %>% 
  filter(disease_type == "CMML") %>% 
  select(SampleID = sample_name, WellID = plate_well_id, 
         plate, sample_id, study_id, collection_date) %>% 
  mutate(collection_date = as.Date(collection_date)) %>% 
  distinct(study_id, collection_date, .keep_all = TRUE)

clinical_data <- clinical_data %>% 
  select(SampleID, mrn, collection_date,
         overallpd : death_status) %>% 
  distinct() %>% 
  mutate(collection_date = case_when(
    !is.na(collection_date)             ~ collection_date,
    is.na(collection_date)              ~ as.Date("2024-07-10")
  )) %>% 
  mutate(mrn = as.character(mrn))

plate_read <- 
  bind_rows(plate1_read, plate2_read, .id = "plate") %>%
  filter(str_detect(SampleID, "CMML_")) %>% 
  filter(AssayQC == "PASS") %>% 
  select(SampleID, plate, WellID, PCNormalizedNPX, Assay)
# a <- plate_read %>% distinct(SampleID, Assay)
rm(plate1_read, plate2_read)

wide_reads <- manifest %>% 
  # Use inner join to filter out the duplicate samples
  inner_join(., clinical_data, 
             by = c("SampleID", "collection_date")) %>% 
  # distinct(SampleID, WellID, plate)
  arrange(mrn, collection_date) %>% 
  select(mrn, collection_date, everything()) %>% 
  # distinct(mrn, SampleID, collection_date) %>% 
  group_by(mrn) %>% 
  mutate(sample_seq_num = row_number(), .after = mrn) %>% 
  mutate(sample_seq_num = paste0("date", sample_seq_num)) %>% 
  ungroup() %>% 
  inner_join(., plate_read,
             by = c("SampleID", "WellID", "plate")) %>% 
  select(-c(SampleID, WellID, plate, sample_id))



wide_reads |>
  dplyr::summarise(n = dplyr::n(), .by = c(mrn, collection_date, Assay)) |>
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
              values_from = -c(mrn, sample_seq_num, study_id, overallpd, overallpddate,
                               enrolldate, os, acuteevent, death, death_status),
              names_vary = "slowest")


write_csv(wide_cmml_data, 
          paste0("data/processed_data",
                 "/MDSProteomics_wide_cmml_data_", 
                 str_remove_all(today(), "-"), ".csv"))

write_csv(wide_cmml_data, 
          paste0(path, "/ProcessedData",
                 "/MDSProteomics_wide_cmml_data_", 
                 str_remove_all(today(), "-"), ".csv"))

















