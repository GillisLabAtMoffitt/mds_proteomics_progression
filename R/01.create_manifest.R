## Import library
library(tidyverse)
library(readxl)
library(janitor)


###################################################################### I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "MDSProteomics")
all_sheets <- excel_sheets(paste0(path, 
                                  "/RawData/MDSProteomics_SampleManifests",
                                  "/MDSProteomics_SampleManifestSubmitted_20250910_ReplicatesIndicated.xlsx"))
all_sheets <- all_sheets[-2]
data_list <- list()
n <- 1
for (sheet in all_sheets) {
  data_list[[sheet]] <- read_excel(paste0(path, 
                                          "/RawData/MDSProteomics_SampleManifests",
                                          "/MDSProteomics_SampleManifestSubmitted_20250910_ReplicatesIndicated.xlsx"),
                                   sheet = sheet, skip = 2) %>% 
    mutate(plate = n)
  n <- n + 1
}
manifest <- bind_rows(data_list) %>% 
  clean_names() %>% 
  select(plate_well_id, plate, sample_name, concentration_ng_ul)

cmml_data <-
  readxl::read_xlsx(paste0(path, 
                           "/RawData/MDSProteomics_SampleManifests",
                           "/MDSProteomics_CMML_SampleLinking_20250909.xlsx")) %>% 
  clean_names()
mds_data <-
  readxl::read_xlsx(paste0(path, 
                           "/RawData/MDSProteomics_NHSClinicalData_20251002",
                           "/Revised_picklist_Padron_Walter_Gillis_pilot_PB_serum_withIDsAndDays.xlsx"), na = "NA") %>% 
  clean_names()

rm(all_sheets, data_list, n)


###################################################################### II ### Merge data
mds_data1 <- mds_data %>% 
  group_by(original_s_sampleid) %>% 
  fill(mdsepid, anatomicsite, sampletypedesc, sample_id_used_for_pilot, collectiondt_days, .direction = "down") %>% 
  mutate(sample_id = paste0("MDS_", sample_id_used_for_pilot)) %>% 
  mutate(disease_type = "MDS") %>% 
  select(sample_id, disease_type, study_id = mdsepid, 
         anatomic_site = anatomicsite, 
         sample_type = sampletypedesc, 
         sample_num, original_s_sampleid,
         sample_id_used_for_pilot, collectiondt_days)

cmml_data1 <- cmml_data %>% 
  filter(!is.na(tcc_id)) %>% 
  mutate(sample_id = paste0("CMML_", sequence_number)) %>% 
  mutate(disease_type = "CMML") %>% 
  select(sample_id, disease_type, study_id = tcc_id, 
         sample_type, controls, 
         collection_date) %>% 
  mutate(study_id = as.character(study_id))

disease_data <- bind_rows(mds_data1, 
                          cmml_data1)


manifest1 <- manifest %>% 
  mutate(sample_id = str_remove(sample_name, " (.*)")) %>% 
  full_join(., 
            disease_data, 
            by = "sample_id") %>% 
  mutate(plate_well_id = case_when(
    !is.na(plate_well_id)         ~ plate_well_id,
    is.na(plate_well_id)          ~ "not analyzed or analyzed using the original id"
  ))
  

write_csv(manifest1, 
          paste0("data/processed_data",
                 "/MDSProteomics_final_manifest_", 
                 str_remove_all(today(), "-"), ".csv"))

write_csv(manifest1, 
          paste0(path, "/ProcessedData",
                 "/MDSProteomics_final_manifest_", 
                 str_remove_all(today(), "-"), ".csv"))


# End Create manifest----

