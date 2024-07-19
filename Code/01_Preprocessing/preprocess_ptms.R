library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))

ubi_dir <- here(external_data_dir, "PTMs", "CPTAC", "Ubiquitylome_CDAP_v1")
output_data_dir <- output_data_base_dir

df_ubi <- read.table(here(ubi_dir, "CPTAC3_Lung_Squamous_Cell_Carcinoma_Ubiquitylome.ubiquitylsite.tmt11.tsv"),
                     sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
  janitor::row_to_names(1)

df_ubi_meta <- read.table(here(external_data_dir, "PDC_biospecimen_manifest.tsv"),
                          sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
    janitor::row_to_names(1)

aliquot_id_mapping <- df_ubi_meta %>%
  rename(Model.ID = "Case Submitter ID",
         Model.AliquotSubmitterID = "Aliquot Submitter ID") %>%
  select(Model.ID, Model.AliquotSubmitterID)

df_ubi_tidy <- df_ubi %>%
  mutate(across(matches("Log Ratio"), as.numeric)) %>%
  pivot_longer(matches("Log Ratio"),
               names_to = "Model.AliquotSubmitterID",
               values_to = "Protein.Ubiquitination.Log2") %>%
  mutate(Model.AliquotSubmitterID = str_split_i(Model.AliquotSubmitterID, " ", 1)) %>%
  left_join(aliquot_id_mapping, by = "Model.AliquotSubmitterID") %>%
  rename(Gene.Symbol = Gene,
         Protein.Ubiquitination.Site = Ubiquitylsite,
         Protein.Ubiquitination.Peptide = Peptide) %>%
  select(Model.ID, Model.AliquotSubmitterID, Gene.Symbol,
         Protein.Ubiquitination.Peptide, Protein.Ubiquitination.Site, Protein.Ubiquitination.Log2)

df_ubi_tidy %>%
  write_parquet(here(output_data_dir, 'ubiquitination_cptac.parquet'),
              version = "2.6")
