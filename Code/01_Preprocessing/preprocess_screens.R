library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)
here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

copynumber_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

df_celllines <- read_csv_arrow(here(copynumber_data_dir, "Model.csv")) %>%
  select(ModelID, CellLineName, SangerModelID) %>%
  rename(CellLine.DepMapModelId = "ModelID",
         CellLine.SangerModelId = "SangerModelID",
         CellLine.Name = "CellLineName")

df_prism <- read_csv_arrow(here(screens_data_dir, "PRISM_Repurposing_Public_23Q2.csv")) %>%
  rename(CellLine.DepMapModelId = 1) %>%
  pivot_longer(everything() & !CellLine.DepMapModelId,
               names_to = "Drug", values_to = "Drug.MFI.Log2FC") %>%
  mutate(Drug.Name = str_extract(Drug, "(.*)\\s+\\((.*)\\)", group = 1),
         Drug.ID = str_extract(Drug, "(.*)\\s+\\((.*)\\)", group = 2)) %>%
  left_join(y = df_celllines, by = "CellLine.DepMapModelId",
            relationship = "many-to-one", na_matches = "never")

write_parquet(df_prism, here(output_data_dir, 'drug_screens.parquet'),
              version = "2.6")