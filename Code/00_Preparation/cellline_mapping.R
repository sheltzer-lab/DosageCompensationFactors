library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

copynumber_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
output_data_dir <- output_data_base_dir

dir.create(output_data_dir, recursive = TRUE)

# Load cell line model file from DepMap and define cell line mapping with custom ID
df_celllines <- read_csv_arrow(here(copynumber_data_dir, "Model.csv")) %>%
  select(ModelID, SangerModelID, CellLineName) %>%
  rename(CellLine.DepMapModelId = "ModelID",
         CellLine.SangerModelId = "SangerModelID",
         CellLine.Name = "CellLineName") %>%
  add_count(CellLine.SangerModelId, name = "SangerReplicates") %>%
  add_count(CellLine.DepMapModelId, name = "DepMapReplicates") %>%
  mutate(SangerReplicates = if_else(is.na(CellLine.SangerModelId), 0, SangerReplicates),
         DepMapReplicates = if_else(is.na(CellLine.DepMapModelId), 0, DepMapReplicates)) %>%
  filter(SangerReplicates < 2 & DepMapReplicates < 2) %>%
  mutate(Model.ID = paste0("CL_", consecutive_id(CellLine.DepMapModelId,
                                                          CellLine.SangerModelId,
                                                          CellLine.Name))) %>%
  select(Model.ID, CellLine.DepMapModelId, CellLine.SangerModelId, CellLine.Name) %>%
  write_parquet(here(output_data_dir, 'celllines.parquet'), version = "2.6")
