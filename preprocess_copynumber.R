library(here)
library(tidyr)
library(dplyr)
library(arrow)

here::i_am("DosageCompensationFactors.Rproj")

source(here("parameters.R"))

copynumber_data_dir <- here(external_data_dir, "CopyNumber")
output_data_dir <- output_data_base_dir
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

df_model <- read_csv_arrow(here(copynumber_data_dir, "Model.csv")) %>%
  select(ModelID, SangerModelID) %>%
  rename(CellLine.DepMapModelId = "ModelID",
         CellLine.SangerModelId = "SangerModelID")

# NOTE! SIDM00018 has a different cell line name in DepMap and expression datasets (KO52 vs K052)!!!
df_cn <-
  read.csv(
    here(copynumber_data_dir, "Copy_Number_Public_23Q2.csv"),
    row.names = 1,
    check.names = FALSE
  ) %>%
  as.matrix() %>%
  as.table() %>%
  as.data.frame() %>%
  setNames(c("CellLine.DepMapModelId", "Gene.Symbol", "Gene.CopyNumber"))  %>%
  inner_join(y = df_model, by = "CellLine.DepMapModelId",
             na_matches = "never", relationship = "many-to-one")

write_parquet(df_cn, here(output_data_dir, 'copy_number.parquet'),
              version = "2.6")