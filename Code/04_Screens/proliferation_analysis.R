library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

screens_data_dir <- here(external_data_dir, "Screens")
output_data_dir <- output_data_base_dir
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet"))
cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_procan.parquet"))

df_proliferation <- cellline_buf_filtered_procan %>%
  inner_join(y = df_growth, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never")

dc_growth_cor <- cor.test(x = df_proliferation$Buffering.CellLine.Ratio.ZScore,
                          y = df_proliferation$CellLine.GrowthRatio,
                          method = "spearman")

cat(capture.output(dc_growth_cor), file = here(reports_dir, "dosage_compensation_proliferation_correlation.txt"),
    append = FALSE, sep = "\n")
