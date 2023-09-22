library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Screens")
output_data_dir <- output_data_base_dir
reports_dir <- reports_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet"))
cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_procan.parquet"))
cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_procan.parquet"))


df_proliferation <- cellline_buf_filtered_procan %>%
  inner_join(y = df_growth, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never")

dc_growth_cor <- cor.test(x = df_proliferation$Buffering.CellLine.Ratio.ZScore,
                          y = df_proliferation$CellLine.GrowthRatio,
                          method = "spearman")

cat(capture.output(dc_growth_cor), file = here(reports_dir, "dosage_compensation_proliferation_correlation.txt"),
    append = FALSE, sep = "\n")


dc_growth_reg <- df_proliferation %>%
  scatter_plot_regression(x = Buffering.CellLine.Ratio.ZScore,
                          y = CellLine.GrowthRatio,
                          formula = CellLine.GrowthRatio ~ Buffering.CellLine.Ratio.ZScore,
                          label_coords = c(0.25, 7.5),
                          title = paste0("ProCan filtered, Correlation: ",
                                         "ρ = ", format(round(dc_growth_cor$estimate[["rho"]], 3), nsmall = 3),
                                         ", p = ", format(round(dc_growth_cor$p.value, 3), nsmall = 3))) %>%
  save_plot("dosage_compensation_proliferation_procan_filtered.png")

