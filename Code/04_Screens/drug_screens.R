library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "DrugScreen")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
drug_screens <- read_parquet(here(output_data_dir, "drug_screens.parquet"))
cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
cellline_buf_filtered_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))

sensitivity_procan <- cellline_buf_filtered_procan %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never")

drug_dc_corr_procan <- sensitivity_procan %>%
  group_by(Drug.ID, Drug.Name) %>%
  summarize(Correlation.Sensitivity_Buffering = cor.test(Buffering.CellLine.Ratio, Drug.MFI.Log2FC,
                                                         method = "pearson")$estimate[["cor"]]) %>%
  arrange(Drug.Name)

sensitivity_depmap <- cellline_buf_filtered_depmap %>%
  inner_join(y = drug_screens, by = "CellLine.Name",
             relationship = "one-to-many", na_matches = "never")

drug_dc_corr_depmap <- sensitivity_depmap %>%
  group_by(Drug.ID, Drug.Name) %>%
  summarize(Correlation.Sensitivity_Buffering = cor.test(Buffering.CellLine.Ratio.ZScore, Drug.MFI.Log2FC,
                                                         method = "pearson")$estimate[["cor"]]) %>%
  arrange(Drug.Name)

cor.test(drug_dc_corr_depmap$Correlation.Sensitivity_Buffering,
         drug_dc_corr_procan$Correlation.Sensitivity_Buffering,
         method = "pearson")

# === Write results ===
write.xlsx(drug_dc_corr_procan, here(tables_base_dir, "sensitivity_correlation_procan_gene_filtered.xlsx"),
           colNames = TRUE)
write.xlsx(drug_dc_corr_depmap, here(tables_base_dir, "sensitivity_correlation_depmap_gene_filtered.xlsx"),
           colNames = TRUE)
