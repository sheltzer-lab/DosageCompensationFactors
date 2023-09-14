library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
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


data_density <- df_celllines_filtered %>%
  summarize_all(~ sum(!is.na(.x)) / nrow(df_celllines_filtered)) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "Density") %>%
  arrange(desc(Density))

sorted_violin_plot <- function(df, x, y) {
  df %>%
  add_count({ { x } }) %>%
  filter(n > 2) %>%
  group_by({ { x } }) %>%
  mutate(Median = median({ { y } })) %>%
  ungroup() %>%
  arrange(Median) %>%
  mutate(cancer_type = factor({ { x } }, levels = unique({ { x } }))) %>%
  violin_plot({ { x } }, { { y } })
}

df_model <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20230801.csv")) %>%
  rename(CellLine.Name = "model_name") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

df_celllines_filtered <- df_model %>%
  inner_join(y = cellline_buf_filtered_procan, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never")

df_metastasis <- df_celllines_filtered %>%
  filter(tissue_status == "Metastasis")

df_tumour <- df_celllines_filtered %>%
  filter(tissue_status == "Tumour")

t.test(df_tumour$Buffering.CellLine.Ratio,
       df_metastasis$Buffering.CellLine.Ratio,
       paired = FALSE, var.equal = TRUE)

df_cancer_type <- df_celllines_filtered %>%
  select(Buffering.CellLine.Ratio, cancer_type) %>%
  group_by(cancer_type) %>%
  skimr::skim()

plot_cancer_type <- df_celllines_filtered %>%
  sorted_violin_plot(cancer_type, Buffering.CellLine.Ratio)

plot_tissue_status <- df_celllines_filtered %>%
  sorted_violin_plot(tissue_status, Buffering.CellLine.Ratio)

plot_msi <- df_celllines_filtered %>%
  sorted_violin_plot(msi_status, Buffering.CellLine.Ratio)

df_msi <- df_celllines_filtered %>%
  filter(msi_status == "MSI")
df_mss <- df_celllines_filtered %>%
  filter(msi_status == "MSS")

plot_smoking <- df_celllines_filtered %>%
  sorted_violin_plot(smoking_status, Buffering.CellLine.Ratio)

plot_gender <- df_celllines_filtered %>%
  sorted_violin_plot(gender, Buffering.CellLine.Ratio)

plot_ethnicity <- df_celllines_filtered %>%
  sorted_violin_plot(ethnicity, Buffering.CellLine.Ratio)

t.test(df_mss$Buffering.CellLine.Ratio,
       df_msi$Buffering.CellLine.Ratio,
       paired = FALSE, var.equal = FALSE)


df_celllines_filtered %>%
  scatter_plot_regression(ploidy_wes, Buffering.CellLine.Ratio,  Buffering.CellLine.Ratio ~ ploidy_wes)

df_celllines_filtered %>%
  scatter_plot_regression(mutational_burden, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ mutational_burden,
                          label_coords = c(500, -0.4))

df_celllines_filtered %>%
  scatter_plot_regression(age_at_sampling, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ age_at_sampling,
                          label_coords = c(50, -0.4))

df_celllines_filtered %>%
  scatter_plot_regression(age_at_sampling, ploidy_wes, ploidy_wes ~ age_at_sampling)