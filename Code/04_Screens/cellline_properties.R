library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Screens", "CellLineProperties")
output_data_dir <- output_data_base_dir
reports_dir <- reports_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_procan.parquet"))
cellline_buf_filtered_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_depmap.parquet"))

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(CellLine.Name, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20230801.csv")) %>%
  rename(CellLine.Name = "model_name") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.Name = "CellLineName") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

aneuploidy_quant <- quantile(copy_number$CellLine.AneuploidyScore, probs = c(0.25, 0.5, 0.75))

# ToDo: Use Cell Line ID instead of cell line name
df_procan <- df_model_procan %>%
  inner_join(y = cellline_buf_filtered_procan, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, TRUE, FALSE),
         `Near-Tetraploid` = if_else(ploidy_wes > 3.5, TRUE, FALSE),
         Aneuploidy = case_when(
           CellLine.AneuploidyScore > aneuploidy_quant["75%"] ~ "Very High",
           CellLine.AneuploidyScore > aneuploidy_quant["50%"] ~ "High",
           CellLine.AneuploidyScore > aneuploidy_quant["25%"] ~ "Low",
           TRUE ~ "Very Low",
         ))


df_depmap <- df_model_depmap %>%
  inner_join(y = cellline_buf_filtered_depmap, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, TRUE, FALSE),
         `Near-Tetraploid` = if_else(CellLine.Ploidy > 3.5, TRUE, FALSE),
         Aneuploidy = case_when(
           CellLine.AneuploidyScore > aneuploidy_quant["75%"] ~ "Very High",
           CellLine.AneuploidyScore > aneuploidy_quant["50%"] ~ "High",
           CellLine.AneuploidyScore > aneuploidy_quant["25%"] ~ "Low",
           TRUE ~ "Very Low",
         ))

plot_categorical_properties <- function (df, cols_categorical) {
  plots <- list()
  for (col in cols_categorical) {
    plots[[col]] <- df %>%
      sorted_violin_plot(col, Buffering.CellLine.Ratio)
  }
  return(plots)
}

data_density_procan <- df_procan %>%
  data_density()

data_density_depmap <- df_procan %>%
  data_density()

cancers_depmap <- df_depmap %>%
  distinct(OncotreeLineage, OncotreePrimaryDisease, OncotreeSubtype, OncotreeCode)

# Plot distributions for categorical variables
cols_procan <- c("tissue_status", "cancer_type", "msi_status",
                 "smoking_status", "gender", "ethnicity",
                 "WGD", "Near-Tetraploid", "Aneuploidy")

cols_depmap <- c("PrimaryOrMetastasis", "OncotreeSubtype", "Sex",
                 "WGD", "Near-Tetraploid", "Aneuploidy")

violoin_plots_procan <- df_procan %>%
  plot_categorical_properties(cols_procan)

violoin_plots_depmap <- df_depmap %>%
  plot_categorical_properties(cols_depmap)

## Plot cell line buffering per cancer type
df_procan %>%
  sorted_beeswarm_plot("cancer_type", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_procan.png", width = 300)
df_depmap %>%
  sorted_beeswarm_plot("OncotreeSubtype", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_depmap.png", width = 300)

### Plot aneuploidy score per cancer type
df_procan %>%
  sorted_beeswarm_plot("cancer_type", CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_aneuploidy_procan.png", width = 300)
df_depmap %>%
  sorted_beeswarm_plot("OncotreeSubtype", CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_aneuploidy_depmap.png", width = 300)

### Control for low aneuploidy score when plotting buffering per cancer type
df_procan %>%
  filter(Aneuploidy %in% c("Low", "Very Low")) %>%
  sorted_beeswarm_plot("cancer_type", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_low-aneuploidy_procan.png", width = 300)
df_depmap %>%
  filter(Aneuploidy %in% c("Low", "Very Low")) %>%
  sorted_beeswarm_plot("OncotreeSubtype", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_low-aneuploidy_depmap.png", width = 300)

# Statistical comparisons
## Difference between Leukemia and Non-Leukemia cancer cells
leuk_procan <- c("B-Cell Non-Hodgkin's Lymphoma", "B-Lymphoblastic Leukemia", "Acute Myeloid Leukemia")
leuk_depmap <- c("Diffuse Large B-Cell Lymphoma, NOS", "B-Lymphoblastic Leukemia/Lymphoma", "Acute Myeloid Leukemia")

df_procan %>%
  mutate(Leukemia_Lymphoma = cancer_type %in% leuk_procan) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_procan.png")
df_depmap %>%
  mutate(Lymphoid_Myeloid = OncotreeLineage %in% leuk_depmap) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_depmap.png")

### Control for low aneuploidy score
df_leuk_procan <- df_procan %>%
  filter(cancer_type %in% leuk_procan)
aneuploidy_leuk_max <- round(quantile(df_leuk_procan$CellLine.AneuploidyScore, probs = 0.9))

df_procan %>%
  filter(CellLine.AneuploidyScore <= aneuploidy_leuk_max) %>%
  mutate(Leukemia_Lymphoma = cancer_type %in% leuk_procan) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_low-aneuploidy_procan.png")
df_depmap %>%
  filter(CellLine.AneuploidyScore <= aneuploidy_leuk_max) %>%
  mutate(Leukemia_Lymphoma = OncotreeSubtype %in% leuk_depmap) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_low-aneuploidy_depmap.png")


## Tumor status
df_procan %>%
  signif_violin_plot(tissue_status, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_tumor-status_procan.png")
df_depmap %>%
  filter(PrimaryOrMetastasis != "NA") %>%
  signif_violin_plot(PrimaryOrMetastasis, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_tumor-status_depmap.png")

## Micro-satellite instability
df_procan %>%
  signif_beeswarm_plot(msi_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi_procan.png")
df_procan %>%
  signif_beeswarm_plot(msi_status, CellLine.AneuploidyScore,
                       test = wilcox.test, cex = 1) %>%
  save_plot("cellline_msi-aneuploidy_procan.png")

### Control MSI/MSS distribution for low low aneuploidy score
msi_procan <- df_procan %>% filter(msi_status == "MSI")
fivenum_msi_aneuploidy <- fivenum(msi_procan$CellLine.AneuploidyScore)
max_msi_aneuploidy <- max(msi_procan$CellLine.AneuploidyScore)

df_procan %>%
  filter(CellLine.AneuploidyScore <= max_msi_aneuploidy) %>%
  signif_beeswarm_plot(msi_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi_low-aneuploidy_procan.png")
df_procan %>%
  filter(CellLine.AneuploidyScore <= max_msi_aneuploidy) %>%
  signif_beeswarm_plot(msi_status, CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi-aneuploidy_low-aneuploidy_procan.png")


### Resample, so that Aneuploidy Score distributions between MSI and MSS are equal
df_split <- split(df_procan, df_procan$msi_status)
df_msi_equal <- df_split$MSI %>%
  equalize_distributions(df_split$MSS, CellLine.AneuploidyScore,
                         with_replacement = FALSE, num_buckets = 5)

df_msi_equal %>%
  signif_beeswarm_plot(msi_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi_equal-aneuploidy_procan.png")
df_msi_equal %>%
  signif_beeswarm_plot(msi_status, CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi-aneuploidy_equal-aneuploidy_procan.png")

## Sex/Gender
df_procan %>%
  filter(gender != "Unknown") %>%
  signif_violin_plot(gender, Buffering.CellLine.Ratio,
                     test = t.test) %>%
  save_plot("cellline_gender_procan.png")
df_depmap %>%
  filter(Sex != "Unknown") %>%
  signif_violin_plot(Sex, Buffering.CellLine.Ratio,
                     test = t.test) %>%
  save_plot("cellline_gender_depmap.png")

## Whole-genome doubling
df_procan %>%
  signif_violin_plot(WGD, Buffering.CellLine.Ratio,
                     test = t.test) %>%
  save_plot("cellline_wgd_procan.png")
df_depmap %>%
  signif_violin_plot(WGD, Buffering.CellLine.Ratio,
                     test = t.test) %>%
  save_plot("cellline_wgd_depmap.png")
df_depmap %>%
  filter(CellLine.WGD == 1 | CellLine.WGD == 2) %>%
  mutate(CellLine.WGD = factor(CellLine.WGD, levels = c(1,2))) %>%
  signif_violin_plot(CellLine.WGD, Buffering.CellLine.Ratio,
                     test = t.test) %>%
  save_plot("cellline_wgd_levels_depmap.png")

## Near-tetraploid cell lines
df_procan %>%
  signif_violin_plot(`Near-Tetraploid`, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_tetraploidy_procan.png")
df_depmap %>%
  signif_violin_plot(`Near-Tetraploid`, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_tetraploidy_depmap.png")

## High vs. Low Aneuploidy Score
test <- df_procan %>%
  filter(Aneuploidy == "High" | Aneuploidy == "Low") %>%
  signif_violin_plot(Aneuploidy, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_aneuploidy-class_procan.png")
df_depmap %>%
  filter(Aneuploidy == "High" | Aneuploidy == "Low") %>%
  signif_violin_plot(Aneuploidy, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_aneuploidy-class_depmap.png")

# Regression analysis
## Age
df_procan %>%
  scatter_plot_regression(age_at_sampling, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ age_at_sampling,
                          label_coords = c(50, -0.4), point_size = 1) %>%
  save_plot("cellline_age_procan.png")
df_depmap %>%
  scatter_plot_regression(Age, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ Age,
                          label_coords = c(50, -0.2), point_size = 1) %>%
  save_plot("cellline_age_depmap.png")

## Ploidy
df_procan %>%
  scatter_plot_regression(ploidy_wes, Buffering.CellLine.Ratio,  Buffering.CellLine.Ratio ~ ploidy_wes,
                          label_coords = c(3.5, 0.7), point_size = 1) %>%
  save_plot("cellline_ploidy_procan.png")
df_depmap %>%
  scatter_plot_regression(CellLine.Ploidy, Buffering.CellLine.Ratio,  Buffering.CellLine.Ratio ~ CellLine.Ploidy,
                          label_coords = c(5, 0.2), point_size = 1) %>%
  save_plot("cellline_ploidy_depmap.png")

## Mutational burden
df_procan %>%
  scatter_plot_regression(mutational_burden, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ mutational_burden,
                          label_coords = c(400, -0.4), point_size = 1) %>%
  save_plot("cellline_mutations_procan.png")

## Aneuploidy score
df_procan %>%
  scatter_plot_regression(CellLine.AneuploidyScore, Buffering.CellLine.Ratio,
                          Buffering.CellLine.Ratio ~ CellLine.AneuploidyScore,
                          label_coords = c(25, 0.8), point_size = 1) %>%
  save_plot("cellline_aneuploidy_procan.png")
df_depmap %>%
  scatter_plot_regression(CellLine.AneuploidyScore, Buffering.CellLine.Ratio,
                          Buffering.CellLine.Ratio ~ CellLine.AneuploidyScore,
                          label_coords = c(25, 0.3), point_size = 1) %>%
  save_plot("cellline_aneuploidy_depmap.png")

## Misc
df_procan %>%
  scatter_plot_regression(age_at_sampling, ploidy_wes, ploidy_wes ~ age_at_sampling,
                          label_coords = c(40, 1), point_size = 1)

### Ploidy measurements differ across datasets
df_procan %>%
  scatter_plot_regression(ploidy_wes, CellLine.Ploidy, CellLine.Ploidy ~ ploidy_wes,
                          label_coords = c(3, 6), point_size = 1)
