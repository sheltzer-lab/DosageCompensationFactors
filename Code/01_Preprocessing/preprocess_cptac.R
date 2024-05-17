library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "annotation.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))

expression_data_dir <- here(external_data_dir, "Expression", "CPTAC", "Proteome_BCM_GENCODE_v34_harmonized_v1")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Preprocessing", "Expression", "CPTAC")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
file_list <- list.files(expression_data_dir, include.dirs = FALSE, recursive = FALSE)
file_list <- file_list[grep(".+_proteomics_.+.txt", file_list)]

df_list <- list()

for (filename in file_list) {
  metadata <- sub('\\.txt$', '', filename) %>%
    str_split_fixed("_", 9) %>%
    as.data.frame() %>%
    select(V1, V9)

  df_name <- paste(metadata$V1, metadata$V9, sep = "_")

  df_list[[df_name]] <- read.table(here(expression_data_dir, filename),
                                   sep = "\t", dec = ".", header = FALSE, stringsAsFactors = FALSE) %>%
    janitor::row_to_names(1) %>%
    pivot_longer(everything() & !idx,
                 names_to = "Model.ID", values_to = "Protein.Expression.Log2",
                 names_ptypes = character(), values_ptypes = double(),
                 names_transform = as.character, values_transform = as.numeric) %>%
    mutate(Model.CancerType = metadata$V1,
           Model.SampleType = metadata$V9)
}

# === Combine & Tidy Datasets ===
expr_cptac <- bind_rows(df_list) %>%
  unite("Model.SampleID", c("Model.ID", "Model.SampleType"), sep = '_', remove = FALSE)

# === Preprocess Datasets ===
# TODO: Reprodcibility Score filtering

expr_cptac_processed <- expr_cptac %>%
  remove_noisefloor(Protein.Expression.Log2) %>%
  normalize_samples(Model.SampleID, Protein.Expression.Log2, idx,
                    normalized_colname = "Protein.Expression.Normalized")

# === Annotation ===
# TODO: ID Mapping

# === Evaluation ===
expr_dist <- plot_expr_dist(expr_cptac_processed) %>%
  save_plot("cptac_expression_distribution.png")

sample_subset <- expr_cptac_processed %>%
  distinct(Model.SampleID, Model.ID, Model.CancerType, Model.SampleType) %>%
  group_by(Model.CancerType, Model.SampleType) %>%
  slice_sample(n = 2)

sample_dist_pre <- expr_cptac_processed %>%
  filter(Model.SampleID %in% sample_subset$Model.SampleID) %>%
  mutate(Model.SampleID = fct_reorder(Model.SampleID, Model.CancerType)) %>%
  vertical_box_plot(Protein.Expression.Log2, Model.SampleID, Model.CancerType) %>%
  save_plot("cptac_sample_distribution_non-norm.png")

sample_dist_post <- expr_cptac_processed %>%
  filter(Model.SampleID %in% sample_subset$Model.SampleID) %>%
  mutate(Model.SampleID = fct_reorder(Model.SampleID, Model.CancerType)) %>%
  vertical_box_plot(Protein.Expression.Normalized, Model.SampleID, Model.CancerType) %>%
  save_plot("cptac_sample_distribution_norm.png")

types_dist_pre <- expr_cptac_processed %>%
  unite("Model.CancerSampleType", c(Model.CancerType, Model.SampleType), sep = '_', remove = FALSE) %>%
  vertical_box_plot(Protein.Expression.Log2, Model.CancerSampleType, Model.CancerType) %>%
  save_plot("cptac_cancer-type_distribution_non-norm.png")

types_dist_post <- expr_cptac_processed %>%
  unite("Model.CancerSampleType", c(Model.CancerType, Model.SampleType), sep = '_', remove = FALSE) %>%
  vertical_box_plot(Protein.Expression.Normalized, Model.CancerSampleType, Model.CancerType) %>%
  save_plot("cptac_cancer-type_distribution_norm.png")

pca_pre <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, idx, Protein.Expression.Log2) %>%
  plot_pca(color_col = Model.CancerType, label_col = NULL) %>%
  save_plot("cptac_pca_non-norm.png")

pca_pre_sample <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, idx, Protein.Expression.Log2) %>%
  plot_pca(color_col = Model.SampleType, label_col = NULL) %>%
  save_plot("cptac_pca_sample_non-norm.png")

pca_norm <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, idx, Protein.Expression.Normalized) %>%
  plot_pca(color_col = Model.CancerType, label_col = NULL) %>%
  save_plot("cptac_pca_norm.png")

pca_norm <- expr_cptac_processed %>%
  calculate_pca(Model.SampleID, Model.CancerType, idx, Protein.Expression.Normalized) %>%
  plot_pca(color_col = Model.SampleType, label_col = NULL) %>%
  save_plot("cptac_pca_sample_norm.png")
