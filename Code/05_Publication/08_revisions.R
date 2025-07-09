library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication", "Revisions")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

# === Load Data ===
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

# === Cluster SHAP values ===
cluster_shap <- function(df, value_col, n_clusters = 3, metric = "pearson", linkage = "average") {
  require(amap)

  df %>%
    select(ID, DosageCompensation.Factor, { { value_col } }) %>%
    pivot_wider(names_from = ID, values_from = { { value_col } }, id_cols = DosageCompensation.Factor) %>%
    tibble::column_to_rownames("DosageCompensation.Factor") %>%
    amap::Dist(method = metric) %>%
    hclust(method = linkage) %>%
    cutree(k = n_clusters) %>%
    as.data.frame() %>%
    rename(Cluster = 1) %>%
    mutate(Cluster = as.character(Cluster)) %>%
    tibble::rownames_to_column("DosageCompensation.Factor")
}

## All (pan-cancer) Models
shap_pan_cancer <- shap_results %>%
  filter(Model.BufferingMethod == "BR" & Model.Samples == "Unaveraged") %>%
  filter(Model.Dataset != "Engineered") %>%
  filter(Model.Subset == "All") %>%
  #filter(Model.Level == "Gene") %>%
  #filter(Model.Condition == "Gain") %>%
  filter(Model.ROC.AUC > 0.65) %>%
  unite("ID", ID, Model.Dataset, Model.Subset, Model.Condition, Model.Level) %>%
  mutate(SHAP.Value.Scaled = as.vector(scale(SHAP.Value)), .by = "Model.Filename")

clusters_pan_cancer <- shap_pan_cancer %>%
  cluster_shap(SHAP.Value.Scaled, n_clusters = 3, metric = "pearson", linkage = "average")

pca_shap <- shap_pan_cancer %>%
  left_join(y = clusters_pan_cancer, by = "DosageCompensation.Factor") %>%
  calculate_pca(DosageCompensation.Factor, NULL, ID, SHAP.Value.Scaled) %>%
  plot_pca(color_col = Cluster, label_col = DosageCompensation.Factor) %>%
  save_plot("shap_pca.png", width = 250, height = 250)

## SHAP per dataset
datasets <- c("DepMap", "ProCan", "CPTAC")
for (dataset in datasets) {
  shap_dataset <- shap_results %>%
    filter(Model.BufferingMethod == "BR" & Model.Samples == "Unaveraged") %>%
    filter(Model.Dataset == dataset) %>%
    filter(Model.Subset == "All") %>%
    filter(Model.Level == "Gene") %>%
    mutate(SHAP.Value.Scaled = as.vector(scale(SHAP.Value)), .by = "Model.Filename")

  clusters_dataset <- shap_dataset %>%
    cluster_shap(SHAP.Value.Scaled, n_clusters = 3, metric = "pearson", linkage = "average")

  pca_shap_dataset <-  shap_dataset %>%
    left_join(y = clusters_dataset, by = "DosageCompensation.Factor") %>%
    calculate_pca(DosageCompensation.Factor, NULL, ID, SHAP.Value.Scaled) %>%
    plot_pca(color_col = Cluster, label_col = DosageCompensation.Factor) %>%
    save_plot(paste0("shap_pca_", tolower(dataset), ".png"), width = 250, height = 250)
}

## SHAP per condition
conditions <- c("Gain", "Loss")

for (condition in conditions) {
  shap_cn <- shap_results %>%
    filter(Model.BufferingMethod == "BR" & Model.Samples == "Unaveraged") %>%
    filter(Model.Dataset != "Engineered") %>%
    filter(Model.Subset == "All") %>%
    filter(Model.Level == "Gene") %>%
    filter(Model.Condition == condition) %>%
    mutate(SHAP.Value.Scaled = as.vector(scale(SHAP.Value, scale = FALSE)), .by = "Model.Filename")

  clusters_cn <- shap_cn %>%
    cluster_shap(SHAP.Value.Scaled, n_clusters = 3, metric = "pearson", linkage = "average")

  pca_shap_cn <- shap_cn %>%
    left_join(y = clusters_cn, by = "DosageCompensation.Factor") %>%
    calculate_pca(DosageCompensation.Factor, NULL, ID, SHAP.Value.Scaled) %>%
    plot_pca(color_col = Cluster, label_col = DosageCompensation.Factor) %>%
    save_plot(paste0("shap_pca_", tolower(condition), ".png"), width = 250, height = 250)
}

## Single Model
## Remark: PCA of single models results in one PC with >90% explained variance.
# TODO: Investigate PC with large coverage
for (condition in conditions) {
  shap_cptac <- shap_results %>%
    filter(Model.BufferingMethod == "BR" & Model.Samples == "Unaveraged") %>%
    filter(Model.Dataset == "CPTAC") %>%
    filter(Model.Subset == "All") %>%
    filter(Model.Level == "Gene") %>%
    filter(Model.Condition == condition)

  clusters_cptac <- shap_cptac %>%
    cluster_shap(SHAP.Value, n_clusters = 3, metric = "pearson", linkage = "average")

  pca_shap_cptac <- shap_cptac %>%
    left_join(y = clusters_cptac, by = "DosageCompensation.Factor") %>%
    calculate_pca(DosageCompensation.Factor, NULL, ID, SHAP.Value) %>%
    plot_pca(color_col = Cluster, label_col = DosageCompensation.Factor) %>%
    save_plot(paste0("shap_pca_cptac_", tolower(condition), ".png"), width = 250, height = 250)
}

## Try UMAP
umap_shap <- shap_pan_cancer %>%
  left_join(y = clusters_pan_cancer, by = "DosageCompensation.Factor") %>%
  calculate_umap(DosageCompensation.Factor, NULL, ID, SHAP.Value.Scaled) %>%
  plot_umap(color_col = Cluster, label_col = DosageCompensation.Factor) %>%
  save_plot("shap_umap.png", width = 250, height = 250)
