library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(msigdbr)

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

# === Evaluate GSEA of UPR Gene Set ===
diff_exp_cptac <- read_parquet(here(output_data_dir, "model_buf_diff-exp_cptac.parquet"))
diff_exp_procan <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))
diff_exp_depmap <- read_parquet(here(output_data_dir, "model_buf_diff-exp_depmap.parquet"))
diff_exp_control <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan_adherent.parquet"))

upr_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>%
  rename(Gene.Symbol = "gene_symbol") %>%
  filter(gs_name == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")

upr_de_procan <- diff_exp_procan %>%
  semi_join(y = upr_gene_set, by = "Gene.Symbol")

upr_de_cptac <- diff_exp_cptac %>%
  semi_join(y = upr_gene_set, by = "Gene.Symbol")

## Perform ORA
upr_down_common <- semi_join(
  x = upr_de_procan %>% filter(Log2FC < 0),
  y = upr_de_cptac %>% filter(Log2FC < 0),
  by = "Gene.Symbol"
)

ora_down_common <- upr_down_common %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

upr_up_common <- semi_join(
  x = upr_de_procan %>% filter(Log2FC > 0),
  y = upr_de_cptac %>% filter(Log2FC > 0),
  by = "Gene.Symbol"
)

ora_up_common <- upr_up_common %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

upr_tumor_down_cellline_up <- semi_join(
  x = upr_de_procan %>% filter(Log2FC > 0),
  y = upr_de_cptac %>% filter(Log2FC < 0),
  by = "Gene.Symbol"
)

ora_tumor_down_cellline_up <- upr_tumor_down_cellline_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

upr_tumor_up_cellline_down <- semi_join(
  x = upr_de_procan %>% filter(Log2FC < 0),
  y = upr_de_cptac %>% filter(Log2FC > 0),
  by = "Gene.Symbol"
)

ora_tumor_up_cellline_down <- upr_tumor_up_cellline_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

## Select unique terms per group
ora_down_common_unique <- ora_down_common$result %>%
  anti_join(y = bind_rows(ora_up_common$result, ora_tumor_down_cellline_up$result, ora_tumor_up_cellline_down$result),
            by = "term_id")

ora_up_common_unique <- ora_up_common$result %>%
  anti_join(y = bind_rows(ora_down_common$result, ora_tumor_down_cellline_up$result, ora_tumor_up_cellline_down$result),
            by = "term_id")

ora_tumor_down_cellline_up_unique <- ora_tumor_down_cellline_up$result %>%
  anti_join(y = bind_rows(ora_down_common$result, ora_up_common$result, ora_tumor_up_cellline_down$result),
            by = "term_id")

ora_tumor_up_cellline_down_unique <- ora_tumor_up_cellline_down$result %>%
  anti_join(y = bind_rows(ora_down_common$result, ora_up_common$result, ora_tumor_down_cellline_up$result),
            by = "term_id")

## Plot ORA results
list(result = ora_down_common_unique) %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "REAC"), string_trunc = 70,
                     custom_color = color_palettes$DiffExp["Down"]) %>%
  save_plot("ora_upr_common_down.png", width = 150, height = 150)

list(result = ora_up_common_unique) %>%
  plot_terms_compact(custom_color = color_palettes$DiffExp["Up"], string_trunc = 70) %>%
  save_plot("ora_upr_common_up.png", width = 150, height = 150)

list(result = ora_tumor_down_cellline_up_unique) %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "REAC"), string_trunc = 70,
                     custom_color = color_palettes$Datasets["CPTAC"]) %>%
  save_plot("ora_upr_tumor_down_cellline_up.png", width = 150, height = 150)

list(result = ora_tumor_up_cellline_down_unique) %>%
  plot_terms_compact(custom_color = color_palettes$Datasets["ProCan"], string_trunc = 70) %>%
  save_plot("ora_upr_tumor_up_cellline_down.png", width = 150, height = 150)
