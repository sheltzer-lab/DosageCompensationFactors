library(here)
library(dplyr)
library(tidyr)
library(arrow)
library(ggplot2)
library(forcats)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "FactorAnalysis", "Multivariate")
dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

## SHAP Results
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

# Factor-SHAP correlation heatmaps
heatmap_shap_depmap <- shap_results %>%
  filter(Model.Dataset == "DepMap") %>%
  filter(!grepl("Log2FC", Model.Variant) & !grepl("Average", Model.Variant)) %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, abs(SHAP.Factor.Corr), .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Model.Variant, Model.Subset,
           Model.Level, Model.Condition, SHAP.Factor.Corr, Label) %>%
  nested_shap_heatmap(as.formula(Model.Subset + Model.Level + Model.Condition ~ .)) %>%
  save_plot("shap_heatmap_depmap.png", width = 300)

heatmap_shap_wgd <- shap_results %>%
  filter(Model.Subset %in% c("WGD", "Non-WGD") & Model.Condition == "Loss") %>%
  filter(!grepl("Log2FC", Model.Variant) & !grepl("Average", Model.Variant)) %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, abs(SHAP.Factor.Corr), .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Model.Variant, Model.Subset,
           Model.Level, Model.Condition, SHAP.Factor.Corr, Label) %>%
  nested_shap_heatmap(as.formula(Model.Subset + Model.Level + Model.Condition ~ .)) %>%
  save_plot("shap_heatmap_wgd.png", width = 300, height = 150)

# TODO: Visualize change between WGD and Non-WGD
## Difference SHAP, t-test, normalize SHAP
## Difference Correlation,
## Correlation of Correlation across methods
# TODO: Random choice of SHAP observations could be confounded by Gene and Cell Line -> sample for each cell line

## DepMap vs. ProCan vs. CPTAC
shap_heatmap_datasets <- shap_results %>%
  filter(Model.Level == "Gene" & Model.Subset == "All") %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, -SHAP.Factor.Corr, .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Model.Dataset, Model.Condition, Model.Variant, SHAP.Factor.Corr, Label) %>%
  nested_shap_heatmap(as.formula(Model.Dataset + Model.Condition ~ .)) %>%
  save_plot("shap_heatmap_datasets.png", width = 280, height = 150)

# === Test Area ===
## SHAP-Value Correlation Matrix
### All conditions
shap_results %>%
  select(ID, Model.Variant, SHAP.Value, DosageCompensation.Factor) %>%
  pivot_wider(values_from = SHAP.Value, names_from = DosageCompensation.Factor, id_cols = c(ID, Model.Variant)) %>%
  select(everything(), -ID, -Model.Variant) %>%
  plot_correlation(adjust = "BY")

### Gain
shap_results %>%
  filter(Model.Condition == "Gain") %>%
  select(ID, Model.Variant, SHAP.Value, DosageCompensation.Factor) %>%
  pivot_wider(values_from = SHAP.Value, names_from = DosageCompensation.Factor, id_cols = c(ID, Model.Variant)) %>%
  select(everything(), -ID, -Model.Variant) %>%
  plot_correlation(adjust = "BY")

### Loss
shap_results %>%
  filter(Model.Condition == "Loss") %>%
  select(ID, Model.Variant, SHAP.Value, DosageCompensation.Factor) %>%
  pivot_wider(values_from = SHAP.Value, names_from = DosageCompensation.Factor, id_cols = c(ID, Model.Variant)) %>%
  select(everything(), -ID, -Model.Variant) %>%
  plot_correlation(adjust = "BY")

## SHAP Confounder Plots
shap_confounder <- shap_results %>%
  filter(Model.Variant == "CPTAC_Gene-Level_Filtered_Gain") %>%
  select(ID, DosageCompensation.Factor.ID, Factor.Value.Relative, SHAP.Value) %>%
  filter(DosageCompensation.Factor.ID %in% c("protein_complexes_corum", "mean_gene_dependency")) %>%
  pivot_wider(names_from = DosageCompensation.Factor.ID,
              values_from = c(Factor.Value.Relative, SHAP.Value),
              id_cols = ID)

### Complexes confound Gene Dependency
shap_confounder %>%
  arrange(Factor.Value.Relative_protein_complexes_corum) %>%
  ggplot() +
  aes(x = Factor.Value.Relative_mean_gene_dependency,
      y = SHAP.Value_mean_gene_dependency,
      color = Factor.Value.Relative_protein_complexes_corum) +
    geom_point() +
    scale_color_viridis(option = "C", trans = "log")


shap_confounder %>%
  arrange(Factor.Value.Relative_protein_complexes_corum) %>%
  ggplot() +
  aes(x = SHAP.Value_protein_complexes_corum,
      y = SHAP.Value_mean_gene_dependency,
      color = Factor.Value.Relative_protein_complexes_corum) +
    geom_point() +
    scale_color_viridis(option = "C", trans = "log")

### Gene Dependency confounds Complexes
shap_confounder %>%
  arrange(Factor.Value.Relative_mean_gene_dependency) %>%
  ggplot() +
  aes(x = Factor.Value.Relative_protein_complexes_corum,
      y = SHAP.Value_protein_complexes_corum,
      color = Factor.Value.Relative_mean_gene_dependency) +
    geom_point() +
    scale_color_viridis(option = "C", trans = "log")


shap_confounder %>%
  arrange(Factor.Value.Relative_mean_gene_dependency) %>%
  ggplot() +
  aes(x = SHAP.Value_mean_gene_dependency,
      y = SHAP.Value_protein_complexes_corum,
      color = Factor.Value.Relative_mean_gene_dependency) +
    geom_point() +
    scale_color_viridis(option = "C", trans = "log")
