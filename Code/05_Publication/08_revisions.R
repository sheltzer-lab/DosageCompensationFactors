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
  save_plot("ora_upr_common_down.png")

list(result = ora_up_common_unique) %>%
  plot_terms_compact(custom_color = color_palettes$DiffExp["Up"], string_trunc = 70) %>%
  save_plot("ora_upr_common_up.png")

list(result = ora_tumor_down_cellline_up_unique) %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "REAC"), string_trunc = 70,
                     custom_color = color_palettes$Datasets["CPTAC"]) %>%
  save_plot("ora_upr_tumor_down_cellline_up.png")

list(result = ora_tumor_up_cellline_down_unique) %>%
  plot_terms_compact(custom_color = color_palettes$Datasets["ProCan"], string_trunc = 70) %>%
  save_plot("ora_upr_tumor_up_cellline_down.png")

# === BR-cutoff sensitivity analysis ===
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))

expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
model_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))

## Get quantile based on z-score cutoff
z_cutoff <- abs(qnorm(0.2))
model_buf_quantiles <- bind_rows(model_buf_depmap, model_buf_procan, model_buf_cptac) %>%
  mutate(Model.Buffering.Ratio.Quantile = ecdf(Model.Buffering.Ratio.ZScore)(Model.Buffering.Ratio.ZScore),
         .by = Dataset)

quantile_dist <- model_buf_quantiles %>%
  ggplot() +
  aes(x = Model.Buffering.Ratio.ZScore, y = Model.Buffering.Ratio.Quantile, color = Dataset) +
  geom_line() +
  geom_vline(xintercept = -z_cutoff) +
  geom_vline(xintercept = z_cutoff) +
  scale_x_continuous(breaks = seq(floor(min(model_buf_quantiles$Model.Buffering.Ratio.ZScore)),
                                  ceiling(max(model_buf_quantiles$Model.Buffering.Ratio.ZScore)),
                                  1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_color_manual(values = color_palettes$Datasets) +
  labs(x = "sampleBR (z-score)", y = "sampleBR Percentiles")

save_plot(quantile_dist, "sampleBR_ZScore_Quantiles.png", width = 200)

## Grid search of cutoffs
cutoff_steps <- seq(5, 50, 5)
df_grid <- expand_grid(cutoff_steps, cutoff_steps)

eval_datasets <- list(
  DepMap = list(ModelBR = model_buf_depmap, Proteome = expr_buf_depmap),
  ProCan = list(ModelBR = model_buf_procan, Proteome = expr_buf_procan),
  CPTAC = list(ModelBR = model_buf_cptac, Proteome = expr_buf_cptac)
)

results <- list()
pb <- txtProgressBar(min = 1, max = nrow(df_grid) * 3, style = 3)
for (eval_dataset in names(eval_datasets)) {
 model_buf_current <- eval_datasets[[eval_dataset]]$ModelBR
 expr_buf_current <- eval_datasets[[eval_dataset]]$Proteome

 for (i in seq_len(nrow(df_grid))) {
  cutoff_low <- df_grid[[i, 1]]
  cutoff_high <- 100 - df_grid[[i, 2]]

  suppressWarnings({
   results[[paste0(eval_dataset, i)]] <- model_buf_current %>%
     split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                        quantile_low = paste0(cutoff_low, "%"), quantile_high = paste0(cutoff_high, "%")) %>%
     inner_join(y = expr_buf_current, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
     select(Model.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
     differential_expression(Gene.Symbol, Model.Buffering.Group, Protein.Expression.Normalized,
                             groups = c("Low", "High")) %>%
     # Evaluation metrics
     mutate(Group.Balance = Count_GroupA / (Count_GroupA + Count_GroupB),
            Group.Balance.Norm = 1 - abs(0.5 - Group.Balance) * 2) %>%
     summarize(Dataset = eval_dataset,
               Low = cutoff_low,
               High = cutoff_high,
               Log2FC.Abs.Max = max(abs(Log2FC)),
               p.adj.Max = max(-log10(Test.p.adj)),
               Observations.Min = min(Count_GroupA, Count_GroupB),
               Group.Balance.Min = min(Group.Balance.Norm),
               Significant.Count = sum(!is.na(Significant)),
               Significant.Genes = list(Gene.Symbol[!is.na(Significant)]))
  })
  setTxtProgressBar(pb, pb$getVal() + 1)
 }
}
close(pb)

df_grid_results <- bind_rows(results) %>%
  write_parquet(here(output_data_dir, "diffexp_sensitivity.parquet"))

jaccard_index <- function(set_list) {
  set_list <- lapply(set_list, unique)
  intersection <- length(Reduce(intersect, set_list))
  union <- length(unique(unlist(set_list)))
  return(intersection / union)
}

df_grid_results <- read_parquet(here(output_data_dir, "diffexp_sensitivity.parquet"))

## Number of Significant Hits
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Significant.Count) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

## Maximum Absolute Log2FC
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Log2FC.Abs.Max) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

## Maximum -log10 adjusted p value
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = p.adj.Max) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

## Minimum number of observation
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Observations.Min) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

## Minimum group balance per parameter across genes
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Group.Balance.Min) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

## Inter-dataset robustness
df_robustness_inter <- df_grid_results %>%
  filter(Dataset != "CPTAC") %>%
  summarize(GeneLists = list(Significant.Genes),
            .by = c(Low, High)) %>%
  mutate(Robustness.InterDataset = purrr::map_dbl(GeneLists, jaccard_index))

df_robustness_inter %>%
  ggplot() +
  aes(x = Low, y = High, fill = Robustness.InterDataset) +
  geom_tile() +
  scale_fill_viridis()

## Intra-data robustness
df_robustness_intra <- df_grid_results %>%
  inner_join(df_grid_results, by = "Dataset", suffix = c("", ".neighbor")) %>%
  filter(abs(Low - Low.neighbor) <= 5 & abs(High - High.neighbor) <= 5) %>%
  summarize(GeneLists = list(Significant.Genes.neighbor),
            .by = c(Low, High, Dataset)) %>%
  mutate(Robustness.IntraDataset = purrr::map_dbl(GeneLists, jaccard_index),
         Robustness.IntraDataset = if_else(is.nan(Robustness.IntraDataset), 0, Robustness.IntraDataset))

df_robustness_intra %>%
  filter(Dataset == "ProCan") %>%
  ggplot() +
  aes(x = Low, y = High, fill = Robustness.IntraDataset) +
  geom_tile() +
  scale_fill_viridis()

## Joint analysis
normalize_min_max <- function(x) (x - min(x)) / (max(x) - min(x))
norm_cols <- c("Log2FC.Abs.Max", "p.adj.Max", "Significant.Count", "Observations.Min",
               "Robustness.InterDataset", "Robustness.IntraDataset", "Group.Balance.Min")

df_sensitivity <- df_grid_results %>%
  inner_join(df_robustness_intra %>% select(-GeneLists), by = c("High", "Low", "Dataset")) %>%
  inner_join(df_robustness_inter %>% select(-GeneLists), by = c("High", "Low")) %>%
  mutate(across(norm_cols, normalize_min_max, .names = "{.col}_norm")) %>%
  mutate(SensitivityScore = (
    0.30 * Significant.Count_norm +    # prioritize number of hits
    0.20 * Robustness.InterDataset +   # inter-dataset consistency
    0.10 * Robustness.IntraDataset +   # local robustness
    0.10 * p.adj.Max_norm +            # strength of statistical evidence
    0.05 * Log2FC.Abs.Max_norm +       # effect size magnitude
    0.10 * Group.Balance.Min +         # balance of observations between compared groups
    0.15 * Observations.Min_norm       # minimum number of observations in either group
  )) %>%
  mutate(Penalty = if_else(Significant.Count < 10 | (Robustness.InterDataset < 0.01 & Observations.Min < 4 & Group.Balance.Min < 0.3), 0.05, 0),
         PenalizedSensitivity = SensitivityScore - Penalty) %>%
  write_parquet(here(output_data_dir, "diffexp_sensitivity_scores.parquet"))

df_sensitivity %>%
  ggplot() +
  aes(x = Low, y = High, fill = PenalizedSensitivity) +
  geom_tile() +
  scale_fill_viridis(limits = c(0, 0.4), oob = scales::squish) +
  facet_wrap(~Dataset)

df_sensitivity %>%
  ggplot() +
  aes(x = Low, y = High, fill = Significant.Count_norm) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

df_sensitivity %>%
  ggplot() +
  aes(x = Low, y = High, fill = p.adj.Max_norm) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~Dataset)

df_sensitivity %>%
  summarize(SensitivityScore.Median = median(SensitivityScore), .by = c(Low, High)) %>%
  ggplot() +
  aes(x = Low, y = High, fill = SensitivityScore.Median) +
  geom_tile() +
  scale_fill_viridis()

df_sensitivity %>%
  summarize(SensitivityScore.SD = sd(SensitivityScore), .by = c(Low, High)) %>%
  ggplot() +
  aes(x = Low, y = High, fill = -log10(SensitivityScore.SD)) +
  geom_tile() +
  scale_fill_viridis()
