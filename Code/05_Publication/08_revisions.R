library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(msigdbr)
library(openxlsx)
library(readxl)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "annotation.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "02_DosageCompensation", "baseline.R"))

plots_dir <- here(plots_base_dir, "Publication", "Revisions")
tables_dir <- here(tables_base_dir, "Publication", "Revisions")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

# === Load Data ===
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

diff_exp_cptac <- read_parquet(here(output_data_dir, "model_buf_diff-exp_cptac.parquet"))
diff_exp_procan <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))
diff_exp_depmap <- read_parquet(here(output_data_dir, "model_buf_diff-exp_depmap.parquet"))
diff_exp_control <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan_adherent.parquet"))

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))

tumor_types <- mskcc.oncotree::get_tumor_types()

df_model_depmap <- read_csv_arrow(here(external_data_dir, "CopyNumber", "DepMap", "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId") %>%
  select(-CellLine.Name) %>%
  mutate(Model.CancerType = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 3)))

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet')) %>%
  # Oncotree uses different cancer type abbrieviations than CPTAC
  mutate(Model.CancerType = str_replace_all(Model.CancerType, c(HNSCC = "HNSC", LSCC = "LUSC", PDAC = "PAAD"))) %>%
  mutate(Model.CancerType = sapply(Model.CancerType, \(x) get_oncotree_parent(tumor_types, x, target_level = 3)))

model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet")) %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, Model.CancerType), by = "Model.ID")
model_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet")) %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, Model.CancerType), by = "Model.ID")
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet")) %>%
  inner_join(y = df_model_cptac %>% select(Model.ID, Model.CancerType), by = "Model.ID")

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
upr_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>%
  rename(Gene.Symbol = "gene_symbol") %>%
  filter(gs_name == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")

upr_de_procan <- diff_exp_procan %>%
  semi_join(y = upr_gene_set, by = "Gene.Symbol")

upr_de_cptac <- diff_exp_cptac %>%
  semi_join(y = upr_gene_set, by = "Gene.Symbol")

## Perform ORA
upr_down_common <- semi_join(
  x = upr_de_procan %>% filter(Log2FC < 0 & Test.p.adj < p_threshold),
  y = upr_de_cptac %>% filter(Log2FC < 0),
  by = "Gene.Symbol"
)

ora_down_common <- upr_down_common %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

upr_up_common <- semi_join(
  x = upr_de_procan %>% filter(Log2FC > 0 & Test.p.adj < p_threshold),
  y = upr_de_cptac %>% filter(Log2FC > 0),
  by = "Gene.Symbol"
)

ora_up_common <- upr_up_common %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

upr_tumor_down_cellline_up <- semi_join(
  x = upr_de_procan %>% filter(Log2FC > 0 & Test.p.adj < p_threshold),
  y = upr_de_cptac %>% filter(Log2FC < 0),
  by = "Gene.Symbol"
)

ora_tumor_down_cellline_up <- upr_tumor_down_cellline_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = FALSE)

upr_tumor_up_cellline_down <- semi_join(
  x = upr_de_procan %>% filter(Log2FC < 0 & Test.p.adj < p_threshold),
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

list(result = ora_tumor_up_cellline_down$result) %>%
  plot_terms_compact(custom_color = color_palettes$Datasets["ProCan"], string_trunc = 70) %>%
  save_plot("ora_upr_tumor_up_cellline_down.png")

## Check for overlap between ORA gene sets in Figure 5 and UPR hallmark gene set
### Check: CCT complex (CORUM), unfolded protein binding (GO:MF), {UGGT1, CANX, TAPBP}
genes_down <- diff_exp_procan %>%
  filter(Significant == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_genes <- ora_down$result %>%
  filter(grepl("CCT", term_name) & source == "CORUM" | grepl("unfolded", term_name) & source == "GO:MF") %>%
  separate_longer_delim(intersection, ",") %>%
  rename(Gene.Symbol = "intersection")

intersect(upr_gene_set$Gene.Symbol, ora_genes$Gene.Symbol)
intersect(upr_gene_set$Gene.Symbol, c("UGGT1", "CANX", "TAPBP"))

### Conclusion: Gene sets identified via ORA and STRING are disjoint to UPR gene set from mSigDB

# === BR-cutoff sensitivity analysis ===
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
      df_split <- model_buf_current %>%
        split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                           quantile_low = paste0(cutoff_low, "%"), quantile_high = paste0(cutoff_high, "%"))

      cancertype_count <- length(unique(df_split$Model.CancerType))

      results[[paste0(eval_dataset, i)]] <- df_split %>%
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
                  Observations.Median = median(c(Count_GroupA, Count_GroupB)),
                  Group.Balance.Min = min(Group.Balance.Norm),
                  Significant.Count = sum(!is.na(Significant)),
                  Significant.Genes = list(Gene.Symbol[!is.na(Significant)]),
                  CancerTypes.Count = cancertype_count)
    })
    setTxtProgressBar(pb, pb$getVal() + 1)
  }
}
close(pb)

bind_rows(results) %>%
  write_parquet(here(output_data_dir, "diffexp_sensitivity.parquet"))

jaccard_index <- function(set_list) {
  set_list <- lapply(set_list, unique)
  intersection <- length(Reduce(intersect, set_list))
  union <- length(unique(unlist(set_list)))
  return(intersection / union)
}

df_grid_results <- read_parquet(here(output_data_dir, "diffexp_sensitivity.parquet"))

## Number of Cancer Types
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = CancerTypes.Count, label = CancerTypes.Count) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  labs(fill = "Cancer Types (Primary Disease)") +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_cancertypes.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Number of Significant Hits
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Significant.Count, label = Significant.Count) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  labs(fill = "Significant Hits") +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_hits.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Maximum Absolute Log2FC
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Log2FC.Abs.Max, label = round(Log2FC.Abs.Max, digits = 1)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  facet_wrap(~Dataset) +
  labs(fill = "max(abs(Log2FC))") +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_Log2FC.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Maximum -log10 adjusted p value
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = p.adj.Max, label = round(p.adj.Max, digits = 1)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  facet_wrap(~Dataset) +
  labs(fill = "max(-log10(p.adj))") +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_p.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Minimum number of observations
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Observations.Min, label = Observations.Min) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  labs(fill = "Min. Obseravtions") +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_observations.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Median number of observation
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Observations.Median, label = round(Observations.Median)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  labs(fill = "Median Obseravtions") +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_observations_median.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Minimum group balance per parameter across genes
df_grid_results %>%
  ggplot() +
  aes(x = Low, y = High, fill = Group.Balance.Min, label = round(Group.Balance.Min, digits = 1)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  labs(fill = "Min. Group Balance") +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_balance.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Inter-dataset robustness
df_robustness_inter <- df_grid_results %>%
  filter(Dataset != "CPTAC") %>%
  summarize(GeneLists = list(Significant.Genes),
            .by = c(Low, High)) %>%
  mutate(Robustness.InterDataset = purrr::map_dbl(GeneLists, jaccard_index))

df_robustness_inter %>%
  ggplot() +
  aes(x = Low, y = High, fill = Robustness.InterDataset, label = round(Robustness.InterDataset, digits = 2)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_robustness_inter.png"), width = 150, height = 180, units = "mm", dpi = 300)

## Intra-data robustness
df_robustness_intra <- df_grid_results %>%
  inner_join(df_grid_results, by = "Dataset", suffix = c("", ".neighbor")) %>%
  filter(abs(Low - Low.neighbor) <= 5 & abs(High - High.neighbor) <= 5) %>%
  summarize(GeneLists = list(Significant.Genes.neighbor),
            .by = c(Low, High, Dataset)) %>%
  mutate(Robustness.IntraDataset = purrr::map_dbl(GeneLists, jaccard_index),
         Robustness.IntraDataset = if_else(is.nan(Robustness.IntraDataset), 0, Robustness.IntraDataset))

df_robustness_intra %>%
  ggplot() +
  aes(x = Low, y = High, fill = Robustness.IntraDataset, label = round(Robustness.IntraDataset, digits = 2)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_robustness_intra.png"), width = 300, height = 150, units = "mm", dpi = 300)

## Joint analysis
replace_nan <- function(x, replacement = 0) replace(x, is.nan(x), replacement)
normalize_min_max <- function(x) replace_nan((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
norm_cols <- c("Log2FC.Abs.Max", "p.adj.Max", "Significant.Count", "Observations.Min", "CancerTypes.Count",
               "Robustness.InterDataset", "Robustness.IntraDataset", "Group.Balance.Min")

df_sensitivity <- df_grid_results %>%
  inner_join(df_robustness_intra %>% select(-GeneLists), by = c("High", "Low", "Dataset")) %>%
  inner_join(df_robustness_inter %>% select(-GeneLists), by = c("High", "Low")) %>%
  mutate(across(all_of(norm_cols), normalize_min_max, .names = "{.col}_norm"), .by = Dataset) %>%
  mutate(SensitivityScore = (
    0.30 * Significant.Count_norm +    # prioritize number of hits
    0.10 * Robustness.InterDataset +   # inter-dataset consistency
    0.05 * Robustness.IntraDataset +   # local robustness
    0.10 * p.adj.Max_norm +            # strength of statistical evidence
    0.05 * Log2FC.Abs.Max_norm +       # effect size magnitude
    0.10 * Group.Balance.Min +         # balance of observations between compared groups
    0.15 * Observations.Min_norm +     # minimum number of observations in either group
    0.15 * CancerTypes.Count_norm      # number of cancer types covered
  )) %>%
  mutate(Penalty = if_else(Significant.Count < 10 | (Robustness.InterDataset < 0.01 & Observations.Min < 4 & Group.Balance.Min < 0.3), 0.05, 0),
         PenalizedSensitivity = SensitivityScore - Penalty) %>%
  write_parquet(here(output_data_dir, "diffexp_sensitivity_scores.parquet"))

df_sensitivity %>%
  select(Dataset, Low, High, all_of(norm_cols), SensitivityScore, Penalty, PenalizedSensitivity) %>%
  write.xlsx(here(tables_dir, "diffexp_sensitivity.xlsx"), asTable = TRUE)

df_sensitivity %>%
  ggplot() +
  aes(x = Low, y = High, fill = SensitivityScore, label = round(SensitivityScore, digits = 2)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis(limits = c(0, 0.6), oob = scales::squish) +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_score.png"), width = 300, height = 150, units = "mm", dpi = 300)

df_sensitivity %>%
  ggplot() +
  aes(x = Low, y = High, fill = PenalizedSensitivity, label = round(PenalizedSensitivity, digits = 2)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis(limits = c(0, 0.5), oob = scales::squish) +
  facet_wrap(~Dataset) +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_score_penalized.png"), width = 300, height = 150, units = "mm", dpi = 300)

df_sensitivity %>%
  summarize(SensitivityScore.Median = median(SensitivityScore), .by = c(Low, High)) %>%
  ggplot() +
  aes(x = Low, y = High, fill = SensitivityScore.Median, label = round(SensitivityScore.Median, digits = 2)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_score_median.png"), width = 150, height = 180, units = "mm", dpi = 300)

df_sensitivity %>%
  summarize(SensitivityScore.SD = sd(SensitivityScore), .by = c(Low, High)) %>%
  ggplot() +
  aes(x = Low, y = High, fill = -log10(SensitivityScore.SD), label = round(-log10(SensitivityScore.SD), digits = 2)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis() +
  theme(legend.position = "top",
        legend.key.width = unit(0.05, 'npc'))

ggsave(here(plots_dir, "sensitivity_score_sd.png"), width = 150, height = 180, units = "mm", dpi = 300)

# === Multiple Myeloma ===
## Identify multiple myeloma cell lines in DepMap CCLE
mm_celllines <- read_excel(here(external_data_dir, "MultipleMyeloma_CellLines.xlsx"),
                           sheet = "Cell Line Patient Correlations") %>%
  mutate(StrippedCellLineName = toupper(str_replace(`Cell Line`, "_Sus", ""))) %>%
  mutate(StrippedCellLineName = str_replace(StrippedCellLineName, "U266", "U266B1")) # Synonymous accordign to Cellosaurus

df_model_mm <- df_model_depmap %>%
  filter(OncotreeCode == "PCM") %>%
  semi_join(mm_celllines, by = "StrippedCellLineName")

## Analyze sample BR of multiple myeloma
df_depmap_mm <- model_buf_depmap %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, OncotreeCode), by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(Model.CancerType = if_else(OncotreeCode == "PCM", "Multiple\nMyeloma", "Other"),
         Dataset = "DepMap")

df_procan_mm <- model_buf_procan %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, OncotreeCode), by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(Model.CancerType = if_else(OncotreeCode == "PCM", "Multiple\nMyeloma", "Other"),
         Dataset = "ProCan")

### MM against other cancer types
plot_mm_br <- bind_rows(df_depmap_mm, df_procan_mbn) %>%
  signif_beeswarm_plot(Model.CancerType, Model.Buffering.Ratio,
                       color_col = CellLine.AneuploidyScore, viridis_color_pal = color_palettes$AneuploidyScore,
                       color_lims = c(0, 0.8), cex = 1, test = wilcox.test) +
  facet_wrap(~Dataset) +
  labs(x = "Cancer Type", y = "Sample Buffering Ratio", color = "Aneuploidy\nScore")

save_plot(plot_mm_br, "buffering_mm.png")

### MM BR per-protein
df_depmap_mm_prot <- expr_buf_depmap %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, OncotreeCode), by = "Model.ID",
             relationship = "many-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(Model.CancerType = if_else(OncotreeCode == "PCM", "Multiple\nMyeloma", "Other"),
         Dataset = "DepMap")

df_procan_mm_prot <- expr_buf_procan %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, OncotreeCode), by = "Model.ID",
             relationship = "many-to-one", na_matches = "never") %>%
  mutate(Model.CancerType = if_else(OncotreeCode == "PCM", "Multiple\nMyeloma", "Other"),
         Dataset = "ProCan")

plot_mm_br_prot <- bind_rows(df_depmap_mm_prot, df_procan_mm_prot) %>%
  drop_na(Model.CancerType, Buffering.GeneLevel.Ratio) %>%
  signif_boxplot(Model.CancerType, Buffering.GeneLevel.Ratio,
                 test = wilcox.test) +
  facet_wrap(~Dataset) +
  labs(x = "Cancer Type", y = "Buffering Ratio")

save_plot(plot_mm_br_prot, "buffering_mm_protein.png")

### MM by aneuploidy score
plot_mm_br_as <- bind_rows(df_depmap_mm, df_procan_mm) %>%
  filter(OncotreeCode == "PCM") %>%
  ggscatter(
    x = "CellLine.AneuploidyScore", y = "Model.Buffering.Ratio.ZScore",
    color = default_color, size = 3,
    add = "reg.line", add.params = list(color = highlight_colors[2]),
    conf.int = TRUE, cor.coef = TRUE,
    cor.coeff.args = list(method = "spearman", label.sep = "\n", cor.coef.name = "rho")
  ) +
  facet_wrap(~Dataset) +
  labs(x = "Aneuploidy Score", y = "Sample Buffering Ratio (z-score)")

save_plot(plot_mm_br_as, "buffering_mm_aneuploidy.png", width = 200)

mm_ids_depmap <- df_depmap_mm %>%
  filter(OncotreeCode == "PCM") %>%
  distinct(Model.ID) %>%
  pull(Model.ID)

mm_ids_procan <- df_procan_mm %>%
  filter(OncotreeCode == "PCM") %>%
  distinct(Model.ID) %>%
  pull(Model.ID)

length(unique(c(mm_ids_depmap, mm_ids_procan)))
length(mm_ids_depmap)
length(mm_ids_procan)

# === Check difference in Gene CN and ChrArm CN distribution to explain differences in BR ===
## Check CN range
fivenum(expr_buf_procan$Gene.CopyNumber)
fivenum(expr_buf_procan$ChromosomeArm.CopyNumber)

cn_range <- bind_rows(expr_buf_depmap, expr_buf_procan) %>%
  summarize(Gene.CN.Min = min(Gene.CopyNumber, na.rm = TRUE),
            Gene.CN.Max = max(Gene.CopyNumber, na.rm = TRUE),
            ChrArm.CN.Min = min(ChromosomeArm.CopyNumber, na.rm = TRUE),
            ChrArm.CN.Max = max(ChromosomeArm.CopyNumber, na.rm = TRUE),
            .by = Dataset)


## Compare CN distribution
cn_compare <- bind_rows(expr_buf_depmap, expr_buf_procan) %>%
  select(Dataset, Gene.CopyNumber, ChromosomeArm.CopyNumber,
         Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio) %>%
  pivot_longer(c(Gene.CopyNumber, ChromosomeArm.CopyNumber),
               names_to = "Level", values_to = "CopyNumber") %>%
  mutate(Level = str_replace(Level, ".CopyNumber", ""))

### All CNs
cn_compare %>%
  ggplot() +
  aes(x = CopyNumber, fill = Level) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 1000000, 100000)) +
  facet_wrap(~Dataset)

ggsave(here(plots_dir, "cn_compare.png"))

### Gain only
cn_compare %>%
  filter(CopyNumber > 3) %>%
  ggplot() +
  aes(x = CopyNumber, fill = Level) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  facet_wrap(~Dataset)

ggsave(here(plots_dir, "cn_compare_gain.png"))

### All CNs, used for valid BRs
cn_compare %>%
  drop_na() %>% # Removes observations where BR (either Gene- or ChrArm-derived) was not calculated
  ggplot() +
  aes(x = CopyNumber, fill = Level) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 1000000, 100000)) +
  facet_wrap(~Dataset)

ggsave(here(plots_dir, "cn_compare_br.png"))

### Gain only, used for valid BRs
cn_compare %>%
  filter(CopyNumber > 3) %>%
  drop_na() %>% # Removes observations where BR (either Gene- or ChrArm-derived) was not calculated
  ggplot() +
  aes(x = CopyNumber, fill = Level) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  facet_wrap(~Dataset)

ggsave(here(plots_dir, "cn_compare_gain_br.png"))

## Check if mean CN significantly differs
cn_compare %>%
  filter(Dataset == "ProCan") %>%
  t.test(CopyNumber ~ Level, data = .)

cn_compare %>%
  filter(Dataset == "ProCan") %>%
  filter(CopyNumber > 2) %>%
  t.test(CopyNumber ~ Level, data = .)

### BR-only
cn_compare %>%
  filter(Dataset == "ProCan") %>%
  drop_na() %>%
  t.test(CopyNumber ~ Level, data = .)

cn_compare %>%
  filter(Dataset == "ProCan") %>%
  filter(CopyNumber > 2) %>%
  drop_na() %>%
  t.test(CopyNumber ~ Level, data = .)

## Check observations with CN discrepancy
cn_diff <- bind_rows(expr_buf_depmap, expr_buf_procan) %>%
  select(Dataset, Gene.CopyNumber, ChromosomeArm.CopyNumber,
         Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio) %>%
  mutate(CN_Diff = Gene.CopyNumber - ChromosomeArm.CopyNumber)

### All copy numbers
cn_diff %>%
  filter(Dataset == "ProCan") %>%
  ggplot() +
  aes(x = CN_Diff) +
  scale_x_continuous(breaks = seq(-7, 7, 1), limits = c(-8, 8), expand = c(0, 0)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                labels = c("1", "10", "100", "1000", "10000", "100000", "1000000")) +
  geom_bar()

### Copy numbers that were used for BRs
cn_diff %>%
  filter(Dataset == "ProCan") %>%
  drop_na() %>% # Removes observations where BR (either Gene- or ChrArm-derived) was not calculated
  ggplot() +
  aes(x = CN_Diff) +
  scale_x_continuous(breaks = seq(-7, 7, 1), limits = c(-8, 8), expand = c(0, 0)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                labels = c("1", "10", "100", "1000", "10000", "100000", "1000000")) +
  geom_bar()

## Compare CN Ratio against baseline
cn_compare_ratio <- bind_rows(expr_buf_depmap, expr_buf_procan) %>%
  select(Dataset, Gene.CopyNumber, Gene.CopyNumber.Baseline,
         ChromosomeArm.CopyNumber, ChromosomeArm.CopyNumber.Baseline,
         Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio) %>%
  mutate(Gene.CNRatio = Gene.CopyNumber / Gene.CopyNumber.Baseline,
         ChrArm.CNRatio = ChromosomeArm.CopyNumber / ChromosomeArm.CopyNumber.Baseline) %>%
  pivot_longer(c(Gene.CNRatio, ChrArm.CNRatio),
               names_to = "Level", values_to = "CopyNumberRatio") %>%
  mutate(Level = str_replace(Level, ".CNRatio", ""))

cn_compare_ratio %>%
  filter(Dataset == "ProCan") %>%
  ggplot() +
  aes(x = CopyNumberRatio, fill = Level) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = seq(0, 4, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1000000, 100000))

ggsave(here(plots_dir, "cn-ratio_compare.png"))

cn_compare_ratio %>%
  filter(Dataset == "ProCan") %>%
  drop_na() %>%
  ggplot() +
  aes(x = CopyNumberRatio, fill = Level) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = seq(0, 4, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1000000, 10000))

ggsave(here(plots_dir, "cn-ratio_compare_br.png"))

# === Compare Buffering Results with Muenzner et al. ===
chr_map <- setNames(as.character(1:16), LETTERS[1:16])

## Load Datasets
muenzner_strains <- read_excel(here(external_data_dir, "Muenzner2024_Yeast_Tables.xlsx"),
                               sheet = 7, skip = 4) %>%
  rename(strain = "Standardized_name")
muenzner_expr <- read_excel(here(external_data_dir, "Muenzner2024_Yeast_Tables.xlsx"),
                            sheet = 3, skip = 3)
muenzner_dc <- read_excel(here(external_data_dir, "Muenzner2024_Yeast_Tables.xlsx"),
                          sheet = 8, skip = 4)
muenzner_turnover <- read_excel(here(external_data_dir, "Muenzner2024_Yeast_Tables.xlsx"),
                                sheet = 17, skip = 4)
## Tidy Datasets
muenzner_cn_aneuploid <- muenzner_strains %>%
  separate_longer_delim(Aneuploidies, ";") %>%
  filter(Aneuploidies != "aneu" & Aneuploidies != "" & Aneuploidies != "euploid") %>%
  separate_wider_delim(Aneuploidies, "*", names = c("ChromosomeArm.CNA", "Chromosome")) %>%
  mutate_at(c("ChromosomeArm.CNA", "Chromosome"), as.integer)

muenzner_cn_euploid <- muenzner_strains %>%
  select(-Aneuploidies) %>%
  mutate(Chromosome = paste(chr_map, collapse = ";")) %>%
  separate_longer_delim(Chromosome, ";") %>%
  mutate(Chromosome = as.integer(Chromosome),
         ChromosomeArm.CNA = 0L) %>%
  anti_join(y = muenzner_cn_aneuploid, by = c("strain", "Chromosome"))

muenzner_cn <- bind_rows(muenzner_cn_aneuploid, muenzner_cn_euploid) %>%
  mutate(ChromosomeArm.CopyNumber = Ploidy + ChromosomeArm.CNA,
         ChromosomeArm.CopyNumber.Baseline = Ploidy)

muenzner_expr_tidy <- muenzner_expr %>%
  reshape2::melt(value.name = "Protein.Expression", id.vars = "Protein") %>%
  rename(Protein = 1,
         strain = "variable") %>%
  mutate(chr_code = str_sub(Protein, 2, 2),
         Chromosome = as.integer(recode(chr_code, !!!chr_map)),
         Protein.Expression = as.numeric(Protein.Expression),
         Protein.Expression.Log2 = log2(Protein.Expression))

## Calculate BR
muenzner_br <- muenzner_expr_tidy %>%
  inner_join(y = muenzner_cn, by = c("strain", "Chromosome")) %>%
  # Calculate protein baseline as median protein abundance on euploid chromosomes
  # TODO: Consider calculating separate baseline for each ploidy
  calculate_baseline(Protein, ChromosomeArm.CNA, Protein.Expression.Log2,
                     distance_col = Ploidy, target_colname = "Protein.Expression.Baseline",
                     weighted = FALSE, summ_func = median) %>%
  mutate(Buffering.ChrArmLevel.Ratio = buffering_ratio(2^Protein.Expression.Baseline, 2^Protein.Expression.Log2,
                                                       ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber),
         Buffering.ChrArmLevel.Class = buffering_class(Buffering.ChrArmLevel.Ratio,
                                                       2^Protein.Expression.Baseline, 2^Protein.Expression.Log2,
                                                       ChromosomeArm.CopyNumber.Baseline, ChromosomeArm.CopyNumber))

muenzner_br_trimmed <- muenzner_br %>%
  drop_na(Buffering.ChrArmLevel.Ratio)

## Calculate average gene-wise BR
muenzner_br_avg <- muenzner_br_trimmed %>%
  summarize(BR.Avg = mean(Buffering.ChrArmLevel.Ratio),
            Attenuated = BR.Avg > br_cutoffs$Buffered,
            .by = Protein)

gene_vec <- AnnotationDbi::mapIds(
  org.Sc.sgd.db::org.Sc.sgd.db,
  keys      = unique(muenzner_br_trimmed$Protein),
  keytype   = "ORF",
  column    = "GENENAME",
  multiVals = "first"  # or "list" if you want all
)

muenzner_br_avg <- muenzner_br_avg %>%
  mutate(gene_name = if_else(is.na(gene_vec), Protein, gene_vec))

## Compare attenuation calls from Muenzner et al. method with BR
muenzner_dc_trimmed <- muenzner_dc %>%
  select(gene_name, `attenuation slope (protein)`) %>%
  mutate(Attenuated = `attenuation slope (protein)` < 0.85)

### Limit BR-based dataset to genes identified in Muenzner DC table
muenzner_br_avg_trimmed <- muenzner_br_avg %>%
  semi_join(y = muenzner_dc_trimmed, by = "gene_name")

nrow(muenzner_dc_trimmed)
nrow(muenzner_br_avg)
nrow(muenzner_br_avg_trimmed)

frac_atten_muenzner <- sum(muenzner_dc_trimmed$Attenuated) / nrow(muenzner_dc_trimmed)
frac_atten_br <- sum(muenzner_br_avg$Attenuated) / nrow(muenzner_br_avg)
frac_atten_br_trimmed <- sum(muenzner_br_avg_trimmed$Attenuated) / nrow(muenzner_br_avg_trimmed)

dc_gene_list <- list(
  Muenzner = muenzner_dc_trimmed$gene_name[muenzner_dc_trimmed$Attenuated],
  BR = muenzner_br_avg$gene_name[muenzner_br_avg$Attenuated],
  BR_trimmed = muenzner_br_avg_trimmed$gene_name[muenzner_br_avg_trimmed$Attenuated]
)

length(dc_gene_list$Muenzner)
length(dc_gene_list$BR)
length(dc_gene_list$BR_trimmed)

jaccard_index(dc_gene_list)
ggvenn::ggvenn(dc_gene_list, c("BR", "Muenzner"))
ggvenn::ggvenn(dc_gene_list, c("BR_trimmed", "Muenzner"))
ggvenn::ggvenn(dc_gene_list, c("BR", "BR_trimmed", "Muenzner"))

## Check correlation between BR and turnover rates to confirm results of Muenzner et al.
muenzner_br_hl <- muenzner_br_trimmed %>%
  rename(systematic_name = Protein) %>%
  inner_join(y = muenzner_turnover, by = c("strain", "systematic_name"))

cor.test(muenzner_br_hl$Buffering.ChrArmLevel.Ratio, muenzner_br_hl$HL)
cor.test(muenzner_br_hl$Buffering.ChrArmLevel.Ratio, muenzner_br_hl$kdp_value)

### Sample median
muenzner_turnover_strain <- muenzner_turnover %>%
  summarize(MedianTurnover = median(kdp_value, na.rm = TRUE),
            MedianHalfLife = median(HL, na.rm = TRUE),
            .by = strain)

muenzner_br_strain <- muenzner_br_trimmed %>%
  summarize(MedianBR = median(Buffering.ChrArmLevel.Ratio, na.rm = TRUE),
            .by = strain)

muenzner_br_hl_strain <- inner_join(x = muenzner_br_strain, y = muenzner_turnover_strain, by = "strain")

cor.test(muenzner_br_hl_strain$MedianBR, muenzner_br_hl_strain$MedianHalfLife)
cor.test(muenzner_br_hl_strain$MedianBR, muenzner_br_hl_strain$MedianTurnover)