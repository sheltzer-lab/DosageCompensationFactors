library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(openxlsx)
library(forcats)
library(magrittr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Screens", "CRISPR")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))
crispr_screens <- read_parquet(here(output_data_dir, "crispr_screens.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
metadata_depmap <- read_csv_arrow(here(external_data_dir, "CopyNumber", "DepMap", "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID")
tp53_mut <- read_parquet(here(output_data_dir, "damaging_mutations_depmap.parquet")) %>%
  filter(Gene.Symbol == "TP53") %>%
  select(-Gene.Symbol) %>%
  rename(TP53.Mutated = MutationStatus)

# === Combine Datasets ===
df_crispr_buf <- crispr_screens %>%
  inner_join(y = expr_buf_depmap %>% select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name),
             by = c("Model.ID", "Protein.Uniprot.Accession", "Gene.Symbol")) %>%
  select(Model.ID, CellLine.Name, Protein.Uniprot.Accession, Gene.Symbol, Gene.CopyNumber, Gene.CopyNumber.Baseline,
         Buffering.GeneLevel.Ratio, CRISPR.EffectScore, CRISPR.DependencyScore)

## Check distributions
df_crispr_buf %>%
  ggplot() +
  aes(Buffering.GeneLevel.Ratio) +
  geom_density(na.rm = TRUE)

df_crispr_buf %>%
  ggplot() +
  aes(CRISPR.EffectScore) +
  geom_density(na.rm = TRUE)

df_crispr_buf %>%
  ggplot() +
  aes(CRISPR.DependencyScore) +
  geom_density(na.rm = TRUE)

# === Analyze Gene-wise Essentiality-Buffering-Correlation ===
df_gene_corr <- df_crispr_buf %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  group_by(Gene.Symbol) %>%
  mutate(
    CRISPR.EffectScore.Average = mean(CRISPR.EffectScore, na.rm = TRUE),
    CRISPR.EffectScore.Corr = cor.test(Buffering.GeneLevel.Ratio, CRISPR.EffectScore, method = "spearman")$estimate[["rho"]],
    CRISPR.EffectScore.Corr.p = cor.test(Buffering.GeneLevel.Ratio, CRISPR.EffectScore, method = "spearman")$p.value,
    CRISPR.DependencyScore.Average = mean(CRISPR.DependencyScore, na.rm = TRUE),
    CRISPR.DependencyScore.Corr = cor.test(Buffering.GeneLevel.Ratio, CRISPR.DependencyScore, method = "spearman")$estimate[["rho"]],
    CRISPR.DependencyScore.Corr.p = cor.test(Buffering.GeneLevel.Ratio, CRISPR.DependencyScore, method = "spearman")$p.value,
    Gene.CopyNumber.GainLossRatio = mean(Gene.CopyNumber - Gene.CopyNumber.Baseline, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  write_parquet(here(output_data_dir, "buffering_gene_dependency_depmap.parquet")) %T>%
  write.xlsx(here(tables_base_dir, "buffering_gene_dependency_depmap.xlsx"), colNames = TRUE)

color_mapping_effect <- scale_color_viridis_c(option = "F", direction = 1, begin = 0.1, end = 0.8)
color_mapping_dep <- scale_color_viridis_c(option = "G", direction = -1, begin = 0.1, end = 0.8)
color_mapping_driver <- scale_color_viridis_d(option = "plasma", na.value = color_palettes$Missing, end = 0.9)

df_gene_corr %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  arrange(desc(CRISPR.EffectScore.Average)) %>%
  mutate(Label = if_else(abs(CRISPR.EffectScore.Corr) > 0.2 &
                           CRISPR.EffectScore.Corr.p < p_threshold &
                           CRISPR.EffectScore.Average < -0.1,
                         Gene.Symbol, NA)) %>%
  plot_volcano(CRISPR.EffectScore.Corr, CRISPR.EffectScore.Corr.p, Label, CRISPR.EffectScore.Average,
               color_mapping = color_mapping_effect, value_threshold = 0.2) %>%
  save_plot("ko-effect_buffering_correlation_volcano.png", width = 300, height = 250)

df_gene_corr %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  arrange(CRISPR.DependencyScore.Average) %>%
  mutate(Label = if_else(abs(CRISPR.DependencyScore.Corr) > 0.2 &
                           CRISPR.DependencyScore.Corr.p < p_threshold &
                           CRISPR.DependencyScore.Average > 0.1,
                         Gene.Symbol, NA)) %>%
  plot_volcano(CRISPR.DependencyScore.Corr, CRISPR.DependencyScore.Corr.p, Label, CRISPR.DependencyScore.Average,
               color_mapping = color_mapping_dep, value_threshold = 0.2) %>%
  save_plot("dependency_buffering_correlation_volcano.png", width = 300, height = 250)

df_gene_corr %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(abs(CRISPR.DependencyScore.Corr) > 0.2 &
                           CRISPR.DependencyScore.Corr.p < p_threshold &
                           !is.na(CancerDriverMode),
                         Gene.Symbol, NA)) %>%
  arrange(!is.na(CancerDriverMode)) %>%
  plot_volcano(CRISPR.DependencyScore.Corr, CRISPR.DependencyScore.Corr.p, Label, CancerDriverMode,
               color_mapping = color_mapping_driver, value_threshold = 0.2) %>%
  save_plot("dependency_buffering_correlation_volcano_cancer-genes.png", width = 300, height = 250)

bot_corr <- df_gene_corr %>%
  filter(CRISPR.DependencyScore.Corr.p < p_threshold) %>%
  filter(CRISPR.DependencyScore.Corr < -0.2) %>%
  filter(CRISPR.DependencyScore.Average > 0.02) %>%
  distinct(Gene.Symbol, .keep_all = FALSE)

top_corr <- df_gene_corr %>%
  filter(CRISPR.DependencyScore.Corr.p < p_threshold) %>%
  filter(CRISPR.DependencyScore.Corr > 0.2) %>%
  filter(CRISPR.DependencyScore.Average > 0.02) %>%
  distinct(Gene.Symbol, .keep_all = FALSE)

df_gene_corr %>%
  semi_join(y = bot_corr, by = "Gene.Symbol") %>%
  drop_na() %>%
  mutate(Label = paste0(Gene.Symbol, " (", print_corr(CRISPR.DependencyScore.Corr), ")")) %>%
  mutate(Label = fct_reorder(Label, desc(CRISPR.DependencyScore.Corr))) %>%
  arrange(Buffering.GeneLevel.Ratio) %>%
  jittered_boxplot(Label, CRISPR.DependencyScore, Buffering.GeneLevel.Ratio,
                   alpha = 3/4, jitter_width = 0.2) %>%
  save_plot("dependency_buffering_bot-correlation.png", height = 220, width = 180)

df_gene_corr %>%
  semi_join(y = top_corr, by = "Gene.Symbol") %>%
  drop_na() %>%
  mutate(Label = paste0(Gene.Symbol, " (", print_corr(CRISPR.DependencyScore.Corr), ")")) %>%
  mutate(Label = fct_reorder(Label, CRISPR.DependencyScore.Corr)) %>%
  arrange(Buffering.GeneLevel.Ratio) %>%
  jittered_boxplot(Label, CRISPR.DependencyScore, Buffering.GeneLevel.Ratio,
                   alpha = 3/4, jitter_width = 0.2) %>%
  save_plot("dependency_buffering_top-correlation.png", height = 220, width = 180)

## Look at selected genes in more detail
selected_genes <- c("KRAS", "EGFR", "TP53", "CDK6", "ELMO2", "CRKL")

for (gene in selected_genes) {
  df_gene_corr %>%
    filter(Gene.Symbol == gene) %>%
    scatter_plot_reg_corr(CRISPR.DependencyScore, Buffering.GeneLevel.Ratio,
                          point_size = 2, title_prefix = gene) %>%
    save_plot(paste0("dependency_buffering_", gene, "_scatterplot.png"))
}

# === Analyze Cell Line Sensitivity to Buffering ===
## High Buffering vs. Low Buffering
df_crispr_model_buf <- model_buf_agg %>%
  split_by_quantiles(Model.Buffering.MeanNormRank, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = crispr_screens, by = "Model.ID") %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, CRISPR.DependencyScore,
                          groups = c("Low", "High"), test = wilcox.test, centr = mean, log2fc_thresh = 0.10) %>%
  write_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap.parquet")) %T>%
  write.xlsx(here(tables_base_dir, "model_buffering_gene_dependency_depmap.xlsx"), colNames = TRUE)

df_crispr_model_buf %>%
  mutate(Label = if_else(!is.na(Significant), Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, value_threshold = 0.1) %>%
  save_plot("model_buffering_dependency_volcano.png", width = 300, height = 250)

df_crispr_model_buf %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA)) %>%
  arrange(!is.na(CancerDriverMode)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, CancerDriverMode,
               value_threshold = 0.1, color_mapping = color_mapping_driver) %>%
  save_plot("model_buffering_dependency_volcano_cancer-genes.png", width = 300, height = 250)

### Apply over representation analysis
ora_up <- df_crispr_model_buf %>%
  filter(Significant == "Up") %>%
  arrange(desc(Log2FC)) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_up %>%
  plot_terms() %>%
  save_plot("ora_terms_up.png", height = 300)

ora_down <- df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  arrange(Log2FC) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_down %>%
  plot_terms() %>%
  save_plot("ora_terms_down.png", height = 300)

### Check if model-wise and gene-wise results match
df_crispr_model_buf %>%
  filter(Significant == "Up") %>%
  pull(Gene.Symbol) %>%
  intersect(top_corr$Gene.Symbol)

df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  pull(Gene.Symbol) %>%
  intersect(bot_corr$Gene.Symbol)

### Get Genes that are essential in Buffering cells only
essential_buf_only <- df_crispr_model_buf %>%
  filter(Significant == "Up") %>%
  filter(Mean_GroupA < 0.5 & Mean_GroupB > 0.5)

essential_buf_only %>%
  arrange(desc(Log2FC)) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = TRUE) %>%
  plot_terms() %>%
  save_plot("ora_terms_essential_buffering_only.png", height = 300)

### Genes that are essential in Scaling cells only
essential_scaling_only <- df_crispr_model_buf %>%
  filter(Significant == "Down") %>%
  filter(Mean_GroupA > 0.5 & Mean_GroupB < 0.5)

essential_scaling_only %>%
  arrange(Log2FC) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = TRUE) %>%
  plot_terms() %>%
  save_plot("ora_terms_essential_scaling_only.png", height = 300)

# ToDo: Calculate correlation of effect score with cell line ranking (per gene)

# === Control: Control for suspension/adhesion growth method and TP53 status
df_crispr_model_buf_adherent <- model_buf_agg %>%
  split_by_quantiles(Model.Buffering.MeanNormRank, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = crispr_screens, by = "Model.ID") %>%
  inner_join(y = metadata_depmap, by = "CellLine.DepMapModelId") %>%
  filter(GrowthPattern == "Adherent") %>%
  select(Model.ID, Gene.Symbol, Model.Buffering.Group, CRISPR.DependencyScore) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, CRISPR.DependencyScore,
                          groups = c("Low", "High"), test = wilcox.test, centr = mean, log2fc_thresh = 0.10) %>%
  write_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap_adherent.parquet")) %T>%
  write.xlsx(here(tables_base_dir, "model_buffering_gene_dependency_depmap_adherent.xlsx"), colNames = TRUE)

df_crispr_model_buf_control <- model_buf_agg %>%
  split_by_quantiles(Model.Buffering.MeanNormRank, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = crispr_screens, by = "Model.ID") %>%
  inner_join(y = metadata_depmap, by = "CellLine.DepMapModelId") %>%
  inner_join(y = tp53_mut, by = "Model.ID") %>%
  filter(GrowthPattern == "Adherent" & !TP53.Mutated) %>%
  select(Model.ID, Gene.Symbol, Model.Buffering.Group, CRISPR.DependencyScore) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, CRISPR.DependencyScore,
                          groups = c("Low", "High"), test = wilcox.test, centr = mean, log2fc_thresh = 0.10) %>%
  write_parquet(here(output_data_dir, "model_buffering_gene_dependency_depmap_control.parquet")) %T>%
  write.xlsx(here(tables_base_dir, "model_buffering_gene_dependency_depmap_control.xlsx"), colNames = TRUE)

essential_scaling_only_adherent <- df_crispr_model_buf_adherent %>%
  semi_join(y = essential_scaling_only, by = c("Gene.Symbol", "Significant"))
