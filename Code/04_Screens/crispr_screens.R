library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(openxlsx)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
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
crispr_screens <- read_parquet(here(output_data_dir, "crispr_screens.parquet"))
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))


# === Combine Datasets ===

df_crispr_buf <- crispr_screens %>%
  inner_join(y = expr_buf_procan %>% select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name),
             by = c("CellLine.CustomId", "Protein.Uniprot.Accession", "Gene.Symbol")) %>%
  select(CellLine.CustomId, CellLine.Name, Protein.Uniprot.Accession, Gene.Symbol,
         Gene.ChromosomeArm, ChromosomeArm.CNA, Buffering.GeneLevel.Ratio, CRISPR.EffectScore, CRISPR.DependencyScore)

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
  group_by(Gene.Symbol) %>%
  mutate(
    # ToDo: Replace with rstatix cor_test
    CRISPR.EffectScore.Average = mean(CRISPR.EffectScore, na.rm = TRUE),
    CRISPR.EffectScore.Corr = cor.test(Buffering.GeneLevel.Ratio, CRISPR.EffectScore, method = "spearman")$estimate[["rho"]],
    CRISPR.EffectScore.Corr.p = cor.test(Buffering.GeneLevel.Ratio, CRISPR.EffectScore, method = "spearman")$p.value,
    CRISPR.DependencyScore.Average = mean(CRISPR.DependencyScore, na.rm = TRUE),
    CRISPR.DependencyScore.Corr = cor.test(Buffering.GeneLevel.Ratio, CRISPR.DependencyScore, method = "spearman")$estimate[["rho"]],
    CRISPR.DependencyScore.Corr.p = cor.test(Buffering.GeneLevel.Ratio, CRISPR.DependencyScore, method = "spearman")$p.value,
    ChromosomeArm.GainLossRatio = mean(ChromosomeArm.CNA, na.rm = TRUE)
  ) %>%
  ungroup()

color_mapping_effect <- scale_color_viridis_c(option = "F", direction = 1, begin = 0.1, end = 0.8)
color_mapping_dep <- scale_color_viridis_c(option = "G", direction = -1, begin = 0.1, end = 0.8)

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
selected_genes <- c("KRAS", "EGFR")

for (gene in selected_genes) {
  df_gene_corr %>%
    filter(Gene.Symbol == gene) %>%
    scatter_plot_reg_corr(CRISPR.DependencyScore, Buffering.GeneLevel.Ratio,
                          point_size = 2, title_prefix = gene) %>%
    save_plot(paste0("dependency_buffering_", gene, "_scatterplot.png"))
}

# === Analyze Cell Line Sensitivity to Buffering ===
# ToDo: Calculate correlation of effect score with cell line ranking (per gene)
# ToDo: Discretize Cell Lines into Low-/Mid-/High-Buffering groups and check gene essentiality difference between groups
