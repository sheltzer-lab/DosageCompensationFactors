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
dir.ceate(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
crispr_screens <- read_parquet(here(output_data_dir, "crispr_screens.parquet"))
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))


# === Combine Datasets ===

df_crispr_buf <- crispr_screens %>%
  inner_join(y = expr_buf_procan %>% select(-CellLine.DepMapModelId, -CellLine.SangerModelId, -CellLine.Name),
             by = c("CellLine.CustomId", "Protein.Uniprot.Accession", "Gene.Symbol")) %>%
  select(CellLine.CustomId, CellLine.Name, Protein.Uniprot.Accession, Gene.Symbol,
         Buffering.GeneLevel.Ratio, CRISPR.EffectScore)

## Check distributions
df_crispr_buf %>%
  ggplot() +
  aes(Buffering.GeneLevel.Ratio) +
  geom_density(na.rm = TRUE)

df_crispr_buf %>%
  ggplot() +
  aes(CRISPR.EffectScore) +
  geom_density(na.rm = TRUE)

# === Analyze Gene-wise Essentiality-Buffering-Correlation ===

df_gene_corr <- df_crispr_buf %>%
  group_by(Protein.Uniprot.Accession, Gene.Symbol) %>%
  mutate(
    # ToDo: Avoid calculating correlation twice
    CRISPR.EffectScore.Average = mean(CRISPR.EffectScore, na.rm = TRUE),
    Corr = cor.test(Buffering.GeneLevel.Ratio, CRISPR.EffectScore, method = "spearman")$estimate[["rho"]],
    Corr.p = cor.test(Buffering.GeneLevel.Ratio, CRISPR.EffectScore, method = "spearman")$p.value
  ) %>%
  arrange(Corr) %>%
  ungroup()

color_mapping <- scale_color_viridis_c(option = "D")

df_gene_corr %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  plot_volcano(Corr, Corr.p, Gene.Symbol, CRISPR.EffectScore.Average, color_mapping = color_mapping) %>%
  save_plot("ko-effect_buffering_correlation_volcano.png")

bot_corr <- df_gene_corr %>%
  filter(Corr.p < p_threshold) %>%
  filter(Corr < -0.2) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol)

top_corr <- df_gene_corr %>%
  filter(Corr.p < p_threshold) %>%
  filter(Corr > 0.2) %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  select(Gene.Symbol)

df_gene_corr %>%
  semi_join(y = bot_corr, by = "Gene.Symbol") %>%
  drop_na() %>%
  mutate(Label = paste0(Gene.Symbol, " (ρ = ", format(round(Corr, 3),
                                                      nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, desc(Corr))) %>%
  arrange(Buffering.GeneLevel.Ratio) %>%
  jittered_boxplot(Label, CRISPR.EffectScore, Buffering.GeneLevel.Ratio,
                   alpha = 1/2, jitter_width = 0.25) %>%
  save_plot("ko-effect_buffering_bot-correlation.png", height = 220, width = 180)

df_gene_corr %>%
  semi_join(y = top_corr, by = "Gene.Symbol") %>%
  drop_na() %>%
  mutate(Label = paste0(Gene.Symbol, " (ρ = ", format(round(Corr, 3),
                                                      nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, Corr)) %>%
  arrange(Buffering.GeneLevel.Ratio) %>%
  jittered_boxplot(Label, CRISPR.EffectScore, Buffering.GeneLevel.Ratio,
                   alpha = 1/2, jitter_width = 0.25) %>%
  save_plot("ko-effect_buffering_top-correlation.png", height = 220, width = 180)

# ToDo: Calculate correlation of effect score with cell line ranking (per gene)

# === Analyze Cell Line Sensitivity to Buffering ===
