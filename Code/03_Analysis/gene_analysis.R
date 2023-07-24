library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "03_Analysis", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Gene")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))

# === Determine Genes that are significantly Buffered on average ===

signif_buf_genes <- function(df, buffering_col, gene_col) {
  df %>%
  select({ { gene_col } }, { { buffering_col } }) %>%
  drop_na() %>%
  group_by({ { gene_col } }) %>%
  # Avoid having not enough samples for t-test
  add_count({ { gene_col } }) %>%
  filter(n > 1) %>%
  summarize(TTest.p = t.test({ { buffering_col } }, mu = 0)$p.value,
            Buffering.Ratio.Average = mean({ { buffering_col } })) %>%
  mutate(TTest.p.adjusted = p.adjust(TTest.p, method = "BY"))
}

# ToDo: Repeat analysis with Log2FC (Gene Level & Chromosome Level)

test_all <- expr_buf_goncalves %>%
  signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)

plot_all <- test_all %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
  plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
  save_plot("volcano_buffering_all.png")

test_gain <- expr_buf_goncalves %>%
  filter_cn_gain() %>%
  signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)

plot_gain <- test_gain %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
  plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
  save_plot("volcano_buffering_cn-gain.png")

test_loss <- expr_buf_goncalves %>%
  filter_cn_loss() %>%
  signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)

plot_loss <- test_loss %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
  plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
  save_plot("volcano_buffering_cn-loss.png")

test_filtered <- expr_buf_goncalves %>%
  filter_cn_diff_quantiles() %>%
  signif_buf_genes(Buffering.GeneLevel.Ratio, Gene.Symbol)

plot_filtered <- test_filtered %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio.Average)) %>%
  plot_volcano_buffered(Buffering.Ratio.Average, TTest.p.adjusted, Gene.Symbol, Buffering.Class) %>%
  save_plot("volcano_buffering_cn-filtered.png")