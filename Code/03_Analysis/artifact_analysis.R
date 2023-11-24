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
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Artifact")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))

# Plot regression between Buffering Ratio, Protein Expression, etc. to uncover non-linearities and compression artifacts
# for Dosage Compensation of heavily amplified genes

## Buffering Ratio vs. Baseline Expression
expr_buf_procan %>%
  scatter_plot_regression(Protein.Expression.Baseline, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Protein.Expression.Baseline,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_br_base-expr.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Protein.Expression.Baseline, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Protein.Expression.Baseline,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_br_base-expr_cn-gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Protein.Expression.Baseline, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Protein.Expression.Baseline,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_br_base-expr_cn-loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Protein.Expression.Baseline = mean(Protein.Expression.Baseline, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Baseline, Buffering.GeneLevel.Ratio.Average,
                          Buffering.GeneLevel.Ratio.Average ~ Protein.Expression.Baseline,
                          label_coords = c(9, 1.5)) %>%
  save_plot("regression_br_base-expr_avg.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Buffering.GeneLevel.Ratio.SD = sd(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Protein.Expression.Baseline = mean(Protein.Expression.Baseline, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Baseline, Buffering.GeneLevel.Ratio.SD,
                          Buffering.GeneLevel.Ratio.SD ~ Protein.Expression.Baseline,
                          label_coords = c(9, 5)) %>%
  save_plot("regression_br_base-expr_sd.png")

## Log2FC vs. Baseline Expression
expr_buf_procan %>%
  scatter_plot_regression(Protein.Expression.Baseline, Log2FC,
                          Log2FC ~ Protein.Expression.Baseline,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_log2fc_base-expr.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Protein.Expression.Baseline, Log2FC,
                          Log2FC ~ Protein.Expression.Baseline,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_log2fc_base-expr_cn-gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Protein.Expression.Baseline, Log2FC,
                          Log2FC ~ Protein.Expression.Baseline,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_log2fc_base-expr_cn-loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Log2FC.Average = mean(Log2FC, na.rm = TRUE),
            Protein.Expression.Baseline = mean(Protein.Expression.Baseline, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Baseline, Log2FC.Average,
                          Log2FC.Average ~ Protein.Expression.Baseline,
                          label_coords = c(9, 1.5)) %>%
  save_plot("regression_log2fc_base-expr_avg.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Log2FC.SD = sd(Log2FC, na.rm = TRUE),
            Protein.Expression.Baseline = mean(Protein.Expression.Baseline, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Baseline, Log2FC.SD,
                          Log2FC.SD ~ Protein.Expression.Baseline,
                          label_coords = c(9, 5)) %>%
  save_plot("regression_log2fc_base-expr_sd.png")

## Buffering Ratio vs. Protein Expression
expr_buf_procan %>%
  scatter_plot_regression(Protein.Expression.Normalized, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Protein.Expression.Normalized,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_br_expr.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Protein.Expression.Normalized, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Protein.Expression.Normalized,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_br_expr_cn-gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Protein.Expression.Normalized, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Protein.Expression.Normalized,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_br_expr_cn-loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Protein.Expression.Normalized.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Normalized.Average, Buffering.GeneLevel.Ratio.Average,
                          Buffering.GeneLevel.Ratio.Average ~ Protein.Expression.Normalized.Average,
                          label_coords = c(9, 1.5)) %>%
  save_plot("regression_br_expr_avg.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Buffering.GeneLevel.Ratio.SD = sd(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Protein.Expression.Normalized.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Normalized.Average, Buffering.GeneLevel.Ratio.SD,
                          Buffering.GeneLevel.Ratio.SD ~ Protein.Expression.Normalized.Average,
                          label_coords = c(9, 5)) %>%
  save_plot("regression_br_expr_sd.png")

## Log2FC vs. Protein Expression
expr_buf_procan %>%
  scatter_plot_regression(Protein.Expression.Normalized, Log2FC,
                          Log2FC ~ Protein.Expression.Normalized,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_log2fc_expr.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Protein.Expression.Normalized, Log2FC,
                          Log2FC ~ Protein.Expression.Normalized,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_log2fc_expr_cn-gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Protein.Expression.Normalized, Log2FC,
                          Log2FC ~ Protein.Expression.Normalized,
                          label_coords = c(9, 10)) %>%
  save_plot("regression_log2fc_expr_cn-loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Log2FC.Average = mean(Log2FC, na.rm = TRUE),
            Protein.Expression.Normalized.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Normalized.Average, Log2FC.Average,
                          Log2FC.Average ~ Protein.Expression.Normalized.Average,
                          label_coords = c(9, 1.5)) %>%
  save_plot("regression_log2fc_expr_avg.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Log2FC.SD = sd(Log2FC, na.rm = TRUE),
            Protein.Expression.Normalized.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
  scatter_plot_regression(Protein.Expression.Normalized.Average, Log2FC.SD,
                          Log2FC.SD ~ Protein.Expression.Normalized.Average,
                          label_coords = c(9, 5)) %>%
  save_plot("regression_log2fc_expr_sd.png")

## Buffering Ratio vs Copy Number
expr_buf_procan %>%
  scatter_plot_regression(Gene.CopyNumber, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Gene.CopyNumber,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_br_cn.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Gene.CopyNumber, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Gene.CopyNumber,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_br_cn_gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Gene.CopyNumber, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Gene.CopyNumber,
                          label_coords = c(0.5, 10)) %>%
  save_plot("regression_br_cn_loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Gene.CopyNumber.Average = mean(Gene.CopyNumber, na.rm = TRUE)) %>%
  scatter_plot_regression(Gene.CopyNumber.Average, Buffering.GeneLevel.Ratio.Average,
                          Buffering.GeneLevel.Ratio.Average ~ Gene.CopyNumber.Average,
                          label_coords = c(1, 2)) %>%
  save_plot("regression_br_cn_avg.png")

## Log2FC vs Copy Number
expr_buf_procan %>%
  scatter_plot_regression(Gene.CopyNumber, Log2FC,
                          Log2FC ~ Gene.CopyNumber,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_log2fc_cn.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Gene.CopyNumber, Log2FC,
                          Log2FC ~ Gene.CopyNumber,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_log2fc_cn_gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Gene.CopyNumber, Log2FC,
                          Log2FC ~ Gene.CopyNumber,
                          label_coords = c(0.5, 10)) %>%
  save_plot("regression_log2fc_cn_loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Log2FC.Average = mean(Log2FC, na.rm = TRUE),
            Gene.CopyNumber.Average = mean(Gene.CopyNumber, na.rm = TRUE)) %>%
  scatter_plot_regression(Gene.CopyNumber.Average, Log2FC.Average,
                          Log2FC.Average ~ Gene.CopyNumber.Average,
                          label_coords = c(1, 2)) %>%
  save_plot("regression_log2fc_cn_avg.png")

## Buffering Ratio vs. Copy Number Diff
expr_buf_procan %>%
  mutate(Log2FC.CopyNumber = Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  scatter_plot_regression(Log2FC.CopyNumber, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Log2FC.CopyNumber,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_br_cn-diff.png")

expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  mutate(Log2FC.CopyNumber = Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  scatter_plot_regression(Log2FC.CopyNumber, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Log2FC.CopyNumber,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_br_cn-diff_gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  mutate(Log2FC.CopyNumber = Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  scatter_plot_regression(Log2FC.CopyNumber, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Log2FC.CopyNumber,
                          label_coords = c(-0.5, 10)) %>%
  save_plot("regression_br_cn-diff_loss.png")

expr_buf_procan %>%
  mutate(Log2FC.CopyNumber = Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  group_by(Gene.Symbol) %>%
  summarize(Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Log2FC.CopyNumber.Average = mean(Log2FC.CopyNumber, na.rm = TRUE)) %>%
  scatter_plot_regression(Log2FC.CopyNumber.Average, Buffering.GeneLevel.Ratio.Average,
                          Buffering.GeneLevel.Ratio.Average ~ Log2FC.CopyNumber.Average,
                          label_coords = c(0, 2)) %>%
  save_plot("regression_br_cn-diff_avg.png")


## Buffering Ratio vs. Log2FC
expr_buf_procan %>%
  scatter_plot_regression(Log2FC, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Log2FC,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_log2fc_br.png")


expr_buf_procan %>%
  filter_cn_gain(remove_below = "50%") %>%
  scatter_plot_regression(Log2FC, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Log2FC,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_log2fc_br_cn-gain.png")

expr_buf_procan %>%
  filter_cn_loss(remove_above = "50%") %>%
  scatter_plot_regression(Log2FC, Buffering.GeneLevel.Ratio,
                          Buffering.GeneLevel.Ratio ~ Log2FC,
                          label_coords = c(4, 10)) %>%
  save_plot("regression_log2fc_br_cn-loss.png")

expr_buf_procan %>%
  group_by(Gene.Symbol) %>%
  summarize(Log2FC.Average = mean(Log2FC, na.rm = TRUE),
            Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE)) %>%
  scatter_plot_regression(Log2FC.Average, Buffering.GeneLevel.Ratio.Average,
                          Buffering.GeneLevel.Ratio.Average ~ Log2FC.Average,
                          label_coords = c(1, 2)) %>%
  save_plot("regression_log2fc_br_avg.png")

## bonus
expr_buf_procan %>%
  select(Protein.Expression.Normalized, Buffering.GeneLevel.Ratio) %>%
  slice_sample(n = 800000) %>%
  drop_na() %>%
  ggplot() +
  aes(x = Protein.Expression.Normalized, y = Buffering.GeneLevel.Ratio) +
  geom_density_2d_filled() +
  geom_point(alpha = 0.05, size = 0.05, color = "white") +
  geom_density_2d_filled(alpha = 0.4) +
  xlab(NULL) +
  ylab(NULL) +
  theme_void() +
  theme(legend.position = "none")
