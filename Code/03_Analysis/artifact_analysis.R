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
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

# Plot regression between Buffering Ratio, Protein Expression, etc. to uncover non-linearities and compression artifacts
# for Dosage Compensation of heavily amplified genes

## Add additional fields (mean, sd)
add_fields <- function(df_buf) {
  df_buf %>%
    group_by(Gene.Symbol) %>%
    mutate(Gene.CopyNumber.Average = mean(Gene.CopyNumber, na.rm = TRUE),
           Log2FC.CopyNumber = Gene.CopyNumber - Gene.CopyNumber.Baseline,
           Log2FC.CopyNumber.Average = mean(Log2FC.CopyNumber, na.rm = TRUE),
           Log2FC.Average = mean(Log2FC, na.rm = TRUE),
           Log2FC.SD = sd(Log2FC, na.rm = TRUE),
           Buffering.GeneLevel.Ratio.Average = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
           Buffering.GeneLevel.Ratio.SD = sd(Buffering.GeneLevel.Ratio, na.rm = TRUE),
           Protein.Expression.Normalized.Average = mean(Protein.Expression.Normalized, na.rm = TRUE)) %>%
    ungroup()
}

expr_buf_procan <- expr_buf_procan %>%
  add_fields()
expr_buf_depmap <- expr_buf_depmap %>%
  add_fields()

## Define datasets to be processed
datasets <- list(
  list(dataset = expr_buf_procan, name = "ProCan"),
  list(dataset = expr_buf_procan %>% filter_cn_gain(remove_below = "50%"), name = "ProCan_CN-Gain"),
  list(dataset = expr_buf_procan %>% filter_cn_loss(remove_above = "50%"), name = "ProCan_CN-Loss"),
  list(dataset = expr_buf_depmap, name = "DepMap"),
  list(dataset = expr_buf_depmap %>% filter_cn_gain(remove_below = "50%"), name = "DepMap_CN-Gain"),
  list(dataset = expr_buf_depmap %>% filter_cn_loss(remove_above = "50%"), name = "DepMap_CN-Loss")
)

## Define conditions to compare
comparisons <- list(
  ## Buffering Ratio vs. Baseline Expression
  list(x = "Protein.Expression.Baseline", y = "Buffering.GeneLevel.Ratio",
       name = "br_base-expr", label_coords = c(9, 10)),
  ## Mean Buffering Ratio vs. Baseline Expression
  list(x = "Protein.Expression.Baseline", y = "Buffering.GeneLevel.Ratio.Average",
       name = "br-mean_base-expr", label_coords = c(9, 1.5)),
  ## Standard Deviation of Buffering Ratio vs. Baseline Expression
  list(x = "Protein.Expression.Baseline", y = "Buffering.GeneLevel.Ratio.SD",
       name = "br-sd_base-expr", label_coords = c(9, 5)),
  ## Log2FC vs. Baseline Expression
  list(x = "Protein.Expression.Baseline", y = "Log2FC",
       name = "log2fc_base-expr", label_coords = c(9, 10)),
  ## Mean Log2FC vs. Baseline Expression
  list(x = "Protein.Expression.Baseline", y = "Log2FC.Average",
       name = "log2fc-mean_base-expr", label_coords = c(9, 1.5)),
  ## Standard Deviation of Log2FC vs. Baseline Expression
  list(x = "Protein.Expression.Baseline", y = "Log2FC.SD",
       name = "log2fc-sd_base-expr", label_coords = c(9, 5)),
  ## Buffering Ratio vs. Protein Expression
  list(x = "Protein.Expression.Normalized", y = "Buffering.GeneLevel.Ratio",
       name = "br_expr", label_coords = c(9, 10)),
  ## Mean Buffering Ratio vs. Mean Protein Expression
  list(x = "Protein.Expression.Normalized.Average", y = "Buffering.GeneLevel.Ratio.Average",
       name = "br-mean_expr-mean", label_coords = c(9, 1.5)),
  ## Standard Deviation of Buffering Ratio vs. Mean Protein Expression
  list(x = "Protein.Expression.Normalized.Average", y = "Buffering.GeneLevel.Ratio.SD",
       name = "br-sd_expr-mean", label_coords = c(9, 5)),
  ## Log2FC vs. Protein Expression
  list(x = "Protein.Expression.Normalized", y = "Log2FC",
       name = "log2fc_expr", label_coords = c(9, 10)),
  ## Mean Log2FC vs. Mean Protein Expression
  list(x = "Protein.Expression.Normalized.Average", y = "Log2FC.Average",
       name = "log2fc-mean_expr-mean", label_coords = c(9, 1.5)),
  ## Standard Deviation of Log2FC vs. Mean Protein Expression
  list(x = "Protein.Expression.Normalized.Average", y = "Log2FC.SD",
       name = "log2fc-sd_expr-mean", label_coords = c(9, 5)),
  ## Buffering Ratio vs. Copy Number
  list(x = "Gene.CopyNumber", y = "Buffering.GeneLevel.Ratio",
       name = "br_cn", label_coords = c(4, 10)),
  ## Mean Buffering Ratio vs. Mean Copy Number
  list(x = "Gene.CopyNumber.Average", y = "Buffering.GeneLevel.Ratio.Average",
       name = "br-mean_cn-mean", label_coords = c(0.5, 10)),
  ## Log2FC vs. Copy Number
  list(x = "Gene.CopyNumber", y = "Log2FC",
       name = "br_cn", label_coords = c(0.5, 10)),
  ## Mean Log2FC vs. Mean Copy Number
  list(x = "Gene.CopyNumber.Average", y = "Log2FC.Average",
       name = "br-mean_cn-mean", label_coords = c(1, 2)),
  ## Buffering Ratio vs. Copy Number Diff
  list(x = "Log2FC.CopyNumber", y = "Buffering.GeneLevel.Ratio",
       name = "br_cn-diff", label_coords = c(4, 10)),
  ## Mean Buffering Ratio vs. Mean Copy Number Diff
  list(x = "Log2FC.CopyNumber.Average", y = "Buffering.GeneLevel.Ratio.Average",
       name = "br-mean_cn-diff-mean", label_coords = c(0, 2)),
  ## Buffering Ratio vs. Log2FC
  list(x = "Log2FC", y = "Buffering.GeneLevel.Ratio",
       name = "br_log2fc", label_coords = c(4, 10)),
  ## Mean Buffering Ratio vs. Mean Log2FC
  list(x = "Log2FC.Average", y = "Buffering.GeneLevel.Ratio.Average",
       name = "br-mean_log2fc-mean", label_coords = c(1, 2))
)

## Compare Variables
compare_variables <- function(df, var_x, var_y, cor_method = "spearman", x_lab = NULL, y_lab = NULL,
                              cor_symbol = utf8_rho, label_coords = NULL) {
  df_ <- df %>%
    select({ { var_x } }, { { var_y } }) %>%
    rename(x = { { var_x } }, y = { { var_y } }) %>%
    distinct()
  cor_result <- df_ %>%
    rstatix::cor_test(x, y, method = cor_method)
  plot_title <- paste0("Correlation (", cor_method , "): ",
                       print_corr(cor_result$cor, p.value = cor_result$p, signif = TRUE, estimate_symbol = cor_symbol))

  df_ %>%
    scatter_plot_regression(x, y, y ~ x, x_lab = x_lab, y_lab = y_lab,
                            label_coords = label_coords, title = plot_title)
}

### Generate plots
print("Analyzing...")
pb <- txtProgressBar(min = 0, max = length(datasets) * length(comparisons), style = 3)
plots <- lapply(datasets, \(dataset) {
  lapply(comparisons, \(comparison) {
    name <- paste(comparison$name, dataset$name, sep = "_")
    filename <- paste0("regression_", comparison$name, ".png")
    sub_dir <- dataset$name

    plot <- dataset$dataset %>%
      rename(x = comparison$x, y = comparison$y) %>%
      compare_variables(var_x = x, var_y = y, x_lab = comparison$x, y_lab = comparison$y)

    setTxtProgressBar(pb, pb$getVal() + 1)
    return(list(name = name, filename = filename, sub_dir = sub_dir, plot = plot))
  })
})
close(pb)
plots <- purrr::list_flatten(plots)

unique_sub_dirs <- unique(unlist(lapply(plots, function(x) x$sub_dir)))
for (sub_dir in unique_sub_dirs) {
  dir.create(here(plots_dir, sub_dir), recursive = TRUE)
}

### Save plots
print("Saving plots to disk...")
pb <- txtProgressBar(min = 0, max = length(plots), style = 3)
plots_output <- lapply(plots, \(plot) {
  suppressMessages({
    plot$plot %>%
      save_plot(plot$filename, dir = here(plots_dir, plot$sub_dir))
  })
  setTxtProgressBar(pb, pb$getVal() + 1)
  return(plot)
})
close(pb)

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
