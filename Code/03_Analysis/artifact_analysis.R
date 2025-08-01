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
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))
# No CN-dependent BRs for P0211 and Chunduri

# Plot regression between Buffering Ratio, Protein Expression, etc. to uncover non-linearities and compression artifacts
# for Dosage Compensation of heavily amplified genes

## Plots may be too large if all values are used
downsample_n <- 50000

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
           Protein.Expression.Normalized.Average = mean(Protein.Expression.Normalized, na.rm = TRUE),
           Protein.Expression.CV = sd(2^Protein.Expression.Normalized, na.rm = TRUE) / mean(2^Protein.Expression.Normalized, na.rm = TRUE)
    ) %>%
    ungroup()
}

expr_buf_procan <- expr_buf_procan %>%
  add_fields()
expr_buf_depmap <- expr_buf_depmap %>%
  add_fields()
expr_buf_cptac <- expr_buf_cptac %>%
  add_fields()

## Define datasets to be processed
datasets <- list(
  list(dataset = expr_buf_procan, name = "ProCan"),
  list(dataset = expr_buf_procan %>% filter_cn_gain(remove_below = "50%"), name = "ProCan_CN-Gain"),
  list(dataset = expr_buf_procan %>% filter_cn_loss(remove_above = "50%"), name = "ProCan_CN-Loss"),
  list(dataset = expr_buf_depmap, name = "DepMap"),
  list(dataset = expr_buf_depmap %>% filter_cn_gain(remove_below = "50%"), name = "DepMap_CN-Gain"),
  list(dataset = expr_buf_depmap %>% filter_cn_loss(remove_above = "50%"), name = "DepMap_CN-Loss"),
  list(dataset = expr_buf_cptac, name = "CPTAC"),
  list(dataset = expr_buf_cptac %>% filter_cn_gain(remove_below = "50%"), name = "CPTAC_CN-Gain"),
  list(dataset = expr_buf_cptac %>% filter_cn_loss(remove_above = "50%"), name = "CPTAC_CN-Loss")
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
       name = "log2fc_cn", label_coords = c(0.5, 10)),
  ## Mean Log2FC vs. Mean Copy Number
  list(x = "Gene.CopyNumber.Average", y = "Log2FC.Average",
       name = "log2fc-mean_cn-mean", label_coords = c(1, 2)),
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
       name = "br-mean_log2fc-mean", label_coords = c(1, 2)),
  ## Buffering Ratio vs. Protein Neutral CV
  list(x = "Protein Neutral CV", y = "Buffering.GeneLevel.Ratio",
       name = "protein-cv-neutral_br", label_coords = c(1, 10)),
  ## Mean Buffering Ratio vs. Protein Neutral CV
  list(x = "Protein Neutral CV", y = "Buffering.GeneLevel.Ratio.Average",
       name = "protein-cv-neutral_br-mean", label_coords = c(1, 2)),
  ## Mean Log2FC vs. Protein Neutral CV
  list(x = "Protein Neutral CV", y = "Log2FC.Average",
       name = "protein-cv-neutral_log2fc-mean", label_coords = c(1, 2)),
  ## Buffering Ratio vs. Protein CV
  list(x = "Protein.Expression.CV", y = "Buffering.GeneLevel.Ratio",
       name = "protein-cv_br", label_coords = c(1, 10)),
  ## Mean Buffering Ratio vs. Protein CV
  list(x = "Protein.Expression.CV", y = "Buffering.GeneLevel.Ratio.Average",
       name = "protein-cv_br-mean", label_coords = c(1, 2)),
  ## Mean Log2FC vs. Protein CV
  list(x = "Protein.Expression.CV", y = "Log2FC.Average",
       name = "protein-cv_log2fc-mean", label_coords = c(1, 2))
)

### Generate plots
print("Analyzing...")
pb <- txtProgressBar(min = 0, max = length(datasets) * length(comparisons), style = 3)
plots <- lapply(datasets, \(dataset) {
  lapply(comparisons, \(comparison) {
    name <- paste(comparison$name, dataset$name, sep = "_")
    filename <- paste0("regression_", comparison$name, ".png")
    sub_dir <- dataset$name

    set.seed(42)

    plot <- dataset$dataset %>%
      rename(x = comparison$x, y = comparison$y) %>%
      slice_sample(n = downsample_n) %>%
      distinct(x, y) %>%
      scatter_plot_reg_corr(x_col = x, y_col = y, x_lab = comparison$x, y_lab = comparison$y)

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

# === Combine Plots for publishing ===
plot1 <- purrr::list_flatten(plots[sapply(plots, \(x) x$name == "br_base-expr_ProCan_CN-Gain")])
plot2 <- purrr::list_flatten(plots[sapply(plots, \(x) x$name == "br_base-expr_ProCan_CN-Loss")])
plot3 <- purrr::list_flatten(plots[sapply(plots, \(x) x$name == "br_cn-diff_DepMap")])
plot4 <- purrr::list_flatten(plots[sapply(plots, \(x) x$name == "br_cn-diff_ProCan_CN-Gain")])

plot_publish <- cowplot::plot_grid(plot1$plot, plot2$plot, plot3$plot, plot4$plot,
                                   nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

cairo_pdf(here(plots_dir, "artifact_publish.pdf"), width = 11)
plot_publish
dev.off()

# === Bonus ===
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
