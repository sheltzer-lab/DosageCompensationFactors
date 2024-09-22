library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
models_base_dir <- here("Output", "Models")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

# === Univariate Ranks Panel ===
rank_gain <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain.parquet"))
rank_loss <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss.parquet"))
rank_loss_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss_WGD.parquet"))
rank_loss_no_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss_NoWGD.parquet"))
rank_gain_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain_WGD.parquet"))
rank_gain_no_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain_NoWGD.parquet"))

rank_univariate <- bind_rows(rank_gain %>% mutate(Condition = "Gain"),
                             rank_loss %>% mutate(Condition = "Loss"))
rank_loss <- bind_rows(rank_loss_wgd %>% mutate(Condition = "Loss, WGD"),
                       rank_loss_no_wgd %>% mutate(Condition = "Loss, Non-WGD"))
rank_gain <- bind_rows(rank_gain_wgd %>% mutate(Condition = "Gain, WGD"),
                       rank_gain_no_wgd %>% mutate(Condition = "Gain, Non-WGD"))

rank_heatmap <- rank_univariate %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  theme(legend.position = "right", legend.direction = "vertical")
rank_loss_heatmap <- rank_loss %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  theme(legend.position = "right", legend.direction = "vertical")
rank_gain_heatmap <- rank_gain %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  theme(legend.position = "right", legend.direction = "vertical")

panel_rank <- cowplot::plot_grid(rank_heatmap + coord_flip() + theme(legend.position = "none"),
                                 rank_gain_heatmap + coord_flip() + theme(legend.position = "none"),
                                 rank_loss_heatmap + coord_flip(),
                                 labels = c("A", "B", "C"),
                                 nrow = 1, ncol = 3, rel_widths = c(0.8, 0.8, 1))

# === Univariate Bootstrap Panel ===
bootstrap_results <- read_parquet(here(output_data_dir, 'bootstrap_univariate.parquet'))

bootstrap_chr_gain <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Chromosome Arm" & Event == "Gain")
bootstrap_chr_loss <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Chromosome Arm" & Event == "Loss")
bootstrap_cn_gain <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Gene Copy Number" & Event == "Gain")
bootstrap_cn_loss <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Gene Copy Number" & Event == "Loss")

results_chrgain_chrloss <- compare_conditions(bootstrap_chr_gain, bootstrap_chr_loss)
results_cngain_cnloss <- compare_conditions(bootstrap_cn_gain, bootstrap_cn_loss)
results_chrgain_cngain <- compare_conditions(bootstrap_chr_gain, bootstrap_cn_gain)
results_chrloss_cnloss <- compare_conditions(bootstrap_chr_loss, bootstrap_cn_loss)

bootstrap_auc <- bind_rows(bootstrap_chr_gain, bootstrap_chr_loss, bootstrap_cn_gain, bootstrap_cn_loss)

rank_tests <- list(
  comparisons = list(c("Chromosome Arm Gain", "Chromosome Arm Loss"),
                     c("Gene Copy Number Gain", "Gene Copy Number Loss"),
                     c("Chromosome Arm Gain", "Gene Copy Number Gain"),
                     c("Chromosome Arm Loss", "Gene Copy Number Loss")),
  annotations = c(print_corr_obj(results_chrgain_chrloss$rank_test, estimate_symbol = utf8_tau, map_p = TRUE),
                 print_corr_obj(results_cngain_cnloss$rank_test, estimate_symbol = utf8_tau, map_p = TRUE),
                 print_corr_obj(results_chrgain_cngain$rank_test, estimate_symbol = utf8_tau, map_p = TRUE),
                 print_corr_obj(results_chrloss_cnloss$rank_test, estimate_symbol = utf8_tau, map_p = TRUE)),
  y_position = c(1, 1, 1.5, 2)
)

panel_bootstrap <- bootstrap_auc %>%
  roc_auc_heatmap(rank_tests)

# === Combine Panels into Figure ===
figure3 <- cowplot::plot_grid(panel_rank, panel_bootstrap,
                              labels = c("", "D"),
                              nrow = 2, ncol = 1, rel_heights = c(1, 0.5))

cairo_pdf(here(plots_dir, "figure03.pdf"), width = 16, height = 15)
figure3
dev.off()