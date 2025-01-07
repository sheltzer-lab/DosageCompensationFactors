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
rank_loss <- bind_rows(rank_loss_wgd %>% mutate(Condition = "WGD"),
                       rank_loss_no_wgd %>% mutate(Condition = "Non-WGD"))
rank_gain <- bind_rows(rank_gain_wgd %>% mutate(Condition = "WGD"),
                       rank_gain_no_wgd %>% mutate(Condition = "Non-WGD"))

# To supplement
rank_heatmap <- rank_univariate %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  theme(legend.position = "right", legend.direction = "vertical")
# To main figure
rank_loss_heatmap <- rank_loss %>%
  # Shorten list of factors
  group_by(DosageCompensation.Factor) %>%
  mutate(MaxMNR = max(MeanNormRank)) %>%
  ungroup() %>%
  filter(MaxMNR > 0.5) %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  ggtitle("Loss") +
  labs(fill = "MNR(ROC AUC)") +
  theme(legend.position = "right", legend.direction = "vertical", axis.text.x = element_text(angle = 0, hjust = 0.5))
rank_gain_heatmap <- rank_gain %>%
  # Shorten list of factors
  group_by(DosageCompensation.Factor) %>%
  mutate(MaxMNR = max(MeanNormRank)) %>%
  ungroup() %>%
  filter(MaxMNR > 0.5) %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  ggtitle("Gain") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

panel_rank <- cowplot::plot_grid(rank_gain_heatmap + coord_flip(),
                                 rank_loss_heatmap + coord_flip(),
                                 labels = c("A", "B"),
                                 nrow = 1, ncol = 2, rel_widths = c(0.8, 1))

# === Univariate Bootstrap Panel ===
bootstrap_results <- read_parquet(here(output_data_dir, 'bootstrap_univariate.parquet')) %>%
  filter(Dataset == "ProCan") %>%
  group_by(DosageCompensation.Factor) %>%
  mutate(MedianAUC = median(DosageCompensation.Factor.ROC.AUC)) %>%
  ungroup() %>%
  filter(MedianAUC > 0.51) %>%
  select(-MedianAUC)

excluded_factors <- c("Protein Abundance", "mRNA Abundance",
                      "Transcription Factors (Repressor)", "Transcription Factors (Activator)")

bootstrap_chr_gain <- bootstrap_results %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Chromosome Arm" & Event == "Gain")
bootstrap_chr_loss <- bootstrap_results %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Chromosome Arm" & Event == "Loss")
bootstrap_cn_gain <- bootstrap_results %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Gene Copy Number" & Event == "Gain")
bootstrap_cn_loss <- bootstrap_results %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Gene Copy Number" & Event == "Loss")

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
  roc_auc_heatmap(rank_tests, color_lab = "Median ROC AUC")

# === Combine Panels into Figure ===
figure3 <- cowplot::plot_grid(panel_rank, panel_bootstrap,
                              labels = c("", "C"),
                              nrow = 2, ncol = 1, rel_heights = c(1, 0.8))

cairo_pdf(here(plots_dir, "figure03.pdf"), width = 11, height = 11)
figure3
dev.off()
