library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
models_base_dir <- here("Output", "Models")
tables_dir <- here(tables_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

# === Univariate Ranks Panel ===
rank_gain_all <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain.parquet"))
rank_loss_all <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss.parquet"))
rank_loss_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss_WGD.parquet"))
rank_loss_no_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss_NoWGD.parquet"))
rank_gain_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain_WGD.parquet"))
rank_gain_no_wgd <- read_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain_NoWGD.parquet"))

rank_univariate <- bind_rows(rank_gain_all %>% mutate(Condition = "Gain"),
                             rank_loss_all %>% mutate(Condition = "Loss"))
rank_loss <- bind_rows(rank_loss_wgd %>% mutate(Condition = "WGD"),
                       rank_loss_no_wgd %>% mutate(Condition = "Non-WGD"))
rank_gain <- bind_rows(rank_gain_wgd %>% mutate(Condition = "WGD"),
                       rank_gain_no_wgd %>% mutate(Condition = "Non-WGD"))

rank_loss_heatmap <- rank_loss %>%
  # Shorten list of factors
  group_by(DosageCompensation.Factor) %>%
  mutate(MaxMNR = max(MeanNormRank)) %>%
  ungroup() %>%
  filter(MaxMNR > 0.5) %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  scale_fill_viridis_c(option = "magma", direction = 1, end = 0.9, oob = scales::squish, limits = c(0.5, 1),
                       breaks = seq(0.5, 1, 0.1), labels = c("\u22640.5", "0.6", "0.7", "0.8", "0.9", "1.0")) +
  coord_flip() +
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
  coord_flip() +
  ggtitle("Gain") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

panel_rank <- cowplot::plot_grid(rank_gain_heatmap, rank_loss_heatmap,
                                 labels = c("A", "B"), nrow = 1, ncol = 2, rel_widths = c(0.8, 1))

# === Univariate Bootstrap Panel ===
bootstrap_results <- read_parquet(here(output_data_dir, 'bootstrap_univariate.parquet'))

bootstrap_results_filtered <- bootstrap_results %>%
  filter(Dataset == "ProCan") %>%
  group_by(DosageCompensation.Factor) %>%
  mutate(MedianAUC = median(DosageCompensation.Factor.ROC.AUC)) %>%
  ungroup() %>%
  filter(MedianAUC > 0.51) %>%
  select(-MedianAUC)

excluded_factors <- c("Protein Abundance", "mRNA Abundance",
                      "Transcription Factors (Repressor)", "Transcription Factors (Activator)")

bootstrap_chr_gain <- bootstrap_results_filtered %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Chromosome Arm" & Event == "Gain")
bootstrap_chr_loss <- bootstrap_results_filtered %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Chromosome Arm" & Event == "Loss")
bootstrap_cn_gain <- bootstrap_results_filtered %>%
  filter(!(DosageCompensation.Factor %in% excluded_factors) & Level == "Gene Copy Number" & Event == "Gain")
bootstrap_cn_loss <- bootstrap_results_filtered %>%
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

# === Tables ===
br_factor_cor <- read_parquet(here(output_data_dir, 'factor_correlation_univariate.parquet'))

bootstrap_summary <- bootstrap_results %>%
  summarize(ROC.AUC.Median = median(DosageCompensation.Factor.ROC.AUC, na.rm = TRUE),
            ROC.AUC.Min = min(DosageCompensation.Factor.ROC.AUC, na.rm = TRUE),
            ROC.AUC.Max = max(DosageCompensation.Factor.ROC.AUC, na.rm = TRUE),
            Observations.Median = median(DosageCompensation.Factor.Observations, na.rm = TRUE),
            Observations.Min = min(DosageCompensation.Factor.Observations, na.rm = TRUE),
            Observations.Max = max(DosageCompensation.Factor.Observations, na.rm = TRUE),
            Bootstrap.Samples = n(),
            .by = c("DosageCompensation.Factor", "Dataset", "Level", "Event"))

wb <- createWorkbook()
sheet_rank_gain <- addWorksheet(wb, "Gain")
sheet_rank_loss <- addWorksheet(wb, "Loss")
sheet_rank_gain_wgd <- addWorksheet(wb, "Gain, WGD")
sheet_rank_loss_wgd <- addWorksheet(wb, "Loss, WGD")
sheet_rank_gain_nonwgd <- addWorksheet(wb, "Gain, Non-WGD")
sheet_rank_loss_nonwgd <- addWorksheet(wb, "Loss, Non-WGD")
sheet_bootstrap_summary <- addWorksheet(wb, "Bootstrapped ROC AUC (Summary)")
sheet_corr <- addWorksheet(wb, "BR-Factor Correlation")
writeDataTable(wb = wb, sheet = sheet_rank_gain, x = rank_gain_all)
writeDataTable(wb = wb, sheet = sheet_rank_loss, x = rank_loss_all)
writeDataTable(wb = wb, sheet = sheet_rank_gain_wgd, x = rank_gain_wgd)
writeDataTable(wb = wb, sheet = sheet_rank_loss_wgd, x = rank_loss_wgd)
writeDataTable(wb = wb, sheet = sheet_rank_gain_nonwgd, x = rank_gain_no_wgd)
writeDataTable(wb = wb, sheet = sheet_rank_loss_nonwgd, x = rank_loss_no_wgd)
writeDataTable(wb = wb, sheet = sheet_bootstrap_summary, x = bootstrap_summary)
writeDataTable(wb = wb, sheet = sheet_corr, x = br_factor_cor)
saveWorkbook(wb, here(tables_dir, "supplementary_table4.xlsx"), overwrite = TRUE)

sup_table5 <- write_parquet(bootstrap_auc, here(tables_dir, "supplementary_table5.parquet"), version = "2.4")

# TODO: Add field descriptions

# === Supplemental Figures ===
## ROC AUC ranks (all)
rank_heatmap <- rank_univariate %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  coord_flip() +
  labs(fill = "MNR\n(ROC AUC)") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

## BR-Factor correlation
panel_br_factor_cor <- br_factor_cor %>%
  mutate(Label = map_signif(p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, cor, .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Dataset, Gene.CNV, cor, Label) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = "", fill = cor, label = Label) +
  geom_raster() +
  geom_text(color = "black") +
  ggh4x::facet_nested(Dataset + Gene.CNV ~ ., switch = "y") +
  scale_fill_gradientn(colors = c(bidirectional_color_pal_viridis[2], "white", bidirectional_color_pal_viridis[4]),
                         space = "Lab", limits = c(-0.25, 0.25), oob = scales::squish) +
  labs(x = "Dosage Compensation Factor", y = "", fill = "Buffering Ratio Correlation") +
  cowplot::theme_minimal_grid() +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color = "lightgrey"),
        axis.ticks.y = element_blank(),
        legend.key.width = unit(24, "points"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 1, vjust = 0.9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 15, 5, 15), "mm"))

## ROC AUC Correlation Plots
auc_corr_plots <- list()

auc_corr_plots$ChrGain$data <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Chromosome Arm" & Event == "Gain")
auc_corr_plots$ChrLoss$data <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Chromosome Arm" & Event == "Loss")
auc_corr_plots$CNGain$data <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Gene Copy Number" & Event == "Gain")
auc_corr_plots$CNLoss$data <- filter(bootstrap_results, Dataset == "ProCan" & Level == "Gene Copy Number" & Event == "Loss")

for (condition in names(auc_corr_plots)) {
  auc_corr_plots[[condition]]$plot <- auc_corr_plots[[condition]]$data %>%
    select(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC, Condition, Bootstrap.Sample) %>%
    pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    select(-Condition, -Bootstrap.Sample) %>%
    plot_correlation_matrix() +
    theme(legend.position = "none",
          text = element_text(size = 10))
}

## Factor Correlation
panel_factor_cor <- read_parquet(here(output_data_dir, 'dosage_compensation_factors_filtered.parquet')) %>%
  select(all_of(dc_factor_cols)) %>%
  plot_correlation_matrix() +
  theme(legend.key.width = unit(24, "points"),
        text = element_text(size = 11),
        legend.title = element_text(size = base_size),
        legend.text = element_text(size = 12))

## Combine figures
panel_roc_cor <- cowplot::plot_grid(
  #auc_corr_plots$ChrGain$plot + ggtitle("Chromosome Arm Gain"), auc_corr_plots$ChrLoss$plot  + ggtitle("Chromosome Arm Loss"),
  auc_corr_plots$CNGain$plot  + ggtitle("Bootstrapped ROC AUCs, Gene CN Gain"),
  auc_corr_plots$CNLoss$plot  + ggtitle("Bootstrapped ROC AUCs, Gene CN Loss"),
  labels = c("C", "D")
)

figure_s3_sub1 <- cowplot::plot_grid(rank_heatmap, panel_factor_cor,
                                     rel_widths = c(0.85, 1), labels = c("A", "B"))
figure_s3 <- cowplot::plot_grid(figure_s3_sub1, panel_roc_cor, panel_br_factor_cor,
                                nrow = 3, rel_heights = c(1, 0.9, 0.8), labels = c("", "", "E"))

cairo_pdf(here(plots_dir, "figure_s3.pdf"), width = 13, height = 19)
figure_s3
dev.off()
