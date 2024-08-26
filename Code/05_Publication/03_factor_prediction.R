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

rank_univariate <- bind_rows(rank_gain %>% mutate(Condition = "Gain"),
                             rank_loss %>% mutate(Condition = "Loss"))
rank_wgd <- bind_rows(rank_loss_wgd %>% mutate(Condition = "Loss, WGD"),
                      rank_loss_no_wgd %>% mutate(Condition = "Loss, Non-WGD"))

rank_heatmap <- rank_univariate %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  theme(legend.position = "right", legend.direction = "vertical")
rank_wgd_heatmap <- rank_wgd %>%
  unidirectional_heatmap(DosageCompensation.Factor, Condition, MeanNormRank) +
  theme(legend.position = "right", legend.direction = "vertical")

panel_rank <- cowplot::plot_grid(rank_heatmap + coord_flip() + theme(legend.position = "none"),
                                 rank_wgd_heatmap + coord_flip(),
                                 labels = c("A", "B"),
                                 nrow = 1, ncol = 2, rel_widths = c(0.8, 1))

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

# === Multivariate Model Performance Panel ===
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

get_rocs <- function(df_models) {
  df_models %>%
    group_by(Model.Variant) %>%
    group_map(\(entry, group) {
      result <- list()
      result[[group$Model.Variant]] <- readRDS(here(models_base_dir, entry$Model.Filename))$evaluation$roc
      return(result)
    }) %>%
    purrr::list_flatten() %>%
    rocs_to_df() %>%
    mutate(Model.Variant = Name) %>%
    left_join(y = df_models, by = "Model.Variant",
              unmatched = "error", relationship = "many-to-one")
}

plot_rocs_publish <- function(df_rocs) {
  df_label <- df_rocs %>%
    distinct(Model.Variant, Model.Condition, AUC) %>%
    arrange(desc(Model.Variant)) %>%
    mutate(Model.Level = case_when(grepl("Chr", Model.Variant) & grepl("Average", Model.Variant) ~ "ChrArm (avg.)",
                                   grepl("Chr", Model.Variant) & !grepl("Average", Model.Variant) ~ "ChrArm",
                                   grepl("Gene", Model.Variant) & !grepl("Average", Model.Variant) ~ "Gene CN"),
           x = if_else(Model.Condition == "Gain", 0.25, 0),
           y = case_when(Model.Level == "ChrArm" ~ 0,
                         Model.Level == "ChrArm (avg.)" ~ 0.1,
                         Model.Level == "Gene CN" ~ 0.2),
           AUC = paste0("AUC = ", format(round(AUC, 3), nsmall = 3)))

  df_legend <- data.frame(
    Label = c("ChrArm", "ChrArm (avg.)", "Gene CN", "Gain", "Loss"),
    x = c(0.5, 0.5, 0.5, 0.25, 0),
    y = c(0.01, 0.11, 0.21, 0.28, 0.28)
  )

  df_rocs %>%
    arrange(Sensitivity) %>%
    ggplot() +
    aes(x = Specificity, y = Sensitivity, color = Model.Variant) +
    geom_abline(slope = 1, intercept = 1, color = "grey") +
    geom_line() +
    geom_label(data = df_label,
               mapping = aes(color = Model.Variant, label = AUC, x = x, y = y),
               hjust = 1, vjust = 0) +
    geom_text(data = df_legend,
              mapping = aes(label = Label, x = x, y = y, color = NULL),
              hjust = 1, vjust = 0, fontface = "bold") +
    scale_x_reverse(limits = c(1, 0)) +
    labs(x = "Specificity", y = "Sensitivity", color = "Model") +
    theme(legend.position = "none")
}

models_depmap <- shap_results %>%
  select(starts_with("Model")) %>%
  distinct() %>%
  filter(Model.Dataset == "DepMap" & Model.Subset == "All" & !grepl("Log2FC", Model.Variant)) %>%
  get_rocs() %>%
  plot_rocs_publish()

models_procan <- shap_results %>%
  select(starts_with("Model")) %>%
  distinct() %>%
  filter(Model.Dataset == "ProCan" & Model.Subset == "All" & !grepl("Log2FC", Model.Variant)) %>%
  get_rocs() %>%
  plot_rocs_publish()

models_cptac <- shap_results %>%
  select(starts_with("Model")) %>%
  distinct() %>%
  filter(Model.Dataset == "CPTAC") %>%
  get_rocs() %>%
  plot_rocs_publish()

panel_roc <- cowplot::plot_grid(models_depmap, models_procan, models_cptac,
                                labels = c("D", "E", "F"),
                                nrow = 1, ncol = 3)

# === Combine Panels into Figure ===
figure3 <- cowplot::plot_grid(panel_rank, panel_bootstrap, panel_roc,
                              labels = c("", "C", ""),
                              nrow = 3, ncol = 1, rel_heights = c(1, 0.5, 0.5))

cairo_pdf(here(plots_dir, "figure03.pdf"), width = 16, height = 20)
figure3
dev.off()