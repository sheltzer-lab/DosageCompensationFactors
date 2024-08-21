library(here)
library(dplyr)
library(arrow)
library(ggplot2)
library(forcats)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "visualization.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "analysis.R"))
source(here("Code", "evaluation.R"))

models_base_dir <- here("Output", "Models")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "FactorAnalysis", "Multivariate")
dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

## SHAP Results
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

### Select models for visualization
shap_gain <- shap_results %>%
  filter(Model.Variant == "ProCan_ChromosomeArm-Level_Gain_Log2FC")
shap_loss <- shap_results %>%
  filter(Model.Variant == "ProCan_ChromosomeArm-Level_Loss_Log2FC")

shap_arrows_plot_gain <- shap_plot_arrows(shap_gain, category_lab = "Factor", show_legend = FALSE, title = "Chromosome Gain")
shap_arrows_plot_loss <- shap_plot_arrows(shap_loss, category_lab = NULL, color_lab = "Factor Value", title = "Chromosome Loss")

shap_raw_plots <- cowplot::plot_grid(shap_arrows_plot_gain, shap_arrows_plot_loss,
                                     nrow = 1, ncol = 2, align = "h", axis = "lr",
                                     rel_widths = c(0.79, 1), labels = c("A", ""))

shap_imp_plots <- cowplot::plot_grid(shap_corr_importance_plot(shap_gain, show_legend = FALSE),
                                     shap_corr_importance_plot(shap_loss, category_lab = NULL),
                                     nrow = 1, ncol = 2, align = "h", axis = "lr",
                                     rel_widths = c(0.75, 1), labels = c("B", ""))

shap_gain_plots <- cowplot::plot_grid(shap_arrows_plot_gain, shap_arrows_plot_loss,
                                     nrow = 1, ncol = 2, align = "h", axis = "lr",
                                     rel_widths = c(0.79, 1), labels = c("A", ""))

plot_publish_shap <- cowplot::plot_grid(shap_raw_plots, shap_imp_plots,
                                   nrow = 2, ncol = 1, rel_heights = c(1, 1))

cairo_pdf(here(plots_dir, "multivariate_shap_publish.pdf"), height = 15, width = 12)
plot_publish_shap
dev.off()

shap_heatmap <- bind_rows(shap_gain, shap_loss) %>%
  mutate(Model.Variant = str_remove(Model.Variant, "ProCan_")) %>%
  mutate(Model.Variant = str_remove(Model.Variant, "-Level")) %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, abs(SHAP.Factor.Corr), .desc = TRUE)) %>%
  simple_heatmap(DosageCompensation.Factor, Model.Variant, SHAP.Factor.Corr, Label,
                 x_lab = "Feature", y_lab = "Model", legend_lab = "SHAP-Value-Feature-Correlation")

plot_publish_shap_alt <- cowplot::plot_grid(shap_raw_plots, shap_heatmap + labs(x = NULL, y = NULL),
                                   nrow = 2, ncol = 1, rel_heights = c(1, 0.4), labels = c("", "B"))

cairo_pdf(here(plots_dir, "multivariate_shap_publish_alt.pdf"), height = 12, width = 11)
plot_publish_shap_alt
dev.off()

### Poster
cairo_pdf(here(plots_dir, "multivariate_shap_raw_poster.pdf"), height = 7, width = 11)
shap_raw_plots
dev.off()
cairo_pdf(here(plots_dir, "multivariate_shap_heatmap_poster.pdf"), height = 3.3, width = 12)
shap_heatmap + labs(x = NULL, y = NULL)
dev.off()

## Model Performance
model_rocs <- shap_results %>%
  distinct(Model.Filename, Model.Variant) %>%
  filter(grepl("ProCan", Model.Filename)) %>%
  mutate(Model.Variant = str_remove(Model.Variant, "ProCan_")) %>%
  mutate(Model.Variant = str_remove(Model.Variant, "-Level")) %>%
  group_by(Model.Variant) %>%
  group_map(\(entry, group) {
    result <- list()
    model <- readRDS(here(models_base_dir, entry$Model.Filename))
    result[[group$Model.Variant]] <- evaluate_model(model, model$datasets$test, plots_dir)$roc
    return(result)
  })

rocs_summary_xgbLinear <- model_rocs %>%
  flatten() %>%
  rocs_to_df() %>%
  plot_rocs(legend_position = "bottom", legend_rows = 5)

rf_gain_importance <- readRDS(here(models_base_dir, "model_rf_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds")) %>%
  explain_model(plots_dir)

plot_model <- cowplot::plot_grid(rocs_summary_xgbLinear, rf_gain_importance,
                                 nrow = 1, ncol = 2, labels = c("A", "B"), rel_widths = c(1, 0.98))

cairo_pdf(here(plots_dir, "multivariate_model_publish.pdf"), width = 12)
plot_model
dev.off()

## Model Performance (Poster)
### Selected Analysis Variants
selected_models <- c("ChromosomeArm_Gain_Log2FC", "ChromosomeArm_Loss_Log2FC",
                     "ChromosomeArm_Gain_Average", "ChromosomeArm_Loss_Average",
                     "Gene_Filtered_Gain", "Gene_Filtered_Loss")

df_model_rocs_selected <- model_rocs %>%
  flatten() %>%
  rocs_to_df() %>%
  filter(Name %in% selected_models) %>%
  mutate(Name = factor(Name, levels = selected_models))

rocs_summary_xgbLinear_selected <- df_model_rocs_selected %>%
  plot_rocs(legend_position = "bottom", legend_rows = 3, label_padding = 0.1)

### Model Architectures
df_gain_models <- data.frame(
  Model.Filename = c("model_rf_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds",
                     "model_pcaNNet_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds",
                     "model_xgbLinear_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds"),
  Model.Variant = c("Random Forests", "pcaNNet", "XGBLinear")
)

model_rocs_gain <- df_gain_models %>%
  group_by(Model.Variant) %>%
  group_map(\(entry, group) {
    result <- list()
    model <- readRDS(here(models_base_dir, entry$Model.Filename))
    result[[group$Model.Variant]] <- evaluate_model(model, model$datasets$test, plots_dir)$roc
    return(result)
  }) %>%
  flatten() %>%
  rocs_to_df() %>%
  mutate(Name = fct_reorder(Name, AUC, .desc = TRUE))

rocs_summary_gain <- model_rocs_gain %>%
  plot_rocs(legend_position = "bottom", legend_rows = 1, label_padding = 0.1)

### Export
plot_model_poster <- cowplot::plot_grid(rocs_summary_xgbLinear_selected +
                                          scale_color_manual(values = rev(categorical_color_pal[1:6])),
                                        rocs_summary_gain +
                                          scale_color_brewer(palette = "Dark2"),
                                        nrow = 1, ncol = 2, align = "hv", axis = "tblr")

cairo_pdf(here(plots_dir, "multivariate_model_publish_poster.pdf"), width = 10, height = 6)
plot_model_poster
dev.off()
