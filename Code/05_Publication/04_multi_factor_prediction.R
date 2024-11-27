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
output_data_dir <- output_data_base_dir
models_base_dir <- here("Output", "Models")


dir.create(plots_dir, recursive = TRUE)

# Load SHAP Results
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

# === Multivariate Model Performance Panel ===
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

plot_rocs_publish <- function(df_rocs, auc_prefix = "AUC = ", title = "") {
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
           AUC = paste0(auc_prefix, format(round(AUC, 3), nsmall = 3)))

  df_legend <- data.frame(
    Label = c("ChrArm", "ChrArm (avg.)", "Gene CN", "Gain", "Loss"),
    x = c(0.5, 0.5, 0.5, 0.25, 0),
    y = c(0.01, 0.11, 0.21, 0.31, 0.31)
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
              hjust = 1, vjust = 0) +
    annotate("text", x = 1, y = 1, hjust = 0, size = 5, label = title) +
    scale_x_reverse(limits = c(1, 0)) +
    labs(x = "Specificity", y = "Sensitivity", color = "Model") +
    theme(legend.position = "none")
}

models_depmap <- shap_results %>%
  select(starts_with("Model")) %>%
  distinct() %>%
  filter(Model.Dataset == "DepMap" & Model.Subset == "All" & !grepl("Log2FC", Model.Variant)) %>%
  get_rocs() %>%
  plot_rocs_publish(auc_prefix = "", title = "DepMap")

models_procan <- shap_results %>%
  select(starts_with("Model")) %>%
  distinct() %>%
  filter(Model.Dataset == "ProCan" & Model.Subset == "All" & !grepl("Log2FC", Model.Variant)) %>%
  get_rocs() %>%
  plot_rocs_publish(auc_prefix = "", title = "ProCan")

models_cptac <- shap_results %>%
  select(starts_with("Model")) %>%
  distinct() %>%
  filter(Model.Dataset == "CPTAC") %>%
  get_rocs() %>%
  plot_rocs_publish(auc_prefix = "", title = "CPTAC")

panel_roc <- cowplot::plot_grid(models_depmap, models_procan, models_cptac,
                                labels = c("A", "B", "C"),
                                nrow = 1, ncol = 3)

# === SHAP Value Panel ===
### Select models for visualization
shap_gain <- shap_results %>%
  filter(Model.Variant == "Gene-Level_Filtered_Gain" & Model.Dataset == "CPTAC") %>%
  group_by(DosageCompensation.Factor) %>%
  slice_sample(n = 200) %>%     # Simplify visualization
  ungroup()
shap_loss <- shap_results %>%
  filter(Model.Variant == "Gene-Level_Filtered_Loss" & Model.Dataset == "CPTAC") %>%
  group_by(DosageCompensation.Factor) %>%
  slice_sample(n = 200) %>%     # Simplify visualization
  ungroup()

shap_arrows_plot_gain <- shap_gain %>%
  shap_plot_arrows(category_lab = "Factor", show_legend = FALSE, title = "Gene CN Gain")
shap_arrows_plot_loss <- shap_loss %>%
  shap_plot_arrows(category_lab = NULL, color_lab = "Factor Value", title = "Gene CN Loss")

shap_raw_plots <- cowplot::plot_grid(shap_arrows_plot_gain, shap_arrows_plot_loss,
                                     nrow = 1, ncol = 2, align = "h", axis = "lr",
                                     labels = c("D", "E"),
                                     rel_widths = c(0.79, 1))

# === SHAP Correlation Heatmap Panel (per Dataset) ===
shap_heatmap_datasets <- shap_results %>%
  filter(Model.Level == "Gene" & Model.Subset == "All") %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, -SHAP.Factor.Corr, .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Model.Dataset, Model.Condition, Model.Variant, SHAP.Factor.Corr, Label) %>%
  nested_shap_heatmap(as.formula(Model.Dataset + Model.Condition ~ .))

# === Combine Panels into Figure ===
figure4 <- cowplot::plot_grid(panel_roc, shap_raw_plots, shap_heatmap_datasets,
                              labels = c("", "", "F"),
                              nrow = 3, ncol = 1, rel_heights = c(0.5, 1, 0.666))

cairo_pdf(here(plots_dir, "figure04.pdf"), width = 12, height = 17)
figure4
dev.off()
