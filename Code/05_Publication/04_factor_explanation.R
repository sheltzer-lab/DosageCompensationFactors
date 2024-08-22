library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

# Load SHAP Results
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

# === SHAP Value Panel ===
### Select models for visualization
shap_gain <- shap_results %>%
  filter(Model.Variant == "CPTAC_Gene-Level_Filtered_Gain")
shap_loss <- shap_results %>%
  filter(Model.Variant == "CPTAC_Gene-Level_Filtered_Loss")

shap_arrows_plot_gain <- shap_plot_arrows(shap_gain, category_lab = "Factor", show_legend = FALSE, title = "Gene CN Gain")
shap_arrows_plot_loss <- shap_plot_arrows(shap_loss, category_lab = NULL, color_lab = "Factor Value", title = "Gene CN Loss")

shap_raw_plots <- cowplot::plot_grid(shap_arrows_plot_gain, shap_arrows_plot_loss,
                                     nrow = 1, ncol = 2, align = "h", axis = "lr",
                                     rel_widths = c(0.79, 1))

# === SHAP Correlation Heatmap Panel (per Dataset) ===
shap_heatmap_datasets <- shap_results %>%
  filter(Model.Level == "Gene" & Model.Subset == "All") %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, -SHAP.Factor.Corr, .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Model.Dataset, Model.Condition, Model.Variant, SHAP.Factor.Corr, Label) %>%
  nested_shap_heatmap(as.formula(Model.Dataset + Model.Condition ~ .))

