library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "publication.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir
models_base_dir <- here("Output", "Models")
tables_dir <- here(tables_base_dir, "Publication")

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)


# Load SHAP Results
model_results <- read_parquet(here(output_data_dir, 'multivariate_model_results.parquet'))
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))
oos_summary <- read_parquet(here(output_data_dir, "out-of-sample-evaluation_summary.parquet"))

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
           x = if_else(Model.Condition == "Gain", 0.2, 0),
           y = case_when(Model.Level == "Gene CN" ~ 0.2,
                         Model.Level == "ChrArm" ~ 0.1,
                         Model.Level == "ChrArm (avg.)" ~ 0),
           AUC = paste0(auc_prefix, format(round(AUC, 3), nsmall = 3)))

  df_legend <- data.frame(
    Label = c("ChrArm (avg.)", "ChrArm", "Gene CN", "Gain", "Loss"),
    x = c(0.4, 0.4, 0.4, 0.2, 0),
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
    scale_color_manual(values = color_palettes$ModelVariants) +
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

# === SHAP Correlation Heatmap Panel (per Dataset) ===
shap_heatmap_datasets <- shap_results %>%
  filter(Model.Level == "Gene" & Model.Subset == "All") %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, -SHAP.Factor.Corr, .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Model.Dataset, Model.Condition, Model.Variant, SHAP.Factor.Corr, Label) %>%
  nested_shap_heatmap(as.formula(Model.Dataset + Model.Condition ~ .))

# === Combine Panels into Figure ===
figure4 <- cowplot::plot_grid(panel_roc, shap_heatmap_datasets,
                              labels = c("", "D"),
                              nrow = 2, ncol = 1, rel_heights = c(0.5, 0.666))

cairo_pdf(here(plots_dir, "figure04.pdf"), width = 12, height = 9)
figure4
dev.off()

# === Tables ===
shap_corr <- shap_results %>%
  distinct(DosageCompensation.Factor, Model.Dataset, Model.Condition, Model.Variant,
           SHAP.Factor.Corr, SHAP.Factor.Corr.p, SHAP.Factor.Corr.p.adj, SHAP.Median.Absolute)

oos_summary_renamed <- oos_summary %>%
  rename_with(~str_replace(.x, "Model.", "ModelA."), starts_with("Model.")) %>%
  rename_with(~str_replace(.x, "Dataset.", "ModelB."), starts_with("Dataset."))

t6_field_descriptions <- c(
  "=== TABLES ===" = "",
  "Model Performance" = "Contains the performance (ROC AUC) of each multifactorial model in predicting protein buffering. Different models have been trained for different conditions (gene copy number gain, chromosome arm loss in WGD cells, etc.).",
  "Out-of-Sample Prediction" = "Contains the performance (ROC AUC) of ModelA in predicting the dataset (training + test set) associated to ModelB.",
  "SHAP-Correlation" = "Correlation between SHAP-values of a model for a given factor and set of observations and the corresponding values of the factor.",
  "=== COLUMNS ===" = "",
  "Model.Filename" = "Filename of the trained and stored model. Serves as an identifier of the model.",
  "Model.ROC.AUC" = "Predictive power of a model (ROC AUC) in classifying a protein as Buffered or Scaling.",
  "Model.BaseModel" = "Model architecture used for training the model. Either xgbLinear (linear models with XGBoost), rf (Random Forests), or pcaNNet (Neural Network with prior PCA feature transformation and reduction).",
  "Model.Variant" = "Describes the copy number variant used for filtering the dataset (e.g., ChromosomeArm-Level_Gain). A combination of the columns Model.Level and Model.Condition.",
  "Model.Dataset" = "Proteomics dataset used as a source for training and evaluating the model after preprocessing (variant filtering, training-test-split).",
  "Model.Subset" = "Indicates whether all samples, whole genome doubled (WGD) samples, or non-WGD samples have been used for training and evaluating the model.",
  "Model.Condition" = "Copy number condition considered for the model (either Gain or Loss).",
  "Model.Level" = "Level of copy number data used for the model (either Chromosome Arm or Gene Copy Number).",
  "Model.Samples" = "Describes whether the samples used for the model have been averaged prior to determining buffering status (as in the ChrArm (avg.) analysis variant using average Log2FC values) or not (GeneCN, ChrArm analysis variants).",
  "Model.BufferingMethod" = "Quantities and thresholds applied for determining buffering classes used in the model.",
  "SHAP.Factor.Corr" = "Spearman's correlation coefficient between SHAP-value of a feature in a model and the associated value of the feature. Negative correlations indicate a higher confidence of the model in predicting a protein as Buffered if the value of the feature is high (e.g., high gene essentiality).",
  "SHAP.Factor.Corr.p" = "P-value of the correlation coefficient.",
  "SHAP.Factor.Corr.p.adj" = "Benjamini-Hochberg-adjusted p-value of the correlation coefficient.",
  "SHAP.Median.Absolute" = "Median absolute SHAP-value of a factor in a model."
)

df_t6_fields <- data.frame(Column = names(t6_field_descriptions), Description = unname(t6_field_descriptions))

wb <- createWorkbook()
sheet_readme <- addWorksheet(wb, "README")
sheet_model_results <- addWorksheet(wb, "Model Performance")
sheet_oos <- addWorksheet(wb, "Out-of-Sample Prediction")
sheet_shap <- addWorksheet(wb, "SHAP-Correlation")
writeDataTable(wb = wb, sheet = sheet_readme, x = df_t6_fields)
writeDataTable(wb = wb, sheet = sheet_model_results, x = model_results)
writeDataTable(wb = wb, sheet = sheet_oos, x = oos_summary)
writeDataTable(wb = wb, sheet = sheet_shap, x = shap_corr)
saveWorkbook(wb, here(tables_dir, "supplementary_table6.xlsx"), overwrite = TRUE)

t7_description <- c("This table contains the raw SHAP-values and derived quantities generated using trained multifactorial models.",
                    "SHAP-values are provided for each model, factor of a model, and observation (protein in a sample).",
                    "Observations are a random subset of the test set associated with a model. Only observations with correct model predictions are considered.",
                    "See manuscript for method details.")

t7_field_descriptions <- c(
  "ID" = "ID of the protein (observation) used for calculating the SHAP value. Consists of a sample ID (specific to the dataset; e.g., DepMap/Sanger cell model ID) and a Uniprot accession.",
  "DosageCompensation.Factor" = "Factor used in a model for predicting protein buffering.",
  "Factor.Value" = "Value of the factor.",
  "Factor.Value.Relative" = "Min-max-scaled value of the factor.",
  "SHAP.Value" = "SHAP-value of a factor in a model for a given protein in a sample. Indicates the confidence of the model in predicting the protein as Buffered (negative value) or Scaling (positive) using the associated relative feature value.",
  "SHAP.p25.Absolute" = "First quartile of absolute SHAP-values of a factor in a model.",
  "SHAP.Median.Absolute" = "Median absolute SHAP-value of a factor in a model.",
  "SHAP.p75.Absolute" = "Third quartile of absolute SHAP-values of a factor in a model.",
  "SHAP.Factor.Corr" = "Spearman's correlation coefficient between SHAP-value of a feature in a model and the associated value of the feature. Negative correlations indicate a higher confidence of the model in predicting a protein as Buffered if the value of the feature is high (e.g., high gene essentiality).",
  "SHAP.Factor.Corr.p" = "P-value of the correlation coefficient.",
  "SHAP.Factor.Corr.p.adj" = "Benjamini-Hochberg-adjusted p-value of the correlation coefficient.",
  "Model.Filename" = "Filename of the trained and stored model. Serves as an identifier of the model.",
  "Model.ROC.AUC" = "Predictive power of a model (ROC AUC) in classifying a protein as Buffered or Scaling.",
  "Model.BaseModel" = "Model architecture used for training the model. Either xgbLinear (linear models with XGBoost), rf (Random Forests), or pcaNNet (Neural Network with prior PCA feature transformation and reduction).",
  "Model.Variant" = "Describes the copy number variant used for filtering the dataset (e.g., ChromosomeArm-Level_Gain). A combination of the columns Model.Level and Model.Condition.",
  "Model.Dataset" = "Proteomics dataset used as a source for training and evaluating the model after preprocessing (variant filtering, training-test-split).",
  "Model.Subset" = "Indicates whether all samples, whole genome doubled (WGD) samples, or non-WGD samples have been used for training and evaluating the model.",
  "Model.Condition" = "Copy number condition considered for the model (either Gain or Loss).",
  "Model.Level" = "Level of copy number data used for the model (either Chromosome Arm or Gene Copy Number).",
  "Model.Samples" = "Describes whether the samples used for the model have been averaged prior to determining buffering status (as in the ChrArm (avg.) analysis variant using average Log2FC values) or not (GeneCN, ChrArm analysis variants).",
  "Model.BufferingMethod" = "Quantities and thresholds applied for determining buffering classes used in the model."
)

sup_table7 <- shap_results %>%
  select(all_of(names(t7_field_descriptions))) %>%
  write_parquet(here(tables_dir, "supplementary_table7.parquet"), version = "2.4")

createParquetReadme(t7_description, t7_field_descriptions, title = "Supplementary Table 7 - README",
                    readme_path = here(tables_dir, "README_supplementary_table7.md"),
                    file_path = "supplementary_table7.parquet", parquet_version = "2.4")

# === Supplemental Figures ===
## Model performance heatmap
### Get model order (xgbLinear outperforms other models on all pan cancer datasets wrt. mean and max ROC AUC)
model_order <- model_results %>%
  filter(!(Model.Level == "Chromosome Arm" & Model.BufferingMethod == "Log2FC" & Model.Samples == "Unaveraged")) %>%
  filter(Model.Dataset != "Engineered") %>%
  summarize(maxROCAUC = max(Model.ROC.AUC), .by = Model.BaseModel) %>%
  arrange(maxROCAUC) %>%
  pull(Model.BaseModel)

panel_model_perf <- model_results %>%
  filter(!(Model.Level == "Chromosome Arm" & Model.BufferingMethod == "Log2FC" & Model.Samples == "Unaveraged")) %>%
  mutate(Model.Dataset = str_replace(Model.Dataset, "Engineered", "Engin."),
         Model.Dataset = factor(Model.Dataset, levels = dataset_order),
         Model.BaseModel = factor(Model.BaseModel, levels = model_order),
         Model.AnalysisVariant = case_when(
           Model.Level == "Gene" ~ "GeneCN",
           Model.Level == "Chromosome Arm" & Model.Samples == "Unaveraged" ~ "ChrArm",
           Model.Level == "Chromosome Arm" & Model.Samples == "Averaged" ~ "ChrArmAvg",
           TRUE ~ NA,
         )) %>%
  ggplot() +
  aes(x = "", y = Model.BaseModel, fill = Model.ROC.AUC) +
  geom_raster() +
  ggh4x::facet_nested(. ~ Model.Dataset + Model.Subset + Model.AnalysisVariant + Model.Condition) +
  # scale_x_discrete(expand = c(0, 0)) +
  # scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(option = "magma", direction = 1, end = 0.9, limits = c(0.5, 0.9), oob = scales::squish,
                       breaks = seq(0.5, 0.9, 0.1), labels = c("\u22640.5", "0.6", "0.7", "0.8", "0.9")) +
  labs(x = NULL, y = NULL, fill = "ROC AUC") +
  cowplot::theme_minimal_grid() +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color = "lightgrey"),
        axis.ticks.y = element_blank(),
        legend.key.width = unit(24, "points"),
        legend.position = "top",
        legend.justification.top = "right",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 1, vjust = 0.9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.margin = margin(0, 8, 0, 8))

## Out-of-sample prediction
panel_oos <- oos_summary %>%
  filter(Model.Condition == Dataset.Condition,
         Model.BaseModel == Dataset.BaseModel,    # Datasets should be equal across base models; avoid duplicates
         Model.Subset == Dataset.Subset,
         Model.Level == Dataset.Level,
         Model.Samples == Dataset.Samples,
         Model.BufferingMethod == Dataset.BufferingMethod,
         Model.Dataset != "Engineered",
         Dataset.Dataset != "Engineered",
         Model.Samples == "Unaveraged",
         Model.BufferingMethod == "BR",
         Model.Subset == "All",
         Model.BaseModel == "xgbLinear") %>%
  mutate(Model.Level = str_replace(Model.Level, "Gene", "Gene CN"),
         Dataset.Dataset = factor(Dataset.Dataset, levels = dataset_order),
         Model.Dataset = factor(Model.Dataset, levels = dataset_order)) %>%
  ggplot() +
  aes(y = Model.Dataset, x = Dataset.Dataset, fill = OOS.ROC.AUC,
      label = format(round(OOS.ROC.AUC, 2), nsmall = 2, scientific = FALSE)) +
  geom_tile() +
  geom_text(color = "white") +
  scale_fill_viridis_c(option = "magma", direction = 1, end = 0.9, limits = c(0.5, 0.9), oob = scales::squish,
                       breaks = seq(0.5, 0.9, 0.1), labels = c("\u22640.5", "0.6", "0.7", "0.8", "0.9")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(Model.Level ~ Model.Condition) +
  labs(fill = "ROC AUC", x = "Evaluation Dataset (training + test)", y = "Model Training Dataset") +
  theme(legend.position = "none")

## Raw SHAP Value Panel
### Select models for visualization
shap_gain <- shap_results %>%
  filter(Model.Variant == "Gene-Level_Filtered_Gain" & Model.Dataset == "CPTAC") %>%
  group_by(DosageCompensation.Factor) %>%
  #slice_sample(n = 200) %>%     # Simplify visualization
  ungroup()
shap_loss <- shap_results %>%
  filter(Model.Variant == "Gene-Level_Filtered_Loss" & Model.Dataset == "CPTAC") %>%
  group_by(DosageCompensation.Factor) %>%
  #slice_sample(n = 200) %>%     # Simplify visualization
  ungroup()

shap_plot_gain <- shap_gain %>%
  shap_plot(category_lab = "Factor", show_legend = FALSE, title = "Gene CN Gain")
shap_plot_loss <- shap_loss %>%
  shap_plot(category_lab = NULL, color_lab = "Factor\nValue", title = "Gene CN Loss")

shap_raw_plots <- cowplot::plot_grid(shap_plot_gain, shap_plot_loss,
                                     nrow = 1, ncol = 2, align = "h", axis = "lr",
                                     labels = c("D", "E"),
                                     rel_widths = c(0.85, 1))


## Combine figures
figure_s4_sub1 <- cowplot::plot_grid(NULL, panel_oos,
                                     ncol = 2, rel_widths = c(1, 0.65), labels = c("A", "C"))
figure_s4 <- cowplot::plot_grid(figure_s4_sub1, panel_model_perf, shap_raw_plots,
                                nrow = 3, rel_heights = c(1, 0.65, 1.5), labels = c("", "B", ""))

cairo_pdf(here(plots_dir, "figure_s4.pdf"), width = 14, height = 16)
figure_s4
dev.off()
