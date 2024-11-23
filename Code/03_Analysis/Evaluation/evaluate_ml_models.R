library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(arrow)
library(limma)
library(pROC)
library(openxlsx)
library(shadowtext)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "evaluation.R"))

models_base_dir <- here("Output", "Models")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "FactorAnalysis", "Multivariate", "OOS-Evaluation")
tables_dir <- tables_base_dir
dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

# Load datasets
model_results <- read_parquet(here(output_data_dir, 'multivariate_model_results.parquet'))

# === Create summary statistics to evaluate performance of ML pipeline ===
model_statistics <- model_results %>%
  filter(Model.Dataset %in% c("DepMap", "ProCan", "CPTAC")) %>%
  group_by(Model.BaseModel, Model.Level, Model.Condition) %>%
  skimr::skim(Model.ROC.AUC) %>%
  write_parquet(here(output_data_dir, "multivariate_model_statistics.parquet")) %T>%
  write.xlsx(here(tables_dir, "multivariate_model_statistics.xlsx"), colNames = TRUE)

# === Out of Sample evaluation of trained models ===
oos_results <- data.frame()

# TODO: P0211 dataset has no 3'-/5'-UTR factor data required for running the ProCan model
# TODO: Limit datasets to chr13 for P0211 models

pb <- txtProgressBar(min = 0, max = nrow(model_results)^2, style = 3)
for (model_filename in model_results$Model.Filename) {
  model <- readRDS(here(models_base_dir, model_filename))
  for (dataset_filename in model_results$Model.Filename) {
    dataset_model <- readRDS(here(models_base_dir, dataset_filename))
    dataset <- bind_rows(dataset_model$datasets)

    tryCatch({
      suppressMessages({ suppressWarnings({ eval <- evaluate_model(model, dataset) }) })
    }, error = function(e) {
      warning(paste("Unable to evaluate model:", model_filename, "on dataset:", dataset_filename))
      next
    })

    oos_results <- list(oos = eval$roc) %>%
      rocs_to_df() %>%
      mutate(Model.Filename = model_filename,
             Dataset.Filename = dataset_filename) %>%
      bind_rows(oos_results)

    setTxtProgressBar(pb, pb$getVal() + 1)
  }
}
close(pb)

dataset_results <- model_results %>%
  rename_all(~str_replace(., "^Model", "Dataset"))

oos_results <- oos_results %>%
  write_parquet(here(output_data_dir, "out-of-sample-evaluation.parquet"))

oos_summary <- oos_results %>%
  rename(OOS.ROC.AUC = AUC) %>%
  select(OOS.ROC.AUC, Model.Filename, Dataset.Filename) %>%
  distinct() %>%
  left_join(y = model_results, by = "Model.Filename") %>%
  left_join(y = dataset_results, by = "Dataset.Filename") %>%
  write_parquet(here(output_data_dir, "out-of-sample-evaluation_summary.parquet")) %T>%
  write.xlsx(here(tables_dir, "out-of-sample-evaluation_summary.xlsx"), colNames = TRUE)

# === Plot results ===
oos_summary <- read_parquet(here(output_data_dir, "out-of-sample-evaluation_summary.parquet"))

oos_pan_cancer <- oos_summary %>%
  filter(Model.Condition == Dataset.Condition,
         Model.BaseModel == Dataset.BaseModel,    # Datasets should be equal across base models; avoid duplicates
         Model.Subset == Dataset.Subset,
         Model.Level == Dataset.Level,
         Model.Samples == Dataset.Samples,
         Model.BufferingMethod == Dataset.BufferingMethod,
         Model.Dataset != "Engin.",
         Dataset.Dataset != "Engin.",
         Model.Samples == "Unaveraged",
         Model.BufferingMethod == "BR",
         Model.Subset == "All",
         Model.BaseModel == "xgbLinear") %>%
  ggplot() +
  aes(y = Model.Dataset, x = Dataset.Dataset, fill = OOS.ROC.AUC,
      label = format(round(OOS.ROC.AUC, 2), nsmall = 2, scientific = FALSE)) +
  geom_tile() +
  geom_shadowtext(color = "white", bg.colour = default_color) +
  scale_fill_viridis_c(option = "magma", direction = 1,
                       limits = c(0.5, 1), oob = scales::squish) +
  facet_grid(Model.Level ~ Model.Condition)

oos_pan_cancer %>%
  save_plot("oos_pan-cancer_xgbLinear.png")
