library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)
library(pROC)
library(caret)
library(randomForest)
library(xgboost)
library(shapr)
library(forcats)
library(openxlsx)


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
tables_dir <- tables_base_dir
dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

# === Load Datasets ===
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors_filtered.parquet"))

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_matched_renorm <- read_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'))
buf_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_wgd.parquet"))
buf_no_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_no-wgd.parquet"))
expr_buf_p0211 <- read_parquet(here(output_data_dir, 'expression_buffering_p0211.parquet'))
expr_buf_chunduri <- read_parquet(here(output_data_dir, 'expression_buffering_chunduri.parquet'))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))

expr_buf_eng <- bind_rows(expr_buf_p0211, expr_buf_chunduri) %>%
  mutate(Dataset = "Engin.",
         Gene.Chromosome = as.character(Gene.Chromosome))

# === Define Functions ===
prepare_datasets <- function(dataset, buffering_class_col, filter_func,
                             df_factors = dc_factors, factor_cols = dc_factor_cols,
                             training_set_ratio = 0.8, target_balance = 0.7, cv_eval = FALSE) {
  set.seed(42)

  # Filter and clean dataset before training
  df_prep <- dataset %>%
    add_factors(df_factors, factor_cols = factor_cols) %>%
    # ToDo: Only use imputation on training set
    # ToDo: Add shadow matrix - https://cran.r-project.org/web/packages/naniar/vignettes/getting-started-w-naniar.html
    impute_na(factor_cols = factor_cols) %>%
    normalize_features(method = "min-max", factor_cols = factor_cols) %>%
    filter_func() %>%
    clean_data({ { buffering_class_col } }, factor_cols = factor_cols) %>%
    select(where(~!all(is.na(.x)))) %>% # Remove empty factors
    # rebalance_binary(buffered, target_balance = target_balance)
    shuffle_rows()

  # Return NA if no training data remains
  if (floor(nrow(df_prep) * training_set_ratio) == 0) return(NA)

  # Don't create a test set if evaluation is done entirely using cross validation
  if (cv_eval == TRUE) return(list(training = df_prep, test = NULL))

  # Split training & test data
  return(split_train_test(df_prep, training_set_ratio))
}

run_analysis <- function(dataset, buffering_class_col, filter_func, model_name, train_control,
                         plots_dir, models_dir, sub_dir, prog_bar, factor_cols = dc_factor_cols,
                         training_set_ratio = 0.8, target_balance = 0.7, cv_eval = FALSE) {
  plots_dir <- here(plots_dir, sub_dir)
  dir.create(plots_dir, recursive = TRUE)
  dir.create(models_dir, recursive = TRUE)

  model_filename <- paste0("model_", model_name, "_", paste0(sub_dir, collapse = "_"), ".rds")

  # Prepare dataset for training and split into training and test set
  datasets <- prepare_datasets(dataset, { { buffering_class_col } }, filter_func,
                               factor_cols = factor_cols, training_set_ratio = training_set_ratio,
                               target_balance = target_balance, cv_eval = cv_eval)
  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  # Check if datasets are empty
  if (!is.list(datasets)) {
    setTxtProgressBar(prog_bar, prog_bar$getVal() + 2)
    return(NA)
  }

  # Train model
  # TODO: Consider stratified cross validation
  model <- caret::train(buffered ~ .,
                        data = datasets$training,
                        method = model_name,
                        trControl = train_control,
                        metric = "ROC")

  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  # Evaluate Model
  explain_model(model, plots_dir, filename = paste0("importance_", model_name, ".png"))
  model_eval <- evaluate_model(model, datasets$test, plots_dir,
                               filename = paste0("ROC-Curve_", model_name, ".png"), cv_eval = cv_eval)
  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  # Save Model
  model$datasets <- datasets  # Add datasets to model
  model$evaluation <- model_eval
  saveRDS(model, here(models_dir, model_filename))

  return(list(modelFileName = model_filename, roc = model_eval$roc))
}

# === Build Models ===
## Define training parameters
tc_base <- trainControl(method = "cv",
                        number = 3,
                        savePredictions = TRUE,
                        classProbs = TRUE,
                        verboseIter = FALSE,
                        returnData = FALSE,
                        allowParallel = TRUE,
                        summaryFunction = twoClassSummary)

tc_nn <- tc_base # ToDo: Increase training iterations

tc_p0211 <- tc_base
tc_p0211["number"] <- 10

## Define models to be trained
models <- list(
  list(modelName = "xgbLinear", tc = tc_base),
  list(modelName = "rf", tc = tc_base),
  list(modelName = "pcaNNet", tc = tc_nn)
)

## Define datasets to train models on
datasets <- list(
  list(dataset = expr_buf_procan, name = "ProCan", sourceDataset = "ProCan"),
  list(dataset = expr_buf_depmap, name = "DepMap", sourceDataset = "DepMap"),
  # list(dataset = expr_buf_matched_renorm, name = "MatchedRenorm", sourceDataset = "DepMap+ProCan"),
  list(dataset = buf_wgd, name = "DepMap-WGD", sourceDataset = "DepMap"),
  list(dataset = buf_no_wgd, name = "DepMap-NoWGD", sourceDataset = "DepMap"),
  # list(dataset = expr_buf_p0211, name = "P0211", cv_eval = TRUE, tc = tc_p0211, sourceDataset = "P0211"),
  list(dataset = expr_buf_eng, name = "Engineered", cv_eval = TRUE, tc = tc_p0211, sourceDataset = "Engineered"),
  list(dataset = expr_buf_cptac, name = "CPTAC", sourceDataset = "CPTAC")
)

## Define training data conditions
analysis_conditions <- list(
  # list(buffering = "Buffering.GeneLevel.Class", filter = identity, sub_dir =  list("Gene-Level", "Unfiltered")),
  # list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff, sub_dir =  list("Gene-Level", "Filtered")),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain, sub_dir =  list("Gene-Level", "Filtered_Gain")),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss, sub_dir =  list("Gene-Level", "Filtered_Loss")),
  list(buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain, sub_dir =  list("ChromosomeArm-Level", "Gain")),
  list(buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss, sub_dir =  list("ChromosomeArm-Level", "Loss")),
  list(buffering = "Buffering.ChrArmLevel.Log2FC.Class", filter = filter_arm_gain, sub_dir = list("ChromosomeArm-Level", "Gain_Log2FC")),
  list(buffering = "Buffering.ChrArmLevel.Log2FC.Class", filter = filter_arm_loss, sub_dir = list("ChromosomeArm-Level", "Loss_Log2FC")),
  list(buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain_gene_avg, sub_dir = list("ChromosomeArm-Level", "Gain_Average")),
  list(buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss_gene_avg, sub_dir = list("ChromosomeArm-Level", "Loss_Average"))
)

## Training & Evaluation Loop
### Test robustness of models by excluding certain factors
### ToDo: Test this separately
excluded_factors <- c("Homology Score", "Protein Abundance", "mRNA Abundance",
                      "Transcription Factors (Repressor)", "Transcription Factors (Activator)",
                      "Triplosensitivity")

analysis_summary <- data.frame()

pb <- txtProgressBar(min = 0,
                     max = length(models) * length(datasets) * length(analysis_conditions) * 3,
                     style = 3)
for (model in models) {
  for (dataset in datasets) {
    analysis_results <- list()
    for (analysis in analysis_conditions) {
      tc <- model$tc
      if (is.list(dataset$tc)) tc <- dataset$tc

      suppressMessages({
        suppressWarnings({
          # ToDo: Allow adjustment of training test ratio for dataset
          results <- run_analysis(dataset = dataset$dataset,
                                  buffering_class_col = get(analysis$buffering),
                                  factor_cols = dc_factor_cols[!dc_factor_cols %in% excluded_factors],
                                  filter_func = analysis$filter,
                                  model_name = model$modelName,
                                  train_control = tc,
                                  plots_dir = plots_dir,
                                  sub_dir = append(dataset$name, analysis$sub_dir),
                                  models_dir = models_base_dir,
                                  prog_bar = pb,
                                  cv_eval = isTRUE(dataset$cv_eval))
        }) })
      result_id <- paste0(analysis$sub_dir, collapse = "_")

      if (is.list(results)) {
        roc_auc <- NA
        if (!is.null(results$roc)) roc_auc <- as.numeric(auc(results$roc))

        analysis_results[[result_id]] <- results
        analysis_summary <- analysis_summary %>%
          bind_rows(data.frame(
            Model.Filename = results$modelFileName,
            Model.ROC.AUC = roc_auc,
            Model.BaseModel = model$modelName,
            Model.Variant = result_id,
            Model.Dataset = dataset$sourceDataset,
            Model.Subset = case_when(grepl("-WGD", dataset$name) ~ "WGD",
                                     grepl("-NoWGD", dataset$name) ~ "Non-WGD",
                                     TRUE ~ "All"),
            Model.Condition = if_else(grepl("Gain", result_id), "Gain", "Loss"),        # TODO: Rename to CNEvent
            Model.Level = if_else(grepl("Gene", result_id), "Gene", "Chromosome Arm"),  # TODO: Rename to CNSource
            Model.Samples = if_else(grepl("Average", result_id), "Averaged", "Unaveraged"),
            Model.BufferingMethod = if_else(grepl("Log2FC", result_id) | grepl("Average", result_id),
                                            "Log2FC", "BR")
          ))
      }
      else
        warning("No results for analysis! Model: ", model$modelName,
                ", Dataset: ", dataset$name, ", Analysis: ", result_id)
    }
    df_rocs <- analysis_results %>%
      lapply(\(x) x$roc) %>%
      rocs_to_df()

    df_rocs %>%
      plot_rocs() %>%
      save_plot(paste0(paste("roc-summary", model$modelName, dataset$name, sep = "_"), ".png"),
                width = 250)
  }
}
close(pb)

analysis_summary %>%
  write_parquet(here(output_data_dir, 'multivariate_model_results.parquet')) %>%
  write.xlsx(here(tables_dir, 'multivariate_model_results.xlsx'), colNames = TRUE)

## ToDo: Train a "universal" model (multiple conditions, or multiple datasets, or both)

# === SHAP Analysis ===
analysis_summary <- read_parquet(here(output_data_dir, 'multivariate_model_results.parquet'))

## Explain XGBoost models using SHAP values
shap_results <- list()
pb <- txtProgressBar(min = 0, max = length(datasets) * length(analysis_conditions), style = 3)
for (dataset in datasets) {
  for (analysis in analysis_conditions) {
    # shapr currently only supports xgboost of all the trained models
    model_name <- "xgbLinear"
    sub_dir <- append(dataset$name, analysis$sub_dir)
    model_filename <- paste0("model_", model_name, "_", paste0(sub_dir, collapse = "_"), ".rds")

    if (!file.exists(here(models_base_dir, model_filename))) {
      warning("Model ", model_filename, " does not exist!")
      setTxtProgressBar(pb, pb$getVal() + 1)
      next
    }

    shap_explanation <- readRDS(here(models_base_dir, model_filename)) %>%
      estimate_shap(n_samples = 300, n_combinations = 300)

    df_explanation <- shap_explanation %>%
      shap2df() %>%
      mutate(Model.Filename = model_filename)

    df_explanation %>%
      shap_plot() %>%
      save_plot(paste0("shap-explanation_", model_name, ".png"), dir = here(plots_dir, sub_dir))
    df_explanation %>%
      shap_corr_importance_plot() %>%
      save_plot(paste0("shap-corr-importance_", model_name, ".png"), dir = here(plots_dir, sub_dir))

    shap_results[[model_filename]] <- df_explanation
    setTxtProgressBar(pb, pb$getVal() + 1)
  }
}
close(pb)

### Create mapping for cleaned names
dc_factor_cols_mapping <- dc_factor_cols
names(dc_factor_cols_mapping) <- janitor::make_clean_names(dc_factor_cols)

### Write results to disk
shap_results <- shap_results %>%
  bind_rows() %>%
  left_join(y = analysis_summary, by = "Model.Filename",
            relationship = "many-to-one") %>%
  rename(DosageCompensation.Factor.ID = DosageCompensation.Factor) %>%
  mutate(DosageCompensation.Factor = sapply(DosageCompensation.Factor.ID, \(x) dc_factor_cols_mapping[[x]])) %>%
  write_parquet(here(output_data_dir, 'shap-analysis.parquet'), version = "2.6") %T>%
  write.xlsx(here(tables_base_dir, "shap-analysis.xlsx"), colNames = TRUE)
