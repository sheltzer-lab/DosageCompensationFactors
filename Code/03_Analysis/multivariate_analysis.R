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


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "visualization.R"))
source(here("Code", "03_Analysis", "analysis.R"))

models_base_dir <- here("Output", "Models")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Multivariate")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_matched_renorm <- read_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'))
buf_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_wgd.parquet"))
buf_no_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_no-wgd.parquet"))

# === Define Functions ===

shuffle_rows <- function(df) {
  df[sample(nrow(df), replace = FALSE),]
}

# ToDo: Explore further imputation methods
impute_na <- function(df) {
  df %>%
    mutate_if(is.numeric, \(x) replace_na(x, median(x, na.rm = TRUE)))
}

balanced_sample <- function(df, class_col, n_per_class) {
  df %>%
    group_by({ { class_col } }) %>%
    slice_sample(n = n_per_class) %>%
    ungroup()
}

# TODO: Explore further techniques
# Oversampling, Undersampling, slice_sample(weight_by=..., replace=TRUE)
rebalance_binary <- function(df, class_col, target_balance = 0.5) {
  total_size <- nrow(df)
  new_samples <- df %>%
    count({ { class_col } }) %>%
    mutate(Ratio = n/total_size) %>%
    mutate(Samples = c(ceiling((target_balance * min(n)) / (1 - target_balance)), min(n))) %>%
    select(-n, -Ratio)

  df %>%
    left_join(y = new_samples, by = quo_name(enquo(class_col))) %>%
    group_by({ { class_col } }) %>%
    group_map(~slice_sample(.x, n = unique(.x$Samples)), .keep = TRUE) %>%
    bind_rows() %>%
    select(-Samples)
}

explain_model <- function(model, dir, filename = "importance.png") {
  # See randomForest::importance() and xgboost::xgb.importance() (rescaled, 0-100)
  ## rf: Mean Decrease in Gini Impurity
  ## xgbLinear: Weight of linear coefficients of feature
  ## xgbTree: fractional contribution of each feature
  if (model$method %in% c("rf", "xgbLinear", "xgbTree")) {
    caret::varImp(model)$importance %>%
      rename(Importance = "Overall") %>%
      tibble::rownames_to_column(var = "Factor") %>%
      mutate(Factor = factor(Factor, levels = Factor[order(Importance)])) %>%
      arrange(Factor) %>%
      vertical_bar_chart(Factor, Importance,
                         value_range = c(0, 100), bar_label_shift = 8,
                         line_intercept = 0, break_steps = 10,
                         category_lab = "Factor", value_lab = "Relative Importance",
                         title = paste0("Model Importance (", model$method, ")")) %>%
      save_plot(filename, dir = dir)
  }
}

evaluate_model <- function(model, test_set, dir, filename = "ROC-Curve.png") {
  test_predicted_prob <- predict(model, test_set, type = "prob")
  model_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob[, "Buffered"]), na.rm = TRUE)
  png(here(dir, filename),
      width = 200, height = 200, units = "mm", res = 200)
  plot(model_roc, print.thres = "best", print.thres.best.method = "closest.topleft",
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.1)
  dev.off()

  return(model_roc)
}

clean_data <- function (dataset, buffering_class_col, factor_cols = dc_factor_cols) {
  dataset %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = factor({ { buffering_class_col } }, levels = c("Scaling", "Buffered"))) %>%
    mutate_if(is.numeric, \(x) as.numeric(x)) %>%
    drop_na(Buffered) %>%
    select(Buffered, all_of(factor_cols)) %>%
    janitor::clean_names()
}

prepare_datasets <- function(dataset, buffering_class_col, filter_func,
                             df_factors = dc_factors, factor_cols = dc_factor_cols,
                             training_set_ratio = 0.8, target_balance = 0.7) {
  set.seed(42)

  # Filter and clean dataset before training
  df_prep <- dataset %>%
    filter_func() %>%
    add_factors(df_factors, factor_cols = factor_cols) %>%
    clean_data({ { buffering_class_col } }, factor_cols = factor_cols) %>%
    # ToDo: Only use imputation on training set
    impute_na() %>%
    # rebalance_binary(buffered, target_balance = target_balance) %>%
    shuffle_rows()

  # Split training & test data
  train_size <- ceiling(nrow(df_prep) * training_set_ratio)
  test_size <- nrow(df_prep) - train_size
  datasets <- list(training = df_prep[1:train_size,],
                   test = df_prep[(train_size + 1):(train_size + test_size),])

  return(datasets)
}

run_analysis <- function(dataset, buffering_class_col, filter_func, model_name, train_control,
                         plots_dir, models_dir, sub_dir, prog_bar,
                         factor_cols = dc_factor_cols, training_set_ratio = 0.8, target_balance = 0.7) {
  plots_dir <- here(plots_dir, sub_dir)
  dir.create(plots_dir, recursive = TRUE)
  dir.create(models_dir, recursive = TRUE)

  # Prepare dataset for training and split into training and test set
  datasets <- prepare_datasets(dataset, { { buffering_class_col } }, filter_func,
                               factor_cols = factor_cols, training_set_ratio = 0.8, target_balance = 0.7)
  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  # Train model
  model <- caret::train(buffered ~ .,
                        data = datasets$training,
                        method = model_name,
                        trControl = train_control,
                        metric = "ROC")

  # Add datasets to model
  model$datasets <- datasets

  # Save Model
  model_filename <- paste0("model_", model_name, "_", paste0(sub_dir, collapse = "_"), ".rds")
  saveRDS(model, here(models_dir, model_filename))
  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  # Evaluate Model
  explain_model(model, plots_dir, filename = paste0("importance_", model_name, ".png"))
  model_roc <- evaluate_model(model, datasets$test, plots_dir, filename = paste0("ROC-Curve_", model_name, ".png"))
  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  return(list(modelFileName = model_filename, roc = model_roc))
}

# === Build Models ===
## Define training parameters
tc_base <- trainControl(method = "cv",
                        number = 2,
                        savePredictions = TRUE,
                        classProbs = TRUE,
                        verboseIter = FALSE,
                        returnData = FALSE,
                        allowParallel = TRUE,
                        summaryFunction = twoClassSummary)

tc_nn <- tc_base # ToDo: Increase training iterations
tc_nn["number"] <- 3

## Define models to be trained
models <- list(
  list(modelName = "xgbLinear", tc = tc_base),
  list(modelName = "rf", tc = tc_base),
  list(modelName = "pcaNNet", tc = tc_nn)
)

## Define datasets to train models on
datasets <- list(
  list(dataset = expr_buf_procan, name = "ProCan"),
  list(dataset = expr_buf_depmap, name = "DepMap"),
  list(dataset = expr_buf_matched_renorm, name = "MatchedRenorm"),
  list(dataset = buf_wgd, name = "DepMap-WGD"),
  list(dataset = buf_no_wgd, name = "DepMap-NoWGD")
)

## Define training data conditions
analysis_conditions <- list(
  # list(buffering = "Buffering.GeneLevel.Class", filter = identity, sub_dir =  list("Gene-Level", "Unfiltered")),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles, sub_dir =  list("Gene-Level", "Filtered")),
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
excluded_factors <- c("Homology Score")

pb <- txtProgressBar(min = 0,
                     max = length(models) * length(datasets) * length(analysis_conditions) * 3,
                     style = 3)
i <- 0
for (model in models) {
  for (dataset in datasets) {
    analysis_results <- list()
    for (analysis in analysis_conditions) {
      suppressMessages({
        suppressWarnings({
          results <- run_analysis(dataset = dataset$dataset,
                                  buffering_class_col = get(analysis$buffering),
                                  factor_cols = dc_factor_cols[!dc_factor_cols %in% excluded_factors],
                                  filter_func = analysis$filter,
                                  model_name = model$modelName,
                                  train_control = model$tc,
                                  plots_dir = plots_dir,
                                  sub_dir = append(dataset$name, analysis$sub_dir),
                                  models_dir = models_base_dir,
                                  prog_bar = pb)
        }) })
      result_id <- paste0(analysis$sub_dir, collapse = "_")
      analysis_results[[result_id]] <- results
    }
    rocs <- lapply(analysis_results, \(x) x$roc)
    rocs_to_df(rocs) %>%
      plot_rocs() %>%
      save_plot(paste0(paste("roc-summary", model$modelName, dataset$name, sep = "_"), ".png"),
                width = 250)
  }
}
close(pb)


## ToDo: Train and compare a set of models on one dataset (or multiple)
## ToDo: Train a "universal" model (multiple conditions, or multiple datasets, or both)

# Use ProCan gain model for evaluating DepMap gene-level gain data
model_procan_arm_gain_rf <- readRDS(here(models_base_dir, "model_rf_ProCan_ChromosomeArm-Level_Gain.rds"))
model_procan_arm_gain_pcaNN <- readRDS(here(models_base_dir, "model_pcaNNet_ProCan_ChromosomeArm-Level_Gain.rds"))
model_procan_arm_gain_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_ProCan_ChromosomeArm-Level_Gain.rds"))
model_matched_renorm_arm_gain_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_MatchedRenorm_ChromosomeArm-Level_Gain.rds"))

test_data_depmap <- expr_buf_depmap %>%
  filter_arm_gain() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Class) %>%
  impute_na() %>%
  shuffle_rows()

eval_dir <- here(plots_dir, "ProCan", "ChromosomeArm-Level", "Gain", "Eval-DepMap")
dir.create(eval_dir, recursive = TRUE)
evaluate_model(model_procan_arm_gain_rf, test_data_depmap, eval_dir, filename = "ROC-Curve_rf.png")
evaluate_model(model_procan_arm_gain_pcaNN, test_data_depmap, eval_dir, filename = "ROC-Curve_pcaNNet.png")
evaluate_model(model_procan_arm_gain_xgb, test_data_depmap, eval_dir, filename = "ROC-Curve_xgbLinear.png")

eval_dir <- here(plots_dir, "MatchedRenorm", "ChromosomeArm-Level", "Gain", "Eval")
dir.create(eval_dir, recursive = TRUE)

test_data_procan <- expr_buf_procan %>%
  filter_arm_gain() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Class) %>%
  impute_na() %>%
  shuffle_rows()

evaluate_model(model_matched_renorm_arm_gain_xgb, test_data_depmap, eval_dir, filename = "ROC-Curve_xgbLinear_DepMap.png")
evaluate_model(model_matched_renorm_arm_gain_xgb, test_data_procan, eval_dir, filename = "ROC-Curve_xgbLinear_ProCan.png")
