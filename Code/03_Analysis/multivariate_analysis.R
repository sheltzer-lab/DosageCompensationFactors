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
source(here("Code", "analysis.R"))

models_base_dir <- here("Output", "Models")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "FactorAnalysis", "Multivariate")
dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_matched_renorm <- read_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'))
buf_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_wgd.parquet"))
buf_no_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_no-wgd.parquet"))
expr_buf_p0211 <- read_parquet(here(output_data_dir, 'expression_buffering_p0211.parquet'))

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

explain_model <- function(model, dir, filename = NULL) {
  # See randomForest::importance() and xgboost::xgb.importance() (rescaled, 0-100)
  ## rf: Mean Decrease in Gini Impurity
  ## xgbLinear: Weight of linear coefficients of feature
  ## xgbTree: fractional contribution of each feature
  if (model$method %in% c("rf", "xgbLinear", "xgbTree")) {
    plot <- caret::varImp(model)$importance %>%
      rename(Importance = "Overall") %>%
      tibble::rownames_to_column(var = "Factor") %>%
      mutate(Factor = factor(Factor, levels = Factor[order(Importance)])) %>%
      arrange(Factor) %>%
      vertical_bar_chart(Factor, Importance,
                         value_range = c(0, 100), bar_label_shift = 1,
                         line_intercept = 0, break_steps = 10,
                         category_lab = "Factor", value_lab = "Relative Importance",
                         title = paste0("Model Importance (", model$method, ")"))

    if (is.null(filename))
      return(plot)

    save_plot(plot, filename, dir = dir)
  }
}

evaluate_model <- function(model, test_set, dir, filename = NULL, cv_eval = FALSE) {
  if (cv_eval == TRUE) {
    # Use results from cross validation generated during training for evaluation
    # Get ROC curve across all folds of k-fold CV for model with best parameters
    cv_best_preds <- model$pred %>%
      semi_join(y = model$bestTune)

    model_roc <- roc(response = cv_best_preds$obs, predictor = cv_best_preds$Buffered, na.rm = TRUE)
  } else {
    # Use separate test set for evaluation
    test_predicted_prob <- predict(model, test_set, type = "prob")
    model_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob[, "Buffered"]), na.rm = TRUE)

    # Add model perfomance metrics
    lvl <- c("Buffered", "Scaling")
    test_predicted_reponse <- test_predicted_prob %>%
      mutate(
        Prediction = factor(if_else(Buffered > 0.5, "Buffered", "Scaling"), levels = lvl),
        Response = factor(test_set$buffered, levels = lvl)
      )
    model_roc["performanceMetrics"] <- list(
      precision = precision(data = test_predicted_reponse$Prediction, reference = test_predicted_reponse$Response),
      recall = recall(data = test_predicted_reponse$Prediction, reference = test_predicted_reponse$Response),
      F1 = F_meas(data = test_predicted_reponse$Prediction, reference = test_predicted_reponse$Response)
    )
  }

  if (is.null(filename))
    return(model_roc)

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

normalize_min_max <- function(x, ...) {
    return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}

normalize_features <- function (df, method = "min-max", factor_cols = dc_factor_cols) {
  if (method == "min-max") {
    df %>%
      mutate_at(factor_cols, \(x) normalize_min_max(x, na.rm = TRUE))
  } else if (method == "z-score") {
    df %>%
      mutate_at(factor_cols, ~(scale(.) %>% as.vector))
  } else {
    warning("Invalid normalization method! Data frame passed as is. Supported methods: min-max, z-score")
    return(df)
  }
}

split_train_test <- function(df_prep, training_set_ratio) {
  df_prep <- df_prep  %>%
    mutate(id = row_number())

  training_set <- df_prep %>%
    slice_sample(by = buffered, prop = training_set_ratio)

  test_set <- df_prep %>%
    anti_join(training_set, by = 'id')

  return(list(training = training_set %>% select(-id),
              test = test_set %>% select(-id)))
}

prepare_datasets <- function(dataset, buffering_class_col, filter_func,
                             df_factors = dc_factors, factor_cols = dc_factor_cols,
                             training_set_ratio = 0.8, target_balance = 0.7, cv_eval = FALSE) {
  set.seed(42)

  # Filter and clean dataset before training
  df_prep <- dataset %>%
    filter_func() %>%
    add_factors(df_factors, factor_cols = factor_cols) %>%
    normalize_features(method = "min-max", factor_cols = factor_cols) %>%
    clean_data({ { buffering_class_col } }, factor_cols = factor_cols) %>%
    # ToDo: Only use imputation on training set
    impute_na() %>%
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

  # Save Model
  model$datasets <- datasets  # Add datasets to model
  saveRDS(model, here(models_dir, model_filename))
  setTxtProgressBar(prog_bar, prog_bar$getVal() + 1)

  # Evaluate Model
  explain_model(model, plots_dir, filename = paste0("importance_", model_name, ".png"))
  model_roc <- evaluate_model(model, datasets$test, plots_dir,
                              filename = paste0("ROC-Curve_", model_name, ".png"), cv_eval = cv_eval)
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
  list(dataset = expr_buf_procan, name = "ProCan"),
  list(dataset = expr_buf_depmap, name = "DepMap"),
  list(dataset = expr_buf_matched_renorm, name = "MatchedRenorm"),
  list(dataset = buf_wgd, name = "DepMap-WGD"),
  list(dataset = buf_no_wgd, name = "DepMap-NoWGD"),
  list(dataset = expr_buf_p0211, name = "P0211", cv_eval = TRUE, tc = tc_p0211)
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
      if (is.list(results))
        analysis_results[[result_id]] <- results
      else
        warning("No results for analysis! Model: ", model$modelName,
                ", Dataset: ", dataset$name, ", Analysis: ", result_id)
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

# === Evaluation ===
model_procan_gain_rf <- readRDS(here(models_base_dir, "model_rf_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_procan_gain_pcaNN <- readRDS(here(models_base_dir, "model_pcaNNet_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_procan_gain_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_ProCan_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_matched_renorm_gain_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_MatchedRenorm_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_p0211_gain_rf <- readRDS(here(models_base_dir, "model_rf_P0211_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_p0211_gain_pcaNN <- readRDS(here(models_base_dir, "model_pcaNNet_P0211_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_p0211_gain_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_P0211_ChromosomeArm-Level_Gain_Log2FC.rds"))
model_procan_loss_rf <- readRDS(here(models_base_dir, "model_rf_ProCan_ChromosomeArm-Level_Loss_Log2FC.rds"))
model_procan_loss_pcaNN <- readRDS(here(models_base_dir, "model_pcaNNet_ProCan_ChromosomeArm-Level_Loss_Log2FC.rds"))
model_procan_loss_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_ProCan_ChromosomeArm-Level_Loss_Log2FC.rds"))
model_matched_renorm_loss_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_MatchedRenorm_ChromosomeArm-Level_Loss_Log2FC.rds"))
model_p0211_loss_rf <- readRDS(here(models_base_dir, "model_rf_P0211_ChromosomeArm-Level_Loss_Log2FC.rds"))
model_p0211_loss_pcaNN <- readRDS(here(models_base_dir, "model_pcaNNet_P0211_ChromosomeArm-Level_Loss_Log2FC.rds"))
model_p0211_loss_xgb <- readRDS(here(models_base_dir, "model_xgbLinear_P0211_ChromosomeArm-Level_Loss_Log2FC.rds"))

## Use ProCan gain model for evaluating DepMap gain data
eval_dir <- here(plots_dir, "ProCan", "ChromosomeArm-Level", "Gain_Log2FC", "Eval-DepMap")
dir.create(eval_dir, recursive = TRUE)

test_data_depmap <- expr_buf_depmap %>%
  filter_arm_gain() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows()

rocs_procan_depmap_gain <- list(
  rf = evaluate_model(model_procan_gain_rf, test_data_depmap, eval_dir, filename = "ROC-Curve_rf.png"),
  pcaNNet = evaluate_model(model_procan_gain_pcaNN, test_data_depmap, eval_dir, filename = "ROC-Curve_pcaNNet.png"),
  xgbLinear = evaluate_model(model_procan_gain_xgb, test_data_depmap, eval_dir, filename = "ROC-Curve_xgbLinear.png")
)

rocs_to_df(rocs_procan_depmap_gain) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

## Use Matched (ProCan + DepMap, renormalized) gain model for evaluating DepMap & ProCan gain data
eval_dir <- here(plots_dir, "MatchedRenorm", "ChromosomeArm-Level", "Gain_Log2FC", "Eval")
dir.create(eval_dir, recursive = TRUE)

test_data_procan <- expr_buf_procan %>%
  filter_arm_gain() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows()

rocs_renorm_eval_gain <- list(
  DepMap = evaluate_model(model_matched_renorm_gain_xgb, test_data_depmap, eval_dir, filename = "ROC-Curve_xgbLinear_DepMap.png"),
  ProCan = evaluate_model(model_matched_renorm_gain_xgb, test_data_procan, eval_dir, filename = "ROC-Curve_xgbLinear_ProCan.png")
)

rocs_to_df(rocs_renorm_eval_gain) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

## Use ProCan gain model to evaluate P0211 gain
eval_dir <- here(plots_dir, "ProCan", "ChromosomeArm-Level", "Gain_Log2FC", "Eval-P0211")
dir.create(eval_dir, recursive = TRUE)

test_data_p0211 <- expr_buf_p0211 %>%
  filter_arm_gain() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows() %>%
  # Hotfix: P0211 has no UTR factor data required for ProCan model
  mutate(mean_3_utr_length = 0,
         mean_5_utr_length = 0)

rocs_procan_p0211_gain <- list(
  rf = evaluate_model(model_procan_gain_rf, test_data_p0211, eval_dir, filename = "ROC-Curve_rf.png"),
  pcaNNet = evaluate_model(model_procan_gain_pcaNN, test_data_p0211, eval_dir, filename = "ROC-Curve_pcaNNet.png"),
  xgbLinear = evaluate_model(model_procan_gain_xgb, test_data_p0211, eval_dir, filename = "ROC-Curve_xgbLinear.png")
)

rocs_to_df(rocs_procan_p0211_gain) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

## Use P0211 gain model to evaluate ProCan gain
eval_dir <- here(plots_dir, "P0211", "ChromosomeArm-Level", "Gain_Log2FC", "Eval-ProCan")
dir.create(eval_dir, recursive = TRUE)

rocs_p0211_procan_gain <- list(
  rf = evaluate_model(model_p0211_gain_rf, test_data_procan, eval_dir, filename = "ROC-Curve_rf.png"),
  pcaNNet = evaluate_model(model_p0211_gain_pcaNN, test_data_procan, eval_dir, filename = "ROC-Curve_pcaNNet.png"),
  xgbLinear = evaluate_model(model_p0211_gain_xgb, test_data_procan, eval_dir, filename = "ROC-Curve_xgbLinear.png")
)

rocs_to_df(rocs_procan_p0211_gain) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

### Only evaluate genes on chromosome 13 (P0211 only has Chromosome CNAs on 13)
test_data_procan_chr13 <- expr_buf_procan %>%
  filter_arm_gain() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  filter(Gene.Chromosome == 13) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows()

rocs_p0211_procan_gain_chr13 <- list(
  rf = evaluate_model(model_p0211_gain_rf, test_data_procan_chr13, eval_dir, filename = "ROC-Curve_rf_chr13.png"),
  pcaNNet = evaluate_model(model_p0211_gain_pcaNN, test_data_procan_chr13, eval_dir, filename = "ROC-Curve_pcaNNet_chr13.png"),
  xgbLinear = evaluate_model(model_p0211_gain_xgb, test_data_procan_chr13, eval_dir, filename = "ROC-Curve_xgbLinear_chr13.png")
)

rocs_to_df(rocs_p0211_procan_gain_chr13) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)


## Use ProCan loss model for evaluating DepMap loss data
eval_dir <- here(plots_dir, "ProCan", "ChromosomeArm-Level", "Loss_Log2FC", "Eval-DepMap")
dir.create(eval_dir, recursive = TRUE)

test_data_depmap <- expr_buf_depmap %>%
  filter_arm_loss() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows()

rocs_procan_depmap_loss <- list(
  rf = evaluate_model(model_procan_loss_rf, test_data_depmap, eval_dir, filename = "ROC-Curve_rf.png"),
  pcaNNet = evaluate_model(model_procan_loss_pcaNN, test_data_depmap, eval_dir, filename = "ROC-Curve_pcaNNet.png"),
  xgbLinear = evaluate_model(model_procan_loss_xgb, test_data_depmap, eval_dir, filename = "ROC-Curve_xgbLinear.png")
)

rocs_to_df(rocs_procan_depmap_loss) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

## Use Matched (ProCan + DepMap, renormalized) loss model for evaluating DepMap & ProCan loss data
eval_dir <- here(plots_dir, "MatchedRenorm", "ChromosomeArm-Level", "Loss_Log2FC", "Eval")
dir.create(eval_dir, recursive = TRUE)

test_data_procan <- expr_buf_procan %>%
  filter_arm_loss() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows()

rocs_renorm_eval_loss <- list(
  DepMap = evaluate_model(model_matched_renorm_loss_xgb, test_data_depmap, eval_dir, filename = "ROC-Curve_xgbLinear_DepMap.png"),
  ProCan = evaluate_model(model_matched_renorm_loss_xgb, test_data_procan, eval_dir, filename = "ROC-Curve_xgbLinear_ProCan.png")
)

rocs_to_df(rocs_renorm_eval_loss) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

## Use ProCan loss model to evaluate P0211 loss
eval_dir <- here(plots_dir, "ProCan", "ChromosomeArm-Level", "Loss_Log2FC", "Eval-P0211")
dir.create(eval_dir, recursive = TRUE)

test_data_p0211 <- expr_buf_p0211 %>%
  filter_arm_loss() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows() %>%
  # Hotfix: P0211 has no UTR factor data required for ProCan model
  mutate(mean_3_utr_length = 0,
         mean_5_utr_length = 0)

rocs_procan_p0211_loss <- list(
  rf = evaluate_model(model_procan_loss_rf, test_data_p0211, eval_dir, filename = "ROC-Curve_rf.png"),
  pcaNNet = evaluate_model(model_procan_loss_pcaNN, test_data_p0211, eval_dir, filename = "ROC-Curve_pcaNNet.png"),
  xgbLinear = evaluate_model(model_procan_loss_xgb, test_data_p0211, eval_dir, filename = "ROC-Curve_xgbLinear.png")
)

rocs_to_df(rocs_procan_p0211_loss) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

## Use P0211 loss model to evaluate ProCan loss
eval_dir <- here(plots_dir, "P0211", "ChromosomeArm-Level", "Loss_Log2FC", "Eval-ProCan")
dir.create(eval_dir, recursive = TRUE)

rocs_p0211_procan_loss <- list(
  rf = evaluate_model(model_p0211_loss_rf, test_data_procan, eval_dir, filename = "ROC-Curve_rf.png"),
  pcaNNet = evaluate_model(model_p0211_loss_pcaNN, test_data_procan, eval_dir, filename = "ROC-Curve_pcaNNet.png"),
  xgbLinear = evaluate_model(model_p0211_loss_xgb, test_data_procan, eval_dir, filename = "ROC-Curve_xgbLinear.png")
)

rocs_to_df(rocs_procan_p0211_loss) %>%
  plot_rocs() %>%
  save_plot("roc-summary.png", dir = eval_dir, width = 250)

### Only evaluate genes on chromosome 13 (P0211 only has Chromosome CNAs on 13)
test_data_procan_chr13 <- expr_buf_procan %>%
  filter_arm_loss() %>%
  add_factors(dc_factors, factor_cols = dc_factor_cols) %>%
  normalize_features(method = "min-max", factor_cols = dc_factor_cols) %>%
  filter(Gene.Chromosome == 13) %>%
  clean_data(Buffering.ChrArmLevel.Log2FC.Class) %>%
  impute_na() %>%
  shuffle_rows()

rocs_p0211_procan_loss_chr13 <- list(
  rf = evaluate_model(model_p0211_loss_rf, test_data_procan_chr13, eval_dir, filename = "ROC-Curve_rf_chr13.png"),
  pcaNNet = evaluate_model(model_p0211_loss_pcaNN, test_data_procan_chr13, eval_dir, filename = "ROC-Curve_pcaNNet_chr13.png"),
  xgbLinear = evaluate_model(model_p0211_loss_xgb, test_data_procan_chr13, eval_dir, filename = "ROC-Curve_xgbLinear_chr13.png")
)

rocs_to_df(rocs_p0211_procan_loss_chr13) %>%
  plot_rocs() %>%
  save_plot("roc-summary_chr13.png", dir = eval_dir, width = 250)

## Explain XGBoost models using SHAP values
estimate_shap <- function(model, n_samples = 100, n_combinations = 250, method = "ctree") {
  set.seed(42)

  # Prepare the data for explanation
  if (is.null(model$datasets$test)) {
    x_small <- model$datasets$training %>%
      group_by(buffered) %>%
      slice_sample(n = 2 * n_samples) %>%
      ungroup() %>%
      split_train_test(training_set_ratio = 0.5)

    training_x_small <- x_small$training
    test_x_small <- x_small$test
  } else {
    training_x_small <- model$datasets$training %>%
      group_by(buffered) %>%
      slice_sample(n = n_samples) %>%
      ungroup() %>%
      select(-buffered)
    test_x_small <- model$datasets$test %>%
      group_by(buffered) %>%
      slice_sample(n = n_samples) %>%
      ungroup() %>%
      select(-buffered)
  }

  # https://cran.r-project.org/web/packages/shapr/vignettes/understanding_shapr.html
  explainer <- shapr(training_x_small, model$finalModel, n_combinations = n_combinations)
  # Specifying the phi_0, i.e. the expected prediction without any features
  # Note: buffered column is a factor: 1 = "Buffered", 2 = "Scaling"
  # Lower (SHAP) values means prediction is closer to Buffered
  p <- mean(as.numeric(model$datasets$training$buffered))

  # Default method: ctree, as counted values are numerical but not continuous
  explanation <- explain(
    test_x_small,
    approach = method,
    explainer = explainer,
    prediction_zero = p
  )

  return(explanation)
}

shap2df <- function(explanation) {
  factor_values <- explanation$x_test %>%
    pivot_longer(everything(), names_to = "DosageCompensation.Factor", values_to = "DosageCompensation.Factor.Value")

  df_explanation <- explanation$dt %>%
    select(-none) %>%
    pivot_longer(everything(), names_to = "DosageCompensation.Factor", values_to = "SHAP.Value") %>%
    mutate(Factor.Value = factor_values$DosageCompensation.Factor.Value) %>%
    group_by(DosageCompensation.Factor) %>%
    mutate(Factor.Value.Relative = normalize_min_max(Factor.Value, na.rm = TRUE),
           SHAP.p25.Absolute = quantile(abs(SHAP.Value), probs = 0.25)[["25%"]],
           SHAP.Median.Absolute = median(abs(SHAP.Value)),
           SHAP.p75.Absolute = quantile(abs(SHAP.Value), probs = 0.75)[["75%"]]) %>%
    ungroup()

  df_corr <- df_explanation %>%
    group_by(DosageCompensation.Factor) %>%
    rstatix::cor_test(Factor.Value, SHAP.Value, method = "spearman", use = "na.or.complete") %>%
    select(DosageCompensation.Factor, cor, p) %>%
    rename(SHAP.Factor.Corr = cor, SHAP.Factor.Corr.p = p) %>%
    mutate(SHAP.Factor.Corr.p.adj = p.adjust(SHAP.Factor.Corr.p, method = "BY"))

  df_explanation <- df_explanation %>%
    left_join(y = df_corr, by = "DosageCompensation.Factor",
              relationship = "many-to-one", unmatched = "error", na_matches = "never")

  return(df_explanation)
}

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

    df_explanation <- readRDS(here(models_base_dir, model_filename)) %>%
      estimate_shap(n_samples = 250, n_combinations = 400) %>%  # Only override defaults if 128GB RAM available
      shap2df() %>%
      mutate(Model.Filename = model_filename,
             Model.BaseModel = model_name,
             Model.Variant = paste0(sub_dir, collapse = "_"))
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

### Write results to disk
shap_results %>%
  bind_rows() %>%
  write_parquet(here(output_data_dir, 'shap-analysis.parquet'), version = "2.6") %>%
  write.xlsx(here(tables_base_dir, "shap-analysis.xlsx"), colNames = TRUE)


# === Combine Plots for publishing ===
## SHAP Results
shap_results <- read_parquet(here(output_data_dir, 'shap-analysis.parquet'))

shap_gain <- shap_results %>%
  filter(Model.Variant == "ProCan_ChromosomeArm-Level_Gain_Log2FC")
shap_loss <- shap_results %>%
  filter(Model.Variant == "ProCan_ChromosomeArm-Level_Loss_Log2FC")

shap_arrows_plot_gain <- shap_plot_arrows(shap_gain, show_legend = FALSE, title = "Chromosome Gain")
shap_arrows_plot_loss <- shap_plot_arrows(shap_loss, category_lab = NULL, title = "Chromosome Loss")

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
    result[[group$Model.Variant]] <- evaluate_model(model, model$datasets$test, plots_dir)
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
    result[[group$Model.Variant]] <- evaluate_model(model, model$datasets$test, plots_dir)
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

# === Draft Area ===
# TODO: Cleanup, set as metadata in model file
# TODO: Control for WGD in ProCan
model_metadata <- data.frame(
  Model.Variant = c("DepMap_ChromosomeArm-Level_Gain_Log2FC", "DepMap_ChromosomeArm-Level_Loss_Log2FC",
                    "DepMap_Gene-Level_Filtered_Gain", "DepMap_Gene-Level_Filtered_Loss",
                    "DepMap-WGD_ChromosomeArm-Level_Loss", "DepMap-WGD_Gene-Level_Filtered_Loss",
                    "DepMap-NoWGD_ChromosomeArm-Level_Loss", "DepMap-NoWGD_Gene-Level_Filtered_Loss"),
  Model.Level = c("Chr", "Chr", "Gene", "Gene", "Chr", "Gene", "Chr", "Gene"),
  Model.Condition = c("Gain", "Loss", "Gain", "Loss", "Loss", "Loss", "Loss", "Loss"),
  Model.Control = c("All", "All", "All", "All", "WGD", "WGD", "Non-WGD", "Non-WGD")
)

model_metadata_wgd <- data.frame(
  Model.Variant = c("DepMap-WGD_ChromosomeArm-Level_Loss", "DepMap-WGD_Gene-Level_Filtered_Loss",
                    "DepMap-NoWGD_ChromosomeArm-Level_Loss", "DepMap-NoWGD_Gene-Level_Filtered_Loss"),
  Model.Level = c("Chr", "Gene", "Chr", "Gene"),
  Model.Condition = c("Loss", "Loss", "Loss", "Loss"),
  Model.Control = c("WGD", "WGD", "Non-WGD", "Non-WGD")
)

shap_results_selected <- shap_results %>%
  inner_join(y = model_metadata, by = "Model.Variant") %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, abs(SHAP.Factor.Corr), .desc = TRUE))

shap_results_selected_wgd <- shap_results %>%
  inner_join(y = model_metadata_wgd, by = "Model.Variant") %>%
  mutate(Label = map_signif(SHAP.Factor.Corr.p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, abs(SHAP.Factor.Corr), .desc = TRUE))

heatmap_shap_selected <- shap_results_selected %>%
  distinct(DosageCompensation.Factor, Model.Variant, Model.Control,
           Model.Level, Model.Condition, SHAP.Factor.Corr, Label) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = "", fill = SHAP.Factor.Corr, label = Label) +
  geom_raster() +
  geom_text(color = "black") +
  ggh4x::facet_nested(Model.Control + Model.Level + Model.Condition ~ ., switch = "y") +
  scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab",
                       limits = c(-1, 1), oob = scales::squish) +
  labs(x = "Feature", y = "Model", fill = "SHAP-Value-Feature-Correlation") +
  cowplot::theme_minimal_grid() +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color = "lightgrey"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(16, "points"),
        legend.key.width = unit(24, "points"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

heatmap_shap_selected %>%
  save_plot("shap_heatmap_selected.png", width = 300)

heatmap_shap_selected_wgd <- shap_results_selected_wgd %>%
  distinct(DosageCompensation.Factor, Model.Variant, Model.Control,
           Model.Level, Model.Condition, SHAP.Factor.Corr, Label) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = "", fill = SHAP.Factor.Corr, label = Label) +
  geom_raster() +
  geom_text(color = "black") +
  ggh4x::facet_nested(Model.Control + Model.Level + Model.Condition ~ ., switch = "y") +
  scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab",
                       limits = c(-1, 1), oob = scales::squish) +
  labs(x = "Feature", y = "Model", fill = "SHAP-Value-Feature-Correlation") +
  cowplot::theme_minimal_grid() +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color = "lightgrey"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(16, "points"),
        legend.key.width = unit(24, "points"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

heatmap_shap_selected_wgd %>%
  save_plot("shap_heatmap_selected_wgd.png", width = 300, height = 150)

# TODO: Visualize change between WGD and Non-WGD
## Difference SHAP, t-test, normalize SHAP
## Difference Correlation,
## Correlation of Correlation across methods
# TODO: Add Observation IDs to SHAP data (Gene Symbol, Cell Line ID, etc.)
# TODO: Random choice of SHAP observations could be confounded by Gene and Cell Line -> sample for each cell line

shap_cor_wgd_diff <- shap_results_selected_wgd %>%
  distinct(DosageCompensation.Factor, Model.Variant, Model.Control,
           Model.Level, Model.Condition, SHAP.Factor.Corr) %>%
  group_by(DosageCompensation.Factor, Model.Level, Model.Condition, Model.Control) %>%
  summarize(SHAP.Factor.Corr = first(SHAP.Factor.Corr), .groups = "drop") %>%
  pivot_wider(names_from = Model.Control, values_from = SHAP.Factor.Corr) %>%
  mutate(SHAP.Factor.Corr.Diff = `WGD` - `Non-WGD`)

shap_cor_wgd_diff %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = "", fill = SHAP.Factor.Corr.Diff) +
  geom_raster() +
  ggh4x::facet_nested(Model.Level + Model.Condition ~ ., switch = "y") +
  scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab",
                       limits = c(-1, 1), oob = scales::squish) +
  labs(x = "Feature", y = "Model", fill = "SHAP-Value-Feature-Correlation-Difference") +
  cowplot::theme_minimal_grid() +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color = "lightgrey"),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(16, "points"),
        legend.key.width = unit(24, "points"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
