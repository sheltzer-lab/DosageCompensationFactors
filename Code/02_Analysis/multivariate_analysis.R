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


here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "02_Analysis", "analysis.R"))

models_base_dir <- here("Output", "Models")
output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Multivariate")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

shuffle_rows <- function(df) {
  df[sample(nrow(df), replace = FALSE),]
}

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

explain_model <- function(model, dir) {
  png(here(dir, paste0("importance", ".png")),
        width = 200, height = 200, units = "mm", res = 200)
  varImpPlot(model$finalModel)
  dev.off()
}

evaluate_model <- function(model, test_set, dir) {
  test_predicted_prob <- predict(model, test_set, type = "prob")
  model_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob[, "Buffered"]), na.rm = TRUE)
  png(here(dir, paste0("ROC-Curve", ".png")),
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
    drop_na(Buffered) %>%
    select(Buffered, all_of(factor_cols)) %>%
    janitor::clean_names()
}

run_analysis <- function(dataset, buffering_class_col, filter_func, model_name, train_control,
                         plots_dir, models_dir, sub_dir,
                         factor_cols = dc_factor_cols, training_set_ratio = 0.7, target_balance = 0.7) {
  set.seed(42)
  plots_dir <- here(plots_dir, sub_dir)
  dir.create(plots_dir, recursive = TRUE)
  dir.create(models_dir, recursive = TRUE)

  # Filter and clean dataset before training
  df_prep <- dataset %>%
    filter_func() %>%
    clean_data({{ buffering_class_col }}, factor_cols = factor_cols) %>%
    impute_na() %>%
    rebalance_binary(buffered, target_balance = target_balance) %>%
    shuffle_rows()

  # Split training & test data
  train_size <- ceiling(nrow(df_prep) * training_set_ratio)
  test_size <- nrow(df_prep) - train_size

  df_train <- df_prep[1:train_size,]
  df_test <- df_prep[(train_size + 1):(train_size + test_size),]

  # Train model
  model <- caret::train(buffered ~ .,
                        data = df_train,
                        method = model_name,
                        trControl = train_control,
                        metric = "ROC")

  # Save Model
  saveRDS(model, here(models_dir, paste0("model_", model_name, "_", paste0(sub_dir, collapse = ""), ".rds")))

  # Evaluate Model
  explain_model(model, plots_dir)
  model_roc <- evaluate_model(model, df_test, plots_dir)

  return(list(model = model, roc = model_roc))
}

# === Build Models ===
tc_rf <- trainControl(method = "cv",
                      number = 2,
                      savePredictions = TRUE,
                      classProbs = TRUE,
                      verboseIter = TRUE,
                      summaryFunction = twoClassSummary)

tc_nn <- trainControl(method = "cv",
                      number = 3,
                      savePredictions = TRUE,
                      classProbs = TRUE,
                      verboseIter = TRUE,
                      summaryFunction = twoClassSummary)

analysis_list <- list(
  # list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = identity,
  #      sub_dir =  list("Goncalves", "Gene-Level", "Unfiltered")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles,
       sub_dir =  list("Goncalves", "Gene-Level", "Filtered")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain,
       sub_dir =  list("Goncalves", "Gene-Level", "FilteredGain")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss,
       sub_dir =  list("Goncalves", "Gene-Level", "FilteredLoss")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain,
       sub_dir =  list("Goncalves", "ChromosomeArm-Level", "Gain")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss,
       sub_dir =  list("Goncalves", "ChromosomeArm-Level", "Loss")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain_gene_avg,
       sub_dir =  list("Goncalves", "ChromosomeArm-Level", "GainAverage")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss_gene_avg,
       sub_dir =  list("Goncalves", "ChromosomeArm-Level", "LossAverage")),

  # list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = identity,
  #      sub_dir =  list("DepMap", "Gene-Level", "Unfiltered")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles,
       sub_dir =  list("DepMap", "Gene-Level", "Filtered")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain,
       sub_dir =  list("DepMap", "Gene-Level", "FilteredGain")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss,
       sub_dir =  list("DepMap", "Gene-Level", "FilteredLoss")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain,
       sub_dir =  list("DepMap", "ChromosomeArm-Level", "Gain")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss,
       sub_dir =  list("DepMap", "ChromosomeArm-Level", "Loss")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain_gene_avg,
       sub_dir =  list("DepMap", "ChromosomeArm-Level", "GainAverage")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss_gene_avg,
       sub_dir =  list("DepMap", "ChromosomeArm-Level", "LossAverage"))
)


## Random Forest
analysis_results_rf <- list()
for (analysis in analysis_list) {
  results <- run_analysis(dataset = analysis$dataset,
                          buffering_class_col = get(analysis$buffering),
                          filter_func = analysis$filter,
                          model_name = "rf",
                          train_control = tc_rf,
                          plots_dir = plots_dir,
                          sub_dir = analysis$sub_dir,
                          models_dir = models_base_dir
  )
  result_id <- paste0(analysis$sub_dir, collapse = "_")
  analysis_results_rf[[result_id]] <- results
}

### Use goncalves gain model for evaluating depmap gene-level gain data
model_goncalves_filtered <- readRDS(here(models_base_dir, "model_rf_GoncalvesGene-LevelFiltered.rds"))

test_data <- expr_buf_depmap %>%
  filter_cn_diff_quantiles() %>%
  clean_data(Buffering.GeneLevel.Class) %>%
  impute_na() %>%
  shuffle_rows()

eval_dir <- here(plots_dir, "Goncalves", "Gene-Level", "Filtered", "Eval-DepMap")
dir.create(eval_dir, recursive = TRUE)
evaluate_model(model_goncalves_filtered, test_data, eval_dir)


## Neural Network
for (analysis in analysis_list) {
  run_analysis(dataset = analysis$dataset,
               buffering_class_col = get(analysis$buffering),
               filter_func = analysis$filter,
               model_name = "pcaNNet",
               train_control = tc_nn,
               plots_dir = plots_dir,
               sub_dir = analysis$sub_dir,
               models_dir = models_base_dir
  )
}