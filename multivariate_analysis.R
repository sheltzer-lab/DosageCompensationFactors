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
library(MLeval)


here::i_am("DosageCompensationFactors.Rproj")

source(here("parameters.R"))
source(here("buffering_ratio.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Multivariate")

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

dc_factor_cols <- c(
  "Protein-Protein Interactions", "Protein Half-Life", "Protein Complexes (CORUM)",
  "Mean 3'-UTR Length", "Mean 5'-UTR Length",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites",
  "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate",
  "Translation Rate", "Protein Length", "mRNA Length",
  "Intrinsic Protein Disorder", "Low Complexity Score", "Homology Score",
  "Loops In Protein Score", "Protein Polyampholyte Score", "Protein Polarity",
  "Non-Exponential Decay Delta", "Mean mRNA Decay Rate", "Aggregation Score",
  ## Dataset-Specific Factors
  "Protein Neutral CV"
)

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

filter_cn_diff_quantiles <- function(df, remove_between = c("5%", "95%")) {
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01))
  df %>%
    filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_between[1]] |
             Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_between[2]])
}

filter_arm_gain <- function(df) {
  df %>% filter(ChromosomeArm.CNA > 0)
}

filter_arm_loss <- function(df) {
  df %>% filter(ChromosomeArm.CNA < 0)
}

explain_model <- function(model, dir) {
  # forest.model$finalModel
  # forest.model.eval <- evalm(forest.model)
  # forest.model.importance <- varImp(forest.model)
  png(here(dir, paste0("importance", ".png")),
        width = 200, height = 200, units = "mm", res = 200)
  varImpPlot(model$finalModel)
  dev.off()
}

evaluate_model <- function(model, test_set, dir) {
  test_predicted_prob <- predict(model, test_set, type = "prob")
  model_roc <- roc(response = test_set$buffered, predictor = as.numeric(test_predicted_prob[, "Buffered"]))
  png(here(dir, paste0("ROC-Curve", ".png")),
      width = 200, height = 200, units = "mm", res = 200)
  plot(model_roc, print.thres = "best", print.thres.best.method = "closest.topleft",
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.1)
  dev.off()
  # rf_coords <- coords(rf_roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))

  return(model_roc)
}

run_analysis <- function(dataset, buffering_class_col, filter_func, model_name, train_control, dir,
                         factor_cols = dc_factor_cols, training_set_ratio = 0.7, target_balance = 0.7) {
  set.seed(42)
  dir.create(dir, recursive = TRUE)

  # Filter and clean dataset before training
  df_prep <- dataset %>%
    filter_func() %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = factor({ { buffering_class_col } }, levels = c("Scaling", "Buffered"))) %>%
    drop_na(Buffered) %>%
    select(Buffered, all_of(factor_cols)) %>%
    impute_na() %>%
    rebalance_binary(Buffered, target_balance = target_balance) %>%
    shuffle_rows() %>%
    janitor::clean_names()


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

  explain_model(model, dir)
  model_roc <- evaluate_model(model, df_test, dir)

  return(list(model = model, test_set = df_test, train_set = df_train, roc = model_roc))
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
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = identity,
       dir = here(plots_dir, "Goncalves", "Gene-Level", "Unfiltered")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles,
       dir = here(plots_dir, "Goncalves", "Gene-Level", "Filtered")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "Gain")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "Loss")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "GainAverage")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "LossAverage")),

  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = identity,
       dir = here(plots_dir, "DepMap", "Gene-Level", "Unfiltered")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles,
       dir = here(plots_dir, "DepMap", "Gene-Level", "Filtered")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "Gain")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "Loss")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "GainAverage")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "LossAverage"))
)


## Random Forest
for (analysis in analysis_list) {
  run_analysis(dataset = analysis$dataset,
               buffering_class_col = get(analysis$buffering),
               filter_func = analysis$filter,
               model_name = "rf",
               train_control = tc_rf,
               dir = analysis$dir
  )
}

## Neural Network
for (analysis in analysis_list) {
  run_analysis(dataset = analysis$dataset,
               buffering_class_col = get(analysis$buffering),
               filter_func = analysis$filter,
               model_name = "pcaNNet",
               train_control = tc_nn,
               dir = analysis$dir
  )
}
