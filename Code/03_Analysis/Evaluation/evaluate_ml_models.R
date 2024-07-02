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
dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# Load datasets
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_matched_renorm <- read_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'))
buf_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_wgd.parquet"))
buf_no_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_no-wgd.parquet"))
expr_buf_p0211 <- read_parquet(here(output_data_dir, 'expression_buffering_p0211.parquet'))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))

# Load Models
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
