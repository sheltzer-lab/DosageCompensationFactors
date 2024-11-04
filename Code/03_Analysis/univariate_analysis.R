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
library(skimr)
library(cowplot)
library(broom)
library(corrplot)
library(openxlsx)
library(parallel)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "FactorAnalysis", "Univariate")
comparison_plots_dir <- here(plots_dir, "Comparison")
procan_plots_dir <- here(plots_dir, "ProCan")
procan_chr_plots_dir <- here(procan_plots_dir, "ChromosomeArm-Level")
procan_gene_plots_dir <- here(procan_plots_dir, "Gene-Level")
depmap_plots_dir <- here(plots_dir, "DepMap")
depmap_chr_plots_dir <- here(depmap_plots_dir, "ChromosomeArm-Level")
depmap_gene_plots_dir <- here(depmap_plots_dir, "Gene-Level")

dir.create(output_data_dir, recursive = TRUE)
dir.create(procan_plots_dir, recursive = TRUE)
dir.create(procan_chr_plots_dir, recursive = TRUE)
dir.create(procan_gene_plots_dir, recursive = TRUE)
dir.create(depmap_plots_dir, recursive = TRUE)
dir.create(depmap_chr_plots_dir, recursive = TRUE)
dir.create(depmap_gene_plots_dir, recursive = TRUE)

# === Load Datasets ===
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_matched_renorm <- read_parquet(here(output_data_dir, 'expression_buffering_matched_renorm.parquet'))
buf_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_wgd.parquet"))
buf_no_wgd <- read_parquet(here(output_data_dir, "expression_buffering_depmap_no-wgd.parquet"))
expr_buf_p0211 <- read_parquet(here(output_data_dir, 'expression_buffering_p0211.parquet'))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))

# === Define Processing Functions ===
establish_binary_classification <- function(df, buffering_class_col) {
  require(dplyr)

  df %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = ifelse({ { buffering_class_col } } == "Buffered", 1, 0)) %>%
    mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
    drop_na(Buffered)
}

reshape_factors <- function(df, factor_cols = dc_factor_cols, id_col = "UniqueId") {
  require(dplyr)
  require(tidyr)

  df %>%
    select(all_of(id_col), Buffered, all_of(factor_cols)) %>%
    pivot_longer(all_of(factor_cols),
                 names_to = "DosageCompensation.Factor",
                 values_to = "DosageCompensation.Factor.Value")
}

roc_silent <- function(response, predictor, factor = "") {
  require(pROC)

  if (length(response) == 0) return(NA)
  if (length(predictor) == 0) return(NA)
  if (all(is.na(response))) return(NA)
  if (all(is.na(predictor))) return(NA)

  tryCatch({
    return(roc(response, predictor, na.rm = TRUE, quiet = TRUE))
  }, error = function(e) {
    warning("An error occured when analyzing factor ", factor, "\nError message: ", e)
    return(NA)
  })
}

auc_na <- function (roc) {
  require(pROC)

  if (!is.list(roc)) return(NA)
  return(auc(roc, partial.auc.correct = TRUE))
}

determine_rocs <- function(df) {
  require(dplyr)

  df %>%
    group_by(DosageCompensation.Factor) %>%
    group_map(~list(factor = .y$DosageCompensation.Factor,
                    roc = roc_silent(.x$Buffered, .x$DosageCompensation.Factor.Value,
                                     factor = .y$DosageCompensation.Factor),
                    observations = sum(!is.na(.x$DosageCompensation.Factor.Value))),
              .keep = TRUE)
}

summarize_roc_auc <- function(factor_rocs) {
  require(dplyr)

  if (length(factor_rocs) == 0) return(data.frame(DosageCompensation.Factor = factor(),
                                                  DosageCompensation.Factor.Observations = integer(),
                                                  DosageCompensation.Factor.ROC.AUC = double()))

  data.frame(t(sapply(factor_rocs,
                      \(x) list(DosageCompensation.Factor = x$factor,
                                DosageCompensation.Factor.Observations = x$observations,
                                DosageCompensation.Factor.ROC.AUC = auc_na(x$roc)) %>% unlist()))) %>%
    mutate(DosageCompensation.Factor.ROC.AUC = as.numeric(DosageCompensation.Factor.ROC.AUC),
           DosageCompensation.Factor.Observations = as.integer(DosageCompensation.Factor.Observations),
           DosageCompensation.Factor = factor(DosageCompensation.Factor,
                                              levels = DosageCompensation.Factor[order(DosageCompensation.Factor.ROC.AUC)])) %>%
    arrange(DosageCompensation.Factor.ROC.AUC)
}

plot_roc_auc_summary <- function(roc_auc_summary, plots_dir, filename) {
  require(dplyr)

  roc_auc_summary %>%
    drop_na() %>%
    vertical_bar_chart(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC,
                       color_col = DosageCompensation.Factor.Observations,
                       value_lab = "ROC AUC", color_lab = "Observations") %>%
    save_plot(filename = filename, dir = plots_dir, height = 200, width = 180)

  return(roc_auc_summary)
}

roc_auc_summary_score <- function(df) {
  mean(abs(df$DosageCompensation.Factor.ROC.AUC - 0.5), na.rm = TRUE)
}

plot_roc_curves <- function(factor_rocs, dir) {
  require(pROC)

  dir <- here(dir, "ROC-Curves")
  dir.create(dir, recursive = TRUE)

  for (factor_roc in factor_rocs) {
    if (is.na(factor_roc$roc)) next

    png(here(dir, paste0(factor_roc$factor, ".png")),
        width = 200, height = 200, units = "mm", res = 200)
    plot(factor_roc$roc, main = factor_roc$factor,
         print.thres = "best", print.thres.best.method = "closest.topleft",
         print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.1)
    dev.off()
  }

  return(factor_rocs)
}

factor_roc_auc_prep_pipeline <- function(df, buffering_class_col, filter_func, df_factors, factor_cols) {
  require(dplyr)

  df %>%
    filter_func() %>%
    add_factors(df_factors) %>%
    establish_binary_classification({ { buffering_class_col } }) %>%
    select(UniqueId, Buffered, all_of(factor_cols))
}

factor_roc_auc_base_pipeline <- function(df) {
  require(dplyr)

  df %>%
    reshape_factors() %>%
    determine_rocs() %>%
    summarize_roc_auc()
}

run_analysis <- function(dataset, buffering_class_col, filter_func,
                         df_factors = dc_factors, factor_cols = dc_factor_cols) {
  require(dplyr)

  dataset %>%
    factor_roc_auc_prep_pipeline({ { buffering_class_col } }, filter_func, df_factors, factor_cols) %>%
    factor_roc_auc_base_pipeline()
}

# === Calculate ROC for all factors in all datasets ===
## Define datasets to train models on
datasets <- list(
  list(dataset = expr_buf_procan, name = "ProCan"),
  list(dataset = expr_buf_depmap, name = "DepMap"),
  list(dataset = expr_buf_matched_renorm, name = "MatchedRenorm"),
  list(dataset = buf_wgd, name = "DepMap-WGD"),
  list(dataset = buf_no_wgd, name = "DepMap-NoWGD"),
  list(dataset = expr_buf_p0211, name = "P0211"),
  list(dataset = expr_buf_cptac, name = "CPTAC")
)

## Define training data conditions
analysis_conditions <- list(
  list(buffering = "Buffering.GeneLevel.Class", filter = identity, sub_dir =  list("Gene-Level", "Unfiltered")),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff, sub_dir =  list("Gene-Level", "Filtered")),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain_abs, sub_dir =  list("Gene-Level", "Filtered_Gain")),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss_abs, sub_dir =  list("Gene-Level", "Filtered_Loss")),
  list(buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain, sub_dir =  list("ChromosomeArm-Level", "Gain")),
  list(buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss, sub_dir =  list("ChromosomeArm-Level", "Loss")),
  list(buffering = "Buffering.ChrArmLevel.Log2FC.Class", filter = filter_arm_gain, sub_dir = list("ChromosomeArm-Level", "Gain_Log2FC")),
  list(buffering = "Buffering.ChrArmLevel.Log2FC.Class", filter = filter_arm_loss, sub_dir = list("ChromosomeArm-Level", "Loss_Log2FC")),
  list(buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain_gene_avg, sub_dir = list("ChromosomeArm-Level", "Gain_Average")),
  list(buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss_gene_avg, sub_dir = list("ChromosomeArm-Level", "Loss_Average"))
)

## Run analysis
analysis_results <- data.frame(AnalysisID = character(),
                               DosageCompensation.Factor = character(),
                               DosageCompensation.Factor.ROC.AUC = numeric(),
                               DosageCompensation.Factor.Observations = integer())

for (dataset in datasets) {
  for (analysis in analysis_conditions) {
    target_dir <- here(plots_dir, append(dataset$name, analysis$sub_dir))
    analysis_id <- paste0(append(dataset$name, analysis$sub_dir), collapse = "_")

    suppressWarnings({ dir.create(target_dir, recursive = TRUE) })
    message("Univariate ROC AUC Analysis: ", analysis_id)

    analysis_results <- dataset$dataset %>%
      rename(BufferingClass = analysis$buffering) %>%
      run_analysis(buffering_class_col = BufferingClass,
                   filter_func = analysis$filter) %>%
      plot_roc_auc_summary(target_dir, "buffering-factors_roc-auc.png") %>%
      mutate(AnalysisID = analysis_id) %>%
      bind_rows(analysis_results)
  }
}

## Create aggregated ranking between methods
### Gain
rank_gain <- analysis_results %>%
  filter(grepl("Gain", AnalysisID) & !grepl("WGD", AnalysisID) & !grepl("P0211", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor)
rank_gain %>%
  vertical_bar_chart(DosageCompensation.Factor, MeanNormRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_gain.png", height = 200, width = 180)

### Loss
rank_loss <- analysis_results %>%
  filter(grepl("Loss", AnalysisID) & !grepl("WGD", AnalysisID) & !grepl("P0211", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor)
rank_loss %>%
  vertical_bar_chart(DosageCompensation.Factor, MeanNormRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_loss.png", height = 200, width = 180)

### WGD, Loss
rank_loss_wgd <- analysis_results %>%
  filter(grepl("Loss", AnalysisID) & grepl("DepMap-WGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor)
rank_loss_wgd %>%
  vertical_bar_chart(DosageCompensation.Factor, MeanNormRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_loss_wgd.png", height = 200, width = 180)

### No-WGD, Loss
rank_loss_no_wgd <- analysis_results %>%
  filter(grepl("Loss", AnalysisID) & grepl("DepMap-NoWGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor)
rank_loss_no_wgd %>%
  vertical_bar_chart(DosageCompensation.Factor, MeanNormRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_loss_no-wgd.png", height = 200, width = 180)

### WGD, Gain
rank_gain_wgd <- analysis_results %>%
  filter(grepl("Gain", AnalysisID) & grepl("DepMap-WGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor)
rank_gain_wgd %>%
  vertical_bar_chart(DosageCompensation.Factor, MeanNormRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_gain_wgd.png", height = 200, width = 180)

### No-WGD, Gain
rank_gain_no_wgd <- analysis_results %>%
  filter(grepl("Gain", AnalysisID) & grepl("DepMap-NoWGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor)
rank_gain_no_wgd %>%
  vertical_bar_chart(DosageCompensation.Factor, MeanNormRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_gain_no-wgd.png", height = 200, width = 180)

## Save results
analysis_results %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate.xlsx"),
             colNames = TRUE)
rank_gain %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate_aggregated_gain.xlsx"),
             colNames = TRUE)
rank_loss %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate_aggregated_loss.xlsx"),
             colNames = TRUE)
rank_loss_wgd %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss_WGD.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate_aggregated_loss_WGD.xlsx"),
             colNames = TRUE)
rank_loss_no_wgd %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_loss_NoWGD.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate_aggregated_loss_NoWGD.xlsx"),
             colNames = TRUE)
rank_gain_wgd %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain_WGD.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate_aggregated_gain_WGD.xlsx"),
             colNames = TRUE)
rank_gain_no_wgd %>%
  write_parquet(here(output_data_dir, "dosage_compensation_factors_univariate_aggregated_gain_NoWGD.parquet")) %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_factors_univariate_aggregated_gain_NoWGD.xlsx"),
             colNames = TRUE)

# === Statistically compare results ===
run_bootstrapped_analysis <- function(dataset, buffering_class_col, filter_func, n, sample_prop, cluster,
                                      df_factors = dc_factors, factor_cols = dc_factor_cols) {
  require(dplyr)

  duration <- system.time({
    results <- dataset %>%
      factor_roc_auc_prep_pipeline({ { buffering_class_col } }, filter_func, df_factors, factor_cols) %>%
      bootstrap_dataframe(cluster, n, sample_prop, factor_roc_auc_base_pipeline) %>%
      mutate(DosageCompensation.Factor = factor(DosageCompensation.Factor,
                                                levels = sort(unique(DosageCompensation.Factor))))
  })
  print(duration)

  return(results)
}

plot_comparison <- function(comparison_results) {
  require(dplyr)
  require(tidyr)
  require(assertr)
  require(ggplot2)
  require(cowplot)

  if(is.null(comparison_results)) return(NULL)

  conditions <- unique(comparison_results$stat_summary$Condition)

  plot1 <- comparison_results$stat_summary %>%
    assertr::verify(length(unique(Condition)) == 2) %>%
    filter(Condition == conditions[1]) %>%
    arrange(DosageCompensation.Factor) %>%
    vertical_bar_chart(DosageCompensation.Factor, numeric.p50,
                       error_low_col = numeric.p25, error_high_col = numeric.p75,
                       title = conditions[1], value_lab = "Median ROC AUC") +
    theme(axis.text.y = element_blank())

  plot2 <- comparison_results$stat_summary %>%
    filter(Condition == conditions[2]) %>%
    arrange(DosageCompensation.Factor) %>%
    vertical_bar_chart(DosageCompensation.Factor, numeric.p50,
                       error_low_col = numeric.p25, error_high_col = numeric.p75,
                       title = conditions[2], value_lab = "Median ROC AUC") +
    theme(axis.text.y = element_blank())

  plot_factor_signif <- comparison_results$factor_test %>%
    arrange(DosageCompensation.Factor) %>%
    plot_text_col(DosageCompensation.Factor, Wilcoxon.significant)

  plot_labels <- comparison_results$factor_test %>%
    arrange(DosageCompensation.Factor) %>%
    plot_text_col(DosageCompensation.Factor, DosageCompensation.Factor, align = "right")

  plot_bracket <- plot_corr_bracket(comparison_results$rank_test,
                                    estimate_symbol = utf8_tau, shift = 2)

  plot_stack1 <- cowplot::plot_grid(plot_labels, plot1, plot_factor_signif, plot2,
                                    nrow = 1, ncol = 4, align = "h", axis = "l",
                                    rel_widths = c(0.75, 1, 0.1, 1))
  plot_stack2 <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                    nrow = 2, ncol = 1,
                                    rel_heights = c(0.1, 1))

  return(plot_stack2)
}

datasets_bootstrap <- list(
  list(dataset = expr_buf_procan, name = "ProCan"),
  list(dataset = expr_buf_depmap, name = "DepMap"),
  list(dataset = expr_buf_cptac, name = "CPTAC")
)

conditions_bootstrap <- list(
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain, level = "Gene Copy Number", event = "Gain"),
  list(buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss, level = "Gene Copy Number", event = "Loss"),
  list(buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain, level = "Chromosome Arm", event = "Gain"),
  list(buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss, level = "Chromosome Arm", event = "Loss")
)

## Create cluster & run bootstrap analysis
cl <- makeCluster(getOption("cl.cores", max(1, detectCores() - 2)))
clusterExport(cl = cl, c("reshape_factors", "determine_rocs", "summarize_roc_auc", "dc_factor_cols",
                         "roc_silent", "auc_na", "factor_roc_auc_base_pipeline", "sample_func"),
              envir = environment())

bootstrap_results <- NULL
pb <- txtProgressBar(min = 0, max = length(datasets_bootstrap) * length(conditions_bootstrap), style = 3)
for (dataset in datasets_bootstrap) {
  for (condition in conditions_bootstrap) {
    result <- dataset$dataset %>%
      run_bootstrapped_analysis(buffering_class_col = get(condition$buffering),
                                filter_func = condition$filter, cluster = cl,
                                n = bootstrap_n, sample_prop = bootstrap_sample_prop) %>%
      mutate(Condition = paste(condition$level, condition$event),
             Level = condition$level,
             Event = condition$event,
             Dataset = dataset$name)

    bootstrap_results <- bind_rows(bootstrap_results, result)
    setTxtProgressBar(pb, pb$getVal() + 1)
  }
}
close(pb)
stopCluster(cl)

## Checkpoint: Save and load bootstrapped results before continuing
bootstrap_results %>%
  write_parquet(here(output_data_dir, 'bootstrap_univariate.parquet'),
                version = "2.6") %>%
  write.xlsx(here(tables_base_dir, "dosage_compensation_univariate_bootstrap.xlsx"),
             colNames = TRUE)

## Compare statistical results between conditions
bootstrap_results <- read_parquet(here(output_data_dir, 'bootstrap_univariate.parquet'))

for (dataset in datasets_bootstrap) {
  current_plots_dir <- here(comparison_plots_dir, dataset$name)
  dir.create(current_plots_dir, recursive = TRUE)

  bootstrap_chr_gain <- filter(bootstrap_results, Dataset == dataset$name & Level == "Chromosome Arm" & Event == "Gain")
  bootstrap_chr_loss <- filter(bootstrap_results, Dataset == dataset$name & Level == "Chromosome Arm" & Event == "Loss")
  bootstrap_cn_gain <- filter(bootstrap_results, Dataset == dataset$name & Level == "Gene Copy Number" & Event == "Gain")
  bootstrap_cn_loss <- filter(bootstrap_results, Dataset == dataset$name & Level == "Gene Copy Number" & Event == "Loss")

  ### Chr Gain vs. Chr Loss
  results_chrgain_chrloss <- compare_conditions(bootstrap_chr_gain, bootstrap_chr_loss)
  plot_chrgain_chrloss <- plot_comparison(results_chrgain_chrloss)
  ggsave(here(current_plots_dir, "roc-auc_comparison_chrgain_chrloss.png"),
         plot = plot_chrgain_chrloss,
         height = 200, width = 320, units = "mm", dpi = 300)

  ### CN gain vs. CN loss
  results_cngain_cnloss <- compare_conditions(bootstrap_cn_gain, bootstrap_cn_loss)
  plot_cngain_cnloss <- plot_comparison(results_cngain_cnloss)
  ggsave(here(current_plots_dir, "roc-auc_comparison_cngain_cnloss.png"),
         plot = plot_cngain_cnloss,
         height = 200, width = 320, units = "mm", dpi = 300)

  ### Chr gain vs. CN gain
  results_chrgain_cngain <- compare_conditions(bootstrap_chr_gain, bootstrap_cn_gain)
  plot_chrgain_cngain <- plot_comparison(results_chrgain_cngain)
  ggsave(here(current_plots_dir, "roc-auc_comparison_chrgain_cngain.png"),
         plot = plot_chrgain_cngain,
         height = 200, width = 320, units = "mm", dpi = 300)

  ### Chr loss vs. CN loss
  results_chrloss_cnloss <- compare_conditions(bootstrap_chr_loss, bootstrap_cn_loss)
  plot_chrloss_cnloss <- plot_comparison(results_chrloss_cnloss)
  ggsave(here(current_plots_dir, "roc-auc_comparison_chrloss_cnloss.png"),
         plot = plot_chrloss_cnloss,
         height = 200, width = 320, units = "mm", dpi = 300)

  ## Investigate correlation of ROC AUC of factors
  png(here(current_plots_dir, "corrplot_chrgain.png"), width = 300, height = 300, units = "mm", res = 200)
  corrplot_chr_gain <- bootstrap_chr_gain %>%
    select(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC, Condition, Bootstrap.Sample) %>%
    pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    select(-Condition, -Bootstrap.Sample) %>%
    plot_correlation()
  dev.off()

  png(here(current_plots_dir, "corrplot_chrloss.png"), width = 300, height = 300, units = "mm", res = 200)
  corrplot_chr_loss <- bootstrap_chr_loss %>%
    select(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC, Condition, Bootstrap.Sample) %>%
    pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    select(-Condition, -Bootstrap.Sample) %>%
    plot_correlation()
  dev.off()

  png(here(current_plots_dir, "corrplot_cngain.png"), width = 300, height = 300, units = "mm", res = 200)
  corrplot_cn_gain <- bootstrap_cn_gain %>%
    select(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC, Condition, Bootstrap.Sample) %>%
    pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    select(-Condition, -Bootstrap.Sample) %>%
    plot_correlation()
  dev.off()

  png(here(current_plots_dir, "corrplot_cnloss.png"), width = 300, height = 300, units = "mm", res = 200)
  corrplot_cn_loss <- bootstrap_cn_loss %>%
    select(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC, Condition, Bootstrap.Sample) %>%
    pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    select(-Condition, -Bootstrap.Sample) %>%
    plot_correlation()
  dev.off()

  ## Plot distribution of ROC AUC of factors
  dist_chr_gain <- bootstrap_chr_gain %>%
    violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
    save_plot("roc-auc_distribution_chrgain.png", current_plots_dir)
  dist_chr_loss <- bootstrap_chr_loss %>%
    violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
    save_plot("roc-auc_distribution_chrloss.png", current_plots_dir)
  dist_cn_gain <- bootstrap_cn_gain %>%
    violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
    save_plot("roc-auc_distribution_cngain.png", current_plots_dir)
  dist_cn_loss <- bootstrap_cn_loss %>%
    violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
    save_plot("roc-auc_distribution_cnloss.png", current_plots_dir)

  ## Plot Heatmap to compare conditions in a condensed way
  results <- list(results_chrgain_chrloss, results_cngain_cnloss, results_chrgain_cngain, results_chrloss_cnloss)
  if (any(sapply(results, is.null))) next

  bootstrap_auc <- bind_rows(bootstrap_chr_gain, bootstrap_chr_loss, bootstrap_cn_gain, bootstrap_cn_loss)

  rank_tests <- list(
    comparisons = list(c("Chromosome Arm Gain", "Chromosome Arm Loss"),
                       c("Gene Copy Number Gain", "Gene Copy Number Loss"),
                       c("Chromosome Arm Gain", "Gene Copy Number Gain"),
                       c("Chromosome Arm Loss", "Gene Copy Number Loss")),
    annotations = c(print_corr_obj(results_chrgain_chrloss$rank_test, estimate_symbol = utf8_tau, map_p = TRUE),
                    print_corr_obj(results_cngain_cnloss$rank_test, estimate_symbol = utf8_tau, map_p = TRUE),
                    print_corr_obj(results_chrgain_cngain$rank_test, estimate_symbol = utf8_tau, map_p = TRUE),
                    print_corr_obj(results_chrloss_cnloss$rank_test, estimate_symbol = utf8_tau, map_p = TRUE)),
    y_position = c(1, 1, 1.5, 2)
  )

  bootstrap_heatmap <- bootstrap_auc %>%
    roc_auc_heatmap(rank_tests) %>%
    save_plot("roc-auc_comparison_heatmap.png", dir = current_plots_dir, width = 400)
}

# === Correlation Analysis between Buffering Ratio and Factors ===
br_factor_cor <- bind_rows(expr_buf_procan, expr_buf_depmap, expr_buf_cptac) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  filter(Buffering.GeneLevel.Class != "Anti-Scaling") %>%
  mutate(Gene.CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  add_factors(dc_factors) %>%
  select(UniqueId, Gene.CNV, Dataset, Buffering.GeneLevel.Ratio, all_of(dc_factor_cols)) %>%
  pivot_longer(all_of(dc_factor_cols),
               names_to = "DosageCompensation.Factor",
               values_to = "DosageCompensation.Factor.Value") %>%
  group_by(Gene.CNV, Dataset, DosageCompensation.Factor) %>%
  rstatix::cor_test(Buffering.GeneLevel.Ratio, DosageCompensation.Factor.Value, method = "spearman") %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BY"))

br_factor_cor %>%
  write_parquet(here(output_data_dir, 'factor_correlation_univariate.parquet'),
                version = "2.6") %>%
  write.xlsx(here(tables_base_dir, "factor_correlation_univariate.xlsx"),
             colNames = TRUE)

br_factor_cor_heatmap <- br_factor_cor %>%
  mutate(Label = map_signif(p.adj),
         DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, cor, .desc = TRUE)) %>%
  distinct(DosageCompensation.Factor, Dataset, Gene.CNV, cor, Label) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = "", fill = cor, label = Label) +
  geom_raster() +
  geom_text(color = "black") +
  ggh4x::facet_nested(Dataset + Gene.CNV ~ ., switch = "y") +
  scale_fill_gradientn(colors = rev(bidirectional_color_pal), space = "Lab",
                       limits = c(-0.3, 0.3), oob = scales::squish) +
  labs(x = "Dosage Compensation Factor", y = "", fill = "Buffering Ratio Correlation") +
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
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 15, 5, 15), "mm"))

save_plot(br_factor_cor_heatmap, "buffering_ratio_factor_correlation.png", width = 280, height = 150)

max_abs_cor <- max(abs(br_factor_cor$cor))
