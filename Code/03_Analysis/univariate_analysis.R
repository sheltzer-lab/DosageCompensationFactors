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

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "03_Analysis", "analysis.R"))

output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Univariate")
procan_plots_dir <- here(plots_base_dir, "Univariate", "ProCan")
procan_chr_plots_dir <- here(procan_plots_dir, "ChromosomeArm-Level")
procan_gene_plots_dir <- here(procan_plots_dir, "Gene-Level")
procan_comparison_plots_dir <- here(procan_plots_dir, "Comparison")
depmap_plots_dir <- here(plots_base_dir, "Univariate", "DepMap")
depmap_chr_plots_dir <- here(depmap_plots_dir, "ChromosomeArm-Level")
depmap_gene_plots_dir <- here(depmap_plots_dir, "Gene-Level")

dir.create(output_data_dir, recursive = TRUE)
dir.create(procan_plots_dir, recursive = TRUE)
dir.create(procan_chr_plots_dir, recursive = TRUE)
dir.create(procan_gene_plots_dir, recursive = TRUE)
dir.create(procan_comparison_plots_dir, recursive = TRUE)
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

# === Define Processing Functions ===
reshape_factors <- function(df, buffering_class_col, factor_cols = dc_factor_cols, id_col = "UniqueId") {
  df %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = ifelse({ { buffering_class_col } } == "Buffered", 1, 0)) %>%
    mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
    drop_na(Buffered) %>%
    select(all_of(id_col), Buffered, all_of(factor_cols)) %>%
    pivot_longer(all_of(factor_cols),
                 names_to = "DosageCompensation.Factor",
                 values_to = "DosageCompensation.Factor.Value")
}

determine_rocs <- function(df) {
  df %>%
    group_by(DosageCompensation.Factor) %>%
    group_map(~list(factor = .y$DosageCompensation.Factor,
                    roc = roc(.x$Buffered, .x$DosageCompensation.Factor.Value, na.rm = TRUE)), .keep = TRUE)
}

summarize_roc_auc <- function(factor_rocs) {
  data.frame(t(sapply(factor_rocs,
                                         \(x) list(DosageCompensation.Factor = x$factor,
                                                   DosageCompensation.Factor.ROC.AUC = auc(x$roc)) %>% unlist()))) %>%
    mutate(DosageCompensation.Factor.ROC.AUC = as.numeric(DosageCompensation.Factor.ROC.AUC),
           DosageCompensation.Factor = factor(DosageCompensation.Factor,
                                              levels = DosageCompensation.Factor[order(DosageCompensation.Factor.ROC.AUC)])) %>%
    arrange(DosageCompensation.Factor.ROC.AUC)
}

plot_roc_auc_summary <- function(roc_auc_summary, plots_dir, filename) {
  roc_auc_summary %>%
    vertical_bar_chart(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC,
                       value_lab = "ROC AUC") %>%
    save_plot(filename = filename, dir = plots_dir, height = 200, width = 180)

  return(roc_auc_summary)
}

roc_auc_summary_score <- function(df) {
  mean(abs(df$DosageCompensation.Factor.ROC.AUC - 0.5))
}

plot_roc_curves <- function(factor_rocs, dir) {
  dir <- here(dir, "ROC-Curves")
  dir.create(dir, recursive = TRUE)

  for (factor_roc in factor_rocs) {
    png(here(dir, paste0(factor_roc$factor, ".png")),
        width = 200, height = 200, units = "mm", res = 200)
    plot(factor_roc$roc, main = factor_roc$factor,
         print.thres = "best", print.thres.best.method = "closest.topleft",
         print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.1)
    dev.off()
  }

  return(factor_rocs)
}

run_analysis <- function(dataset, buffering_class_col, filter_func, df_factors = dc_factors) {
  dataset %>%
      filter_func() %>%
      add_factors(df_factors) %>%
      reshape_factors({ { buffering_class_col } }) %>%
      determine_rocs() %>%
      summarize_roc_auc()
}

# === Calculate ROC for all factors in all datasets ===
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
  list(buffering = "Buffering.GeneLevel.Class", filter = identity, sub_dir =  list("Gene-Level", "Unfiltered")),
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

## Run analysis
analysis_results <- data.frame(AnalysisID = character(),
                               DosageCompensation.Factor = character(),
                               DosageCompensation.Factor.ROC.AUC = numeric())
for (dataset in datasets) {
  for (analysis in analysis_conditions) {
    target_dir <- here(plots_dir, append(dataset$name, analysis$sub_dir))
    dir.create(target_dir, recursive = TRUE)

    analysis_results <- dataset$dataset %>%
      run_analysis(buffering_class_col = get(analysis$buffering),
                   filter_func = analysis$filter) %>%
      plot_roc_auc_summary(target_dir, "buffering-factors_roc-auc.png") %>%
      mutate(AnalysisID = paste0(append(dataset$name, analysis$sub_dir), collapse = "_")) %>%
      bind_rows(analysis_results)
  }
}

## Create aggregated ranking between methods
rank_gain <- analysis_results %>%
  filter(grepl("Gain", AnalysisID) & !grepl("WGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor) %>%
  vertical_bar_chart(DosageCompensation.Factor, AggregatedRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_gain.png", height = 200, width = 180)


rank_loss <- analysis_results %>%
  filter(grepl("Loss", AnalysisID) & !grepl("WGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor) %>%
  vertical_bar_chart(DosageCompensation.Factor, AggregatedRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_loss.png", height = 200, width = 180)

rank_loss_wgd <- analysis_results %>%
  filter(grepl("Loss", AnalysisID) & grepl("DepMap-WGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor) %>%
  vertical_bar_chart(DosageCompensation.Factor, AggregatedRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_loss_wgd.png", height = 200, width = 180)

rank_loss_no_wgd <- analysis_results %>%
  filter(grepl("Loss", AnalysisID) & grepl("DepMap-NoWGD", AnalysisID)) %>%
  mean_norm_rank(DosageCompensation.Factor.ROC.AUC,
                 AnalysisID, DosageCompensation.Factor) %>%
  vertical_bar_chart(DosageCompensation.Factor, AggregatedRank,
                     value_range = c(0, 1), break_steps = 0.1, value_lab = "Aggregated Rank",
                     bar_label_shift = 0.07, line_intercept = 0) %>%
  save_plot("buffering-factors_rank_loss_no-wgd.png", height = 200, width = 180)

# === Statistically compare results ===

run_bootstrapped_analysis <- function(dataset, buffering_class_col, filter_func,
                                      n, sample_prop, df_factors = dc_factors) {
  set.seed(42)

  dataset <- dataset %>%
    filter_func() %>%
    add_factors(df_factors)

  results <- data.frame(DosageCompensation.Factor = character(),
                        DosageCompensation.Factor.ROC.AUC = numeric(),
                        Bootstrap.Sample = integer())
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n) {
    suppressMessages({
      results <- dataset %>%
        slice_sample(prop = sample_prop, replace = TRUE) %>%
        reshape_factors({ { buffering_class_col } }) %>%
        determine_rocs() %>%
        summarize_roc_auc() %>%
        mutate(Bootstrap.Sample = i) %>%
        bind_rows(results)
    })
    setTxtProgressBar(pb, i)
  }
  results <- results %>%
    mutate(DosageCompensation.Factor = factor(DosageCompensation.Factor,
                                              levels = sort(unique(DosageCompensation.Factor))))
  close(pb)
  return(results)
}

# Notes for Statistical Tests and Multiple Testing Correction
# * Bootstrapping values independently samples values from the same distribution (with replacement) => iid.
# ** Proof: Show that sampling is iid. (linear map from set of indices to set of objects)
# * Check : ROC AUC of factors may be correlated -> Not independent!

compare_conditions <- function(df_condition1, df_condition2) {
  # Merge dataframes for unified handling
  df_merged <- df_condition1 %>%
    bind_rows(df_condition2) %>%
    assertr::verify(length(unique(Condition)) == 2)

  conditions <- unique(df_merged$Condition)

  # Compare statistical significance between each factor
  df_factor_test <- df_merged %>%
    pivot_wider(id_cols = c(DosageCompensation.Factor, Bootstrap.Sample),
                names_from = Condition, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    group_by(DosageCompensation.Factor) %>%
    summarize(Wilcoxon.p.value = wilcox.test(get(conditions[1]), get(conditions[2]), paired = TRUE)$p.value) %>%
    mutate(Wilcoxon.p.adjusted = p.adjust(Wilcoxon.p.value, method = "BY")) %>%
    mutate(Wilcoxon.significant = map_signif(Wilcoxon.p.adjusted))

  # Calculate summary statistics for each factor in each condition
  df_stat <- df_merged %>%
    pivot_wider(id_cols = c(DosageCompensation.Factor, Bootstrap.Sample),
                names_from = Condition, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    group_by(DosageCompensation.Factor) %>%
    skimr::skim(conditions[1], conditions[2]) %>%
    rename(Condition = skim_variable) %>%
    ungroup()

  # Compare statistical significance between ranks of median values of factors
  df_rank_test <- df_stat %>%
    # Introduce pertubation to avoid equal ranks, otherwise p-value can't be calculated accurately
    mutate(RankValue = numeric.p50 + runif(length(numeric.p50), min = -1e-10, max = 1e-10)) %>%
    group_by(Condition) %>%
    mutate(DosageCompensation.Factor.Rank = as.integer(rank(-RankValue))) %>%
    select(Condition, DosageCompensation.Factor, DosageCompensation.Factor.Rank) %>%
    pivot_wider(id_cols = DosageCompensation.Factor,
                names_from = Condition, values_from = DosageCompensation.Factor.Rank)

  df_rank_test <- cor.test(df_rank_test[[conditions[1]]],
                           df_rank_test[[conditions[2]]],
                           method = "kendall")

  return(list(factor_test = df_factor_test,
              rank_test = df_rank_test,
              stat_summary = df_stat))
}

plot_comparison <- function(comparison_results) {
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

  plot_bracket <- broom::tidy(comparison_results$rank_test) %>%
    mutate(Label = paste0(print_signif(p.value),
                          ", τ = ", format(round(estimate, 3), nsmall = 3))) %>%
    ggplot() +
    aes(x = 0, y = 0, label = Label) +
    geom_segment(aes(x = 4, y = 1, xend = 8, yend = 1)) +
    geom_text(x = 6, color = "black", y = 2) +
    xlab(NULL) +
    ylab(NULL) +
    xlim(c(0, 10)) +
    ylim(c(0, 3)) +
    theme_void() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())


  plot_stack1 <- cowplot::plot_grid(plot_labels, plot1, plot_factor_signif, plot2,
                                    nrow = 1, ncol = 4, align = "h", axis = "l",
                                    rel_widths = c(0.75, 1, 0.1, 1))
  plot_stack2 <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                    nrow = 2, ncol = 1,
                                    rel_heights = c(0.1, 1))

  return(plot_stack2)
}

n <- 1000
sample_prop <- 0.9

bootstrap_chr_gain <- expr_buf_procan %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.ChrArmLevel.Log2FC.Class,
                            filter_func = filter_arm_gain,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Chromosome Arm Gain")

bootstrap_chr_loss <- expr_buf_procan %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.ChrArmLevel.Log2FC.Class,
                            filter_func = filter_arm_loss,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Chromosome Arm Loss")

bootstrap_cn_gain <- expr_buf_procan %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.GeneLevel.Class,
                            filter_func = filter_cn_gain,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Gene Copy Number Gain")

bootstrap_cn_loss <- expr_buf_procan %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.GeneLevel.Class,
                            filter_func = filter_cn_loss,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Gene Copy Number Loss")

## Checkpoint: Save and load bootstrapped results before continuing
write_parquet(bootstrap_chr_gain, here(output_data_dir, 'bootstrap_univariate_procan_chrgain.parquet'),
              version = "2.6")
write_parquet(bootstrap_chr_loss, here(output_data_dir, 'bootstrap_univariate_procan_chrloss.parquet'),
              version = "2.6")
write_parquet(bootstrap_cn_gain, here(output_data_dir, 'bootstrap_univariate_procan_cngain.parquet'),
              version = "2.6")
write_parquet(bootstrap_cn_loss, here(output_data_dir, 'bootstrap_univariate_procan_cnloss.parquet'),
              version = "2.6")
bootstrap_chr_gain <- read_parquet(here(output_data_dir, "bootstrap_univariate_procan_chrgain.parquet"))
bootstrap_chr_loss <- read_parquet(here(output_data_dir, "bootstrap_univariate_procan_chrloss.parquet"))
bootstrap_cn_gain <- read_parquet(here(output_data_dir, "bootstrap_univariate_procan_cngain.parquet"))
bootstrap_cn_loss <- read_parquet(here(output_data_dir, "bootstrap_univariate_procan_cnloss.parquet"))

## Compare statistical results between conditions
### Chr Gain vs. Chr Loss
results_chrgain_chrloss <- compare_conditions(bootstrap_chr_gain, bootstrap_chr_loss)
plot_chrgain_chrloss <- plot_comparison(results_chrgain_chrloss)
ggsave(here(procan_comparison_plots_dir, "roc-auc_comparison_chrgain_chrloss.png"),
       plot = plot_chrgain_chrloss,
       height = 200, width = 320, units = "mm", dpi = 300)

### CN gain vs. CN loss
results_cngain_cnloss <- compare_conditions(bootstrap_cn_gain, bootstrap_cn_loss)
plot_cngain_cnloss <- plot_comparison(results_cngain_cnloss)
ggsave(here(procan_comparison_plots_dir, "roc-auc_comparison_cngain_cnloss.png"),
       plot = plot_cngain_cnloss,
       height = 200, width = 320, units = "mm", dpi = 300)

### Chr gain vs. CN gain
results_chrgain_cngain <- compare_conditions(bootstrap_chr_gain, bootstrap_cn_gain)
plot_chrgain_cngain <- plot_comparison(results_chrgain_cngain)
ggsave(here(procan_comparison_plots_dir, "roc-auc_comparison_chrgain_cngain.png"),
       plot = plot_chrgain_cngain,
       height = 200, width = 320, units = "mm", dpi = 300)

### Chr loss vs. CN loss
results_chrloss_cnloss <- compare_conditions(bootstrap_chr_loss, bootstrap_cn_loss)
plot_chrloss_cnloss <- plot_comparison(results_chrloss_cnloss)
ggsave(here(procan_comparison_plots_dir, "roc-auc_comparison_chrloss_cnloss.png"),
       plot = plot_chrloss_cnloss,
       height = 200, width = 320, units = "mm", dpi = 300)

## Investigate correlation of ROC AUC of factors
png(here(procan_comparison_plots_dir, "corrplot_chrgain.png"), width = 300, height = 300, units = "mm", res = 200)
corrplot_chr_gain <- bootstrap_chr_gain %>%
  pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
  select(-Condition, -Bootstrap.Sample) %>%
  plot_correlation()
dev.off()

png(here(procan_comparison_plots_dir, "corrplot_chrloss.png"), width = 300, height = 300, units = "mm", res = 200)
corrplot_chr_loss <- bootstrap_chr_loss %>%
  pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
  select(-Condition, -Bootstrap.Sample) %>%
  plot_correlation()
dev.off()

png(here(procan_comparison_plots_dir, "corrplot_cngain.png"), width = 300, height = 300, units = "mm", res = 200)
corrplot_cn_gain <- bootstrap_cn_gain %>%
  pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
  select(-Condition, -Bootstrap.Sample) %>%
  plot_correlation()
dev.off()

png(here(procan_comparison_plots_dir, "corrplot_cnloss.png"), width = 300, height = 300, units = "mm", res = 200)
corrplot_cn_loss <- bootstrap_cn_loss %>%
  pivot_wider(names_from = DosageCompensation.Factor, values_from = DosageCompensation.Factor.ROC.AUC) %>%
  select(-Condition, -Bootstrap.Sample) %>%
  plot_correlation()
dev.off()

## Plot distribution of ROC AUC of factors
dist_chr_gain <- bootstrap_chr_gain %>%
  violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
  save_plot("roc-auc_distribution_chrgain.png", procan_comparison_plots_dir)
dist_chr_loss <- bootstrap_chr_loss %>%
  violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
  save_plot("roc-auc_distribution_chrloss.png", procan_comparison_plots_dir)
dist_cn_gain <- bootstrap_cn_gain %>%
  violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
  save_plot("roc-auc_distribution_cngain.png", procan_comparison_plots_dir)
dist_cn_loss <- bootstrap_cn_loss %>%
  violin_plot(DosageCompensation.Factor, DosageCompensation.Factor.ROC.AUC) %>%
  save_plot("roc-auc_distribution_cnloss.png", procan_comparison_plots_dir)
