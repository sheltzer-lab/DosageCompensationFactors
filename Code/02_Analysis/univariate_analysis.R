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

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "02_Analysis", "analysis.R"))


output_data_dir <- output_data_base_dir
plots_dir <- here(plots_base_dir, "Univariate")
goncalves_plots_dir <- here(plots_base_dir, "Univariate", "Goncalves")
goncalves_chr_plots_dir <- here(goncalves_plots_dir, "ChromosomeArm-Level")
goncalves_gene_plots_dir <- here(goncalves_plots_dir, "Gene-Level")
depmap_plots_dir <- here(plots_base_dir, "Univariate", "DepMap")
depmap_chr_plots_dir <- here(depmap_plots_dir, "ChromosomeArm-Level")
depmap_gene_plots_dir <- here(depmap_plots_dir, "Gene-Level")

dir.create(output_data_dir, recursive = TRUE)
dir.create(goncalves_plots_dir, recursive = TRUE)
dir.create(goncalves_chr_plots_dir, recursive = TRUE)
dir.create(goncalves_gene_plots_dir, recursive = TRUE)
dir.create(depmap_plots_dir, recursive = TRUE)
dir.create(depmap_chr_plots_dir, recursive = TRUE)
dir.create(depmap_gene_plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_goncalves <- read_parquet(here(output_data_dir, "expression_buffering_goncalves.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))

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
           ROC.AUC.Label = format(round(DosageCompensation.Factor.ROC.AUC, 3), nsmall = 3),
           DosageCompensation.Factor = factor(DosageCompensation.Factor,
                                              levels = DosageCompensation.Factor[order(DosageCompensation.Factor.ROC.AUC)])) %>%
    arrange(DosageCompensation.Factor.ROC.AUC)
}

plot_roc_auc_summary <- function(roc_auc_summary, plots_dir, filename) {
  roc_auc_summary_plot <- roc_auc_summary %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = DosageCompensation.Factor.ROC.AUC, label = ROC.AUC.Label) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.5) +
    geom_text(color = "white", nudge_y = -0.01) +
    scale_y_continuous(breaks = seq(0.45, 0.65, 0.05)) +
    xlab("") +
    ylab("ROC AUC") +
    coord_flip(ylim = c(0.45, 0.65)) +
    theme_light()

  dir.create(plots_dir)
  ggsave(here(plots_dir, filename), plot = roc_auc_summary_plot,
         height = 200, width = 180, units = "mm", dpi = 300)

  return(roc_auc_summary %>% select(-ROC.AUC.Label))
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

run_analysis <- function(dataset, buffering_class_col, filter_func, dir = NULL) {
  if (!is.null(dir)) {
    dir.create(dir, recursive = TRUE)

    roc_auc_summary <- dataset %>%
    filter_func() %>%
    reshape_factors({ { buffering_class_col } }) %>%
    determine_rocs() %>%
    plot_roc_curves(dir) %>%
    summarize_roc_auc() %>%
    plot_roc_auc_summary(dir, "buffering-factors_roc-auc.png")

    return(roc_auc_summary)
  } else {
    roc_auc_summary <- dataset %>%
    filter_func() %>%
    reshape_factors({ { buffering_class_col } }) %>%
    determine_rocs() %>%
    summarize_roc_auc()

    return(roc_auc_summary)
  }
}

# === Calculate ROC for all factors in all datasets ===
analysis_list <- list(
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = identity,
       dir = here(plots_dir, "Goncalves", "Gene-Level", "Unfiltered")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles,
       dir = here(plots_dir, "Goncalves", "Gene-Level", "Filtered")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain,
       dir = here(plots_dir, "Goncalves", "Gene-Level", "FilteredGain")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss,
       dir = here(plots_dir, "Goncalves", "Gene-Level", "FilteredLoss")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "Gain")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "Loss")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain_gene_avg,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "GainAverage")),
  list(dataset = expr_buf_goncalves, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss_gene_avg,
       dir = here(plots_dir, "Goncalves", "ChromosomeArm-Level", "LossAverage")),

  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = identity,
       dir = here(plots_dir, "DepMap", "Gene-Level", "Unfiltered")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_diff_quantiles,
       dir = here(plots_dir, "DepMap", "Gene-Level", "Filtered")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_gain,
       dir = here(plots_dir, "DepMap", "Gene-Level", "FilteredGain")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.GeneLevel.Class", filter = filter_cn_loss,
       dir = here(plots_dir, "DepMap", "Gene-Level", "FilteredLoss")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_gain,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "Gain")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Class", filter = filter_arm_loss,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "Loss")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_gain_gene_avg,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "GainAverage")),
  list(dataset = expr_buf_depmap, buffering = "Buffering.ChrArmLevel.Average.Class", filter = filter_arm_loss_gene_avg,
       dir = here(plots_dir, "DepMap", "ChromosomeArm-Level", "LossAverage"))
)

for (analysis in analysis_list) {
  run_analysis(dataset = analysis$dataset,
               buffering_class_col = get(analysis$buffering),
               filter_func = analysis$filter,
               dir = analysis$dir
  )
}

# === Statistically compare results ===

run_bootstrapped_analysis <- function(dataset, buffering_class_col, filter_func, n, sample_prop) {
  set.seed(42)
  dataset <- dataset %>%
    filter_func()
  results <- data.frame(DosageCompensation.Factor = character(),
                        DosageCompensation.Factor.ROC.AUC = numeric(),
                        Bootstrap.Sample = integer())
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n) {
    suppressMessages({
      results <- dataset %>%
        slice_sample(prop = sample_prop, replace = TRUE) %>%
        run_analysis(buffering_class_col = { { buffering_class_col } },
                     filter_func = identity) %>%
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

compare_conditions <- function(df_condition1, df_condition2) {
  # Merge dataframes for unified handling
  results_merged <- df_condition1 %>%
    bind_rows(df_condition2) %>%
    assertr::verify(length(unique(Condition)) == 2)

  conditions <- unique(results_merged$Condition)

  # Compare statistical significance between each factor
  results_factor_test <- results_merged %>%
    pivot_wider(id_cols = c(DosageCompensation.Factor, Bootstrap.Sample),
                names_from = Condition, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    group_by(DosageCompensation.Factor) %>%
    summarize(Wilcoxon.p.value = wilcox.test(get(conditions[1]), get(conditions[2]), paired = TRUE)$p.value) %>%
    mutate(Wilcoxon.significant = case_when(
      Wilcoxon.p.value < 0.0001 ~ "***",
      Wilcoxon.p.value < 0.001 ~ "**",
      Wilcoxon.p.value < 0.01 ~ "*",
      TRUE ~ "N.S."
    ))

  # Calculate summary statistics for each factor in each condition
  results_stat <- results_merged %>%
    pivot_wider(id_cols = c(DosageCompensation.Factor, Bootstrap.Sample),
                names_from = Condition, values_from = DosageCompensation.Factor.ROC.AUC) %>%
    group_by(DosageCompensation.Factor) %>%
    skimr::skim(conditions[1], conditions[2]) %>%
    rename(Condition = skim_variable) %>%
    ungroup()

  # Compare statistical significance between ranks of median values of factors
  results_rank_test <- results_stat %>%
    # Introduce pertubation to avoid equal ranks, otherwise p-value can't be calculated accurately
    mutate(RankValue = numeric.p50 + runif(length(numeric.p50), min = -1e-10, max = 1e-10)) %>%
    group_by(Condition) %>%
    mutate(DosageCompensation.Factor.Rank = as.integer(rank(-RankValue))) %>%
    select(Condition, DosageCompensation.Factor, DosageCompensation.Factor.Rank) %>%
    pivot_wider(id_cols = DosageCompensation.Factor,
                names_from = Condition, values_from = DosageCompensation.Factor.Rank)

  results_rank_test <- cor.test(results_rank_test[[conditions[1]]],
                                results_rank_test[[conditions[2]]],
                                method = "kendall")

  return(list(factor_test = results_factor_test,
              rank_test = results_rank_test,
              stat_summary = results_stat))
}

n <- 20
sample_prop <- 0.9

results_chr_gain <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.ChrArmLevel.Class,
                            filter_func = filter_arm_gain,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Chromosome Arm Gain")

results_chr_loss <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.ChrArmLevel.Class,
                            filter_func = filter_arm_loss,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Chromosome Arm Loss")

results_cn_gain <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.GeneLevel.Class,
                            filter_func = filter_cn_gain,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Gene Copy Number Gain")

results_cn_loss <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.GeneLevel.Class,
                            filter_func = filter_cn_loss,
                            n = n, sample_prop = sample_prop) %>%
  mutate(Condition = "Gene Copy Number Loss")

## Chr Gain vs. Chr Loss
results_chrgain_chrloss <- compare_conditions(results_chr_gain, results_chr_loss)
## CN gain vs. CN loss
results_cngain_cnloss <- compare_conditions(results_cn_gain, results_cn_loss)
## Chr gain vs. CN gain
results_chrgain_cngain <- compare_conditions(results_chr_gain, results_cn_gain)
## Chr loss vs. CN loss
results_chrloss_cnloss <- compare_conditions(results_chr_loss, results_cn_loss)

conditions <- unique(results_chrgain_chrloss$stat_summary$Condition)

plot1 <- results_chrgain_chrloss$stat_summary %>%
  assertr::verify(length(unique(Condition)) == 2) %>%
  filter(Condition == conditions[1]) %>%
  arrange(DosageCompensation.Factor) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = numeric.p50, label = format(round(numeric.p50, 3), nsmall = 3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.5) +
  geom_text(color = "white", nudge_y = -0.015) +
  scale_y_continuous(breaks = seq(0.45, 0.60, 0.05)) +
  ggtitle(conditions[1]) +
  xlab("") +
  ylab("Median ROC AUC") +
  coord_flip(ylim = c(0.45, 0.60)) +
  theme_minimal()  +
  theme(axis.text.y = element_blank())

plot2 <- results_chrgain_chrloss$stat_summary %>%
  assertr::verify(length(unique(Condition)) == 2) %>%
  filter(Condition == conditions[2]) %>%
  arrange(DosageCompensation.Factor) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = numeric.p50, label = format(round(numeric.p50, 3), nsmall = 3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.5) +
  geom_text(color = "white", nudge_y = -0.015) +
  scale_y_continuous(breaks = seq(0.45, 0.60, 0.05)) +
  ggtitle(conditions[2]) +
  xlab("") +
  ylab("Median ROC AUC") +
  coord_flip(ylim = c(0.45, 0.60)) +
  theme_minimal() +
  theme(axis.text.y = element_blank())

plot_factor_signif <- results_chrgain_chrloss$factor_test %>%
  arrange(DosageCompensation.Factor) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = Wilcoxon.p.value,
      label = Wilcoxon.significant) +
  geom_text(y = 1, color = "black") +
  xlab("") +
  ylab("") +
  coord_flip(ylim = c(0, 2)) +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot_labels <- results_chrgain_chrloss$factor_test %>%
  arrange(DosageCompensation.Factor) %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = Wilcoxon.p.value,
      label = DosageCompensation.Factor) +
  geom_text(y = 2, color = "black", hjust = 1) +
  xlab("") +
  ylab("") +
  coord_flip(ylim = c(0, 2)) +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ToDo: Use geom_bracket from ggpubr
plot_bracket <- broom::tidy(results_chrgain_chrloss$rank_test) %>%
  mutate(Label = paste0("p = ", format(round(p.value, 3), nsmall = 3),
                        ", τ = ", format(round(estimate, 3), nsmall = 3))) %>%
  ggplot() +
  aes(x = 0, y = 0, label = Label) +
  geom_segment(aes(x = 4, y = 1, xend = 8, yend = 1)) +
  geom_text(x = 6, color = "black", y = 2) +
  xlab("") +
  ylab("") +
  xlim(c(0, 10)) +
  ylim(c(0, 3)) +
  theme_void() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


plot_stack1 <- cowplot::plot_grid(plot_labels, plot1, plot_factor_signif, plot2,
                                  nrow = 1, ncol = 4, align = "h", axis = "l",
                                  rel_widths = c(0.8, 1, 0.1, 1))

plot_stack2 <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                  nrow = 2, ncol = 1,
                                  rel_heights = c(0.1, 1))

results_chr_gain %>%
  ggplot() +
  aes(x = DosageCompensation.Factor, y = DosageCompensation.Factor.ROC.AUC) +
  geom_violin(trim = FALSE)