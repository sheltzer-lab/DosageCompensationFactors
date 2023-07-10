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
  results <- data.frame(DosageCompensation.Factor = factor(),
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
  close(pb)
  return(results)
}

n <- 100
sample_prop <- 0.9

results_chr_gain <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.ChrArmLevel.Average.Class,
                            filter_func = filter_arm_gain_gene_avg,
                            n = n, sample_prop = sample_prop)

results_chr_gain_ranked <- results_chr_gain %>%
  group_by(Bootstrap.Sample) %>%
  mutate(DosageCompensation.Factor.Rank = as.integer(rank(-DosageCompensation.Factor.ROC.AUC))) %>%
  ungroup() %>%
  arrange(DosageCompensation.Factor, Bootstrap.Sample) %>%
  select(Bootstrap.Sample, DosageCompensation.Factor.Rank) %>%
  split(~ Bootstrap.Sample) %>%
  sapply(\(df) list(df$DosageCompensation.Factor.Rank))

results_chr_loss <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.ChrArmLevel.Average.Class,
                            filter_func = filter_arm_loss_gene_avg,
                            n = n, sample_prop = sample_prop)

results_chr_loss_ranked <- results_chr_loss %>%
  group_by(Bootstrap.Sample) %>%
  mutate(DosageCompensation.Factor.Rank = as.integer(rank(-DosageCompensation.Factor.ROC.AUC))) %>%
  ungroup() %>%
  arrange(DosageCompensation.Factor, Bootstrap.Sample) %>%
  select(Bootstrap.Sample, DosageCompensation.Factor.Rank) %>%
  split(~ Bootstrap.Sample) %>%
  sapply(\(df) list(df$DosageCompensation.Factor.Rank))



results_cn_gain <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.GeneLevel.Class,
                            filter_func = filter_cn_gain,
                            n = n, sample_prop = sample_prop)

results_cn_gain_ranked <- results_cn_gain %>%
  group_by(Bootstrap.Sample) %>%
  mutate(DosageCompensation.Factor.Rank = as.integer(rank(-DosageCompensation.Factor.ROC.AUC))) %>%
  ungroup() %>%
  arrange(DosageCompensation.Factor, Bootstrap.Sample) %>%
  select(Bootstrap.Sample, DosageCompensation.Factor.Rank) %>%
  split(~ Bootstrap.Sample) %>%
  sapply(\(df) list(df$DosageCompensation.Factor.Rank))

results_cn_loss <- expr_buf_goncalves %>%
  run_bootstrapped_analysis(buffering_class_col = Buffering.GeneLevel.Class,
                            filter_func = filter_cn_loss,
                            n = n, sample_prop = sample_prop)

results_cn_loss_ranked <- results_cn_loss%>%
  group_by(Bootstrap.Sample) %>%
  mutate(DosageCompensation.Factor.Rank = as.integer(rank(-DosageCompensation.Factor.ROC.AUC))) %>%
  ungroup() %>%
  arrange(DosageCompensation.Factor, Bootstrap.Sample) %>%
  select(Bootstrap.Sample, DosageCompensation.Factor.Rank) %>%
  split(~ Bootstrap.Sample) %>%
  sapply(\(df) list(df$DosageCompensation.Factor.Rank))

# ToDo: Problematic, unlisting creates a single vector of all observations across all bootstrap samples
# Results in each bootstrap sample can be very different
cor.test(unlist(results_chr_gain_ranked),
         unlist(results_chr_loss_ranked),
         method = "kendall")

cor.test(unlist(results_cn_gain_ranked),
         unlist(results_cn_loss_ranked),
         method = "kendall")

cor.test(unlist(results_chr_gain_ranked),
         unlist(results_cn_gain_ranked),
         method = "kendall")

cor.test(unlist(results_chr_loss_ranked),
         unlist(results_cn_loss_ranked),
         method = "kendall")

## Averaged tests

results_chr_gain_avg <- results_chr_gain %>%
  group_by(DosageCompensation.Factor) %>%
  summarize(DosageCompensation.Factor.ROC.AUC = mean(DosageCompensation.Factor.ROC.AUC)) %>%
  mutate(DosageCompensation.Factor.Rank = as.integer(rank(-DosageCompensation.Factor.ROC.AUC))) %>%
  arrange(DosageCompensation.Factor)

results_cn_gain_avg <- results_cn_gain %>%
  group_by(DosageCompensation.Factor) %>%
  summarize(DosageCompensation.Factor.ROC.AUC = mean(DosageCompensation.Factor.ROC.AUC)) %>%
  mutate(DosageCompensation.Factor.Rank = as.integer(rank(-DosageCompensation.Factor.ROC.AUC))) %>%
  arrange(DosageCompensation.Factor)

cor.test(unlist(results_chr_gain_avg$DosageCompensation.Factor.Rank),
         unlist(results_cn_gain_avg$DosageCompensation.Factor.Rank),
         method = "kendall")