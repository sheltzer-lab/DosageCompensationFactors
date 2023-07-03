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

source(here("parameters.R"))
source(here("buffering_ratio.R"))

output_data_dir <- output_data_base_dir
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


# === Calculate ROC AUC of factors ===
dc_factor_cols <- c(
  "Protein-Protein Interactions", "Protein Half-Life", "Protein Complexes (CORUM)",
  "Mean 3'-UTR Length", "Mean 5'-UTR Length",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites",
  "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate",
  "Translation Rate", "Protein Length", "mRNA Length",
  "Intrinsic Protein Disorder", "Low Complexity Score", "Homology Score",
  "Loops In Protein Score", "Protein Polyampholyte Score", "Protein Polarity",
  "Non-Exponential Decay Delta", "Mean mRNA Decay Rate"
)

analyze_roc_auc <- function(df, buffering_class_col, factor_cols = dc_factor_cols, id_col = "UniqueId") {
  df %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = ifelse({ { buffering_class_col } } == "Buffered", 1, 0)) %>%
    mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
    drop_na(Buffered) %>%
    select(all_of(id_col), Buffered, all_of(factor_cols)) %>%
    pivot_longer(all_of(factor_cols),
                 names_to = "DosageCompensation.Factor",
                 values_to = "DosageCompensation.Factor.Value") %>%
    group_by(DosageCompensation.Factor) %>%
    summarize(DosageCompensation.Factor.ROC.AUC =
                list(auc(Buffered, DosageCompensation.Factor.Value, na.rm = TRUE)) %>% unlist()) %>%
    mutate(DosageCompensation.Factor =
             factor(DosageCompensation.Factor,
                    levels = DosageCompensation.Factor[order(DosageCompensation.Factor.ROC.AUC)])) %>%
    arrange(DosageCompensation.Factor.ROC.AUC)
}

plot_roc_auc_summary <- function(factors_roc_auc, plots_dir, filename) {
  roc_auc_summary_plot <- factors_roc_auc %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = DosageCompensation.Factor.ROC.AUC,
        label = format(round(DosageCompensation.Factor.ROC.AUC, 3), nsmall = 3)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.5) +
    geom_text(color = "white", nudge_y = -0.007) +
    scale_y_continuous(breaks = seq(0.45, 0.60, 0.05)) +
    xlab("") +
    ylab("ROC AUC") +
    coord_flip(ylim = c(0.45, 0.60)) +
    theme_light()

  ggsave(here(plots_dir, filename), plot = roc_auc_summary_plot,
         height = 200, width = 180, units = "mm", dpi = 300)
}

roc_auc_summary_score <- function(df) {
  mean(abs(df$DosageCompensation.Factor.ROC.AUC - 0.5))
}

plot_roc_curves <- function(df, buffering_class_col, dir, factor_cols = dc_factor_cols) {
  dir <- here(dir, "ROC-Curves")
  dir.create(dir, recursive = TRUE)

  plot_roc_curve <- function(df_dc_factors, factor) {
    df <- df_dc_factors %>%
      filter(DosageCompensation.Factor == factor)
    roc <- roc(df$Buffered, df$DosageCompensation.Factor.Value, na.rm = TRUE)
    roc_plot <- plot(roc, main = factor, print.thres = "best", print.thres.best.method = "closest.topleft", print.auc = TRUE)
  }

  df_dc_factors <- df %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = ifelse({ { buffering_class_col } } == "Buffered", 1, 0)) %>%
    mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
    drop_na(Buffered) %>%
    select(Buffered, all_of(factor_cols)) %>%
    pivot_longer(all_of(factor_cols),
                 names_to = "DosageCompensation.Factor",
                 values_to = "DosageCompensation.Factor.Value")

  for (factor in factor_cols) {
    png(here(dir, paste0(factor, ".png")),
        width = 200, height = 200, units = "mm", res = 200)
    plot_roc_curve(df_dc_factors, factor)
    dev.off()
  }
}

# === Calculate ROC AUCs for data from Goncalves et al. ===

## Gene-level Dosage Compensation
### Unfiltered
factors_roc_auc_gene <- expr_buf_goncalves %>%
  analyze_roc_auc(Buffering.GeneLevel.Class)

plot_roc_auc_summary(factors_roc_auc_gene, here(goncalves_gene_plots_dir, "Unfiltered"), "buffering-factors_roc-auc_gene-level.png")
plot_roc_curves(expr_buf_goncalves, Buffering.GeneLevel.Class, here(goncalves_gene_plots_dir, "Unfiltered"),
                factor_cols = union(dc_factor_cols, c("Buffering.GeneLevel.Ratio", "Log2FC")))

auc_score_gene <- roc_auc_summary_score(factors_roc_auc_gene)

### Filtered by copy number difference
cn_diff_quantiles <- quantile(expr_buf_goncalves$Gene.CopyNumber - expr_buf_goncalves$Gene.CopyNumber.Baseline,
                              probs = seq(0, 1, 0.01))

factors_roc_auc_gene_filtered <- expr_buf_goncalves %>%
  filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles["5%"] |
           Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles["95%"]) %>%
  analyze_roc_auc(Buffering.GeneLevel.Class)

plot_roc_auc_summary(factors_roc_auc_gene_filtered, here(goncalves_gene_plots_dir, "Filtered"), "buffering-factors_roc-auc_gene-level_filtered.png")
plot_roc_curves(expr_buf_goncalves %>%
                  filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles["5%"] |
                           Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles["95%"]),
                Buffering.GeneLevel.Class, here(goncalves_gene_plots_dir, "Filtered"),
                factor_cols = union(dc_factor_cols, c("Buffering.GeneLevel.Ratio", "Log2FC")))

auc_score_gene_filtered <- roc_auc_summary_score(factors_roc_auc_gene_filtered)

## Chromosome-arm-level Dosage Compensation
### Chromosome Arm Gain
factors_roc_auc_chr_gain <- expr_buf_goncalves %>%
  filter(ChromosomeArm.CNA > 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Class)

plot_roc_auc_summary(factors_roc_auc_chr_gain, here(goncalves_chr_plots_dir, "Gain"), "buffering-factors_roc-auc_chr-level_gain.png")
plot_roc_curves(expr_buf_goncalves %>% filter(ChromosomeArm.CNA > 0),
                Buffering.ChrArmLevel.Class, here(goncalves_chr_plots_dir, "Gain"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Ratio", "Log2FC")))

auc_score_chr_gain <- roc_auc_summary_score(factors_roc_auc_chr_gain)

### Chromosome Arm Loss
factors_roc_auc_chr_loss <- expr_buf_goncalves %>%
  filter(ChromosomeArm.CNA < 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Class)

plot_roc_auc_summary(factors_roc_auc_chr_loss, here(goncalves_chr_plots_dir, "Loss"), "buffering-factors_roc-auc_chr-level_loss.png")
plot_roc_curves(expr_buf_goncalves %>% filter(ChromosomeArm.CNA < 0),
                Buffering.ChrArmLevel.Class, here(goncalves_chr_plots_dir, "Loss"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Ratio", "Log2FC")))

auc_score_chr_loss <- roc_auc_summary_score(factors_roc_auc_chr_loss)

## Chromosome-arm-level Dosage Compensation (Close to reference method of Schukken & Sheltzer)
### Chromosome Arm Gain
factors_roc_auc_chr_gain_avg <- expr_buf_goncalves %>%
  filter(ChromosomeArm.CNA > 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Average.Class)

plot_roc_auc_summary(factors_roc_auc_chr_gain_avg, here(goncalves_chr_plots_dir, "AverageGain"), "buffering-factors_roc-auc_chr-level_gain_averaged.png")
plot_roc_curves(expr_buf_goncalves %>% filter(ChromosomeArm.CNA > 0),
                Buffering.ChrArmLevel.Average.Class, here(goncalves_chr_plots_dir, "AverageGain"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Average.Ratio", "Log2FC.Average")))

auc_score_chr_gain_avg <- roc_auc_summary_score(factors_roc_auc_chr_gain_avg)

### Chromosome Arm Loss
factors_roc_auc_chr_loss_avg <- expr_buf_goncalves %>%
  filter(ChromosomeArm.CNA < 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Average.Class)

plot_roc_auc_summary(factors_roc_auc_chr_loss_avg, here(goncalves_chr_plots_dir, "AverageLoss"), "buffering-factors_roc-auc_chr-level_loss_averaged.png")
plot_roc_curves(expr_buf_goncalves %>% filter(ChromosomeArm.CNA < 0),
                Buffering.ChrArmLevel.Average.Class, here(goncalves_chr_plots_dir, "AverageLoss"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Average.Ratio", "Log2FC.Average")))

auc_score_chr_loss_avg <- roc_auc_summary_score(factors_roc_auc_chr_loss_avg)


# === Calculate ROC AUCs for data from DepMap ===

## Gene-level Dosage Compensation
### Unfiltered
factors_roc_auc_gene <- expr_buf_depmap %>%
  analyze_roc_auc(Buffering.GeneLevel.Class)

plot_roc_auc_summary(factors_roc_auc_gene, here(depmap_gene_plots_dir, "Unfiltered"), "buffering-factors_roc-auc_gene-level.png")
plot_roc_curves(expr_buf_depmap, Buffering.GeneLevel.Class, here(depmap_gene_plots_dir, "Unfiltered"),
                factor_cols = union(dc_factor_cols, c("Buffering.GeneLevel.Ratio", "Log2FC")))

auc_score_gene <- roc_auc_summary_score(factors_roc_auc_gene)

### Filtered by copy number difference
cn_diff_quantiles <- quantile(expr_buf_depmap$Gene.CopyNumber - expr_buf_depmap$Gene.CopyNumber.Baseline,
                              probs = seq(0, 1, 0.01))

factors_roc_auc_gene_filtered <- expr_buf_depmap %>%
  filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles["5%"] |
           Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles["95%"]) %>%
  analyze_roc_auc(Buffering.GeneLevel.Class)

plot_roc_auc_summary(factors_roc_auc_gene_filtered, here(depmap_gene_plots_dir, "Filtered"), "buffering-factors_roc-auc_gene-level_filtered.png")
plot_roc_curves(expr_buf_depmap %>%
                  filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles["5%"] |
                           Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles["95%"]),
                Buffering.GeneLevel.Class, here(depmap_gene_plots_dir, "Filtered"),
                factor_cols = union(dc_factor_cols, c("Buffering.GeneLevel.Ratio", "Log2FC")))

auc_score_gene_filtered <- roc_auc_summary_score(factors_roc_auc_gene_filtered)

## Chromosome-arm-level Dosage Compensation
### Chromosome Arm Gain
factors_roc_auc_chr_gain <- expr_buf_depmap %>%
  filter(ChromosomeArm.CNA > 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Class)

plot_roc_auc_summary(factors_roc_auc_chr_gain, here(depmap_chr_plots_dir, "Gain"), "buffering-factors_roc-auc_chr-level_gain.png")
plot_roc_curves(expr_buf_depmap %>% filter(ChromosomeArm.CNA > 0),
                Buffering.ChrArmLevel.Class, here(depmap_chr_plots_dir, "Gain"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Ratio", "Log2FC")))

auc_score_chr_gain <- roc_auc_summary_score(factors_roc_auc_chr_gain)

### Chromosome Arm Loss
factors_roc_auc_chr_loss <- expr_buf_depmap %>%
  filter(ChromosomeArm.CNA < 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Class)

plot_roc_auc_summary(factors_roc_auc_chr_loss, here(depmap_chr_plots_dir, "Loss"), "buffering-factors_roc-auc_chr-level_loss.png")
plot_roc_curves(expr_buf_depmap %>% filter(ChromosomeArm.CNA < 0),
                Buffering.ChrArmLevel.Class, here(depmap_chr_plots_dir, "Loss"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Ratio", "Log2FC")))

auc_score_chr_loss <- roc_auc_summary_score(factors_roc_auc_chr_loss)

## Chromosome-arm-level Dosage Compensation (Close to reference method of Schukken & Sheltzer)
### Chromosome Arm Gain
factors_roc_auc_chr_gain_avg <- expr_buf_depmap %>%
  filter(ChromosomeArm.CNA > 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Average.Class)

plot_roc_auc_summary(factors_roc_auc_chr_gain_avg, here(depmap_chr_plots_dir, "AverageGain"), "buffering-factors_roc-auc_chr-level_gain_averaged.png")
plot_roc_curves(expr_buf_depmap %>% filter(ChromosomeArm.CNA > 0),
                Buffering.ChrArmLevel.Average.Class, here(depmap_chr_plots_dir, "AverageGain"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Average.Ratio", "Log2FC.Average")))

auc_score_chr_gain_avg <- roc_auc_summary_score(factors_roc_auc_chr_gain_avg)

### Chromosome Arm Loss
factors_roc_auc_chr_loss_avg <- expr_buf_depmap %>%
  filter(ChromosomeArm.CNA < 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Average.Class)

plot_roc_auc_summary(factors_roc_auc_chr_loss_avg, here(depmap_chr_plots_dir, "AverageLoss"), "buffering-factors_roc-auc_chr-level_loss_averaged.png")
plot_roc_curves(expr_buf_depmap %>% filter(ChromosomeArm.CNA < 0),
                Buffering.ChrArmLevel.Average.Class, here(depmap_chr_plots_dir, "AverageLoss"),
                factor_cols = union(dc_factor_cols, c("Buffering.ChrArmLevel.Average.Ratio", "Log2FC.Average")))

auc_score_chr_loss_avg <- roc_auc_summary_score(factors_roc_auc_chr_loss_avg)