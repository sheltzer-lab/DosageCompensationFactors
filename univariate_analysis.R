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
plots_dir <- plots_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_avg <- read_parquet(here(output_data_dir, "expression_average.parquet"))
dc_factors <- read_parquet(here(output_data_dir, "dosage_compensation_factors.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet"))

# === Calculate Buffering & Dosage Compensation ===

calculate_baseline_copynumber <- function(df, gene_col, copynumber_col) {
  df %>%
    group_by({{ gene_col }}) %>%
    mutate(Gene.CopyNumber.Baseline = median({{ copynumber_col }}, na.rm = TRUE)) %>%
    ungroup()
}

calculate_baseline_expression <- function(df, gene_col, chr_arm_cna_col, expr_col) {
  baseline_expr <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { expr_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    summarize(Protein.Expression.Baseline = median({ { expr_col } }, na.rm = TRUE)) %>%
    ungroup() %>%
    select({ { gene_col } }, Protein.Expression.Baseline)

  df %>%
    inner_join(y = baseline_expr, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one")
}

expr_avg_dc <- expr_avg %>%
  inner_join(y = copy_number, by = c("CellLine.SangerModelId", "Gene.Symbol"),
             na_matches = "never", relationship = "one-to-one") %>%
  calculate_baseline_copynumber(Gene.Symbol, Gene.CopyNumber) %>%
  calculate_baseline_expression(Gene.Symbol, ChromosomeArm.CNA, Protein.Expression.Normalized) %>%
  mutate(Log2FC = log2(Protein.Expression.Normalized) - log2(Protein.Expression.Baseline)) %>%
  mutate(Buffering.GeneLevel.Ratio = buffering_ratio(Protein.Expression.Baseline, Protein.Expression.Normalized,
                                                     Gene.CopyNumber.Baseline, Gene.CopyNumber)) %>%
  mutate(Buffering.GeneLevel.Class = buffering_class(Buffering.GeneLevel.Ratio)) %>%
  mutate(Buffering.ChrArmLevel.Class = buffering_class_log2fc(Log2FC,
                                                              cn_base = ChromosomeArm.CopyNumber.Baseline,
                                                              cn_var = ChromosomeArm.CopyNumber)) %>%
  left_join(y = dc_factors, by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one")

# === Calculate ROC AUC of factors ===
dc_factor_cols <- c(
  "Protein-Protein Interactions", "Protein Half-Life",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites",
  "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate",
  "Translation Rate", "Protein Length", "mRNA Length"
)

analyze_roc_auc <- function(df, buffering_class_col, factor_cols = dc_factor_cols, id_col = "UniqueId") {
  df %>%
    mutate(Buffered = ifelse({ { buffering_class_col } } == "Buffered", 1, 0)) %>%
    mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
    drop_na(Buffered) %>%
    select(all_of(id_col), Buffered, all_of(factor_cols)) %>%
    pivot_longer(all_of(factor_cols),
                 names_to = "DosageCompensation.Factor",
                 values_to = "DosageCompensation.Factor.Value") %>%
    group_by(DosageCompensation.Factor) %>%
    summarize(DosageCompensation.Factor.ROC_AUC =
                list(auc(Buffered, DosageCompensation.Factor.Value, na.rm = TRUE)) %>% unlist()) %>%
    mutate(DosageCompensation.Factor =
             factor(DosageCompensation.Factor,
                    levels = DosageCompensation.Factor[order(DosageCompensation.Factor.ROC_AUC)])) %>%
    arrange(DosageCompensation.Factor.ROC_AUC)
}

plot_roc_auc_summary <- function(factors_roc_auc, filename) {
  roc_auc_summary_plot <- factors_roc_auc %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = DosageCompensation.Factor.ROC_AUC,
        label = format(round(DosageCompensation.Factor.ROC_AUC, 3), nsmall = 3)) +
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
  mean(abs(df$DosageCompensation.Factor.ROC_AUC - 0.5))
}

## Gene-level Dosage Compensation
### Unfiltered
factors_roc_auc_gene <- expr_avg_dc %>%
  analyze_roc_auc(Buffering.GeneLevel.Class)

plot_roc_auc_summary(factors_roc_auc_gene, "buffering-factors_roc-auc_gene-level.png")

auc_score_gene <- roc_auc_summary_score(factors_roc_auc_gene)

### Filtered by copy number difference
cn_diff_quantiles <- quantile(expr_avg_dc$Gene.CopyNumber - expr_avg_dc$Gene.CopyNumber.Baseline,
                              probs = seq(0, 1, 0.05))

factors_roc_auc_gene_filtered <- expr_avg_dc %>%
  filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles["5%"] |
           Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles["95%"]) %>%
  analyze_roc_auc(Buffering.GeneLevel.Class)

plot_roc_auc_summary(factors_roc_auc_gene_filtered, "buffering-factors_roc-auc_gene-level_filtered.png")

auc_score_gene_filtered <- roc_auc_summary_score(factors_roc_auc_gene_filtered)

## Chromosome-arm-level Dosage Compensation
### Chromosome Arm Gain
factors_roc_auc_chr_gain <- expr_avg_dc %>%
  filter(ChromosomeArm.CNA > 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Class)

plot_roc_auc_summary(factors_roc_auc_chr_gain, "buffering-factors_roc-auc_chr-level_gain.png")

auc_score_chr_gain <- roc_auc_summary_score(factors_roc_auc_chr_gain)

### Chromosome Arm Loss
factors_roc_auc_chr_loss <- expr_avg_dc %>%
  filter(ChromosomeArm.CNA < 0) %>%
  analyze_roc_auc(Buffering.ChrArmLevel.Class)

plot_roc_auc_summary(factors_roc_auc_chr_loss, "buffering-factors_roc-auc_chr-level_loss.png")

auc_score_chr_loss <- roc_auc_summary_score(factors_roc_auc_chr_loss)


# === Write Processed dataset to disk ===
write_parquet(expr_avg_dc, here(output_data_dir, 'expression_average_buffering.parquet'),
              version = "2.6")