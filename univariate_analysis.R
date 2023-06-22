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

calculate_baseline_expression <- function(df, gene_col, expr_col) {
  df %>%
    group_by({{ gene_col }}) %>%
    mutate(Protein.Expression.Baseline = median({{ expr_col }}, na.rm = TRUE)) %>%
    ungroup()
}

expr_avg_dc <- expr_avg %>%
  inner_join(y = copy_number, by = c("CellLine.SangerModelId", "Gene.Symbol"),
             na_matches = "never", relationship = "one-to-one") %>%
  calculate_baseline_copynumber(Gene.Symbol, Gene.CopyNumber) %>%
  calculate_baseline_expression(Gene.Symbol, Protein.Expression.Normalized) %>%
  mutate(Buffering.Ratio = buffering_ratio(Protein.Expression.Baseline, Protein.Expression.Normalized,
                                           Gene.CopyNumber.Baseline, Gene.CopyNumber)) %>%
  mutate(Buffering.Class = buffering_class(Buffering.Ratio)) %>%
  mutate(Buffered = ifelse(Buffering.Class == "Buffered", 1, 0)) %>%
  mutate(Buffered = factor(Buffered, levels = c(0, 1))) %>%
  left_join(y = dc_factors, by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
            na_matches = "never", relationship = "many-to-one")

# === Calculate ROC AUC of factors
dc_factor_names <- c(
  "Protein-Protein Interactions", "Protein Half-Life",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites", "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate", "Translation Rate", "Protein Length", "mRNA Length"
)

factors_roc_auc <- expr_avg_dc %>%
  drop_na(Buffered) %>%
  select(UniqueId, Buffered, all_of(dc_factor_names)) %>%
  pivot_longer(all_of(dc_factor_names),
               names_to = "DosageCompensation.Factor",
               values_to = "DosageCompensation.Factor.Value") %>%
  group_by(DosageCompensation.Factor) %>%
  summarize(DosageCompensation.Factor.ROC_AUC = list(auc(Buffered, DosageCompensation.Factor.Value, na.rm = TRUE)) %>% unlist()) %>%
  mutate(DosageCompensation.Factor = factor(DosageCompensation.Factor,
                                            levels = DosageCompensation.Factor[order(DosageCompensation.Factor.ROC_AUC)])) %>%
  arrange(DosageCompensation.Factor.ROC_AUC)

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

ggsave(here(plots_dir, "buffering_factors_roc_auc.png"), plot = roc_auc_summary_plot,
       height = 200, width = 180, units = "mm", dpi = 300)

avg_abs_auc <- factors_roc_auc %>%
  summarize(Score = mean(abs(DosageCompensation.Factor.ROC_AUC - 0.5)))

# === Write Processed dataset to disk ===

write_parquet(expr_avg_dc, here(output_data_dir, 'expression_average_dc.parquet'),
              version = "2.6")