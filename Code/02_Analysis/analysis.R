# Dosage Compensation Factors
dc_factor_cols <- c(
  "Protein-Protein Interactions", "Protein Half-Life", "Protein Complexes (CORUM)",
  "Mean 3'-UTR Length", "Mean 5'-UTR Length",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites",
  "Methylation Sites", "Acetylation Sites", "Regulatory Sites",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate",
  "Translation Rate", "Protein Length", "mRNA Length",
  "Intrinsic Protein Disorder", "Low Complexity Score", "Homology Score",
  "Loops In Protein Score", "Protein Polyampholyte Score", "Protein Polarity",
  "Non-Exponential Decay Delta", "Mean mRNA Decay Rate", "Aggregation Score"
)

## Dataset-Specific Dosage Compensation Factors
dc_factor_cols_specific <- c("Protein Neutral CV")

# Dataset Filter Functions
filter_cn_diff_quantiles <- function(df, remove_between = c("5%", "95%")) {
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01))
  df %>%
    filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_between[1]] |
             Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_between[2]])
}

filter_cn_diff <- function(df, remove_between = c(-0.01, 0.01)) {
  df %>%
    filter(Gene.CopyNumber - Gene.CopyNumber.Baseline < remove_between[1] |
             Gene.CopyNumber - Gene.CopyNumber.Baseline > remove_between[2])
}

filter_cn_gain <- function(df, remove_below = "90%") {
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01))
  df %>%
    filter(Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_below])
}

filter_cn_loss <- function(df, remove_above = "10%") {
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01))
  df %>%
    filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_above])
}

filter_arm_gain <- function(df) {
  df %>% filter(ChromosomeArm.CNA > 0)
}

filter_arm_loss <- function(df) {
  df %>% filter(ChromosomeArm.CNA < 0)
}

filter_arm_gain_gene_avg <- function(df) {
  df %>%
    filter(ChromosomeArm.CNA > 0) %>%
    distinct(Gene.Symbol, .keep_all = TRUE)
}

filter_arm_loss_gene_avg <- function(df) {
  df %>%
    filter(ChromosomeArm.CNA < 0) %>%
    distinct(Gene.Symbol, .keep_all = TRUE)
}

# Parameter search
grid_search <- function(func, param_range) {
  values <- NULL
  pb <- txtProgressBar(min = min(param_range), max = max(param_range), style = 3)

  for (param in param_range) {
    setTxtProgressBar(pb, param)
    new_value <- func(param)
    values <- c(values, new_value)
  }

  names(values) <- param_range
  close(pb)

  return(values)
}