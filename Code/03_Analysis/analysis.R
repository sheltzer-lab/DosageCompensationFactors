library(dplyr)
library(tidyr)
library(here)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

# Add dosage compensation factors to a dataframe
add_factors <- function(df, df_factors, factor_cols = dc_factor_cols) {
  df %>%
    left_join(y = df_factors %>% select(Protein.Uniprot.Accession, Gene.Symbol, all_of(factor_cols)),
              by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
              na_matches = "never", relationship = "many-to-one")
}

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

# Convert list of ROC objects to data frame
rocs_to_df <- function(rocs) {
  df_rocs <- data.frame()
  # Extract relevant information from ROC objects into data frame
  for (name in names(rocs)) {
    current_roc <- rocs[[name]]
    df_current_roc <- data.frame(
      Name = rep(name, length(current_roc$specificities)),
      Specificity = current_roc$specificities,
      Sensitivity = current_roc$sensitivities,
      AUC = rep(auc(current_roc), length(current_roc$specificities))
    )
    df_rocs <- rbind(df_rocs, df_current_roc)
  }
  return(df_rocs)
}

# === Aggregation & Consensus Methods ===

mean_norm_rank <- function(df, value_col, group_col, id_col) {
  total_groups <- length(unique(df[[quo_name(enquo(group_col))]]))

  df %>%
    select({ { value_col } }, { { group_col } }, { { id_col } }) %>%
    # Calculate Normalized Ranks
    group_by({ { group_col } }) %>%
    add_count(name = "n.Features") %>%
    mutate(Rank = as.integer(rank({ { value_col } })) / n.Features) %>%
    ungroup() %>%
    # Remove entries that don't have values for more than half of the samples
    group_by({ { id_col } }) %>%
    add_count(name = "n.Samples") %>%
    filter(n.Samples >= round(total_groups / 2)) %>%
    # Create NAs for missing samples
    select(Rank, { { group_col } }, { { id_col } }) %>%
    pivot_wider(values_from = Rank,
                names_from = quo_name(enquo(group_col)),
                id_cols = quo_name(enquo(id_col))) %>%
    pivot_longer(-{ { id_col } },
                 values_to = "Rank",
                 names_to = quo_name(enquo(group_col))) %>%
    # Impute missing data
    mutate(Rank = replace_na(Rank, 0.5)) %>%
    # Aggregate normalized ranks
    summarize(AggregatedRank = mean(Rank)) %>%
    mutate(AggregatedRank = as.numeric(AggregatedRank),
           { { id_col } } := factor({ { id_col } },
                                    levels = { { id_col } }[order(AggregatedRank)])) %>%
    arrange(desc(AggregatedRank))
}