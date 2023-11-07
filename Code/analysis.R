library(dplyr)
library(tidyr)
library(here)
library(broom)

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
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01), na.rm = TRUE)
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
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01), na.rm = TRUE)
  df %>%
    filter(Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_below])
}

filter_cn_loss <- function(df, remove_above = "10%") {
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01), na.rm = TRUE)
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

# Calculate correlation between datasets (cell lines as separate samples)
dataset_correlation <- function(df, dataset_col, value_col, comparison_name, method = "pearson") {
  test <- df %>%
    pivot_wider(id_cols = c("CellLine.CustomId", "Gene.Symbol", "Protein.Uniprot.Accession"),
                names_from = quo_name(enquo(dataset_col)),
                values_from = quo_name(enquo(value_col))) %>%
    rename(DatasetA = 4,
           DatasetB = 5) %>%
    group_by(CellLine.CustomId) %>%
    summarize(Correlation = cor(DatasetA, DatasetB,
                                method = method, use = "na.or.complete")) %>%
    mutate(Comparison = comparison_name)
}

calculate_pca <- function(df, sample_col, sample_group_col, value_group_col, value_col) {
  pca_fit <- df %>%
    select({ { sample_col } }, { { sample_group_col } }, { { value_group_col } }, { { value_col } }) %>%
    pivot_wider(names_from = { { sample_col } }, values_from = { { value_col } },
                id_cols = { { value_group_col } }) %>%
    drop_na() %>%
    select(where(is.numeric)) %>%
    tibble::rownames_to_column() %>%
    pivot_longer(-rowname) %>%
    pivot_wider(names_from = rowname, values_from = value) %>%
    tibble::column_to_rownames(var = "name") %>%
    scale() %>%
    prcomp()

  sample_metadata <- df %>%
    distinct({ { sample_col } }, { { sample_group_col } })

  df_pca <- pca_fit %>%
    augment(sample_metadata) %>%
    select(-.rownames)

  eigenvalues <- pca_fit %>%
    tidy(matrix = "eigenvalues")

  return(list(pca = pca_fit, df_pca = df_pca, eigenvalues = eigenvalues))
}

# Resample the target distribution in a stratified way using a reference distribution,
# so that both distributions approximately the same
equalize_distributions <- function(df_reference, df_target, value_col,
                                   with_replacement = TRUE, num_buckets = 10) {
  set.seed(42)

  df_resample <- function(df, group_key) {
    df %>%
      slice_sample(prop = first(df$Value.Prop.Adj), replace = with_replacement)
  }

  buckets <- unique(quantile(df_reference[[quo_name(enquo(value_col))]],
                             probs = seq(0, 1, 1 / num_buckets)))

  df_ref_buckets <- df_reference %>%
    mutate(Value.Bucket = cut({ { value_col } }, breaks = buckets, include.lowest = TRUE)) %>%
    count(Value.Bucket, name = "Value.Count") %>%
    mutate(Value.Prop = Value.Count / nrow(df_reference))

  df_joined <- df_target %>%
    mutate(Value.Bucket = cut({ { value_col } }, breaks = buckets, include.lowest = TRUE)) %>%
    inner_join(y = df_ref_buckets, by = "Value.Bucket")

  df_adjusted <- df_joined %>%
    add_count(Value.Bucket, name = "Value.Count.New") %>%
    mutate(Value.Prop.New = Value.Count.New / nrow(df_joined)) %>%
    mutate(Value.Prop.Adj = Value.Prop / Value.Prop.New) %>%
    group_by(Value.Bucket) %>%
    group_modify(df_resample) %>%
    ungroup()

  df_equal <- bind_rows(df_reference, df_adjusted) %>%
    select(-Value.Bucket, -Value.Prop, -Value.Prop.New, -Value.Prop.Adj,
           -Value.Count, -Value.Count.New)

  return(df_equal)
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

standardized_mean <- function (df, value_col, group_col, id_col) {
  df %>%
    group_by({{group_col}}) %>%
    mutate(Scaled = scale({ { value_col } })[,1]) %>%
    ungroup() %>%
    group_by({{id_col}}) %>%
    summarize(StandardizedMean = mean(Buffering.CellLine.Ratio, na.rm = TRUE))
}