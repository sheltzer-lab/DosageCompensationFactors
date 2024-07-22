library(dplyr)
library(tidyr)
library(here)
library(broom)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))

# Add dosage compensation factors to a dataframe
add_factors <- function(df, df_factors, factor_cols = dc_factor_cols) {
  df %>%
    inner_join(y = df_factors %>% select(Protein.Uniprot.Accession, Gene.Symbol, all_of(factor_cols)),
              by = c("Protein.Uniprot.Accession", "Gene.Symbol"),
              na_matches = "never", relationship = "many-to-one")
}

# === Filtering ===
## Quantile Filtering
remove_outer_quantiles <- function(df, col, keep_between = c("5%", "95%")) {
  quantiles <- quantile(df[[quo_name(enquo(col))]], probs = seq(0, 1, 0.01), na.rm = TRUE)
  df %>%
    filter({ { col } } > quantiles[keep_between[1]] & { { col } } < quantiles[keep_between[2]])
}

split_by_quantiles <- function(df, col, quantile_low = "20%", quantile_high = "80%",
                               target_group_col = "SplitGroup") {
  quantiles <- quantile(df[[quo_name(enquo(col))]], probs = seq(0, 1, 0.01), na.rm = TRUE)
  df %>%
    filter({ { col } } < quantiles[quantile_low] | { { col } } > quantiles[quantile_high]) %>%
    mutate(!!target_group_col := if_else({ { col } } < quantiles[quantile_low], "Low", "High"))
}

split_by_3_quantiles <- function(df, col, quantile_range = 20L, target_group_col = "SplitGroup") {
  quantiles <- quantile(df[[quo_name(enquo(col))]], probs = seq(0, 1, 0.01), na.rm = TRUE)
  low_quantile <- quantiles[paste0(quantile_range, "%")]
  high_quantile <- quantiles[paste0(100 - quantile_range, "%")]
  center_quantiles <- c(quantiles[50L - ceiling(quantile_range / 2)],
                        quantiles[50L + floor(quantile_range / 2)])

  df %>%
    mutate(Group_ = if_else({ { col } } < low_quantile, "Low",
                                         if_else({ { col } } > high_quantile, "High",
                                                 if_else({ { col } } > center_quantiles[1] & { { col } } < center_quantiles[2], "Center",
                                                         NA)))) %>%
    drop_na(Group_) %>%
    rename(!!target_group_col := Group_)
}

## Gene / Chromosome Copy Number Filtering
filter_cn_diff_quantiles <- function(df, remove_between = c("5%", "95%")) {
  cn_diff_quantiles <- quantile(df$Gene.CopyNumber - df$Gene.CopyNumber.Baseline, probs = seq(0, 1, 0.01), na.rm = TRUE)
  df %>%
    filter(abs(Gene.CopyNumber - Gene.CopyNumber.Baseline) > 0) %>%
    filter(Gene.CopyNumber < Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_between[1]] |
             Gene.CopyNumber > Gene.CopyNumber.Baseline + cn_diff_quantiles[remove_between[2]])
}

filter_cn_diff <- function(df, remove_between = c(-0.02, 0.02)) {
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

filter_cn_gain_abs <- function(df) {
  df %>% filter(Gene.CopyNumber > Gene.CopyNumber.Baseline)
}

filter_cn_loss_abs <- function(df) {
  df %>% filter(Gene.CopyNumber < Gene.CopyNumber.Baseline)
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
    pivot_wider(id_cols = c("Model.ID", "Gene.Symbol", "Protein.Uniprot.Accession"),
                names_from = quo_name(enquo(dataset_col)),
                values_from = quo_name(enquo(value_col))) %>%
    rename(DatasetA = 4,
           DatasetB = 5) %>%
    group_by(Model.ID) %>%
    summarize(Correlation = cor(DatasetA, DatasetB,
                                method = method, use = "na.or.complete"),
              Observations = min(c(sum(!is.na(DatasetA)), sum(!is.na(DatasetB))))
    ) %>%
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
    distinct({ { sample_col } }, { { sample_group_col } }, .keep_all = TRUE)

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

# TODO: Move to evaluation
data_density <- function(df) {
  df %>%
    summarize_all(~sum(!is.na(.x)) / nrow(df)) %>%
    pivot_longer(everything(), names_to = "Column", values_to = "Density") %>%
    arrange(desc(Density))
}

z_score <- function(values) {
  return((values - mean(values, na.rm = TRUE)) / sd(values, na.rm = TRUE))
}

differential_expression <- function(df, id_col, group_col, expr_col,
                                    paired = FALSE, p.adj.method = "BH", groups = NULL) {

  if (length(groups) != 2) {
    warning("'groups' parameter does not contain exactly two values.
    Groups will be defined from unique values in 'group_col'.")
    groups <- unique(df[[quo_name(enquo(expr_col))]])
  }

  df %>%
    drop_na() %>%
    group_by({ { id_col } }) %>%
    add_count({ { group_col } }) %>%                    # Count observations per group and protein
    filter(n > 2) %>%                                   # Remove proteins with insufficient observations
    filter(length(unique({ { group_col } })) > 1) %>%   # Remove leftover groups
    summarise(
      GroupA = groups[1],
      GroupB = groups[2],
      Mean_GroupA = mean({ { expr_col } }[{ { group_col } } == groups[1]]),
      Mean_GroupB = mean({ { expr_col } }[{ { group_col } } == groups[2]]),
      Log2FC = Mean_GroupB - Mean_GroupA,
      TTest.p = t.test(as.formula(paste(quo_name(enquo(expr_col)), "~", quo_name(enquo(group_col)))),
                       paired = paired)$p.value,
      .groups = 'drop'
    ) %>%
    mutate(TTest.p.adj = p.adjust(TTest.p, method = p.adj.method))
}

# === (Parallelized) Bootstrapping Methods ===

sample_func <- function(i, df, sample_prop, func, ...) {
  require(dplyr)

  set.seed(i)

  df %>%
    slice_sample(prop = sample_prop, replace = TRUE) %>%
    func(...) %>%
    mutate(Bootstrap.Sample = as.integer(i))
}

bootstrap_dataframe <- function (df, cluster, n, sample_prop, func, ...) {
  require(dplyr)

  results <- list()

  # Non-parallelized with progress bar
  if (is.null(cluster)) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
      results[[i]] <- sample_func(i, df, sample_prop, func, ...)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  # Parallelized without progress bar
  else {
    results <- parLapply(cluster, 1:n, sample_func, df, sample_prop, func, ...)
  }

  return(bind_rows(results))
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
    summarize(MeanNormRank = mean(Rank)) %>%
    mutate(MeanNormRank = as.numeric(MeanNormRank),
           { { id_col } } := factor({ { id_col } },
                                    levels = { { id_col } }[order(MeanNormRank)])) %>%
    arrange(desc(MeanNormRank))
}

standardized_mean <- function (df, value_col, group_col, id_col) {
  df %>%
    group_by({{group_col}}) %>%
    mutate(Scaled = scale({ { value_col } })[,1]) %>%
    ungroup() %>%
    group_by({{id_col}}) %>%
    summarize(StandardizedMean = mean(Model.Buffering.Ratio, na.rm = TRUE))
}

# === SHAP Analysis ===
# Only override defaults if enough RAM available
estimate_shap <- function(model, n_samples = 100, n_combinations = 250, method = "ctree") {
  require(shapr)
  require(dplyr)

  set.seed(42)
  # Prepare the data for explanation
  if (is.null(model$datasets$test)) {
    x_small <- model$datasets$training %>%
      tibble::rownames_to_column("ID") %>%
      group_by(buffered) %>%
      slice_sample(n = 2 * n_samples) %>%
      ungroup() %>%
      tibble::column_to_rownames("ID") %>%
      split_train_test(training_set_ratio = 0.5)  # From preprocessing.R

    training_x_small <- x_small$training
    test_x_small <- x_small$test
  } else {
    training_x_small <- model$datasets$training %>%
      group_by(buffered) %>%
      slice_sample(n = n_samples) %>%
      ungroup() %>%
      select(-buffered)
    test_x_small <- model$datasets$test %>%
      tibble::rownames_to_column("ID") %>%
      group_by(buffered) %>%
      slice_sample(n = n_samples) %>%
      ungroup() %>%
      tibble::column_to_rownames("ID") %>%
      select(-buffered)
  }

  # https://cran.r-project.org/web/packages/shapr/vignettes/understanding_shapr.html
  explainer <- shapr(training_x_small, model$finalModel, n_combinations = n_combinations)
  # Specifying the phi_0, i.e. the expected prediction without any features
  # Note: buffered column is a factor: 1 = "Buffered", 2 = "Scaling"
  # Lower (SHAP) values means prediction is closer to Buffered
  p <- mean(as.numeric(model$datasets$training$buffered))

  # Default method: ctree, as counted values are numerical but not continuous
  explanation <- explain(
    test_x_small,
    approach = method,
    explainer = explainer,
    prediction_zero = p
  )
  explanation$x_test <- test_x_small

  return(explanation)
}

shap2df <- function(explanation) {
  require(dplyr)
  require(tidyr)

  factor_values <- explanation$x_test %>%
    tibble::rownames_to_column("ID") %>%
    pivot_longer(c(everything(), -ID),
                 names_to = "DosageCompensation.Factor",
                 values_to = "DosageCompensation.Factor.Value")

  df_explanation <- explanation$dt %>%
    select(-none) %>%
    pivot_longer(everything(), names_to = "DosageCompensation.Factor", values_to = "SHAP.Value") %>%
    mutate(ID = factor_values$ID,
           Factor.Value = factor_values$DosageCompensation.Factor.Value) %>%
    group_by(DosageCompensation.Factor) %>%
    mutate(Factor.Value.Relative = normalize_min_max(Factor.Value, na.rm = TRUE),
           SHAP.p25.Absolute = quantile(abs(SHAP.Value), probs = 0.25)[["25%"]],
           SHAP.Median.Absolute = median(abs(SHAP.Value)),
           SHAP.p75.Absolute = quantile(abs(SHAP.Value), probs = 0.75)[["75%"]]) %>%
    ungroup()

  df_corr <- df_explanation %>%
    group_by(DosageCompensation.Factor) %>%
    rstatix::cor_test(Factor.Value, SHAP.Value, method = "spearman", use = "na.or.complete") %>%
    select(DosageCompensation.Factor, cor, p) %>%
    rename(SHAP.Factor.Corr = cor, SHAP.Factor.Corr.p = p) %>%
    mutate(SHAP.Factor.Corr.p.adj = p.adjust(SHAP.Factor.Corr.p, method = "BY"))

  df_explanation <- df_explanation %>%
    left_join(y = df_corr, by = "DosageCompensation.Factor",
              relationship = "many-to-one", unmatched = "error", na_matches = "never") %>%
    select(ID, DosageCompensation.Factor, starts_with("Factor"), starts_with("SHAP"))

  return(df_explanation)
}