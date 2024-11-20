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
  require(pROC)

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
                                    test = t.test, paired = FALSE, centr = mean, p.adj.method = "BH", groups = NULL,
                                    log2fc_thresh = log2fc_threshold, p_thresh = p_threshold) {

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
      Mean_GroupA = centr({ { expr_col } }[{ { group_col } } == groups[1]]),
      Mean_GroupB = centr({ { expr_col } }[{ { group_col } } == groups[2]]),
      Count_GroupA = sum(!is.na({ { expr_col } }[{ { group_col } } == groups[1]])),
      Count_GroupB = sum(!is.na({ { expr_col } }[{ { group_col } } == groups[2]])),
      Log2FC = Mean_GroupB - Mean_GroupA,
      Test.p = test(as.formula(paste(quo_name(enquo(expr_col)), "~", quo_name(enquo(group_col)))),
                    paired = paired)$p.value,
      .groups = 'drop'
    ) %>%
    mutate(Test.p.adj = p.adjust(Test.p, method = p.adj.method),
           Significant = case_when(Log2FC > log2fc_thresh & Test.p.adj < p_thresh ~ "Up",
                                   Log2FC < -log2fc_thresh & Test.p.adj < p_thresh ~ "Down",
                                   TRUE ~ NA))
}

analyze_low_br_variance <- function (df_expr_buf) {
  source(here("Code", "buffering_ratio.R"))

  # Per gene summarize Mean, SD, and share of samples buffered across dataset
  mean_var_buf <- df_expr_buf %>%
    select(Dataset, Gene.Symbol, Buffering.GeneLevel.Ratio, Buffering.GeneLevel.Class) %>%
    drop_na(Buffering.GeneLevel.Ratio) %>%
    group_by(Dataset, Gene.Symbol) %>%
    summarize(Gene.BR.Mean = mean(Buffering.GeneLevel.Ratio),
              Gene.BR.SD = sd(Buffering.GeneLevel.Ratio),
              Gene.Buffered.Count = sum(Buffering.GeneLevel.Class == "Buffered"),
              Samples = sum(!is.na(Buffering.GeneLevel.Class)),
              Gene.Buffered.Share = Gene.Buffered.Count / Samples,
              .groups = "drop")

  # Filter genes by number of samples, share of buffered samples, and mean buffering ratio and rank by SD of BR
  buf_rank <- mean_var_buf %>%
    filter(Samples > 20) %>%
    filter(Gene.Buffered.Share > 1 / 3 & Gene.BR.Mean > br_cutoffs$Buffered & Gene.BR.SD < 2) %>%
    add_count(Gene.Symbol) %>%
    filter(n >= length(unique(df_expr_buf$Dataset))) %>%
    mean_norm_rank(Gene.BR.SD, Dataset, Gene.Symbol) %>%
    rename(Gene.BR.SD.MeanNormRank = MeanNormRank)

  # Rank genes by lowest BR variance
  mean_var_buf %>%
    left_join(y = buf_rank, by = "Gene.Symbol") %>%
    mutate(Top50 = Gene.Symbol %in% slice_min(buf_rank, Gene.BR.SD.MeanNormRank, n = 50)$Gene.Symbol,
           Top10 = Gene.Symbol %in% slice_min(buf_rank, Gene.BR.SD.MeanNormRank, n = 10)$Gene.Symbol) %>%
    arrange(Gene.BR.SD.MeanNormRank) %>%
    return()
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

# === Statistically compare bootstrapped results between conditions ===
# Notes for Statistical Tests and Multiple Testing Correction
# * Bootstrapping values independently samples values from the same distribution (with replacement) => iid.
# * However: ROC AUC of factors may be correlated -> Not independent!
compare_conditions <- function(df_condition1, df_condition2) {
  require(dplyr)
  require(tidyr)
  require(skimr)
  require(assertr)

  if (nrow(df_condition1) == 0 || nrow(df_condition2) == 0) return(NULL)

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
    # Introduce perturbation to avoid equal ranks, otherwise p-value can't be calculated accurately
    mutate(RankValue = numeric.p50 + runif(length(numeric.p50), min = -1e-10, max = 1e-10)) %>%
    group_by(Condition) %>%
    mutate(DosageCompensation.Factor.Rank = as.integer(rank(-RankValue))) %>%
    select(Condition, DosageCompensation.Factor, DosageCompensation.Factor.Rank) %>%
    pivot_wider(id_cols = DosageCompensation.Factor,
                names_from = Condition, values_from = DosageCompensation.Factor.Rank)

  rank_test <- cor.test(df_rank_test[[conditions[1]]],
                        df_rank_test[[conditions[2]]],
                        method = "kendall")

  return(list(factor_test = df_factor_test,
              rank_test = rank_test,
              stat_summary = df_stat))
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
estimate_shap <- function(model, n_samples = 100, n_combinations = 250, method = "ctree", remove_incorrect = TRUE) {
  require(shapr)
  require(dplyr)

  test_set <- model$datasets$test
  training_set <- model$datasets$training

  # Remove incorrect predictions from test set
  if (!is.null(test_set) && remove_incorrect) {
    test_set <- test_set %>%
      mutate(Prediction = model$evaluation$predictedResponse$Prediction) %>%
      filter(buffered == Prediction) %>%
      select(-Prediction)
  }

  # Prepare the data for explanation
  set.seed(42)
  if (is.null(test_set)) {
    x_small <- training_set %>%
      tibble::rownames_to_column("ID") %>%
      group_by(buffered) %>%
      slice_sample(n = 2 * n_samples) %>%
      ungroup() %>%
      tibble::column_to_rownames("ID") %>%
      split_train_test(training_set_ratio = 0.5)  # From preprocessing.R
    training_x_small <- x_small$training
    test_x_small <- x_small$test
  } else {
    training_x_small <- training_set %>%
      group_by(buffered) %>%
      slice_sample(n = n_samples) %>%
      ungroup()
    test_x_small <- test_set %>%
      tibble::rownames_to_column("ID") %>%
      group_by(buffered) %>%
      slice_sample(n = n_samples) %>%
      ungroup() %>%
      tibble::column_to_rownames("ID")
  }
  training_x_small <- training_x_small %>% select(-buffered)
  test_x_small <- test_x_small %>% select(-buffered)

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

# === Enrichment Analysis ===
single_gene_set_enrichment <- function(df, gene_set, score_col, group_col, target_group, sample_col,
                                       gene_col = Gene.Symbol, ...) {
  # https://rdrr.io/bioc/limma/man/roast.html
  require(limma)

  data <- df %>%
    select({ { gene_col } }, { { sample_col } }, { { score_col } }) %>%
    drop_na() %>%
    group_by({ { gene_col } }, { { sample_col } }) %>%
    summarize(Score = mean({ { score_col } }), .groups = "drop") %>%
    pivot_wider(names_from = { { sample_col } }, values_from = Score, id_cols = { { gene_col } }) %>%
    na.omit()

  design <- df %>%
    distinct({ { sample_col } }, { { group_col } }) %>%
    filter({ { sample_col } } %in% colnames(data)) %>%
    mutate(Group = as.integer({ { group_col } } == target_group),
           Intercept = 1) %>%
    arrange({ { sample_col } }) %>%
    select(Intercept, Group) %>%
    as.matrix()

  index <- data %>%
    mutate(InSet = { { gene_col } } %in% gene_set) %$%
    which(InSet)

  data_mat <- data %>%
    select(-{ { gene_col } }) %>%
    as.matrix()

  roast(data_mat, index, design, contrast = 2, ...)
}

ssGSEA <- function(df, gene_sets, expr_col, sample_col,
                   gene_col = Gene.Symbol, method = "ssGSEA") {
  require(GSVA)
  require(dplyr)
  # https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html

  paramFunc <- switch(method,
                      "ssGSEA" = GSVA::ssgseaParam,
                      "GSVA" = GSVA::gsvaParam,
                      "PLAGE" = GSVA::plageParam,
                      "zscore" = GSVA::zscoreParam)

  data_mat <- df %>%
    select({ { gene_col } }, { { sample_col } }, { { expr_col } }) %>%
    drop_na() %>%
    group_by({ { gene_col } }, { { sample_col } }) %>%
    summarize(Score = mean({ { expr_col } }), .groups = "drop") %>%
    pivot_wider(names_from = { { sample_col } }, values_from = Score, id_cols = { { gene_col } }) %>%
    na.omit() %>%
    tibble::column_to_rownames(quo_name(enquo(gene_col))) %>%
    as.matrix()

  gseaPar <- paramFunc(data_mat, gene_sets)
  gsea_result <- gsva(gseaPar)
}

create_string_network <- function(df, gene_col, logfc_col, string_db, min_score = 700) {
  require(STRINGdb)
  require(dplyr)
  require(magrittr)
  require(scales)

  max_val <- df %>% pull({ { logfc_col } }) %>% abs() %>% max()
  domain <- c(-max_val, max_val)
  color_func <- scales::col_numeric(palette = bidirectional_color_pal, domain = domain)

  df_mapped <- df %>%
    select({ { gene_col } }, { { logfc_col } }) %>%
    as.data.frame() %>%
    string_db$map(quo_name(enquo(gene_col)), removeUnmappedRows = TRUE) %>%
    mutate(Color = color_func({ { logfc_col } }))

  payload <- df_mapped %$%
    string_db$post_payload(STRING_id, colors = Color)

  df_mapped %$%
    return(list(df_mapped = df_mapped,
                payload = payload,
                network = string_db$get_subnetwork(STRING_id),
                link = string_db$get_link(STRING_id, payload_id = payload, required_score = min_score)))
}

overrepresentation_analysis <- function(genes, ordered = TRUE, p_thresh = p_threshold, ref_background = NULL,
                                        databases = c("GO:BP", "GO:MF", "KEGG", "REAC", "WP", "CORUM")) {
  require(gprofiler2)
  gost(query = genes,
       organism = "hsapiens", ordered_query = ordered,
       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
       measure_underrepresentation = FALSE, evcodes = FALSE,
       user_threshold = p_thresh, correction_method = "g_SCS",
       domain_scope = "annotated", custom_bg = ref_background,
       numeric_ns = "", sources = databases, as_short_link = FALSE, highlight = TRUE)
}

ora_webgestalt <- function(genes, ref_background, p_thresh = p_threshold,
                           databases = c("pathway_KEGG", "pathway_Reactome", "pathway_Wikipathway_cancer",
                                         "geneontology_Biological_Process", "geneontology_Molecular_Function")) {
  require(WebGestaltR)
  WebGestaltR(enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = databases,
              interestGene = as.vector(genes), interestGeneType = "genesymbol",
              referenceGene = as.vector(ref_background), referenceGeneType = "genesymbol",
              fdrThr = p_thresh, isOutput = FALSE)
}

plot_terms <- function(ora, selected_sources = c("CORUM", "KEGG", "REAC", "WP", "GO:MF", "GO:BP"),
                       terms_per_source = 5, p_thresh = p_threshold, string_trunc = 50) {
  require(forcats)
  require(stringr)

  if (is.null(ora)) return(NULL)

  ora$result %>%
    filter(p_value < p_thresh) %>%
    filter(source %in% selected_sources) %>%
    group_by(source) %>%
    slice_min(p_value, n = terms_per_source, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(`-log10(p)` = -log10(p_value),
           term_name = fct_reorder(str_trunc(term_name, string_trunc), p_value, .desc = TRUE)) %>%
    vertical_bar_chart(term_name, `-log10(p)`, color_col = source, color_discrete = TRUE,
                       value_range = c(1, max(.$`-log10(p)`)),
                       line_intercept = 0, bar_label_shift = 0.18, break_steps = 2,
                       category_lab = "Enriched Term", value_lab = "-log10(p)") +
    facet_wrap(~source, scales = "free_y", ncol = 1) +
    theme(legend.position = "none")
}

plot_terms_compact <- function(ora, selected_sources = c("CORUM", "KEGG", "REAC", "WP", "GO:MF", "GO:BP"),
                               terms_per_source = 5, p_thresh = p_threshold, string_trunc = 50,
                               custom_color = NULL, luminance_shift = +1) {
  require(forcats)
  require(stringr)
  require(ggplot2)
  require(shadowtext)
  require(gt)

  if (is.null(ora)) return(NULL)
  if (!is.null(custom_color)) custom_color <- adjust_luminance(custom_color, luminance_shift)

  df <- ora$result %>%
    filter(p_value < p_thresh) %>%
    filter(source %in% selected_sources) %>%
    group_by(source) %>%
    slice_min(p_value, n = terms_per_source, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(`-log10(p)` = -log10(p_value),
           term_name = fct_reorder(str_trunc(term_name, string_trunc), p_value, .desc = TRUE))

  df %>%
    ggplot() +
    aes(x = term_name, y = `-log10(p)`,
        label = term_name) +
    { if (!is.null(custom_color)) geom_bar(fill = custom_color, stat = "identity")
      else geom_bar(aes(fill = source), stat = "identity") } +
    geom_text(color = default_color, y = 1, hjust = 0) +
    scale_y_continuous(breaks = seq(1, max(df$`-log10(p)`), 2)) +
    scale_fill_viridis(option = "G", direction = -1, begin = 0.2, end = 0.8, discrete = TRUE) +
    labs(x = "Enriched Term", y = "-log10(p)") +
    coord_flip(ylim = c(1, max(df$`-log10(p)`))) +
    facet_wrap(~source, scales = "free_y", ncol = 1) +
    theme(legend.position = "none",
          axis.text.y = element_blank())
}

# === Mutual Exclusivity ===
df2contingency <- function(df, x_col, y_col, count_col = n) {
  df %>%
    select({ { x_col } }, { { y_col } }, { { count_col } }) %>%
    drop_na() %>%
    arrange({ { x_col } }, { { y_col } }) %>%
    pivot_wider(names_from = quo_name(enquo(y_col)), values_from = quo_name(enquo(count_col))) %>%
    tibble::column_to_rownames(quo_name(enquo(x_col))) %>%
    as.matrix()
}

mutex_score <- function (contingency_matrix, ignore_bad_matrix = TRUE) {
  if (nrow(contingency_matrix) != 2 || ncol(contingency_matrix) != 2 || any(is.na(contingency_matrix))) {
    if (ignore_bad_matrix) return(NA)
    stop("Contingency matrix must have exactly 2 rows and columns and no NA values")
  }

  # Source: Girish, et al. “Oncogene-like Addiction to Aneuploidy in Human Cancers.” Preprint. Cancer Biology, January 10, 2023. https://doi.org/10.1101/2023.01.09.523344.
  test <- fisher.test(contingency_matrix)
  score <- abs(qnorm(test$p.value)) * sign(-log2(test$estimate[[1]]))
  return(list(MutEx = score, p = test$p.value, OR = test$estimate[[1]]))
}
