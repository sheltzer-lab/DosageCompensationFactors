library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(limma)

normalize_samples <- function(df, sample_col, value_col, group_col,
                              normalized_colname = "Normalized", batch_col = NULL) {
  require(limma)
  df_pre <- data.frame(df)
  batch <- NULL

  if (!quo_is_null(enquo(batch_col)))
    batch <- df %>% distinct({{ sample_col }}, {{ batch_col }}) %>% pull({{ batch_col }}, name = {{ sample_col }})

  df %>%
    assert_rows(col_concat, is_uniq, {{ sample_col }}, {{ group_col }}) %>%
    select({{ sample_col }}, {{ value_col }}, {{ group_col }}) %>%
    pivot_wider(names_from = {{ sample_col }}, values_from = {{ value_col }}) %>%
    tibble::column_to_rownames(var = quo_name(enquo(group_col))) %>%
    select(where(is.numeric)) %>%
    normalizeBetweenArrays(method = "cyclicloess", cyclic.method = "fast") %>%
  { if (!is.null(batch)) removeBatchEffect(., batch) else . } %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = quo_name(enquo(group_col))) %>%
    pivot_longer(everything() & !{{ group_col }},
                 names_to = quo_name(enquo(sample_col)), values_to = normalized_colname) %>%
    inner_join(y = df_pre,
               by = c(quo_name(enquo(sample_col)), quo_name(enquo(group_col))),
               relationship = "one-to-one", na_matches = "never")
}

remove_noisefloor <- function(df, value_col, percentile_cutoff = 0.0001) {
  # ToDo: Determine Noise Floor by Signal-To-Noise-Ratio
  noise_floor <- quantile(df[[quo_name(enquo(value_col))]],
                          probs = percentile_cutoff, na.rm = TRUE)[[1]]

  df %>%
    mutate(!!enquo(value_col) := replace({{ value_col }}, {{ value_col }} < noise_floor, NA))
}

# Combine two datasets by matching cell lines and genes for each matched cell line
match_datasets <- function(df_dataset1, df_dataset2) {
  matched_set1 <- df_dataset1 %>%
    semi_join(y = df_dataset2,
            by = c("Model.ID", "Gene.Symbol", "Protein.Uniprot.Accession"),
            na_matches = "never")
  matched_set2 <- df_dataset2 %>%
    semi_join(y = df_dataset1,
            by = c("Model.ID", "Gene.Symbol", "Protein.Uniprot.Accession"),
            na_matches = "never")

  return(bind_rows(matched_set1, matched_set2))
}

plot_expr_dist <- function(df) {
  df %>%
    ggplot() +
    geom_density(na.rm = TRUE, aes(Protein.Expression.Log2, color = "Non-Normalized"), show.legend = TRUE) +
    geom_density(na.rm = TRUE, aes(Protein.Expression.Normalized, color = "Normalized"), show.legend = TRUE) +
    xlab("Protein Expression")
}

# Multivariate
shuffle_rows <- function(df) {
  df[sample(nrow(df), replace = FALSE),]
}

# ToDo: Explore further imputation methods
impute_na <- function(df, factor_cols = NULL) {
  imputation_func <- \(x) replace_na(x, median(x, na.rm = TRUE))

  if (is.null(factor_cols)) {
    df %>%
      mutate_if(is.numeric, imputation_func) %>%
      return()
  }

  df %>%
    mutate_at(factor_cols, imputation_func) %>%
    return()
}

impute_na_long <- function(df, group_col, value_col) {
  df %>%
    group_by({ { group_col } }) %>%
    mutate(!!target_group_col := replace_na({ { value_col } }, median({ { value_col } }, na.rm = TRUE)),
           IsImputed = is.na({ { value_col } })) %>%
    ungroup()
}

balanced_sample <- function(df, class_col, n_per_class) {
  df %>%
    group_by({ { class_col } }) %>%
    slice_sample(n = n_per_class) %>%
    ungroup()
}

# TODO: Explore further techniques
# Oversampling, Undersampling, slice_sample(weight_by=..., replace=TRUE)
rebalance_binary <- function(df, class_col, target_balance = 0.5) {
  total_size <- nrow(df)
  new_samples <- df %>%
    count({ { class_col } }) %>%
    mutate(Ratio = n/total_size) %>%
    mutate(Samples = c(ceiling((target_balance * min(n)) / (1 - target_balance)), min(n))) %>%
    select(-n, -Ratio)

  df %>%
    left_join(y = new_samples, by = quo_name(enquo(class_col))) %>%
    group_by({ { class_col } }) %>%
    group_map(~slice_sample(.x, n = unique(.x$Samples)), .keep = TRUE) %>%
    bind_rows() %>%
    select(-Samples)
}

clean_data <- function (dataset, buffering_class_col, factor_cols = dc_factor_cols, id_col = "UniqueId") {
  dataset %>%
    filter({ { buffering_class_col } } == "Buffered" | { { buffering_class_col } } == "Scaling") %>%
    mutate(Buffered = factor({ { buffering_class_col } }, levels = c("Scaling", "Buffered"))) %>%
    mutate_if(is.numeric, \(x) as.numeric(x)) %>%
    drop_na(Buffered) %>%
    tibble::column_to_rownames(id_col) %>%
    select(Buffered, all_of(factor_cols)) %>%
    janitor::clean_names()
}

normalize_min_max <- function(x, ...) {
    return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}

normalize_features <- function (df, method = "min-max", factor_cols = dc_factor_cols) {
  if (method == "min-max") {
    df %>%
      mutate_at(factor_cols, \(x) normalize_min_max(x, na.rm = TRUE))
  } else if (method == "z-score") {
    df %>%
      mutate_at(factor_cols, ~(scale(.) %>% as.vector))
  } else {
    warning("Invalid normalization method! Data frame passed as is. Supported methods: min-max, z-score")
    return(df)
  }
}

split_train_test <- function(df_prep, training_set_ratio) {
  df_prep <- df_prep  %>%
    mutate(id = row_number())

  training_set <- df_prep %>%
    slice_sample(by = buffered, prop = training_set_ratio)

  test_set <- df_prep %>%
    anti_join(training_set, by = 'id')

  return(list(training = training_set %>% select(-id),
              test = test_set %>% select(-id)))
}