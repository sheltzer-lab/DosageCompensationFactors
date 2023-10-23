library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(limma)

normalize_samples <- function(df, sample_col, value_col, group_col, normalized_colname = "Normalized") {
  df_pre <- data.frame(df)

  df %>%
    assert_rows(col_concat, is_uniq, {{ sample_col }}, {{ group_col }}) %>%
    select({{ sample_col }}, {{ value_col }}, {{ group_col }}) %>%
    pivot_wider(names_from = {{ sample_col }}, values_from = {{ value_col }}) %>%
    tibble::column_to_rownames(var = quo_name(enquo(group_col))) %>%
    select(where(is.numeric)) %>%
    normalizeBetweenArrays(method = "cyclicloess", cyclic.method = "fast") %>%
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
            by = c("CellLine.CustomId", "Gene.Symbol", "Protein.Uniprot.Accession"),
            na_matches = "never")
  matched_set2 <- df_dataset2 %>%
    semi_join(y = df_dataset1,
            by = c("CellLine.CustomId", "Gene.Symbol", "Protein.Uniprot.Accession"),
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