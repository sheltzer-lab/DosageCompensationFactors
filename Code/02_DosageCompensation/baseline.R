library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(readxl)
library(assertr)
library(ggplot2)
library(limma)

here::i_am("DosageCompensationFactors.Rproj")

calculate_weights <- function(distances) {
  weights <- (1/(1+distances)) / sum(1/(1+distances), na.rm = TRUE)
  return(weights)
}

calculate_baseline <- function(df, gene_col, chr_arm_cna_col, value_col, summ_func = mean,
                               target_colname = "Baseline", distance_col = NULL, weighted = FALSE) {
  baseline_expr <- df %>%
    select({ { gene_col } }, { { chr_arm_cna_col } }, { { value_col } }, { { distance_col } }) %>%
    filter({ { chr_arm_cna_col } } == 0) %>%
    group_by({ { gene_col } }) %>%
    mutate(Weights = calculate_weights({ { distance_col } })) %>%
    summarize(!!target_colname := if_else(weighted,
                                          sum({ { value_col } } * Weights, na.rm = TRUE),
                                          summ_func({ { value_col } }, na.rm = TRUE))) %>%
    ungroup()

  df %>%
    inner_join(y = baseline_expr, by = quo_name(enquo(gene_col)),
               unmatched = "error", na_matches = "never", relationship = "many-to-one")
}

calculate_median_baseline <- function(df, gene_col, value_col, target_colname = "Baseline") {
  df %>%
    group_by({ { gene_col } }) %>%
    mutate(!!target_colname := median({ { value_col } }, na.rm = TRUE)) %>%
    ungroup()
}
