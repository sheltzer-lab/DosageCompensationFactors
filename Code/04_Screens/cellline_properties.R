library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)
library(ggsignif)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Screens")
output_data_dir <- output_data_base_dir
reports_dir <- reports_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

cellline_buf_filtered_procan <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_procan.parquet"))
cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_procan.parquet"))
cellline_buf_filtered_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_filtered_depmap.parquet"))
cellline_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_depmap.parquet"))

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20230801.csv")) %>%
  rename(CellLine.Name = "model_name") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.Name = "CellLineName") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

df_procan <- df_model_procan %>%
  inner_join(y = cellline_buf_filtered_procan, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never")

df_depmap <- df_model_depmap %>%
  inner_join(y = cellline_buf_filtered_depmap, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never")


sorted_violin_plot <- function(df, x, y) {
  plot <- df %>%
    add_count(get(x)) %>%
    filter(n > 2) %>%
    group_by(get(x)) %>%
    mutate(Median = median({ { y } }),
           Label = paste0(get(x), " (n=", n, ")")) %>%
    ungroup() %>%
    arrange(Median) %>%
    mutate(Label = factor(Label, levels = unique(Label))) %>%
    violin_plot(Label, { { y } })
  plot <- plot +
    xlab(x)

  return(plot)
}

signif_violin_plot <- function(df, x, y, test = t.test, test.args = NULL) {
  df <- df %>%
    add_count({ { x } }) %>%
    filter(n > 2) %>%
    group_by({ { x } }) %>%
    mutate(Median = median({ { y } }),
           Label = paste0({ { x } }, " (n=", n, ")")) %>%
    ungroup() %>%
    arrange(Median) %>%
    mutate(Label = factor(Label, levels = unique(Label)))

  plot <- df %>%
    ggplot() +
    aes(x = Label, y = { { y } }) +
    geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75),
                color = "#4080DB") +
    geom_signif(
      comparisons = list(levels(df$Label)),
      map_signif_level = FALSE, tip_length = 0,
      test = test, test.args = test.args
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(as_name(enquo(x)))

  return(plot)
}

plot_categorical_properties <- function (df, cols_categorical) {
  plots <- list()
  for (col in cols_categorical) {
    plots[[col]] <- df %>%
      sorted_violin_plot(col, Buffering.CellLine.Ratio)
  }
  return(plots)
}

data_density <- function(df) {
  df %>%
    summarize_all(~sum(!is.na(.x)) / nrow(df)) %>%
    pivot_longer(everything(), names_to = "Column", values_to = "Density") %>%
    arrange(desc(Density))
}

data_density_procan <- df_procan %>%
  data_density()

data_density_depmap <- df_procan %>%
  data_density()

cols_procan <- c("tissue_status", "cancer_type", "msi_status",
                 "smoking_status", "gender", "ethnicity")

cols_depmap <- c("PrimaryOrMetastasis", "OncotreeSubtype", "Sex")

violoin_plots_procan <- df_procan %>%
  plot_categorical_properties(cols_procan)

violoin_plots_depmap <- df_depmap %>%
  plot_categorical_properties(cols_depmap)

df_procan %>%
  signif_violin_plot(tissue_status, Buffering.CellLine.Ratio,
                     test.args = list(var.equal = TRUE))

df_procan %>%
  signif_violin_plot(msi_status, Buffering.CellLine.Ratio)

df_procan %>%
  filter(gender != "Unknown") %>%
  signif_violin_plot(gender, Buffering.CellLine.Ratio,
                     test.args = list(var.equal = TRUE))

df_depmap %>%
  filter(Sex != "Unknown") %>%
  signif_violin_plot(Sex, Buffering.CellLine.Ratio,
                     test.args = list(var.equal = TRUE))

df_depmap %>%
  filter(PrimaryOrMetastasis != "NA") %>%
  signif_violin_plot(PrimaryOrMetastasis, Buffering.CellLine.Ratio,
                     test.args = list(var.equal = TRUE))


df_procan %>%
  scatter_plot_regression(ploidy_wes, Buffering.CellLine.Ratio,  Buffering.CellLine.Ratio ~ ploidy_wes)

df_procan %>%
  scatter_plot_regression(mutational_burden, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ mutational_burden,
                          label_coords = c(500, -0.4))

df_procan %>%
  scatter_plot_regression(age_at_sampling, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ age_at_sampling,
                          label_coords = c(50, -0.4))

df_procan %>%
  scatter_plot_regression(age_at_sampling, ploidy_wes, ploidy_wes ~ age_at_sampling)


df_depmap %>%
  scatter_plot_regression(Age, Buffering.CellLine.Ratio, Buffering.CellLine.Ratio ~ Age,
                          label_coords = c(50, -0.4))
