library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Screens", "Proliferation")
output_data_dir <- output_data_base_dir
reports_dir <- here(reports_base_dir, "Proliferation")

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet"))
df_growth_chronos <- read_parquet(here(output_data_dir, "cellline_growth_chronos.parquet")) %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)
cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
cellline_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))
cellline_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
# Load copy number dataset to obtain metadata
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet"))

# ToDo: Use Cell Line ID instead of cell line name

merge_datasets <- function(cellline_dc_dataset, copy_number_dataset, prolif_dataset) {
  cellline_cn_metadata <- copy_number_dataset %>%
    distinct(CellLine.Name, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

  df_prolif <- cellline_dc_dataset %>%
    inner_join(y = prolif_dataset, by = "CellLine.Name",
               relationship = "one-to-one", na_matches = "never") %>%
    inner_join(y = cellline_cn_metadata, by = "CellLine.Name",
               relationship = "one-to-one", na_matches = "never") %>%
    mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"))

  return(df_prolif)
}

create_plots <- function (df_prolif, color_col = NULL, dataset_name = NULL, condition_name = NULL) {
  filename_suffix <- tolower(paste(c(dataset_name, condition_name), collapse = "_"))
  title <- paste(c(dataset_name, condition_name), collapse = " ")

  # Scatter plot
  prolif_plot <- df_prolif %>%
    scatter_plot_reg_corr(Model.Buffering.Ratio, CellLine.GrowthRatio,
                          point_size = 2, color_col = { { color_col } },
                          title_prefix = title) %>%
    save_plot(paste0("dosage_compensation_proliferation_", filename_suffix, ".png"))

  ## Correlation Report
  report_file <- here(reports_dir, paste0("dosage_compensation_proliferation_correlation_", filename_suffix, ".txt"))
  dc_prolif_corr <- cor.test(x = df_prolif$Model.Buffering.Ratio, y = df_prolif$CellLine.GrowthRatio,
                             method = "spearman")
  dc_prolif_corr %>%
    capture.output() %>%
    writeLines(con = report_file)

  # Statistical comparison between cell line groups
  dc_growth_violin <- df_prolif %>%
    split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
    signif_violin_plot(Model.Buffering.Group, CellLine.GrowthRatio,
                       test = wilcox.test, title = title) %>%
    save_plot(paste0("dosage_compensation_proliferation_test_", filename_suffix, ".png"))

  ## Ratio between cell line groups
  df_median <- df_prolif %>%
    split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
    group_by(Model.Buffering.Group) %>%
    summarize(MedianGrowthRatio = median(CellLine.GrowthRatio, na.rm = TRUE)) %>%
    pivot_wider(names_from = Model.Buffering.Group, values_from = MedianGrowthRatio)

  growth_change <- df_median$High / df_median$Low

  return(list(scatter_plot = prolif_plot, corr = dc_prolif_corr,
              violin_plot = dc_growth_violin, growth_change = growth_change))
}

proliferation_analysis <- function (cellline_dc_dataset, copy_number_dataset, prolif_dataset, dataset_name = NULL) {
  df_prolif <- merge_datasets(cellline_dc_dataset, copy_number_dataset, prolif_dataset)

  # All cells
  results_all <- df_prolif %>%
    create_plots(dataset_name, color_col = WGD)

  # WGD-only cells
  results_wgd <- df_prolif %>%
    filter(CellLine.WGD > 0) %>%
    create_plots(dataset_name, color_col = CellLine.AneuploidyScore, condition_name = "WGD")

  # Non-WGD cells
  results_nonwgd <- df_prolif %>%
    filter(CellLine.WGD == 0) %>%
    create_plots(dataset_name, color_col = CellLine.AneuploidyScore, condition_name = "Non-WGD")

  return(list(All = results_all, WGD = results_wgd, NonWGD = results_nonwgd))
}

results_procan <- cellline_buf_procan %>%
  proliferation_analysis(copy_number, df_growth, dataset_name = "ProCan")

results_depmap <- cellline_buf_depmap %>%
  proliferation_analysis(copy_number, df_growth, dataset_name = "DepMap")

results_agg <- cellline_buf_agg %>%
  mutate(Model.Buffering.Ratio = Model.Buffering.MeanNormRank) %>%
  proliferation_analysis(copy_number, df_growth, dataset_name = "Aggregated")


# === Combine Plots for publishing ===
plot_publish <- cowplot::plot_grid(results_procan$All$scatter_plot, results_procan$NonWGD$scatter_plot,
                                   results_procan$All$violin_plot, results_procan$NonWGD$violin_plot,
                                   ncol = 2, nrow = 2, labels = c("A", "", "B", ""))

cairo_pdf(here(plots_dir, "proliferation_publish.pdf"), height = 10, width = 12)
plot_publish
dev.off()

## Poster
df_prolif <- merge_datasets(cellline_buf_procan, copy_number, df_growth) %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"))

growth_poster <- (df_prolif %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, CellLine.GrowthRatio,
                          point_size = 2, color_col = WGD,
                          title_prefix = "ProCan")) +
  theme_light(base_size = 20) +
  labs(x = "Mean Buffering Ratio", y = "Cell Line Growth Rate") +
  scale_color_manual(values = two_class_color_pal) +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        legend.background = element_blank())

growth_poster_nowgd <- (df_prolif %>%
  filter(WGD == "Non-WGD") %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, CellLine.GrowthRatio,
                          point_size = 2, color_col = WGD,
                          title_prefix = "ProCan, Non-WGD")) +
  theme_light(base_size = 20) +
  labs(x = "Mean Buffering Ratio", y = "Cell Line Growth Rate") +
  scale_color_manual(values = two_class_color_pal) +
  theme(legend.position = "none")


plot_poster <- cowplot::plot_grid(growth_poster, growth_poster_nowgd,
                                   ncol = 1, nrow = 2, align = "h", axis = "tblr")

cairo_pdf(here(plots_dir, "proliferation_poster.pdf"), height = 12)
plot_poster
dev.off()
