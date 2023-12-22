library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Screens", "Proliferation")
output_data_dir <- output_data_base_dir
reports_dir <- reports_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet"))
cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_procan.parquet"))
cellline_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(CellLine.Name, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

# ToDo: Use Cell Line ID instead of cell line name
df_prolif <- cellline_buf_procan %>%
  inner_join(y = df_growth, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never")

dc_growth_cor <- cor.test(x = df_prolif$Buffering.CellLine.Ratio,
                          y = df_prolif$CellLine.GrowthRatio,
                          method = "spearman")

cat(capture.output(dc_growth_cor), file = here(reports_dir, "dosage_compensation_proliferation_correlation.txt"),
    append = FALSE, sep = "\n")

dc_growth_reg <- df_prolif %>%
  scatter_plot_reg_corr(Buffering.CellLine.Ratio, CellLine.GrowthRatio,
                          point_size = 1.5, label_coords = c(0, 0),
                          title_prefix = "ProCan") %>%
  save_plot("dosage_compensation_proliferation_procan.png")

# Exclude whole-genome doubling
df_prolif_nowgd <- df_prolif %>%
  filter(CellLine.WGD == 0)

dc_growth_cor_nowgd <- cor.test(x = df_prolif_nowgd$Buffering.CellLine.Ratio,
                                y = df_prolif_nowgd$CellLine.GrowthRatio,
                                method = "spearman")

cat(capture.output(dc_growth_cor_nowgd),
    file = here(reports_dir, "dosage_compensation_proliferation_correlation_no-wgd.txt"),
    append = FALSE, sep = "\n")

dc_growth_reg_nowgd <- df_prolif_nowgd %>%
  scatter_plot_reg_corr(Buffering.CellLine.Ratio, CellLine.GrowthRatio,
                          point_size = 1.5, label_coords = c(0, 0),
                          title_prefix = "ProCan, No-WGD") %>%
  save_plot("dosage_compensation_proliferation_procan_nowgd.png")

# === Combine Plots for publishing ===
plot_publish <- cowplot::plot_grid(dc_growth_reg, dc_growth_reg_nowgd)

cairo_pdf(here(plots_dir, "proliferation_publish.pdf"), height = 5, width = 12)
plot_publish
dev.off()