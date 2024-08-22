library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(viridisLite)
library(mskcc.oncotree)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

tumor_types <- get_tumor_types()

df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))

cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
cellline_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))
cellline_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))
cellline_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20240110.csv")) %>%
  rename(CellLine.SangerModelId = "model_id") %>%
  inner_join(y = df_celllines, by = "CellLine.SangerModelId") %>%
  mutate(OncotreeCode = nci_to_oncotree(cancer_type_ncit_id, expand = TRUE)$oncotree_code) %>%
  select(-CellLine.Name)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId") %>%
  select(-CellLine.Name)

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

intersect(unique(df_model_depmap$OncotreeCode), unique(df_model_cptac$Model.CancerType))
intersect(unique(tumor_types$oncotree_code), unique(df_model_cptac$Model.CancerType))

df_depmap <- cellline_buf_depmap %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

df_procan <- cellline_buf_procan %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

df_cptac <- cellline_buf_cptac %>%
  inner_join(y = df_model_cptac, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  rename(OncotreeCode = Model.CancerType) %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

df_agg <- cellline_buf_agg %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

get_oncotree_parent <- function(df_tumor_types = NULL, start_code = NULL, target_level = 3) {
  require(mskcc.oncotree)

  if (is.null(start_code) || is.na(start_code)) return(NA)

  if (is.null(df_tumor_types))
    df_tumor_types <- mskcc.oncotree::get_tumor_types()

  current <- df_tumor_types %>%
    filter(oncotree_code == start_code)

  if (nrow(current) == 0) return(NA)

  if (current$level <= target_level) return(current$oncotree_code)
  else get_oncotree_parent(df_tumor_types, current$parent, target_level)
}

# === Cancer Types Panel ===

# TODO: Control for aneuploidy score
df_cancer_heatmap <- bind_rows(df_depmap, df_procan, df_cptac) %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 1))) %>%
  group_by(Dataset, OncotreeCode) %>%
  summarize(Mean.BR = mean(Model.Buffering.Ratio, na.rm = TRUE),
            Mean.AS = mean(CellLine.AneuploidyScore, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(OncotreeCode = fct_reorder(OncotreeCode, Mean.BR)) %>%
  drop_na(OncotreeCode) %>%
  complete(Dataset, OncotreeCode)

cancer_heatmap_br <- df_cancer_heatmap %>%
  ggplot() +
  aes(y = Dataset, x = OncotreeCode, fill = Mean.BR) +
  geom_tile(aes(color = Mean.AS), alpha = 0) +
  geom_tile() +
  scale_color_viridis_c(na.value = color_palettes$Missing, option = color_palettes$AneuploidyScore) +
  scale_fill_viridis_c(na.value = color_palettes$Missing, option = color_palettes$BufferingRatio) +
  scale_x_discrete(position = "top") +
  theme_void() +
  labs(x = NULL, fill = "Mean Buffering Ratio", color = "Mean Aneuploidy") +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_colourbar(order = 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        legend.key.size = unit(16, "points"),
        legend.box = "horizontal")

cancer_heatmap_as <- df_cancer_heatmap %>%
  group_by(OncotreeCode) %>%
  summarize(Mean.AS = mean(Mean.AS, na.rm = TRUE),
            Dataset = "Mean Aneuploidy") %>%
  ggplot() +
  aes(y = Dataset, x = OncotreeCode, fill = Mean.AS) +
  geom_tile() +
  scale_fill_viridis_c(na.value = color_palettes$missing, option = color_palettes$AneuploidyScore) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

cowplot::plot_grid(cancer_heatmap_br, cancer_heatmap_as,
                   nrow = 2, ncol = 1, align = "v", axis = "tb",
                   rel_heights = c(1, 0.2))

# === WGD & Aneuploidy Score ===
print_corr(cor_as_test$estimate)

df_procan_wgd <- df_procan %>%
  distinct(Model.ID, Model.Buffering.Ratio, CellLine.WGD, CellLine.AneuploidyScore) %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"))

cor_as_test <- cor.test(df_procan$CellLine.AneuploidyScore,
                        df_procan$Model.Buffering.Ratio, method = "spearman")

wgd_base_panel <- df_procan_wgd %>%
  ggplot() +
  aes(x = CellLine.AneuploidyScore, y = Model.Buffering.Ratio, color = WGD) +
  geom_point(size = 2) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0, y = 1.8, hjust = 0, size = 5,
           label = paste0(print_corr(cor_as_test$estimate), ", p = ", formatC(cor_as_test$p.value, format = "e", digits = 2))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "Aneuploidy Score", y = "Model Buffering Ratio") +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.90),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size))

wgd_br_test <- wilcox.test(Model.Buffering.Ratio ~ WGD, data = df_procan_wgd)
wgd_as_test <- wilcox.test(CellLine.AneuploidyScore ~ WGD, data = df_procan_wgd)

wgd_median <- df_procan_wgd %>%
  group_by(WGD) %>%
  summarize_if(is.numeric, median)

as_median_range <- as.matrix(dist(wgd_median$CellLine.AneuploidyScore))[1,2]
br_median_range <- as.matrix(dist(wgd_median$Model.Buffering.Ratio))[1,2]

as_density_panel <- df_procan_wgd %>%
  ggplot() +
  aes(x = CellLine.AneuploidyScore, color = WGD, fill = WGD) +
  geom_density(alpha = 1 / 4) +
  geom_vline(data = wgd_median, aes(xintercept = CellLine.AneuploidyScore, color = WGD), linetype = "dashed") +
  geom_label(data = data.frame(x = mean(wgd_median$CellLine.AneuploidyScore), y = 0.08), aes(x = x, y = y),
             label = paste0("p = ", formatC(wgd_as_test$p.value, format = "e", digits = 2)),
             label.size = NA, fill = "white", color = "black", size = 4, alpha = 0.5, vjust = 0) +
  geom_segment(aes(x = wgd_median$CellLine.AneuploidyScore[1] + as_median_range * 0.1,
                   xend = wgd_median$CellLine.AneuploidyScore[2] - as_median_range * 0.1,
                   y = 0.08, yend = 0.08), color = "black", linewidth = 0.3) +
  scale_color_manual(values = color_palettes$WGD) +
  ylim(c(0, 0.12)) +
  theme_void() +
  theme(legend.position = "none")

br_density_panel <- df_procan_wgd %>%
  ggplot() +
  aes(x = Model.Buffering.Ratio, color = WGD, fill = WGD) +
  geom_density(alpha = 1 / 4) +
  geom_vline(data = wgd_median, aes(xintercept = Model.Buffering.Ratio, color = WGD), linetype = "dashed") +
  geom_label(data = data.frame(x = mean(wgd_median$Model.Buffering.Ratio), y = 3), aes(x = x, y = y),
             label = paste0("p = ", formatC(wgd_br_test$p.value, format = "e", digits = 2)),
             label.size = NA, fill = "white", color = "black", size = 4, alpha = 0.5, angle = -90, vjust = 0) +
  geom_segment(aes(x = wgd_median$Model.Buffering.Ratio[1] + br_median_range * 0.1,
                   xend = wgd_median$Model.Buffering.Ratio[2] - br_median_range * 0.1,
                   y = 3, yend = 3), color = "black", linewidth = 0.2) +
  scale_color_manual(values = color_palettes$WGD) +
  ylim(c(0, 4.5)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip()

wgd_panel_pre <- cowplot::insert_xaxis_grob(wgd_base_panel, as_density_panel, position = "top")
wgd_panel <- cowplot::insert_yaxis_grob(wgd_panel_pre, br_density_panel, position = "right")

cowplot::ggdraw(wgd_panel)

# === Proliferation ===
df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet")) %>%
  select(Model.ID, CellLine.GrowthRatio)

df_prolif <- df_agg %>%
  select(Model.ID, Model.Buffering.MeanNormRank, CellLine.WGD) %>%
  inner_join(y = df_growth, by = "Model.ID") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"))

cor_prolif <- cor.test(df_prolif$Model.Buffering.MeanNormRank,
                       df_prolif$CellLine.GrowthRatio, method = "spearman")

prolif_base_panel <- df_prolif %>%
  ggplot() +
  aes(x = Model.Buffering.MeanNormRank, y = CellLine.GrowthRatio, color = WGD) +
  geom_point(size = 2) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0, y = 15, hjust = 0, size = 5,
           label = paste0(print_corr(cor_prolif$estimate), ", ", print_signif(cor_prolif$p.value))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "Mean Normalized Buffering Ranks", y = "Growth Ratio (Day4/Day1)") +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.90),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size))

# TODO: Mention adjusted BR

# TODO: Why are more samples in df_agg than in df_procan and df_depmap?