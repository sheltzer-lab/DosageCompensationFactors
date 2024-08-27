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
source(here("Code", "analysis.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

tumor_types <- get_tumor_types()

df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))

model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
model_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))

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

df_depmap <- model_buf_depmap %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

df_procan <- model_buf_procan %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

df_cptac <- model_buf_cptac %>%
  inner_join(y = df_model_cptac, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  rename(OncotreeCode = Model.CancerType) %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

df_agg <- model_buf_agg %>%
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

# === Cell Lines Panel ===
waterfall_procan <- model_buf_procan %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name)
waterfall_depmap <- model_buf_depmap %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name)

model_buf_celllines <- model_buf_procan %>%
  bind_rows(model_buf_depmap) %>%
  arrange(CellLine.Name) %>%
  pivot_wider(names_from = "Dataset", values_from = "Model.Buffering.Ratio", id_cols = "Model.ID")

cellline_pearson_gene <- cor.test(x = model_buf_celllines$ProCan,
                                  y = model_buf_celllines$DepMap,
                                  method = "pearson")

plot_agg_top <- model_buf_agg %>%
  slice_max(Model.Buffering.MeanNormRank, n = 10) %>%
  vertical_bar_chart(CellLine.Name, Model.Buffering.MeanNormRank,
                     default_fill_color = head(bidirectional_color_pal, n = 1),
                     text_color = "black", bar_label_shift = 0.1, break_steps = 0.25,
                     value_range = c(0,1), line_intercept = 0.5, value_lab = "Mean Normalized Rank")

plot_agg_bot <- model_buf_agg %>%
  slice_min(Model.Buffering.MeanNormRank, n = 10) %>%
  vertical_bar_chart(CellLine.Name, Model.Buffering.MeanNormRank,
                     default_fill_color = tail(bidirectional_color_pal, n = 1),
                     text_color = "black", bar_label_shift = 0.1, break_steps = 0.25,
                     value_range = c(0,1), line_intercept = 0.5, value_lab = "Mean Normalized Rank")

plot_bracket <- plot_corr_bracket(cellline_pearson_gene)
plot_stack1 <- cowplot::plot_grid(waterfall_procan, waterfall_depmap + ylab(NULL),
                                  nrow = 1, ncol = 2, align = "hv", axis = "tblr", labels = c("ProCan", "DepMap"),
                                  label_y = 0.98, label_x = 0.1, rel_widths = c(1, 1))
panel_celllines <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                      nrow = 2, ncol = 1,
                                      rel_heights = c(0.1, 1))

panel_agg <- cowplot::plot_grid(plot_agg_top + ylab(NULL) + cowplot::theme_minimal_vgrid()
                                 + theme(axis.text.x = element_blank()),
                                plot_agg_bot + cowplot::theme_minimal_vgrid(),
                                nrow = 2, ncol = 1, align = "v", axis = "lr")

panel_celllines_agg <- cowplot::plot_grid(panel_celllines, panel_agg,
                                          nrow = 1, ncol = 2, labels = c("A", "B"),
                                          rel_widths = c(1.75, 1))

# === Cancer Types Panel ===
df_cancer_heatmap <- bind_rows(df_depmap, df_procan, df_cptac) %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2))) %>%
  group_by(Dataset, OncotreeCode) %>%
  summarize(Mean.BR = mean(Model.Buffering.Ratio, na.rm = TRUE),
            Mean.AS = mean(CellLine.AneuploidyScore, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(OncotreeCode = fct_reorder(OncotreeCode, Mean.BR)) %>%
  drop_na(OncotreeCode) %>%
  complete(Dataset, OncotreeCode) %>%
  mutate_all(~if_else(is.nan(.), NA, .))

cancer_heatmap_br <- df_cancer_heatmap %>%
  ggplot() +
  aes(x = Dataset, y = OncotreeCode, fill = Mean.BR) +
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
        legend.box = "vertical",
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5, position = "top"))

cancer_heatmap_as <- df_cancer_heatmap %>%
  group_by(OncotreeCode) %>%
  summarize(Mean.AS = mean(Mean.AS, na.rm = TRUE),
            Dataset = "Mean Aneuploidy") %>%
  ggplot() +
  aes(x = Dataset, y = OncotreeCode, fill = Mean.AS) +
  geom_tile() +
  scale_fill_viridis_c(na.value = color_palettes$Missing, option = color_palettes$AneuploidyScore) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

panel_types <- cowplot::plot_grid(cancer_heatmap_br, cancer_heatmap_as,
                                  nrow = 1, ncol = 2, align = "h", axis = "lr",
                                  rel_widths = c(1, 0.2))

# == Lymphoma Leukemia Panel ===
leukemia_codes <- c("MNM", "LNM")

df_leuk <- df_depmap %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2))) %>%
  filter(OncotreeCode %in% leukemia_codes)

min_aneuploidy <- min(df_depmap$CellLine.AneuploidyScore)
max_aneuploidy <- max(df_depmap$CellLine.AneuploidyScore)
max_aneuploidy_leuk <- round(quantile(df_leuk$CellLine.AneuploidyScore, probs = 0.9))

leuk_plot <- df_depmap %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2))) %>%
  mutate(`Myeloid / Lymphoid` = OncotreeCode %in% leukemia_codes) %>%
  signif_beeswarm_plot(`Myeloid / Lymphoid`, Model.Buffering.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test, signif_label = \(p) paste0("p = ", formatC(p, format = "e", digits = 2)))

leuk_plot_low <- df_depmap %>%
  filter(CellLine.AneuploidyScore <= max_aneuploidy_leuk) %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2))) %>%
  mutate(`Myeloid / Lymphoid` = OncotreeCode %in% leukemia_codes) %>%
  signif_beeswarm_plot(`Myeloid / Lymphoid`, Model.Buffering.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test, signif_label = \(p) paste0("p = ", formatC(p, format = "e", digits = 2)))

leuk_poster <- leuk_plot +
  theme_light(base_size = 20) +
  theme(legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top") +
  guides(color = guide_colourbar(ticks.linewidth = 1)) +
  labs(color = "\nAneuploidy Score", y = "Mean Buffering Ratio") +
  scale_colour_viridis_c(option = "D", direction = 1,
                         limits = c(min_aneuploidy, max_aneuploidy),
                         breaks = c(min_aneuploidy, max_aneuploidy))

leuk_poster_low <- leuk_plot_low +
  theme_light(base_size = 20) +
  theme(legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top") +
  guides(color = guide_colourbar(ticks.linewidth = 1)) +
  labs(color = "Aneuploidy Score\n(Limited Range)", y = NULL) +
  scale_colour_viridis_c(option = "D", direction = 1,
                         limits = c(min_aneuploidy, max_aneuploidy),
                         breaks = c(min_aneuploidy, max_aneuploidy_leuk[["90%"]], max_aneuploidy),
                         labels = c(min_aneuploidy, max_aneuploidy_leuk[["90%"]], max_aneuploidy))

panel_leuk <- cowplot::plot_grid(leuk_poster,
                                 leuk_poster_low,
                                 nrow = 1, ncol = 2)

# === WGD & Aneuploidy Score ===
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
panel_wgd <- cowplot::insert_yaxis_grob(wgd_panel_pre, br_density_panel, position = "right")

cowplot::ggdraw(panel_wgd)

# === Proliferation ===
df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet")) %>%
  select(Model.ID, CellLine.GrowthRatio)

df_prolif <- df_agg %>%
  select(Model.ID, Model.Buffering.MeanNormRank, CellLine.WGD) %>%
  inner_join(y = df_growth, by = "Model.ID") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"))

cor_prolif <- cor.test(df_prolif$Model.Buffering.MeanNormRank,
                       df_prolif$CellLine.GrowthRatio, method = "spearman")

panel_prolif_base <- df_prolif %>%
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

# === Proteotoxic Stress Panel ===
library(msigdbr)
hallmark_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>%
  rename(Gene.Symbol = "gene_symbol")

gene_sets <- hallmark_gene_set %>%
  group_by(gs_name) %>%
  group_map(~list(.x$Gene.Symbol)) %>%
  purrr::list_flatten()
names(gene_sets) <- unique(hallmark_gene_set$gs_name)

gsea_cptac <- expr_buf_cptac %>%
  ssGSEA(gene_sets, Protein.Expression.Normalized, Model.ID) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model.ID") %>%
  inner_join(y = model_buf_cptac, by = "Model.ID")

panel_proteotox <- gsea_cptac %>%
  ggplot() +
  aes(x = Model.Buffering.Ratio, y = HALLMARK_UNFOLDED_PROTEIN_RESPONSE) +
  geom_point(size = 2) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0.2, y = 0.3, hjust = 0, size = 5,
           label = paste0(print_corr(cor_as_test$estimate), ", p = ", formatC(cor_as_test$p.value, format = "e", digits = 2))) +
  labs(x = "Model Buffering Ratio", y = "Unfolded Protein Response") +
  theme(legend.position = "none")

# === Combine Panels into Figure ===
figure2_sub1 <- cowplot::plot_grid(panel_wgd, panel_leuk,
                                   panel_prolif_base, panel_proteotox,
                                   nrow = 2, ncol = 2, labels = c("D", "E", "F", "G"))

figure2_sub2 <- cowplot::plot_grid(panel_types, figure2_sub1,
                                   nrow = 1, ncol = 2, rel_widths = c(0.2, 1),
                                   labels = c("C", ""))

figure2 <- cowplot::plot_grid(panel_celllines_agg, figure2_sub2,
                              nrow = 2, ncol = 1, rel_heights = c(0.5, 1))

cairo_pdf(here(plots_dir, "figure02.pdf"), width = 15, height = 20)
figure2
dev.off()
