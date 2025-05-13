library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(viridisLite)
library(openxlsx)
library(mskcc.oncotree)
library(survival)
library(survminer)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "publication.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Publication")
tables_dir <- here(tables_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

tumor_types <- get_tumor_types()

df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))

model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
model_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

tp53_mut <- read_parquet(here(output_data_dir, "damaging_mutations_depmap.parquet")) %>%
  filter(Gene.Symbol == "TP53") %>%
  select(-Gene.Symbol) %>%
  rename(TP53.Mutated = MutationStatus)

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20240110.csv")) %>%
  rename(CellLine.SangerModelId = "model_id") %>%
  inner_join(y = df_celllines, by = "CellLine.SangerModelId") %>%
  mutate(OncotreeCode = nci_to_oncotree(cancer_type_ncit_id, expand = TRUE)$oncotree_code) %>%
  select(-CellLine.Name)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId") %>%
  select(-CellLine.Name)

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet')) %>%
  # Oncotree uses different cancer type abbrieviations than CPTAC
  mutate(Model.CancerType = str_replace_all(Model.CancerType, c(HNSCC = "HNSC", LSCC = "LUSC", PDAC = "PAAD")))

intersect(unique(df_model_depmap$OncotreeCode), unique(df_model_cptac$Model.CancerType))
intersect(unique(tumor_types$oncotree_code), unique(df_model_cptac$Model.CancerType))

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

df_depmap <- model_buf_depmap %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2)))

df_procan <- model_buf_procan %>%
  inner_join(y = df_model_depmap %>% select(Model.ID, OncotreeCode), by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = df_model_procan %>% select(Model.ID, growth_properties), by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2))) %>%
  rename(GrowthPattern = growth_properties)

df_cptac <- model_buf_cptac %>%
  inner_join(y = df_model_cptac, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  rename(OncotreeCode = Model.CancerType) %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(Model.Ploidy = if_else(is.na(WGS_ploidy), WES_ploidy, WGS_ploidy),
         OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2)))

df_agg <- model_buf_agg %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  mutate(OncotreeCode = sapply(OncotreeCode, \(x) get_oncotree_parent(tumor_types, x, target_level = 2)))

# === Cell Lines Panel ===
waterfall_procan <- model_buf_procan %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name, centrality_measure = median,
                 color_low = color_palettes$BufferingClasses["Scaling"],
                 color_high = color_palettes$BufferingClasses["Buffered"])
waterfall_depmap <- model_buf_depmap %>%
  waterfall_plot(Model.Buffering.Ratio.ZScore, Rank, CellLine.Name, centrality_measure = median,
                 color_low = color_palettes$BufferingClasses["Scaling"],
                 color_high = color_palettes$BufferingClasses["Buffered"]) +
  labs(y = "Sample Buffering Ratio (z-score)")

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
                     default_fill_color = color_palettes$BufferingClasses["Buffered"], text_color = "black",
                     bar_label_shift = 0.1, break_steps = 0.25,
                     value_range = c(0,1), line_intercept = 0.5, value_lab = "Mean Normalized Rank")

plot_agg_bot <- model_buf_agg %>%
  slice_min(Model.Buffering.MeanNormRank, n = 10) %>%
  vertical_bar_chart(CellLine.Name, Model.Buffering.MeanNormRank,
                     default_fill_color = color_palettes$BufferingClasses["Scaling"], text_color = "black",
                     bar_label_shift = 0.1, break_steps = 0.25,
                     value_range = c(0,1), line_intercept = 0.5, value_lab = "Mean Normalized Rank")

plot_bracket <- plot_corr_bracket(cellline_pearson_gene)
plot_stack1 <- cowplot::plot_grid(waterfall_depmap + ylim(-2, 6),
                                  waterfall_procan + ylim(-2, 6) + ylab(NULL),
                                  nrow = 1, ncol = 2, align = "v", axis = "tb", labels = c("DepMap", "ProCan"),
                                  label_y = 0.98, label_x = 0.1, rel_widths = c(1, 1))
panel_celllines <- cowplot::plot_grid(plot_bracket, plot_stack1,
                                      nrow = 2, ncol = 1,
                                      rel_heights = c(0.1, 1), align = "h", axis = "lr")

panel_agg <- cowplot::plot_grid(plot_agg_top + ylab(NULL) + cowplot::theme_minimal_vgrid()
                                 + theme(axis.text.x = element_blank()),
                                plot_agg_bot + cowplot::theme_minimal_vgrid(),
                                nrow = 2, ncol = 1, align = "v", axis = "lr")

panel_celllines_agg <- cowplot::plot_grid(panel_celllines, panel_agg,
                                          nrow = 1, ncol = 2, labels = c("A", "B"),
                                          rel_widths = c(1.8, 1))

# === Cancer Types Panel ===
df_cancer_heatmap <- bind_rows(df_depmap, df_procan, df_cptac) %>%
  mutate(Suspension = GrowthPattern == "Suspension") %>%
  group_by(Dataset, OncotreeCode) %>%
  summarize(Mean.BR = mean(Model.Buffering.Ratio, na.rm = TRUE),
            Mean.AS = mean(CellLine.AneuploidyScore, na.rm = TRUE),
            Mean.Suspension = mean(Suspension, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(OncotreeCode = fct_reorder(OncotreeCode, Mean.BR)) %>%
  drop_na(OncotreeCode) %>%
  complete(Dataset, OncotreeCode) %>%
  mutate_all(~if_else(is.nan(.), NA, .)) %>%
  mutate(Dataset = factor(Dataset, levels = rev(dataset_order)))

cancer_heatmap_br <- df_cancer_heatmap %>%
  mutate(Suspension = case_when(Mean.Suspension >= 0.5 ~ "High",
                                Mean.Suspension >= 0.1 ~ "Low",
                                TRUE ~ "None")) %>%
  ggplot() +
  aes(x = OncotreeCode, y = Dataset, fill = Mean.BR) +
  geom_tile(aes(color = Mean.AS), alpha = 0) +
  geom_tile() +
  geom_point(aes(shape = Suspension), color = "white", size = 3) + # Shape 18
  scale_color_viridis_c(na.value = color_palettes$Missing, option = color_palettes$AneuploidyScore, end = 0.8) +
  scale_fill_viridis_c(na.value = color_palettes$Missing, option = color_palettes$BufferingRatio, end = 0.8, begin = 0.1) +
  scale_shape_manual(values = c(High = 18, Low = 5, None = NA), labels = c("\u226550%", "\u226510%", "")) +
  #scale_alpha_continuous(range = c(0, 1)) +
  scale_x_discrete(position = "top") +
  theme_void() +
  labs(x = NULL, fill = "Mean\nBuffering Ratio", color = "Mean\nAneuploidy", shape = "Suspension\nCulture") +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_colourbar(order = 2),
         shape = guide_legend(order = 3, override.aes = list(color = "black"))
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        legend.key.size = unit(16, "points"),
        legend.box = "horizontal",
        legend.position = "right",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0))

cancer_heatmap_as <- df_cancer_heatmap %>%
  drop_na() %>%
  summarize(Mean.AS = mean(Mean.AS, na.rm = TRUE),
            Dataset = "Mean Aneuploidy", .by = "OncotreeCode") %>%
  ggplot() +
  aes(x = OncotreeCode, y = Dataset, fill = Mean.AS) +
  geom_tile() +
  scale_fill_viridis_c(na.value = color_palettes$Missing, option = color_palettes$AneuploidyScore, end = 0.8) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

panel_types <- cowplot::plot_grid(cancer_heatmap_br, cancer_heatmap_as,
                                  nrow = 2, ncol = 1, align = "v", axis = "tb",
                                  rel_heights = c(1, 0.2))

# == Lymphoma Leukemia Panel ===
leukemia_codes <- c("MNM", "LNM")

df_leuk <- df_depmap %>% filter(OncotreeCode %in% leukemia_codes)

min_aneuploidy <- min(df_depmap$CellLine.AneuploidyScore)
max_aneuploidy <- max(df_depmap$CellLine.AneuploidyScore)
max_aneuploidy_leuk <- round(quantile(df_leuk$CellLine.AneuploidyScore, probs = 0.9))

leuk_plot <- df_depmap %>%
  mutate(`Myeloid / Lymphoid` = if_else(OncotreeCode %in% leukemia_codes, "Myeloid /\nLymphoid", "Other")) %>%
  signif_beeswarm_plot(`Myeloid / Lymphoid`, Model.Buffering.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test, count_y = 0)

leuk_plot_low_as <- df_depmap %>%
  filter(CellLine.AneuploidyScore <= max_aneuploidy_leuk) %>%
  mutate(`Myeloid / Lymphoid` = if_else(OncotreeCode %in% leukemia_codes, "Myeloid /\nLymphoid", "Other")) %>%
  signif_beeswarm_plot(`Myeloid / Lymphoid`, Model.Buffering.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test, count_y = 0)

leuk_plot <- leuk_plot +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.margin = margin(0,0,0,0, unit = 'cm')) +
  guides(color = guide_colourbar(ticks.linewidth = 1)) +
  labs(color = "\nAneuploidy Score", y = "Sample Buffering Ratio", x = NULL) +
  scale_colour_viridis_c(option = color_palettes$AneuploidyScore, direction = 1, end = 0.8,
                         limits = c(min_aneuploidy, max_aneuploidy),
                         breaks = c(min_aneuploidy, max_aneuploidy))

leuk_plot_low_as <- leuk_plot_low_as +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.margin = margin(0,0,0,0, unit = 'cm')) +
  guides(color = guide_colourbar(ticks.linewidth = 1)) +
  labs(color = "Aneuploidy Score\n(Limited Range)", y = NULL, x = NULL) +
  scale_colour_viridis_c(option = color_palettes$AneuploidyScore, direction = 1, end = 0.8,
                         limits = c(min_aneuploidy, max_aneuploidy),
                         breaks = c(min_aneuploidy, max_aneuploidy_leuk[["90%"]], max_aneuploidy),
                         labels = c(min_aneuploidy, max_aneuploidy_leuk[["90%"]], max_aneuploidy))

panel_leuk <- cowplot::plot_grid(leuk_plot,
                                 leuk_plot_low_as,
                                 nrow = 1, ncol = 2)

# === Growth Pattern Panel ===
df_low_as <- df_procan %>%
  filter(GrowthPattern %in% c("Adherent", "Suspension")) %>%
  filter(CellLine.AneuploidyScore <= min(
    max(df_procan[df_procan$GrowthPattern == "Adherent",]$CellLine.AneuploidyScore),
    max(df_procan[df_procan$GrowthPattern == "Suspension",]$CellLine.AneuploidyScore)
  )) %>%
  mutate(Condition = "Low Aneuploidy")
df_split_growth <- split(df_procan, df_procan$GrowthPattern)
df_equal_as <- df_split_growth$Suspension %>%
  equalize_distributions(df_split_growth$Adherent, CellLine.AneuploidyScore,
                         with_replacement = FALSE, num_buckets = 6) %>%
  mutate(Condition = "Equal Aneuploidy")

panel_growth <- bind_rows(df_procan %>% mutate(Condition = "Uncontrolled"), df_low_as, df_equal_as) %>%
  mutate(Condition = factor(Condition, levels = c("Uncontrolled", "Low Aneuploidy", "Equal Aneuploidy")),
         GrowthPattern = str_replace(GrowthPattern, "Semi-Adherent", "Mixed")) %>%
  filter(GrowthPattern != "Unknown") %>%
  ggplot() +
  aes(x = GrowthPattern, y = Model.Buffering.Ratio, color = GrowthPattern) +
  geom_boxplot(outliers = FALSE, size = 1, alpha = 0) +
  stat_summary(aes(y = 0.2), fun.data = show.n,
               geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("Adherent", "Suspension")),
              map_signif_level = print_signif, y_position = 1, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c(Adherent = discrete_color_pal2_bright[3],
                                Suspension = discrete_color_pal2_bright[2],
                                Mixed = default_color), guide = NULL) +
  labs(x = "Growth Pattern", y = "Sample Buffering Ratio") +
  facet_grid(~Condition, scales = "free_x", space = "free_x")

## Validate that aneuploidy score distributions are equal
wilcox.test(CellLine.AneuploidyScore ~ GrowthPattern, data = df_equal_as)

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
           label = paste0(print_corr(cor_as_test$estimate), ", ", print_signif(cor_as_test$p.value))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "Aneuploidy Score", y = "Sample Buffering Ratio") +
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
             label = print_signif(wgd_as_test$p.value),
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
             label = print_signif(wgd_br_test$p.value),
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

# === Combine Panels into Figure ===
figure2_sub1 <- cowplot::plot_grid(panel_wgd, panel_growth,
                                   nrow = 1, ncol = 2, rel_widths = c(0.7, 1),
                                   labels = c("D", "E"))

figure2 <- cowplot::plot_grid(panel_celllines_agg, panel_types, figure2_sub1,
                              labels = c("", "C", ""),
                              nrow = 3, ncol = 1, rel_heights = c(1, 0.4, 1))

cairo_pdf(here(plots_dir, "figure02.pdf"), width = 12, height = 12)
figure2
dev.off()

# === Tables ===
t3_field_descriptions <- c(
  "Dataset" = "Proteomics dataset used for analysis (e.g., DepMap, ProCan, CPTAC, etc.).",
  "Model.ID" = "Unique identifier for cell lines and tumor samples.",
  "CellLine.DepMapModelId" = "Unique identifier for cell lines; Provided by DepMap CCLE.",
  "CellLine.SangerModelId" = "Unique identifier for cell lines; Provided by Sanger Cell Model Passports.",
  "CellLine.Name" = "Name of the cell line.",
  "OncotreeCode" = "Cancer type classification using OncoTree codes (hierarchy level 2).",
  "Model.Buffering.Ratio" = "Sample buffering ratio (sample BR) calculated as the mean gene copy number-derived buffering ratio of (filtered) proteins in a cell line or tumor sample. Proteins missing in at least one of the datasets were removed prior to calculation.",
  "Model.Buffering.Ratio.ZScore" = "Standardized sample BR calculated as the z-score of Model.Buffering.Ratio.",
  "Observations" = "Number of valid gene copy number-derived buffering ratios in a cell line or tumor sample after filtering.",
  "SD" = "Standard deviation of gene copy number-derived buffering ratios of (filtered) proteins in a cell line or tumor sample.",
  "Rank" = "Rank of the sample BR (ascending).",
  "Model.Buffering.MeanNormRank" = "Mean normalized rank of sample BRs of a cell line across datasets (DepMap, ProCan).",
  "Model.Buffering.StandardizedMean" = "Mean of the standardized sample BRs of a cell line across datasets (DepMap, ProCan).",
  "Model.AneuploidyScore.Estimate" = "Estimated aneuploidy score quantifying the degree of aneuploidy in tumor samples. See manuscript for details."
)

fields <- names(t3_field_descriptions)
df_t3_fields <- data.frame(Column = names(t3_field_descriptions), Description = unname(t3_field_descriptions))

wb <- createWorkbook()
sheet_readme <- addWorksheet(wb, "README")
sheet_depmap <- addWorksheet(wb, "DepMap")
sheet_procan <- addWorksheet(wb, "ProCan")
sheet_cptac <- addWorksheet(wb, "CPTAC")
sheet_agg  <- addWorksheet(wb, "Cell Lines (Aggregated)")
writeDataTable(wb = wb, sheet = sheet_readme, x = df_t3_fields)
writeDataTable(wb = wb, sheet = sheet_depmap, x = df_depmap %>% select(any_of(fields)))
writeDataTable(wb = wb, sheet = sheet_procan, x = df_procan %>% select(any_of(fields)))
writeDataTable(wb = wb, sheet = sheet_cptac, x = df_cptac %>% select(any_of(fields)))
writeDataTable(wb = wb, sheet = sheet_agg, x = df_agg %>% select(any_of(fields)))
saveWorkbook(wb, here(tables_dir, "supplementary_table3.xlsx"), overwrite = TRUE)

# === Supplemental Figures ===
## Kaplan-Meyer Plots
df_cptac_split <- df_cptac %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Buffering")

surv_os <- survfit(Surv(OS_days, OS_event) ~ Buffering, data = df_cptac_split) %>%
  ggsurvplot(data = df_cptac_split, pval = TRUE)
surv_pfs <- survfit(Surv(PFS_days, PFS_event) ~ Buffering, data = df_cptac_split) %>%
  ggsurvplot(data = df_cptac_split, pval = TRUE)

panel_surv <- cowplot::plot_grid(surv_os$plot + labs(x = NULL, y = "Overall\nsurvival probability") + theme(legend.position = "none"),
                                 surv_pfs$plot + labs(y = "Progression-free\nsurvival probability") + theme(legend.margin = margin(0, 0, 0, 0, unit = 'cm')),
                                 nrow = 2, rel_heights = c(0.75, 1))

## Age
cor_age <- cor.test(df_agg$Model.Buffering.MeanNormRank,
                       df_agg$Age, method = "spearman")

panel_age_cor <- df_agg %>%
  ggplot() +
  aes(y = Model.Buffering.MeanNormRank, x = Age) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = lm, color = highlight_color) +
  annotate("text", x = 0, y = 1.1, hjust = 0, size = 5,
           label = paste0(print_corr(cor_age$estimate), ", ", print_signif(cor_age$p.value))) +
  labs(y = "Sample BR MNR", x = "Age")

panel_age_cat <- df_agg %>%
  filter(AgeCategory %in% c("Pediatric", "Adult")) %>%
  signif_beeswarm_plot(AgeCategory, Model.Buffering.MeanNormRank,
                       color_col = CellLine.AneuploidyScore, viridis_color_pal = color_palettes$AneuploidyScore,
                       color_lims = c(0, 0.8), cex = 1, test = wilcox.test) +
  labs(x = "Age Group", y = "Sample BR MNR", color = "Aneuploidy Score") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.margin = margin(0,0,0,0, unit = 'cm'))

## TP53
df_p53 <- bind_rows(df_depmap, df_procan) %>%
  left_join(y = tp53_mut, by = "Model.ID") %>%
  bind_rows(df_cptac %>% mutate(TP53.Mutated = as.logical(TP53_mutation))) %>%
  drop_na(TP53.Mutated) %>%
  mutate(Aneuploidy.Estimate = if_else(!is.na(CellLine.AneuploidyScore),
                                       CellLine.AneuploidyScore, Model.AneuploidyScore.Estimate),
         Dataset = factor(Dataset, levels = dataset_order))

panel_p53 <- df_p53 %>%
  signif_beeswarm_plot(TP53.Mutated, Model.Buffering.Ratio,
                       facet_col = Dataset, color_col = Aneuploidy.Estimate,
                       viridis_color_pal = color_palettes$AneuploidyScore, cex = 1, color_lims = c(0, 0.8)) +
  labs(x = "Damaging TP53 Mutation", y = "Sample Buffering Ratio", color = "(Estimated) Aneuploidy Score") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.margin = margin(0,0,0,0, unit = 'cm'))

## Ploidy
cor_ploidy <- cor.test(df_procan$Model.Buffering.Ratio,
                       df_procan$CellLine.Ploidy, method = "spearman")

panel_ploidy <- df_procan %>%
  ggplot() +
  aes(y = Model.Buffering.Ratio, x = CellLine.Ploidy) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = loess, color = highlight_color) +
  annotate("text", x = 0, y = 2, hjust = 0, size = 5,
           label = paste0(print_corr(cor_ploidy$estimate), ", ", print_signif(cor_ploidy$p.value))) +
  labs(y = "Sample Buffering Ratio", x = "Cell Line Ploidy")

## Aneuploidy Score (all)
panel_as_all <- bind_rows(df_depmap, df_procan, df_cptac) %>%
  mutate(Aneuploidy.Estimate = if_else(!is.na(CellLine.AneuploidyScore),
                                       CellLine.AneuploidyScore, Model.AneuploidyScore.Estimate),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  distinct(Dataset, Model.ID, Model.Buffering.Ratio, CellLine.WGD, Aneuploidy.Estimate) %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD")) %>%
  ggplot() +
  aes(x = Aneuploidy.Estimate, y = Model.Buffering.Ratio, color = WGD) +
  geom_point(size = 2) +
  stat_smooth(method = lm, color = "dimgrey") +
  #annotate("text", x = 0, y = 1.8, hjust = 0, size = 5, label = paste0(print_corr(cor_as_test$estimate), ", ", print_signif(cor_as_test$p.value))) +
  stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE,
           p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "(Estimated) Aneuploidy Score", y = "Sample Buffering Ratio") +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.85),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size)) +
  facet_grid(~Dataset, scales = "free_x")

## Growth Pattern (DepMap)
df_low_as_depmap <- df_depmap %>%
  filter(GrowthPattern %in% c("Adherent", "Suspension")) %>%
  filter(CellLine.AneuploidyScore <= min(
    max(df_depmap[df_depmap$GrowthPattern == "Adherent",]$CellLine.AneuploidyScore),
    max(df_depmap[df_depmap$GrowthPattern == "Suspension",]$CellLine.AneuploidyScore)
  )) %>%
  mutate(Condition = "Low Aneuploidy")
df_split_growth_depmap <- split(df_depmap, df_depmap$GrowthPattern)
df_equal_as_depmap <- df_split_growth_depmap$Suspension %>%
  equalize_distributions(df_split_growth_depmap$Adherent, CellLine.AneuploidyScore,
                         with_replacement = FALSE, num_buckets = 6) %>%
  mutate(Condition = "Equal Aneuploidy")

panel_growth_depmap <- bind_rows(df_depmap %>% mutate(Condition = "Uncontrolled"), df_low_as_depmap, df_equal_as_depmap) %>%
  mutate(Condition = factor(Condition, levels = c("Uncontrolled", "Low Aneuploidy", "Equal Aneuploidy"))) %>%
  filter(GrowthPattern %in% c("Adherent", "Suspension")) %>%
  ggplot() +
  aes(x = GrowthPattern, y = Model.Buffering.Ratio, color = GrowthPattern) +
  geom_boxplot(outliers = FALSE, size = 1, alpha = 0) +
  stat_summary(aes(y = 0.2), fun.data = show.n,
               geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("Adherent", "Suspension")),
              map_signif_level = print_signif, y_position = 0.9, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c(Adherent = discrete_color_pal2_bright[3],
                                Suspension = discrete_color_pal2_bright[2]),
                                guide = NULL) +
  labs(x = "Growth Pattern", y = "Sample Buffering Ratio") +
  facet_grid(~Condition, scales = "free_x", space = "free_x")

## Aneuploidy Score per growth pattern
panel_growth_as <- bind_rows(df_depmap %>% mutate(Condition = "Uncontrolled"), df_low_as_depmap, df_equal_as_depmap) %>%
  mutate(Condition = factor(Condition, levels = c("Uncontrolled", "Low Aneuploidy", "Equal Aneuploidy"))) %>%
  filter(GrowthPattern %in% c("Adherent", "Suspension")) %>%
  ggplot() +
  aes(x = GrowthPattern, y = CellLine.AneuploidyScore, color = GrowthPattern) +
  geom_boxplot(outliers = FALSE, size = 1, alpha = 0) +
  stat_summary(aes(y = 0), fun.data = show.n,
               geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("Adherent", "Suspension")),
              map_signif_level = print_signif, y_position = 30, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c(Adherent = discrete_color_pal2_bright[3],
                                Suspension = discrete_color_pal2_bright[2]),
                                guide = NULL) +
  labs(x = "Growth Pattern", y = "Aneuploidy Score") +
  facet_grid(~Condition)

## Combine figures
figure_s2_sub1 <- cowplot::plot_grid(panel_as_all, panel_ploidy,
                                     ncol = 2, rel_widths = c(1, 0.4), labels = c("A", "B"))
figure_s2_sub2 <- cowplot::plot_grid(panel_surv, panel_age_cat, panel_p53,
                                     ncol = 3, rel_widths = c(0.8, 0.5, 1), labels = c("C", "D", "E"))
figure_s2_sub3 <- cowplot::plot_grid(panel_growth_as, panel_growth_depmap,
                                     ncol = 2, rel_widths = c(1, 1), labels = c("F", "G"))

figure_s2 <- cowplot::plot_grid(figure_s2_sub1, figure_s2_sub2, figure_s2_sub3, nrow = 3)

cairo_pdf(here(plots_dir, "figure_s2.pdf"), width = 12, height = 12)
figure_s2
dev.off()
