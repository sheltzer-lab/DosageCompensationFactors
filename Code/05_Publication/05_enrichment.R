library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(STRINGdb)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
temp_dir <- temp_base_dir
output_data_dir <- output_data_base_dir
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
tables_dir <- here(tables_base_dir, "Publication")

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))
expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))
metadata_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

diff_exp_cptac <- read_parquet(here(output_data_dir, "model_buf_diff-exp_cptac.parquet"))
diff_exp_procan <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))
diff_exp_depmap <- read_parquet(here(output_data_dir, "model_buf_diff-exp_depmap.parquet"))
diff_exp_control <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan_adherent.parquet"))

ssgsea_cptac <- read_parquet(here(output_data_dir, "ssgsea_unfolded_cptac.parquet"))
gsea_all <- read_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer.parquet"))

df_growth <- read_parquet(here(output_data_dir, "cellline_growth.parquet")) %>%
  select(Model.ID, CellLine.GrowthRatio)

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.DepMapModelId = "ModelID") %>%
  inner_join(y = df_celllines, by = "CellLine.DepMapModelId") %>%
  select(-CellLine.Name)

df_agg <- model_buf_agg %>%
  inner_join(y = df_model_depmap, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  left_join(y = copy_number, by = "Model.ID", relationship = "many-to-one", na_matches = "never")

# === Volcano Plot Panel ===
color_mapping <- scale_color_manual(values = color_palettes$DiffExp,
                                    na.value = color_palettes$Missing)

top_diff <- diff_exp_procan %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  filter(Test.p.adj < p_threshold) %>%
  slice_max(abs(Log2FC), n = 10)

diff_exp_procan <- diff_exp_procan  %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & (!is.na(CancerDriverMode) | Gene.Symbol %in% top_diff$Gene.Symbol),
                         Gene.Symbol, NA))

panel_volcano <- diff_exp_procan %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping) +
  theme(legend.position = "none") +
  labs(y = "-log10(p.adj)", x = "Protein Log2FC (High - Low Buffering)")

# === ORA Panel ===
genes_up <- diff_exp_procan %>%
  filter(Significant == "Up") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

genes_down <- diff_exp_procan %>%
  filter(Significant == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_up <- genes_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("REAC", "GO:MF", "CORUM"),
                     custom_color = color_palettes$DiffExp["Up"])

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("KEGG", "GO:MF", "CORUM"),
                     custom_color = color_palettes$DiffExp["Down"])

# === STRING Panel ===
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                          network_type = "full", input_directory = temp_dir)

string_up <- genes_up %>%
  create_string_network(Gene.Symbol, Log2FC, string_db, min_score = 900)

string_down <- genes_down %>%
  create_string_network(Gene.Symbol, Log2FC, string_db, min_score = 900)

# Use link in string_up & string_down to create plot manually

# TODO: Mention adjusted BR
# TODO: Why are more samples in df_agg than in df_procan and df_depmap?

# === 2D Enrichment Panel ===
max_abs_nes_procan <- gsea_all %>% filter(Dataset == "ProCan") %>% pull(NES) %>% abs() %>% max() %>% round()
max_abs_nes_cptac <- gsea_all %>% filter(Dataset == "CPTAC") %>% pull(NES) %>% abs() %>% max() %>% round()

selected_pathways <- c("DNA_REPAIR", "APOPTOSIS", "G2M_CHECKPOINT", "MYC_TARGETS_V2", "IL2_STAT5_SIGNALING",
                       "UNFOLDED_PROTEIN_RESPONSE", "EPITHELIAL_MESENCHYMAL_TRANSITION", "APICAL_JUNCTION")

panel_2d_enrichment <- gsea_all %>%
  filter(Dataset %in% c("CPTAC", "ProCan")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_ProCan < p_threshold & padj_CPTAC < p_threshold ~ "Both",
    padj_ProCan >= p_threshold & padj_CPTAC < p_threshold ~ "CPTAC",
    padj_ProCan < p_threshold & padj_CPTAC >= p_threshold ~ "ProCan",
    TRUE ~ "None"),
         Significant = factor(Significant, levels = c("ProCan", "CPTAC", "Both", "None")),
         Label = str_replace(pathway, "HALLMARK_", ""),
         Label = if_else(Significant == "Both" | Label %in% selected_pathways,
                         Label, NA)) %>%
  ggplot() +
  aes(x = NES_ProCan, y = NES_CPTAC, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(xintercept = 0, yintercept = 0, slope = 1, color = default_color) +
  geom_point(size = 2) +
  geom_label_repel(size = ceiling(base_size / 4), force = 30, min.segment.length = 0.01) +
  scale_color_manual(values = c(Both = highlight_colors[1], color_palettes$Datasets,
                                None = default_color)) +
  lims(x = c(-max_abs_nes_procan, max_abs_nes_procan), y = c(-max_abs_nes_cptac, max_abs_nes_cptac)) +
  labs(x = "Normalized Enrichment Score (ProCan)",
       y = "Normalized Enrichment Score (CPTAC)") +
  theme(legend.position = "top", legend.direction = "horizontal")

# === Proteotoxic Stress Panel ===
quantile(model_buf_cptac$Observations, probs = seq(0,1,0.1))

ssgsea_cptac_filtered <- ssgsea_cptac %>%
  inner_join(y = model_buf_cptac, by = "Model.ID") %>%
  filter(Observations > 500) %>%
  left_join(y = metadata_cptac %>% select(-HALLMARK_UNFOLDED_PROTEIN_RESPONSE), by = "Model.ID")

cor_proteotox_test <- cor.test(ssgsea_cptac_filtered$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                               ssgsea_cptac_filtered$Model.Buffering.Ratio, method = "spearman")

panel_proteotox <- ssgsea_cptac_filtered %>%
  mutate(Ploidy = if_else(Model.Ploidy >= 3, "\u2265 3", "< 3")) %>%
  ggplot() +
  aes(x = Model.Buffering.Ratio, y = HALLMARK_UNFOLDED_PROTEIN_RESPONSE, color = Model.AneuploidyScore.Estimate) +
  geom_point(aes(shape = Ploidy), size = 3) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0.35, y = 0.275, hjust = 0, size = 5,
           label = paste0(print_corr(cor_proteotox_test$estimate), ", ", print_signif(cor_proteotox_test$p.value))) +
  labs(x = "Sample Buffering Ratio", y = "Unfolded Protein Response", color = "eAS") +
  ylim(c(0.085, 0.275)) +
  scale_color_viridis_c(option = color_palettes$AneuploidyScore, end = 0.8) +
  theme(legend.position = c("left", "bottom"),
        legend.position.inside = c(0, 0),
        legend.justification = c("left", "bottom"),
        legend.box = "horizontal",
        legend.text = element_text(size = base_size),
        legend.background = element_rect(fill = alpha('white', 2/3))) +
  guides(colour = guide_colourbar(order = 1),
         shape = guide_legend(order = 2))

## Check if tumor purity influences UPR
ssgsea_cptac_filtered %>%
  scatter_plot_reg_corr(HALLMARK_UNFOLDED_PROTEIN_RESPONSE, Model.TumorPurity)

ssgsea_cptac_filtered %>%
  count(HALLMARK_UNFOLDED_PROTEIN_RESPONSE > 0.2, Model.TumorPurity > 0.7) %>%
  pivot_wider(names_from = "Model.TumorPurity > 0.7", values_from = "n") %>%
  tibble::column_to_rownames("HALLMARK_UNFOLDED_PROTEIN_RESPONSE > 0.2") %>%
  as.matrix() %>%
  fisher.test()

# === Proliferation ===
df_prolif <- df_agg %>%
  select(Model.ID, Model.Buffering.MeanNormRank, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy) %>%
  inner_join(y = df_growth, by = "Model.ID") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"),
         Model.Buffering.Adjusted = Model.Buffering.MeanNormRank / (1+CellLine.AneuploidyScore/max(CellLine.AneuploidyScore)),
         Aneuploidy = factor(if_else(CellLine.AneuploidyScore > median(CellLine.AneuploidyScore),
                                     "High Aneuploidy", "Low Aneuploidy"), levels = c("Low Aneuploidy", "High Aneuploidy")))

cor_prolif <- cor.test(df_prolif$Model.Buffering.MeanNormRank,
                       df_prolif$CellLine.GrowthRatio, method = "spearman")

cor_prolif_adj <- cor.test(df_prolif$Model.Buffering.Adjusted,
                           df_prolif$CellLine.GrowthRatio, method = "spearman")

cor_prolif_as <- cor.test(df_prolif$CellLine.AneuploidyScore,
                          df_prolif$CellLine.GrowthRatio, method = "spearman")

df_prolif %>%
  group_by(WGD) %>%
  rstatix::cor_test(Model.Buffering.MeanNormRank, CellLine.GrowthRatio, method = "spearman")

df_prolif %>%
  group_by(WGD) %>%
  rstatix::cor_test(CellLine.AneuploidyScore, CellLine.GrowthRatio, method = "spearman")

panel_prolif_buf <- df_prolif %>%
  ggplot() +
  aes(x = Model.Buffering.MeanNormRank, y = CellLine.GrowthRatio, color = WGD) +
  geom_point(aes(shape = Aneuploidy), size = 3) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0, y = 14, hjust = 0, size = 5,
           label = paste0(print_corr(cor_prolif$estimate), ", ", print_signif(cor_prolif$p.value))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "Mean Normalized Sample Buffering Ranks", y = "Growth Ratio (Day4/Day1)") +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.90),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        legend.background = element_rect(fill = alpha('white', 2/3)),
        legend.direction = "horizontal",
        legend.spacing.y = unit(-10, "pt"))

panel_prolif_as <- df_prolif %>%
  ggplot() +
  aes(x = CellLine.AneuploidyScore, y = CellLine.GrowthRatio, color = WGD) +
  geom_point(size = 3) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0, y = 14, hjust = 0, size = 5,
           label = paste0(print_corr(cor_prolif_as$estimate), ", ", print_signif(cor_prolif_as$p.value))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "Aneuploidy Score", y = "Growth Ratio (Day4/Day1)") +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.90),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        legend.background = element_rect(fill = alpha('white', 2/3)))

## Multiple Linear Regression with variable standardization
reg_as_br <- lm(CellLine.GrowthRatio ~ scale(Model.Buffering.MeanNormRank) +
  scale(CellLine.AneuploidyScore) + scale(CellLine.WGD),
                data = df_prolif)
summary(reg_as_br)

# === Combine Panels into Figure ===
figure5_01 <- cowplot::plot_grid(panel_volcano, ora_down, ora_up,
                              labels = c("A", "B", "C"),
                              rel_widths = c(1, 1, 1),
                              nrow = 1, ncol = 3)

cairo_pdf(here(plots_dir, "figure05_01.pdf"), width = 13, height = 5)
figure5_01
dev.off()

# D & E reserved for STRING clusters

fig5_02_sub1 <- cowplot::plot_grid(panel_proteotox, panel_prolif_buf,
                                   nrow = 2, ncol = 1, labels = c("G", "H"))

figure5_02 <- cowplot::plot_grid(panel_2d_enrichment, fig5_02_sub1,
                              labels = c("F", ""),
                              rel_widths = c(1, 0.7),
                              nrow = 1, ncol = 2)

cairo_pdf(here(plots_dir, "figure05_02.pdf"), width = 12, height = 8)
figure5_02
dev.off()

# === Tables ===
diff_exp_all <- bind_rows(diff_exp_depmap %>% mutate(Dataset = "DepMap"),
                          diff_exp_procan %>% mutate(Dataset = "ProCan"),
                          diff_exp_control %>% mutate(Dataset = "ProCan (adherent control)"),
                          diff_exp_cptac %>% mutate(Dataset = "CPTAC")) %>%
  select(-CancerDriverMode, -Occurrences, -Label) %>%
  mutate(GroupA = "Low Buffering", GroupB = "High Buffering")

diff_exp_common <- diff_exp_all %>%
  drop_na(Significant) %>%
  filter(Dataset != "CPTAC") %>%
  add_count(Gene.Symbol, Significant) %>%
  filter(n == 3) %>%
  arrange(Gene.Symbol) %>%
  select(-n)

ssgsea_cptac_export <- ssgsea_cptac %>%
  left_join(y = metadata_cptac %>% select(starts_with("Model.")), by = "Model.ID") %>%
  mutate(Dataset = "CPTAC")

wb <- createWorkbook()
sheet_diffexp <- addWorksheet(wb, "Differential Expression")
sheet_diffexp_common <- addWorksheet(wb, "DiffExp (common genes)")
sheet_gsea <- addWorksheet(wb, "GSEA")
sheet_ssgsea <- addWorksheet(wb, "ssGSEA")
writeDataTable(wb = wb, sheet = sheet_diffexp, x = diff_exp_all)
writeDataTable(wb = wb, sheet = sheet_diffexp_common, x = diff_exp_common)
writeDataTable(wb = wb, sheet = sheet_gsea, x = gsea_all)
writeDataTable(wb = wb, sheet = sheet_ssgsea, x = ssgsea_cptac_export)
saveWorkbook(wb, here(tables_dir, "supplementary_table8.xlsx"), overwrite = TRUE)

# === Supplemental Figures ===
string_common <- diff_exp_common %>%
  summarize(Log2FC.Mean = mean(Log2FC), .by = c("Gene.Symbol", "Significant")) %>%
  #filter(Significant == "Up") %>%
  create_string_network(Gene.Symbol, Log2FC.Mean, string_db, min_score = 900)

