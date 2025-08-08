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
source(here("Code", "publication.R"))

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

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))

ssgsea_cptac <- read_parquet(here(output_data_dir, "ssgsea_unfolded_cptac.parquet"))
gsea_all <- read_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer.parquet"))
gsea_all_as <- read_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer_aneuploidy.parquet"))

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
  slice_max(abs(Log2FC), n = 10) %>%
  pull(Gene.Symbol)

common_genes <- bind_rows(diff_exp_depmap %>% mutate(Dataset = "DepMap"),
                             diff_exp_procan %>% mutate(Dataset = "ProCan"),
                             diff_exp_control %>% mutate(Dataset = "ProCan (adherent control)"),
                             diff_exp_cptac %>% mutate(Dataset = "CPTAC")) %>%
  drop_na(Significant) %>%
  filter(Dataset != "CPTAC") %>%
  add_count(Gene.Symbol, Significant) %>%
  filter(n == 3) %>%
  distinct(Gene.Symbol) %>%
  pull(Gene.Symbol)

selected_common_genes <- c("CTSA", "UBE2N", "ITGA3", "ITGAV")

panel_volcano <- diff_exp_procan %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(IsCommon = Gene.Symbol %in% common_genes,
         Font = if_else(IsCommon, "bold", "plain"),
         Label = if_else(!is.na(Significant) & (CancerDriverMode %in% c("OG", "TSG", "OG/TSG") | Gene.Symbol %in% top_diff | Gene.Symbol %in% selected_common_genes),
                         Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, font_col = Font, color_mapping = color_mapping) +
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
         Label = str_replace_all(str_replace(pathway, "HALLMARK_", ""), "_", " "),
         Label = if_else(Significant != "None" | Label %in% selected_pathways,
                         Label, NA)) %>%
  ggplot() +
  aes(x = NES_ProCan, y = NES_CPTAC, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(xintercept = 0, yintercept = 0, slope = 1, color = default_color) +
  geom_point(size = 3) +
  geom_label_repel(size = ceiling(base_size / 4), force = 40, min.segment.length = 0.01) +
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
  labs(x = "Sample Buffering Ratio", y = "Unfolded Protein Response", color = "AS") +
  ylim(c(0.085, 0.275)) +
  scale_color_viridis_c(option = color_palettes$AneuploidyScore, end = 0.8) +
  theme(legend.position = c("left", "bottom"),
        legend.position.inside = c(0, 0),
        legend.justification = c("left", "bottom"),
        legend.box = "horizontal",
        legend.box.just = "bottom",
        legend.spacing.x = unit(0, "pt"),
        legend.text = element_text(size = base_size),
        legend.background = element_rect(fill = alpha('white', 2/3))) +
  guides(colour = guide_colourbar(order = 1),
         shape = guide_legend(order = 2))

## Check if tumor purity influences UPR
cor_upr_purity <- cor.test(ssgsea_cptac_filtered$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                           ssgsea_cptac_filtered$Model.TumorPurity)

panel_upr_purity <- ssgsea_cptac_filtered %>%
  ggplot() +
  aes(y = Model.TumorPurity, x = HALLMARK_UNFOLDED_PROTEIN_RESPONSE) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = lm, color = highlight_color) +
  annotate("text", x = 0.1, y = 1, hjust = 0, size = 5,
           label = paste0(print_corr(cor_upr_purity$estimate), ", ", print_signif(cor_upr_purity$p.value))) +
  labs(y = "Tumor Purity", x = "Unfolded Protein Response")

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
         Model.Buffering.Adjusted = Model.Buffering.MeanNormRank / (1+CellLine.AneuploidyScore/max(CellLine.AneuploidyScore)))

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
  aes(x = Model.Buffering.MeanNormRank, y = CellLine.GrowthRatio, color = CellLine.AneuploidyScore) +
  geom_point(aes(shape = WGD), size = 3) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0, y = 14, hjust = 0, size = 5,
           label = paste0(print_corr(cor_prolif$estimate), ", ", print_signif(cor_prolif$p.value))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_viridis_c(option = color_palettes$AneuploidyScore, end = 0.8) +
  labs(x = "Mean Normalized Sample Buffering Ranks", y = "Growth Ratio (Day4/Day1)", color = "AS") +
  guides(shape = guide_legend(title = NULL, order = 1),
         color = guide_colorbar(title.vjust = 0.9, order = 2)) +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.92),
        legend.justification = c("left", "top"),
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

gsea_export <- bind_rows(
  gsea_all %>% mutate(Comparison = "Buffering (High - Low)"),
  gsea_all_as %>% mutate(Comparison = "Aneuploidy (High - Low)")
)

t8_field_descriptions <- c(
  "=== TABLES ===" = "",
  "Differential Expression" = "Contains the differential expression profiles of genes in high buffering samples (cell lines or tumor samples; >80% sample buffering ratio) against low buffering samples (<20% sample BR) across different datasets (DepMap, ProCan, CPTAC, etc.).",
  "DiffExp (common genes)" = "Contains the differential expression profiles of genes that are commonly up- or downregulated in all cell line datasets.",
  "GSEA" = "Gene Set Enrichment Analysis results for the HALLMARK gene sets from MSigDB; Generated using the ranks of the differential expression profiles. This table was generated using the fGSEA R-package. See documentation of fgsea::fgseaMultilevel() for further details.",
  "ssGSEA" = "Single-sample Gene Set Enrichment Analysis results for the HALLMARK gene sets from MSigDB; ssGSEA scores are generated for each HALLMARK gene set and each tumor sample in CPTAC using the normalized protein abundance (supplementary table 1). This table was generated using the GSVA R-package. See documentation of GSVA::gsva() for further details.",
  "=== COLUMNS ===" = "",
  "Gene.Symbol" = "HGNC gene symbol; updated using HGNChelper.",
  "GroupA" = "First group in the statistical analysis. Here: Samples with low average buffering (<20% of sample buffering ratios of a dataset).",
  "GroupB" = "Second group in the statistical analysis. Here: Samples with high average buffering (>80% of sample buffering ratios of a dataset).",
  "Mean_GroupA" = "Mean normalized log2 protein abundance of GroupA.",
  "Mean_GroupB" = "Mean normalized log2 protein abundance of GroupB.",
  "Count_GroupA" = "Number of valid protein abundance values in GroupA.",
  "Count_GroupB" = "Number of valid protein abundance values in GroupB.",
  "Log2FC" = "Protein abundance log2 fold-change between high and low buffering samples for each gene. Calculated as Mean_GroupB - Mean_GroupA.",
  "Test.p" = "P-value of the two-sided Welch's unequal variances t-test used to determine whether the difference between GroupA and GroupB is significant.",
  "Test.p.adj" = "Benjamini-Hochberg-adjusted p-value.",
  "Significant" = "Indicates whether a protein was significantly up- or downregulated in high buffering samples (Test.p.adj < 0.05, |Log2FC| > 0.5). Either Up, Down, or empty.",
  "Dataset" = "Proteomics dataset used for analysis (e.g., DepMap, ProCan, CPTAC, etc.).",
  "pathway" = "Gene set used for the enrichment analysis.",
  "pval" = "P-value of the Gene Set Enrichment Analysis (GSEA).",
  "padj" = "Benjamini-Hochberg-adjusted p-value of the Gene Set Enrichment Analysis (GSEA).",
  "log2err" = "Expected error for the standard deviation of the GSEA p-value.",
  "ES" = "GSEA enrichment score, same as the Broad GSEA implementation. Positive values indicate an upregulation of the gene set.",
  "NES" = "Enrichment score normalized to the mean enrichment of random samples of the same size.",
  "size" = "Size of the gene set after removing genes not present in the differntial expression results.",
  "leadingEdge" = "Vector with leading edge genes that drive the enrichment, see http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading",
  "Comparison" = "Indicates which groups were compared for the enrichment analysis (high vs. low buffering, high vs. low aneuploidy).",
  "Model.ID" = "Unique identifier of the tumor sample.",
  "HALLMARK_*" = "Single-sample GSEA enrichment score for each gene set (columns) and sample (rows).",
  "Model.CancerType" = "CPTAC cancer type associated with te tumor sample.",
  "Model.TumorPurity" = "Purity of the tumor sample reported by CPTAC.",
  "Model.Ploidy" = "Ploidy of the tumor sample reported by CPTAC.",
  "Model.WGD.Estimate" = "Estimated WGD status of the tumor sample; Defined as Model.Ploidy >= 3.",
  "Model.AneuploidyScore.Estimate" = "Estimated aneuploidy score quantifying the degree of aneuploidy in tumor samples. See manuscript for details."
)

df_t8_fields <- data.frame(Column = names(t8_field_descriptions), Description = unname(t8_field_descriptions))

wb <- createWorkbook()
sheet_readme <- addWorksheet(wb, "README")
sheet_diffexp <- addWorksheet(wb, "Differential Expression")
sheet_diffexp_common <- addWorksheet(wb, "DiffExp (common genes)")
sheet_gsea <- addWorksheet(wb, "GSEA")
sheet_ssgsea <- addWorksheet(wb, "ssGSEA")
writeDataTable(wb = wb, sheet = sheet_readme, x = df_t8_fields)
writeDataTable(wb = wb, sheet = sheet_diffexp, x = diff_exp_all)
writeDataTable(wb = wb, sheet = sheet_diffexp_common, x = diff_exp_common)
writeDataTable(wb = wb, sheet = sheet_gsea, x = gsea_export)
writeDataTable(wb = wb, sheet = sheet_ssgsea, x = ssgsea_cptac_export)
saveWorkbook(wb, here(tables_dir, "supplementary_table8.xlsx"), overwrite = TRUE)

# === Supplemental Figures ===
## DiffExp Network for commonly deregulated genes
string_common <- diff_exp_common %>%
  summarize(Log2FC.Mean = mean(Log2FC), .by = c("Gene.Symbol", "Significant")) %>%
  #filter(Significant == "Up") %>%
  create_string_network(Gene.Symbol, Log2FC.Mean, string_db, min_score = 900)

## Volcano Plots (DepMap, CPTAC, Adherent Control)
top_diff_all <- diff_exp_all %>%
  filter(Test.p.adj < p_threshold & !is.na(Significant)) %>%
  filter(Dataset != "ProCan") %>%
  group_by(Dataset, Significant) %>%
  slice_max(abs(Log2FC), n = 5) %>%
  ungroup() %>%
  distinct(Gene.Symbol, .keep_all = TRUE)

panel_volcano_all <- diff_exp_all %>%
  filter(Dataset != "ProCan") %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) &
                           (!is.na(CancerDriverMode) | Gene.Symbol %in% top_diff_all$Gene.Symbol),
                         Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping = color_mapping) +
  theme(legend.position = "none") +
  labs(y = "-log10(p.adj)", x = "Protein Log2FC (High - Low Buffering)") +
  facet_grid(~Dataset)

## ORA (adherent control)
genes_up_control <- diff_exp_control %>%
  filter(Significant == "Up") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

genes_down_control <- diff_exp_control %>%
  filter(Significant == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_up_control <- genes_up_control %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "GO:MF", "REAC"),
                     custom_color = color_palettes$DiffExp["Up"], string_trunc = 45)

ora_down_control <- genes_down_control %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("GO:BP", "KEGG", "REAC"),
                     custom_color = color_palettes$DiffExp["Down"], string_trunc = 45)

## Cancer Driver Buffering
panel_og_tsg <- expr_buf_procan %>%
  inner_join(y = diff_exp_procan, by = "Gene.Symbol") %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(GeneGroup = case_when(CancerDriverMode == "TSG" & Significant == "Down" ~ "TSG (down)",
                               CancerDriverMode == "OG" & Significant == "Up" ~ "OG (up)",
                               TRUE ~ NA),
         CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.GeneLevel.Ratio, GeneGroup, CNV) %>%
  ggplot() +
  aes(x = CNV, y = Buffering.GeneLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  stat_summary(aes(y = -2), fun.data = \(x) show.n(x, prefix = "n="), geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = 2.5,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  facet_grid(~GeneGroup) +
  labs(x = "Gene Copy Number", y = "Buffering Ratio")

## Oncogene Buffering (chr arm)
panel_oncogenes <- expr_buf_procan %>%
  filter(Gene.Symbol %in% c("EGFR", "RRAS", "CTNNB1", "LMNA")) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CNA = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.ChrArmLevel.Ratio) %>%
  signif_beeswarm_plot(CNA, Buffering.ChrArmLevel.Ratio,
                       color_col = CellLine.AneuploidyScore, facet_col = Gene.Symbol,
                       viridis_color_pal = color_palettes$AneuploidyScore, color_lims = c(0, 0.8), cex = 1,
                       n.prefix = "n=", test = wilcox.test) +
  labs(x = "Chromosome Arm", y = "Buffering Ratio", color = "Aneuploidy Score") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(hjust = 1, vjust = 0.9),
        legend.margin = margin(0,0,0,0, unit = 'cm'))

## EGFR Buffering Classes (chr arm)
egfr_classes_chr <- expr_buf_procan %>%
  filter(Gene.Symbol == "EGFR") %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  drop_na(Buffering.ChrArmLevel.Ratio) %>%
  mutate(CNA = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Gain", "Loss"))

panel_egfr <- egfr_classes_chr %>%
  count(Buffering.ChrArmLevel.Class, CNA) %>%
  group_by(CNA) %>%
  mutate(Share = (n / sum(n)) * 100) %>%
  ungroup() %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Class, y = Share, x = CNA) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  labs(x = "Chromosome Arm", y = "Fraction", fill = NULL) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(0,0,0,0, unit = 'cm')) +
  ggtitle("EGFR") # Alternative: As fill label

### Test if fraction of scaling genes is increased upon gain
egfr_classes_chr %>%
  mutate(Scaling = Buffering.ChrArmLevel.Class == "Scaling") %>%
  count(CNA, Scaling) %>%
  pivot_wider(names_from = "Scaling", values_from = "n") %>%
  tibble::column_to_rownames("CNA") %>%
  as.matrix() %>%
  fisher.test()

## All GSEA comparisons
panel_gsea_all <- bind_rows(gsea_all %>% mutate(Comparison = "Buffering"),
                            gsea_all_as %>% mutate(Comparison = "Aneuploidy")) %>%
  ### Estimate difference between Buffering and Aneuploidy enrichment scores
  mutate(Score = mean(NES[Comparison == "Buffering"]) - mean(NES[Comparison == "Aneuploidy"]), .by = pathway) %>%
  mutate(Label = map_signif(padj),
         pathway = str_to_lower(str_replace_all(pathway, c("HALLMARK_" = "", "_" = " "))),
         pathway = fct_reorder(pathway, Score, .desc = FALSE)) %>%
  ### Remove pathways that are insignificant everywhere
  add_count(pathway, Label == "n.s.") %>%
  filter(n < 6) %>%
  ### Create plot
  distinct(pathway, Dataset, Comparison, NES, Label) %>%
  ggplot() +
  aes(x = pathway, y = "", fill = NES, label = Label) +
  geom_raster() +
  geom_text(color = "black") +
  ggh4x::facet_nested(Comparison + Dataset ~ ., switch = "y") +
  scale_fill_gradientn(colors = bidirectional_color_pal,
                       space = "Lab", limits = c(-3, 3), oob = scales::squish) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Cancer Hallmark Pathway", y = "", fill = "Normalized Enrichment Score") +
  cowplot::theme_minimal_grid() +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_rect(color = "lightgrey"),
        axis.ticks.y = element_blank(),
        legend.key.width = unit(24, "points"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(hjust = 1, vjust = 0.9, size = base_size),
        legend.text = element_text(size = 12),
        text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 15, 5, 15), "mm"))

## Combine figures
figure_s5_sub1 <- cowplot::plot_grid(ora_down_control, ora_up_control, panel_upr_purity,
                                     ncol = 3, rel_widths = c(1, 1, 1),
                                     labels = c("B", "C", "D"))
figure_s5_sub2 <- cowplot::plot_grid(panel_og_tsg, panel_oncogenes, panel_egfr,
                                     ncol = 3, rel_widths = c(0.6, 1, 0.6),
                                     labels = c("E", "F", "G"))
figure_s5 <- cowplot::plot_grid(panel_volcano_all, figure_s5_sub1, figure_s5_sub2, panel_gsea_all,
                                nrow = 4, rel_heights = c(0.95, 1, 0.9, 1.35),
                                labels = c("A", "", "", "H"))

cairo_pdf(here(plots_dir, "figure_s5.pdf"), width = 12, height = 18)
figure_s5
dev.off()
