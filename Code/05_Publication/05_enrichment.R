library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(STRINGdb)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

plots_dir <- here(plots_base_dir, "Publication")
temp_dir <- temp_base_dir
output_data_dir <- output_data_base_dir
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")

dir.create(plots_dir, recursive = TRUE)

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
df_celllines <- read_parquet(here(output_data_dir, 'celllines.parquet'))
expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))
metadata_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))


diff_exp_cptac <- read_parquet(here(output_data_dir, "model_buf_diff-exp_cptac.parquet"))
diff_exp_procan <- read_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))

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
  labs(y = "-log10(p.adj)", x = "Log2FC (High - Low Buffering)")

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
  plot_terms_compact(selected_sources = c("REAC", "GO:MF", "CORUM"), custom_color = color_palettes$DiffExp["Up"])

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  plot_terms_compact(selected_sources = c("KEGG", "GO:MF", "CORUM"), custom_color = color_palettes$DiffExp["Down"])

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
hallmark_gene_set <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  rename(Gene.Symbol = "gene_symbol")

gene_sets <- hallmark_gene_set %>%
  group_by(gs_name) %>%
  group_map(~list(.x$Gene.Symbol)) %>%
  purrr::list_flatten()
names(gene_sets) <- unique(hallmark_gene_set$gs_name)

fgsea_procan <- diff_exp_procan %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol) %>%
  fgsea::fgsea(pathways = gene_sets, stats = .)

fgsea_cptac <- diff_exp_cptac %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol) %>%
  fgsea::fgsea(pathways = gene_sets, stats = .)

max_abs_nes_procan <- round(max(abs(fgsea_procan$NES)))
max_abs_nes_cptac <- round(max(abs(fgsea_cptac$NES)))

selected_pathways <- c("DNA_REPAIR", "APOPTOSIS", "G2M_CHECKPOINT", "MYC_TARGETS_V2", "IL2_STAT5_SIGNALING",
                       "UNFOLDED_PROTEIN_RESPONSE", "EPITHELIAL_MESENCHYMAL_TRANSITION", "APICAL_JUNCTION")

panel_2d_enrichment <- bind_rows(fgsea_procan %>% mutate(Dataset = "ProCan"),
                                  fgsea_cptac %>% mutate(Dataset = "CPTAC")) %>%
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
  geom_point() +
  geom_label_repel(force = 30, min.segment.length = 0.01) +
  scale_color_manual(values = c(Both = highlight_colors[1], color_palettes$Datasets,
                                None = default_color)) +
  lims(x = c(-max_abs_nes_procan, max_abs_nes_procan), y = c(-max_abs_nes_cptac, max_abs_nes_cptac)) +
  labs(x = "Normalized Enrichment Score (ProCan)",
       y = "Normalized Enrichment Score (CPTAC)") +
  theme(legend.position = "top", legend.direction = "horizontal")

# === Proteotoxic Stress Panel ===
gsea_cptac_all <- expr_buf_cptac %>%
  ssGSEA(gene_sets, Protein.Expression.Normalized, Model.ID) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model.ID") %>%
  inner_join(y = model_buf_cptac, by = "Model.ID") %>%
  left_join(y = metadata_cptac %>% select(-HALLMARK_UNFOLDED_PROTEIN_RESPONSE), by = "Model.ID")


quantile(gsea_cptac_all$Observations, probs = seq(0,1,0.1))

gsea_cptac <- gsea_cptac_all %>%
  filter(Observations > 500)

cor_proteotox_test <- cor.test(gsea_cptac$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                               gsea_cptac$Model.Buffering.Ratio, method = "spearman")

panel_proteotox <- gsea_cptac %>%
  ggplot() +
  aes(x = Model.Buffering.Ratio, y = HALLMARK_UNFOLDED_PROTEIN_RESPONSE, color = Model.AneuploidyScore.Estimate) +
  geom_point(size = 2) +
  stat_smooth(method = lm, color = "dimgrey") +
  annotate("text", x = 0.35, y = 0.275, hjust = 0, size = 5,
           label = paste0(print_corr(cor_proteotox_test$estimate), ", ", print_signif(cor_proteotox_test$p.value))) +
  labs(x = "Model Buffering Ratio", y = "Unfolded Protein Response", color = "eAS") +
  scale_color_viridis_c(option = color_palettes$AneuploidyScore, end = 0.8) +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0, 0),
        legend.justification = c("left", "bottom"),
        legend.text = element_text(size = base_size),
        legend.background = element_rect(fill = alpha('white', 2/3)))

# === Proliferation ===
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
  annotate("text", x = 0, y = 14, hjust = 0, size = 5,
           label = paste0(print_corr(cor_prolif$estimate), ", ", print_signif(cor_prolif$p.value))) +
  # stat_cor(aes(color = NULL), method = "spearman", show.legend = FALSE, p.accuracy = 0.001, r.accuracy = 0.001, size = 5, cor.coef.name = "rho") +
  scale_color_manual(values = color_palettes$WGD) +
  labs(x = "Mean Normalized Buffering Ranks", y = "Growth Ratio (Day4/Day1)") +
  theme(legend.position = c("left", "top"),
        legend.position.inside = c(0.01, 0.90),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        legend.background = element_rect(fill = alpha('white', 2/3)))

# === Combine Panels into Figure ===
figure5_01 <- cowplot::plot_grid(panel_volcano, ora_down, ora_up,
                              labels = c("A", "B", "C"),
                              rel_widths = c(1, 1, 1),
                              nrow = 1, ncol = 3)

cairo_pdf(here(plots_dir, "figure05_01.pdf"), width = 13, height = 5)
figure5_01
dev.off()

# D & E reserved for STRING clusters

fig5_02_sub1 <- cowplot::plot_grid(panel_proteotox, panel_prolif_base,
                                   nrow = 2, ncol = 1, labels = c("G", "H"))

figure5_02 <- cowplot::plot_grid(panel_2d_enrichment, fig5_02_sub1,
                              labels = c("F", ""),
                              rel_widths = c(1, 0.7),
                              nrow = 1, ncol = 2)

cairo_pdf(here(plots_dir, "figure05_02.pdf"), width = 12, height = 8)
figure5_02
dev.off()
