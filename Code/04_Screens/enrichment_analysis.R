library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)
library(openxlsx)
library(gprofiler2)
library(WebGestaltR)
library(STRINGdb)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Screens", "Enrichment")
reports_dir <- reports_base_dir
temp_dir <- temp_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)
dir.create(temp_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
metadata_procan <- read_csv_arrow(here(external_data_dir, "CopyNumber", "ProCan", "model_list_20240110.csv")) %>%
  rename(CellLine.SangerModelId = "model_id")

expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
model_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))
metadata_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

cancer_genes <- read_parquet(here(output_data_dir, "cancer_genes.parquet"))

# === Identify proteins with significant expression differences between cell lines with high and low buffering ===
diff_exp_procan <- model_buf_procan %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_procan, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))

diff_exp_depmap <- model_buf_depmap %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_depmap, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "model_buf_diff-exp_depmap.parquet"))

diff_exp_cptac <- model_buf_cptac %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_cptac, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "model_buf_diff-exp_cptac.parquet"))

## Volcano Plots
color_mapping <- scale_color_manual(values = color_palettes$DiffExp,
                                    na.value = color_palettes$Missing)

volcano_plot_procan <- diff_exp_procan %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping) %>%
  save_plot("volcano_plot_procan.png")

volcano_plot_depmap <- diff_exp_depmap %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping) %>%
  save_plot("volcano_plot_depmap.png")

volcano_plot_cptac <- diff_exp_cptac %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping) %>%
  save_plot("volcano_plot_cptac.png")

## Skip analysis of CPTAC - no clear enrichment

# === Enrichment Analysis ===
genes_up_procan <- diff_exp_procan %>%
  filter(Significant == "Up") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

genes_down_procan <- diff_exp_procan %>%
  filter(Significant == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_up <- genes_up_procan %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_down <- genes_down_procan %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

# TODO: Use WebGestaltR::affinityPropagation() for reducing redundant terms
# ora_down_reduced <- WebGestaltR::affinityPropagation(ora_down$result$term_name, -log10(ora_down$result$p_value))

 ora_up %>%
   plot_terms() %>%
   save_plot("overrepresentation_terms_up.png", height = 300)

 ora_down %>%
   plot_terms() %>%
   save_plot("overrepresentation_terms_down.png", height = 300)

# === Common Enriched Genes ===
genes_up_common <- diff_exp_depmap %>%
  filter(Significant == "Up") %>%
  pull(Gene.Symbol) %>%
  intersect(genes_up_procan$Gene.Symbol)

genes_down_common <- diff_exp_depmap %>%
  filter(Significant == "Down") %>%
  pull(Gene.Symbol) %>%
  intersect(genes_down_procan$Gene.Symbol)

ora_up_common <- genes_up_common %>%
  overrepresentation_analysis() %>%
  plot_terms() %>%
  save_plot("overrepresentation_terms_up_common.png", height = 300)

ora_down_common <- genes_down_common %>%
  overrepresentation_analysis() %>%
  plot_terms() %>%
  save_plot("overrepresentation_terms_down_common.png", height = 300)

# === Network Analysis ===
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                          network_type = "full", input_directory = temp_dir)

string_up <- genes_up_procan %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_up_png <- string_db$get_png(string_up$df_mapped$STRING_id, payload_id = string_up$payload,
                                   required_score = 700, file = here(plots_dir, "string_up.png"))

string_down <- genes_down_procan %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_down_png <- string_db$get_png(string_down$df_mapped$STRING_id, payload_id = string_down$payload,
                                     required_score = 700, file = here(plots_dir, "string_down.png"))

## Common genes
string_up_common <- data.frame(Gene.Symbol = genes_up_common, Log2FC = 0) %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_up_common_png <- string_db$get_png(string_up_common$df_mapped$STRING_id, payload_id = string_up_common$payload,
                                   required_score = 700, file = here(plots_dir, "string_up_common.png"))

string_down_common <- data.frame(Gene.Symbol = genes_down_common, Log2FC = 0) %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_down_common_png <- string_db$get_png(string_down_common$df_mapped$STRING_id, payload_id = string_down_common$payload,
                                   required_score = 700, file = here(plots_dir, "string_down_common.png"))

# ===  Karyotype & Expression by ChrArm between High & Low buffering groups ===
gene_metadata <- expr_buf_procan %>%
  distinct(Gene.Symbol, Gene.Chromosome, Gene.ChromosomeArm, Gene.StartPosition, Gene.EndPosition) %>%
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome))

diff_exp_procan %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position") %>%
  save_plot("diffexp_by_chromosome_procan.png", height = 150, width = 200)

diff_exp_depmap %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position") %>%
  save_plot("diffexp_by_chromosome_depmap.png", height = 150, width = 200)

diff_exp_cptac %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position") %>%
  save_plot("diffexp_by_chromosome_cptac.png", height = 150, width = 200)

# TODO: Nested facet with chromosome arm

## Copy Number Differences
diff_exp_cn_procan <- model_buf_procan %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_procan, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Gene.CopyNumber) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Gene.CopyNumber,
                          groups = c("Low", "High"))

diff_exp_cn_depmap <- model_buf_depmap %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_depmap, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Gene.CopyNumber) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Gene.CopyNumber,
                          groups = c("Low", "High"))

diff_exp_cn_cptac <- model_buf_cptac %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_cptac, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Gene.CopyNumber) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Gene.CopyNumber,
                          groups = c("Low", "High"))

diff_exp_cn_procan %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position", value_range = c(0, 3)) %>%
  save_plot("diffexp_by_chromosome_procan_copy-number.png", height = 150, width = 200)

diff_exp_cn_depmap %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position", value_range = c(0, 3)) %>%
  save_plot("diffexp_by_chromosome_depmap_copy-number.png", height = 150, width = 200)

diff_exp_cn_cptac %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position", value_range = c(0, 4)) %>%
  save_plot("diffexp_by_chromosome_cptac_copy-number.png", height = 150, width = 200)

# === Cancer Driver Genes Buffering ===
og_dc <- expr_buf_procan %>%
  inner_join(y = cancer_genes, by = "Gene.Symbol") %>%
  filter(CancerDriverMode == "OG") %>%
  filter(Gene.Symbol %in% genes_up_procan$Gene.Symbol) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.GeneLevel.Ratio) %>%
  signif_boxplot(CNV, Buffering.GeneLevel.Ratio) %>%
  save_plot("oncogene_upregulated_buffering.png")

tsg_dc <- expr_buf_procan %>%
  inner_join(y = cancer_genes, by = "Gene.Symbol") %>%
  filter(CancerDriverMode == "TSG") %>%
  filter(Gene.Symbol %in% genes_down_procan$Gene.Symbol) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.GeneLevel.Ratio) %>%
  signif_boxplot(CNV, Buffering.GeneLevel.Ratio) %>%
  save_plot("tumorsupressor_downregulated_buffering.png")

## Common upregulated genes
og_common_dc <- expr_buf_procan %>%
  inner_join(y = cancer_genes, by = "Gene.Symbol") %>%
  filter(CancerDriverMode == "OG") %>%
  filter(Gene.Symbol %in% genes_up_common) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.GeneLevel.Ratio) %>%
  signif_boxplot(CNV, Buffering.GeneLevel.Ratio, facet_col = Gene.Symbol) %>%
  save_plot("oncogene_common_upregulated_buffering.png", width = 400)

og_common_dc <- expr_buf_procan %>%
  inner_join(y = cancer_genes, by = "Gene.Symbol") %>%
  filter(CancerDriverMode == "OG") %>%
  filter(Gene.Symbol %in% genes_up_common) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.GeneLevel.Ratio) %>%
  signif_beeswarm_plot(CNV, Buffering.GeneLevel.Ratio, facet_col = Gene.Symbol, color_col = CellLine.AneuploidyScore) %>%
  save_plot("oncogene_common_upregulated_buffering.png", width = 400)

## EGFR & RRAS (chr arm)
og_common_dc_chr <- expr_buf_procan %>%
  filter(Gene.Symbol %in% c("EGFR", "RRAS")) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CNA = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Gain", "Loss")) %>%
  drop_na(Buffering.ChrArmLevel.Ratio) %>%
  signif_beeswarm_plot(CNA, Buffering.ChrArmLevel.Ratio, facet_col = Gene.Symbol, color_col = CellLine.AneuploidyScore) %>%
  save_plot("oncogene_buffering_selected_chr.png", width = 200)

## EGFR Scaling Frequency Chr Arm Gain & Loss
egfr_classes_chr <- expr_buf_procan %>%
  filter(Gene.Symbol == "EGFR") %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  drop_na(Buffering.ChrArmLevel.Ratio) %>%
  mutate(CNA = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Gain", "Loss"))

(egfr_classes_chr %>%
  count(Buffering.ChrArmLevel.Class, CNA) %>%
  group_by(CNA) %>%
  mutate(Share = (n / sum(n)) * 100) %>%
  ungroup() %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Class, y = Share, x = CNA) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palettes$BufferingClasses)) %>%
  save_plot("egfr_buffering_classes_chr.png", width = 200)

egfr_classes_chr %>%
  mutate(Scaling = Buffering.ChrArmLevel.Class == "Scaling") %>%
  count(CNA, Scaling) %>%
  pivot_wider(names_from = "Scaling", values_from = "n") %>%
  tibble::column_to_rownames("CNA") %>%
  as.matrix() %>%
  fisher.test()

# === Gene Set Enrichment Analysis ===
## All Hallmark Gene Sets (ssGSEA)
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
  write_parquet(here(output_data_dir, "ssgsea_unfolded_cptac.parquet"))

gsea_procan <- expr_buf_procan %>%
  ssGSEA(gene_sets, Protein.Expression.Normalized, Model.ID) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model.ID") %>%
  write_parquet(here(output_data_dir, "ssgsea_unfolded_procan.parquet"))

gsea_depmap <- expr_buf_depmap %>%
  ssGSEA(gene_sets, Protein.Expression.Normalized, Model.ID) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model.ID") %>%
  write_parquet(here(output_data_dir, "ssgsea_unfolded_depmap.parquet"))

## Proteotoxic Stress / Unfolded Proteins / Autophagosome
gsea_cptac %>%
  inner_join(y = model_buf_cptac, by = "Model.ID") %>%
  left_join(y = metadata_cptac %>% select(-HALLMARK_UNFOLDED_PROTEIN_RESPONSE), by = "Model.ID") %>%
  filter(Observations > 500) %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                        color_col = Model.AneuploidyScore.Estimate) %>%
  save_plot("gsea_unfolded_cptac.png")

gsea_procan %>%
  inner_join(y = model_buf_procan, by = "Model.ID") %>%
  filter(Observations > 500) %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, HALLMARK_UNFOLDED_PROTEIN_RESPONSE) %>%
  save_plot("gsea_unfolded_procan.png")

gsea_depmap %>%
  inner_join(y = model_buf_procan, by = "Model.ID") %>%
  filter(Observations > 500) %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, HALLMARK_UNFOLDED_PROTEIN_RESPONSE) %>%
  save_plot("gsea_unfolded_depmap.png")

## Manual Analysis of Unfolded Protein Response
unfolded_gene_set <- hallmark_gene_set %>%
  filter(gs_name == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE") %>%
  pull(Gene.Symbol)

unfolded_protein_response <- function(df_expr_buf, df_model_buf, gene_set = unfolded_gene_set) {
 df_expr_buf %>%
   # filter_cn_gain_abs() %>%
   drop_na(Protein.Expression.Normalized) %>%
   mutate(Unfolded.Gene.Set = Gene.Symbol %in% gene_set) %>%
   group_by(Model.ID) %>%
   summarize(Unfolded.Protein.Mean = mean(Protein.Expression.Normalized[Unfolded.Gene.Set], na.rm = TRUE),
             Background.Protein.Mean = mean(Protein.Expression.Normalized[!Unfolded.Gene.Set], na.rm = TRUE),
             Unfolded.Protein.Response = Unfolded.Protein.Mean - Background.Protein.Mean) %>%
   ungroup() %>%
   inner_join(y = df_model_buf, by = "Model.ID", relationship = "one-to-many", na_matches = "never")
}

### CPTAC
mean_unfolded_cptac <- unfolded_protein_response(expr_buf_cptac, model_buf_cptac)

mean_unfolded_cptac %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_scatter_cptac.png", height = 150, width = 150)

mean_unfolded_cptac %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                     quantile_low = "30%", quantile_high = "70%") %>%
  signif_violin_plot(Model.Buffering.Group, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_split_cptac.png", width = 100, height = 150)

### ProCan
mean_unfolded_procan <- unfolded_protein_response(expr_buf_procan, model_buf_procan)

mean_unfolded_procan %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_scatter_procan.png", height = 150, width = 150)

mean_unfolded_procan %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                     quantile_low = "30%", quantile_high = "70%") %>%
  signif_violin_plot(Model.Buffering.Group, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_split_procan.png", width = 100, height = 150)

### DepMap
mean_unfolded_depmap <- unfolded_protein_response(expr_buf_depmap, model_buf_depmap)

mean_unfolded_depmap %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_scatter_depmap.png", height = 150, width = 150)

mean_unfolded_depmap %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                     quantile_low = "30%", quantile_high = "70%") %>%
  signif_violin_plot(Model.Buffering.Group, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_split_depmap.png", width = 100, height = 150)


## 2D-GSEA
ranks_depmap <- diff_exp_depmap %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol)

ranks_procan <- diff_exp_procan %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol)

ranks_cptac <- diff_exp_cptac %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol)

fgsea_depmap <- fgsea::fgsea(pathways = gene_sets, stats = ranks_depmap)
fgsea_procan <- fgsea::fgsea(pathways = gene_sets, stats = ranks_procan)
fgsea_cptac <- fgsea::fgsea(pathways = gene_sets, stats = ranks_cptac)

fgsea_results <- bind_rows(fgsea_depmap %>% mutate(Dataset = "DepMap"),
                           fgsea_procan %>% mutate(Dataset = "ProCan"),
                           fgsea_cptac %>% mutate(Dataset = "CPTAC")) %>%
  write_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer.parquet"))

### CPTAC vs ProCan
enrichment_2d_procan <- fgsea_results %>%
  filter(Dataset %in% c("CPTAC", "ProCan")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_ProCan < p_threshold & padj_CPTAC < p_threshold ~ "Both",
    padj_ProCan >= p_threshold & padj_CPTAC < p_threshold ~ "CPTAC",
    padj_ProCan < p_threshold & padj_CPTAC >= p_threshold ~ "ProCan",
    TRUE ~ "None"),
         Label = str_replace(pathway, "HALLMARK_", "")) %>%
  ggplot() +
  aes(x = NES_ProCan, y = NES_CPTAC, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(slope = 1, color = default_color) +
  geom_point() +
  geom_label_repel() +
  scale_color_manual(values = c(Both = "purple", CPTAC = "blue", ProCan = "red", None = default_color)) +
  labs(x = "Normalized Enrichment Score (ProCan, High vs. Low Model Buffering)",
       y = "Normalized Enrichment Score (CPTAC, High vs. Low Model Buffering)")

enrichment_2d_procan %>%
  save_plot("gsea_2d_hallmark_cptac_procan.png", height = 250, width = 300)

### CPTAC vs DepMap
enrichment_2d_depmap <- fgsea_results %>%
  filter(Dataset %in% c("CPTAC", "DepMap")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_DepMap < p_threshold & padj_CPTAC < p_threshold ~ "Both",
    padj_DepMap >= p_threshold & padj_CPTAC < p_threshold ~ "CPTAC",
    padj_DepMap < p_threshold & padj_CPTAC >= p_threshold ~ "DepMap",
    TRUE ~ "None"),
         Label = str_replace(pathway, "HALLMARK_", "")) %>%
  ggplot() +
  aes(x = NES_DepMap, y = NES_CPTAC, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(slope = 1, color = default_color) +
  geom_point() +
  geom_label_repel() +
  scale_color_manual(values = c(Both = "purple", CPTAC = "blue", DepMap = "red", None = default_color)) +
  labs(x = "Normalized Enrichment Score (DepMap, High vs. Low Model Buffering)",
       y = "Normalized Enrichment Score (CPTAC, High vs. Low Model Buffering)")

enrichment_2d_depmap %>%
  save_plot("gsea_2d_hallmark_cptac_depmap.png", height = 250, width = 300)

### DepMap vs ProCan
enrichment_2d_procan_depmap <- fgsea_results %>%
  filter(Dataset %in% c("DepMap", "ProCan")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_ProCan < p_threshold & padj_DepMap < p_threshold ~ "Both",
    padj_ProCan >= p_threshold & padj_DepMap < p_threshold ~ "DepMap",
    padj_ProCan < p_threshold & padj_DepMap >= p_threshold ~ "ProCan",
    TRUE ~ "None"),
         Label = str_replace(pathway, "HALLMARK_", "")) %>%
  ggplot() +
  aes(x = NES_ProCan, y = NES_DepMap, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(slope = 1, color = default_color) +
  geom_point() +
  geom_label_repel() +
  scale_color_manual(values = c(Both = "purple", DepMap = "blue", ProCan = "red", None = default_color)) +
  labs(x = "Normalized Enrichment Score (ProCan, High vs. Low Model Buffering)",
       y = "Normalized Enrichment Score (DepMap, High vs. Low Model Buffering)")

enrichment_2d_procan_depmap %>%
  save_plot("gsea_2d_hallmark_procan_depmap.png", height = 250, width = 300)

# === Split datasets by aneuploidy score as control and perform DiffExp & GSEA ===
## DiffExp
diff_exp_as_procan <- expr_buf_procan %>%
  split_by_quantiles(CellLine.AneuploidyScore, target_group_col = "Aneuploidy") %>%
  select(Aneuploidy, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Aneuploidy, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "aneuploidy_diff-exp_procan.parquet"))

diff_exp_as_depmap <- expr_buf_depmap %>%
  split_by_quantiles(CellLine.AneuploidyScore, target_group_col = "Aneuploidy") %>%
  select(Aneuploidy, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Aneuploidy, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "aneuploidy_diff-exp_depmap.parquet"))

diff_exp_as_cptac <- expr_buf_cptac %>%
  inner_join(y = metadata_cptac %>% select(Model.ID, Model.AneuploidyScore.Estimate),
             by = "Model.ID", relationship = "many-to-one", na_matches = "never") %>%
  split_by_quantiles(Model.AneuploidyScore.Estimate, target_group_col = "Aneuploidy") %>%
  select(Aneuploidy, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Aneuploidy, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "aneuploidy_diff-exp_cptac.parquet"))

## GSEA
ranks_as_depmap <- diff_exp_as_depmap %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol)

ranks_as_procan <- diff_exp_as_procan %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol)

ranks_as_cptac <- diff_exp_as_cptac %>%
  arrange(Log2FC) %>%
  pull(Log2FC, name = Gene.Symbol)

fgsea_as_depmap <- fgsea::fgsea(pathways = gene_sets, stats = ranks_as_depmap)
fgsea_as_procan <- fgsea::fgsea(pathways = gene_sets, stats = ranks_as_procan)
fgsea_as_cptac <- fgsea::fgsea(pathways = gene_sets, stats = ranks_as_cptac)

fgsea_results_as <- bind_rows(fgsea_as_depmap %>% mutate(Dataset = "DepMap"),
                           fgsea_as_procan %>% mutate(Dataset = "ProCan"),
                           fgsea_as_cptac %>% mutate(Dataset = "CPTAC")) %>%
  write_parquet(here(output_data_dir, "gsea_hallmark_pan-cancer_aneuploidy.parquet"))

### CPTAC vs. ProCan
enrichment_2d_as_cptac_procan <- fgsea_results_as %>%
  filter(Dataset %in% c("CPTAC", "ProCan")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_ProCan < p_threshold & padj_CPTAC < p_threshold ~ "Both",
    padj_ProCan >= p_threshold & padj_CPTAC < p_threshold ~ "CPTAC",
    padj_ProCan < p_threshold & padj_CPTAC >= p_threshold ~ "ProCan",
    TRUE ~ "None"),
         Label = str_replace(pathway, "HALLMARK_", "")) %>%
  ggplot() +
  aes(x = NES_ProCan, y = NES_CPTAC, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(slope = 1, color = default_color) +
  geom_point() +
  geom_label_repel() +
  scale_color_manual(values = c(Both = "purple", CPTAC = "blue", ProCan = "red", None = default_color)) +
  labs(x = "Normalized Enrichment Score (ProCan, High vs. Low Model Buffering)",
       y = "Normalized Enrichment Score (CPTAC, High vs. Low Model Buffering)")

enrichment_2d_as_cptac_procan %>%
  save_plot("gsea_2d_hallmark_aneuploidy_cptac_procan.png", height = 250, width = 300)

### CPTAC vs. DepMap
enrichment_2d_as_cptac_depmap <- fgsea_results_as %>%
  filter(Dataset %in% c("CPTAC", "DepMap")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_DepMap < p_threshold & padj_CPTAC < p_threshold ~ "Both",
    padj_DepMap >= p_threshold & padj_CPTAC < p_threshold ~ "CPTAC",
    padj_DepMap < p_threshold & padj_CPTAC >= p_threshold ~ "DepMap",
    TRUE ~ "None"),
         Label = str_replace(pathway, "HALLMARK_", "")) %>%
  ggplot() +
  aes(x = NES_DepMap, y = NES_CPTAC, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(slope = 1, color = default_color) +
  geom_point() +
  geom_label_repel() +
  scale_color_manual(values = c(Both = "purple", CPTAC = "blue", DepMap = "red", None = default_color)) +
  labs(x = "Normalized Enrichment Score (DepMap, High vs. Low Model Buffering)",
       y = "Normalized Enrichment Score (CPTAC, High vs. Low Model Buffering)")

enrichment_2d_as_cptac_depmap %>%
  save_plot("gsea_2d_hallmark_aneuploidy_cptac_depmap.png", height = 250, width = 300)

### DepMap vs. ProCan
enrichment_2d_as_procan_depmap <- fgsea_results_as %>%
  filter(Dataset %in% c("DepMap", "ProCan")) %>%
  select(pathway, Dataset, NES, padj) %>%
  pivot_wider(id_cols = "pathway", names_from = "Dataset", values_from = c("padj", "NES")) %>%
  mutate(Significant = case_when(
    padj_ProCan < p_threshold & padj_DepMap < p_threshold ~ "Both",
    padj_ProCan >= p_threshold & padj_DepMap < p_threshold ~ "DepMap",
    padj_ProCan < p_threshold & padj_DepMap >= p_threshold ~ "ProCan",
    TRUE ~ "None"),
         Label = str_replace(pathway, "HALLMARK_", "")) %>%
  ggplot() +
  aes(x = NES_ProCan, y = NES_DepMap, label = Label, color = Significant) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_vline(xintercept = 0, color = default_color) +
  geom_abline(slope = 1, color = default_color) +
  geom_point() +
  geom_label_repel() +
  scale_color_manual(values = c(Both = "purple", DepMap = "blue", ProCan = "red", None = default_color)) +
  labs(x = "Normalized Enrichment Score (ProCan, High vs. Low Model Buffering)",
       y = "Normalized Enrichment Score (DepMap, High vs. Low Model Buffering)")

enrichment_2d_as_procan_depmap %>%
  save_plot("gsea_2d_hallmark_aneuploidy_procan_depmap.png", height = 250, width = 300)

# === Control: Remove cells with suspension growth methods ===
# Rationale: Integrins detected in ORA, Blood cancers have low BR, Blood cancers are grown in suspension

## DiffExp
diff_exp_procan_adherent <- model_buf_procan %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  inner_join(y = expr_buf_procan, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = metadata_procan, by = "CellLine.SangerModelId") %>%
  filter(growth_properties == "Adherent") %>%
  select(Model.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "model_buf_diff-exp_procan_adherrent.parquet"))

volcano_plot_procan_adherent <- diff_exp_procan_adherent %>%
  left_join(y = cancer_genes, by = "Gene.Symbol") %>%
  mutate(Label = if_else(!is.na(Significant) & !is.na(CancerDriverMode), Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, color_mapping) %>%
  save_plot("volcano_plot_procan_adherrent.png")

## ORA
genes_up_procan_adherent <- diff_exp_procan_adherent %>%
  filter(Significant == "Up") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

genes_down_procan_adherent <- diff_exp_procan_adherent %>%
  filter(Significant == "Down") %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_up_adherent <- genes_up_procan_adherent %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  semi_join(x = .$result, y = ora_up$result, by = "term_id") # Filter for terms that were discovered without control

ora_down_adherent <- genes_down_procan_adherent %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis() %>%
  semi_join(x = .$result, y = ora_down$result, by = "term_id") # Filter for terms that were discovered without control

list(result = ora_up_adherent) %>%
  plot_terms() %>%
  save_plot("overrepresentation_terms_up_adherent.png", height = 300)

list(result = ora_down_adherent) %>%
  plot_terms() %>%
  save_plot("overrepresentation_terms_down_adherent.png", height = 300)

## Common genes
genes_up_common_adherent <- intersect(genes_up_common, genes_up_procan_adherent$Gene.Symbol)
genes_down_common_adherent <- intersect(genes_down_common, genes_down_procan_adherent$Gene.Symbol)

genes_up_common_adherent %>%
  overrepresentation_analysis() %>%
  plot_terms() %>%
  save_plot("overrepresentation_terms_up_common_adherent.png", height = 300)

genes_down_common_adherent %>%
  overrepresentation_analysis() %>%
  plot_terms() %>%
  save_plot("overrepresentation_terms_down_common_adherent.png", height = 300)
