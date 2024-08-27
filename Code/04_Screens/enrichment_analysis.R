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

expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac_pure.parquet"))
model_buf_cptac <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_cptac.parquet"))

# === Identify proteins with significant expression differences between cell lines with high and low buffering ===

diff_exp <- model_buf_procan %>%
  split_by_3_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  filter(Model.Buffering.Group != "Center") %>%
  inner_join(y = expr_buf_procan, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Model.Buffering.Group, Gene.Symbol, Protein.Expression.Normalized) %>%
  differential_expression(Gene.Symbol, Model.Buffering.Group, Protein.Expression.Normalized,
                          groups = c("Low", "High")) %>%
  write_parquet(here(output_data_dir, "model_buf_diff-exp_procan.parquet"))

## Volcano Plots
color_mapping <- scale_color_manual(values = c(Down = bidirectional_color_pal[1],
                                               Up = bidirectional_color_pal[5]),
                                    na.value = color_palettes$Missing)

volcano_plot <- diff_exp %>%
  mutate(Regulation = case_when(Log2FC > log2fc_threshold & TTest.p.adj < p_threshold ~ "Up",
                                Log2FC < -log2fc_threshold & TTest.p.adj < p_threshold ~ "Down",
                                TRUE ~ NA),
         Label = if_else(!is.na(Regulation), Gene.Symbol, NA)) %>%
  plot_volcano(Log2FC, TTest.p.adj, Label, Regulation, color_mapping) %>%
  save_plot("volcano_plot.png")

# === Enrichment Analysis ===
ref_gene <- unique(diff_exp$Gene.Symbol)

genes_up <- diff_exp %>%
  filter(Log2FC > log2fc_threshold & TTest.p.adj < p_threshold) %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = TRUE)) %>%
  arrange(Gene.Symbol)

ora_up <- genes_up %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_up_wg <- genes_up %>%
  pull(Gene.Symbol) %>%
  ora_webgestalt(ref_background = ref_gene)

genes_down <- diff_exp %>%
  filter(Log2FC < -log2fc_threshold & TTest.p.adj < p_threshold) %>%
  mutate(Gene.Symbol = fct_reorder(Gene.Symbol, Log2FC, .desc = FALSE)) %>%
  arrange(Gene.Symbol)

ora_down <- genes_down %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_down_wg <- genes_down %>%
  pull(Gene.Symbol) %>%
  ora_webgestalt(ref_background = ref_gene)

#gostplot(ora_up, capped = TRUE, interactive = TRUE)
#gostplot(ora_down, capped = TRUE, interactive = TRUE)

# TODO: Use WebGestaltR::affinityPropagation() for reducing redundant terms
# ora_down_reduced <- WebGestaltR::affinityPropagation(ora_down$result$term_name, -log10(ora_down$result$p_value))

 ora_up %>%
   plot_terms() %>%
   save_plot("overrepresentation_terms_up.png", height = 300)

 ora_down %>%
   plot_terms() %>%
   save_plot("overrepresentation_terms_down.png", height = 300)

# === Network Analysis ===
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                          network_type = "full", input_directory = temp_dir)

string_up <- genes_up %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_up_png <- string_db$get_png(string_up$df_mapped$STRING_id, payload_id = string_up$payload,
                                   required_score = 700, file = here(plots_dir, "string_up.png"))

string_down <- genes_down %>%
  create_string_network(Gene.Symbol, Log2FC, string_db)

string_down_png <- string_db$get_png(string_down$df_mapped$STRING_id, payload_id = string_down$payload,
                                     required_score = 700, file = here(plots_dir, "string_down.png"))

# ===  Karyotype & Expression by ChrArm ===
# TODO: Compare Karyotype & Expression per Chromosome between buffered and scaling cell lines
gene_metadata <- expr_buf_procan %>%
  distinct(Gene.Symbol, Gene.Chromosome, Gene.ChromosomeArm, Gene.StartPosition, Gene.EndPosition) %>%
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome))

# TODO: Nested facet with chromosome arm
diff_exp %>%
  inner_join(y = gene_metadata, by = "Gene.Symbol", relationship = "one-to-one") %>%
  bucketed_scatter_plot(Log2FC, Gene.StartPosition, Gene.Chromosome,
                        x_lab = "Chromosome & Gene Position")

# TODO: Check copy number / chromosome CNA between two cohorts

egfr_dc <- expr_buf_procan %>%
  filter(Gene.Symbol == "EGFR") %>%
  mutate(CNDiff = ChromosomeArm.CopyNumber - ChromosomeArm.CopyNumber.Baseline) %>%
  filter(abs(CNDiff) > 0.5) %>%
  mutate(CNV = if_else(CNDiff > 0, "Gain", "Loss")) %>%
  drop_na(Buffering.ChrArmLevel.Ratio) %>%
  signif_beeswarm_plot(CNV, Buffering.ChrArmLevel.Ratio, color_col = CellLine.AneuploidyScore) %>%
  save_plot("EGFR_DC_ChrArm.png", height = 150, width = 150)


# Gene Sets
## Proteotoxic Stress / Unfolded Proteins / Autophagosome
library(msigdbr)
hallmark_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>%
  rename(Gene.Symbol = "gene_symbol")

expr_buf_procan_hallmark <- model_buf_procan %>%
  split_by_3_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group") %>%
  filter(Model.Buffering.Group != "Center") %>%
  inner_join(y = expr_buf_procan, by = "Model.ID", relationship = "one-to-many", na_matches = "never") %>%
  select(Gene.Symbol, Model.ID, Model.Buffering.Group, Protein.Expression.Normalized) %>%
  right_join(y = hallmark_gene_set, by = "Gene.Symbol", relationship = "many-to-many")

hallmark_tests <- expr_buf_procan_hallmark %>%
  group_by(gs_name) %>%
  group_modify(~tidy(t.test(Protein.Expression.Normalized ~ Model.Buffering.Group, data = .x))) %>%
  ungroup() %>%
  mutate(p.value.adj = p.adjust(p.value, method = "bonferroni"))

# === TEST AREA ===
gene_sets <- hallmark_gene_set %>%
  group_by(gs_name) %>%
  group_map(~list(.x$Gene.Symbol)) %>%
  purrr::list_flatten()
names(gene_sets) <- unique(hallmark_gene_set$gs_name)

gsea_cptac <- expr_buf_cptac %>%
  ssGSEA(gene_sets, Protein.Expression.Normalized, Model.ID)

gsea_cptac %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model.ID") %>%
  inner_join(y = model_buf_cptac, by = "Model.ID") %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, HALLMARK_UNFOLDED_PROTEIN_RESPONSE)

gsea_procan <- expr_buf_procan %>%
  ssGSEA(gene_sets, Protein.Expression.Normalized, Model.ID)

gsea_procan %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model.ID") %>%
  inner_join(y = model_buf_procan, by = "Model.ID") %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, HALLMARK_UNFOLDED_PROTEIN_RESPONSE)

unfolded_gene_set <- hallmark_gene_set %>%
  filter(gs_name == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE") %>%
  pull(Gene.Symbol)

### ProCan
mean_unfolded_procan <- expr_buf_procan %>%
  # filter_cn_gain_abs() %>%
  drop_na(Protein.Expression.Normalized) %>%
  mutate(Unfolded.Gene.Set = Gene.Symbol %in% unfolded_gene_set) %>%
  group_by(Model.ID) %>%
  summarize(Unfolded.Protein.Mean = mean(Protein.Expression.Normalized[Unfolded.Gene.Set], na.rm = TRUE),
            Background.Protein.Mean = mean(Protein.Expression.Normalized[!Unfolded.Gene.Set], na.rm = TRUE),
            Unfolded.Protein.Response = Unfolded.Protein.Mean - Background.Protein.Mean) %>%
  ungroup() %>%
  inner_join(y = model_buf_procan, by = "Model.ID", relationship = "one-to-many", na_matches = "never")

mean_unfolded_procan %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_scatter_procan.png", height = 150, width = 150)

mean_unfolded_procan %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                     quantile_low = "30%", quantile_high = "70%") %>%
  signif_violin_plot(Model.Buffering.Group, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_split_procan.png", width = 100, height = 150)

### CPTAC
mean_unfolded_cptac <- expr_buf_cptac %>%
  # filter_cn_gain_abs() %>%
  drop_na(Protein.Expression.Normalized) %>%
  mutate(Unfolded.Gene.Set = Gene.Symbol %in% unfolded_gene_set) %>%
  group_by(Model.ID) %>%
  summarize(Unfolded.Protein.Mean = mean(Protein.Expression.Normalized[Unfolded.Gene.Set], na.rm = TRUE),
            Background.Protein.Mean = mean(Protein.Expression.Normalized[!Unfolded.Gene.Set], na.rm = TRUE),
            Unfolded.Protein.Response = Unfolded.Protein.Mean - Background.Protein.Mean) %>%
  ungroup() %>%
  inner_join(y = model_buf_cptac, by = "Model.ID", relationship = "one-to-many", na_matches = "never")

mean_unfolded_cptac %>%
  scatter_plot_reg_corr(Model.Buffering.Ratio, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_scatter_cptac.png", height = 150, width = 150)

mean_unfolded_cptac %>%
  split_by_quantiles(Model.Buffering.Ratio, target_group_col = "Model.Buffering.Group",
                     quantile_low = "30%", quantile_high = "70%") %>%
  signif_violin_plot(Model.Buffering.Group, Unfolded.Protein.Response) %>%
  save_plot("unfolded_protein_response_split_cptac.png", width = 100, height = 150)

