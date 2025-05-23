library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))
source(here("Code", "preprocessing.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "publication.R"))

plots_dir <- here(plots_base_dir, "Publication")
tables_dir <- here(tables_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

# === Load Datasets ===
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))
expr_buf_p0211 <- read_parquet(here(output_data_dir, "expression_buffering_p0211.parquet"))
expr_buf_chunduri <- read_parquet(here(output_data_dir, "expression_buffering_chunduri.parquet"))

expr_buf_eng <- bind_rows(expr_buf_p0211, expr_buf_chunduri) %>%
  mutate(Dataset = "Engin.",
         Gene.Chromosome = as.character(Gene.Chromosome))

df_model_cptac <- read_parquet(here(output_data_dir, 'metadata_cptac.parquet'))

low_var_buf <- read_parquet(here(output_data_dir, "frequently_buffered_genes.parquet"))
low_var_buf_gain <- read_parquet(here(output_data_dir, "frequently_buffered_genes_gain.parquet"))
low_var_buf_loss <- read_parquet(here(output_data_dir, "frequently_buffered_genes_loss.parquet"))

# === DC Class Illustration Panel ===
dc_class_line <- expr_buf_depmap %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  filter(Gene.CopyNumber.Baseline == 2) %>%
  filter(Gene.CopyNumber %in% c(1, 3, 4)) %>%
  mutate(Gene.CopyNumber = as.factor(Gene.CopyNumber),
         Buffering.GeneLevel.Class = factor(Buffering.GeneLevel.Class,
                                            levels = c("Scaling", "Buffered", "Anti-Scaling"))) %>%
  drop_na(Buffering.GeneLevel.Class, Protein.Expression.Normalized) %>%
  group_by(Buffering.GeneLevel.Class, Gene.CopyNumber) %>%
  summarize(Mean.Protein.Ratio = (mean(2^Protein.Expression.Normalized / 2^Protein.Expression.Baseline) - 1),
            .groups = "drop") %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Class, y = Mean.Protein.Ratio,
      group = Gene.CopyNumber, color = Gene.CopyNumber, label = Gene.CopyNumber) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_path(linewidth = 1) +
  scale_y_continuous(limits = c(-0.6, 2.5),
                     breaks = seq(-0.5, 2.5, 0.5),
                     labels = paste0(seq(-0.5, 2.5, 0.5) * 100, "%")) +
  scale_color_manual(values = color_palettes$CopyNumbers) +
  labs(x = "Buffering Class", y = "Mean Protein Abundance", color = "Gene Copy Number") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm'))

# === Categorical Distribution Panel ===
df_share_gene <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  select(Buffering.GeneLevel.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.GeneLevel.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n))) %>%
  ungroup()

df_share_chr <- bind_rows(expr_buf_depmap, expr_buf_eng) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "ChrArm Gain", "ChrArm Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  select(Buffering.ChrArmLevel.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.ChrArmLevel.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n))) %>%
  ungroup()

stacked_buf_cn <- df_share_gene %>%
  ggplot() +
  aes(fill = Buffering.GeneLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = NULL, y = "Fraction") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm'))

stacked_buf_cn_chr <- df_share_chr %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = "Buffering Class") +
  theme(legend.position = "none")

# === DC Class Protein Expression Panel ===
prot_exp_dc_class <- expr_buf_depmap %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  drop_na(Buffering.GeneLevel.Class, Protein.Expression.Normalized) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss"),
         Buffering.GeneLevel.Class = factor(Buffering.GeneLevel.Class,
                                            levels = c("Scaling", "Buffered", "Anti-Scaling"))) %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Class, y = Protein.Expression.Normalized, color = Buffering.GeneLevel.Class) +
  geom_hline(yintercept = 0, color = default_color) +
  geom_hline(aes(yintercept = mean(Protein.Expression.Baseline, na.rm = TRUE)), color = highlight_colors[2]) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  stat_summary(aes(y = -2), fun.data = show.n, geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("Anti-Scaling", "Buffered"),
                                 c("Buffered", "Scaling"),
                                 c("Anti-Scaling", "Scaling")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(1.5, 1.5, 1.8),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  facet_grid(~Gene.CopyNumber.Event) +
  labs(x = "Buffering Class", y = "Protein Abundance") +
  scale_color_manual(values = color_palettes$BufferingClasses) +
  theme(legend.position = "none")

# === Tumor vs Cell Lines Panel ===
dc_dataset_dist <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter(abs(Gene.CopyNumber - Gene.CopyNumber.Baseline) > 0.2) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  select(Dataset, Buffering.GeneLevel.Ratio, Gene.CopyNumber.Event, Gene.Symbol) %>%
  drop_na() %>%
  ggplot() +
  aes(x = Dataset, y = Buffering.GeneLevel.Ratio, color = Dataset) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  stat_summary(aes(y = -1), fun.data = show.n, geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("DepMap", "ProCan"),
                                 c("CPTAC", "ProCan"),
                                 c("CPTAC", "DepMap")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(1.8, 1.8, 2.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  theme(legend.position = "") +
  scale_color_manual(values = color_palettes$Datasets) +
  facet_wrap(~Gene.CopyNumber.Event, ncol = 2) +
  labs(x = "Dataset", y = "Buffering Ratio")

# === Buffering Gain vs. Loss Panel ===
br_by_cna <- bind_rows(expr_buf_depmap, expr_buf_eng) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CopyNumber.Event = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Gain", "Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  group_by(Gene.Symbol, CopyNumber.Event, Dataset) %>%
  summarize(Buffering.ChrArmLevel.Ratio = mean(Buffering.ChrArmLevel.Ratio)) %>%
  ggplot() +
  aes(x = CopyNumber.Event, y = Buffering.ChrArmLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  stat_summary(aes(y = -0.3), fun.data = \(x) show.n(x, prefix = "n="),
               geom = "text", color = default_color, size = floor(base_size / 4)) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              map_signif_level = print_signif, y_position = 1.9, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  labs(x = "Chromosome Arm", y = "Mean Buffering Ratio") +
  ylim(c(-0.5, 2.5)) +
  facet_grid(~Dataset)

br_by_cnv <- bind_rows(expr_buf_cptac, expr_buf_depmap) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss"),
         Dataset = factor(Dataset, levels = dataset_order)) %>%
  group_by(Gene.Symbol, Gene.CopyNumber.Event, Dataset) %>%
  summarize(Buffering.GeneLevel.Ratio = mean(Buffering.GeneLevel.Ratio)) %>%
  ggplot() +
  aes(x = Gene.CopyNumber.Event, y = Buffering.GeneLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8, alpha = 0) +
  stat_summary(aes(y = -0.3), fun.data = \(x) show.n(x, prefix = "n="),
               geom = "text", color = default_color, size = floor(base_size / 4)) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              map_signif_level = print_signif, y_position = 1.9, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black") +
  labs(x = "Gene Copy Number", y = "Mean Buffering Ratio") +
  ylim(c(-0.5, 2.5)) +
  facet_grid(~Dataset)

# === Low Variance Buffering Panel ===
scatter_signif_buffered <- low_var_buf %>%
  mutate(Label = if_else(Top10, Gene.Symbol, NA)) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Top50) %>%
  ggplot() +
  aes(x = Gene.BR.SD, y = Gene.BR.Mean, label = Label, color = Top50, alpha = Top50) +
  geom_point() +
  geom_label_repel(na.rm = TRUE, max.overlaps = 50) +
  scale_color_manual(values = c(default_color, highlight_color)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Standard Deviation of Buffering Ratio", y = "Mean Buffering Ratio",
       color = "Frequently Buffered", alpha = "Frequently Buffered") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm')) +
  guides(colour = guide_legend(title.position="top", title.hjust = 0))

ora_buf <- low_var_buf %>%
  filter(Top50) %>%
  filter(Dataset == "CPTAC") %>%
  arrange(Gene.BR.SD.MeanNormRank) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis(ordered = TRUE)

ora_buf_terms <- ora_buf %>%
  plot_terms_compact(selected_sources = c("GO:MF", "KEGG", "WP"), custom_color = highlight_color)

# === Combine Panels into Figure ===
#dc_illus <- cowplot::draw_image(here(illustrations_dir, "dc_illustration.svg"))
buf_classes <- plot_buffering_ratio_classes()
illustrations <- cowplot::plot_grid(NULL, buf_classes, labels = c("A", "B"), rel_widths = c(0.63, 0.37))

stacked_buf_all <- cowplot::plot_grid(stacked_buf_cn, stacked_buf_cn_chr + ylab(NULL),
                                      ncol = 2, align = "h", axis = "tb",
                                      labels = c("E", "F"), rel_widths = c(1, 0.7))

figure1_sub1 <- cowplot::plot_grid(dc_class_line, stacked_buf_all,
                                   rel_widths = c(0.5, 1), labels = c("C", ""), ncol = 2)

figure1_sub2 <- cowplot::plot_grid(prot_exp_dc_class, dc_dataset_dist,
                                   rel_widths = c(1, 1), labels = c("D", "G"), ncol = 2)

br_by_cnv_all <- cowplot::plot_grid(br_by_cnv, br_by_cna + ylab(NULL),
                                    rel_widths = c(1, 0.9), labels = c("H", "I"), ncol = 2,
                                    align = "h", axis = "tb")

figure1_sub3 <- cowplot::plot_grid(br_by_cnv_all, scatter_signif_buffered, ora_buf_terms,
                                   rel_widths = c(1, 0.8, 0.75), labels = c("", "J", "K"), ncol = 3)

figure1 <- cowplot::plot_grid(figure1_sub1, figure1_sub2, figure1_sub3,
                              nrow = 3, rel_heights = c(0.8, 1, 1))

figure1_illustr <- cowplot::plot_grid(illustrations, figure1, nrow = 2, rel_heights = c(0.25, 1))

cairo_pdf(here(plots_dir, "figure01_01.pdf"), width = 13, height = 16)
figure1_illustr
dev.off()


# === Tables ===

## Buffering Classes & Ratios
t1_field_descriptions <- c(
  "Dataset" = "Proteomics dataset used for analysis (e.g., DepMap, ProCan, CPTAC, etc.).",
  "Model.ID" = "Unique identifier for cell lines and tumor samples.",
  "CellLine.DepMapModelId" = "Unique identifier for cell lines; Provided by DepMap CCLE.",
  "CellLine.SangerModelId" = "Unique identifier for cell lines; Provided by Sanger Cell Model Passports.",
  "CellLine.Name" = "Name of the cell line.",
  "CellLine.Replicate" = "Biological replicate of the cell line; Relevant for engineered cell lines.",
  "Model.Type" = "Cell Line or Tumor Sample.",
  "Gene.Symbol" = "HGNC gene symbol; updated using HGNChelper.",
  "Gene.Chromosome" = "Chromosome the gene is encoded on.",
  "Gene.ChromosomeArm" = "Chromosome arm the gene is encoded on.",
  "Gene.ChromosomeBand" = "Band of the chromosome the gene is encoded on.",
  "Gene.CopyNumber" = "Absolute copy number of a gene; Generated using ABSOLUTE for DepMap and ProCan, and using AscatNGS for CPTAC. See manuscript for dataset references.",
  "Gene.CopyNumber.Baseline" = "Gene copy number assumed as baseline for calculating the buffering ratio. For tumor samples: 2; For cell lines: Median of Gene.CopyNumber for samples where ChromosomeArm.CNA = 0, for each Gene.Symbol and Dataset.",
  "ChromosomeArm.CopyNumber" = "Copy number of the chromosome; Determined as ChromosomeArm.CopyNumber.Baseline + ChromosomeArm.CNA.",
  "ChromosomeArm.CopyNumber.Baseline" = "Chromosome copy number assumed as baseline for calculating the buffering ratio. Determined as the rounded average ploidy of a cell line.",
  "ChromosomeArm.CNA" = "Chromosome copy number alterations relative to background ploidy (loss: -1, neutral: 0, gain: +1); Obtained from Cohen-Sharir et al. (2021) through DepMap.",
  "Protein.Uniprot.Accession" = "Uniprot accession number of the protein.",
  "Protein.Expression.Normalized" = "LOESS-normalized log2 protein abundance. See manuscript for dataset references.",
  "Protein.Expression.Baseline.Unweighted" = "Mean neutral log2 protein abundance. Calculated as the mean of Protein.Expression.Normalized for samples where ChromosomeArm.CNA = 0, for each Gene.Symbol and Dataset.",
  "Protein.Expression.Baseline" = "Median neutral log2 protein abundance. Calculated as the median of Protein.Expression.Normalized for samples where ChromosomeArm.CNA = 0, for each Gene.Symbol and Dataset.",
  "Protein.Expression.Average" = "Mean log2 protein abundance upon chromosome arm gain or loss. Calculated as the mean of Protein.Expression.Normalized for samples where ChromosomeArm.CNA is +1 or -1, for each Gene.Symbol and Dataset.",
  "Log2FC" = "Log2 fold-change between observed protein abundance and baseline protein abundance. Calculated as Protein.Expression.Normalized - Protein.Expression.Baseline.",
  "Log2FC.Average" = "Log2 fold-change between mean protein abundance upon chromosome arm gain or loss and mean neutral protein abundance. Calculated as Protein.Expression.Average - Protein.Expression.Baseline.Unweighted.",
  "Buffering.GeneLevel.Ratio" = "Buffering Ratio calculated using gene copy numbers (GeneCN analysis variant). Calculated using Protein.Expression.Normalized, Protein.Expression.Baseline, Gene.CopyNumber, and Gene.CopyNumber.Baseline.",
  "Buffering.GeneLevel.Ratio.Confidence" = "Confidence score of the GeneCN-based buffering ratio. Higher values indicate higher confidence, i.e., less potential influence of noise in the buffering ratio.",
  "Buffering.GeneLevel.Class" = "Buffering class (Scaling, Buffered, Anti-Scaling) for the GeneCN analysis variant; Derived from Buffering.GeneLevel.Ratio and Log2FC.",
  "Buffering.ChrArmLevel.Ratio" = "Buffering Ratio calculated using chromosome arm copy numbers (ChrArm analysis variant). Calculated using Protein.Expression.Normalized, Protein.Expression.Baseline, ChromosomeArm.CopyNumber, and ChromosomeArm.CopyNumber.Baseline.",
  "Buffering.ChrArmLevel.Ratio.Confidence" = "Confidence score of the ChrArm-based buffering ratio. Higher values indicate higher confidence, i.e., less potential influence of noise in the buffering ratio.",
  "Buffering.ChrArmLevel.Class" = "Buffering class (Scaling, Buffered, Anti-Scaling) for the ChrArm analysis variant; Derived from Buffering.ChrArmLevel.Ratio and Log2FC.",
  "Buffering.ChrArmLevel.Average.Class" = "Buffering class (Scaling, Buffered, Anti-Scaling) for the ChrArm (avg.) analysis variant; Derived from ChromosomeArm.CNA and Log2FC.Average using Log2FC-based thresholds.",
  "Buffering.ChrArmLevel.Log2FC.Class" = "Buffering class (Scaling, Buffered, Anti-Scaling) for the ChrArm (Log2FC) analysis variant; Derived from ChromosomeArm.CNA and Log2FC using Log2FC-based thresholds."
)

expr_buf_p0211_publish <- expr_buf_p0211 %>%
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome),
         Dataset = "P0211",
         Model.Type = "Cell Line (Engineered)"
  )

expr_buf_chunduri_publish <- expr_buf_chunduri %>%
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome),
         Dataset = "Chunduri et al., 2021",
         Model.Type = "Cell Line (Engineered)"
  )

sup_table1 <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  mutate(Gene.Chromosome = as.integer(Gene.Chromosome)) %>%
  bind_rows(expr_buf_p0211_publish, expr_buf_chunduri_publish) %>%
  mutate(CellLine.Replicate = as.integer(CellLine.Replicate)) %>%
  select(all_of(names(t1_field_descriptions))) %>%
  write_parquet(here(tables_dir, "supplementary_table1.parquet"), version = "2.4")

### Export README for Table S1
t1_description <- c(
  "This table contains buffering ratios (BR) and buffering classes for all analysis varaints described in the methods section of the manuscript (GeneCN, ChrArm, ChrArm (avg.), ChrArm (Log2FC)).",
  "Protein abundance and copy number baselines required for calculating the BR are also included.",
  "The values of this table were determined for the proteomics datasets of DepMap CCLE, Sanger ProCan, CPTAC, Chunduri et al. (2021), and for P0211 (chromosome-engineered RPE-1 cell lines)."
)

createParquetReadme(t1_description, t1_field_descriptions, title = "Supplementary Table 1 - README",
                    readme_path = here(tables_dir, "README_supplementary_table1.md"),
                    file_path = "supplementary_table1.parquet", parquet_version = "2.4")

## Frequently Buffered Genes
t2_field_descriptions <- c(
  "=== TABLES ===" = "",
  "Frequently Buffered Genes" = "This table contains information about whether a gene was frequently and consistently buffered across multiple datasets, independent of whether the gene was affected by copy number gain or loss.",
  "Freq. Buffered Genes (CN Gain)" = "This table contains information about whether a gene was frequently and consistently buffered upon gene copy number gain across multiple datasets.",
  "Freq. Buffered Genes (CN Loss)" = "This table contains information about whether a gene was frequently and consistently buffered upon gene copy number loss across multiple datasets.",
  "=== COLUMNS ===" = "",
  "Dataset" = "Proteomics dataset used for analysis (e.g., DepMap, ProCan, CPTAC, etc.).",
  "Gene.Symbol" = "HGNC gene symbol; updated using HGNChelper.",
  "Gene.BR.Mean" = "Mean of gene copy number-derived buffering ratios across samples within a dataset.",
  "Gene.BR.SD" = "Standard deviation of gene copy number-derived buffering ratios across samples within a dataset.",
  "Samples" = "Number of samples in a dataset with valid gene copy number-derived buffering ratios.",
  "Gene.Buffered.Count" = "Number of samples within a dataset where the gene was classified as Buffered using gene copy number-derived buffering ratios.",
  "Gene.Buffered.Share" = "Fraction of samples within a dataset where the gene was classified as Buffered using gene copy number-derived buffering ratios.",
  "Gene.Buffered.Share.Median" = "Median fraction of Buffered genes within a dataset; Calculated as median of Gene.Buffered.Share across all genes within a dataset.",
  "Gene.BR.Test" = "Test used to determine whether the distribution of buffering ratios for this gene within is significantly higher than the threshold in Gene.BR.Test.Mu.",
  "Gene.BR.Test.Mu" = "Threshold on the buffering ratio used for Gene.BR.Test.",
  "Gene.BR.Test.p" = "p-value resulting from the test specified in Gene.BR.Test.",
  "Gene.BR.Test.p.adjusted" = "Benjamini-Hochberg adjusted p-value resulting from the test specified in Gene.BR.Test.",
  "Gene.BR.SD.MeanNormRank" = "Mean normalized rank (MNR) of buffering ratio standard deviations of a gene; Calculated for each gene as the MNR of Gene.BR.SD across datasets.",
  "Top50" = "Top 50 frequently and consistently buffered genes (ranked by lowest Gene.BR.SD.MeanNormRank). Genes were considered frequently buffered if Gene.Buffered.Share was at least 1/3 and above Gene.Buffered.Share.Median. Genes were considered consitently buffered, if Gene.BR.SD < 2 and Gene.BR.Test.p.adjusted < 0.05.",
  "Top10" = "Top 10 frequently and consistently buffered genes (ranked by lowest Gene.BR.SD.MeanNormRank)."
)

df_t2_fields <- data.frame(Column = names(t2_field_descriptions), Description = unname(t2_field_descriptions))

wb <- createWorkbook()
sheet_readme <- addWorksheet(wb, "README")
sheet_low_buf <- addWorksheet(wb, "Frequently Buffered Genes")
sheet_low_buf_gain <- addWorksheet(wb, "Freq. Buffered Genes (CN Gain)")
sheet_low_buf_loss <- addWorksheet(wb, "Freq. Buffered Genes (CN Loss)")
writeDataTable(wb = wb, sheet = sheet_readme, x = df_t2_fields)
writeDataTable(wb = wb, sheet = sheet_low_buf, x = low_var_buf)
writeDataTable(wb = wb, sheet = sheet_low_buf_gain, x = low_var_buf_gain)
writeDataTable(wb = wb, sheet = sheet_low_buf_loss, x = low_var_buf_loss)
saveWorkbook(wb, here(tables_dir, "supplementary_table2.xlsx"), overwrite = TRUE)

# === Supplemental Figures ===
## Copy Number Shares
expr_buf_all <- bind_rows(expr_buf_procan, expr_buf_cptac) %>%
  mutate(WGD = "All")

gene_cn_share <- expr_buf_procan %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD")) %>%
  bind_rows(expr_buf_all) %>%
  filter(Gene.CopyNumber.Baseline == 2) %>%
  mutate(GeneCN = as.factor(case_when(Gene.CopyNumber == 1 ~ 1L,
                                      Gene.CopyNumber == 2 ~ 2L,
                                      Gene.CopyNumber == 3 ~ 3L,
                                      Gene.CopyNumber >= 4 ~ 4L))) %>%
  select(GeneCN, Dataset, WGD) %>%
  drop_na() %>%
  count(GeneCN, Dataset, WGD) %>%
  group_by(Dataset, WGD) %>%
  mutate(Share = n / sum(n)) %>%
  ungroup()

panel_gene_cn_share <- gene_cn_share %>%
  mutate(Dataset = factor(Dataset, levels = dataset_order)) %>%
  ggplot() +
  aes(fill = GeneCN, y = Share, x = WGD) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_palettes$CopyNumbers, labels = c("1", "2", "3", "\u22654")) +
  labs(x = NULL, y = "Fraction", fill = "Gene\nCopy Number") +
  guides(fill = guide_legend(nrow = 2,byrow = TRUE)) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = c("left", "top"),
        legend.title = element_text(hjust = 1),
        legend.margin = margin(0,0,0,-1, unit = 'cm'))

### Chr Arm
chr_share <- expr_buf_procan %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD")) %>%
  bind_rows(expr_buf_all) %>%
  mutate(CNA = case_when(ChromosomeArm.CNA > 0 ~ "Gain",
                         ChromosomeArm.CNA == 0 ~ "Neutral",
                         ChromosomeArm.CNA < 0 ~ "Loss"),
         CNA = factor(CNA, levels = c("Gain", "Neutral", "Loss"))) %>%
  distinct(Model.ID, Gene.ChromosomeArm, CNA, WGD) %>%
  drop_na() %>%
  count(CNA, WGD) %>%
  mutate(Share = n / sum(n), .by = WGD)

panel_chr_share <- chr_share %>%
  ggplot() +
  aes(fill = CNA, y = Share, x = WGD) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(Gain = color_palettes$CopyNumbers[["4"]],
                               Neutral = "darkgrey",
                               Loss = color_palettes$CopyNumbers[["1"]])) +
  labs(x = NULL, y = "Fraction", fill = "Chromosome\nArm") +
  guides(fill = guide_legend(nrow = 2,byrow = TRUE)) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = c("left", "top"),
        legend.title = element_text(hjust = 1),
        legend.margin = margin(0,0,0,-1, unit = 'cm'))

## Tumor Purity
mean_br_cptac <- expr_buf_cptac %>%
  #filter(Buffering.GeneLevel.Ratio.Confidence > 0.3) %>%
  summarize(Mean.BR = mean(Buffering.GeneLevel.Ratio, na.rm = TRUE),
            Observations = sum(!is.na(Buffering.GeneLevel.Ratio)),
            .by = Model.ID) %>%
  #filter(Observations > 50) %>%
  left_join(y = df_model_cptac, by = "Model.ID")

cor_purity <- cor.test(mean_br_cptac$Mean.BR,
                       mean_br_cptac$Model.TumorPurity, method = "spearman")

panel_purity <- mean_br_cptac %>%
  ggplot() +
  aes(y = Mean.BR, x = Model.TumorPurity) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = lm, color = highlight_color) +
  annotate("text", x = min(mean_br_cptac$Model.TumorPurity), y = 1.5, hjust = 0, size = 5,
           label = paste0(print_corr(cor_purity$estimate), ", ", print_signif(cor_purity$p.value))) +
  labs(y = "Mean Buffering Ratio", x = "Tumor Purity")

## Immune Scores
cor_immune <- cor.test(mean_br_cptac$Mean.BR,
                       mean_br_cptac$ESTIMATE_ImmuneScore, method = "spearman")
cor_micro <- cor.test(mean_br_cptac$Mean.BR,
                      mean_br_cptac$xCell_microenvironment_score, method = "spearman")
cor_stroma <- cor.test(mean_br_cptac$Mean.BR,
                       mean_br_cptac$xCell_stroma_score, method = "spearman")

panel_immune <- mean_br_cptac %>%
  ggplot() +
  aes(y = Mean.BR, x = ESTIMATE_ImmuneScore) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = lm, color = highlight_color) +
  annotate("text", x = min(mean_br_cptac$ESTIMATE_ImmuneScore), y = 1.5, hjust = 0, size = 5,
           label = paste0(print_corr(cor_immune$estimate), ", ", print_signif(cor_immune$p.value))) +
  labs(y = "Mean Buffering Ratio", x = "Immune Score (ESTIMATE)")

panel_micro <- mean_br_cptac %>%
  ggplot() +
  aes(y = Mean.BR, x = xCell_microenvironment_score) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = lm, color = highlight_color) +
  annotate("text", x = min(mean_br_cptac$xCell_microenvironment_score), y = 1.5, hjust = 0, size = 5,
           label = paste0(print_corr(cor_micro$estimate), ", ", print_signif(cor_micro$p.value))) +
  labs(y = "Mean Buffering Ratio", x = "Microenvironment Score (xCell)")

panel_stroma <- mean_br_cptac %>%
  ggplot() +
  aes(y = Mean.BR, x = xCell_stroma_score) +
  geom_point(size = 2, color = default_color) +
  stat_smooth(method = lm, color = highlight_color) +
  annotate("text", x = min(mean_br_cptac$xCell_stroma_score), y = 1.5, hjust = 0, size = 5,
           label = paste0(print_corr(cor_stroma$estimate), ", ", print_signif(cor_stroma$p.value))) +
  labs(y = "Mean Buffering Ratio", x = "Stroma Score (xCell)")

## Buffering Ratio, WGD-Controlled
panel_br_wgd <- bind_rows(expr_buf_depmap, expr_buf_procan) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss"),
         WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD")) %>%
  select(Buffering.GeneLevel.Ratio, Dataset, CNV, WGD) %>%
  drop_na() %>%
  ggplot() +
  aes(y = Buffering.GeneLevel.Ratio, x = CNV) +
  geom_boxplot(outliers = FALSE) +
  stat_summary(aes(y = -1), fun.data = \(x) show.n(x, prefix = "n="), geom = "text", color = default_color) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = 1,
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  facet_grid(vars(WGD), vars(Dataset)) +
  labs(x = "Gene Copy Number", y = "Buffering Ratio")

## Buffering Classes, WGD-Controlled
df_share_chr_wgd <- bind_rows(expr_buf_depmap, expr_buf_procan) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Chr Arm Gain", "Chr Arm Loss"),
         WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD")) %>%
  select(Buffering.ChrArmLevel.Class, Dataset, CNV, WGD) %>%
  drop_na() %>%
  count(Buffering.ChrArmLevel.Class, Dataset, CNV, WGD) %>%
  group_by(Dataset, CNV, WGD) %>%
  mutate(Share = n / sum(n)) %>%
  ungroup()

panel_buf_chr_wgd <- df_share_chr_wgd %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_grid(vars(WGD), vars(CNV)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  labs(fill = NULL, x = NULL, y = "Fraction") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit = 'cm'))

## Buffering Classes (all)
# TODO: Buffering Classes (all)

## Frequently Buffered Genes (divided by Gain & Loss)
# TODO: Frequently Buffered (divided by Gain & Loss)

## Inter-Dataset Correlation
buf_matched <- match_datasets(expr_buf_procan, expr_buf_depmap)
corr_gene <- dataset_correlation(buf_matched,
                                 Dataset, Buffering.GeneLevel.Ratio,
                                 "GeneCN")
corr_chr <- dataset_correlation(buf_matched,
                                Dataset, Buffering.ChrArmLevel.Ratio,
                                "ChrArm")
corr_chr_avg <- dataset_correlation(buf_matched,
                                    Dataset, Buffering.ChrArmLevel.Average.Ratio,
                                    "ChrArm (avg.)")

panel_corr <- bind_rows(corr_gene, corr_chr) %>%
  jittered_boxplot(Comparison, Correlation, alpha = 1, jitter_width = 0.25, size = 2) +
  labs(y = "BR Correlation\n(DepMap ~ ProCan)", x = "Analysis Variant")

## Inter-Dataset Correlation (frequently buffered genes)
panel_corr_top50 <- read_parquet(here(output_data_dir, "br_correlation_gene.parquet")) %>%
  signif_boxplot(Top50, cor) +
  lims(y = c(0, 1)) +
  labs(x = "Top50 Frequently Buffered Genes", y = "BR Correlation\n(DepMap ~ ProCan)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

## Combine Figures
figure_s1_sub1 <- cowplot::plot_grid(panel_gene_cn_share, panel_br_wgd, panel_buf_chr_wgd,
                                     labels = c("A", "B", "C"), ncol = 3, rel_widths = c(0.9, 1, 1))
figure_s1_sub2 <- cowplot::plot_grid(panel_purity, panel_micro + ylab(NULL) + theme(axis.text.y = element_blank()),
                                     panel_immune + ylab(NULL) + theme(axis.text.y = element_blank()),
                                     panel_stroma + ylab(NULL) + theme(axis.text.y = element_blank()),
                                     rel_widths = c(1, 0.9, 0.85, 0.85), nrow = 1)
figure_s1_sub3 <- cowplot::plot_grid(panel_corr, panel_corr_top50,
                                     rel_widths = c(1, 0.6), labels = c("E", "F"))

figure_s1 <- cowplot::plot_grid(figure_s1_sub1, figure_s1_sub2, figure_s1_sub3,
                                nrow = 3, rel_heights = c(1, 0.75, 0.9), labels = c("", "D", ""))

cairo_pdf(here(plots_dir, "figure_s1.pdf"), width = 12, height = 11)
figure_s1
dev.off()
