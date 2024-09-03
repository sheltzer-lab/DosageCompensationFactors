library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(assertr)
library(ggplot2)
library(ggrepel)
library(skimr)
library(openxlsx)
library(gtsummary)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))
source(here("Code", "analysis.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Distribution")

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===

expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_cptac <- read_parquet(here(output_data_dir, "expression_buffering_cptac.parquet"))
expr_buf_p0211 <- read_parquet(here(output_data_dir, "expression_buffering_p0211.parquet")) %>%
  mutate(Gene.Chromosome = as.character(Gene.Chromosome))

# === Summarize Distribution of Obersavtions ===

summary_tbl_depmap <- expr_buf_depmap %>%
  tbl_summary(include = c(Gene.CopyNumber, Protein.Expression.Log2, Protein.Expression.Normalized, Protein.Expression.Average,
                          Protein.Expression.Baseline, Protein.Expression.Baseline.Unweighted,
                          Log2FC, Log2FC.Average,
                          Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio, Buffering.ChrArmLevel.Average.Ratio,
                          Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class, Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

summary_tbl_procan <- expr_buf_procan %>%
  tbl_summary(include = c(Gene.CopyNumber, Protein.Expression.Log2, Protein.Expression.Normalized, Protein.Expression.Average,
                          Protein.Expression.Baseline, Protein.Expression.Baseline.Unweighted,
                          Log2FC, Log2FC.Average,
                          Buffering.GeneLevel.Ratio, Buffering.ChrArmLevel.Ratio, Buffering.ChrArmLevel.Average.Ratio,
                          Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class, Buffering.ChrArmLevel.Average.Class),
              by = ChromosomeArm.CNA) %>%
  modify_header(label = "**Chromosome Arm CNA**") %>%
  bold_labels()

summary_tbl_merged <- tbl_merge(tbls = list(summary_tbl_depmap, summary_tbl_procan),
                                tab_spanner = c("**DepMap**", "**ProCan**"))
summary_tbl_merged %>%
  as_gt() %>%
  gt::gtsave(filename = here(tables_dir, "summary_table.html"))

# === Plot Copy Number Distribution and Filter Thresholds ===
# TODO: This has been previously defined for relative copy numbers, not absolute copy numbers
cn_diff_quantiles <- quantile(expr_buf_procan$Gene.CopyNumber - expr_buf_procan$Gene.CopyNumber.Baseline,
                              probs = seq(0, 1, 0.01))

cn_loss_thresh <- cn_diff_quantiles["10%"]
cn_gain_thresh <- cn_diff_quantiles["90%"]
breaks <- c(seq(-2, 6, 0.5), signif(cn_gain_thresh, digits = 2), signif(cn_loss_thresh, digits = 2))
labels <- c(seq(-2, 6, 0.5), signif(cn_gain_thresh, digits = 2), signif(cn_loss_thresh, digits = 2))
colors <- c(rep("black", length(seq(-2, 6, 0.5))), "red", "red")

cn_dist <- expr_buf_procan %>%
  mutate(Log2FC.CopyNumber = Gene.CopyNumber - Gene.CopyNumber.Baseline) %>%
  ggplot() +
  aes(x = Log2FC.CopyNumber) +
  geom_density() +
  geom_vline(xintercept = cn_loss_thresh, linetype = "dashed", color = "red") +
  geom_vline(xintercept = cn_gain_thresh, linetype = "dashed", color = "red") +
  annotate(x = cn_loss_thresh, y = -Inf, label = paste0("I\n", signif(cn_loss_thresh, digits = 2)), geom = "text",
           color = "red", lineheight = .8, vjust = 1) +
  annotate(x = cn_gain_thresh, y = -Inf, label = paste0("I\n", signif(cn_gain_thresh, digits = 2)), geom = "text",
           color = "red", lineheight = .8, vjust = 1) +
  scale_x_continuous(breaks = seq(-2, 6, 0.5)) +
  ylab("Density") +
  coord_cartesian(clip = "off")
save_plot(cn_dist, "copy_number_distribution_procan.png", height = 100)

# === Plot categorical distribution of buffering classes by analysis type ===
## Gene CN only (across datasets)
df_share_gene <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gene CN Gain", "Gene CN Loss")) %>%
  select(Buffering.GeneLevel.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.GeneLevel.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n)) * 100) %>%
  ungroup()

stacked_buf_cn <- df_share_gene %>%
  ggplot() +
  aes(fill = Buffering.GeneLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses)

stacked_buf_cn %>%
  save_plot("buffering_class_gene.png", width = 200)

## Chr Arm only (across datasets)
df_share_chr <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_p0211) %>%
  filter(ChromosomeArm.CopyNumber != ChromosomeArm.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(ChromosomeArm.CopyNumber > ChromosomeArm.CopyNumber.Baseline, "Chr Arm Gain", "Chr Arm Loss")) %>%
  select(Buffering.ChrArmLevel.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.ChrArmLevel.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n)) * 100) %>%
  ungroup()

stacked_buf_chr <- df_share_chr %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses)

stacked_buf_chr %>%
  save_plot("buffering_class_chromosome.png", width = 200)

## Chr Arm Average only (across datasets)
df_share_chr_avg <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_p0211) %>%
  distinct(Gene.Symbol, Dataset, Buffering.ChrArmLevel.Average.Class, ChromosomeArm.CNA) %>%
  filter(ChromosomeArm.CNA != 0) %>%
  mutate(CNV = if_else(ChromosomeArm.CNA > 0, "Chr Arm Gain", "Chr Arm Loss")) %>%
  select(Buffering.ChrArmLevel.Average.Class, Dataset, CNV) %>%
  drop_na() %>%
  count(Buffering.ChrArmLevel.Average.Class, Dataset, CNV) %>%
  group_by(Dataset, CNV) %>%
  mutate(Share = (n / sum(n)) * 100) %>%
  ungroup()

stacked_buf_chr_avg <- df_share_chr_avg %>%
  ggplot() +
  aes(fill = Buffering.ChrArmLevel.Average.Class, y = Share, x = Dataset) +
  geom_bar(stat = "identity") +
  facet_wrap(~CNV) +
  scale_fill_manual(values = color_palettes$BufferingClasses)

stacked_buf_chr_avg %>%
  save_plot("buffering_class_chromosome_average.png", width = 200)

## Compare all methods
df_share_all <- bind_rows(df_share_gene, df_share_chr, df_share_chr_avg) %>%
  pivot_longer(c(Buffering.GeneLevel.Class, Buffering.ChrArmLevel.Class, Buffering.ChrArmLevel.Average.Class),
               names_to = "AnalysisVariant", values_to = "BufferingClass") %>%
  drop_na()

dataset_order <- c("DepMap", "ProCan", "CPTAC", "P0211")

stacked_buf_class <- df_share_all %>%
  mutate(CNV = if_else(grepl("Gain", CNV), "Gain", "Loss"),
         Dataset = factor(Dataset, levels = dataset_order),
         AnalysisVariant = case_when(AnalysisVariant == "Buffering.GeneLevel.Class" ~ "Gene CN",
                                     AnalysisVariant == "Buffering.ChrArmLevel.Class" ~ "ChrArm",
                                     AnalysisVariant == "Buffering.ChrArmLevel.Average.Class" ~ "ChrArm (avg.)")) %>%
  ggplot() +
  aes(fill = BufferingClass, y = Share, x = AnalysisVariant) +
  geom_bar(position = "fill", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = color_palettes$BufferingClasses) +
  ggh4x::facet_nested(~Dataset + CNV)

stacked_buf_class %>%
  save_plot("buffering_class_all.png", width = 400)
