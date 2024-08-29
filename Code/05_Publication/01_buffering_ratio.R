library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))

plots_dir <- here(plots_base_dir, "Publication")
output_data_dir <- output_data_base_dir

dir.create(plots_dir, recursive = TRUE)

# === Load Datasets ===
expr_buf_procan <- read_parquet(here(output_data_dir, "expression_buffering_procan.parquet"))
expr_buf_depmap <- read_parquet(here(output_data_dir, "expression_buffering_depmap.parquet"))
expr_buf_cptac <- read_parquet(here(output_data_dir, 'expression_buffering_cptac_pure.parquet'))
expr_buf_p0211 <- read_parquet(here(output_data_dir, "expression_buffering_p0211.parquet"))

# === DC Class Illustration Panel ===
dc_illustration <- expr_buf_depmap %>%
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
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_path(linewidth = 1) +
  scale_y_continuous(limits = c(-0.6, 2.5),
                     breaks = seq(-0.5, 2.5, 0.5),
                     labels = paste0(seq(-0.5, 2.5, 0.5) * 100, "%")) +
  scale_color_manual(values = color_palettes$CopyNumbers) +
  labs(x = "Buffering Class", y = "Mean Protein Expression", color = "Gene Copy Number")

# === DC Class Protein Expression Panel ===
prot_exp_dc_class <- expr_buf_depmap %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  drop_na(Buffering.GeneLevel.Class, Protein.Expression.Normalized) %>%
  mutate(Gene.CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss"),
         Buffering.GeneLevel.Class = factor(Buffering.GeneLevel.Class,
                                            levels = c("Scaling", "Buffered", "Anti-Scaling"))) %>%
  ggplot() +
  aes(x = Buffering.GeneLevel.Class, y = Protein.Expression.Normalized, color = Buffering.GeneLevel.Class) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("Anti-Scaling", "Buffered"),
                                 c("Buffered", "Scaling"),
                                 c("Anti-Scaling", "Scaling")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(1.8, 1.8, 2.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  facet_grid(~Gene.CNV)

# === Tumor vs Cell Lines Panel ===
dc_dataset_dist <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac_pure) %>%
  filter(abs(Gene.CopyNumber - Gene.CopyNumber.Baseline) > 0.2) %>%
  mutate(Gene.CopyNumber.Event = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  select(Dataset, Buffering.GeneLevel.Ratio, Gene.CopyNumber.Event) %>%
  drop_na() %>%
  ggplot() +
  aes(x = Dataset, y = Buffering.GeneLevel.Ratio, color = Dataset) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("ProCan", "DepMap"),
                                 c("CPTAC", "DepMap"),
                                 c("CPTAC", "ProCan")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(1.8, 1.8, 2.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  theme(legend.position = "") +
  scale_color_manual(values = color_palettes$Datasets) +
  facet_wrap(~Gene.CopyNumber.Event, ncol = 2)

# === Buffering Gain vs. Loss Panel ===
expr_buf_p0211 %>%
  filter(CellLine.Name != "RPE1") %>%
  filter(Gene.Chromosome == 13) %>%
  ggplot() +
  aes(x = CellLine.Name, y = Buffering.ChrArmLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("RM13", "Rtr13")),
              map_signif_level = print_signif, y_position = 0.8, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black")

expr_buf_cptac %>%
  filter(Gene.CopyNumber.Baseline == 2) %>%
  filter(Gene.CopyNumber != Gene.CopyNumber.Baseline) %>%
  mutate(CNV = if_else(Gene.CopyNumber > Gene.CopyNumber.Baseline, "Gain", "Loss")) %>%
  ggplot() +
  aes(x = CNV, y = Buffering.GeneLevel.Ratio) +
  geom_boxplot(outliers = FALSE, size = 0.8) +
  geom_signif(comparisons = list(c("Gain", "Loss")),
              map_signif_level = print_signif, y_position = 2, size = 1,
              tip_length = 0, extend_line = -0.05, color = "black")

# === Low Variance Buffering Panel ===
low_var_buf <- bind_rows(expr_buf_depmap, expr_buf_procan, expr_buf_cptac) %>%
  analyze_low_br_variance()

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
       color = "Frequently Buffered", alpha = "Frequently Buffered")

ora_buf <- low_var_buf %>%
  filter(Top50) %>%
  pull(Gene.Symbol) %>%
  overrepresentation_analysis()

ora_buf %>%
  plot_terms(selected_sources = c("CORUM", "GO:MF", "KEGG", "WP"))