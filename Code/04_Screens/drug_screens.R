library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(openxlsx)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "analysis.R"))
source(here("Code", "visualization.R"))
source(here("Code", "buffering_ratio.R"))

output_data_dir <- output_data_base_dir
tables_dir <- tables_base_dir
plots_dir <- here(plots_base_dir, "Screens", "DrugSensitivity")
reports_dir <- reports_base_dir

dir.create(output_data_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)
dir.create(reports_dir, recursive = TRUE)

# === Load Datasets ===
drug_screens <- read_parquet(here(output_data_dir, "drug_screens.parquet")) %>% select(-CellLine.Name)
drug_meta <- read_parquet(here(output_data_dir, "drug_metadata.parquet"))
model_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
model_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))
model_buf_agg <- read_parquet(here(output_data_dir, "cellline_buffering_aggregated.parquet"))
copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(Model.ID, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

# === Calculate Correlation between Drug Sensitivity and Cell Line Buffering ===
analyze_sensitivity <- function(df_model_buf, df_drug_screens, buffering_col, drug_effect_col) {
  df_model_buf %>%
    inner_join(y = df_drug_screens, by = "Model.ID",
               relationship = "one-to-many", na_matches = "never") %>%
    mutate(Buffering.Median = median({ { buffering_col } }, na.rm = TRUE),
           Buffering.Group = if_else({ { buffering_col } } > Buffering.Median, "High", "Low")) %>%
    group_by(Drug.ID, Buffering.Group) %>%
    mutate(Group.Median = median({{drug_effect_col}}, na.rm = TRUE)) %>%
    group_by(Drug.ID, Drug.Name) %>%
    # https://stackoverflow.com/questions/48041504/calculate-pairwise-correlation-in-r-using-dplyrmutate
    summarize(
      # Correlation analysis Buffering vs Sensitivity
      Corr.DrugEffect_Buffering = cor.test({ { buffering_col } }, { { drug_effect_col } },
                                            method = "spearman")$estimate[["rho"]],
      Corr.p = cor.test({ { buffering_col } }, { { drug_effect_col } },
                        method = "spearman")$p.value,
      # Statistical test sensitivity low vs. high buffering
      BufferingGroup.DrugEffect.High = Group.Median[Buffering.Group == "High"][1],
      BufferingGroup.DrugEffect.Low = Group.Median[Buffering.Group == "Low"][1],
      BufferingGroup.DrugEffect.Diff = BufferingGroup.DrugEffect.High - BufferingGroup.DrugEffect.Low,
      BufferingGroup.Diff.Test.p = wilcox.test({{drug_effect_col}}[Buffering.Group == "High"],
                                               {{drug_effect_col}}[Buffering.Group == "Low"])$p.value,
    ) %>%
    ungroup() %>%
    mutate(Corr.p.adj = p.adjust(Corr.p, method = "BH"),
           BufferingGroup.Diff.Test.p.adj = p.adjust(BufferingGroup.Diff.Test.p, method = "BH")) %>%
    arrange(Drug.Name)
}

drug_dc_corr_procan <- model_buf_procan %>%
  analyze_sensitivity(drug_screens, Model.Buffering.Ratio.ZScore, Drug.MFI.Log2FC)
drug_dc_corr_depmap <- model_buf_depmap %>%
  analyze_sensitivity(drug_screens, Model.Buffering.Ratio.ZScore, Drug.MFI.Log2FC)
drug_dc_corr_agg_rank <- model_buf_agg %>%
  analyze_sensitivity(drug_screens, Model.Buffering.MeanNormRank, Drug.MFI.Log2FC)
drug_dc_corr_agg_mean <- model_buf_agg %>%
  analyze_sensitivity(drug_screens, Model.Buffering.StandardizedMean, Drug.MFI.Log2FC)


# === Evaluate Reproducibility of Results ===
cor_datasets <- cor.test(drug_dc_corr_depmap$Corr.DrugEffect_Buffering,
         drug_dc_corr_procan$Corr.DrugEffect_Buffering,
         method = "pearson")

cor_agg_methods <- cor.test(drug_dc_corr_agg_rank$Corr.DrugEffect_Buffering,
         drug_dc_corr_agg_mean$Corr.DrugEffect_Buffering,
         method = "pearson")

cat(capture.output(cor_datasets), file = here(reports_dir, "drug_correlation_datasets.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cor_agg_methods), file = here(reports_dir, "drug_correlation_aggregation-methods.txt"),
    append = FALSE, sep = "\n")


drug_dc_corr_agg_rank %>% mutate(Method = "StandardizedMean") %>%
  bind_rows(drug_dc_corr_agg_mean %>% mutate(Method = "MeanNormRank")) %>%
  ggplot(aes(x = Corr.DrugEffect_Buffering, color = Method)) +
  geom_density()

# === Write results ===
drug_dc_corr_procan %>%
  write_parquet(here(output_data_dir, "sensitivity_correlation_procan_gene_filtered.parquet")) %>%
  write.xlsx(here(tables_base_dir, "sensitivity_correlation_procan_gene_filtered.xlsx"), colNames = TRUE)
drug_dc_corr_depmap %>%
  write_parquet(here(output_data_dir, "sensitivity_correlation_depmap_gene_filtered.parquet")) %>%
  write.xlsx(here(tables_base_dir, "sensitivity_correlation_depmap_gene_filtered.xlsx"), colNames = TRUE)
drug_dc_corr_agg_rank %>%
  write_parquet(here(output_data_dir, "sensitivity_correlation_mean-norm-rank.parquet")) %>%
  write.xlsx(here(tables_base_dir, "sensitivity_correlation_mean-norm-rank.xlsx"), colNames = TRUE)
drug_dc_corr_agg_mean %>%
  write_parquet(here(output_data_dir, "sensitivity_correlation_standardized-mean.parquet")) %>%
  write.xlsx(here(tables_base_dir, "sensitivity_correlation_standardized-mean.xlsx"), colNames = TRUE)

# === Determine Drug Candidates using Buffering-DrugEffect Correlation ===
## Drug candidates for high-buffering cells
top_sensitivity <- drug_dc_corr_agg_rank %>%
  filter(BufferingGroup.Diff.Test.p < p_threshold & Corr.p < p_threshold) %>%
  filter(BufferingGroup.DrugEffect.Diff < 0 & BufferingGroup.DrugEffect.High < 0) %>%
  slice_min(Corr.DrugEffect_Buffering, n = 10)
## Drug candidates for high-buffering cells
bot_sensitivity <- drug_dc_corr_agg_rank %>%
  filter(BufferingGroup.Diff.Test.p < p_threshold & Corr.p < p_threshold) %>%
  filter(BufferingGroup.DrugEffect.Diff > 0 & BufferingGroup.DrugEffect.Low < 0) %>%
  slice_max(Corr.DrugEffect_Buffering, n = 10)

df_sensitivity_agg <- model_buf_agg %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  mutate(Model.Buffering = if_else(Model.Buffering.MeanNormRank > 0.5, "High", "Low"))

## Visualize drug effect of candidates in all cell lines
df_sensitivity_agg %>%
  inner_join(y = top_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(Corr.DrugEffect_Buffering, 3),
                                                                 nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, desc(Corr.DrugEffect_Buffering))) %>%
  arrange(Model.Buffering.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Model.Buffering.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_top.png", width = 300)

df_sensitivity_agg %>%
  inner_join(y = bot_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(Corr.DrugEffect_Buffering, 3),
                                                                 nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, Corr.DrugEffect_Buffering)) %>%
  arrange(Model.Buffering.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Model.Buffering.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_bot.png", width = 300)

## Visualize candidate drug effect difference between high and low buffering cells
df_sensitivity_agg %>%
  semi_join(y = top_sensitivity, by = "Drug.ID") %>%
  signif_violin_plot(Model.Buffering, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("buffering_sensitivity_distributions_top.png"), width = 350)

df_sensitivity_agg %>%
  semi_join(y = bot_sensitivity, by = "Drug.ID") %>%
  signif_violin_plot(Model.Buffering, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("buffering_sensitivity_distributions_bot.png"), width = 350)

## Scatter- & Violin-plot visualization of selected drugs
selected_drugs <- c("REGORAFENIB", "MPS1-IN-5", "SECLIDEMSTAT", "ATIPRIMOD", "INARIGIVIR", "C-021", "G-749",
                    "NIMORAZOLE", "GELDANAMYCIN", "ZOTAROLIMUS", "LERCANIDIPINE", "BAFILOMYCIN-A1")

for (drug in selected_drugs) {
  df_sensitivity_agg %>%
    filter(Drug.Name == drug) %>%
    scatter_plot_reg_corr(Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                          point_size = 2, title_prefix = drug) %>%
    save_plot(paste0("buffering_sensitivity_", drug, "_scatterplot.png"))
}

df_sensitivity_agg %>%
  filter(Drug.Name %in% selected_drugs) %>%
  signif_violin_plot(Model.Buffering, Drug.MFI.Log2FC, Drug.Name) %>%
  save_plot(paste0("buffering_sensitivity_distributions_selected.png"), width = 350)

## Get mechanism of action and target of drug candidates
top_meta <- drug_meta %>%
  semi_join(y = top_sensitivity, by = "Drug.ID")

bot_meta <- drug_meta %>%
  semi_join(y = bot_sensitivity, by = "Drug.ID")

selected_meta <- drug_meta %>%
  filter(Drug.Name %in% selected_drugs)

## Heatmap plotting Drug Sensitivity, Cell Line Buffering and potential confounders
drug_confounder_heatmap <- function(df, desc = FALSE) {
  theme_settings <- theme(axis.line.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.key.size = unit(12, "points"),
                          legend.title = element_text(size = 10),
                          legend.text = element_text(size = 10))

  horizontal_legend <- theme(legend.direction = "horizontal",
                             legend.title = element_blank())

  df <- df %>%
    mutate(Drug.Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(Corr.DrugEffect_Buffering, 3),
                                                                        nsmall = 3, scientific = FALSE), ")")) %>%
    mutate(Drug.Label = fct_reorder(Drug.Label, Corr.DrugEffect_Buffering, .desc = desc),
           CellLine.Name = fct_reorder(CellLine.Name, Model.Buffering.MeanNormRank),
           CellLine.WGD = as.factor(CellLine.WGD))

  heatmap_buffering <- df %>%
    distinct(CellLine.Name, Model.Buffering.MeanNormRank) %>%
    ggplot() +
    aes(x = CellLine.Name, y = "Buffering Score", fill = Model.Buffering.MeanNormRank) +
    geom_raster() +
    scale_fill_viridis_c(option = "mako", labels = c(0, 0.5, 1), breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    cowplot::theme_minimal_hgrid() +
    theme_settings +
    horizontal_legend

  heatmap_drugs <- df %>%
    ggplot() +
    aes(x = CellLine.Name, y = Drug.Label, fill = Drug.MFI.Log2FC) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    cowplot::theme_minimal_hgrid() +
    theme_settings

  heatmap_wgd <- df %>%
    distinct(CellLine.Name, CellLine.WGD) %>%
    ggplot() +
    aes(x = CellLine.Name, y = "WGD Status", fill = CellLine.WGD) +
    geom_raster() +
    scale_fill_viridis_d(option = "viridis") +
    cowplot::theme_minimal_hgrid() +
    theme_settings +
    horizontal_legend

  heatmap_aneuploidy <- df %>%
    distinct(CellLine.Name, CellLine.AneuploidyScore) %>%
    ggplot() +
    aes(x = CellLine.Name, y = "Aneuploidy Score", fill = CellLine.AneuploidyScore) +
    geom_raster() +
    scale_fill_viridis_c(option = "rocket") +
    cowplot::theme_minimal_hgrid() +
    theme_settings +
    horizontal_legend

  cowplot::plot_grid(
    heatmap_buffering, NULL, heatmap_drugs, NULL, heatmap_wgd, NULL, heatmap_aneuploidy,
    labels = NULL, ncol = 1, align = "v", axis = "lr",
    rel_heights = c(2, -0.5, length(selected_drugs), -0.5, 2, -0.8, 2)
  )
}

model_buf_agg %>%
  inner_join(y = copy_number, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = drug_dc_corr_agg_rank, by = c("Drug.ID", "Drug.Name")) %>%
  filter(Drug.Name %in% selected_drugs) %>%
  drug_confounder_heatmap(desc = TRUE) %>%
  save_plot("drug_confounder_heatmap_selected.png")

model_buf_agg %>%
  inner_join(y = copy_number, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = top_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  drug_confounder_heatmap(desc = TRUE) %>%
  save_plot("drug_confounder_heatmap_top.png")

model_buf_agg %>%
  inner_join(y = copy_number, by = "Model.ID",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  inner_join(y = bot_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  drug_confounder_heatmap() %>%
  save_plot("drug_confounder_heatmap_bot.png")

# === Determine Drug Mechanisms & Targets using Buffering-DrugEffect Correlation ===
## Calculate & plot correlation by mechanism of action
moa_corr <- model_buf_agg %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  drop_na(Drug.MOA) %>%
  mutate(Drug.MOA = str_squish(Drug.MOA)) %>%
  group_by(Drug.MOA) %>%
  summarize(
    # ToDo: Avoid calculating correlation twice
    Corr.DrugEffect_Buffering = cor.test(Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                                          method = "spearman")$estimate[["rho"]],
    Corr.p = cor.test(Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                      method = "pearson")$p.value
  ) %>%
  mutate(Corr.p.adj = p.adjust(Corr.p, method = "BH"))

moa_corr %>%
  filter(Corr.p.adj < p_threshold) %>%
  slice_max(abs(Corr.DrugEffect_Buffering), n = 20) %>%
  mutate(Drug.MOA = fct_reorder(Drug.MOA, Corr.DrugEffect_Buffering, .desc = TRUE),
         `-log10(p)` = -log10(Corr.p.adj)) %>%
  vertical_bar_chart(Drug.MOA, Corr.DrugEffect_Buffering, `-log10(p)`,
                     value_range = c(-0.2, 0.2), line_intercept = 0, bar_label_shift = 0.1,
                     category_lab = "Mechanism of Action", value_lab = "Correlation Buffering-DrugEffect",
                     color_lab = "-log10(p)") %>%
  save_plot("mechanism_buffering_sensitivity_top.png", width = 200)

## Calculate & plot correlation by drug target
target_corr <- model_buf_agg %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.Target, delim = ", ") %>%
  mutate(Drug.Target = str_squish(Drug.Target)) %>%
  drop_na(Drug.Target) %>%
  group_by(Drug.Target) %>%
  summarize(
    # ToDo: Avoid calculating correlation twice
    Corr.DrugEffect_Buffering = cor.test(Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                                          method = "spearman")$estimate[["rho"]],
    Corr.p = cor.test(Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                      method = "pearson")$p.value
  ) %>%
  mutate(Corr.p.adj = p.adjust(Corr.p, method = "BH"))

target_corr %>%
  filter(Corr.p.adj < p_threshold) %>%
  slice_max(abs(Corr.DrugEffect_Buffering), n = 20) %>%
  mutate(Drug.Target = fct_reorder(Drug.Target, Corr.DrugEffect_Buffering, .desc = TRUE),
         `-log10(p)` = -log10(Corr.p.adj)) %>%
  vertical_bar_chart(Drug.Target, Corr.DrugEffect_Buffering, `-log10(p)`,
                     value_range = c(0, 0.15), line_intercept = 0, bar_label_shift = 0.005,
                     category_lab = "Drug Target", value_lab = "Correlation Buffering-DrugEffect",
                     color_lab = "-log10(p)") %>%
  save_plot("target_buffering_sensitivity_top.png")

# === Compare MOA & Target between groups of positive and negative drug effect difference ===
## Note: nothing significant when filtering for adjusted p-value

lower_diff <- drug_dc_corr_agg_rank %>%
  filter(BufferingGroup.DrugEffect.Diff < -0.15 & BufferingGroup.DrugEffect.High < 0) %>%
  filter(BufferingGroup.Diff.Test.p < p_threshold) %>%
  left_join(y = drug_meta %>% select(-Drug.Name), by = "Drug.ID")

higher_diff <- drug_dc_corr_agg_rank %>%
  filter(BufferingGroup.DrugEffect.Diff > 0.15 & BufferingGroup.DrugEffect.Low < 0) %>%
  filter(BufferingGroup.Diff.Test.p < p_threshold) %>%
  left_join(y = drug_meta, by = "Drug.ID")

lower_diff_moa <- lower_diff %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  mutate(Drug.MOA = trimws(Drug.MOA)) %>%
  group_by(Drug.MOA) %>%
  add_count() %>%
  summarize(BufferingGroup.Diff.MOA.Median = median(BufferingGroup.DrugEffect.Diff, na.rm = TRUE),
            MOA.Count = first(n))

higher_diff_moa <- higher_diff %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  mutate(Drug.MOA = trimws(Drug.MOA)) %>%
  group_by(Drug.MOA) %>%
  add_count() %>%
  summarize(BufferingGroup.Diff.MOA.Median = median(BufferingGroup.DrugEffect.Diff, na.rm = TRUE),
            MOA.Count = first(n))

lower_diff_target <- lower_diff %>%
  separate_longer_delim(Drug.Target, delim = ", ") %>%
  mutate(Drug.Target = trimws(Drug.Target)) %>%
  group_by(Drug.Target) %>%
  add_count() %>%
  summarize(BufferingGroup.Diff.Target.Median = median(BufferingGroup.DrugEffect.Diff, na.rm = TRUE),
            Target.Count = first(n))

higher_diff_target <- higher_diff %>%
  separate_longer_delim(Drug.Target, delim = ", ") %>%
  mutate(Drug.Target = trimws(Drug.Target)) %>%
  group_by(Drug.Target) %>%
  add_count() %>%
  summarize(BufferingGroup.Diff.Target.Median = median(BufferingGroup.DrugEffect.Diff, na.rm = TRUE),
            Target.Count = first(n))

# === Difference in Drug Effect between High and Low Buffering Groups ===
model_buf_sensitivity <- model_buf_agg %>%
  split_by_quantiles(Model.Buffering.MeanNormRank, target_group_col = "Model.Buffering.Group",
                     quantile_low = "20%", quantile_high = "80%") %>%
  inner_join(y = drug_screens, by = "Model.ID") %>%
  differential_expression(Drug.ID, Model.Buffering.Group, Drug.MFI.Log2FC,
                          groups = c("Low", "High"), test = wilcox.test, centr = mean, log2fc_thresh = 0.2) %>%
  inner_join(y = drug_meta, by = "Drug.ID")

model_buf_sensitivity %>%
  mutate(Label = if_else(!is.na(Significant), Drug.Name, NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, value_threshold = 0.2) %>%
  save_plot("buffering_drug_sensitivity_volcano.png", width = 300, height = 250)

## Mechanism of Action
moa_diff <- model_buf_agg %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  mutate(Drug.MOA = str_squish(Drug.MOA)) %>%
  split_by_quantiles(Model.Buffering.MeanNormRank, target_group_col = "CellLine.Buffering.Group",
                     quantile_low = "20%", quantile_high = "80%") %>%
  differential_expression(Drug.MOA, CellLine.Buffering.Group, Drug.MFI.Log2FC,
                          groups = c("Low", "High"), log2fc_thresh = 0.2)

moa_diff %>%
  mutate(Label = if_else(!is.na(Significant), str_trunc(Drug.MOA, 20), NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, value_threshold = 0.2) %>%
  save_plot("buffering_drug_mechanism_volcano.png", width = 300, height = 250)


## Drug Target
target_diff <- model_buf_agg %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.Target, delim = ", ") %>%
  mutate(Drug.Target = str_squish(Drug.Target)) %>%
  split_by_quantiles(Model.Buffering.MeanNormRank, target_group_col = "CellLine.Buffering.Group",
                     quantile_low = "20%", quantile_high = "80%") %>%
  differential_expression(Drug.Target, CellLine.Buffering.Group, Drug.MFI.Log2FC,
                          groups = c("Low", "High"), log2fc_thresh = 0.2)

target_diff %>%
  mutate(Label = if_else(!is.na(Significant), str_trunc(Drug.Target, 20), NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, value_threshold = 0.2) %>%
  save_plot("buffering_drug_target_volcano.png", width = 300, height = 250)

# === Common Mechanisms and Targets between correlation and fold-change methods ===
## Mechanisms
moa_corr_lowbuf <- moa_corr %>%
  filter(Corr.p.adj < p_threshold & Corr.DrugEffect_Buffering > 0)
moa_corr_highbuf <- moa_corr %>%
  filter(Corr.p.adj < p_threshold & Corr.DrugEffect_Buffering < 0)

common_moa_lowbuf <- moa_diff %>%
  filter(Significant == "Up" & Mean_GroupA < 0) %>%
  #inner_join(higher_diff_moa, by = "Drug.MOA") %>%
  inner_join(moa_corr_lowbuf, by = "Drug.MOA")

common_moa_highbuf <- moa_diff %>%
  filter(Significant == "Down" & Mean_GroupB < 0) %>%
  #inner_join(lower_diff_moa, by = "Drug.MOA") %>%
  inner_join(moa_corr_highbuf, by = "Drug.MOA")

## Targets
target_corr_lowbuf <- target_corr %>%
  filter(Corr.p.adj < p_threshold & Corr.DrugEffect_Buffering > 0)
target_corr_highbuf <- target_corr %>%
  filter(Corr.p.adj < p_threshold & Corr.DrugEffect_Buffering < 0)

common_target_lowbuf <- target_diff %>%
  filter(Significant == "Up" & Mean_GroupA < 0) %>%
  #inner_join(higher_diff_target, by = "Drug.Target") %>%
  inner_join(target_corr_lowbuf, by = "Drug.Target")

common_target_highbuf <- target_diff %>%
  filter(Significant == "Down" & Mean_GroupB < 0) %>%
  #inner_join(lower_diff_target, by = "Drug.Target") %>%
  inner_join(target_corr_highbuf, by = "Drug.Target")

## Save results
bind_rows(common_moa_lowbuf, common_moa_highbuf) %>%
  write_parquet(here(output_data_dir, "drug_mechanisms_common.parquet")) %>%
  write.xlsx(here(tables_base_dir, "drug_mechanisms_common.xlsx"), colNames = TRUE)

bind_rows(common_target_lowbuf, common_target_highbuf) %>%
  write_parquet(here(output_data_dir, "drug_targets_common.parquet")) %>%
  write.xlsx(here(tables_base_dir, "drug_targets_common.xlsx"), colNames = TRUE)

# === Difference in drug sensitivity between high-buffering, drug-responsive and drug-resistant cells ===
median_response <- drug_screens %>%
  #semi_join(model_buf_agg, by = "Model.ID") %>%
  group_by(Model.ID) %>%
  summarize(Drug.Effect.Median = median(Drug.MFI.Log2FC, na.rm = TRUE),
            Drug.Effect.Observations = sum(!is.na(Drug.MFI.Log2FC)),
            .groups = "drop") %>%
  filter(Drug.Effect.Observations > 1000)

high_buf_thresh <- quantile(model_buf_agg$Model.Buffering.MeanNormRank, probs = 0.8)[[1]]
low_buf_thresh <- quantile(model_buf_agg$Model.Buffering.MeanNormRank, probs = 0.2)[[1]]
resistant_thresh <- quantile(median_response$Drug.Effect.Median, probs = 0.8)[[1]]
responsive_thresh <- quantile(median_response$Drug.Effect.Median, probs = 0.2)[[1]]

median_response_buf <- median_response %>%
  left_join(model_buf_agg, by = "Model.ID") %>%
  mutate(Buffering = case_when(Model.Buffering.MeanNormRank > high_buf_thresh ~ "High",
                            Model.Buffering.MeanNormRank < low_buf_thresh ~ "Low",
                            TRUE ~ NA)) %>%
  mutate(DrugStatus = case_when(Drug.Effect.Median > resistant_thresh ~ "Resistant",
                            Drug.Effect.Median < responsive_thresh ~ "Responsive",
                            TRUE ~ NA))

median_response_buf %>%
  scatter_plot_reg_corr(Model.Buffering.MeanNormRank, Drug.Effect.Median) %>%
  save_plot("median_drug_effect_correlation.png")

median_response_plot_buf <- median_response_buf %>%
  # drop_na(Model.Buffering.MeanNormRank) %>%
  mutate(Buffering = replace_na(Buffering, "Other")) %>%
  ggplot() +
  aes(x = Buffering, y = Drug.Effect.Median, color = Buffering) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("High", "Low"),
                                 c("High", "Other")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(0.05, 0.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c("High" = color_palettes$BufferingClasses[["Buffered"]],
                                "Low" = color_palettes$BufferingClasses[["Scaling"]],
                                "Other" = default_color)) +
  theme(legend.position = "none") +
  ylim(c(-0.4, 0.15)) +
  ylab("Median Drug Effect")

median_response_plot_control <- median_response_buf %>%
  drop_na(DrugStatus) %>%
  ggplot() +
  aes(x = DrugStatus, y = Drug.Effect.Median) +
  geom_boxplot() +
  theme(legend.position = "none") +
  xlab("Drug Control") +
  ylim(c(-0.4, 0.15))

cowplot::plot_grid(median_response_plot_buf,
                   median_response_plot_control
                     + ylab(NULL)
                     + theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank()),
                   ncol = 2, rel_widths = c(1, 2/3)) %>%
  save_plot("median_drug_effect.png")

## Check if drug resistance and high buffering are stochastically independent
response_counts <- median_response_buf %>%
  mutate(Buffering = factor(if_else(Model.Buffering.MeanNormRank > median(Model.Buffering.MeanNormRank, na.rm = TRUE),
                             "Buffering", "Other"), levels = c("Buffering", "Other")),
         Resistant = factor(if_else(Drug.Effect.Median > median(Drug.Effect.Median, na.rm = TRUE),
                             "Resistant", "Other"), levels = c("Resistant", "Other"))) %>%
  count(Buffering, Resistant) %>%
  drop_na() %>%
  arrange(Buffering, Resistant) %>%
  pivot_wider(names_from = "Resistant", values_from = "n") %>%
  tibble::column_to_rownames("Buffering") %>%
  as.matrix()

fisher.test(response_counts)

# === Drug Mechanism Groups ===
# drug_meta_long <- drug_meta %>%
#   separate_longer_delim(Drug.MOA, ", ") %>%
#   mutate(Drug.MOA = str_squish(Drug.MOA))
# mechanisms <- drug_meta_long %>%
#   count(Drug.MOA) %>%
#   drop_na()
# mechanism_dist <- outer(mechanisms$Drug.MOA, mechanisms$Drug.MOA,
#                         \(x, y) stringdist::stringdist(x, y, method = "lv"))
# dimnames(mechanism_dist) <- list(mechanisms$Drug.MOA, mechanisms$Drug.MOA)
# mechanism_clust <- hclust(as.dist(mechanism_dist))
# pheatmap::pheatmap(mechanism_dist)