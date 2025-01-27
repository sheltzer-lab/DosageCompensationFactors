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
source(here("Code", "evaluation.R"))

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
analyze_sensitivity <- function(df_model_buf, df_drug_screens, buffering_col, drug_effect_col,
                                id_col = Drug.ID, group_cols = c("Drug.ID", "Drug.Name"),
                                p_thresh = p_threshold, log2fc_thresh = log2fc_threshold) {
  df_merged <- df_model_buf %>%
    inner_join(y = df_drug_screens, by = "Model.ID",
               relationship = "one-to-many", na_matches = "never")

  df_corr <- df_merged %>%
    summarize(Drug.Effect.Median = median({ { drug_effect_col } }, na.rm = TRUE),
              .by = c("Model.ID", group_cols, quo_name(enquo(buffering_col)))) %>%
    summarize(DrugEffect.Buffering.Corr = cor.test({ { buffering_col } }, Drug.Effect.Median,
                                                   method = "spearman")$estimate[["rho"]],
              DrugEffect.Buffering.Corr.p = cor.test({ { buffering_col } }, Drug.Effect.Median,
                                                     method = "spearman")$p.value,
              .by = group_cols) %>%
    mutate(DrugEffect.Buffering.Corr.p.adj = p.adjust(DrugEffect.Buffering.Corr.p, method = "BH"))

  df_diff <- df_merged %>%
    split_by_quantiles({ { buffering_col } }, target_group_col = "Buffering.Group", quantile_low = "20%", quantile_high = "80%") %>%
    differential_expression({ { id_col } }, Buffering.Group, { { drug_effect_col } },
                            groups = c("Low", "High"), test = t.test,
                            log2fc_thresh = log2fc_thresh, p_thresh = p_thresh)

  df_results <- df_corr %>%
    inner_join(y = df_diff, by = quo_name(enquo(id_col))) %>%
    select(-GroupA, -GroupB) %>%
    rename(
      DrugEffect.Buffering.Group.Low.Count = Count_GroupA,
      DrugEffect.Buffering.Group.High.Count = Count_GroupB,
      DrugEffect.Buffering.Group.Low.Mean = Mean_GroupA,
      DrugEffect.Buffering.Group.High.Mean = Mean_GroupB,
      DrugEffect.Buffering.Group.Log2FC = Log2FC,
      DrugEffect.Buffering.Group.Log2FC.p = Test.p,
      DrugEffect.Buffering.Group.Log2FC.p.adj = Test.p.adj,
      DrugEffect.Buffering.Group.Log2FC.Significant = Significant
    ) %>%
    mutate(
      EffectiveIn = case_when(
        DrugEffect.Buffering.Group.Log2FC.Significant == "Down" &
          DrugEffect.Buffering.Group.High.Mean < -log2fc_thresh ~ "High Buffering",
        DrugEffect.Buffering.Group.Log2FC.Significant == "Up" &
          DrugEffect.Buffering.Group.Low.Mean < -log2fc_thresh ~ "Low Buffering",
        TRUE ~ NA
      ),
      CommonEffect = DrugEffect.Buffering.Corr.p.adj < p_thresh &
        DrugEffect.Buffering.Group.Log2FC.p.adj < p_thresh &
        sign(DrugEffect.Buffering.Corr) == sign(DrugEffect.Buffering.Group.Log2FC)
    )

  return(df_results)
}

drug_dc_corr_procan <- model_buf_procan %>%
  analyze_sensitivity(drug_screens, Model.Buffering.Ratio.ZScore, Drug.MFI.Log2FC, log2fc_thresh = 0.2)
drug_dc_corr_depmap <- model_buf_depmap %>%
  analyze_sensitivity(drug_screens, Model.Buffering.Ratio.ZScore, Drug.MFI.Log2FC, log2fc_thresh = 0.2)
drug_dc_corr_agg_rank <- model_buf_agg %>%
  analyze_sensitivity(drug_screens, Model.Buffering.MeanNormRank, Drug.MFI.Log2FC, log2fc_thresh = 0.2)
drug_dc_corr_agg_mean <- model_buf_agg %>%
  analyze_sensitivity(drug_screens, Model.Buffering.StandardizedMean, Drug.MFI.Log2FC, log2fc_thresh = 0.2)


# === Evaluate Reproducibility of Results ===
cor_datasets <- cor.test(drug_dc_corr_depmap$DrugEffect.Buffering.Corr,
         drug_dc_corr_procan$DrugEffect.Buffering.Corr,
         method = "pearson")

cor_agg_methods <- cor.test(drug_dc_corr_agg_rank$DrugEffect.Buffering.Corr,
         drug_dc_corr_agg_mean$DrugEffect.Buffering.Corr,
         method = "pearson")

cat(capture.output(cor_datasets), file = here(reports_dir, "drug_correlation_datasets.txt"),
    append = FALSE, sep = "\n")
cat(capture.output(cor_agg_methods), file = here(reports_dir, "drug_correlation_aggregation-methods.txt"),
    append = FALSE, sep = "\n")


drug_dc_corr_agg_rank %>% mutate(Method = "StandardizedMean") %>%
  bind_rows(drug_dc_corr_agg_mean %>% mutate(Method = "MeanNormRank")) %>%
  ggplot(aes(x = DrugEffect.Buffering.Corr, color = Method)) +
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

# === Violin Plot ===
drug_dc_corr_agg_rank %>%
  mutate(Label = if_else(!is.na(DrugEffect.Buffering.Group.Log2FC.Significant), Drug.Name, NA)) %>%
  plot_volcano(DrugEffect.Buffering.Group.Log2FC, DrugEffect.Buffering.Group.Log2FC.p.adj, Label,
               DrugEffect.Buffering.Group.Log2FC.Significant, value_threshold = 0.2) %>%
  save_plot("buffering_drug_sensitivity_volcano.png", width = 300, height = 250)

# === Determine Drug Candidates using Buffering-DrugEffect Correlation ===
## Drug candidates for high-buffering cells
top_sensitivity <- drug_dc_corr_agg_rank %>%
  #filter(EffectiveIn == "High Buffering") %>%
  filter(DrugEffect.Buffering.Group.Log2FC.p < p_threshold & DrugEffect.Buffering.Group.Log2FC < -0.2 &
           DrugEffect.Buffering.Group.High.Mean < -0.2) %>%
  slice_min(DrugEffect.Buffering.Corr, n = 10)
## Drug candidates for low-buffering cells
bot_sensitivity <- drug_dc_corr_agg_rank %>%
  #filter(EffectiveIn == "Low Buffering") %>%
  filter(DrugEffect.Buffering.Group.Log2FC.p < p_threshold & DrugEffect.Buffering.Group.Log2FC > 0.2 &
           DrugEffect.Buffering.Group.Low.Mean < -0.2) %>%
  slice_max(DrugEffect.Buffering.Corr, n = 10)

df_sensitivity_agg <- model_buf_agg %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  mutate(Model.Buffering = if_else(Model.Buffering.MeanNormRank > 0.5, "High", "Low"))

## Visualize drug effect of candidates in all cell lines
df_sensitivity_agg %>%
  inner_join(y = top_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(DrugEffect.Buffering.Corr, 3),
                                                                 nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, desc(DrugEffect.Buffering.Corr))) %>%
  arrange(Model.Buffering.MeanNormRank) %>%
  jittered_boxplot(Label, Drug.MFI.Log2FC, Model.Buffering.MeanNormRank,
                   alpha = 3/4, jitter_width = 0.25) %>%
  save_plot("correlation_buffering_sensitivity_top.png", width = 300)

df_sensitivity_agg %>%
  inner_join(y = bot_sensitivity, by = c("Drug.ID", "Drug.Name")) %>%
  mutate(Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(DrugEffect.Buffering.Corr, 3),
                                                                 nsmall = 3, scientific = FALSE), ")")) %>%
  mutate(Label = fct_reorder(Label, DrugEffect.Buffering.Corr)) %>%
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
    mutate(Drug.Label = paste0(Drug.Name, " (", utf8_rho, " = ", format(round(DrugEffect.Buffering.Corr, 3),
                                                                        nsmall = 3, scientific = FALSE), ")")) %>%
    mutate(Drug.Label = fct_reorder(Drug.Label, DrugEffect.Buffering.Corr, .desc = desc),
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

# === Analyze Difference & Correlation of Drug Effect by Drug Mechanism ===
drug_screens_moa <- drug_screens %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID", relationship = "many-to-many", na_matches = "never") %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  drop_na(Drug.MOA) %>%
  mutate(Drug.MOA = str_squish(Drug.MOA))

results_moa <- model_buf_agg %>%
  analyze_sensitivity(drug_screens_moa, Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                      id_col = Drug.MOA, group_cols = "Drug.MOA", log2fc_thresh = 0.2) %>%
  write_parquet(here(output_data_dir, "drug_effect_buffering_mechanism.parquet"))

results_moa %>%
  filter(DrugEffect.Buffering.Corr.p.adj < p_threshold) %>%
  slice_max(abs(DrugEffect.Buffering.Corr), n = 20) %>%
  mutate(Drug.MOA = fct_reorder(Drug.MOA, DrugEffect.Buffering.Corr, .desc = TRUE),
         `-log10(p)` = -log10(DrugEffect.Buffering.Corr.p.adj)) %>%
  vertical_bar_chart(Drug.MOA, DrugEffect.Buffering.Corr, `-log10(p)`,
                     value_range = c(-0.2, 0.2), line_intercept = 0, bar_label_shift = 0.1,
                     category_lab = "Mechanism of Action", value_lab = "Correlation Buffering-DrugEffect",
                     color_lab = "-log10(p)") %>%
  save_plot("mechanism_buffering_sensitivity_top.png", width = 200)

results_moa %>%
  mutate(Label = if_else(!is.na(DrugEffect.Buffering.Group.Log2FC.Significant), str_trunc(Drug.MOA, 20), NA)) %>%
  plot_volcano(DrugEffect.Buffering.Group.Log2FC, DrugEffect.Buffering.Group.Log2FC.p.adj, Label,
               DrugEffect.Buffering.Group.Log2FC.Significant, value_threshold = 0.2) %>%
  save_plot("mechanism_buffering_sensitivity_volcano.png", width = 300, height = 250)


# === Analyze Difference & Correlation of Drug Effect by Drug Target ===
drug_screens_target <- drug_screens %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID", relationship = "many-to-many", na_matches = "never") %>%
  separate_longer_delim(Drug.Target, delim = ", ") %>%
  drop_na(Drug.Target) %>%
  mutate(Drug.Target = str_squish(Drug.Target))

results_target <- model_buf_agg %>%
  analyze_sensitivity(drug_screens_target, Model.Buffering.MeanNormRank, Drug.MFI.Log2FC,
                      id_col = Drug.Target, group_cols = "Drug.Target", log2fc_thresh = 0.2) %>%
  write_parquet(here(output_data_dir, "drug_effect_buffering_target.parquet"))

results_target %>%
  filter(DrugEffect.Buffering.Corr.p.adj < p_threshold) %>%
  slice_max(abs(DrugEffect.Buffering.Corr), n = 20) %>%
  mutate(Drug.Target = fct_reorder(Drug.Target, DrugEffect.Buffering.Corr, .desc = TRUE),
         `-log10(p)` = -log10(DrugEffect.Buffering.Corr.p.adj)) %>%
  vertical_bar_chart(Drug.Target, DrugEffect.Buffering.Corr, `-log10(p)`,
                     value_range = c(-0.2, 0.2), line_intercept = 0, bar_label_shift = 0.1,
                     category_lab = "Drug Target", value_lab = "Correlation Buffering-DrugEffect",
                     color_lab = "-log10(p)") %>%
  save_plot("target_buffering_sensitivity_top.png", width = 200)

results_target %>%
  mutate(Label = if_else(!is.na(DrugEffect.Buffering.Group.Log2FC.Significant), str_trunc(Drug.Target, 20), NA)) %>%
  plot_volcano(DrugEffect.Buffering.Group.Log2FC, DrugEffect.Buffering.Group.Log2FC.p.adj, Label,
               DrugEffect.Buffering.Group.Log2FC.Significant, value_threshold = 0.2) %>%
  save_plot("target_buffering_sensitivity_volcano.png", width = 300, height = 250)

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
  left_join(y = copy_number, by = "Model.ID") %>%
  mutate(Buffering = factor(if_else(Model.Buffering.MeanNormRank > median(Model.Buffering.MeanNormRank, na.rm = TRUE),
                             "Buffering", "Other"), levels = c("Buffering", "Other")),
         Resistant = factor(if_else(Drug.Effect.Median > median(Drug.Effect.Median, na.rm = TRUE),
                             "Resistant", "Other"), levels = c("Resistant", "Other")),
         Aneuploidy = factor(if_else(CellLine.AneuploidyScore > median(CellLine.AneuploidyScore, na.rm = TRUE),
                                     "High Aneuploidy", "Low Aneuploidy"), levels = c("High Aneuploidy", "Low Aneuploidy"))) %>%
  count(Buffering, Resistant, Aneuploidy) %>%
  drop_na() %>%
  arrange(Buffering, Resistant, Aneuploidy)

buf_res_plot <- response_counts %>%
  group_by(Buffering) %>%
  mutate(Share = (n / sum(n))) %>%
  ungroup() %>%
  ggplot() +
  aes(x = Buffering, y = Share, fill = Resistant, label = n) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(Resistant = highlight_colors[2],
                               Other = color_palettes$Missing))

buf_res_plot %>%
  save_plot("drug_resistance_buffering_share.png")

buf_res_test <- response_counts %>%
  summarize(n = sum(n), .by = c("Buffering", "Resistant")) %>%
  pivot_wider(names_from = "Resistant", values_from = "n") %>%
  tibble::column_to_rownames("Buffering") %>%
  as.matrix() %>%
  fisher.test()

### Control for aneuploidy
buf_res_plot_aneuploidy <- response_counts %>%
  group_by(Buffering, Aneuploidy) %>%
  mutate(Share = (n / sum(n))) %>%
  ungroup() %>%
  ggplot() +
  aes(x = Buffering, y = Share, fill = Resistant, label = n) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~Aneuploidy) +
  scale_fill_manual(values = c(Resistant = highlight_colors[2],
                               Other = color_palettes$Missing))

buf_res_plot_aneuploidy %>%
  save_plot("drug_resistance_buffering_share_aneuploidy.png")

### Test stochastic independence, high aneuploidy
buf_res_test_high_as <- response_counts %>%
  filter(Aneuploidy == "High Aneuploidy") %>%
  summarize(n = sum(n), .by = c("Buffering", "Resistant")) %>%
  pivot_wider(names_from = "Resistant", values_from = "n") %>%
  tibble::column_to_rownames("Buffering") %>%
  as.matrix() %>%
  fisher.test()

### Test stochastic independence, low aneuploidy
buf_res_test_low_as <- response_counts %>%
  filter(Aneuploidy == "Low Aneuploidy") %>%
  summarize(n = sum(n), .by = c("Buffering", "Resistant")) %>%
  pivot_wider(names_from = "Resistant", values_from = "n") %>%
  tibble::column_to_rownames("Buffering") %>%
  as.matrix() %>%
  fisher.test()

list("No Control" = buf_res_test,
     "High Aneuploidy" = buf_res_test_high_as,
     "Low Aneuploidy" = buf_res_test_low_as) %>%
  write_report(here(reports_base_dir, "drug_resistance_buffering_test.txt"))

# === Control: Compare drug sensitivity and drug mechanisms between high & low aneuploidy groups
## Drug sensitivity by aneuploidy
high_as_thresh <- quantile(copy_number$CellLine.AneuploidyScore, probs = 0.8)[[1]]
low_as_thresh <- quantile(copy_number$CellLine.AneuploidyScore, probs = 0.2)[[1]]

median_response_aneuploidy <- median_response %>%
  left_join(y = copy_number, by = "Model.ID") %>%
  drop_na(CellLine.AneuploidyScore) %>%
  mutate(Aneuploidy = case_when(CellLine.AneuploidyScore > high_as_thresh ~ "High",
                                CellLine.AneuploidyScore < low_as_thresh ~ "Low",
                                TRUE ~ NA)) %>%
  mutate(DrugStatus = case_when(Drug.Effect.Median > resistant_thresh ~ "Resistant",
                                Drug.Effect.Median < responsive_thresh ~ "Responsive",
                                TRUE ~ NA))

median_response_plot_aneuploidy <- median_response_aneuploidy %>%
  mutate(Aneuploidy = replace_na(Aneuploidy, "Other")) %>%
  ggplot() +
  aes(x = Aneuploidy, y = Drug.Effect.Median) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("High", "Low"),
                                 c("High", "Other")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(0.05, 0.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  theme(legend.position = "none") +
  ylim(c(-0.4, 0.15)) +
  ylab("Median Drug Effect")

median_response_plot_aneuploidy %>%
  save_plot("median_drug_effect_aneuploidy.png")

## Drug sensitivity by buffering (filtered by aneuploidy)
median_response_buf_aneuploidy <- median_response_buf %>%
  inner_join(y = copy_number, by = "Model.ID") %>%
  drop_na(CellLine.AneuploidyScore, Buffering) %>%
  mutate(Aneuploidy = if_else(CellLine.AneuploidyScore >= median(CellLine.AneuploidyScore),
                              "High Aneuploidy", "Low Aneuploidy"))

median_response_plot_buf_aneuploidy <- median_response_buf_aneuploidy %>%
  ggplot() +
  aes(x = Buffering, y = Drug.Effect.Median, color = Buffering) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("High", "Low")),
              test = wilcox.test,
              map_signif_level = print_signif, y_position = c(0.05, 0.1),
              size = 0.8, tip_length = 0, extend_line = -0.05, color = "black") +
  scale_color_manual(values = c("High" = color_palettes$BufferingClasses[["Buffered"]],
                                "Low" = color_palettes$BufferingClasses[["Scaling"]])) +
  theme(legend.position = "none") +
  ylim(c(-0.4, 0.15)) +
  ylab("Median Drug Effect") +
  facet_grid(~Aneuploidy)

median_response_plot_buf_aneuploidy %>%
  save_plot("median_drug_effect_by-aneuploidy.png")

## Drug mechanism (Difference) by aneuploidy
moa_diff_aneuploidy <- copy_number %>%
  inner_join(y = drug_screens, by = "Model.ID",
             relationship = "one-to-many", na_matches = "never") %>%
  select(-"Drug.Name") %>%
  left_join(y = drug_meta, by = "Drug.ID") %>%
  separate_longer_delim(Drug.MOA, delim = ", ") %>%
  mutate(Drug.MOA = str_squish(Drug.MOA)) %>%
  split_by_quantiles(CellLine.AneuploidyScore, target_group_col = "Aneuploidy",
                     quantile_low = "20%", quantile_high = "80%") %>%
  differential_expression(Drug.MOA, Aneuploidy, Drug.MFI.Log2FC,
                          groups = c("Low", "High"), log2fc_thresh = 0.2)

moa_diff_aneuploidy %>%
  mutate(Label = if_else(!is.na(Significant), str_trunc(Drug.MOA, 20), NA)) %>%
  plot_volcano(Log2FC, Test.p.adj, Label, Significant, value_threshold = 0.2) %>%
  save_plot("aneuploidy_drug_mechanism_volcano.png", width = 300, height = 250)

### Find drug mechanisms unique in bufffering
aneuploidy_up <- moa_diff_aneuploidy %>% filter(Significant == "Up")

unique_buf_moa <- results_moa %>%
  filter(DrugEffect.Buffering.Group.Log2FC.Significant == "Up") %>%
  filter(!(Drug.MOA %in% aneuploidy_up$Drug.MOA))

common_moa_lowbuf <- results_moa %>%
  filter(CommonEffect & EffectiveIn == "Low Buffering")

unique_buf_moa_common <- common_moa_lowbuf %>%
  filter(!(Drug.MOA %in% aneuploidy_up$Drug.MOA))

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
