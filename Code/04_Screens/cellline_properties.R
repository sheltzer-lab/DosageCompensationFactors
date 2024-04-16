library(here)
library(tidyr)
library(dplyr)
library(arrow)
library(stringr)

here::i_am("DosageCompensationFactors.Rproj")

source(here("Code", "parameters.R"))
source(here("Code", "visualization.R"))
source(here("Code", "analysis.R"))

procan_cn_data_dir <- here(external_data_dir, "CopyNumber", "ProCan")
depmap_cn_data_dir <- here(external_data_dir, "CopyNumber", "DepMap")
screens_data_dir <- here(external_data_dir, "Screens")
plots_dir <- here(plots_base_dir, "Screens", "CellLineProperties")
output_data_dir <- output_data_base_dir
reports_dir <- reports_base_dir

dir.create(plots_dir, recursive = TRUE)
dir.create(output_data_dir, recursive = TRUE)

cellline_buf_procan <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_procan.parquet"))
cellline_buf_depmap <- read_parquet(here(output_data_dir, "cellline_buffering_gene_filtered_depmap.parquet"))

copy_number <- read_parquet(here(output_data_dir, "copy_number.parquet")) %>%
  distinct(CellLine.Name, CellLine.AneuploidyScore, CellLine.WGD, CellLine.Ploidy)

df_model_procan <- read_csv_arrow(here(procan_cn_data_dir, "model_list_20230801.csv")) %>%
  rename(CellLine.Name = "model_name") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

df_model_depmap <- read_csv_arrow(here(depmap_cn_data_dir, "Model.csv")) %>%
  rename(CellLine.Name = "CellLineName") %>%
  add_count(CellLine.Name) %>%
  filter(n == 1) %>%
  select(-n)

aneuploidy_quant <- quantile(copy_number$CellLine.AneuploidyScore, probs = c(0.25, 0.5, 0.75))

# ToDo: Use Cell Line ID instead of cell line name
df_procan <- df_model_procan %>%
  inner_join(y = cellline_buf_procan, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"),
         `Near-Tetraploid` = if_else(ploidy_wes > 3.5, TRUE, FALSE),
         Aneuploidy = if_else(CellLine.AneuploidyScore > aneuploidy_quant["50%"], "High", "Low"),
         `Aneuploidy Quantiles` = case_when(
           CellLine.AneuploidyScore > aneuploidy_quant["75%"] ~ "Very High",
           CellLine.AneuploidyScore > aneuploidy_quant["50%"] ~ "High",
           CellLine.AneuploidyScore > aneuploidy_quant["25%"] ~ "Low",
           TRUE ~ "Very Low",
         ))


df_depmap <- df_model_depmap %>%
  inner_join(y = cellline_buf_depmap, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  inner_join(y = copy_number, by = "CellLine.Name",
             relationship = "one-to-one", na_matches = "never") %>%
  mutate(WGD = if_else(CellLine.WGD > 0, "WGD", "Non-WGD"),
         `Near-Tetraploid` = if_else(CellLine.Ploidy > 3.5, TRUE, FALSE),
         Aneuploidy = if_else(CellLine.AneuploidyScore > aneuploidy_quant["50%"], "High", "Low"),
         `Aneuploidy Quantiles` = case_when(
           CellLine.AneuploidyScore > aneuploidy_quant["75%"] ~ "Very High",
           CellLine.AneuploidyScore > aneuploidy_quant["50%"] ~ "High",
           CellLine.AneuploidyScore > aneuploidy_quant["25%"] ~ "Low",
           TRUE ~ "Very Low",
         ))

plot_categorical_properties <- function (df, cols_categorical) {
  plots <- list()
  for (col in cols_categorical) {
    plots[[col]] <- df %>%
      sorted_violin_plot(col, Buffering.CellLine.Ratio)
  }
  return(plots)
}

data_density_procan <- df_procan %>%
  data_density()

data_density_depmap <- df_procan %>%
  data_density()

cancers_depmap <- df_depmap %>%
  distinct(OncotreeLineage, OncotreePrimaryDisease, OncotreeSubtype, OncotreeCode)

# Plot distributions for categorical variables
cols_procan <- c("tissue_status", "cancer_type", "msi_status",
                 "smoking_status", "gender", "ethnicity",
                 "WGD", "Near-Tetraploid", "Aneuploidy Quantiles")

cols_depmap <- c("PrimaryOrMetastasis", "OncotreeSubtype", "Sex",
                 "WGD", "Near-Tetraploid", "Aneuploidy Quantiles")

violoin_plots_procan <- df_procan %>%
  plot_categorical_properties(cols_procan)

violoin_plots_depmap <- df_depmap %>%
  plot_categorical_properties(cols_depmap)

## Plot cell line buffering per cancer type
df_procan %>%
  sorted_beeswarm_plot("cancer_type", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_procan.png", width = 300)
df_depmap %>%
  sorted_beeswarm_plot("OncotreeSubtype", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_depmap.png", width = 300)

### Plot aneuploidy score per cancer type
df_procan %>%
  sorted_beeswarm_plot("cancer_type", CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_aneuploidy_procan.png", width = 300)
df_depmap %>%
  sorted_beeswarm_plot("OncotreeSubtype", CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_aneuploidy_depmap.png", width = 300)

### Control for low aneuploidy score when plotting buffering per cancer type
df_procan %>%
  filter(Aneuploidy == "Low") %>%
  sorted_beeswarm_plot("cancer_type", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_low-aneuploidy_procan.png", width = 300)
df_depmap %>%
  filter(Aneuploidy == "Low") %>%
  sorted_beeswarm_plot("OncotreeSubtype", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) %>%
  save_plot("cellline_cancer-type_low-aneuploidy_depmap.png", width = 300)

# Statistical comparisons
## Difference between Leukemia and Non-Leukemia cancer cells
leuk_procan <- c("B-Cell Non-Hodgkin's Lymphoma", "B-Lymphoblastic Leukemia", "Acute Myeloid Leukemia")
leuk_depmap <- c("Diffuse Large B-Cell Lymphoma, NOS", "B-Lymphoblastic Leukemia/Lymphoma", "Acute Myeloid Leukemia")

df_procan %>%
  mutate(Leukemia_Lymphoma = cancer_type %in% leuk_procan) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_procan.png")
leuk_plot_depmap <- df_depmap %>%
  mutate(Leukemia_Lymphoma = OncotreeSubtype %in% leuk_depmap) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_depmap.png")

### Control for low aneuploidy score
df_leuk_procan <- df_procan %>%
  filter(cancer_type %in% leuk_procan)
max_aneuploidy_leuk <- round(quantile(df_leuk_procan$CellLine.AneuploidyScore, probs = 0.9))

df_procan %>%
  filter(CellLine.AneuploidyScore <= max_aneuploidy_leuk) %>%
  mutate(Leukemia_Lymphoma = cancer_type %in% leuk_procan) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_low-aneuploidy_procan.png")
leuk_plot_depmap_low <- df_depmap %>%
  filter(CellLine.AneuploidyScore <= max_aneuploidy_leuk) %>%
  mutate(Leukemia_Lymphoma = OncotreeSubtype %in% leuk_depmap) %>%
  signif_beeswarm_plot(Leukemia_Lymphoma, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 1,
                       test = wilcox.test) %>%
  save_plot("cellline_leukemia_low-aneuploidy_depmap.png")


## Tumor status
df_procan %>%
  signif_beeswarm_plot(tissue_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_tumor-status_procan.png")
df_depmap %>%
  filter(PrimaryOrMetastasis != "NA") %>%
  signif_beeswarm_plot(PrimaryOrMetastasis, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_tumor-status_depmap.png")

## Micro-satellite instability
msi_comparison <- df_procan %>%
  signif_beeswarm_plot(msi_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi_procan.png")
df_procan %>%
  signif_beeswarm_plot(msi_status, CellLine.AneuploidyScore,
                       test = wilcox.test, cex = 1) %>%
  save_plot("cellline_msi-aneuploidy_procan.png")

### Control MSI/MSS distribution for low low aneuploidy score
msi_procan <- df_procan %>% filter(msi_status == "MSI")
fivenum_msi_aneuploidy <- fivenum(msi_procan$CellLine.AneuploidyScore)
max_msi_aneuploidy <- max(msi_procan$CellLine.AneuploidyScore)

msi_comparison_low <- df_procan %>%
  filter(CellLine.AneuploidyScore <= max_msi_aneuploidy) %>%
  signif_beeswarm_plot(msi_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi_low-aneuploidy_procan.png")
df_procan %>%
  filter(CellLine.AneuploidyScore <= max_msi_aneuploidy) %>%
  signif_beeswarm_plot(msi_status, CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi-aneuploidy_low-aneuploidy_procan.png")


### Resample, so that Aneuploidy Score distributions between MSI and MSS are equal
df_split <- split(df_procan, df_procan$msi_status)
df_msi_equal <- df_split$MSI %>%
  equalize_distributions(df_split$MSS, CellLine.AneuploidyScore,
                         with_replacement = FALSE, num_buckets = 5)

df_msi_equal %>%
  signif_beeswarm_plot(msi_status, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi_equal-aneuploidy_procan.png")
df_msi_equal %>%
  signif_beeswarm_plot(msi_status, CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_msi-aneuploidy_equal-aneuploidy_procan.png")

## Sex/Gender
df_procan %>%
  filter(gender != "Unknown") %>%
  signif_beeswarm_plot(gender, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_gender_procan.png")
df_depmap %>%
  filter(Sex != "Unknown") %>%
  signif_beeswarm_plot(Sex, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_gender_depmap.png")

## Whole-genome doubling
wgd_comparison <- df_procan %>%
  signif_beeswarm_plot(WGD, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_wgd_procan.png")
df_depmap %>%
  signif_beeswarm_plot(WGD, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_wgd_depmap.png")
df_depmap %>%
  filter(CellLine.WGD == 1 | CellLine.WGD == 2) %>%
  mutate(CellLine.WGD = factor(CellLine.WGD, levels = c(1, 2))) %>%
  signif_beeswarm_plot(CellLine.WGD, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_wgd_levels_depmap.png")

### Check aneuploidy score
df_procan %>%
  signif_beeswarm_plot(WGD, CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_wgd-aneuploidy_procan.png")
df_depmap %>%
  signif_beeswarm_plot(WGD, CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_wgd-aneuploidy_depmap.png")

### Control for low aneuploidy score
max_aneuploidy_nowgd <- round(quantile((df_procan %>% filter(WGD == "Non-WGD"))$CellLine.AneuploidyScore,
                                       probs = 0.9))
wgd_comparison_low <- df_procan %>%
  filter(CellLine.AneuploidyScore <= max_aneuploidy_nowgd) %>%
  signif_beeswarm_plot(WGD, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_wgd_low-aneuploidy_procan.png")
df_depmap %>%
  filter(CellLine.AneuploidyScore <= max_aneuploidy_nowgd) %>%
  signif_beeswarm_plot(WGD, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = t.test) %>%
  save_plot("cellline_wgd_low-aneuploidy_depmap.png")

## Near-tetraploid cell lines
df_procan %>%
  signif_beeswarm_plot(`Near-Tetraploid`, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_tetraploidy_procan.png")
df_depmap %>%
  signif_beeswarm_plot(`Near-Tetraploid`, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_tetraploidy_depmap.png")

## High vs. Low Aneuploidy Score
df_procan %>%
  signif_violin_plot(Aneuploidy, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_aneuploidy-class_procan.png")
df_depmap %>%
  signif_violin_plot(Aneuploidy, Buffering.CellLine.Ratio,
                     test = wilcox.test) %>%
  save_plot("cellline_aneuploidy-class_depmap.png")

### Control for Whole Genome Doubling
df_procan %>%
  signif_violin_plot(Aneuploidy, Buffering.CellLine.Ratio, WGD,
                     test = wilcox.test) %>%
  save_plot("cellline_aneuploidy-class_by-wgd_procan.png")
aneuploidy_comparison_wgd <- df_depmap %>%
  signif_violin_plot(Aneuploidy, Buffering.CellLine.Ratio, WGD,
                     test = wilcox.test) %>%
  save_plot("cellline_aneuploidy-class_by-wgd_depmap.png")


## Mutational Burden
mutational_burden_comparison <- df_procan %>%
  mutate(`Mutational Burden` = if_else(mutational_burden > mean(df_procan$mutational_burden),
                                       "High", "Low")) %>%
  signif_beeswarm_plot(`Mutational Burden`, Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore,
                       test = wilcox.test) %>%
  save_plot("cellline_mutations_comparison_procan.png")

# Regression analysis
## Age
df_procan %>%
  scatter_plot_reg_corr(age_at_sampling, Buffering.CellLine.Ratio,
                        point_size = 2, cor_method = "pearson") %>%
  save_plot("cellline_age_procan.png")
df_depmap %>%
  scatter_plot_reg_corr(Age, Buffering.CellLine.Ratio,
                        point_size = 2, cor_method = "pearson") %>%
  save_plot("cellline_age_depmap.png")

## Ploidy
df_procan %>%
  scatter_plot_reg_corr(ploidy_wes, Buffering.CellLine.Ratio, point_size = 2) %>%
  save_plot("cellline_ploidy_procan.png")
df_depmap %>%
  scatter_plot_reg_corr(CellLine.Ploidy, Buffering.CellLine.Ratio, point_size = 2) %>%
  save_plot("cellline_ploidy_depmap.png")

## Mutational burden
mutational_burden_reg_plot <- df_procan %>%
  scatter_plot_reg_corr(mutational_burden, Buffering.CellLine.Ratio, point_size = 2) %>%
  save_plot("cellline_mutations_procan.png")

## Aneuploidy score
aneuploidy_reg_plot <- df_procan %>%
  scatter_plot_reg_corr(CellLine.AneuploidyScore, Buffering.CellLine.Ratio, color_col = WGD, point_size = 2) %>%
  save_plot("cellline_aneuploidy_procan.png")
df_depmap %>%
  scatter_plot_reg_corr(CellLine.AneuploidyScore, Buffering.CellLine.Ratio, color_col = WGD, point_size = 2) %>%
  save_plot("cellline_aneuploidy_depmap.png")

## Misc
df_procan %>%
  scatter_plot_reg_corr(age_at_sampling, ploidy_wes, point_size = 2)

### Ploidy measurements differ across datasets
df_procan %>%
  scatter_plot_reg_corr(ploidy_wes, CellLine.Ploidy, point_size = 2)

# === Combine Plots for publishing ===
as_legend_label <- str_wrap("Aneuploidy Score", 10)

cancer_type_comparison <- df_depmap %>%
  sorted_beeswarm_plot("OncotreeCode", Buffering.CellLine.Ratio,
                       color_col = CellLine.AneuploidyScore, cex = 0.5) +
  labs(x = "Oncotree Cancer Code", color = as_legend_label)

bottom <- theme(legend.position = "bottom")
legend <- theme(legend.key.size = unit(14, "points"),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10))

plot1 <- cowplot::plot_grid(aneuploidy_reg_plot + legend + labs(x = "Aneuploidy Score") + theme(legend.title = element_blank()),
                            aneuploidy_comparison_wgd,
                            wgd_comparison + legend + labs(color = as_legend_label, x = NULL),
                            wgd_comparison_low + legend + labs(color = as_legend_label, x = NULL),
                            nrow = 2, ncol = 2, labels = c("A", "B", "C", ""))

plot3 <- cowplot::plot_grid(msi_comparison + legend + labs(color = as_legend_label, x = NULL),
                            msi_comparison_low + legend + labs(x = NULL, color = as_legend_label),
                            mutational_burden_reg_plot + legend + labs(x = "Mutational Burden"),
                            mutational_burden_comparison + legend + labs(color = as_legend_label),
                            nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

cairo_pdf(here(plots_dir, "cellline_properties_publish1.pdf"), height = 6, width = 9)
plot1
dev.off()

cairo_pdf(here(plots_dir, "cellline_properties_publish2.pdf"), height = 5, width = 9)
cancer_type_comparison
dev.off()

cairo_pdf(here(plots_dir, "cellline_properties_publish3.pdf"), height = 6, width = 9)
plot3
dev.off()

## Poster
### Confounders
plot_conf <- cowplot::plot_grid(aneuploidy_reg_plot +
                                  legend +
                                  theme_light(base_size = 20) +
                                  labs(x = "Aneuploidy Score", y = "Mean Buffering Ratio") +
                                  scale_color_manual(values = two_class_color_pal) +
                                  theme(legend.position = c(.95, .95),
                                        legend.justification = c("right", "top"),
                                        legend.box.just = "right",
                                        legend.margin = margin(6, 6, 6, 6),
                                        legend.title = element_blank(),
                                        legend.background = element_blank()),
                                aneuploidy_comparison_wgd +
                                  labs(y = "Mean Buffering Ratio") +
                                  theme_light(base_size = 20) +
                                  theme(strip.background = element_rect(fill = "grey90"),
                                        strip.text = element_text(color = "black")),
                                nrow = 2, ncol = 1)

cairo_pdf(here(plots_dir, "cellline_properties_poster_confounders.pdf"), height = 9)
plot_conf
dev.off()

### Leukemia
min_aneuploidy <- min(df_depmap$CellLine.AneuploidyScore)
max_aneuploidy <- max(df_depmap$CellLine.AneuploidyScore)

leuk_poster <- leuk_plot_depmap +
  theme_light(base_size = 20) +
  theme(legend.position = "top") +
  labs(color = NULL, y = "Mean Buffering Ratio") +
  scale_colour_viridis_c(option = "D", direction = 1,
                         limits = c(min_aneuploidy, max_aneuploidy),
                         breaks = c(min_aneuploidy, max_aneuploidy))

leuk_poster_low <- leuk_plot_depmap_low +
  theme_light(base_size = 20) +
  theme(legend.position = "top") +
  labs(color = NULL, y = NULL) +
  scale_colour_viridis_c(option = "D", direction = 1,
                         limits = c(min_aneuploidy, max_aneuploidy),
                         breaks = c(min_aneuploidy, max_aneuploidy_leuk[["90%"]], max_aneuploidy),
                         labels = c(min_aneuploidy, max_aneuploidy_leuk[["90%"]], max_aneuploidy))

plot_leuk <- cowplot::plot_grid(leuk_poster,
                                leuk_poster_low,
                                nrow = 1, ncol = 2)

cairo_pdf(here(plots_dir, "cellline_properties_poster_leukemia.pdf"), width = 8)
plot_leuk
dev.off()
