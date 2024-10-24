library(ggplot2)
library(here)
library(RColorBrewer)
library(viridisLite)

 # === Set UTF-8 Encoding ===
options(encoding = "utf-8")
Sys.setlocale(category="LC_ALL", locale = "English_United States.1252")

## Special UTF-8 characters
utf8_rho <- enc2utf8("\u03C1")    # ρ
utf8_tau <- enc2utf8("\u03C4")    # τ
utf8_delta <- enc2utf8("\u0394")  # Δ

# === Define Thresholds ===
log2fc_threshold <- 0.5     # Log2FC threshold, Default: 1 (= 2-fold change)
p_threshold <- 0.05         # p-value threshold, Default: 0.05

noisefloor_percentile_threshold <- 0.0001   # Percentile threshold to determine value of noise floor, Default: 0.0001 (0.01%)

# === Bootstrap Settings ===
bootstrap_n <- 10000
bootstrap_sample_prop <- 1.0

# === Define Directories ===
input_data_dir <- here("Data")
external_data_dir <- here(input_data_dir, "External")
output_data_base_dir <- here("Output", "Data")
plots_base_dir <- here("Output", "Plots")
tables_base_dir <- here("Output", "Tables")
reports_base_dir <- here("Output", "Reports")
temp_base_dir <- here("Output", "Temp")
downloads_base_dir <- here("Downloads")
illustrations_dir <- here("Illustrations")

# === Theme Defaults ===
default_color <- "#2B2B2B"
highlight_colors <- c("#03A678", "#F27405", "#43BF71")
highlight_color <- highlight_colors[1]
unidirectional_color_pal <- brewer.pal(5, "Greens")
bidirectional_color_pal <- rev(brewer.pal(5, "RdBu"))
bidirectional_color_pal2 <- rev(brewer.pal(5, "PiYG"))
bidirectional_color_pal_viridis <- viridis(n = 5, option = "D", direction = 1)
categorical_color_pal <- brewer.pal(12, "Paired")
categorical_color_pal_dark <- brewer.pal(8, "Dark2")
two_class_color_pal <- c(categorical_color_pal[4], categorical_color_pal[8])
discrete_color_pal1 <- c("#E6B000", "#E6882C", "#2C7FE6", "#465566", "#6B644B")
discrete_color_pal1b <- c("#D9572B", "#D99100", "#2B74D9", "#3D4E66", "#66583D")
discrete_color_pal2 <- c("#33C653", "#C6A433", "#3349C6", "#C63345")
discrete_color_pal2_bright <- c("#47E669", "#E6C047", "#475EE6", "#E6475A")
discrete_color_pal2_dark <- c("#C8A51A", "#8F3635", "#2E3390", "#39724B")
discrete_color_pal3 <- c("#593659", "#803380", "#BF4DBF", "#CC7ACC")
discrete_color_pal4 <- c("#D9AB21", "#99823D", "#59513A", "#33322E")

color_palettes <- list(
  BufferingRatio = "viridis",
  AneuploidyScore = "rocket",
  Missing = "#808080",
  WGD = c("WGD" = categorical_color_pal[4],
          "Non-WGD" = categorical_color_pal[8]),
  CopyNumbers = c("1" = bidirectional_color_pal[1],
                  "3" = bidirectional_color_pal[4],
                  "4" = bidirectional_color_pal[5]),
  Datasets = c("DepMap" = discrete_color_pal3[1],
               "ProCan" = discrete_color_pal3[2],
               "Engin." = discrete_color_pal3[3],
               "CPTAC" = discrete_color_pal4[1]),
  BufferingClasses = c("Scaling" = viridis(2, begin = 0.15, end = 0.7)[1],
                       "Buffered" = viridis(2, begin = 0.15, end = 0.7)[2],
                       "Anti-Scaling" = "dimgrey"),
  DiffExp = c("Up" = bidirectional_color_pal[5],
              "Down" = bidirectional_color_pal[1]),
  DiffExpBackground = c("Up" = rev(brewer.pal(9, "RdBu"))[8],
                        "Down" = rev(brewer.pal(9, "RdBu"))[2])
)

base_size <- 15
default_theme <- theme_light(base_size = base_size) +
  theme(strip.background = element_rect(fill = default_color))
theme_set(default_theme)

# === Dosage Compensation Factors ===
dc_factor_cols <- c(
  "Protein-Protein Interactions", "Protein Half-Life", "Protein Complexes (CORUM)",
  "Mean 3'-UTR Length", "Mean 5'-UTR Length",
  "Phosphorylation Sites", "Ubiquitination Sites", "Sumoylation Sites",
  "Methylation Sites", "Acetylation Sites", "Regulatory Sites", "Kinase Interactions",
  "mRNA Abundance", "Protein Abundance", "Transcription Rate",
  "Translation Rate", "Protein Length", "mRNA Length",
  "Intrinsic Protein Disorder", "Low Complexity Score", "Homology Score",
  "Loops In Protein Score", "Protein Polyampholyte Score", "Protein Polarity",
  "Non-Exponential Decay Delta", "Mean mRNA Decay Rate", "Aggregation Score",
  "Haploinsufficiency", "Triplosensitivity", "Mean Gene Dependency", "Random Allelic Expression",
  "Transcription Factors (Repressor)", "Transcription Factors (Activator)",
  "Transcription Factors", "Mean TF Regulation Mode"
)

## Dataset-Specific Dosage Compensation Factors
dc_factor_cols_specific <- c("Protein Neutral CV")

# === Misc ===
dataset_order <- c("DepMap", "ProCan", "Engin.", "CPTAC")
