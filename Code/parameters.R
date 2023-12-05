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
log2fc_threshold <- 1     # Log2FC threshold, Default: 1 (= 2-fold change)
p_threshold <- 0.05       # p-value threshold, Default: 0.05

noisefloor_percentile_threshold <- 0.0001   # Percentile threshold to determine value of noise floor, Default: 0.0001 (0.01%)

# === Bootstrap Settings ===
bootstrap_n <- 10000
bootstrap_sample_prop <- 0.8

# === Define Directories ===
input_data_dir <- here("Data")
external_data_dir <- here(input_data_dir, "External")
output_data_base_dir <- here("Output", "Data")
plots_base_dir <- here("Output", "Plots")
tables_base_dir <- here("Output", "Tables")
reports_base_dir <- here("Output", "Reports")

# === Theme Defaults ===
theme_set(theme_light())
default_color <- "darkgrey"
highlight_color <- "#66CCB4"
unidirectional_color_pal <- brewer.pal(5, "Greens")
bidirectional_color_pal <- rev(brewer.pal(5, "RdBu"))
bidirectional_color_pal_viridis <- viridis(n = 5, option = "D", direction = 1)
categorical_color_pal <- brewer.pal(12, "Paired")

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
  "Haploinsufficiency Score", "Mean Gene Dependency"
)

## Dataset-Specific Dosage Compensation Factors
dc_factor_cols_specific <- c("Protein Neutral CV")