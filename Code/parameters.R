library(ggplot2)
library(RColorBrewer)

# === Define Thresholds ===
log2fc_threshold <- 1     # Log2FC threshold, Default: 1 (= 2-fold change)
p_threshold <- 0.05       # p-value threshold, Default: 0.05

noisefloor_percentile_threshold <- 0.0001   # Percentile threshold to determine value of noise floor, Default: 0.0001 (0.01%)


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
biderectional_color_pal <- rev(brewer.pal(5, "RdBu"))
categorical_color_pal <- brewer.pal(12, "Paired")