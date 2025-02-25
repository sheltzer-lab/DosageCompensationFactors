library(dplyr)
library(plotly)
library(ggplot2)

# Dosage Compensation Measures

buffering_ratio <- function(expr_base, expr_obs, cn_base = 2, cn_obs = 3) {
  br <- log2(cn_obs / cn_base) - log2(expr_obs / expr_base)
  br <- ifelse(cn_obs > cn_base, br, -br)
  br <- ifelse(cn_obs == cn_base, NA, br)
  ifelse(is.finite(br), br, NA)
}

# TODO: explain why (ratio - 1) is necessary
scaling_factor <- function(expr_base, expr_obs, cn_base = 2, cn_obs = 3) {
  sf <- ((expr_obs / expr_base) - 1) / ((cn_obs / cn_base) - 1)
  ifelse(is.finite(sf), sf, NA)
}

# Expr Log2FC and CN Log2FC are zero if there is no change in expression / copy number
scaling_ratio <- function(expr_base, expr_obs, cn_base = 2, cn_obs = 3, f = log2) {
  sr <- f(expr_obs / expr_base) / f(cn_obs / cn_base)
  ifelse(is.finite(sr) && sr >= 0, sr, NA)
}

# Only quantifies buffering against scaling - anti-scaling not quantified, as log2 of negative numbers is undefined
log_buffering_score <- function(expr_base, expr_obs, cn_base = 2, cn_obs = 3) {
  lbs <- log2(log2(cn_obs / cn_base) / log2(expr_obs / expr_base))
  ifelse(is.finite(lbs), lbs, NA)
}

# Define thresholds on dosage compensation measures for classification
## Scaling:       Change in protein expression is equal or more than the change in DNA copy number
## Buffered:      Change in protein expression is less than the change in DNA copy number
## Anti-Scaling:  Change in protein expression is opposite to the change in DNA copy number

br_cutoffs <- list(Buffered = 0.2, AntiScaling = 0.6)
buffering_class <- function(buffering_ratio, expr_base, expr_obs, cn_base, cn_obs, br_cutoffs_ = br_cutoffs) {
  scaling_direction <- sign(expr_obs - expr_base) * sign(cn_obs - cn_base)
  expr_abs_log2fc <- abs(log2(expr_obs / expr_base))

  ifelse(scaling_direction >= 0,
         ifelse(buffering_ratio > br_cutoffs_$Buffered, "Buffered",
                "Scaling"),
         ifelse(expr_abs_log2fc > 0.1, "Anti-Scaling",
                "Buffered")
  )
}

buffering_class_br_only <- function(buffering_ratio, expr_base, expr_obs, cn_base, cn_obs, br_cutoffs_ = br_cutoffs) {
  scaling_direction <- sign(expr_obs - expr_base) * sign(cn_obs - cn_base)

  ifelse(scaling_direction >= 0,
         ifelse(buffering_ratio > br_cutoffs_$Buffered, "Buffered",
                "Scaling"),
         ifelse(buffering_ratio > br_cutoffs_$AntiScaling, "Anti-Scaling",
                "Buffered")
  )
}

sf_cutoffs <- list(Scaling = 0.3784142, Buffered = -0.133934)
buffering_class_sf <- function(scaling_factor, sf_cutoffs_ = sf_cutoffs) {
  ifelse(scaling_factor > sf_cutoffs_$Scaling, "Scaling",
         ifelse(scaling_factor > sf_cutoffs_$Buffered, "Buffered",
                "Anti-Scaling"))
}

sr_cutoffs <- list(Scaling = 0.8, Buffered = -0.2)
buffering_class_sr <- function(scaling_ratio, sr_cutoffs_ = sr_cutoffs) {
  ifelse(scaling_ratio > sr_cutoffs_$Scaling, "Scaling",
         ifelse(scaling_ratio > sr_cutoffs_$Buffered, "Buffered",
                "Anti-Scaling"))
}

lbs_cutoffs <- list(Buffered = 0.3)
buffering_class_lbs <- function (log_buffering_score, lbs_cutoffs_ = lbs_cutoffs) {
  ifelse(log_buffering_score > lbs_cutoffs_$Buffered, "Buffered", "Scaling")
}

# Values from Schukken & Sheltzer, 2022 (DOI: 10.1101/gr.276378.121) for chromosome arm gain (inverted for arm loss):
# Scaling:        0.25 < Log2FC
# Buffered:      -0.1  < Log2FC <  0.25
# Anti-Scaling:          Log2FC < -0.1
buffering_class_log2fc <- function (log2fc, cn_base = 2, cn_var = 3) {
  ifelse(cn_var > cn_base,
         ifelse(log2fc > 0.25, "Scaling",             # Expression level as expected or higher
                ifelse(log2fc > -0.1, "Buffered",     # Less expression in trisomy than expected
                       "Anti-Scaling")),               # Buffering of trisomy expression below disomy level
         ifelse(log2fc < -0.25, "Scaling",
                ifelse(log2fc < 0.1, "Buffered",
                       "Anti-Scaling")))
}

# Confidence Scores
buffering_ratio_confidence <- function(cn_base, cn_obs, expr_neutral_cv) {
  conf <- abs(log2(cn_obs / cn_base)) * (1 / (1 + expr_neutral_cv))
  conf <- ifelse(cn_obs == cn_base, NA, conf)
  ifelse(is.finite(conf), conf, NA)
}

# Convert BR cutoffs to Log2FC cutoffs
br2logfc_cutoffs <- function(br_cutoff, cn_base, cn_obs) {
  ifelse(cn_obs == cn_base, NA,
         ifelse(cn_obs > cn_base,
                log2(cn_obs / cn_base) - br_cutoff,
                log2(cn_obs / cn_base) + br_cutoff)
  )
}

# === Example Code ===
plot_buffering_ratio <- function(br_func, cnv_lim = c(-1, 1), expr_lim = c(-1, 1)) {
  require(tibble)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(viridisLite)
  require(geomtextpath)

  cn_diff <- seq(cnv_lim[1], cnv_lim[2], by = 0.01)
  cn_base <- rep(2, length(cn_diff))
  expr_diff <- seq(expr_lim[1], expr_lim[2], by = 0.01)
  expr_base <- rep(2, length(expr_diff))

  br <- outer(expr_diff + expr_base, cn_diff + cn_base,
              FUN = \(x,y) br_func(expr_obs = x, cn_obs = y,
                                   expr_base = expr_base, cn_base = cn_base))
  dimnames(br) <- list(expr_diff + expr_base, cn_diff + cn_base)

  br_df <- br %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ExpressionLevel") %>%
    tidyr::pivot_longer(everything() & !ExpressionLevel, names_to = "CopyNumber", values_to = "BufferingRatio") %>%
    mutate_all(as.numeric)

  expr_ticks <- seq(expr_lim[1], expr_lim[2], 0.25) + expr_base[1]
  cn_ticks <- seq(cnv_lim[1], cnv_lim[2], 0.25) + cn_base[1]

  br_df %>%
    ggplot() +
    aes(x = CopyNumber, y = ExpressionLevel, z = BufferingRatio) +
    geom_raster(aes(fill = BufferingRatio)) +
    geom_vline(xintercept = cn_ticks, color = "darkgrey", alpha = 0.5) +
    geom_hline(yintercept = expr_ticks, color = "darkgrey", alpha = 0.5) +
    #geom_contour(color = "white") +
    geom_textcontour(color = "white") +
    scale_fill_viridis_c() +
    scale_x_continuous(limits = cnv_lim + cn_base[1], breaks = cn_ticks) +
    scale_y_continuous(limits = expr_lim + expr_base[1], breaks = expr_ticks) +
    theme_light()
}

plot_buffering_ratio_classes <- function(cnv_lim = c(-1.99, 2), buffered_threshold = br_cutoffs$Buffered) {
  require(tibble)
  require(dplyr)
  require(tidyr)
  require(ggplot2)

  cn_diff <- seq(cnv_lim[1], cnv_lim[2], by = 0.01)
  cn_base <- rep(2, length(cn_diff))

  df_logfc <- data.frame(CN_Base = cn_base, CN_Obs = cn_diff + cn_base) %>%
    mutate(CN_Log2FC = log2(CN_Obs) - log2(CN_Base),
           Expr_Log2FC = br2logfc_cutoffs(buffered_threshold, CN_Base, CN_Obs))

  buf_poly1 <- df_logfc %>%
    select(CN_Log2FC, Expr_Log2FC) %>%
    filter(CN_Log2FC %in% c(-1, max(df_logfc$CN_Log2FC[df_logfc$CN_Log2FC < 0]))) %>%
    bind_rows(data.frame(
      CN_Log2FC = c(max(df_logfc$CN_Log2FC[df_logfc$CN_Log2FC < 0]), -1),
      Expr_Log2FC = c(1, 1)
    )) %>%
    mutate(id = 1)

  buf_poly2 <- df_logfc %>%
    select(CN_Log2FC, Expr_Log2FC) %>%
    filter(CN_Log2FC %in% c(1, min(df_logfc$CN_Log2FC[df_logfc$CN_Log2FC > 0]))) %>%
    bind_rows(data.frame(
      CN_Log2FC = c(1, min(df_logfc$CN_Log2FC[df_logfc$CN_Log2FC > 0])),
      Expr_Log2FC = c(-1, -1)
    )) %>%
    mutate(id = 2)

  # Expr_Log2FC = 0 around abs(CN_Log2FC) = 0.3 for BR = 0.2
  undefined_cn_range <- c(-0.3, 0.3)

  bind_rows(buf_poly1, buf_poly2) %>%
    ggplot() +
    aes(x = CN_Log2FC, y = Expr_Log2FC, group = id) +
    geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = 1, fill = color_palettes$BufferingClasses[["Scaling"]]) +
    geom_polygon(fill = color_palettes$BufferingClasses[["Buffered"]]) +
    geom_rect(xmin = -1, xmax = 0, ymin = 0.1, ymax = 1, fill = color_palettes$BufferingClasses[["Anti-Scaling"]]) +
    geom_rect(xmin = 0, xmax = 1, ymin = -1, ymax = -0.1, fill = color_palettes$BufferingClasses[["Anti-Scaling"]]) +
    geom_rect(xmin = undefined_cn_range[1], xmax = undefined_cn_range[2], ymin = -1, ymax = 1, fill = "white") +
    geom_hline(yintercept = 0, color = default_color) +
    geom_vline(xintercept = 0, color = default_color) +
    #geom_line(color = color_palettes$BufferingClasses[["Buffered"]]) +
    annotate("text", x = -0.95, y = -0.9, label = "Scaling", hjust = 0, vjust = 0, color = "white", size = base_size/3) +
    annotate("text", x = -0.95, y = -0.1, label = "Buffered", hjust = 0, color = "white", size = base_size/3) +
    annotate("text", x = -0.95, y = 0.9, label = "Anti-Scaling", hjust = 0,  vjust = 1, color = "white", size = base_size/3) +
    annotate("text", x = 0.95, y = 0.9, label = "Scaling", hjust = 1, vjust = 1, color = "white", size = base_size/3) +
    annotate("text", x = 0.95, y = 0.1, label = "Buffered", hjust = 1, color = "white", size = base_size/3) +
    annotate("text", x = 0.95, y = -0.9, label = "Anti-Scaling", hjust = 1, vjust = 0, color = "white", size = base_size/3) +
    annotate("label", x = 0, y = 0, label = paste0("BR > ", buffered_threshold),
             hjust = 0.5, vjust = 0.5, color = color_palettes$BufferingClasses[["Buffered"]], size = base_size/3) +
    scale_x_continuous(limits = c(-1,1), breaks = c(seq(-1, 1, 0.5), undefined_cn_range), minor_breaks = NULL) +
    scale_y_continuous(limits = c(-1,1), breaks = seq(-1, 1, 0.5), minor_breaks = NULL) +
    labs(x = "Copy Number Log2FC", y = "Protein Log2FC") +
    theme_minimal(base_size = base_size)
}

plot_buffering_ratio_3d <- function(br_func, cnv_lim = c(-1, 1), expr_lim = c(-1, 1)) {
  cn_diff <- seq(cnv_lim[1], cnv_lim[2], by = 0.01)
  cn_base <- rep(2, length(cn_diff))
  expr_diff <- seq(expr_lim[1], expr_lim[2], by = 0.01)
  expr_base <- rep(2, length(expr_diff))

  br <- outer(expr_diff + expr_base, cn_diff + cn_base,
              FUN = \(x,y) br_func(expr_obs = x, cn_obs = y,
                                   expr_base = expr_base, cn_base = cn_base))
  plot_ly(y = ~expr_diff, x = ~cn_diff, z = ~br, type = 'surface')
}

plot_buffering_ratio_expr <- function(br_func, expr_lim = c(-1, 1), cn_diff = 1,
                                      buffered_threshold = br_cutoffs$Buffered,
                                      anti_scaling_threshold = br_cutoffs$AntiScaling) {
  expr_diff <- seq(expr_lim[1], expr_lim[2], by = 0.01)
  expr_base <- rep(2, length(expr_diff))
  cn_base <- rep(2, length(expr_diff))
  cn_diff <- rep(cn_diff, length(expr_diff))
  br_values <- br_func(expr_obs = expr_diff + expr_base, cn_obs = cn_diff + cn_base,
                       expr_base = expr_base, cn_base = cn_base)

  ggplot() +
    aes(x = expr_diff, y = br_values) +
    geom_line() +
    geom_hline(yintercept = buffered_threshold, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = anti_scaling_threshold, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = expr_lim, breaks = seq(expr_lim[1], expr_lim[2], 0.2)) +
    labs(x = "Expression Difference", y = "Buffering Ratio")
}

plot_buffering_ratio_cn <- function(br_func, expr_diff = 1, cn_lim = c(-1, 1), tick_distance = 0.1,
                                    buffered_threshold = br_cutoffs$Buffered,
                                    anti_scaling_threshold = br_cutoffs$AntiScaling) {
  cn_diff <- seq(cn_lim[1], cn_lim[2], by = 0.01)
  cn_base <- rep(2, length(cn_diff))
  expr_diff <- rep(expr_diff, length(cn_diff))
  expr_base <- rep(2, length(expr_diff))
  br_values <- br_func(expr_obs = expr_diff + expr_base, cn_obs = cn_diff + cn_base,
                       expr_base = expr_base, cn_base = cn_base)

  ggplot() +
    aes(x = cn_diff, y = br_values) +
    geom_line() +
    geom_hline(yintercept = buffered_threshold, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = anti_scaling_threshold, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = cn_lim, breaks = seq(cn_lim[1], cn_lim[2], tick_distance)) +
    labs(x = "Copy Number Difference", y = "Buffering Ratio")
}

buffering_example <- function() {
  eps <- 1e-8

  df1 <- data.frame(ExprBase = c(10, 10, 10, 10), ExprVar = c(10, 15, 20, 30),
                    CNBase = c(2, 2, 2, 2), CNVar = c(3, 3, 3, 3))
  df2 <- data.frame(ExprBase = c(10, 15, 20, 30), ExprVar = c(10, 10, 10, 10),
                    CNBase = c(2, 2, 2, 2), CNVar = c(1, 1, 1, 1))
  df3 <- data.frame(ExprBase = c(10, 15, 20, 30), ExprVar = c(10, 10, 10, 10),
                    CNBase = c(3, 3, 3, 3), CNVar = c(2, 2, 2, 2))
  df4 <- data.frame(ExprBase = c(10, 10, 10, 10), ExprVar = c(10, 15, 20, 40),
                    CNBase = c(2, 2, 2, 2), CNVar = c(4, 4, 4, 4))

  for (df in list(df1, df2, df3, df4)) {
    df %>%
      mutate(ExprVar = ExprVar + runif(length(ExprVar), min = -eps, max = eps)
      ) %>%
      mutate(Log2FC = log2(ExprVar) - log2(ExprBase),
             BR = buffering_ratio(ExprBase, ExprVar, CNBase, CNVar),
             SF = scaling_factor(ExprBase, ExprVar, CNBase, CNVar),
             SR = scaling_ratio(ExprBase, ExprVar, CNBase, CNVar),
             LBS = log_buffering_score(ExprBase, ExprVar, CNBase, CNVar)
      ) %>%
      mutate(BC_log2fc = buffering_class_log2fc(Log2FC, CNBase, CNVar),
             BC_br = buffering_class(BR),
             BC_sf = buffering_class_sf(SF),
             BC_sr = buffering_class_sr(SR),
             BC_lbs = buffering_class_lbs(LBS)
      ) %>%
      mutate_if(is.numeric, ~round(., digits = 5)) %>%
      print(digits = 5)
  }
}

# scaling_ratio(expr_base = 2, expr_obs = 1.149, cn_base = 2, cn_obs = 1) = 0.7996212
# scaling_ratio(expr_base = 1, expr_obs = 1.741, cn_base = 1, cn_obs = 2) = 0.7999162

# Strange bahaviour for scaling ratio & scaling factor
# Measure should show stronger Anti-Scaling if expression change is positive and the loss of copy number is stronger
# Example: Positive expression change (2->3), copy number loss (2->1 vs. 2->0.2)
## scaling_ratio(expr_base = 2, expr_obs = 3, cn_base = 2, cn_obs = 1) = -1
## scaling_ratio(expr_base = 2, expr_obs = 3, cn_base = 2, cn_obs = 0.2) = -0.5555556
### SR = -0.555 > -1, indicating less buffering/anti-scaling)
# Buffering Ratio works as expected
## buffering_ratio(expr_base = 2, expr_obs = 3, cn_base = 2, cn_obs = 1) = 1.584963
## buffering_ratio(expr_base = 2, expr_obs = 3, cn_base = 2, cn_obs = 0.2) = 3.906891
### BR = 3.9 > 1.6 indicating more buffering/anti-scaling

# Conclusion: SR & SF not adequate for quantifying anti-scaling
