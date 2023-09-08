library(dplyr)
library(plotly)
library(ggplot2)

buffering_ratio <- function(expr_base, expr_var, cn_base = 2, cn_var = 3) {
  br <- log2(cn_var / cn_base) - log2(expr_var / expr_base)
  ifelse(cn_var > cn_base, br, -br)
}

buffering_ratio_old <- function(expr_base, expr_var, cn_base = 2, cn_var = 3) {
  br <- 1 - (expr_var / expr_base) * (cn_base / cn_var)
  ifelse(cn_var > cn_base, br, -br)
}

scaling_factor <- function(expr_base, expr_var, cn_base = 2, cn_var = 3) {
  sf <- ((expr_var / expr_base) - 1) / ((cn_var / cn_base) - 1)
  ifelse(is.finite(sf), sf, NA)
}

buffering_class <- function(buffering_ratio) {
  ifelse(buffering_ratio > 0.6849625, "Anti-Scaling",        # Buffering of trisomy expression below disomy level
         ifelse(buffering_ratio > 0.3349625, "Buffered",    # Less expression in trisomy than expected
                "Scaling"))                                  # Expression level as expected or higher
}

buffering_class_old <- function(buffering_ratio) {
    ifelse(buffering_ratio > 0.377978, "Anti-Scaling",        # Buffering of trisomy expression below disomy level
         ifelse(buffering_ratio > 0.2071953, "Buffered",    # Less expression in trisomy than expected
                "Scaling"))                                  # Expression level as expected or higher
}

buffering_class_sf <- function(scaling_factor) {
    ifelse(scaling_factor > 0.3784142, "Scaling",
         ifelse(scaling_factor > -0.133934, "Buffered",
                "Anti-Scaling"))
}

# Values from Schukken & Sheltzer, 2022 (DOI: 10.1101/gr.276378.121) on chromosome arm gain (inverted for arm loss):
# Scaling:        0.25 < Log2FC
# Buffered:      -0.1  < Log2FC <  0.25
# Anti-Scaling:          Log2FC < -0.1
buffering_class_log2fc <- function (log2fc, cn_base = 2, cn_var = 3) {
  ifelse(cn_var > cn_base,
         ifelse(log2fc > 0.25, "Scaling",
                ifelse(log2fc > -0.1, "Buffered",
                       "Anti-Scaling")),
         ifelse(log2fc < -0.25, "Scaling",
                ifelse(log2fc < 0.1, "Buffered",
                       "Anti-Scaling")))
}

# === Example Code ===

plot_buffering_ratio_3d <- function(br_func, cnv_lim = c(-0.99, 0.99), expr_lim = c(-1, 1)) {
  cn_diff <- seq(cnv_lim[1], cnv_lim[2], by = 0.01)
  cn_base <- rep(2, length(cn_diff))
  expr_diff <- seq(expr_lim[1], expr_lim[2], by = 0.01)
  expr_base <- rep(2, length(expr_diff))

  br <- outer(expr_diff + expr_base, cn_diff + cn_base,
              FUN = \(x,y) br_func(expr_var = x, cn_var = y,
                                   expr_base = expr_base, cn_base = cn_base))
  plot_ly(y = ~expr_diff, x = ~cn_diff, z = ~br, type = 'surface')
}

plot_buffering_ratio_expr <- function(br_func, expr_lim = c(-1, 1), cn_diff = 1,
                                 buffered_threshold = 0.3349625, anti_scaling_threshold = 0.6849625) {
  expr_diff <- seq(expr_lim[1], expr_lim[2], by = 0.01)
  expr_base <- rep(2, length(expr_diff))
  cn_base <- rep(2, length(expr_diff))
  cn_diff <- rep(cn_diff, length(expr_diff))
  br_values <- br_func(expr_var = expr_diff + expr_base, cn_var = cn_diff + cn_base,
                       expr_base = expr_base, cn_base = cn_base)

  ggplot() +
    aes(x = expr_diff, y = br_values) +
    geom_line() +
    geom_hline(yintercept = buffered_threshold, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = anti_scaling_threshold, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = expr_lim, breaks = seq(expr_lim[1], expr_lim[2], 0.2))
}

plot_buffering_ratio_cn <- function(br_func, expr_diff = 1, cn_lim = c(-1, 1),
                                    buffered_threshold = 0.3349625, anti_scaling_threshold = 0.6849625) {
  cn_diff <- seq(cn_lim[1], cn_lim[2], by = 0.01)
  cn_base <- rep(2, length(cn_diff))
  expr_diff <- rep(expr_diff, length(cn_diff))
  expr_base <- rep(2, length(expr_diff))
  br_values <- br_func(expr_var = expr_diff + expr_base, cn_var = cn_diff + cn_base,
                       expr_base = expr_base, cn_base = cn_base)

  ggplot() +
    aes(x = cn_diff, y = br_values) +
    geom_line() +
    geom_hline(yintercept = buffered_threshold, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = anti_scaling_threshold, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = cn_lim, breaks = seq(cn_lim[1], cn_lim[2], 0.1))
}

buffering_example <- function() {
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
      mutate(Log2FC = log2(ExprVar) - log2(ExprBase),
             BR = buffering_ratio(ExprBase, ExprVar, CNBase, CNVar),
             BR_old = buffering_ratio_old(ExprBase, ExprVar, CNBase, CNVar),
             SF = scaling_factor(ExprBase, ExprVar, CNBase, CNVar)) %>%
      mutate(Log2FC_BC = buffering_class_log2fc(Log2FC, CNBase, CNVar),
             BC = buffering_class(BR),
             BC_old = buffering_class_old(BR_old),
             BC_sf = buffering_class_sf(SF)) %>%
      print()
  }
}