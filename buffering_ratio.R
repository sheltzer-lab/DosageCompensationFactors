library(dplyr)
library(plotly)

buffering_ratio <- function(expr_base, expr_var, cn_base = 2, cn_var = 3) {
  br <- log2(cn_var / cn_base) - log2(expr_var / expr_base)
  ifelse(cn_var > cn_base, br, -br)
}

buffering_ratio_old <- function(expr_base, expr_var, cn_base = 2, cn_var = 3) {
  br <- 1 - (expr_var / expr_base) * (cn_base / cn_var)
  ifelse(cn_var > cn_base, br, -br)
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

plot_buffering_ratio <- function(br_func, cnv_lim = c(-1, 1), expr_lim = c(-1, 1)) {
  cn_base <- 2
  cn_diff <- seq(cnv_lim[1], cnv_lim[2], by = 0.01)
  expr_base <- 2
  expr_diff <- seq(expr_lim[1], expr_lim[2], by = 0.01)

  # ToDo: Something is wrong here
  br <- outer(expr_diff + expr_base, cn_diff + cn_base,
              FUN = \(x,y) br_func(expr_var = x, cn_var = y,
                                   expr_base = expr_base, cn_base = cn_base))
  plot_ly(x = ~expr_diff, y = ~cn_diff, z = ~br, type = 'surface')
}

buffering_example <- function() {
  df1 <- data.frame(ExprBase = c(10, 10, 10, 10), ExprVar = c(10, 15, 20, 30),
                    CNBase = c(2, 2, 2, 2), CNVar = c(3, 3, 3, 3))
  df2 <- data.frame(ExprBase = c(10, 15, 20, 30), ExprVar = c(10, 10, 10, 10),
                    CNBase = c(3, 3, 3, 3), CNVar = c(2, 2, 2, 2))
  df3 <- data.frame(ExprBase = c(10, 10, 10, 10), ExprVar = c(10, 15, 20, 40),
                    CNBase = c(2, 2, 2, 2), CNVar = c(4, 4, 4, 4))

  for (df in list(df1, df2, df3)) {
    df_buff <- df %>%
      mutate(BR = buffering_ratio(ExprBase, ExprVar, CNBase, CNVar)) %>%
      mutate(BC = buffering_class(BR)) %>%
      mutate(Log2FC = log2(ExprVar) - log2(ExprBase)) %>%
      mutate(Log2FC_BC = buffering_class_log2fc(Log2FC, CNBase, CNVar))
    print(df_buff)
  }
}