library(ggplot2)
library(dplyr)
library(rlang)
library(EnhancedVolcano)
library(psych)
library(corrplot)
library(ggsignif)
library(viridisLite)
library(viridis)
library(ggbeeswarm)
library(pheatmap)
library(forcats)
library(ggrepel)

here::i_am("DosageCompensationFactors.Rproj")

vertical_bar_chart <- function(df, category_col, value_col,
                               color_col = NULL, default_fill_color = default_color, text_color = NULL,
                               error_low_col = NULL, error_high_col = NULL,
                               value_range = c(0.45, 0.65), break_steps = 0.05,
                               line_intercept = 0.5, bar_label_shift = 0.002, color_discrete = FALSE,
                               title = NULL, category_lab = NULL, value_lab = NULL, color_lab = NULL) {
  require(ggplot2)
  require(shadowtext)

  df %>%
    ggplot() +
    aes(x = { { category_col } }, y = { { value_col } },
        label = format(round({ { value_col } }, 3), nsmall = 3)) +
    { if (!quo_is_null(enquo(color_col))) geom_bar(aes(fill = {{ color_col }}), stat = "identity")
      else geom_bar(stat = "identity", fill = default_fill_color) } +
    geom_hline(yintercept = line_intercept) +
    { if (!quo_is_null(enquo(error_low_col)) & !quo_is_null(enquo(error_high_col)))
      geom_pointrange(aes(x = { { category_col } }, y = { { value_col } },
                          ymin = { { error_low_col } }, ymax = { { error_high_col } }),
                      colour = "orange", fatten = 1) } +
    { if (!is.null(text_color)) geom_text(color = text_color, y = value_range[1] + bar_label_shift, hjust = 0)
      else geom_shadowtext(color = "white", y = value_range[1] + bar_label_shift, hjust = 0, bg.colour = default_color) } +
    scale_y_continuous(breaks = seq(value_range[1], value_range[2], break_steps)) +
    scale_fill_viridis(option = "G", direction = -1, begin = 0.2, end = 0.8, discrete = color_discrete) +
    labs(title = title, x = category_lab, y = value_lab, fill = color_lab) +
    coord_flip(ylim = c(value_range[1], value_range[2]))
}

plot_text_col <- function(df, x_col, label_col, align = "center") {
  align_param <- case_when(
    align == "left" ~ 0,
    align == "center" ~ 0.5,
    align == "right" ~ 1,
    TRUE ~ 0.5
  )

  df %>%
    ggplot() +
    aes(x = { { x_col } }, y = align_param, label = { { label_col } }) +
    geom_text(color = "black", hjust = align_param) +
    xlab("") +
    ylab("") +
    coord_flip(ylim = c(0, 1)) +
    theme_void() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

violin_plot <- function(df, x, y) {
  df %>%
    ggplot() +
    aes(x = { { x } }, y = { { y } }) +
    geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75),
                color = "#4080DB") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

beeswarm_plot <- function(df, x, y, color_col = NULL, cex = 2) {
  df %>%
    ggplot() +
    aes(x = { { x } }, y = { { y } }) +
    geom_violin(trim = FALSE, draw_quantiles = 0.5,
                fill = "darkgrey", color = "darkgrey", alpha = 1/3) +
    {
      if (quo_is_null(enquo(color_col))) {
        geom_beeswarm(priority = "density", color = "#4080DB", cex = cex)
      } else {
        geom_beeswarm(aes(color = { { color_col } }), priority = "density", cex = cex)
      }
    } +
    scale_colour_viridis_c(option = "D", direction = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_correlation <- function(df, method = "spearman", adjust = "none") {
  if(nrow(df) == 0) return(NULL)

  cor_matrix <- psych::corr.test(df, method = method, adjust = adjust)

  corrplot(cor_matrix$r, p.mat = cor_matrix$p,
           type = "upper", order = "hclust",
           tl.col = "black", tl.srt = 45,
           pch.cex = 1, pch.col = "darkgrey")
}


plot_volcano_buffered <- function(df, ratio_col, signif_col, label_col, class_col,
                                  color_mapping = NULL,
                                  value_threshold = log2fc_threshold, signif_threshold = p_threshold,
                                  title = NULL, subtitle = NULL) {
  buffering_levels <- c("N.S.", "Scaling", "Buffered", "Anti-Scaling")
  color_palette <- c("darkgrey", "#0E992C", "#E6C45A", "#E64356")
  names(color_palette) <- buffering_levels
  color_mapping <- scale_colour_manual(name = "Buffering Class", values = color_palette)

  df %>%
    mutate(Buffering.Class = if_else(Test.p.adjusted > p_threshold, "N.S.", Buffering.Class)) %>%
    mutate(Buffering.Class = factor({ { class_col } }, levels = buffering_levels)) %>%
    mutate(Label = if_else(Buffering.Class %in% c("Buffered", "Anti-Scaling"),
                           { { label_col } }, NA)) %>%
    plot_volcano({ { ratio_col } }, { { signif_col } }, Label, Buffering.Class,
                 color_mapping = color_mapping,
                 value_threshold = value_threshold, signif_threshold = signif_threshold,
                 title = title, subtitle = subtitle)
}

plot_volcano <- function(df, value_col, signif_col, label_col, color_col,
                         color_mapping = scale_color_manual(values = color_palettes$DiffExp,
                                                            na.value = color_palettes$Missing),
                         value_threshold = log2fc_threshold, signif_threshold = p_threshold,
                         title = NULL, subtitle = NULL) {
  max_abs_value <- df %>% pull({ { value_col } }) %>% abs() %>% max(na.rm = TRUE)

  df %>%
    mutate(`-log10(p)` = -log10({ { signif_col } })) %>%
    ggplot() +
    aes(x = { { value_col } }, y = `-log10(p)`,
        label = { { label_col } }, color = { { color_col } }) +
    geom_point(alpha = 0.5, size = 1) +
    color_mapping +
    geom_hline(yintercept = -log10(signif_threshold),
               linetype = "dashed", color = "black") +
    geom_vline(xintercept = abs(value_threshold),
               linetype = "dashed", color = "black") +
    geom_vline(xintercept = -abs(value_threshold),
               linetype = "dashed", color = "black") +
    geom_label_repel(min.segment.length = 0.01, label.size = 0.15,
                     seed = 42, max.iter = 30000, max.time = 1.5,
                     point.padding = 0.3, label.padding = 0.3, box.padding = 0.3,
                     force = 5, max.overlaps = 20) +
    labs(title = title, subtitle = subtitle) +
    xlim(c(-max_abs_value, max_abs_value))
}

# Get density of points in 2 dimensions.
# https://slowkow.com/notes/ggplot2-color-by-density/
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Source: https://rpubs.com/Bio-Geek/71339, https://slowkow.com/notes/ggplot2-color-by-density/
scatter_plot_regression <- function(df, x_col, y_col, formula, color_col = NULL, x_lab = NULL, y_lab = NULL,
                                    label_coords = NULL, title = NULL, point_size = 1) {
  x_col_name <- quo_name(enquo(x_col))
  y_col_name <- quo_name(enquo(y_col))

  eps <- 1e-8

  df <- df %>%
    select({ { x_col } }, { { y_col } }, { { color_col } }) %>%
    drop_na() %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    mutate(Density = get_density({ { x_col } }, { { y_col } },
                                 # Ensure that bandwidth is stricly positive
                                 h = c(MASS::bandwidth.nrd({ { x_col } }) + eps,
                                       MASS::bandwidth.nrd({ { y_col } }) + eps),
                                 n = 100))

  suppressWarnings({
    regression <- lm(formula, df)
    pred.int <- predict(regression, interval = "prediction")
    regression_summary <- summary(regression)
    df <- cbind(df, pred.int)
    slope <- regression$coefficients[[x_col_name]]
    intercept <- regression$coefficients[["(Intercept)"]]
  })

  if (is.null(x_lab)) x_lab <- x_col_name
  if (is.null(y_lab)) y_lab <- y_col_name

  x_range <- c(min(df[[x_col_name]]), max(df[[x_col_name]]))
  y_range <- c(min(df[[y_col_name]]), max(df[[y_col_name]]))

  if (is.null(label_coords)) {
    label_coords <- list(x = mean(x_range),
                         y = c(y_range[1] + 0.05 * abs(y_range[2] - y_range[1]),
                               y_range[1] + 0.01 * abs(y_range[2] - y_range[1])))
  }

  regression_plot <- df %>%
    ggplot() +
    {
      if (quo_is_null(enquo(color_col))) aes(x = { { x_col } }, y = { { y_col } }, color = Density)
      else aes(x = { { x_col } }, y = { { y_col } }, color = { { color_col } })
    } +
    geom_point(alpha = 0.8, size = point_size) +
    stat_smooth(method = lm, color = "blue") +
    geom_line(aes(y = lwr), color = "red", linetype = "dashed") +
    geom_line(aes(y = upr), color = "red", linetype = "dashed") +
    annotate("text", x = label_coords$x, y = label_coords$y, color = "blue", parse = TRUE,
             label = c(paste0("y == ", signif(slope, 2), " * x + ", signif(intercept, 2)),
                       paste0("R^{2} == ", signif(regression_summary$r.squared, 2)))) +
    {
      if (quo_is_null(enquo(color_col))) scale_colour_viridis_c(option = "D", direction = 1)
    } +
    labs(x = x_lab, y = y_lab, title = title)

  return(regression_plot)
}

scatter_plot_reg_corr <- function(df, x_col, y_col, color_col = NULL, cor_method = "spearman",
                                  cor_symbol = utf8_rho, label_coords = NULL, point_size = 1,
                                  x_lab = NULL, y_lab = NULL, title_prefix = NULL) {
  df_ <- df %>%
    mutate(x = { { x_col } }, y = { { y_col } })

  cor_result <- df_ %>%
    rstatix::cor_test(x, y, method = cor_method)

  if (is.null(title_prefix)) {
    plot_title <- paste0("Correlation: ",
                         print_corr(cor_result$cor, p.value = cor_result$p, signif = TRUE, estimate_symbol = cor_symbol))
  } else {
    plot_title <- paste0(title_prefix, " (",
                         print_corr(cor_result$cor, p.value = cor_result$p, signif = TRUE, estimate_symbol = cor_symbol), ")")
  }

  if (is.null(x_lab)) x_lab <- quo_name(enquo(x_col))
  if (is.null(y_lab)) y_lab <- quo_name(enquo(y_col))

  df_ %>%
    scatter_plot_regression(x, y, y ~ x, color_col = { { color_col } }, x_lab = x_lab, y_lab = y_lab,
                            label_coords = label_coords, title = plot_title, point_size = point_size)
}

waterfall_plot <- function(df, value_col, rank_col, label_col, n = 5,
                           color_low = tail(bidirectional_color_pal, n = 1),
                           color_high = head(bidirectional_color_pal, n = 1),
                           centrality_measure = mean,
                           y_margin = 0.02, font_size = theme_get()$text$size / 4) {
  rank_col_name <- quo_name(enquo(rank_col))
  value_col_name <- quo_name(enquo(value_col))

  xlim <- c(0, max(df[[rank_col_name]]))
  label_nudge_x <- floor(xlim[2] / 4)

  y_range <- max(df[[value_col_name]]) - min(df[[value_col_name]])
  ylim1 <- c(min(df[[value_col_name]] - y_range * y_margin),
             (df %>% filter({ { rank_col } } == label_nudge_x))[[value_col_name]])
  ylim2 <- c((df %>% filter({ { rank_col } } == xlim[2] - label_nudge_x))[[value_col_name]],
             max(df[[value_col_name]] + y_range * y_margin))

  central_value <- centrality_measure(df[[value_col_name]], na.rm = TRUE)

  df %>%
    ggplot() +
    aes(x = { { rank_col } }, y = { { value_col } }, label = { { label_col } }) +
    geom_hline(yintercept = 0, color = default_color) +
    geom_point(size = 0.3) +
    geom_hline(yintercept = central_value, color = "red") +
    geom_text_repel(data = df %>% slice_min({ { rank_col } }, n = n),
                    xlim = xlim, ylim = ylim1, direction = "y", nudge_x = label_nudge_x,
                    seed = 42, force = 2.2, color = color_low, size = font_size, hjust = 0) +
    geom_text_repel(data = df %>% slice_max({ { rank_col } }, n = n),
                    xlim = xlim, ylim = ylim2, direction = "y", nudge_x = -label_nudge_x,
                    seed = 42, force = 2.2, color = color_high, size = font_size, hjust = 1) +
    lims(y = c(ylim1[1], ylim2[2]))
}

save_plot <- function(plot, filename, dir = plots_dir,
                      height = 150, width = 150, dpi = 300) {
  ggsave(here(dir, filename), plot = plot,
       height = height, width = width, units = "mm", dpi = dpi)
  return(plot)
}

plot_rocs <- function(df_rocs, legend_position = "right", legend_rows = 10, label_padding = 0.05) {
  df_label <- df_rocs %>%
    distinct(Name, AUC) %>%
    arrange(desc(Name)) %>%
    mutate(x = 0.1,
           y = seq(0.05, 0.75, label_padding)[seq_along(unique(df_rocs$Name))],
           AUC = paste0("AUC = ", format(round(AUC, 3), nsmall = 3)))

  plot <- df_rocs %>%
    arrange(Sensitivity) %>%
    ggplot() +
    aes(x = Specificity, y = Sensitivity, color = Name) +
    geom_abline(slope = 1, intercept = 1, color = "grey") +
    geom_line() +
    geom_label(data = df_label, mapping = aes(color = Name, label = AUC, x = x, y = y)) +
    scale_x_reverse(limits = c(1, 0)) +
    labs(x = "Specificity", y = "Sensitivity", color = "Model") +
    theme(legend.position = legend_position,
          legend.text = element_text(size = 10)) +
    guides(color = guide_legend(nrow = legend_rows))

  return(plot)
}

print_signif <- function(p, digits = 3, scientific = TRUE, digits_scientific = 1) {
  ifelse(!is.na(p) & p == 0,
         "p < 1.0e-99",
         ifelse(scientific,
                paste0("p = ", formatC(p, format = "e", digits = digits_scientific)),
                paste0("p ", if_else(p < 10^(-digits),
                                     paste0("< ", format(10^(-digits), nsmall = digits, scientific = FALSE)),
                                     paste0("= ", format(round(p, digits), nsmall = digits, scientific = FALSE))))
         )
  )
}

print_corr <- function(corr, p.value = NULL, estimate_symbol = utf8_rho, signif = FALSE, digits = 3, map_p = FALSE) {
  corr_str <- paste(estimate_symbol, "=", format(round(corr, digits), nsmall = digits, scientific = FALSE))
  if (signif == FALSE) return(corr_str)
  else {
    p_str <- print_signif(p.value, digits = digits)
    if (map_p == TRUE) p_str <- map_signif(p.value)
    return(paste0(corr_str, ", ", p_str))
  }
}

print_corr_obj <- function (corr, estimate_symbol = utf8_rho, signif = TRUE, digits = 3, map_p = FALSE) {
  return(print_corr(corr$estimate, p.value = corr$p.value,
                    estimate_symbol = estimate_symbol, signif = signif, digits = digits, map_p = map_p))
}

map_signif <- function (p, thresholds = c(0.01, 0.001, 0.0001)) {
  case_when(
      p < thresholds[3] ~ "***",
      p < thresholds[2] ~ "**",
      p < thresholds[1] ~ "*",
      TRUE ~ "n.s."
    )
}

plot_corr_bracket <- function(corr, estimate_symbol = utf8_rho, signif = TRUE, digits = 3, map_p = FALSE,
                              margin = 2, shift = 0, size = default_theme$text$size / 2) {
  x_range <- c(0,10)
  y_range <- c(0,3)
  x_range_adj <- c(x_range[1] + margin + shift, x_range[2] - margin)

  data.frame(Label = print_corr_obj(corr, estimate_symbol, signif, digits, map_p)) %>%
    ggplot() +
    aes(x = x_range[1], y = y_range[1], label = Label) +
    geom_segment(aes(x = x_range_adj[1], y = 1,
                     xend = x_range_adj[2], yend = 1)) +
    geom_text(x = mean(x_range_adj), color = "black", y = 2, size = size) +
    lims(x = x_range, y = y_range) +
    cowplot::theme_nothing()
}

vertical_box_plot <- function(df, value_col, id_col, group_col) {
  df %>%
    arrange({ { id_col } }) %>%
    ggplot() +
    aes(x = { { value_col } }, y = { { id_col } }, fill = { { group_col } }) +
    geom_boxplot()
}

jittered_boxplot <- function(df, group_col, value_col, color_col = NULL, alpha = 0.5, jitter_width = 0.15) {
  df %>%
    ggplot() +
    aes(x = { { group_col } }, y = { { value_col } }) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    {
      if (quo_is_null(enquo(color_col))) {
        geom_quasirandom(fill = "darkgrey", color = "white",
                         shape = 21, alpha = alpha, width = jitter_width)
      } else {
        geom_quasirandom(aes(color = { { color_col } }), alpha = alpha, width = jitter_width)
      }
    } +
    coord_flip() +
    scale_colour_viridis_c(option = "D", direction = 1)
    #scale_color_gradientn(colors = biderectional_color_pal, space = "Lab")
}

plot_pca <- function(pca, color_col = CellLine.Name, label_col = Sample.Name) {
  eigenvalues <- pca$eigenvalues
  df_pca <- pca$df_pca

  df_pca %>%
    ggplot() +
    aes(.fittedPC1, .fittedPC2, color = { { color_col } }, label = { { label_col } }) +
    geom_point(size = 1.5) +
    {
      if (!quo_is_null(enquo(label_col))) geom_label_repel(force = 10, seed = 123)
    } +
    xlab(sprintf("PC1 (%0.1f%% Variance Explained)", eigenvalues[eigenvalues$PC == 1,]$percent * 100)) +
    ylab(sprintf("PC2 (%0.1f%% Variance Explained)", eigenvalues[eigenvalues$PC == 2,]$percent * 100))
}

scree_plot <- function(pca) {
  pca$eigenvalues %>%
    ggplot(aes(PC, percent)) +
    geom_col() +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = expansion(mult = c(0, 0.01))
    ) +
    xlab("Principal Component") +
    ylab("Explained Variance")
}

bidirectional_heatmap <- function(df, value_col, sample_col, group_col, text_col = NULL,
                                  cluster_rows = FALSE, cluster_cols = FALSE,
                                  show_rownames = TRUE, show_colnames = TRUE,
                                  transpose = FALSE, palette_length = 100, color_pal = bidirectional_color_pal) {
  mat <- df %>%
    group_by({ { sample_col } }, { { group_col } }) %>%
    summarize(Value = mean({ { value_col } }, na.rm = TRUE), .groups = "drop") %>%
    arrange({ { sample_col } }, { { group_col } }) %>%
    pivot_wider(names_from = { { sample_col } }, values_from = Value, id_cols = { { group_col } }) %>%
    tibble::column_to_rownames(var = quo_name(enquo(group_col)))

  # Add text labels
  text_mat <- matrix(data = "", nrow = nrow(mat), ncol = ncol(mat))
  if (!quo_is_null(enquo(text_col))) {
    text_mat <- df %>%
      select({ { sample_col } }, { { group_col } }, { { text_col } }) %>%
      distinct() %>%
      arrange({ { sample_col } }, { { group_col } }) %>%
      pivot_wider(names_from = { { sample_col } }, values_from = { { text_col } }, id_cols = { { group_col } }) %>%
      tibble::column_to_rownames(var = quo_name(enquo(group_col)))
  }

  if (transpose) {
    mat <- t(mat)
    text_mat <- t(text_mat)
  }

  # https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
  color <- colorRampPalette(color_pal, space = "Lab")(palette_length)
  # use floor and ceiling to deal with even/odd length pallettelengths
  breaks <- c(seq(min(mat, 0 - 10e-10, na.rm = T), 0,
                  length.out = ceiling(palette_length / 2) + 1),
              seq(max(mat, 0 + 10e-10, na.rm = T) / palette_length, max(mat, na.rm = T),
                  length.out = floor(palette_length / 2)))

  pheatmap(mat, na_col = "black",
           show_colnames = show_colnames, show_rownames = show_rownames,
           cluster_cols = cluster_cols, cluster_rows = cluster_rows,
           display_numbers = text_mat, breaks = breaks, color = color)
}

simple_heatmap <- function(df, x_col, y_col, color_col, label_col,
                           x_lab = NULL, y_lab = NULL, legend_lab = NULL) {
  theme_settings <- theme(legend.key.size = unit(16, "points"),
                          legend.key.width = unit(24, "points"),
                          legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10),
                          legend.position = "top",
                          legend.direction = "horizontal",
                          axis.text.x = element_text(angle = 45, hjust = 1))

  df %>%
    distinct({ { x_col } }, { { y_col } }, { { color_col } }, { { label_col } }) %>%
    ggplot() +
    aes(x = { { x_col } }, y = { { y_col } }, fill = { { color_col } }, label = { { label_col } }) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab",
                         limits = c(-1, 1), oob = scales::squish) +
    labs(x = x_lab, y = y_lab, fill = legend_lab) +
    cowplot::theme_minimal_grid() +
    theme_settings
}

unidirectional_heatmap <- function(df, x_col, y_col, color_col, order_desc = FALSE, color_lab = NULL) {
  theme_settings <- theme(legend.key.size = unit(16, "points"),
                        legend.key.width = unit(24, "points"),
                        legend.title = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.position = "top",
                        legend.direction = "horizontal",
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())

  color_col_name <- quo_name(enquo(color_col))
  x_col_name <- quo_name(enquo(x_col))

  if (is.null(color_lab)) color_lab <- color_col_name

  max_value <- round(max(df[[color_col_name]]), digits = 1)

  df %>%
    group_by({ { x_col } }, { { y_col } }) %>%
    summarize(!!color_col_name := median({ { color_col } }, na.rm = TRUE), .groups = "drop") %>%
    mutate(!!x_col_name := fct_reorder({ { x_col } }, { { color_col } }, .desc = order_desc)) %>%
    ggplot() +
    aes(x = { { x_col } }, y = { { y_col } },
        fill = { { color_col } }, label = format(round({ { color_col } }, 2), nsmall = 2, scientific = FALSE)) +
    geom_tile() +
    geom_text(color = "white") +
    scale_fill_viridis_c(option = "magma", direction = 1, end = 0.9,
                         limits = c(0.5, max_value), oob = scales::squish) +
    cowplot::theme_minimal_grid() +
    labs(fill = color_lab) +
    theme_settings
}

roc_auc_heatmap <- function(df, rank_tests, color_lab = NULL) {
  heatmap_boot_auc <- df %>%
    unidirectional_heatmap(DosageCompensation.Factor, Condition, DosageCompensation.Factor.ROC.AUC,
                           order_desc = TRUE, color_lab = color_lab)

  signif_bars <- data.frame(Condition = unique(df$Condition), y = c(2, 2, 2, 2)) %>%
    ggplot() +
    aes(x = Condition, y = y) +
    geom_signif(comparisons = rank_tests$comparisons,
                annotations = rank_tests$annotations,
                y_position = rank_tests$y_position,
                size = 1, textsize = 4, vjust = -0.5
    ) +
    cowplot::theme_nothing() +
    coord_flip(ylim = c(1, 2.5))

  cowplot::plot_grid(
    heatmap_boot_auc, signif_bars,
    labels = NULL, nrow = 1, align = "h", axis = "tb",
    rel_widths = c(12, 1)
  )
}

bucketed_scatter_plot <- function(df, value_col, x_value_col, bucket_col,
                                  threshold_low = NULL, threshold_high = NULL, value_range = c(-2, 2),
                                  highlight_buckets = NULL, x_lab = NULL, title = element_blank()) {
  df %>%
    group_by({ { bucket_col } }) %>%
    mutate(Position = ({ { x_value_col } } - min({ { x_value_col } }, na.rm = TRUE)) /
      max({ { x_value_col } }, na.rm = TRUE),
           Value.Average = mean({ { value_col } }, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Bucket = { { bucket_col } },
           Color = if_else({ { bucket_col } } %in% highlight_buckets,
                           highlight_color, default_color)) %>%
    ggplot() +
    aes(x = Position, y = { { value_col } }, color = Color) +
    geom_hline(yintercept = threshold_low, color = bidirectional_color_pal[1], linetype="dashed", linewidth = 1/3) +
    geom_hline(yintercept = threshold_high, color = "orange", linetype="dashed", linewidth = 1/3) +
    geom_hline(yintercept = 0, color = default_color) +
    geom_point(alpha = 1/2, size = 1/3) +
    geom_hline(aes(yintercept = Value.Average), color = "red") +
    facet_grid(~Bucket) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) +
    scale_colour_identity() +
    scale_y_continuous(limits = value_range, breaks = seq(c(value_range, 1))) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0.25, "mm")) +
    xlab(x_lab) +
    ggtitle(title)
}

sorted_violin_plot <- function(df, x, y) {
  plot <- df %>%
    add_count(get(x)) %>%
    filter(n > 2) %>%
    group_by(get(x)) %>%
    mutate(Median = median({ { y } }),
           Label = paste0(get(x), " (n=", n, ")")) %>%
    ungroup() %>%
    arrange(Median) %>%
    mutate(Label = factor(Label, levels = unique(Label))) %>%
    violin_plot(Label, { { y } })
  plot <- plot +
    xlab(x)

  return(plot)
}

sorted_beeswarm_plot <- function(df, x, y, color_col = NULL, cex = 2) {
  plot <- df %>%
    add_count(get(x)) %>%
    filter(n > 2) %>%
    group_by(get(x)) %>%
    mutate(Median = median({ { y } }),
           Label = paste0(get(x), " (n=", n, ")")) %>%
    ungroup() %>%
    arrange(Median) %>%
    mutate(Label = factor(Label, levels = unique(Label))) %>%
    beeswarm_plot(Label, { { y } }, color_col = { { color_col } }, cex = cex)
  plot <- plot +
    xlab(x)

  return(plot)
}

prep_signif <- function(df, x, facet_col = NULL) {
  df_prep <- df %>%
    {
      if (!quo_is_null(enquo(facet_col))) {
        mutate(., Bucket = { { facet_col } }) %>%
          group_by(Bucket)
      } else {
        .
      }
    } %>%
    add_count({ { x } }) %>%
    filter(n > 2) %>%
    ungroup() %>%
    { if (!quo_is_null(enquo(facet_col))) group_by(., Bucket, { { x } })
    else group_by(., { { x } }) } %>%
    mutate(Label = factor(paste0({ { x } }, " (n=", n, ")")),
           x = factor({ { x } })) %>%
    ungroup()

  labels <- df_prep %>%
    mutate(Levels = { { x } }) %>%
    { if (!quo_is_null(enquo(facet_col))) group_by(., Bucket)
    else . } %>%
    distinct(Levels, Label)

  return(list(df = df_prep, labels = labels))
}

show.n <- function (x, y_quantile = 0.01, prefix = "n = ") {
  return(data.frame(y = quantile(x, probs = y_quantile)[[1]], label = paste0(prefix, length(x))))
}

signif_boxplot <- function(df, x, y, facet_col = NULL,
                           test = wilcox.test, test.args = NULL,
                           signif_label = print_signif, title = NULL) {
  prep <- df %>%
    prep_signif({ { x } }, { { facet_col } })

  df_n <- prep$df %>%
    group_by({ { x } }, { { facet_col } }) %>%
    summarize(n = first(n),
              Q1 = quantile({ { y } }, probs = 0.25)[[1]],
              Q3 = quantile({ { y } }, probs = 0.75)[[1]],
              IQR = IQR({ { y } }),
              y_min = Q1 - IQR * 1.7,
              y_max = Q3 + IQR * 1.7,
              .groups = "drop") %>%
    mutate(y_min = min(y_min),
           y_max = max(y_max))

  plot <- prep$df %>%
    ggplot() +
    aes(x = { { x } }, y = { { y } }) +
    geom_boxplot(outliers = FALSE) +
    geom_signif(
      comparisons = list(levels(prep$df$x)),
      map_signif_level = signif_label,
      tip_length = 0, extend_line = -0.05, y_position = max(df_n$y_max),
      test = test, test.args = test.args
    ) +
    geom_text(data = df_n, aes(y = y_min, label = paste0("n = ", n))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  { if (!quo_is_null(enquo(facet_col))) facet_grid(~Bucket) } +
    xlab(as_name(enquo(x))) +
    # ToDo: Use name of test as subtitle
    ggtitle(title, subtitle = NULL)

  return(plot)
}

signif_violin_plot <- function(df, x, y, facet_col = NULL,
                               test = wilcox.test, test.args = NULL,
                               signif_label = print_signif, title = NULL,
                               count_y = NULL, count_size = base_size / 4) {

  prep <- df %>%
    prep_signif({ { x } }, { { facet_col } })

  plot <- prep$df %>%
    ggplot() +
    aes(x = { { x } }, y = { { y } }) +
    geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75),
                color = "#4080DB") +
    geom_signif(
      comparisons = list(levels(prep$df$x)),
      map_signif_level = signif_label,
      tip_length = 0, extend_line = -0.05,
      test = test, test.args = test.args
    ) +
    { if (is.null(count_y)) stat_summary(aes(y = min({ { y } }) - 0.2 * abs(max({ { y } }) - min({ { y } }))),
      fun.data = show.n, geom = "text", color = default_color, size = count_size)
      else stat_summary(aes(y = count_y), fun.data = show.n, geom = "text", color = default_color, size = count_size) } +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    {if (!quo_is_null(enquo(facet_col))) facet_grid(~Bucket)} +
    xlab(as_name(enquo(x))) +
    # ToDo: Use name of test as subtitle
    ggtitle(title, subtitle = NULL)

  return(plot)
}

signif_beeswarm_plot <- function(df, x, y, facet_col = NULL, color_col = NULL,
                                 test = wilcox.test, test.args = NULL, cex = 2,
                                 signif_label = print_signif, title = NULL,
                                 count_y = NULL, count_size = base_size / 4) {

  prep <- df %>%
    prep_signif({ { x } }, { { facet_col } })

  plot <- prep$df %>%
    ggplot() +
    aes(x = { { x } }, y = { { y } }) +
    geom_violin(trim = FALSE, draw_quantiles = 0.5,
                fill = "darkgrey", color = "darkgrey", alpha = 1/3) +
    {
      if (quo_is_null(enquo(color_col))) {
        geom_beeswarm(priority = "density", color = "#4080DB", cex = cex)
      } else {
        geom_beeswarm(aes(color = { { color_col } }), priority = "density", cex = cex)
      }
    } +
    geom_signif(
      comparisons = list(levels(prep$df$x)),
      map_signif_level = signif_label,
      tip_length = 0, extend_line = -0.05,
      test = test, test.args = test.args
    ) +
    { if (is.null(count_y)) stat_summary(aes(y = min({ { y } }) - 0.2 * abs(max({ { y } }) - min({ { y } }))),
      fun.data = show.n, geom = "text", color = default_color, size = count_size)
      else stat_summary(aes(y = count_y), fun.data = show.n, geom = "text", color = default_color, size = count_size) } +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_viridis_c(option = "D", direction = 1) +
    { if (!quo_is_null(enquo(facet_col))) facet_grid(~Bucket) } +
    xlab(as_name(enquo(x))) +
    # ToDo: Use name of test as subtitle
    ggtitle(title, subtitle = NULL)

  return(plot)
}

inverse_distance_weighting <- function(distances, values, p = 1) {
  if (any(distances == 0, na.rm = TRUE)) {
    return(values[distances == 0][1])
  }
  else {
    weights <- 1 / distances^p
    return(sum(weights * values, na.rm = TRUE) / sum(weights, na.rm = TRUE))
  }
}

# https://stackoverflow.com/questions/41087157/gradient-fill-violin-plots-using-ggplot2/78316904#78316904
gradient_violin_plot <- function(df, group_col, value_col, color_col, p = 2) {
  require(ggplot2)
  require(dplyr)

  df_input <- df %>%
    select({ { group_col } }, { { color_col } }, { { value_col } }) %>%
    rename(group = { { group_col } },
           color = { { color_col } },
           value = { { value_col } })

  scale_violin_width <- .45

  p <- ggplot(df_input, aes(x = group, y = value)) + geom_violin()

  breaks <- unique(as.integer(df_input$group))
  labels <- unique(df_input$group)

  vl_fill <- data.frame(ggplot_build(p)$data) %>%
    mutate(xnew = x - scale_violin_width * violinwidth, xend = x + scale_violin_width * violinwidth) %>%
    group_by(group) %>%
    mutate(color = sapply(y, \(y_item) inverse_distance_weighting(
      abs(df_input$value - y_item) * if_else(as.integer(df_input$group) == first(group), 1, Inf),
      df_input$color, p = 2))
    )

  ggplot() +
    geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y,
                                     color = color)) +
    geom_violin(data = df_input, aes(x = as.integer(group), y = value, fill = group),
                color = "grey", alpha = 0, draw_quantiles = c(0.25, 0.5, 0.75),
                show.legend = FALSE) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    scale_color_viridis_c() +
    labs(x = quo_name(enquo(group_col)),
         y = quo_name(enquo(value_col)),
         color = quo_name(enquo(color_col))) +
    theme(legend.position = "bottom")
}

# === SHAP-Plots ===

shap_plot <- function(df_explanation, alpha = 0.75, jitter_width = 0.25, show_legend = TRUE, title = NULL,
                      category_lab = "Feature", value_lab = "SHAP-Value", color_lab = "Feature Value") {
  max_abs_shap <- round(max(abs(df_explanation$SHAP.Value)), digits = 2)

  df_explanation %>%
    mutate(SHAP.Factor.Corr.Absolute = abs(SHAP.Factor.Corr)) %>%
    mutate(DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, SHAP.Factor.Corr.Absolute)) %>%
    arrange(Factor.Value.Relative) %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = SHAP.Value) +
    geom_hline(yintercept = 0) +
    # geom_boxplot(outlier.shape = NA, color = "black") +
    geom_quasirandom(aes(color = Factor.Value.Relative),
                     show.legend = show_legend, alpha = alpha, width = jitter_width) +
    scale_colour_viridis_c(option = "C", direction = 1, end = 0.95,
                           limits = c(0,1), breaks = c(0,1), labels = c("Low", "High")) +
    labs(title = title, x = category_lab, y = value_lab, color = color_lab) +
    coord_flip(ylim = c(-max_abs_shap, max_abs_shap))
}

shap_importance_plot <- function(df_explanation, bar_label_shift = 0.002, show_legend = TRUE, title = NULL,
                                 category_lab = "Feature", value_lab = "Median Absolute SHAP-Value",
                                 color_lab = "Feature Correlation") {
  df_explanation %>%
    select(-SHAP.Value) %>%
    distinct(DosageCompensation.Factor, .keep_all = TRUE) %>%
    mutate(DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, SHAP.Median.Absolute)) %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = SHAP.Median.Absolute,
        label = format(round(SHAP.Median.Absolute, 3), nsmall = 3)) +
    geom_hline(yintercept = 0) +
    geom_bar(aes(fill = SHAP.Factor.Corr), show.legend = show_legend, color = "black", stat = "identity") +
    geom_pointrange(aes(x = DosageCompensation.Factor, y = SHAP.Median.Absolute,
                        ymin = SHAP.p25.Absolute, ymax = SHAP.p75.Absolute),
                    colour = "orange", fatten = 1) +
    geom_text(color = "black", y = 0 + bar_label_shift, hjust = 0) +
    scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab", limits  = c(-1, 1)) +
    labs(title = title, x = category_lab, y = value_lab, fill = color_lab) +
    coord_flip(ylim = c(0, max(df_explanation$SHAP.p75.Absolute)))
}

shap_corr_importance_plot <- function(df_explanation, bar_label_shift = 0.02, show_legend = TRUE, title = NULL,
                                      category_lab = "Feature", value_lab = "Absolute SHAP-Value-Feature Correlation",
                                      color_lab = "Feature Correlation") {
  df_explanation %>%
    select(-SHAP.Value) %>%
    distinct(DosageCompensation.Factor, .keep_all = TRUE) %>%
    mutate(SHAP.Factor.Corr.Absolute = abs(SHAP.Factor.Corr)) %>%
    mutate(DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, SHAP.Factor.Corr.Absolute)) %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = SHAP.Factor.Corr.Absolute,
        label = print_corr(SHAP.Factor.Corr, SHAP.Factor.Corr.p.adj, signif = TRUE, map_p = TRUE)) +
    geom_hline(yintercept = 0) +
    geom_bar(aes(fill = SHAP.Factor.Corr), show.legend = show_legend, color = "black", stat = "identity") +
    geom_text(color = "black", y = 0 + bar_label_shift, hjust = 0) +
    scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab", limits = c(-1, 1)) +
    labs(title = title, x = category_lab, y = value_lab, fill = color_lab) +
    coord_flip(ylim = c(0, max(abs(df_explanation$SHAP.Factor.Corr))))
}

shap_plot_arrows <- function(df_shap, show_legend = TRUE, title = NULL,
                             category_lab = "Feature", value_lab = "SHAP-Value", color_lab = "Feature Value") {
  arrow_element <- arrow(length = unit(0.30, "cm"), ends = "last", type = "closed")
  max_abs_shap <- round(max(abs(df_shap$SHAP.Value)), digits = 2)

  plot_arrows <- df_shap %>%
    slice_sample(n = 1) %>%
    ggplot() +
    aes(x = SHAP.Value) +
    geom_segment(aes(x = -max_abs_shap * 0.05, y = 2, xend = -max_abs_shap, yend = 2),
                 arrow = arrow_element, color = color_palettes$BufferingClasses[["Buffered"]]) +
    geom_segment(aes(x = max_abs_shap * 0.05, y = 2, xend = max_abs_shap, yend = 2),
                 arrow = arrow_element, color = color_palettes$BufferingClasses[["Scaling"]]) +
    geom_text(aes(x = -max_abs_shap * 0.5, y = 1, label = "Buffering"),
              color = color_palettes$BufferingClasses[["Buffered"]]) +
    geom_text(aes(x = max_abs_shap * 0.5, y = 1, label = "Scaling"),
              color = color_palettes$BufferingClasses[["Scaling"]]) +
    ylim(c(0, 3)) +
    xlim(c(-max_abs_shap, max_abs_shap)) +
    cowplot::theme_nothing()

  plot_shap <- shap_plot(df_shap, show_legend = show_legend, title = title,
                         category_lab = category_lab, value_lab = value_lab, color_lab = color_lab)

  cowplot::plot_grid(plot_shap, plot_arrows,
                     labels = NULL, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.05))
}

nested_shap_heatmap <- function(df, nesting_formula) {
  df %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = "", fill = SHAP.Factor.Corr, label = Label) +
    geom_raster() +
    geom_text(color = "black") +
    ggh4x::facet_nested(nesting_formula, switch = "y") +
    scale_fill_gradientn(colors = c(bidirectional_color_pal_viridis[4], "white", bidirectional_color_pal_viridis[2]),
                         space = "Lab", limits = c(-1, 1), oob = scales::squish) +
    labs(x = "Feature", y = "Model", fill = "SHAP-Value Correlation") +
    cowplot::theme_minimal_grid() +
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_rect(color = "lightgrey"),
          axis.ticks.y = element_blank(),
          legend.key.size = unit(16, "points"),
          legend.key.width = unit(24, "points"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = "top",
          legend.direction = "horizontal",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(5, 15, 5, 15), "mm"))
}
