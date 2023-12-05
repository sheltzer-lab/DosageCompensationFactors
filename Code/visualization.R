library(ggplot2)
library(dplyr)
library(rlang)
library(EnhancedVolcano)
library(psych)
library(corrplot)
library(ggsignif)
library(viridisLite)
library(ggbeeswarm)
library(pheatmap)
library(forcats)

here::i_am("DosageCompensationFactors.Rproj")

vertical_bar_chart <- function(df, category_col, value_col,
                               color_col = NULL,
                               error_low_col = NULL, error_high_col = NULL,
                               value_range = c(0.45, 0.65), break_steps = 0.05,
                               line_intercept = 0.5, bar_label_shift = 0.002,
                               title = NULL, category_lab = NULL, value_lab = NULL, color_lab = NULL) {
  df %>%
    ggplot() +
    aes(x = { { category_col } }, y = { { value_col } },
        label = format(round({ { value_col } }, 3), nsmall = 3)) +
    { if (!quo_is_null(enquo(color_col))) geom_bar(aes(fill = {{ color_col }}), stat = "identity")
      else geom_bar(stat = "identity") } +
    geom_hline(yintercept = line_intercept) +
    { if (!quo_is_null(enquo(error_low_col)) & !quo_is_null(enquo(error_high_col)))
      geom_pointrange(aes(x = { { category_col } }, y = { { value_col } },
                          ymin = { { error_low_col } }, ymax = { { error_high_col } }),
                      colour = "orange", fatten = 1) } +
    geom_text(color = "white", y = value_range[1] + bar_label_shift, hjust = 0) +
    scale_y_continuous(breaks = seq(value_range[1], value_range[2], break_steps)) +
    scale_fill_viridis_c(option = "G", direction = -1, begin = 0.2, end = 0.8) +
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

plot_correlation <- function(df, method = "spearman") {
  cor_matrix <- psych::corr.test(df, method = method,
                                 adjust = "none")

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
    mutate(Buffering.Class = if_else(TTest.p.adjusted > p_threshold, "N.S.", Buffering.Class)) %>%
    mutate(Buffering.Class = factor({ { class_col } }, levels = buffering_levels)) %>%
    mutate(Label = if_else(Buffering.Class %in% c("Buffered", "Anti-Scaling"),
                           { { label_col } }, NA)) %>%
    plot_volcano({ { ratio_col } }, { { signif_col } }, Label, Buffering.Class,
                 color_mapping = color_mapping,
                 value_threshold = value_threshold, signif_threshold = signif_threshold,
                 title = title, subtitle = subtitle)
}

plot_volcano <- function(df, value_col, signif_col, label_col, color_col,
                         color_mapping = NULL,
                         value_threshold = log2fc_threshold, signif_threshold = p_threshold,
                         title = NULL, subtitle = NULL) {
  df %>%
    mutate(`-Log10(p)` = -log10({ { signif_col } })) %>%
    ggplot() +
    aes(x = { { value_col } }, y = `-Log10(p)`,
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
                     force = 2, max.overlaps = 20)
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

scatter_plot_regression <- function(df, x_col, y_col, formula, color_col = NULL,
                                    label_coords = c(0, 0), title = NULL, point_size = 0.3) {
  df <- df %>%
    select({ { x_col } }, { { y_col } }, { { color_col } }) %>%
    drop_na() %>%
    mutate(Density = get_density({ { x_col } }, { { y_col } }, n = 100))

  regression <- lm(formula, df)
  pred.int <- predict(regression, interval = "prediction")
  regression_summary <- summary(regression)
  df <- cbind(df, pred.int)
  slope <- regression$coefficients[[quo_name(enquo(x_col))]]
  intercept <- regression$coefficients[["(Intercept)"]]

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
    ggplot2::annotate("text", x = label_coords[1], y = label_coords[2], color = "blue",
                      label = paste("y =", format(round(slope, 5), nsmall = 5),
                                    "* x +", format(round(intercept, 5), nsmall = 5),
                                    ", R² = ", format(round(regression_summary$r.squared, 5), nsmall = 5)
                      )) +
    {
      if (quo_is_null(enquo(color_col))) scale_colour_viridis_c(option = "D", direction = 1)
    } +
    xlab(quo_name(enquo(x_col))) +
    ylab(quo_name(enquo(y_col))) +
    ggtitle(title)

  return(regression_plot)
}

waterfall_plot <- function(df, value_col, rank_col, label_col) {
  xlim <- c(0, max(df[[quo_name(enquo(rank_col))]]))
  label_nudge_x <- floor(xlim[2] / 5)

  ylim1 <- c(min(df[[quo_name(enquo(value_col))]]),
             (df %>% filter({ { rank_col } } == label_nudge_x))[[quo_name(enquo(value_col))]])
  ylim2 <- c((df %>% filter({ { rank_col } } == xlim[2] - label_nudge_x))[[quo_name(enquo(value_col))]],
             max(df[[quo_name(enquo(value_col))]]))

  df %>%
    ggplot() +
    aes(x = { { rank_col } }, y = { { value_col } }, label = { { label_col } }) +
    geom_hline(yintercept = 0, color = "red") +
    geom_point(size = 0.3) +
    geom_text_repel(data = df %>% slice_min({ { rank_col } }, n = 5),
                    xlim = xlim, ylim = ylim1, direction = "y", nudge_x = label_nudge_x,
                    seed = 42, color = "darkblue") +
    geom_text_repel(data = df %>% slice_max({ { rank_col } }, n = 5),
                    xlim = xlim, ylim = ylim2, direction = "y", nudge_x = -label_nudge_x,
                    seed = 42, color = "darkred")
}

save_plot <- function(plot, filename, dir = plots_dir,
                      height = 200, width = 200, dpi = 300) {
  ggsave(here(dir, filename), plot = plot,
       height = height, width = width, units = "mm", dpi = dpi)
  return(plot)
}

plot_rocs <- function(df_rocs) {
  df_label <- df_rocs %>%
    distinct(Name, AUC) %>%
    arrange(desc(Name)) %>%
    mutate(x = 0.1,
           y = seq(0.05, 0.75, 0.05)[seq_along(unique(df_rocs$Name))],
           AUC = paste0("AUC = ", format(round(AUC, 3), nsmall = 3)))

  plot <- df_rocs %>%
    arrange(Sensitivity) %>%
    ggplot() +
    aes(x = Specificity, y = Sensitivity, color = Name) +
    geom_abline(slope = 1, intercept = 1, color = "grey") +
    geom_line() +
    geom_label(data = df_label, mapping = aes(color = Name, label = AUC, x = x, y = y)) +
    scale_x_reverse(limits = c(1, 0)) +
    labs(x = "Specificity", y = "Sensitivity", color = "Model")

  return(plot)
}

print_signif <- function(p, digits = 4) {
  paste0("p ", if_else(p < 10^(-digits),
                      paste0("< ", format(10^(-digits), nsmall = digits, scientific = FALSE)),
                      paste0("= ", format(round(p, digits), nsmall = digits, scientific = FALSE))))
}

map_signif <- function (p) {
  case_when(
      p < 0.0001 ~ "***",
      p < 0.001 ~ "**",
      p < 0.01 ~ "*",
      TRUE ~ "N.S."
    )
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

plot_pca <- function(pca) {
  eigenvalues <- pca$eigenvalues
  df_pca <- pca$df_pca

  df_pca %>%
    ggplot() +
    aes(.fittedPC1, .fittedPC2, color = CellLine.Name, label = Sample.Name) +
    geom_point(size = 1.5) +
    geom_label_repel(force = 10, seed = 123) +
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

bidirectional_heatmap <- function(df, value_col, sample_col, group_col,
                                  cluster_rows = FALSE, cluster_cols = FALSE,
                                  show_rownames = TRUE, show_colnames = TRUE,
                                  transpose = FALSE, palette_length = 100, color_pal = bidirectional_color_pal) {
  mat <- df %>%
    group_by({ { sample_col } }, { { group_col } }) %>%
    summarize(Value = mean({ { value_col } }, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = { { sample_col } }, values_from = Value, id_cols = { { group_col } }) %>%
    arrange({ { group_col } }) %>%
    tibble::column_to_rownames(var = quo_name(enquo(group_col)))

  if (transpose) {
    mat <- t(mat)
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
           breaks = breaks, color = color)
}


bucketed_scatter_plot <- function(df, value_col, x_value_col, bucket_col,
                                  threshold_low = NULL, threshold_high = NULL,
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
    scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 1)) +
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

signif_violin_plot <- function(df, x, y, facet_col = NULL,
                               test = wilcox.test, test.args = NULL,
                               signif_label = print_signif, title = NULL) {

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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = prep$labels$Label, breaks = prep$labels$Levels) +
    {if (!quo_is_null(enquo(facet_col))) facet_grid(~Bucket)} +
    xlab(as_name(enquo(x))) +
    # ToDo: Use name of test as subtitle
    ggtitle(title, subtitle = NULL)

  return(plot)
}

signif_beeswarm_plot <- function(df, x, y, facet_col = NULL, color_col = NULL,
                                 test = wilcox.test, test.args = NULL, cex = 2,
                                 signif_label = print_signif, title = NULL) {

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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = prep$labels$Label, breaks = prep$labels$Levels) +
    scale_colour_viridis_c(option = "D", direction = 1) +
    { if (!quo_is_null(enquo(facet_col))) facet_grid(~Bucket) } +
    xlab(as_name(enquo(x))) +
    # ToDo: Use name of test as subtitle
    ggtitle(title, subtitle = NULL)

  return(plot)
}


shap_plot <- function(df_explanation, alpha = 0.75, jitter_width = 0.15) {
  max_abs_shap <- round(max(abs(df_explanation$SHAP.Value)), digits = 2)

  df_explanation %>%
    mutate(DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, SHAP.Median.Absolute)) %>%
    arrange(Factor.Value.Relative) %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = SHAP.Value) +
    geom_hline(yintercept = 0) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_quasirandom(aes(color = Factor.Value.Relative), alpha = alpha, width = jitter_width) +
    coord_flip(ylim = c(-max_abs_shap, max_abs_shap)) +
    scale_colour_viridis_c(option = "C", direction = 1, end = 0.95)
}

shap_importance_plot <- function(df_explanation,
                                 bar_label_shift = 0.002,
                                 title = NULL, category_lab = "Feature", value_lab = "Median Absolute SHAP-Value",
                                 color_lab = "Feature Correlation") {
  df_explanation %>%
    select(-SHAP.Value) %>%
    distinct(DosageCompensation.Factor, .keep_all = TRUE) %>%
    mutate(DosageCompensation.Factor = fct_reorder(DosageCompensation.Factor, SHAP.Median.Absolute)) %>%
    ggplot() +
    aes(x = DosageCompensation.Factor, y = SHAP.Median.Absolute,
        label = format(round(SHAP.Median.Absolute, 3), nsmall = 3)) +
    geom_hline(yintercept = 0) +
    geom_bar(aes(fill = SHAP.Factor.Corr), color = "black", stat = "identity") +
    geom_pointrange(aes(x = DosageCompensation.Factor, y = SHAP.Median.Absolute,
                        ymin = SHAP.p25.Absolute, ymax = SHAP.p75.Absolute),
                    colour = "orange", fatten = 1) +
    geom_text(color = "black", y = 0 + bar_label_shift, hjust = 0) +
    scale_fill_gradientn(colors = bidirectional_color_pal, space = "Lab", limits  = c(-1, 1)) +
    labs(title = title, x = category_lab, y = value_lab, fill = color_lab) +
    coord_flip(ylim = c(0, max(df_explanation$SHAP.p75.Absolute)))
}
