#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
})

# ggtree releases before the ggplot2 4.0 transition call the former internal
# helper name. Keep the shim local to this renderer until the installed ggtree
# release adopts ggplot2::is_waiver().
is.waive <- ggplot2:::is_waiver

ANALYSIS_LEVELS <- c(
  "true_tree_fixed_lengths",
  "true_topology_ml_lengths",
  "ml_tree_unadjusted",
  "ml_tree_bonferroni"
)
ANALYSIS_LABELS <- c(
  true_tree_fixed_lengths = "true tree +\ntrue branch\nlengths",
  true_topology_ml_lengths = "true tree +\nML-estimated\nbranch lengths",
  ml_tree_unadjusted = "ML-tree",
  ml_tree_bonferroni = "ML-tree &\nBonferroni\ncorrection"
)
ANALYSIS_COLORS <- c(
  true_tree_fixed_lengths = "#D95F8D",
  true_topology_ml_lengths = "#7B3294",
  ml_tree_unadjusted = "#2C7FB8",
  ml_tree_bonferroni = "#1A9850"
)
SITE_STYLES <- c(`100` = "22", `1000` = "solid", `10000` = "42")
SITE_ALPHA <- c(`100` = 0.48, `1000` = 1.0, `10000` = 0.42)
TREE_COLORS <- c(A = "#3B82F6", B = "#F97316", `0` = "#222222")

parse_args <- function(arguments) {
  result <- list()
  index <- 1L
  while (index <= length(arguments)) {
    key <- arguments[[index]]
    if (!startsWith(key, "--") || index == length(arguments)) {
      stop("Arguments must use --name value syntax; unexpected argument: ", key)
    }
    result[[substring(key, 3L)]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  result
}

balanced_newick <- function(labels) {
  if (length(labels) == 1L) return(labels)
  midpoint <- length(labels) %/% 2L
  paste0(
    "(",
    balanced_newick(labels[seq_len(midpoint)]),
    ",",
    balanced_newick(labels[(midpoint + 1L):length(labels)]),
    ")"
  )
}

make_tree <- function(internal) {
  if (internal) {
    newick <- paste0(
      "(",
      balanced_newick(paste0("A", seq_len(8L))),
      ",",
      balanced_newick(paste0("B", seq_len(8L))),
      ");"
    )
  } else {
    newick <- paste0("(", balanced_newick(paste0("A", seq_len(4L))), ",B1);")
  }
  tree <- read.tree(text = newick)
  groupOTU(
    tree,
    list(
      A = grep("^A", tree$tip.label, value = TRUE),
      B = grep("^B", tree$tip.label, value = TRUE)
    ),
    group_name = "clade"
  )
}

tree_inset <- function(internal) {
  tree <- make_tree(internal)
  plot <- ggtree(tree, layout = "unrooted", branch.length = "none", aes(color = clade), size = 1.05) +
    geom_tippoint(aes(color = clade), size = 1.5) +
    scale_color_manual(values = TREE_COLORS, na.value = "#222222", guide = "none") +
    theme_tree() +
    theme(
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
  plot_data <- plot$data
  tip_data <- plot_data[plot_data$isTip %in% TRUE, , drop = FALSE]
  a_tips <- tip_data[grepl("^A", tip_data$label), , drop = FALSE]
  b_tips <- tip_data[grepl("^B", tip_data$label), , drop = FALSE]
  center <- c(mean(tip_data$x), mean(tip_data$y))
  a_center <- c(mean(a_tips$x), mean(a_tips$y))
  b_center <- c(mean(b_tips$x), mean(b_tips$y))
  span <- max(diff(range(tip_data$x)), diff(range(tip_data$y)))
  outward <- function(point) {
    direction <- point - center
    length <- sqrt(sum(direction^2))
    point + if (length > 0) direction / length * span * 0.10 else c(0, 0)
  }
  a_label <- outward(a_center)
  b_label <- outward(b_center)
  question_label <- (a_center + b_center) / 2
  plot +
    annotate(
      "text",
      x = a_label[[1L]],
      y = a_label[[2L]],
      label = "A",
      color = TREE_COLORS[["A"]],
      fontface = "bold",
      size = 5.0
    ) +
    annotate(
      "text",
      x = b_label[[1L]],
      y = b_label[[2L]],
      label = "B",
      color = TREE_COLORS[["B"]],
      fontface = "bold",
      size = 5.0
    ) +
    annotate(
      "text",
      x = question_label[[1L]],
      y = question_label[[2L]],
      label = "?",
      color = "#111111",
      fontface = "bold",
      size = 4.8
    ) +
    coord_equal(clip = "off")
}

prepare_rows <- function(summary_path) {
  rows <- read.delim(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
  rows <- rows[
    rows$simulation_model == "JC" &
      rows$evaluation_model == "JC" &
      rows$formula == "dominant",
    ,
    drop = FALSE
  ]
  rows$analysis <- NA_character_
  rows$analysis[
    rows$scenario == "true_tree_fixed_lengths" & rows$decision_rule == "unadjusted"
  ] <- "true_tree_fixed_lengths"
  rows$analysis[
    rows$scenario == "true_topology_ml_lengths" & rows$decision_rule == "unadjusted"
  ] <- "true_topology_ml_lengths"
  rows$analysis[
    rows$scenario == "ml_tree" & rows$decision_rule == "unadjusted"
  ] <- "ml_tree_unadjusted"
  rows$analysis[
    rows$scenario == "ml_tree" & rows$decision_rule == "taxon_bonferroni"
  ] <- "ml_tree_bonferroni"
  rows <- rows[!is.na(rows$analysis), , drop = FALSE]
  rows$analysis <- factor(rows$analysis, levels = ANALYSIS_LEVELS)
  rows$nsites <- factor(as.integer(rows$nsites), levels = c(100L, 1000L, 10000L))
  rows$percent_informative <- 100 * as.numeric(rows$fraction_informative)
  rows$branch_length <- as.numeric(rows$branch_length)
  rows
}

power_panel <- function(rows, tree_case, panel_label, internal) {
  panel_rows <- rows[rows$tree_case == tree_case, , drop = FALSE]
  plot <- ggplot(
    panel_rows,
    aes(
      x = branch_length,
      y = percent_informative,
      color = analysis,
      linetype = nsites,
      alpha = nsites,
      group = interaction(analysis, nsites)
    )
  ) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 5, fill = "#BDBDBD", alpha = 0.72) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.3, stroke = 0) +
    annotate("label", x = -0.48, y = 103, label = panel_label, size = 9.5, fontface = "bold", fill = "#E5E5E5", linewidth = 0) +
    annotate("text", x = 0.62, y = 70, label = "100", color = "#8A8A8A", fontface = "bold", size = 3.4) +
    annotate("text", x = 1.40, y = 99, label = "1,000", color = "#111111", fontface = "bold", size = 3.4) +
    annotate("text", x = 3.20, y = 99, label = "10,000", color = "#AAAAAA", fontface = "bold", size = 3.4) +
    scale_color_manual(values = ANALYSIS_COLORS, labels = ANALYSIS_LABELS, name = NULL) +
    scale_linetype_manual(values = SITE_STYLES, guide = "none") +
    scale_alpha_manual(values = SITE_ALPHA, guide = "none") +
    scale_x_continuous(breaks = 0:10, limits = c(-0.65, 10.15), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 5, 25, 50, 75, 100), limits = c(-1, 105), expand = c(0, 0)) +
    labs(x = "Branch length", y = "Informative instances (%)") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1, linewidth = 0.9))) +
    theme_minimal(base_family = "sans", base_size = 9.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.35),
      axis.line = element_line(color = "#777777", linewidth = 0.35),
      axis.title = element_text(size = 10.5),
      legend.text = element_text(size = 8.2, lineheight = 0.9),
      legend.key.width = unit(1.8, "lines"),
      plot.margin = margin(8, 8, 4, 22)
    ) +
    coord_cartesian(clip = "off")
  plot + inset_element(
    tree_inset(internal),
    left = 0.53,
    bottom = 0.55,
    right = 0.96,
    top = 0.91,
    align_to = "panel",
    clip = FALSE
  )
}

main <- function() {
  file_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_argument) != 1L) stop("Could not resolve renderer path")
  script_path <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
  project <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
  arguments <- parse_args(commandArgs(trailingOnly = TRUE))
  summary_path <- if (is.null(arguments$summary)) {
    file.path(
      project,
      "artifacts/cluster/lisc/paper-fig2-jc-reps1000-v4-217bd244e77-20260715",
      "merged/head_to_head_summary.tsv"
    )
  } else {
    arguments$summary
  }
  output_path <- if (is.null(arguments$output)) {
    file.path(project, "manuscript/figures/figure_jc_baseline_power.pdf")
  } else {
    arguments$output
  }
  if (!file.exists(summary_path)) stop("Summary does not exist: ", summary_path)

  rows <- prepare_rows(summary_path)
  panel_a <- power_panel(rows, "five_external", "a", FALSE)
  panel_b <- power_panel(rows, "sixteen_internal", "b", TRUE)
  composition <- panel_a / guide_area() / panel_b +
    plot_layout(heights = c(1, 0.15, 1), guides = "collect") &
    theme(legend.position = "bottom")

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    output_path,
    composition,
    width = 7.0866,
    height = 11.6929,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  cat(output_path, "\n")
}

main()
