#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(plot3D)
  library(plotly)
})

FORMULAS <- c("dominant", "eigenvalue_weighted")
FORMULA_LABELS <- c(
  dominant = "Published dominant coefficient",
  eigenvalue_weighted = "Eigenvalue-weighted coefficient"
)
RULES <- c("unadjusted", "taxon_bonferroni", "by_fdr")
RULE_LABELS <- c(
  unadjusted = "Unadjusted",
  taxon_bonferroni = "Taxon Bonferroni",
  by_fdr = "Benjamini–Yekutieli"
)
TREE_LABELS <- c(
  five_external = "Five-taxon external branch",
  sixteen_internal = "16-taxon internal branch"
)
MODEL_LABELS <- c(JC = "JC", GTR_PF06346 = "GTR PF06346")
SITE_COLORS <- c("100" = "#0072B2", "1000" = "#E69F00")
SITE_SHAPES <- c("100" = 21, "1000" = 24)
SITE_LINES <- c("100" = "solid", "1000" = "22")

parse_args <- function(values) {
  result <- list()
  index <- 1L
  while (index <= length(values)) {
    key <- values[[index]]
    if (!startsWith(key, "--") || index == length(values)) {
      stop("Arguments must use --name value syntax; unexpected argument: ", key)
    }
    result[[substring(key, 3L)]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required argument --", name)
  value
}

safe_name <- function(value) {
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  gsub("^_|_$", "", value)
}

read_strict <- function(path, required, label) {
  if (!file.exists(path)) stop(label, " file does not exist: ", path)
  rows <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(rows))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  if (!nrow(rows)) stop(label, " contains no data rows")
  rows
}

wilson_interval <- function(successes, total, z = qnorm(0.975)) {
  if (total <= 0) return(c(NA_real_, NA_real_))
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half <- z * sqrt(proportion * (1 - proportion) / total + z^2 / (4 * total^2)) / denominator
  c(max(0, center - half), min(1, center + half))
}

bootstrap_mean_interval <- function(values, replicates) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (!length(values)) return(c(NA_real_, NA_real_, NA_real_))
  estimate <- mean(values)
  if (length(values) == 1L || all(values == values[[1L]])) {
    return(c(estimate, estimate, estimate))
  }
  frequencies <- table(values)
  support <- as.numeric(names(frequencies))
  draws <- rmultinom(replicates, length(values), prob = as.numeric(frequencies))
  bootstrap_means <- as.numeric(crossprod(support, draws)) / length(values)
  limits <- quantile(bootstrap_means, c(0.025, 0.975), names = FALSE, type = 6)
  c(estimate, limits[[1L]], limits[[2L]])
}

validate_design_grid <- function(rows, expected_replicates, expected_samples, label) {
  expected_models <- c("JC", "GTR_PF06346")
  expected_trees <- c("five_external", "sixteen_internal")
  expected_sites <- c(100L, 1000L)
  expected_branches <- c(4, 8)
  checks <- list(
    formula = c(FORMULAS), decision_rule = c(RULES), simulation_model = expected_models,
    evaluation_model = expected_models, tree_case = expected_trees,
    nsites = expected_sites, branch_length = expected_branches, sample = expected_samples
  )
  for (column in names(checks)) {
    observed <- sort(unique(rows[[column]]))
    expected <- sort(checks[[column]])
    if (!identical(as.character(observed), as.character(expected))) {
      stop(label, " has unexpected ", column, ": ", paste(observed, collapse = ", "))
    }
  }
  if (any(rows$simulation_model != rows$evaluation_model)) {
    stop(label, " must contain matched simulation/evaluation models only")
  }
  key_columns <- c(
    "design", "sample", "tree_case", "simulation_model", "evaluation_model",
    "nsites", "branch_length", "formula", "decision_rule", "replicate"
  )
  if (anyDuplicated(rows[key_columns])) stop(label, " contains duplicate replicate rows")
  cell_columns <- setdiff(key_columns, "replicate")
  expected_cells <- expand.grid(
    sample = expected_samples,
    tree_case = expected_trees,
    simulation_model = expected_models,
    nsites = expected_sites,
    branch_length = expected_branches,
    formula = FORMULAS,
    decision_rule = RULES,
    stringsAsFactors = FALSE
  )
  expected_cells$design <- unique(rows$design)
  expected_cells$evaluation_model <- expected_cells$simulation_model
  expected_cells <- expected_cells[cell_columns]
  observed_cells <- unique(rows[cell_columns])
  cell_key <- function(values) do.call(paste, c(values, sep = "\r"))
  if (!setequal(cell_key(observed_cells), cell_key(expected_cells))) {
    stop(label, " does not contain the complete prespecified cell grid")
  }
  for (sample_name in expected_samples) {
    sample_rows <- rows[rows$sample == sample_name, , drop = FALSE]
    sample_groups <- split(
      sample_rows,
      interaction(sample_rows[cell_columns], drop = TRUE, lex.order = TRUE)
    )
    expected_ids <- expected_replicates[[sample_name]]
    expected_count <- length(expected_ids)
    sample_counts <- vapply(sample_groups, nrow, integer(1))
    if (!length(sample_counts) || any(sample_counts != expected_count)) {
      stop(label, " has incomplete ", sample_name, " cells; expected ", expected_count, " replicates")
    }
    complete_ids <- vapply(
      sample_groups,
      function(group) identical(sort(as.integer(group$replicate)), expected_ids),
      logical(1)
    )
    if (!all(complete_ids)) {
      stop(label, " has incorrect replicate identifiers in at least one ", sample_name, " cell")
    }
  }
  invisible(TRUE)
}

summarise_metrics <- function(rows, bootstrap_reps) {
  key_columns <- c(
    "design", "sample", "tree_case", "simulation_model", "evaluation_model",
    "nsites", "branch_length", "formula", "decision_rule"
  )
  group_id <- interaction(rows[key_columns], drop = TRUE, lex.order = TRUE)
  summaries <- lapply(split(rows, group_id), function(group) {
    metric_values <- list(
      empirical_fdr = as.numeric(group$fdp),
      power = as.numeric(group$power),
      null_rejection_rate = ifelse(
        as.integer(group$null_branches) > 0L,
        as.integer(group$false_discoveries) / as.integer(group$null_branches),
        NA_real_
      )
    )
    do.call(rbind, lapply(names(metric_values), function(metric) {
      values <- metric_values[[metric]]
      valid <- values[is.finite(values)]
      if (metric == "null_rejection_rate" && length(valid)) {
        successes <- sum(valid)
        interval <- wilson_interval(successes, length(valid))
        estimate <- successes / length(valid)
      } else {
        interval3 <- bootstrap_mean_interval(valid, bootstrap_reps)
        estimate <- interval3[[1L]]
        interval <- interval3[2:3]
      }
      data.frame(
        group[1L, key_columns, drop = FALSE], metric = metric,
        replicates = length(valid), estimate = estimate,
        ci_low = interval[[1L]], ci_high = interval[[2L]],
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, summaries)
}

summarise_paired_differences <- function(rows, bootstrap_reps) {
  pair_columns <- c(
    "design", "sample", "tree_case", "simulation_model", "evaluation_model",
    "nsites", "branch_length", "decision_rule", "replicate", "seed"
  )
  dominant <- rows[rows$formula == "dominant", , drop = FALSE]
  weighted <- rows[rows$formula == "eigenvalue_weighted", , drop = FALSE]
  if (anyDuplicated(dominant[pair_columns]) || anyDuplicated(weighted[pair_columns])) {
    stop("Duplicate formula rows prevent one-to-one calibration pairing")
  }
  paired <- merge(
    dominant, weighted, by = pair_columns, all = FALSE, sort = FALSE,
    suffixes = c("_dominant", "_weighted")
  )
  if (nrow(paired) != nrow(dominant) || nrow(paired) != nrow(weighted)) {
    stop("Calibration formula rows are not complete one-to-one pairs")
  }
  group_columns <- setdiff(pair_columns, c("replicate", "seed"))
  group_id <- interaction(paired[group_columns], drop = TRUE, lex.order = TRUE)
  summaries <- lapply(split(paired, group_id), function(group) {
    deltas <- list(
      empirical_fdr = as.numeric(group$fdp_weighted) - as.numeric(group$fdp_dominant),
      power = as.numeric(group$power_weighted) - as.numeric(group$power_dominant),
      null_rejection_rate = ifelse(
        as.integer(group$null_branches_dominant) > 0L &
          as.integer(group$null_branches_weighted) > 0L,
        as.integer(group$false_discoveries_weighted) /
          as.integer(group$null_branches_weighted) -
          as.integer(group$false_discoveries_dominant) /
          as.integer(group$null_branches_dominant),
        NA_real_
      )
    )
    do.call(rbind, lapply(names(deltas), function(metric) {
      values <- deltas[[metric]]
      values <- values[is.finite(values)]
      interval <- bootstrap_mean_interval(values, bootstrap_reps)
      discordant <- values[values != 0]
      sign_p <- if (!length(discordant)) 1 else {
        binom.test(sum(discordant > 0), length(discordant), p = 0.5)$p.value
      }
      data.frame(
        group[1L, group_columns, drop = FALSE], metric = metric,
        pairs = length(values), estimate = interval[[1L]],
        ci_low = interval[[2L]], ci_high = interval[[3L]],
        weighted_higher = sum(values > 0), dominant_higher = sum(values < 0),
        equal = sum(values == 0), paired_sign_p = sign_p,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, summaries)
}

decorate <- function(rows) {
  rows$model_label <- unname(MODEL_LABELS[rows$simulation_model])
  rows$tree_label <- unname(TREE_LABELS[rows$tree_case])
  rows$rule_label <- unname(RULE_LABELS[rows$decision_rule])
  rows$site_label <- paste(
    format(rows$nsites, big.mark = ",", scientific = FALSE, trim = TRUE),
    "sites"
  )
  if ("metric" %in% names(rows)) {
    rows$metric_label <- unname(c(
      empirical_fdr = "Empirical FDR", power = "Power",
      null_rejection_rate = "Null rejection rate"
    )[rows$metric])
  }
  rows$tree_label <- factor(rows$tree_label, levels = unname(TREE_LABELS))
  rows$rule_label <- factor(rows$rule_label, levels = unname(RULE_LABELS))
  rows$site_label <- factor(rows$site_label, levels = c("100 sites", "1,000 sites"))
  if ("metric_label" %in% names(rows)) {
    rows$metric_label <- factor(
      rows$metric_label,
      levels = c("Empirical FDR", "Power", "Null rejection rate")
    )
  }
  rows
}

common_theme <- theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#444444", margin = margin(b = 8)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 9.5),
    legend.position = "bottom",
    plot.caption = element_text(color = "#555555", hjust = 0)
  )

save_static <- function(plot, base_path, width, height) {
  ggsave(paste0(base_path, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggsave(
    paste0(base_path, ".svg"), plot, width = width, height = height,
    device = svglite::svglite
  )
  ggsave(
    paste0(base_path, ".png"), plot, width = width, height = height,
    dpi = 450, device = ragg::agg_png, background = "white"
  )
}

operating_plot <- function(rows, title, difference = FALSE) {
  rows <- rows[rows$metric %in% c("empirical_fdr", "power"), , drop = FALSE]
  rows$panel_label <- factor(
    paste(as.character(rows$tree_label), as.character(rows$metric_label), sep = "\n"),
    levels = as.vector(t(outer(
      unname(TREE_LABELS), c("Empirical FDR", "Power"),
      function(tree, metric) paste(tree, metric, sep = "\n")
    )))
  )
  reference <- unique(rows[
    rows$metric == "empirical_fdr", c("panel_label", "rule_label"), drop = FALSE
  ])
  limit <- if (difference) {
    max(0.05, max(abs(c(rows$ci_low, rows$ci_high)), na.rm = TRUE) * 1.08)
  } else {
    1
  }
  plot <- ggplot(
    rows,
    aes(
      x = branch_length, y = estimate, color = site_label,
      linetype = site_label, shape = site_label, group = site_label
    )
  ) +
    geom_hline(
      data = reference, aes(yintercept = if (difference) 0 else 0.05),
      inherit.aes = FALSE, color = "#555555", linewidth = 0.45
    ) +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), linewidth = 0.5, alpha = 0.55) +
    geom_line(linewidth = 0.9) +
    geom_point(fill = "white", size = 2.4, stroke = 0.8) +
    facet_grid(rows = vars(panel_label), cols = vars(rule_label), scales = "free_y") +
    scale_x_continuous(breaks = c(4, 8)) +
    scale_color_manual(values = c("100 sites" = SITE_COLORS[["100"]], "1,000 sites" = SITE_COLORS[["1000"]])) +
    scale_linetype_manual(values = c("100 sites" = SITE_LINES[["100"]], "1,000 sites" = SITE_LINES[["1000"]])) +
    scale_shape_manual(values = c("100 sites" = SITE_SHAPES[["100"]], "1,000 sites" = SITE_SHAPES[["1000"]])) +
    labs(
      title = title,
      subtitle = if (difference) {
        "Eigenvalue-weighted minus dominant; 95% paired bootstrap intervals"
      } else {
        "Prespecified 50:50 target-null/target-alternative mixture; 95% intervals"
      },
      x = "Target branch length", y = if (difference) "Paired difference" else "Probability",
      color = "Alignment length", linetype = "Alignment length", shape = "Alignment length",
      caption = if (difference) {
        "The horizontal reference is zero. Lines connect only the two prespecified branch lengths."
      } else {
        "The empirical-FDR panels include the nominal 0.05 reference. Lines connect only the two prespecified branch lengths."
      }
    ) + common_theme +
    theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 9))
  if (difference) {
    plot <- plot + scale_y_continuous(limits = c(-limit, limit), labels = scales::label_percent(accuracy = 1))
  } else {
    plot <- plot + scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1))
  }
  plot
}

null_calibration_plot <- function(rows, title) {
  ggplot(
    rows,
    aes(
      x = branch_length, y = estimate, color = calibration_type,
      linetype = site_label, shape = site_label,
      group = interaction(calibration_type, site_label)
    )
  ) +
    geom_hline(yintercept = 0.05, color = "#555555", linewidth = 0.45) +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), linewidth = 0.5, alpha = 0.55) +
    geom_line(linewidth = 0.9) +
    geom_point(fill = "white", size = 2.4, stroke = 0.8) +
    facet_grid(rows = vars(tree_label), cols = vars(rule_label)) +
    scale_x_continuous(breaks = c(4, 8)) +
    scale_y_continuous(
      limits = c(0, max(0.12, rows$ci_high, na.rm = TRUE) * 1.05),
      labels = scales::label_percent(accuracy = 1)
    ) +
    scale_color_manual(values = c("Native α = 0.05" = "#0072B2", "Training-derived threshold" = "#E69F00")) +
    scale_linetype_manual(values = c("100 sites" = SITE_LINES[["100"]], "1,000 sites" = SITE_LINES[["1000"]])) +
    scale_shape_manual(values = c("100 sites" = SITE_SHAPES[["100"]], "1,000 sites" = SITE_SHAPES[["1000"]])) +
    labs(
      title = title,
      subtitle = "Exact-null held-out target rejection with 95% Wilson intervals",
      x = "Target branch length", y = "Held-out null rejection rate",
      color = "Threshold", linetype = "Alignment length", shape = "Alignment length",
      caption = "The training-derived threshold is estimated from 2,500 independent training replicates and evaluated on 2,500 held-out replicates."
    ) + common_theme
}

static_3d <- function(rows, title, pdf_path, difference = FALSE) {
  tree_values <- levels(rows$tree_label)
  rule_values <- levels(rows$rule_label)
  site_values <- sort(unique(rows$nsites))
  panel_ids <- matrix(
    seq_len(length(tree_values) * length(rule_values)),
    nrow = length(tree_values), ncol = length(rule_values), byrow = TRUE
  )
  legend_id <- max(panel_ids) + 1L
  layout_matrix <- rbind(panel_ids, rep(legend_id, ncol(panel_ids)))
  if (difference) {
    z_limit <- max(0.05, max(abs(c(rows$ci_low, rows$ci_high)), na.rm = TRUE) * 1.08)
    z_range <- c(-z_limit, z_limit)
  } else if (unique(rows$metric) == "empirical_fdr") {
    z_range <- c(0, max(0.1, max(rows$ci_high, na.rm = TRUE) * 1.08))
  } else {
    z_range <- c(0, 1)
  }
  grDevices::cairo_pdf(pdf_path, width = 15.5, height = 8.8, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(layout_matrix, heights = c(rep(1, length(tree_values)), 0.22))
  par(
    oma = c(2.0, 0.5, 4.4, 0.5), mar = c(1.0, 1.0, 2.8, 0.4),
    family = "sans", fg = "#333333", col.axis = "#444444", col.lab = "#333333"
  )
  for (tree_value in tree_values) {
    for (rule_value in rule_values) {
      selected <- rows[
        as.character(rows$tree_label) == tree_value &
          as.character(rows$rule_label) == rule_value,
        , drop = FALSE
      ]
      first <- TRUE
      for (site_value in site_values) {
        trace <- selected[selected$nsites == site_value, , drop = FALSE]
        trace <- trace[order(trace$branch_length), , drop = FALSE]
        if (!nrow(trace)) next
        plot3D::scatter3D(
          x = trace$branch_length, y = log10(trace$nsites), z = trace$estimate,
          type = "o", add = !first, colvar = NULL,
          col = SITE_COLORS[[as.character(site_value)]], colkey = FALSE,
          lty = SITE_LINES[[as.character(site_value)]], lwd = 2.1,
          pch = SITE_SHAPES[[as.character(site_value)]], cex = 0.78, bg = "white",
          xlim = c(4, 8), ylim = c(2, 3), zlim = z_range,
          xlab = "Branch length", ylab = "log10 sites",
          zlab = if (difference) "Weighted − dominant" else as.character(unique(trace$metric_label)),
          theta = 42, phi = 23, d = 2.35, expand = 0.8,
          bty = "u", ticktype = "detailed", nticks = 4,
          col.panel = "#FAFAFA", col.grid = "#DDDDDD", col.axis = "#555555",
          lwd.panel = 0.65, lwd.grid = 0.55,
          main = if (first) paste(tree_value, rule_value, sep = "\n") else NULL,
          cex.main = 0.86
        )
        first <- FALSE
      }
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center", legend = paste(format(site_values, big.mark = ","), "sites"),
    col = unname(SITE_COLORS[as.character(site_values)]),
    lty = unname(SITE_LINES[as.character(site_values)]),
    pch = unname(SITE_SHAPES[as.character(site_values)]), pt.bg = "white",
    lwd = 2.1, pt.cex = 0.9, horiz = TRUE, bty = "n", cex = 0.9
  )
  mtext(title, side = 3, outer = TRUE, line = 2.65, font = 2, cex = 1.2)
  mtext(
    "Static 3D view of measured alignment-length trajectories",
    side = 3, outer = TRUE, line = 1.2, cex = 0.9, col = "#555555"
  )
  mtext(
    "Only the two prespecified branch lengths and two alignment lengths are connected; no surface interpolation.",
    side = 1, outer = TRUE, line = 0.35, cex = 0.8, col = "#555555"
  )
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

interactive_3d <- function(rows, title, html_path, difference = FALSE) {
  selections <- unique(rows[c("tree_label", "rule_label")])
  selections <- selections[order(selections$tree_label, selections$rule_label), , drop = FALSE]
  site_values <- sort(unique(rows$nsites))
  if (difference) {
    z_limit <- max(0.05, max(abs(c(rows$ci_low, rows$ci_high)), na.rm = TRUE) * 1.08)
    z_range <- c(-z_limit, z_limit)
  } else if (unique(rows$metric) == "empirical_fdr") {
    z_range <- c(0, max(0.1, max(rows$ci_high, na.rm = TRUE) * 1.08))
  } else {
    z_range <- c(0, 1)
  }
  figure <- plot_ly(width = 1280, height = 860)
  for (selection_index in seq_len(nrow(selections))) {
    selection <- selections[selection_index, , drop = FALSE]
    selected <- rows[
      as.character(rows$tree_label) == as.character(selection$tree_label) &
        as.character(rows$rule_label) == as.character(selection$rule_label),
      , drop = FALSE
    ]
    for (site_value in site_values) {
      trace <- selected[selected$nsites == site_value, , drop = FALSE]
      trace <- trace[order(trace$branch_length), , drop = FALSE]
      hover <- paste0(
        "<b>", if (difference) "Weighted − dominant" else as.character(trace$metric_label), "</b>",
        "<br>Branch length: ", trace$branch_length,
        "<br>Alignment length: ", format(trace$nsites, big.mark = ","),
        "<br>Estimate: ", sprintf("%+.4f", trace$estimate),
        "<br>95% CI: [", sprintf("%+.4f", trace$ci_low), ", ", sprintf("%+.4f", trace$ci_high), "]",
        "<br>Replicates/pairs: ", if ("pairs" %in% names(trace)) trace$pairs else trace$replicates
      )
      figure <- add_trace(
        figure, x = trace$branch_length, y = trace$nsites, z = trace$estimate,
        type = "scatter3d", mode = "lines+markers",
        name = paste(format(site_value, big.mark = ","), "sites"),
        legendgroup = as.character(site_value), showlegend = selection_index == 1L,
        visible = selection_index == 1L,
        line = list(color = SITE_COLORS[[as.character(site_value)]], width = 6),
        marker = list(
          color = SITE_COLORS[[as.character(site_value)]], size = 5,
          symbol = "circle", line = list(color = "white", width = 0.8)
        ),
        text = hover, hoverinfo = "text"
      )
    }
  }
  buttons <- lapply(seq_len(nrow(selections)), function(selection_index) {
    visible <- rep(FALSE, nrow(selections) * length(site_values))
    first <- (selection_index - 1L) * length(site_values) + 1L
    visible[first:(first + length(site_values) - 1L)] <- TRUE
    selection <- selections[selection_index, , drop = FALSE]
    selected_title <- paste0(
      title, "<br><sup>", selection$tree_label, " · ", selection$rule_label, "</sup>"
    )
    list(
      method = "update", args = list(list(visible = visible), list(title = list(text = selected_title))),
      label = paste(selection$tree_label, selection$rule_label, sep = " · ")
    )
  })
  initial <- selections[1L, , drop = FALSE]
  figure <- layout(
    figure,
    title = list(text = paste0(title, "<br><sup>", initial$tree_label, " · ", initial$rule_label, "</sup>"), x = 0.5),
    scene = list(
      xaxis = list(title = "Target branch length", showspikes = FALSE),
      yaxis = list(title = "Alignment length (sites)", type = "log", showspikes = FALSE),
      zaxis = list(
        title = if (difference) "Weighted − dominant" else as.character(unique(rows$metric_label)),
        range = z_range, zeroline = TRUE, zerolinecolor = "#444444", showspikes = FALSE
      ),
      aspectmode = "manual", aspectratio = list(x = 1.4, y = 1.0, z = 0.9),
      camera = list(eye = list(x = 1.55, y = -1.65, z = 1.15))
    ),
    legend = list(orientation = "h", x = 0, y = -0.12),
    margin = list(l = 20, r = 20, b = 100, t = 110),
    updatemenus = list(list(
      type = "dropdown", direction = "down", x = 0.01, y = 1.08,
      xanchor = "left", yanchor = "top", buttons = buttons
    )),
    annotations = list(list(
      text = "Measured grid points only; no surface interpolation.",
      x = 0.5, y = -0.17, xref = "paper", yref = "paper", showarrow = FALSE,
      font = list(size = 12, color = "#555555")
    ))
  )
  figure <- config(figure, responsive = TRUE, displaylogo = FALSE, scrollZoom = TRUE)
  htmlwidgets::saveWidget(figure, html_path, selfcontained = TRUE, title = title)
  dependency_dir <- sub("\\.html$", "_files", html_path)
  if (dir.exists(dependency_dir)) {
    html_text <- paste(readLines(html_path, warn = FALSE), collapse = "")
    if (grepl(basename(dependency_dir), html_text, fixed = TRUE)) {
      stop("Interactive 3D output is not self-contained: ", html_path)
    }
    unlink(dependency_dir, recursive = TRUE)
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
exact_path <- required_arg(args, "exact-replicate")
mixed_path <- required_arg(args, "mixed-replicate")
heldout_path <- required_arg(args, "heldout")
outdir <- required_arg(args, "outdir")
expected_exact <- if (is.null(args[["expected-exact"]])) 5000L else as.integer(args[["expected-exact"]])
training_reps <- if (is.null(args[["training-reps"]])) 2500L else as.integer(args[["training-reps"]])
expected_mixed <- if (is.null(args[["expected-mixed"]])) 2000L else as.integer(args[["expected-mixed"]])
bootstrap_reps <- if (is.null(args[["bootstrap-reps"]])) 5000L else as.integer(args[["bootstrap-reps"]])
if (bootstrap_reps < 1000L) stop("--bootstrap-reps must be at least 1000")
if (training_reps <= 0L || training_reps >= expected_exact) {
  stop("--training-reps must split the exact-null replicates into non-empty samples")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

replicate_required <- c(
  "schema_version", "design", "sample", "tree_case", "simulation_model",
  "evaluation_model", "nsites", "branch_length", "replicate", "seed", "formula",
  "decision_rule", "family_size", "null_branches", "alternative_branches",
  "discoveries", "false_discoveries", "true_discoveries", "fdp", "power"
)
exact <- read_strict(exact_path, replicate_required, "Exact-null replicate table")
mixed <- read_strict(mixed_path, replicate_required, "Mixed-null replicate table")
if (any(exact$schema_version != 2L) || any(mixed$schema_version != 2L)) {
  stop("Calibration renderer requires schema version 2")
}
if (!all(exact$design == "exact_null") || !all(mixed$design == "mixed_null")) {
  stop("Calibration design labels do not match their input files")
}
for (column in c("nsites", "replicate", "family_size", "null_branches", "alternative_branches", "discoveries", "false_discoveries", "true_discoveries")) {
  exact[[column]] <- as.integer(exact[[column]])
  mixed[[column]] <- as.integer(mixed[[column]])
}
for (column in c("branch_length", "fdp", "power")) {
  exact[[column]] <- as.numeric(exact[[column]])
  mixed[[column]] <- as.numeric(mixed[[column]])
}
validate_design_grid(
  exact,
  list(
    training = seq_len(training_reps),
    holdout = seq.int(training_reps + 1L, expected_exact)
  ),
  c("training", "holdout"), "Exact-null replicate table"
)
validate_design_grid(
  mixed, list(operating = seq_len(expected_mixed)), c("operating"),
  "Mixed-null replicate table"
)

set.seed(20260715)
exact_summary <- summarise_metrics(exact, bootstrap_reps)
mixed_summary <- summarise_metrics(mixed, bootstrap_reps)
exact_paired <- summarise_paired_differences(exact, bootstrap_reps)
mixed_paired <- summarise_paired_differences(mixed, bootstrap_reps)
write.table(
  rbind(exact_summary, mixed_summary), file.path(outdir, "fdr_operating_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  rbind(exact_paired, mixed_paired), file.path(outdir, "fdr_paired_difference_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

heldout_required <- c(
  "schema_version", "design", "tree_case", "simulation_model", "evaluation_model",
  "nsites", "branch_length", "formula", "decision_rule", "training_nulls",
  "empirical_alpha", "empirical_threshold", "holdout_nulls", "holdout_rejections",
  "holdout_rejection_rate"
)
heldout <- read_strict(heldout_path, heldout_required, "Held-out calibration table")
heldout$nsites <- as.integer(heldout$nsites)
heldout$branch_length <- as.numeric(heldout$branch_length)
heldout$holdout_nulls <- as.integer(heldout$holdout_nulls)
heldout$holdout_rejections <- as.integer(heldout$holdout_rejections)
heldout_key <- c(
  "design", "tree_case", "simulation_model", "evaluation_model", "nsites",
  "branch_length", "formula", "decision_rule"
)
if (any(heldout$schema_version != 2L) || !all(heldout$design == "exact_null")) {
  stop("Held-out calibration table must use schema version 2 and design exact_null")
}
if (anyDuplicated(heldout[heldout_key])) stop("Held-out calibration table contains duplicate cells")
if (
  !setequal(heldout$simulation_model, c("JC", "GTR_PF06346")) ||
    !setequal(heldout$tree_case, c("five_external", "sixteen_internal")) ||
    !setequal(heldout$nsites, c(100L, 1000L)) ||
    !setequal(heldout$branch_length, c(4, 8)) ||
    !setequal(heldout$formula, FORMULAS) ||
    !setequal(heldout$decision_rule, RULES) ||
    any(heldout$simulation_model != heldout$evaluation_model)
) {
  stop("Held-out calibration table does not match the prespecified model/tree/site/branch/rule grid")
}
if (nrow(heldout) != 2L * 2L * 2L * 2L * length(FORMULAS) * length(RULES)) {
  stop("Held-out calibration table does not contain exactly one row per prespecified cell")
}
expected_heldout <- expand.grid(
  tree_case = c("five_external", "sixteen_internal"),
  simulation_model = c("JC", "GTR_PF06346"),
  nsites = c(100L, 1000L), branch_length = c(4, 8),
  formula = FORMULAS, decision_rule = RULES,
  stringsAsFactors = FALSE
)
expected_heldout$design <- "exact_null"
expected_heldout$evaluation_model <- expected_heldout$simulation_model
expected_heldout <- expected_heldout[heldout_key]
heldout_key_text <- function(values) do.call(paste, c(values, sep = "\r"))
if (!setequal(heldout_key_text(heldout[heldout_key]), heldout_key_text(expected_heldout))) {
  stop("Held-out calibration table is missing at least one prespecified cell")
}
if (
  any(heldout$training_nulls != training_reps) ||
    any(heldout$holdout_nulls != expected_exact - training_reps) ||
    any(abs(as.numeric(heldout$empirical_alpha) - 0.05) > 1e-12)
) {
  stop("Held-out calibration counts or nominal alpha do not match the requested design")
}
heldout$estimate <- heldout$holdout_rejections / heldout$holdout_nulls
heldout_intervals <- t(mapply(wilson_interval, heldout$holdout_rejections, heldout$holdout_nulls))
heldout$ci_low <- heldout_intervals[, 1L]
heldout$ci_high <- heldout_intervals[, 2L]
heldout$calibration_type <- "Training-derived threshold"

native_null <- exact_summary[
  exact_summary$sample == "holdout" & exact_summary$metric == "null_rejection_rate",
  , drop = FALSE
]
native_null$training_nulls <- training_reps
native_null$empirical_alpha <- 0.05
native_null$empirical_threshold <- 0.05
native_null$holdout_nulls <- native_null$replicates
native_null$holdout_rejections <- round(native_null$estimate * native_null$replicates)
native_null$holdout_rejection_rate <- native_null$estimate
native_null$calibration_type <- "Native α = 0.05"
heldout_columns <- intersect(names(heldout), names(native_null))
null_summary <- rbind(heldout[heldout_columns], native_null[heldout_columns])
null_summary <- decorate(null_summary)
write.table(
  null_summary, file.path(outdir, "heldout_null_plot_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

mixed_summary <- decorate(mixed_summary)
mixed_paired <- decorate(mixed_paired)
for (model in c("JC", "GTR_PF06346")) {
  model_label <- MODEL_LABELS[[model]]
  for (formula in FORMULAS) {
    method <- mixed_summary[
      mixed_summary$simulation_model == model & mixed_summary$formula == formula,
      , drop = FALSE
    ]
    method_title <- paste0(FORMULA_LABELS[[formula]], ": mixed-null FDR and power under ", model_label)
    base <- file.path(
      outdir,
      paste0("figure_", formula, "_fdr_power_2d_model_", safe_name(model))
    )
    save_static(operating_plot(method, method_title), base, width = 11.8, height = 10.5)

    null_method <- null_summary[
      null_summary$simulation_model == model & null_summary$formula == formula,
      , drop = FALSE
    ]
    null_title <- paste0(FORMULA_LABELS[[formula]], ": held-out exact-null calibration under ", model_label)
    null_base <- file.path(
      outdir,
      paste0("figure_", formula, "_null_calibration_2d_model_", safe_name(model))
    )
    save_static(null_calibration_plot(null_method, null_title), null_base, width = 11.8, height = 7.8)

    for (metric in c("empirical_fdr", "power")) {
      metric_rows <- method[method$metric == metric, , drop = FALSE]
      metric_name <- if (metric == "empirical_fdr") "Empirical FDR" else "Power"
      title <- paste0(FORMULA_LABELS[[formula]], ": ", metric_name, " under ", model_label)
      stem <- paste0("figure_", formula, "_", metric, "_3d_model_", safe_name(model))
      static_3d(metric_rows, title, file.path(outdir, paste0(stem, ".pdf")))
      interactive_3d(metric_rows, title, file.path(outdir, paste0(stem, ".html")))
    }
  }

  difference <- mixed_paired[mixed_paired$simulation_model == model, , drop = FALSE]
  difference_title <- paste0("Paired method difference: mixed-null FDR and power under ", model_label)
  difference_base <- file.path(
    outdir, paste0("figure_paired_difference_fdr_power_2d_model_", safe_name(model))
  )
  save_static(
    operating_plot(difference, difference_title, difference = TRUE),
    difference_base, width = 11.8, height = 10.5
  )
  for (metric in c("empirical_fdr", "power")) {
    metric_rows <- difference[difference$metric == metric, , drop = FALSE]
    metric_name <- if (metric == "empirical_fdr") "Empirical FDR" else "Power"
    title <- paste0("Paired weighted − dominant ", metric_name, " under ", model_label)
    stem <- paste0("figure_paired_difference_", metric, "_3d_model_", safe_name(model))
    static_3d(
      metric_rows, title, file.path(outdir, paste0(stem, ".pdf")), difference = TRUE
    )
    interactive_3d(
      metric_rows, title, file.path(outdir, paste0(stem, ".html")), difference = TRUE
    )
  }
}

cat("FDR calibration figures: ", normalizePath(outdir), "\n", sep = "")
