#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(plot3D)
  library(plotly)
})

FORMULA_LEVELS <- c("dominant", "eigenvalue_weighted")
FORMULA_LABELS <- c(
  dominant = "Published dominant coefficient",
  eigenvalue_weighted = "Eigenvalue-weighted coefficient"
)
FORMULA_COLORS <- c(
  dominant = "#0072B2",
  eigenvalue_weighted = "#D55E00"
)
FORMULA_SYMBOLS <- c(
  dominant = "circle",
  eigenvalue_weighted = "diamond"
)
SCENARIO_LABELS <- c(
  true_tree_fixed_lengths = "True tree, fixed lengths",
  true_topology_ml_lengths = "True topology, ML lengths",
  ml_tree = "ML tree"
)
TREE_LABELS <- c(
  five_external = "Five-taxon external branch",
  sixteen_internal = "16-taxon internal branch"
)

parse_args <- function(args) {
  result <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must use --name value syntax; unexpected argument: ", key)
    }
    result[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  result
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required argument --", name)
  value
}

resolve_alias <- function(value, available) {
  if (value %in% available) return(value)
  patterns <- c(
    JC = "^JC$",
    GTR_PF06346 = "^GTR\\{0\\.6676,3\\.7807,4\\.2833,0\\.5354,0\\.8718,1\\.0\\}\\+F\\{0\\.125,0\\.436,0\\.191,0\\.245\\}$",
    LG = "^LG$", WAG = "^WAG$", JTT = "^JTT$", Q.PFAM = "^Q\\.pfam$",
    LG_G4 = "^LG\\+G4\\{", WAG_G4 = "^WAG\\+G4\\{", JTT_G4 = "^JTT\\+G4\\{",
    Q.PFAM_G4 = "^Q\\.pfam\\+G4\\{"
  )
  pattern <- patterns[[value]]
  matches <- if (is.null(pattern)) character() else grep(pattern, available, value = TRUE)
  if (length(matches) != 1L) {
    stop("Could not resolve model '", value, "'. Available values: ", paste(available, collapse = "; "))
  }
  matches[[1L]]
}

short_model_label <- function(model) {
  if (model == "JC") return("JC")
  if (grepl("^GTR\\{0\\.6676,", model)) return("GTR_PF06346")
  if (grepl("^Q\\.pfam", model)) return(sub("^Q\\.pfam", "Q.PFAM", model))
  sub("\\{.*$", "", model)
}

safe_name <- function(value) {
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  gsub("^_|_$", "", value)
}

wilson_interval <- function(successes, total, z = qnorm(0.975)) {
  result <- matrix(NA_real_, nrow = length(total), ncol = 2L)
  valid <- total > 0
  p <- successes[valid] / total[valid]
  denominator <- 1 + z^2 / total[valid]
  center <- (p + z^2 / (2 * total[valid])) / denominator
  half <- z * sqrt(p * (1 - p) / total[valid] + z^2 / (4 * total[valid]^2)) / denominator
  result[valid, 1L] <- pmax(0, center - half)
  result[valid, 2L] <- pmin(1, center + half)
  result
}

site_dashes <- function(site_values) {
  candidates <- c("solid", "22", "42", "13", "73", "F2")
  setNames(rep(candidates, length.out = length(site_values)), site_values)
}

plotly_dashes <- function(site_values) {
  candidates <- c("solid", "dash", "dot", "dashdot", "longdash", "longdashdot")
  setNames(rep(candidates, length.out = length(site_values)), site_values)
}

site_point_shapes <- function(site_values) {
  candidates <- c(21, 24, 22, 25, 23, 8)
  setNames(rep(candidates, length.out = length(site_values)), site_values)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
summary_path <- required_arg(args, "summary")
detail_path <- required_arg(args, "detail")
outdir <- required_arg(args, "outdir")
simulation_arg <- args[["simulation-model"]]
evaluation_arg <- args[["evaluation-model"]]
expected_reps <- if (is.null(args[["expected-reps"]])) NA_integer_ else as.integer(args[["expected-reps"]])
bootstrap_reps <- if (is.null(args[["bootstrap-reps"]])) 5000L else as.integer(args[["bootstrap-reps"]])
decision_rule <- if (is.null(args[["decision-rule"]])) "unadjusted" else args[["decision-rule"]]
decision_columns <- list(
  unadjusted = "decision_unadjusted",
  taxon_bonferroni = "decision_taxon_bonf",
  by_fdr = "decision_fdr"
)
if (!decision_rule %in% names(decision_columns)) {
  stop("--decision-rule must be one of: ", paste(names(decision_columns), collapse = ", "))
}
decision_column <- decision_columns[[decision_rule]]

if (!file.exists(summary_path)) stop("Summary file does not exist: ", summary_path)
if (!file.exists(detail_path)) stop("Detail file does not exist: ", detail_path)
if (bootstrap_reps < 1000L) stop("--bootstrap-reps must be at least 1000")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rows <- read.delim(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "tree_case", "simulation_model", "evaluation_model", "nsites", "branch_length",
  "scenario", "decision_rule", "formula", "evaluated", "informative", "missing_split", "fraction_informative"
)
missing_columns <- setdiff(required_columns, names(rows))
if (length(missing_columns)) stop("Summary is missing columns: ", paste(missing_columns, collapse = ", "))
rows <- rows[rows$decision_rule == decision_rule, , drop = FALSE]
if (!nrow(rows)) stop("No rows remain for decision rule: ", decision_rule)

available_sim <- unique(rows$simulation_model)
available_eval <- unique(rows$evaluation_model)
if (is.null(simulation_arg) || is.null(evaluation_arg)) {
  pairs <- unique(rows[c("simulation_model", "evaluation_model")])
  if (nrow(pairs) != 1L) stop("Pass --simulation-model and --evaluation-model for a multi-model summary")
  simulation_model <- pairs$simulation_model[[1L]]
  evaluation_model <- pairs$evaluation_model[[1L]]
} else {
  simulation_model <- resolve_alias(simulation_arg, available_sim)
  evaluation_model <- resolve_alias(evaluation_arg, available_eval)
}

rows <- rows[
  rows$simulation_model == simulation_model & rows$evaluation_model == evaluation_model,
  , drop = FALSE
]
rows <- rows[rows$formula %in% FORMULA_LEVELS, , drop = FALSE]
if (!nrow(rows)) stop("No rows remain after model/formula filtering")
if (!setequal(unique(rows$formula), FORMULA_LEVELS)) {
  stop("The comparison requires exactly dominant and eigenvalue_weighted rows")
}
if (!is.na(expected_reps)) {
  totals <- rows$evaluated + rows$missing_split
  if (any(totals != expected_reps)) {
    stop("At least one grid cell has evaluated + missing_split != --expected-reps")
  }
}

rows$nsites <- as.integer(rows$nsites)
rows$branch_length <- as.numeric(rows$branch_length)
rows$evaluated <- as.integer(rows$evaluated)
rows$informative <- as.integer(rows$informative)
rows$fraction_informative <- as.numeric(rows$fraction_informative)
intervals <- wilson_interval(rows$informative, rows$evaluated)
rows$ci_low <- intervals[, 1L]
rows$ci_high <- intervals[, 2L]
rows$site_label <- factor(
  paste(format(rows$nsites, big.mark = ",", scientific = FALSE), "sites"),
  levels = paste(format(sort(unique(rows$nsites)), big.mark = ",", scientific = FALSE), "sites")
)
rows$tree_label <- unname(TREE_LABELS[rows$tree_case])
rows$tree_label[is.na(rows$tree_label)] <- rows$tree_case[is.na(rows$tree_label)]
rows$scenario_label <- unname(SCENARIO_LABELS[rows$scenario])
rows$scenario_label[is.na(rows$scenario_label)] <- rows$scenario[is.na(rows$scenario_label)]
rows$tree_label <- factor(rows$tree_label, levels = unique(rows$tree_label))
rows$scenario_label <- factor(rows$scenario_label, levels = unique(rows$scenario_label))

sim_label <- short_model_label(simulation_model)
eval_label <- short_model_label(evaluation_model)
suffix <- paste0("rule_", safe_name(decision_rule), "__sim_", safe_name(sim_label), "__eval_", safe_name(eval_label))
site_styles <- site_dashes(levels(rows$site_label))
n_scenarios <- length(unique(rows$scenario))
n_trees <- length(unique(rows$tree_case))
width_in <- max(9.0, 4.1 * n_scenarios)
height_in <- max(6.5, 3.5 * n_trees + 2.0)

common_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#444444"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 10.5),
    legend.position = "bottom"
  )

save_static <- function(plot, base_path) {
  ggsave(paste0(base_path, ".pdf"), plot, width = width_in, height = height_in, units = "in", device = cairo_pdf)
  ggsave(paste0(base_path, ".svg"), plot, width = width_in, height = height_in, units = "in", device = svglite::svglite)
  ggsave(
    paste0(base_path, ".png"), plot, width = width_in, height = height_in, units = "in",
    dpi = 450, device = ragg::agg_png, background = "white"
  )
  cat("2D PDF: ", paste0(base_path, ".pdf"), "\n", sep = "")
  cat("2D SVG: ", paste0(base_path, ".svg"), "\n", sep = "")
  cat("2D PNG: ", paste0(base_path, ".png"), "\n", sep = "")
}

for (formula_name in FORMULA_LEVELS) {
  method_rows <- rows[rows$formula == formula_name, , drop = FALSE]
  method_label <- FORMULA_LABELS[[formula_name]]
  method_plot <- ggplot(
    method_rows,
    aes(x = branch_length, y = fraction_informative, linetype = site_label, group = site_label)
  ) +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), color = FORMULA_COLORS[[formula_name]], alpha = 0.3, linewidth = 0.35) +
    geom_line(color = FORMULA_COLORS[[formula_name]], linewidth = 0.95) +
    geom_point(shape = if (formula_name == "dominant") 21 else 23, color = FORMULA_COLORS[[formula_name]], fill = "white", size = 2.3, stroke = 0.75) +
    facet_grid(rows = vars(tree_label), cols = vars(scenario_label)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.03))) +
    scale_linetype_manual(values = site_styles) +
    labs(
      title = paste0(method_label, ": sim ", sim_label, ", eval ", eval_label),
      subtitle = paste0("Decision rule: ", decision_rule, "; measured informative fraction with 95% Wilson intervals"),
      x = "Target branch length (log scale)", y = "Fraction classified informative",
      linetype = "Alignment length"
    ) + common_theme
  save_static(method_plot, file.path(outdir, paste0("figure_", formula_name, "_2d_", suffix)))
}

# Pair decisions on the same alignment before estimating the method difference.
detail <- read.delim(detail_path, check.names = FALSE, stringsAsFactors = FALSE)
detail_required <- c(
  "tree_case", "simulation_model", "evaluation_model", "nsites", "branch_length", "replicate",
  "seed", "target_split", "scenario", "formula", "target_found", decision_column
)
detail_missing <- setdiff(detail_required, names(detail))
if (length(detail_missing)) stop("Detail is missing columns: ", paste(detail_missing, collapse = ", "))
detail <- detail[
  detail$simulation_model == simulation_model & detail$evaluation_model == evaluation_model &
    detail$formula %in% FORMULA_LEVELS,
  , drop = FALSE
]
pair_keys <- c(
  "tree_case", "simulation_model", "evaluation_model", "nsites", "branch_length", "replicate",
  "seed", "target_split", "scenario"
)
dominant <- detail[detail$formula == "dominant", c(pair_keys, "target_found", decision_column)]
weighted <- detail[detail$formula == "eigenvalue_weighted", c(pair_keys, "target_found", decision_column)]
names(dominant)[(length(pair_keys) + 1L):(length(pair_keys) + 2L)] <- c("target_found_dominant", "decision_dominant")
names(weighted)[(length(pair_keys) + 1L):(length(pair_keys) + 2L)] <- c("target_found_weighted", "decision_weighted")
if (anyDuplicated(dominant[pair_keys]) || anyDuplicated(weighted[pair_keys])) stop("Duplicate formula rows prevent one-to-one pairing")
paired <- merge(dominant, weighted, by = pair_keys, all = FALSE, sort = FALSE)
if (nrow(paired) != nrow(dominant) || nrow(paired) != nrow(weighted)) stop("Formula rows are not a complete one-to-one pair set")
paired$valid_pair <- paired$target_found_dominant == 1L & paired$target_found_weighted == 1L
paired$delta <- ifelse(
  paired$valid_pair,
  as.integer(paired$decision_weighted == "informative") - as.integer(paired$decision_dominant == "informative"),
  NA_integer_
)

group_keys <- c("tree_case", "nsites", "branch_length", "scenario")
group_id <- interaction(paired[group_keys], drop = TRUE, lex.order = TRUE)
set.seed(20260714)
paired_summary <- do.call(rbind, lapply(split(paired, group_id), function(group) {
  valid_delta <- group$delta[!is.na(group$delta)]
  n <- length(valid_delta)
  counts <- table(factor(valid_delta, levels = c(-1L, 0L, 1L)))
  if (n == 0L) {
    estimate <- ci_low <- ci_high <- NA_real_
  } else {
    estimate <- mean(valid_delta)
    if (counts[[1L]] + counts[[3L]] == 0L) {
      ci_low <- ci_high <- 0
    } else {
      bootstrap_counts <- rmultinom(bootstrap_reps, n, prob = as.numeric(counts) / n)
      bootstrap_delta <- (bootstrap_counts[3L, ] - bootstrap_counts[1L, ]) / n
      limits <- quantile(bootstrap_delta, c(0.025, 0.975), names = FALSE, type = 6)
      ci_low <- limits[[1L]]
      ci_high <- limits[[2L]]
    }
  }
  discordant <- counts[[1L]] + counts[[3L]]
  p_value <- if (discordant == 0L) 1 else binom.test(counts[[3L]], discordant, p = 0.5)$p.value
  data.frame(
    tree_case = group$tree_case[[1L]], nsites = as.integer(group$nsites[[1L]]),
    branch_length = as.numeric(group$branch_length[[1L]]), scenario = group$scenario[[1L]],
    total_pairs = nrow(group), evaluated_pairs = n, missing_pairs = nrow(group) - n,
    dominant_only_informative = as.integer(counts[[1L]]), concordant = as.integer(counts[[2L]]),
    weighted_only_informative = as.integer(counts[[3L]]), delta_fraction_informative = estimate,
    ci_low = ci_low, ci_high = ci_high, mcnemar_exact_p = p_value,
    stringsAsFactors = FALSE
  )
}))
paired_summary <- paired_summary[order(
  paired_summary$tree_case, paired_summary$scenario, paired_summary$nsites, paired_summary$branch_length
), ]
paired_summary$tree_label <- unname(TREE_LABELS[paired_summary$tree_case])
paired_summary$tree_label[is.na(paired_summary$tree_label)] <- paired_summary$tree_case[is.na(paired_summary$tree_label)]
paired_summary$scenario_label <- unname(SCENARIO_LABELS[paired_summary$scenario])
paired_summary$scenario_label[is.na(paired_summary$scenario_label)] <- paired_summary$scenario[is.na(paired_summary$scenario_label)]
paired_summary$tree_label <- factor(paired_summary$tree_label, levels = levels(rows$tree_label))
paired_summary$scenario_label <- factor(paired_summary$scenario_label, levels = levels(rows$scenario_label))
paired_summary$site_label <- factor(
  paste(format(paired_summary$nsites, big.mark = ",", scientific = FALSE), "sites"),
  levels = levels(rows$site_label)
)
write.table(
  paired_summary[setdiff(names(paired_summary), c("tree_label", "scenario_label", "site_label"))],
  file.path(outdir, paste0("paired_formula_difference_", suffix, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

delta_limit <- max(0.05, max(abs(c(paired_summary$ci_low, paired_summary$ci_high)), na.rm = TRUE) * 1.08)
delta_limit <- min(1, delta_limit)
difference_plot <- ggplot(
  paired_summary,
  aes(
    x = branch_length, y = delta_fraction_informative, color = site_label,
    linetype = site_label, group = site_label
  )
) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = 0.45) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), alpha = 0.38, linewidth = 0.45) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.2) +
  facet_grid(rows = vars(tree_label), cols = vars(scenario_label)) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-delta_limit, delta_limit), labels = scales::label_percent(accuracy = 1)) +
  scale_color_viridis_d(option = "D", end = 0.85) +
  scale_linetype_manual(values = site_styles) +
  labs(
    title = paste0("Paired formula difference: sim ", sim_label, ", eval ", eval_label),
    subtitle = paste0("Decision rule: ", decision_rule, "; eigenvalue-weighted minus dominant; 95% paired bootstrap intervals"),
    x = "Target branch length (log scale)", y = "Difference in informative fraction",
    color = "Alignment length", linetype = "Alignment length"
  ) + common_theme
save_static(difference_plot, file.path(outdir, paste0("figure_paired_difference_2d_", suffix)))

# True 3D line plots: x = branch length, y = alignment length, z = the observed
# metric. Each trajectory is a measured alignment-length slice; no surface is
# interpolated between site lengths.
build_3d <- function(data, title, z_column, z_title, z_range, line_color, html_path, difference = FALSE) {
  selection_keys <- unique(data[c("tree_case", "tree_label", "scenario", "scenario_label")])
  selection_keys <- selection_keys[order(selection_keys$tree_case, selection_keys$scenario), , drop = FALSE]
  site_values <- sort(unique(data$nsites))
  dash_styles <- plotly_dashes(as.character(site_values))
  traces_per_selection <- length(site_values)
  figure <- plot_ly(width = 1280, height = 860)
  site_colors <- setNames(viridisLite::viridis(length(site_values), option = "D", end = 0.85), site_values)

  for (selection_index in seq_len(nrow(selection_keys))) {
    selection <- selection_keys[selection_index, , drop = FALSE]
    selected <- data[data$tree_case == selection$tree_case & data$scenario == selection$scenario, , drop = FALSE]
    for (site_value in site_values) {
      trace <- selected[selected$nsites == site_value, , drop = FALSE]
      trace <- trace[order(trace$branch_length), , drop = FALSE]
      z <- trace[[z_column]]
      if (difference) {
        hover <- paste0(
          "<b>Eigenvalue-weighted − dominant</b>",
          "<br>Branch length: ", format(trace$branch_length, trim = TRUE),
          "<br>Alignment length: ", format(site_value, big.mark = ","),
          "<br>Difference: ", sprintf("%+.3f", z),
          "<br>95% paired CI: [", sprintf("%+.3f", trace$ci_low), ", ", sprintf("%+.3f", trace$ci_high), "]",
          "<br>Weighted only: ", trace$weighted_only_informative,
          "<br>Dominant only: ", trace$dominant_only_informative,
          "<br>Evaluated pairs: ", trace$evaluated_pairs
        )
      } else {
        hover <- paste0(
          "<br>Branch length: ", format(trace$branch_length, trim = TRUE),
          "<br>Alignment length: ", format(site_value, big.mark = ","),
          "<br>Informative: ", sprintf("%.3f", z),
          "<br>95% CI: [", sprintf("%.3f", trace$ci_low), ", ", sprintf("%.3f", trace$ci_high), "]",
          "<br>Evaluated: ", trace$evaluated
        )
      }
      trace_color <- if (difference) site_colors[[as.character(site_value)]] else line_color
      figure <- add_trace(
        figure, x = trace$branch_length, y = trace$nsites, z = z,
        type = "scatter3d", mode = "lines+markers",
        name = paste0(format(site_value, big.mark = ","), " sites"),
        legendgroup = as.character(site_value), showlegend = selection_index == 1L,
        visible = selection_index == 1L,
        line = list(color = trace_color, width = 6, dash = dash_styles[[as.character(site_value)]]),
        marker = list(color = trace_color, size = 4.5, symbol = "circle", line = list(color = "white", width = 0.8)),
        text = hover, hoverinfo = "text"
      )
    }
  }

  buttons <- vector("list", nrow(selection_keys))
  total_traces <- nrow(selection_keys) * traces_per_selection
  for (selection_index in seq_len(nrow(selection_keys))) {
    visible <- rep(FALSE, total_traces)
    first <- (selection_index - 1L) * traces_per_selection + 1L
    visible[first:(first + traces_per_selection - 1L)] <- TRUE
    selection <- selection_keys[selection_index, , drop = FALSE]
    selected_title <- paste0(title, "<br><sup>", selection$tree_label, " · ", selection$scenario_label, "</sup>")
    buttons[[selection_index]] <- list(
      method = "update", args = list(list(visible = visible), list(title = list(text = selected_title))),
      label = paste(selection$tree_label, selection$scenario_label, sep = " · ")
    )
  }

  initial <- selection_keys[1L, , drop = FALSE]
  initial_title <- paste0(title, "<br><sup>", initial$tree_label, " · ", initial$scenario_label, "</sup>")
  figure <- layout(
    figure,
    title = list(text = initial_title, x = 0.5),
    scene = list(
      xaxis = list(title = "Target branch length", type = "log", showspikes = FALSE),
      yaxis = list(title = "Alignment length (sites)", type = "log", showspikes = FALSE),
      zaxis = list(title = z_title, range = z_range, zeroline = TRUE, zerolinecolor = "#444444", showspikes = FALSE),
      aspectmode = "manual", aspectratio = list(x = 1.45, y = 1.0, z = 0.9),
      camera = list(eye = list(x = 1.55, y = -1.65, z = 1.15))
    ),
    legend = list(orientation = "h", x = 0, y = -0.12),
    margin = list(l = 20, r = 20, b = 100, t = 110),
    updatemenus = list(list(
      type = "dropdown", direction = "down", x = 0.01, y = 1.08,
      xanchor = "left", yanchor = "top", buttons = buttons
    )),
    annotations = list(list(
      text = "Observed alignment-length slices only; no surface interpolation.",
      x = 0.5, y = -0.17, xref = "paper", yref = "paper", showarrow = FALSE,
      font = list(size = 12, color = "#555555")
    ))
  )
  figure <- config(figure, responsive = TRUE, displaylogo = FALSE, scrollZoom = TRUE)
  htmlwidgets::saveWidget(figure, html_path, selfcontained = TRUE, title = title)
  dependency_dir <- sub("\\.html$", "_files", html_path)
  if (dir.exists(dependency_dir) && !grepl(basename(dependency_dir), paste(readLines(html_path, warn = FALSE), collapse = ""), fixed = TRUE)) {
    unlink(dependency_dir, recursive = TRUE)
  }
  cat("3D HTML: ", html_path, "\n", sep = "")
}

# A vector, publication-oriented companion to the interactive Plotly figure.
# Every tree/scenario selection is shown as a panel on the same PDF page, so no
# information is hidden behind the HTML dropdown. The x and y coordinates use
# the same log10 transformations as the interactive plot's logarithmic axes.
build_static_3d <- function(data, title, z_column, z_title, z_range, line_color, pdf_path, difference = FALSE) {
  tree_values <- unique(as.character(data$tree_label))
  scenario_values <- unique(as.character(data$scenario_label))
  site_values <- sort(unique(data$nsites))
  line_styles <- site_dashes(as.character(site_values))
  point_shapes <- site_point_shapes(as.character(site_values))
  site_colors <- setNames(viridisLite::viridis(length(site_values), option = "D", end = 0.85), site_values)
  safe_log_range <- function(values) {
    result <- range(log10(values), finite = TRUE)
    if (!all(is.finite(result))) return(c(-1, 1))
    if (result[[1L]] == result[[2L]]) result <- result + c(-0.5, 0.5)
    result
  }
  x_range <- safe_log_range(data$branch_length)
  y_range <- safe_log_range(data$nsites)

  panel_ids <- matrix(
    seq_len(length(tree_values) * length(scenario_values)),
    nrow = length(tree_values), ncol = length(scenario_values), byrow = TRUE
  )
  legend_id <- max(panel_ids) + 1L
  layout_matrix <- rbind(panel_ids, rep(legend_id, ncol(panel_ids)))

  grDevices::cairo_pdf(
    pdf_path,
    width = max(10.5, 4.25 * length(scenario_values)),
    height = max(7.5, 3.45 * length(tree_values) + 1.45),
    onefile = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(layout_matrix, heights = c(rep(1, length(tree_values)), 0.24))
  par(
    oma = c(2.15, 0.6, 4.6, 0.6), mar = c(1.0, 1.0, 3.0, 0.45),
    family = "sans", fg = "#333333", col.axis = "#444444", col.lab = "#333333"
  )

  for (tree_label_value in tree_values) {
    for (scenario_label_value in scenario_values) {
      selected <- data[
        as.character(data$tree_label) == tree_label_value &
          as.character(data$scenario_label) == scenario_label_value,
        , drop = FALSE
      ]
      if (!nrow(selected)) {
        plot.new()
        title(main = paste(tree_label_value, scenario_label_value, sep = "\n"), cex.main = 0.9)
        text(0.5, 0.5, "No observations", col = "#666666")
        next
      }
      if (length(unique(selected$branch_length)) < 2L || length(unique(selected$nsites)) < 2L) {
        plot.new()
        title(main = paste(tree_label_value, scenario_label_value, sep = "\n"), cex.main = 0.9)
        text(0.5, 0.5, "At least two branch and site values are required for a 3D projection", col = "#666666", cex = 0.8)
        next
      }

      first_trace <- TRUE
      for (site_value in site_values) {
        trace <- selected[selected$nsites == site_value, , drop = FALSE]
        trace <- trace[order(trace$branch_length), , drop = FALSE]
        if (!nrow(trace)) next
        trace_color <- if (difference) site_colors[[as.character(site_value)]] else line_color
        plot3D::scatter3D(
          x = log10(trace$branch_length), y = log10(trace$nsites), z = trace[[z_column]],
          type = "o", add = !first_trace, colvar = NULL, col = trace_color,
          colkey = FALSE, lty = line_styles[[as.character(site_value)]], lwd = 2.0,
          pch = point_shapes[[as.character(site_value)]], cex = 0.72, bg = "white",
          xlim = x_range,
          ylim = y_range, zlim = z_range,
          xlab = "log10 branch length", ylab = "log10 sites", zlab = z_title,
          theta = 42, phi = 23, d = 2.35, expand = 0.78,
          bty = "u", ticktype = "detailed", nticks = 4,
          col.panel = "#FAFAFA", col.grid = "#DDDDDD", col.axis = "#555555",
          lwd.panel = 0.65, lwd.grid = 0.55,
          main = if (first_trace) paste(tree_label_value, scenario_label_value, sep = "\n") else NULL,
          cex.main = 0.86
        )
        first_trace <- FALSE
      }
    }
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend_colors <- if (difference) unname(site_colors[as.character(site_values)]) else rep(line_color, length(site_values))
  legend(
    "center",
    legend = paste(format(site_values, big.mark = ",", scientific = FALSE), "sites"),
    col = legend_colors, lty = unname(line_styles[as.character(site_values)]),
    pch = unname(point_shapes[as.character(site_values)]), pt.bg = "white",
    lwd = 2.0, pt.cex = 0.9, horiz = TRUE, bty = "n", cex = 0.88,
    x.intersp = 0.75, y.intersp = 0.8
  )
  mtext(title, side = 3, outer = TRUE, line = 2.75, font = 2, cex = 1.22, col = "#222222")
  mtext(
    "Static 3D view of measured alignment-length trajectories",
    side = 3, outer = TRUE, line = 1.25, cex = 0.92, col = "#555555"
  )
  mtext(
    "Points are observed means; axes x and y are log10-scaled; no surface interpolation.",
    side = 1, outer = TRUE, line = 0.45, cex = 0.82, col = "#555555"
  )
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  cat("3D PDF: ", pdf_path, "\n", sep = "")
}

for (formula_name in FORMULA_LEVELS) {
  method_rows <- rows[rows$formula == formula_name, , drop = FALSE]
  static_title <- paste0(FORMULA_LABELS[[formula_name]], ": sim ", sim_label, ", eval ", eval_label)
  build_3d(
    method_rows,
    static_title,
    "fraction_informative", "Fraction informative", c(0, 1), FORMULA_COLORS[[formula_name]],
    file.path(outdir, paste0("figure_", formula_name, "_3d_", suffix, ".html"))
  )
  build_static_3d(
    method_rows, static_title,
    "fraction_informative", "Fraction informative", c(0, 1), FORMULA_COLORS[[formula_name]],
    file.path(outdir, paste0("figure_", formula_name, "_3d_", suffix, ".pdf"))
  )
}
static_difference_title <- paste0("Paired formula difference (", decision_rule, "): sim ", sim_label, ", eval ", eval_label)
build_3d(
  paired_summary,
  static_difference_title,
  "delta_fraction_informative", "Weighted − dominant", c(-delta_limit, delta_limit), "#6A3D9A",
  file.path(outdir, paste0("figure_paired_difference_3d_", suffix, ".html")), difference = TRUE
)
build_static_3d(
  paired_summary, static_difference_title,
  "delta_fraction_informative", "Weighted − dominant", c(-delta_limit, delta_limit), "#6A3D9A",
  file.path(outdir, paste0("figure_paired_difference_3d_", suffix, ".pdf")), difference = TRUE
)
