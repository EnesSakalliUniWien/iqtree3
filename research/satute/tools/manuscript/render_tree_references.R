#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) {
  normalizePath(sub("^--file=", "", file_arg))
} else {
  normalizePath("research/satute/tools/manuscript/render_tree_references.R", mustWork = FALSE)
}

project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
out_dir <- file.path(project_dir, "artifacts", "releases", "generated-figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

protein_tree_path <- file.path(
  project_dir,
  "artifacts/releases/biological/independent_enhanced_sliding_window_20260619_151135",
  "trees/protein_based_2D_ToL.treefile"
)
rrna_tree_path <- file.path(
  project_dir,
  "artifacts/releases/biological/independent_enhanced_sliding_window_20260619_151135",
  "trees/rRNA_based_3D_ToL.treefile"
)

blue <- "#153E81"
light_blue <- "#58767C"
orange <- "#DC840C"
green <- "#8D2F14"
brown <- "#303030"
dark <- "#303030"
mixed <- "grey55"

draw_to_pdf_and_png <- function(stem, width, height, draw_fun) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))

  if (capabilities("cairo")) {
    cairo_pdf(pdf_path, width = width, height = height, family = "Helvetica")
  } else {
    pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  }
  draw_fun()
  dev.off()

  png(
    png_path,
    width = ceiling(width * 260),
    height = ceiling(height * 260),
    res = 260,
    bg = "white"
  )
  draw_fun()
  dev.off()

  invisible(c(pdf = pdf_path, png = png_path))
}

tip_group_from_prefix <- function(labels) {
  first <- sub("_.*", "", labels)
  ifelse(
    first %in% c("Bacteria", "BAC"),
    "Bacteria",
    ifelse(
      first == "Archaea",
      "Archaea",
      ifelse(first %in% c("Eukaryota", "Archaeplastida"), "Eukaryota", "Other")
    )
  )
}

assign_edge_groups <- function(tree, tip_group) {
  n_tip <- Ntip(tree)
  node_group <- rep(NA_character_, n_tip + tree$Nnode)
  node_group[seq_len(n_tip)] <- tip_group
  children <- split(tree$edge[, 2], tree$edge[, 1])

  resolve_group <- function(node) {
    if (!is.na(node_group[node])) {
      return(node_group[node])
    }
    child_nodes <- children[[as.character(node)]]
    groups <- unique(vapply(child_nodes, resolve_group, character(1)))
    groups <- groups[groups != "Other"]
    node_group[node] <<- if (length(groups) == 1) groups else "Mixed"
    node_group[node]
  }

  resolve_group(n_tip + 1)
  node_group[tree$edge[, 2]]
}

descendant_tip_prefix_group <- function(tree) {
  n_tip <- Ntip(tree)
  tip_prefix <- sub("[0-9]+$", "", tree$tip.label)
  tip_prefix[tip_prefix == ""] <- tree$tip.label[tip_prefix == ""]
  assign_edge_groups(tree, tip_prefix)
}

balanced_four <- function(prefix) {
  sprintf(
    "((%s1:0.2,%s2:0.2):0.1,(%s3:0.2,%s4:0.2):0.1)",
    prefix, prefix, prefix, prefix
  )
}

balanced_eight <- function(prefix) {
  sprintf(
    "(((%s1:0.2,%s2:0.2):0.1,(%s3:0.2,%s4:0.2):0.1):0.1,((%s5:0.2,%s6:0.2):0.1,(%s7:0.2,%s8:0.2):0.1):0.1)",
    prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix
  )
}

draw_simulation_panel <- function(
  tree,
  title,
  target_edge,
  node_labels = NULL,
  node_ids = NULL,
  tip_labels = NULL,
  tip_ids = NULL
) {
  edge_group <- descendant_tip_prefix_group(tree)
  edge_color <- ifelse(edge_group == "A", blue, ifelse(edge_group == "B", orange, "grey65"))
  edge_width <- rep(2.2, nrow(tree$edge))
  edge_lty <- rep(1, nrow(tree$edge))
  edge_color[target_edge] <- dark
  edge_width[target_edge] <- 2.8
  edge_lty[target_edge] <- 3

  plot.phylo(
    tree,
    type = "unrooted",
    show.tip.label = FALSE,
    edge.color = edge_color,
    edge.width = edge_width,
    edge.lty = edge_lty,
    use.edge.length = TRUE,
    no.margin = FALSE
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  title(title, adj = 0, font.main = 2, cex.main = 0.98, line = 0.7)

  if (!is.null(node_labels) && !is.null(node_ids)) {
    nodelabels(
      node_labels,
      node = node_ids,
      frame = "circle",
      bg = "white",
      col = dark,
      cex = 0.88,
      font = 2
    )
  }
  if (!is.null(tip_labels) && !is.null(tip_ids)) {
    tiplabels(
      tip_labels,
      tip = tip_ids,
      frame = "circle",
      bg = "white",
      col = dark,
      cex = 0.88,
      font = 2
    )
  }
  edgelabels("varied AB", edge = target_edge, frame = "none", bg = "white", cex = 0.72)
}

draw_simulation_reference <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(0.6, 0.4, 2.3, 0.4), oma = c(0.1, 0, 0.2, 0), family = "Helvetica", xpd = NA)

  five <- read.tree(text = paste0("(", balanced_four("A"), ":0.0,B:1.0);"))
  five_root <- Ntip(five) + 1
  five_target <- which(five$edge[, 2] == match("B", five$tip.label))
  draw_simulation_panel(
    five,
    "(a) Five-taxon external branch",
    five_target,
    node_labels = "A",
    node_ids = five_root,
    tip_labels = "B",
    tip_ids = match("B", five$tip.label)
  )

  a8 <- balanced_eight("A")
  b8 <- balanced_eight("B")
  sixteen <- read.tree(text = paste0("(", a8, ":0.12,", b8, "B:1.0);"))
  sixteen_root <- Ntip(sixteen) + 1
  b_node <- Ntip(sixteen) + which(sixteen$node.label == "B")
  sixteen_target <- which(sixteen$edge[, 2] == b_node)
  draw_simulation_panel(
    sixteen,
    "(b) Sixteen-taxon internal branch",
    sixteen_target,
    node_labels = c("A", "B"),
    node_ids = c(sixteen_root, b_node)
  )
}

draw_domain_hulls <- function(pp, tip_group) {
  fills <- c(
    Bacteria = adjustcolor(blue, alpha.f = 0.08),
    Archaea = adjustcolor(light_blue, alpha.f = 0.11),
    Eukaryota = adjustcolor(orange, alpha.f = 0.14)
  )
  for (group in c("Bacteria", "Archaea", "Eukaryota")) {
    idx <- which(tip_group == group)
    if (length(idx) > 2) {
      hull <- chull(pp$xx[idx], pp$yy[idx])
      polygon(pp$xx[idx][hull], pp$yy[idx][hull], col = fills[group], border = NA)
    }
  }
}

draw_tree_of_life_panel <- function(tree_path, title, show_legend = FALSE) {
  tree <- reorder.phylo(read.tree(tree_path), "postorder")
  tip_group <- tip_group_from_prefix(tree$tip.label)
  edge_group <- assign_edge_groups(tree, tip_group)
  palette <- c(Bacteria = blue, Archaea = light_blue, Eukaryota = orange, Other = "grey75", Mixed = mixed)
  edge_color <- palette[edge_group]
  edge_color[is.na(edge_color)] <- "grey70"
  edge_width <- rep(if (Ntip(tree) > 2500) 0.30 else 0.38, nrow(tree$edge))

  euk_tips <- tree$tip.label[tip_group == "Eukaryota"]
  euk_node <- getMRCA(tree, euk_tips)
  euk_edge <- which(tree$edge[, 2] == euk_node)
  edge_color[euk_edge] <- green
  edge_width[euk_edge] <- 3.4

  yeast_label <- grep(
    "Saccharomyces_cerevisiae_baker_s_yeast|Saccharomyces_cerevisiae$",
    tree$tip.label,
    value = TRUE
  )[1]
  yeast_tip <- match(yeast_label, tree$tip.label)
  yeast_edge <- which(tree$edge[, 2] == yeast_tip)
  edge_color[yeast_edge] <- brown
  edge_width[yeast_edge] <- 2.4

  plot.phylo(
    tree,
    type = "unrooted",
    show.tip.label = FALSE,
    edge.color = edge_color,
    edge.width = edge_width,
    use.edge.length = TRUE,
    no.margin = FALSE
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  draw_domain_hulls(pp, tip_group)

  for (edge_id in c(euk_edge, yeast_edge)) {
    parent <- tree$edge[edge_id, 1]
    child <- tree$edge[edge_id, 2]
    segments(
      pp$xx[parent],
      pp$yy[parent],
      pp$xx[child],
      pp$yy[child],
      col = edge_color[edge_id],
      lwd = edge_width[edge_id] * 1.65,
      lend = "round"
    )
  }

  points(pp$xx[yeast_tip], pp$yy[yeast_tip], pch = 21, bg = brown, col = brown, cex = 1.05)

  euk_parent <- tree$edge[euk_edge, 1]
  euk_child <- tree$edge[euk_edge, 2]
  points(
    mean(pp$xx[c(euk_parent, euk_child)]),
    mean(pp$yy[c(euk_parent, euk_child)]),
    pch = 21,
    bg = green,
    col = green,
    cex = 1.05
  )

  title(title, adj = 0, font.main = 2, cex.main = 1.0, line = 0.7)

  if (show_legend) {
    legend(
      "bottomleft",
      inset = c(0, -0.08),
      bty = "n",
      horiz = TRUE,
      legend = c("Bacteria", "Archaea", "Eukaryota", "branch to Eukaryota", "yeast branch"),
      col = c(blue, light_blue, orange, green, brown),
      lwd = c(2.2, 2.2, 2.2, 3.2, 2.4),
      cex = 0.72
    )
  }
}

draw_tree_of_life_reference <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 1), mar = c(1.0, 0.3, 2.0, 0.3), oma = c(0.25, 0, 0.05, 0), family = "Helvetica", xpd = NA)
  draw_tree_of_life_panel(
    protein_tree_path,
    "(a) 2D Tree of Life: ribosomal protein alignment",
    show_legend = FALSE
  )
  draw_tree_of_life_panel(
    rrna_tree_path,
    "(b) 3D Tree of Life: 16S rRNA alignment",
    show_legend = TRUE
  )
}

draw_to_pdf_and_png("figure_simulation_topology_reference", 7.2, 3.4, draw_simulation_reference)
draw_to_pdf_and_png("figure_tree_of_life_reference", 7.2, 8.0, draw_tree_of_life_reference)

cat("Wrote tree reference figures to ", out_dir, "\n", sep = "")
