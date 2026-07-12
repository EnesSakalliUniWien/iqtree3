#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) {
  normalizePath(sub("^--file=", "", file_arg))
} else {
  normalizePath("research/satute/tools/manuscript/render_workflow.R", mustWork = FALSE)
}

project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
out_dir <- file.path(project_dir, "artifacts", "releases", "generated-figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blue <- "#153E81"
blue_fill <- adjustcolor(blue, alpha.f = 0.13)
orange <- "#DC840C"
orange_fill <- adjustcolor(orange, alpha.f = 0.15)
brown <- "#8D2F14"
green <- "#58767C"
red <- "#8D2F14"
dark <- "#202020"
mid <- "#777777"
light <- "#EEF1F4"
panel_bg <- "#FBFCFD"
panel_border <- "#D7DEE8"

panel_label <- function(label, x, y) {
  grid.circle(x, y, r = unit(0.021, "npc"),
              gp = gpar(fill = "#263238", col = NA))
  grid.text(label, x, y - 0.001,
            gp = gpar(fontface = "bold", fontsize = 16, col = "white"))
}

panel_box <- function(x, y, w, h) {
  grid.roundrect(x, y, width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(0.014, "npc"),
                 gp = gpar(fill = panel_bg, col = panel_border, lwd = 1))
}

txt <- function(label, x, y, size = 10, col = dark, fontface = "plain", just = "centre",
                family = "Helvetica", rot = 0, lineheight = 0.95) {
  grid.text(label, x, y, just = just, rot = rot,
            gp = gpar(fontsize = size, col = col, fontface = fontface,
                      fontfamily = family, lineheight = lineheight))
}

expr_txt <- function(expr, x, y, size = 12, col = dark, just = "centre") {
  grid.text(expr, x, y, just = just, gp = gpar(fontsize = size, col = col))
}

seg <- function(x0, y0, x1, y1, col = dark, lwd = 1.2, arrow = FALSE, lty = 1) {
  grid.segments(x0, y0, x1, y1,
                gp = gpar(col = col, lwd = lwd, lty = lty,
                          lineend = "round"),
                arrow = if (arrow) grid::arrow(length = unit(0.018, "npc"),
                                               type = "closed") else NULL)
}

curve_arrow <- function(x0, y0, x1, y1, id, col = "black") {
  xs <- seq(x0, x1, length.out = 80)
  if (id == "top") {
    ys <- y0 + (y1 - y0) * seq(0, 1, length.out = 80) + 0.025 * sin(pi * seq(0, 1, length.out = 80))
  } else if (id == "right") {
    ys <- seq(y0, y1, length.out = 80)
    xs <- x0 + (x1 - x0) * seq(0, 1, length.out = 80) + 0.05 * sin(pi * seq(0, 1, length.out = 80))
  } else if (id == "bottom") {
    ys <- y0 + (y1 - y0) * seq(0, 1, length.out = 80) - 0.04 * sin(pi * seq(0, 1, length.out = 80))
  } else {
    ys <- seq(y0, y1, length.out = 80)
    xs <- x0 + (x1 - x0) * seq(0, 1, length.out = 80) - 0.04 * sin(pi * seq(0, 1, length.out = 80))
  }
  grid.lines(xs, ys, gp = gpar(col = col, lwd = 1.5, lineend = "round"))
  n <- length(xs)
  grid.segments(xs[n - 4], ys[n - 4], xs[n], ys[n],
                gp = gpar(col = col, lwd = 1.5, lineend = "round"),
                arrow = grid::arrow(length = unit(0.018, "npc"), type = "closed"))
}

draw_alignment_panel <- function() {
  panel_label("a", 0.082, 0.840)
  txt("Alignment and site pattern", 0.115, 0.842, size = 12.5, fontface = "bold", just = "left")
  txt("taxa", 0.082, 0.690, size = 9, fontface = "bold", rot = 90)
  site_x <- 0.255
  left_end_x <- site_x - 0.025
  right_start_x <- site_x + 0.026

  txt(expression(partialdiff), site_x, 0.785, size = 18)

  taxa <- paste0("s", 1:5, ":")
  left <- c("... ATTAC", "... CTTAT", "... AATAG", "... ATTAG", "... ATTAG")
  site <- c("C", "T", "T", "C", "T")
  right <- c("TTAC ATTAC ...", "TTAT CTTAC ...", "ATAG ATTAC ...",
             "TTAC ATTAC ...", "TCAC GTCAC ...")
  cols <- c(blue, orange, blue, blue, orange)
  y <- seq(0.754, 0.650, length.out = 5)
  for (i in seq_along(taxa)) {
    txt(taxa[i], 0.126, y[i], size = 10, fontface = "bold", just = "right", family = "Courier")
    txt(left[i], left_end_x, y[i], size = 10.5, col = cols[i], fontface = "bold",
        just = "right", family = "Courier")
    txt(site[i], site_x, y[i], size = 11, col = cols[i], fontface = "bold",
        family = "Courier")
    txt(right[i], right_start_x, y[i], size = 10.5, col = cols[i], fontface = "bold",
        just = "left", family = "Courier")
  }
  grid.roundrect(site_x, mean(range(y)), width = unit(0.022, "npc"), height = unit(0.135, "npc"),
                 r = unit(0.015, "npc"), gp = gpar(fill = NA, col = "black", lwd = 0.8))
  txt("site pattern", site_x, 0.616, size = 8.5, col = mid)
}

draw_likelihood_panel <- function() {
  panel_label("b", 0.515, 0.840)
  txt("Branch likelihood at AB", 0.548, 0.842, size = 12.5, fontface = "bold", just = "left")

  # subpatterns
  grid.roundrect(0.494, 0.618, width = unit(0.032, "npc"), height = unit(0.145, "npc"),
                 r = unit(0.015, "npc"), gp = gpar(fill = "white", col = "black", lwd = 0.7))
  txt("C", 0.494, 0.663, size = 14, col = blue, fontface = "bold")
  txt("T", 0.494, 0.618, size = 14, col = blue, fontface = "bold")
  txt("C", 0.494, 0.573, size = 14, col = blue, fontface = "bold")
  txt(expression(partialdiff[A]), 0.489, 0.704, size = 13)

  grid.roundrect(0.882, 0.618, width = unit(0.032, "npc"), height = unit(0.145, "npc"),
                 r = unit(0.015, "npc"), gp = gpar(fill = "white", col = "black", lwd = 0.7))
  txt("T", 0.882, 0.655, size = 14, col = orange, fontface = "bold")
  txt("T", 0.882, 0.585, size = 14, col = orange, fontface = "bold")
  txt(expression(partialdiff[B]), 0.882, 0.704, size = 13)

  # subtree wedges
  grid.polygon(c(0.514, 0.660, 0.514), c(0.714, 0.625, 0.536),
               gp = gpar(fill = blue_fill, col = NA))
  grid.polygon(c(0.862, 0.716, 0.862), c(0.714, 0.625, 0.536),
               gp = gpar(fill = orange_fill, col = NA))
  txt(expression(T[A]), 0.600, 0.716, size = 15)
  txt(expression(T[B]), 0.820, 0.716, size = 15)

  left_tips <- c("s1", "s3", "s4")
  left_y <- c(0.690, 0.626, 0.562)
  for (i in seq_along(left_y)) {
    txt(left_tips[i], 0.540, left_y[i], size = 9, fontface = "bold")
    seg(0.560, left_y[i], 0.675, 0.625, col = blue, lwd = 1.2)
  }
  right_tips <- c("s2", "s5")
  right_y <- c(0.682, 0.570)
  for (i in seq_along(right_y)) {
    txt(right_tips[i], 0.846, right_y[i], size = 9, fontface = "bold")
    seg(0.828, right_y[i], 0.718, 0.625, col = orange, lwd = 1.2)
  }

  grid.circle(0.675, 0.625, r = unit(0.012, "npc"),
              gp = gpar(fill = "white", col = blue, lwd = 0.8))
  grid.circle(0.718, 0.625, r = unit(0.012, "npc"),
              gp = gpar(fill = "white", col = orange, lwd = 0.8))
  txt("A", 0.675, 0.651, size = 15, col = blue, fontface = "bold")
  txt("B", 0.718, 0.651, size = 15, col = orange, fontface = "bold")
  seg(0.687, 0.625, 0.706, 0.625, col = dark, lwd = 1.1, lty = 3)
  txt("AB", 0.697, 0.603, size = 8.5, col = mid)

  # likelihood boxes
  txt(expression(L(partialdiff[A])), 0.642, 0.774, size = 10)
  grid.rect(0.642, 0.733, width = unit(0.055, "npc"), height = unit(0.052, "npc"),
            gp = gpar(fill = "white", col = "black", lwd = 0.5))
  states <- c("A", "C", "G", "T")
  vals <- c("0.0620", "0.0068", "0.0620", "0.0067")
  for (i in 1:4) {
    txt(states[i], 0.626, 0.751 - (i - 1) * 0.012, size = 5.7, fontface = "bold")
    txt(vals[i], 0.652, 0.751 - (i - 1) * 0.012, size = 5.3, family = "Courier")
  }
  txt(expression(L(partialdiff[B])), 0.743, 0.774, size = 10)
  grid.rect(0.743, 0.733, width = unit(0.055, "npc"), height = unit(0.052, "npc"),
            gp = gpar(fill = "white", col = "black", lwd = 0.5))
  vals_b <- c("0.0620", "0.0068", "0.0620", "0.0067")
  for (i in 1:4) {
    txt(states[i], 0.727, 0.751 - (i - 1) * 0.012, size = 5.7, fontface = "bold")
    txt(vals_b[i], 0.753, 0.751 - (i - 1) * 0.012, size = 5.3, family = "Courier")
  }

  txt("coherence coefficients", 0.696, 0.540, size = 10)
  expr_txt(expression(C[partialdiff*i] ==
                        group("(", list(frac(L(partialdiff[A]), P(partialdiff[A])), h[i]), ")") %.%
                        group("(", list(frac(L(partialdiff[B]), P(partialdiff[B])), h[i]), ")")),
           0.696, 0.505, size = 8.4)
}

draw_statistic_panel <- function() {
  panel_label("c", 0.515, 0.400)
  txt("Weighted SatuTe statistic", 0.548, 0.402, size = 12.5, fontface = "bold", just = "left")
  txt("branch-length weights", 0.610, 0.350, size = 7.0, col = mid)
  expr_txt(expression(w[i](t) == exp((lambda[i] - lambda[1]) * t)), 0.660, 0.322, size = 9.8)
  txt("site-pattern coefficient", 0.615, 0.286, size = 7.0, col = mid)
  expr_txt(expression(C[partialdiff]^lambda(t) == sum(w[i](t) * C[partialdiff*i], i == 1, d)),
           0.660, 0.254, size = 8.2)
  txt("alignment average", 0.607, 0.214, size = 7.0, col = mid)
  expr_txt(expression(hat(C)[lambda] == sum(frac(n[partialdiff], n) * C[partialdiff]^lambda(t), partialdiff)),
           0.660, 0.183, size = 8.0)
  expr_txt(expression(Z[lambda] == frac(sqrt(n) * hat(C)[lambda], hat(sigma)[lambda])),
           0.660, 0.121, size = 9.0)
  expr_txt(expression(H[0] * ":" ~~ hat(C)[lambda] %~% N(0, hat(sigma)[lambda]^2 / n)),
           0.660, 0.078, size = 7.0)

  # mode persistence inset
  txt("relative persistence\nover AB", 0.895, 0.354, size = 7.6, col = mid, lineheight = 0.9)
  x0 <- 0.815
  for (i in 1:3) {
    y <- 0.303 - (i - 1) * 0.030
    len <- c(0.105, 0.072, 0.043)[i]
    txt(paste0("mode ", i), x0, y, size = 7.5, just = "right", col = mid)
    grid.roundrect(x0 + 0.006 + len / 2, y, width = unit(len, "npc"),
                   height = unit(0.010, "npc"), r = unit(0.004, "npc"),
                   gp = gpar(fill = c(blue, adjustcolor(blue, alpha.f = 0.62),
                                       adjustcolor(blue, alpha.f = 0.34))[i], col = NA))
  }
}

draw_decision_panel <- function() {
  panel_label("d", 0.082, 0.400)
  txt("Decision under the null", 0.115, 0.402, size = 12.5, fontface = "bold", just = "left")
  x <- seq(-4, 4, length.out = 240)
  y <- dnorm(x)
  x0 <- 0.185
  y0 <- 0.148
  w <- 0.260
  h <- 0.150
  sx <- function(v) x0 + (v + 4) / 8 * w
  sy <- function(v) y0 + v / max(y) * h
  grid.lines(sx(x), sy(y), gp = gpar(col = "black", lwd = 0.8))
  seg(sx(-4.1), y0, sx(4.1), y0, col = "black", lwd = 0.7)
  for (tick in -4:4) {
    seg(sx(tick), y0, sx(tick), y0 - 0.007, col = "black", lwd = 0.5)
    txt(as.character(tick), sx(tick), y0 - 0.022, size = 7.5, col = mid)
  }
  z_alpha <- 2.33
  right <- x >= z_alpha
  grid.polygon(c(sx(z_alpha), sx(x[right]), sx(4), sx(z_alpha)),
               c(y0, sy(y[right]), y0, y0),
               gp = gpar(fill = green, col = NA))
  seg(sx(0), y0, sx(0), sy(dnorm(0)), col = mid, lwd = 0.5, lty = 3)
  seg(sx(z_alpha), y0, sx(z_alpha), sy(dnorm(z_alpha)), col = green, lwd = 0.7)
  txt("H0: independent\nsubtrees", 0.102, y0 + 0.143, size = 7.3,
      just = "left", col = mid, lineheight = 0.95)
  txt("Z under subtree\nindependence", 0.102, y0 + 0.106, size = 7.7,
      just = "left", lineheight = 0.88)
  txt(expression(z[alpha]), sx(z_alpha) + 0.008, y0 - 0.022, size = 8, col = mid)
  txt(expression(alpha == 0.01), sx(3.10), y0 + 0.030, size = 8, col = green)
  txt("saturated", sx(-2.15), y0 + h + 0.040, size = 9.5, col = red)
  txt("phylogenetically\ninformative", sx(2.92), y0 + h + 0.044, size = 9.5, col = green)
  seg(sx(2.34), y0 + h + 0.020, sx(4.10), y0 + h + 0.020, col = green, lwd = 1.1, arrow = TRUE)
  seg(sx(2.24), y0 + h + 0.020, sx(-3.05), y0 + h + 0.020, col = red, lwd = 1.1, arrow = TRUE)
  txt("reject independence when", 0.248, 0.072, size = 9.5)
  expr_txt(expression(Z[lambda] > z[alpha]), 0.420, 0.072, size = 12)
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  txt("SatuTe in a nutshell", 0.500, 0.954, size = 18, fontface = "bold", col = brown)
  panel_box(0.260, 0.675, 0.420, 0.405)
  panel_box(0.735, 0.675, 0.440, 0.405)
  panel_box(0.735, 0.250, 0.440, 0.395)
  panel_box(0.260, 0.250, 0.420, 0.395)

  seg(0.460, 0.865, 0.505, 0.865, col = "#53636F", lwd = 1.1, arrow = TRUE)
  txt("ML-tree inference", 0.482, 0.887, size = 7.5, col = "#53636F")
  seg(0.938, 0.505, 0.938, 0.445, col = "#53636F", lwd = 1.1, arrow = TRUE)
  txt("SatuTe", 0.912, 0.478, size = 8.5, col = brown, fontface = "bold", rot = 90)
  seg(0.505, 0.105, 0.460, 0.105, col = "#53636F", lwd = 1.1, arrow = TRUE)
  seg(0.035, 0.445, 0.035, 0.505, col = "#53636F", lwd = 1.1, arrow = TRUE)
  txt("data\nre-evaluation", 0.064, 0.475, size = 8.2, col = brown,
      fontface = "bold", rot = 90, lineheight = 0.9)

  draw_alignment_panel()
  draw_likelihood_panel()
  draw_statistic_panel()
  draw_decision_panel()
}

write_outputs <- function() {
  pdf_path <- file.path(out_dir, "figure1_satute_weighted_workflow.pdf")
  png_path <- file.path(out_dir, "figure1_satute_weighted_workflow.png")
  svg_path <- file.path(out_dir, "figure1_satute_weighted_workflow.svg")

  cairo_pdf(pdf_path, width = 7.1, height = 7.0, family = "Helvetica")
  draw_figure()
  dev.off()

  png(png_path, width = 7.1 * 320, height = 7.0 * 320, res = 320, bg = "white")
  draw_figure()
  dev.off()

  svg(svg_path, width = 7.1, height = 7.0, family = "Helvetica")
  draw_figure()
  dev.off()

  message("Wrote ", pdf_path)
  message("Wrote ", png_path)
  message("Wrote ", svg_path)
}

write_outputs()
