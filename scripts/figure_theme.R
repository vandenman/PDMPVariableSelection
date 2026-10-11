# Shared settings for the manuscript figures.
#
# Every figure is saved at the width at which it is printed (the text width of
# the manuscript), so `base_size` is the font size in the printed paper. All
# figures use ggplot2; only network panels are drawn with qgraph.

library(ggplot2)

figure_width_cm <- 15.24                    # \textwidth of proposal/main.tex
figure_dir <- "../proposal/figures"         # relative to the project root

theme_paper <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          legend.position = "bottom")
}

sampler_labels <- c(AdaptiveBoomerang = "Boomerang", PreconditionedBPS = "Bouncy Particle",
                    PreconditionedZigZag = "Zig-Zag", nimble = "NIMBLE")

# Save a figure as PDF (for the manuscript) and PNG (for quick viewing).
save_figure <- function(plot, name, height_cm, width_cm = figure_width_cm,
                        dir = figure_dir) {
  # cairo_pdf keeps hyphens as hyphens (the pdf device draws them as minus signs).
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = width_cm,
         height = height_cm, units = "cm", device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = width_cm,
         height = height_cm, units = "cm", dpi = 300)
}

# Median and interquartile range over replications, for
# stat_summary(fun.data = median_iqr): the points and bars of the figures.
median_iqr <- function(x) {
  x <- x[is.finite(x)]
  data.frame(y = stats::median(x), ymin = stats::quantile(x, 0.25),
             ymax = stats::quantile(x, 0.75))
}

# Log-scale axes labelled at 1, 2, 5, 10, ...: breaks_125() gives the labels
# within a panel's range, and blank_125() extends every panel to the next such
# value below and above what it shows (the interquartile ranges over `groups`,
# or with `iqr = FALSE` the values themselves), so that no point lies beyond
# the outermost label. `panels` are the facet variables.
values_125 <- sort(as.vector(outer(c(1, 2, 5), 10^(-4:6))))
labels_125 <- function(breaks) format(breaks, drop0trailing = TRUE, trim = TRUE, scientific = FALSE)
breaks_125 <- function(limits)
  values_125[values_125 >= limits[1] / 1.001 & values_125 <= limits[2] * 1.001]
blank_125 <- function(d, y, panels, groups = character(), iqr = TRUE) {
  d <- d[is.finite(d[[y]]) & d[[y]] > 0, ]
  keys <- interaction(d[c(panels, groups)], drop = TRUE)
  ends <- do.call(rbind, lapply(split(d, keys), function(g) {
    v <- if (iqr) stats::quantile(g[[y]], c(0.25, 0.75)) else range(g[[y]])
    cbind(g[1, panels, drop = FALSE], lo = v[1], hi = v[2])
  }))
  ends <- do.call(rbind, lapply(split(ends, interaction(ends[panels], drop = TRUE)), function(g)
    cbind(g[c(1, 1), panels, drop = FALSE],
          y = c(max(values_125[values_125 <= min(g$lo)]), min(values_125[values_125 >= max(g$hi)])))))
  geom_blank(data = ends, aes(y = .data$y), inherit.aes = FALSE)
}
