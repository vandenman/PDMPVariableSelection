rm(list = ls())
library(ggplot2)
library(patchwork)
library(dplyr)
library(PDMPSamplersR)
# PDMPSamplersR::use_local_pdmpsamplers("/home/don/github/juliapkgs/PDMPSamplers.jl")

dbanana <- function(X, Y, mu = 0, a = 2, b = 0.5,
                    rho = 0.5, curvature = 1) {

  sx2 <- 1 / (2 * a)
  C <- 2 * curvature^2 * sx2^2 + 1 / (2 * b)
  c <- sign(rho) * sqrt(rho^2 * C / (sx2 * (1 - rho^2)))
  u <- X - mu
  meanY <- curvature * u^2 + c * u

  dnorm(X, mu, sqrt(sx2)) * dnorm(Y, meanY, 1 / sqrt(2 * b))
}


grad_log_banana <- function(X, Y, mu = 0, a = 2, b = 0.5,
                            rho = 0.5, curvature = 1) {

  sx2 <- 1 / (2 * a)
  C <- 2 * curvature^2 * sx2^2 + 1 / (2 * b)
  c <- sign(rho) * sqrt(rho^2 * C / (sx2 * (1 - rho^2)))
  u <- X - mu
  r <- Y - curvature * u^2 - c * u

  c(X = -2 * a * u + 2 * b * r * (2 * curvature * u + c),
    Y = -2 * b * r
  )
}


# Suggested parameters
mu <- -1
a <- 1 / 2
b <- 2

rho <- 0.5        # try -0.7, 0, 0.7
curvature <- 0.5    # try 0.5, 1, 2


# Grid for plotting
np <- 2^8

x <- seq(-3.5, 3, length.out = np)
y <- seq(-2, 6, length.out = np)

df <- data.frame(
  x = rep(x, length(y)),
  y = rep(y, each = length(x))
)

df$z <- dbanana(
  df$x,
  df$y,
  mu = mu,
  a = a,
  b = b,
  rho = rho,
  curvature = curvature
)

ggplot(df, aes(x, y, z = z)) +
  geom_contour(linewidth = 0.25) +
  theme_bw()

x0 <- c(-1, 0)
grad_log_banana2 <- function(x) {
  -grad_log_banana(x[1], x[2], mu, a, b, rho = rho, curvature = curvature)
}


#
# dbanana <- function(X, Y, mu, a, b) {
#   pX <- dnorm(X, mu, 1 / sqrt(2 * a))
#   pY <- dnorm(Y, X^2, 1/sqrt(2 * b))
#   return(pX * pY)
# }
#
# grad_log_banana <- function(X, Y, mu, a, b) {
#
#   r <- Y - X^2
#   return(c(
#     X = -2 * a * (X - mu) + 4 * b * X * r,
#     Y = -2 * b * r
#   ))
# }
# mu <- 0
# a <- 1 / 2
# b <- 2
#
# np <- 2^8
# x <- seq(-3, 3, length.out = np)
# y <- seq(-1, 5, length.out = np)
#
# df <- data.frame(
#   x = rep(x, length(y)),
#   y = rep(y, each = length(x))
# )
# df$z <- dbanana(df$x, df$y, mu, a, b)
#
# ggplot(df, aes(x, y, z = z)) +
#   # geom_contour_filled() +
#   geom_contour(linewidth = .25) +
#   scale_color_discrete(palette = "grey")+
#   scale_fill_discrete(palette = "grey") +
#   theme_bw()
#   # geom_raster() +
# x0 <- c(-1, 1)
# grad_log_banana2 <- function(x) {
#   -grad_log_banana(x[1], x[2], mu, a, b)
# }

set.seed(44)
trace_zz   <- pdmp_sample(grad_log_banana2, 2, flow = "ZigZag",            algorithm = "GridThinningStrategy", x0 = x0, T = 5000)
trace_bps  <- pdmp_sample(grad_log_banana2, 2, flow = "BouncyParticle",    algorithm = "GridThinningStrategy", x0 = x0, T = 5000)

sx2 <- 1 / (2 * a)
C <- 2 * curvature^2 * sx2^2 + 1 / (2 * b)
sy2 <- C / (1 - rho^2)

ref_mean <- c(mu, curvature * sx2)

ref_cov <- matrix(c(
  sx2,                     rho * sqrt(sx2 * sy2),
  rho * sqrt(sx2 * sy2),   sy2
), byrow = TRUE, 2, 2)
trace_boom <- pdmp_sample(grad_log_banana2, 2, flow = "Boomerang", algorithm = "GridThinningStrategy", x0 = x0, T = 5000, flow_mean = ref_mean, flow_cov = ref_cov)

# Discretize the trace for plotting purposes, for estimators of posterior quantities, something like mean(trace_zz) or cov(trace_zz) is preferred.
# as this uses continuous time estimators based on analytic results
samples_zz   <- discretize(trace_zz, 0.01)
samples_bps  <- discretize(trace_bps, 0.01)
samples_boom <- discretize(trace_boom, 0.01)

to_tib <- function(samples, nmax = nrow(samples)) {
  tibble(t = (0:(nmax - 1)) * 0.01, x = samples[1:nmax, 1], y = samples[1:nmax, 2])
}
nmax <- 3000
df_zz   <- to_tib(samples_zz,   nmax = nmax)
df_bps  <- to_tib(samples_bps,  nmax = nmax)
df_boom <- to_tib(samples_boom, nmax = nmax)

# df_zz   <- to_tib(samples_zz)
# df_bps  <- to_tib(samples_bps)
# df_boom <- to_tib(samples_boom)

add_trace <- function(df) {
  list(
    geom_path(data = df, mapping = aes(x = x, y = y), inherit.aes = FALSE,
              alpha = .9, color = "gold3"),
    geom_point(data = df[1, ], mapping = aes(x = x, y = y), inherit.aes = FALSE,
               alpha = .9, color = "white", fill = "grey", shape = 21, size = 2.5)
  )
}

# base_plot <- ggplot(df, aes(x, y, z = z)) +
#   geom_raster(aes(fill = z), alpha = .95, show.legend = FALSE)

base_plot <- ggplot(df, aes(x, y, z = z)) +
  geom_contour(linewidth = .2, bins = 15)

p1 <- base_plot + add_trace(df_zz)   + ggtitle("ZigZag")# + theme(title = element_text(hjust = .5))
p2 <- base_plot + add_trace(df_bps)  + ggtitle("BPS") # + theme(title = element_text(hjust = .5))
p3 <- base_plot + add_trace(df_boom) + ggtitle("Boomerang") #+ theme(title = element_text(hjust = .5))

p123 <- (p1 + p2 + p3) &
  scale_x_continuous(breaks = seq(-3, 2, 1), labels = NULL, limits = c(-3, 1)) &
  scale_y_continuous(breaks = seq(-2, 4, 2), labels = NULL, limits = c(-2, 4)) &
  labs(x = NULL, y = NULL) & theme_bw(base_size = 18) &
  theme(aspect.ratio = 1, axis.ticks.length = unit(1, "mm"),plot.title = element_text(hjust = .5)) &
  theme(margins = unit(c(0, 0, 0, 0), "cm"))

p123

h <- 9
w <- 16
# ggsave(file.path("figures", "dynamics.pdf"), p123, height = h, width = w, units = "cm")
# ggsave(file.path("/home/don/hdd/surfdrive/Postdoc/Psychometrika_PDMP_Variable_Selection/proposal/figures", "dynamics.pdf"), p123, height = h, width = w, units = "cm")


colors <- viridisLite::viridis(3)
# scales::show_col(colors[1:2])


add_trace <- function(df) {
  list(
    geom_path(
      data = df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = .6,
      alpha = .85,
      color = colors[2]
    ),
    geom_point(
      data = df[1, ],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = .7,
      color = "grey20",
      fill = "white"
    )
  )
}

base_plot <- ggplot(df, aes(x, y, z = z)) + geom_contour(bins = 9, linewidth = .2, color = colors[1])

p1 <- base_plot + add_trace(df_zz, head_len_mm = 2.0, arrow_len = 0.05) + ggtitle("ZigZag")
p2 <- base_plot + add_trace(df_bps, head_len_mm = 2.0, arrow_len = 0.05) + ggtitle("BPS")
p3 <- base_plot + add_trace(df_boom, arrow_len = 0, head_len_mm = 0.5) + ggtitle("Boomerang")

p123 <- (p1 + p2 + p3) &
  scale_x_continuous(breaks = seq(-3, 2, 1), labels = NULL, limits = c(-3.5, 1.5), expand = expansion(mult = 0)) &
  scale_y_continuous(breaks = seq(-2, 4, 2), labels = NULL, limits = c(-2, 4),     expand = expansion(mult = 0)) &
  labs(x = NULL, y = NULL) &
  coord_fixed() & theme_bw(base_size = 18) &
  theme(
    axis.ticks.length = unit(1, "mm"),
    axis.line = element_line(linewidth = .4, color = "grey30"),
    # plot.title = element_text(hjust = .5),
    plot.title=element_text(hjust = .1, margin=margin(t=40,b=-30)),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid.major = element_line(color = "#EBEBEBFF", linewidth = .2),
    panel.grid.minor = element_blank()
  )

p123
h <- 9
w <- 16
ggsave(file.path("/home/don/surfdrive/Postdoc/Psychometrika_PDMP_Variable_Selection/PDMPVariableSelection/figures", "dynamics_new.pdf"), p123, height = h, width = w, units = "cm")



mu <- -1
a <- 1 / 2
b <- 1

rho <- 0.5        # try -0.7, 0, 0.7
curvature <- .5    # try 0.5, 1, 2


# Grid for plotting
np <- 2^8

x <- seq(-3.5, 3, length.out = np)
y <- seq(-2, 6, length.out = np)

df <- data.frame(x = rep(x, length(y)), y = rep(y, each = length(x)))

df$z <- dbanana(df$x, df$y, mu = mu, a = a, b = b,rho = rho, curvature = curvature)

ggplot(df, aes(x, y, z = z)) +
  geom_contour(linewidth = 0.25) +
  theme_bw()

x0 <- c(-1, 0)
grad_log_banana2 <- function(x) {
  -grad_log_banana(x[1], x[2], mu, a, b, rho = rho, curvature = curvature)
}


sx2 <- 1 / (2 * a)
C <- 2 * curvature^2 * sx2^2 + 1 / (2 * b)
cc <- sign(rho) * sqrt(rho^2 * C / (sx2 * (1 - rho^2)))

fx0 <- dnorm(0, mean = mu, sd = sqrt(sx2))

fy0 <- integrate(
  function(u) { dnorm(u) * dnorm(0, mean = curvature * sx2 * u^2 + cc * sqrt(sx2) * u, sd = 1 / sqrt(2 * b)) },
  lower = -Inf, upper = Inf
)$value

density_at_0 <- c(x = fx0, y = fy0)

set.seed(44)
trace_zz_vs   <- pdmp_sample(grad_log_banana2, 2, flow = "ZigZag",            algorithm = "GridThinningStrategy", x0 = x0, T = 5000, sticky = TRUE,
                             model_prior = bernoulli(.3), parameter_prior = density_at_0, can_stick = c(TRUE, TRUE))
samples_zz_vs   <- discretize(trace_zz_vs, 0.01)
nmax <- 10000
df_zz_vs   <- to_tib(samples_zz_vs,   nmax = nmax)


colors <- viridisLite::viridis(3)
# scales::show_col(colors[1:2])

head(df_zz_vs)

add_trace <- function(df) {
  list(
    geom_path(
      data = df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = .6,
      alpha = .85,
      color = colors[2]
    ),
    geom_point(
      data = df[1, ],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = .7,
      color = "grey20",
      fill = "white"
    )
  )
}

add_trace <- function(df, tol = 1e-6, arrow_frac = 0.6) {

  sgn <- function(x)
    ifelse(abs(x) < tol, 0L, sign(x))

  direction_changes <- df |>
    dplyr::mutate(
      x_prev = dplyr::lag(x),
      y_prev = dplyr::lag(y),
      x_next = dplyr::lead(x),
      y_next = dplyr::lead(y),

      dx = x - x_prev,
      dy = y - y_prev,

      sx = sgn(dx),
      sy = sgn(dy),

      sx_prev = dplyr::lag(sx),
      sy_prev = dplyr::lag(sy),

      changed = sx != sx_prev | sy != sy_prev,

      # Exclude sticky events: neither component should be
      # stationary before or after the direction change
      nonsticky =
        sx      != 0 & sy      != 0 &
        sx_prev != 0 & sy_prev != 0,

      # Short outgoing arrow, pointing in the new direction
      x_arrow = x + arrow_frac * (x_next - x),
      y_arrow = y + arrow_frac * (y_next - y)
    ) |>
    dplyr::filter(
      changed,
      nonsticky,
      !is.na(x_next),
      !is.na(y_next)
    )

  list(
    geom_path(
      data = df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = .6,
      alpha = .85,
      color = colors[2]
    ),

    geom_segment(
      data = direction_changes,
      aes(x = x, y = y, xend = x_arrow, yend = y_arrow),
      inherit.aes = FALSE,
      linewidth = .6,
      color = colors[2],
      # alpha = .75,
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(2, "mm")
      )
    ),

    geom_point(
      data = df[1, ],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = .7,
      color = "grey20",
      fill = "white"
    )
  )
}

add_trace <- function(df, tol = 1e-6, arrow_len = 0.35, head_len_mm = 1.5) {

  sgn <- function(x) {
    ifelse(abs(x) < tol, 0L, sign(x))
  }

  direction_changes <- df |>
    dplyr::mutate(
      x_prev = dplyr::lag(x),
      y_prev = dplyr::lag(y),
      x_next = dplyr::lead(x),
      y_next = dplyr::lead(y),

      dx = x - x_prev,
      dy = y - y_prev,

      sx = sgn(dx),
      sy = sgn(dy),

      sx_prev = dplyr::lag(sx),
      sy_prev = dplyr::lag(sy),

      changed = (sx != sx_prev | sy != sy_prev),

      # Exclude sticky events
      nonsticky =
        sx      != 0 & sy      != 0 &
        sx_prev != 0 & sy_prev != 0,

      # Outgoing direction
      vx = x_next - x,
      vy = y_next - y,
      vnorm = sqrt(vx^2 + vy^2),

      ux = vx / vnorm,
      uy = vy / vnorm,

      # Arrow endpoint: start at turning point, point in new direction
      x_arrow = x + arrow_len * ux,
      y_arrow = y + arrow_len * uy
    ) |>
    dplyr::filter(
      changed,
      nonsticky,
      !is.na(x_next),
      !is.na(y_next),
      is.finite(vnorm),
      vnorm > 0
    )

  list(
    geom_path(
      data = df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = .6,
      alpha = .85,
      color = colors[2]
    ),

    geom_segment(
      data = direction_changes,
      aes(
        x = x,
        y = y,
        xend = x_arrow,
        yend = y_arrow
      ),
      inherit.aes = FALSE,
      linewidth = .6,
      color = colors[2],
      lineend = "round",
      arrow = grid::arrow(
        ends = "last",
        type = "closed",
        length = grid::unit(head_len_mm, "mm")
      )
    ),

    geom_point(
      data = df[1, ],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = .7,
      color = "grey20",
      fill = "white"
    )
  )
}

"grey75"
"#EBEBEBFF"
p1_vs <- ggplot(df, aes(x, y, z = z)) +
  geom_hline(yintercept = 0, linewidth = 0.75, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.75, color = "grey80") +
  geom_contour(bins = 9, linewidth = .2, color = colors[1]) +
  add_trace(df_zz_vs, head_len_mm = 2.0, arrow_len = 0.05) +
  scale_x_continuous(breaks = seq(-3, 2, 1), limits = c(-3.5, 1.5), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(-2, 4, 2), limits = c(-2, 4),     expand = expansion(mult = 0)) +
  labs(x = NULL, y = NULL, title = "Sticky ZigZag") +
  coord_fixed() + theme_bw(base_size = 18) +
  theme(
    axis.ticks.length = unit(1, "mm"),
    axis.line = element_line(linewidth = .4, color = "grey30"),
    # plot.title = element_text(hjust = .5),
    plot.title=element_text(hjust = .1, margin=margin(t=40,b=-30)),
    plot.margin = margin(2, 2, 2, 2),
    panel.grid.major = element_line(color = "#EBEBEBFF", linewidth = .2),
    panel.grid.minor = element_blank()
  )
p1_vs
# p123
h <- 16
w <- 16
ggsave(file.path("/home/don/surfdrive/Postdoc/Psychometrika_PDMP_Variable_Selection/PDMPVariableSelection/figures", "dynamics_sticky.pdf"), p1_vs, height = h, width = w, units = "cm")


# estimators
base_plot <- ggplot(df, aes(x, y, z = z)) +
  geom_contour(linewidth = .2, bins = 15)

p1 <- base_plot + add_trace(df_zz)   + ggtitle("ZigZag")# + theme(title = element_text(hjust = .5))

trace_zz$skeleton[[1]]$event_positions
head(trace_zz)

trace_zz$skeleton[[1]]$initial_velocity
p0 <- trace_zz$skeleton[[1]]$initial_position
trace_zz$skeleton[[1]]$sparse
trace_zz$skeleton[[1]]$event_indices[1:3]
PDMPSamplersR:::
PDMPSamplersR:::

df_zz
direction_changes <- df_zz |>
  dplyr::mutate(
    x_prev = dplyr::lag(x),
    y_prev = dplyr::lag(y),
    x_next = dplyr::lead(x),
    y_next = dplyr::lead(y),

    dx = x - x_prev,
    dy = y - y_prev,

    sx = sgn(dx),
    sy = sgn(dy),

    sx_prev = dplyr::lag(sx),
    sy_prev = dplyr::lag(sy),

    changed = (sx != sx_prev | sy != sy_prev),

    # Exclude sticky events
    nonsticky =
      sx      != 0 & sy      != 0 &
      sx_prev != 0 & sy_prev != 0,

    # Outgoing direction
    vx = x_next - x,
    vy = y_next - y,
    vnorm = sqrt(vx^2 + vy^2),

    ux = vx / vnorm,
    uy = vy / vnorm,

    # Arrow endpoint: start at turning point, point in new direction
    x_arrow = x + arrow_len * ux,
    y_arrow = y + arrow_len * uy
  ) |>
  dplyr::filter(
    changed,
    nonsticky,
    !is.na(x_next),
    !is.na(y_next),
    is.finite(vnorm),
    vnorm > 0
  )



library(ggplot2)
library(dplyr)
library(patchwork)

library(ggplot2)
library(dplyr)
library(patchwork)

reconstruct_factorized_segments <- function(skel, t_max) {
  stopifnot(isTRUE(skel$sparse))
  stopifnot(length(skel$initial_position) == 2)
  stopifnot(length(skel$initial_velocity) == 2)
  stopifnot(t_max > skel$initial_time)

  pos <- skel$initial_position
  vel <- skel$initial_velocity
  t_cur <- skel$initial_time

  out <- vector("list", length(skel$event_times) + 1L)
  k <- 0L

  for (i in seq_along(skel$event_times)) {
    t_evt <- skel$event_times[i]

    if (t_evt > t_max)
      break

    dt <- t_evt - t_cur
    pos_evt <- pos + vel * dt

    # Sparse/factorized update: only one coordinate is stored
    j <- skel$event_indices[i]
    pos_evt[j] <- skel$event_positions[i]

    k <- k + 1L
    out[[k]] <- tibble(
      segment = k,
      t0 = t_cur,
      t1 = t_evt,
      dt = dt,
      x0 = pos[1],
      y0 = pos[2],
      x1 = pos_evt[1],
      y1 = pos_evt[2],
      vx = vel[1],
      vy = vel[2]
    )

    pos <- pos_evt
    vel[j] <- skel$event_velocities[i]
    t_cur <- t_evt
  }

  # Final partial segment up to t_max
  if (t_cur < t_max) {
    dt <- t_max - t_cur
    pos_end <- pos + vel * dt

    k <- k + 1L
    out[[k]] <- tibble(
      segment = k,
      t0 = t_cur,
      t1 = t_max,
      dt = dt,
      x0 = pos[1],
      y0 = pos[2],
      x1 = pos_end[1],
      y1 = pos_end[2],
      vx = vel[1],
      vy = vel[2]
    )
  }

  bind_rows(out[seq_len(k)])
}


make_marginal_layers <- function(segments, axis = c("x", "y"), t_max) {
  axis <- match.arg(axis)

  if (axis == "x") {
    lo <- pmin(segments$x0, segments$x1)
    hi <- pmax(segments$x0, segments$x1)
  } else {
    lo <- pmin(segments$y0, segments$y1)
    hi <- pmax(segments$y0, segments$y1)
  }

  # All points where the marginal density can change
  breaks <- sort(unique(c(lo, hi)))

  atoms <- tibble(
    start = head(breaks, -1),
    end   = tail(breaks, -1)
  ) |>
    mutate(mid = (start + end) / 2)

  # For ordinary ZigZag |v| = 1, so each segment contributes
  # exactly 1 / t_max density wherever it passes.
  out <- lapply(seq_len(nrow(atoms)), function(i) {
    a <- atoms$start[i]
    b <- atoms$end[i]
    m <- atoms$mid[i]

    active <- which(lo <= m & hi >= m)

    if (length(active) == 0)
      return(NULL)

    tibble(
      start = a,
      end = b,
      source_segment = active,
      layer = seq_along(active),
      height = seq_along(active) / t_max
    )
  })

  bind_rows(out)
}


make_marginal_envelope <- function(layers) {
  layers |>
    group_by(start, end) |>
    summarise(
      density = max(height),
      .groups = "drop"
    )
}


plot_zigzag_density_build_up <- function(
    skel,
    t_max,
    contour_df,
    xlim = c(-3.5, 1.5),
    ylim = c(-2, 4),
    bins = 9,
    base_size = 16,
    trace_color = colors[2],
    contour_color = colors[1]) {

  segments <- reconstruct_factorized_segments(
    skel = skel,
    t_max = t_max
  )

  mx <- make_marginal_layers(
    segments,
    axis = "x",
    t_max = t_max
  )

  my <- make_marginal_layers(
    segments,
    axis = "y",
    t_max = t_max
  )

  mx_env <- make_marginal_envelope(mx)
  my_env <- make_marginal_envelope(my)


  # Main trace panel
  p_main <- ggplot() +
    geom_hline(
      yintercept = 0,
      linewidth = 0.5,
      color = "#EBEBEBFF"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.5,
      color = "#EBEBEBFF"
    ) +
    geom_contour(
      data = contour_df,
      aes(x = x, y = y, z = z),
      bins = bins,
      linewidth = 0.3,
      color = contour_color
    ) +
    geom_segment(
      data = segments,
      aes(
        x = x0,
        y = y0,
        xend = x1,
        yend = y1
      ),
      linewidth = 0.6,
      alpha = 0.85,
      color = trace_color
    ) +
    geom_point(
      data = segments[1, ],
      aes(x = x0, y = y0),
      shape = 21,
      size = 3,
      stroke = 0.7,
      color = "grey20",
      fill = "white"
    ) +
    scale_x_continuous(
      limits = xlim,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = ylim,
      expand = expansion(mult = 0)
    ) +
    coord_fixed(clip = "off") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.grid.major = element_line(
        color = "#EBEBEBFF",
        linewidth = 0.2
      ),
      panel.grid.minor = element_blank()
    )


  # Top marginal:
  # each horizontal dashed line is one trajectory contribution.
  p_top <- ggplot() +
    geom_segment(
      data = mx,
      aes(
        x = start,
        xend = end,
        y = height,
        yend = height
      ),
      linewidth = 0.5,
      linetype = "dashed",
      color = trace_color
    ) +
    geom_segment(
      data = mx_env,
      aes(
        x = start,
        xend = end,
        y = density,
        yend = density
      ),
      linewidth = 0.4,
      color = "grey25"
    ) +
    scale_x_continuous(
      limits = xlim,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = NULL,
      y = "Density"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.margin = margin(2, 0, 2, 0),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )


  # Right marginal:
  # each vertical dashed line is one trajectory contribution.
  p_right <- ggplot() +
    geom_segment(
      data = my,
      aes(
        x = height,
        xend = height,
        y = start,
        yend = end
      ),
      linewidth = 0.5,
      linetype = "dashed",
      color = trace_color
    ) +
    geom_segment(
      data = my_env,
      aes(
        x = density,
        xend = density,
        y = start,
        yend = end
      ),
      linewidth = 0.4,
      color = "grey25"
    ) +
    scale_y_continuous(
      limits = ylim,
      expand = expansion(mult = 0)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Density",
      y = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.margin = margin(0, 2, 0, 2),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )


  # Compose
  (p_top + plot_spacer()) /
    (p_main + p_right) +
    plot_layout(
      widths = c(4, 1),
      heights = c(1, 4)
    )
}


# Example
plot_zigzag_density_build_up(
  skel = trace_zz$skeleton[[1]],
  t_max = 12,
  contour_df = df,
  xlim = c(-3.5, 1.5),
  ylim = c(-2, 4)
)
