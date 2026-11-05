## =====================================================================
## 06-schematic.R -- a diagram of the model itself.
## Everything is driven by the config, so editing rates updates the figure.
## =====================================================================
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

## node layout: A and B adjacent (strong exchange), C and D peripheral
## circle as a polygon, so we don't need ggforce
circle_df <- function(x0, y0, r, id, n = 72) {
  th <- seq(0, 2*pi, length.out = n)
  data.frame(x = x0 + r*cos(th), y = y0 + r*sin(th), id = id)
}

NODE_XY <- tibble::tribble(
  ~region, ~x,    ~y,
  "A",     0.00,  0.62,
  "B",     1.30,  0.62,
  "C",     0.00, -0.52,
  "D",     1.30, -0.52
)

p_regions <- function(mig = MIG, pop = POP, r0 = R0, seed_region = SEED_REGION) {

  nodes <- NODE_XY |>
    mutate(pop = pop[region], R0 = r0[region],
           rad = 0.085 + 0.095 * sqrt(pop / max(pop)),
           seed = region == seed_region)

  ## one row per directed flow, trimmed so arrows stop at the circle edge
  ed <- expand.grid(from = REGIONS, to = REGIONS, stringsAsFactors = FALSE) |>
    filter(from != to) |>
    mutate(rate = mapply(function(i, j) mig[i, j], from, to)) |>
    filter(rate > 0) |>
    left_join(nodes |> select(from = region, x0 = x, y0 = y, r0rad = rad), by = "from") |>
    left_join(nodes |> select(to   = region, x1 = x, y1 = y, r1rad = rad), by = "to") |>
    mutate(dx = x1 - x0, dy = y1 - y0, L = sqrt(dx^2 + dy^2),
           xs = x0 + dx / L * (r0rad + 0.02), ys = y0 + dy / L * (r0rad + 0.02),
           xe = x1 - dx / L * (r1rad + 0.055), ye = y1 - dy / L * (r1rad + 0.055))

  ggplot() +
    ## migration arrows: width and colour both encode the rate
    geom_curve(data = ed,
               aes(xs, ys, xend = xe, yend = ye, linewidth = rate,
                   colour = rate, alpha = rate),
               curvature = 0.22, angle = 90, ncp = 12,
               arrow = arrow(length = unit(0.016, "npc"), type = "closed"),
               lineend = "round") +
    scale_linewidth(range = c(0.3, 2.8), guide = "none") +
    scale_alpha(range = c(0.42, 1), guide = "none") +
    scale_colour_gradient(low = "#AEBECE", high = "#1F3B5B",
                          breaks = c(0.002, 0.02, 0.05),
                          labels = c("0.002\nvery low", "0.02\nmedium", "0.05\nvery high"),
                          name = "migration rate (per lineage per week)",
                          guide = guide_colourbar(title.position = "top", title.hjust = .5,
                                                  barwidth = unit(11, "lines"),
                                                  barheight = unit(0.5, "lines"))) +
    ## regions
    geom_polygon(data = do.call(rbind, Map(circle_df, nodes$x, nodes$y, nodes$rad, nodes$region)),
                 aes(x, y, group = id, fill = id), colour = "white", linewidth = 1.1) +
    geom_path(data = { sd <- filter(nodes, seed)
                       do.call(rbind, Map(circle_df, sd$x, sd$y, sd$rad + 0.035, sd$region)) },
              aes(x, y, group = id), colour = "#C44E52",
              linetype = "21", linewidth = .7, inherit.aes = FALSE) +
    scale_fill_manual(values = REG_COL, guide = "none") +
    geom_text(data = nodes, aes(x, y, label = region),
              colour = "white", fontface = "bold", size = 7) +
    geom_text(data = nodes,
              aes(x, y - rad - 0.075,
                  label = sprintf("N = %s   R₀ = %.1f", format(pop, big.mark = ",", trim = TRUE), R0)),
              size = 3.1, colour = "grey25") +
    geom_text(data = filter(nodes, seed), aes(x, y + rad + 0.135, label = "index case"),
              size = 3, colour = "#C44E52", fontface = "italic") +
    coord_fixed(clip = "off") +
    labs(title = "Four regions, asymmetric migration") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 12),
          legend.position = "bottom",
          legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 7.5, lineheight = .95),
          plot.margin = margin(6, 6, 6, 6))
}

p_compartments <- function() {
  bx <- tibble::tribble(
    ~id,    ~x,   ~y,   ~lab,            ~fill,
    "S",  0.00,  0.00, "S",             "#9FB4C7",
    "I",  1.00,  0.00, "I",             "#C44E52",
    "R",  2.10,  0.62, "recovered",     "#E8E8E8",
    "C",  2.10,  0.00, "case",          "#F0C9A0",
    "Q",  2.10, -0.62, "sequenced",     "#A8C8A0"
  )
  w <- 0.30; h <- 0.19
  bx <- bx |> mutate(w = ifelse(id %in% c("S","I"), 0.22, 0.46), h = h)

  seg <- tibble::tribble(
    ~x0,   ~y0,   ~x1,   ~y1,   ~lab,                ~lx,   ~ly,  ~col,
    0.11,  0.00,  0.88,  0.00, "beta[i]*S[i]*I[i]",  0.495, 0.085, "grey20",
    1.12,  0.06,  1.86,  0.56, "mu",                 1.44,  0.40,  "grey45",
    1.12,  0.00,  1.86,  0.00, "psi[case]",          1.49,  0.085, "grey20",
    1.12, -0.06,  1.86, -0.56, "psi[seq]",           1.44, -0.40,  "grey20"
  )

  ggplot() +
    geom_segment(data = seg, aes(x0, y0, xend = x1, yend = y1, colour = col),
                 arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
                 linewidth = .55, lineend = "round") +
    scale_colour_identity() +
    geom_text(data = seg, aes(lx, ly, label = lab), parse = TRUE,
              size = 3.1, colour = "grey20") +
    ## migration in/out of the infectious class
    geom_curve(aes(x = 0.80, y = 0.40, xend = 0.97, yend = 0.11),
               curvature = -0.35, linewidth = .5, colour = "#1F3B5B",
               arrow = arrow(length = unit(0.018, "npc"), type = "closed")) +
    geom_curve(aes(x = 1.03, y = 0.11, xend = 1.20, yend = 0.40),
               curvature = -0.35, linewidth = .5, colour = "#1F3B5B",
               arrow = arrow(length = unit(0.018, "npc"), type = "closed")) +
    annotate("text", x = 1.0, y = 0.50, label = "migration~(m[ij])", parse = TRUE,
             size = 3.1, colour = "#1F3B5B") +
    geom_tile(data = bx, aes(x, y, width = w, height = h, fill = fill),
              colour = "white", linewidth = 1) +
    scale_fill_identity() +
    geom_text(data = bx, aes(x, y, label = lab),
              size = ifelse(bx$id %in% c("S","I"), 4.6, 3.2),
              fontface = ifelse(bx$id %in% c("S","I"), "bold", "plain"),
              colour = ifelse(bx$id == "I", "white", "grey15")) +
    annotate("text", x = 1.05, y = -0.90, size = 3,
             colour = "grey35", hjust = .5, fontface = "italic",
             label = "detection removes the individual; sequenced tips are a 20% thinning of cases") +
    coord_fixed(clip = "off", xlim = c(-0.2, 2.45), ylim = c(-1.0, 0.62)) +
    labs(title = "Within each region") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 12),
          plot.margin = margin(6, 6, 6, 6))
}

p_model <- function() {
  (p_regions() | p_compartments()) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      title = "Multi-region birth-death model as simulated in ReMASTER",
      subtitle = "time in weeks · 1-week infectious period · 40% of infections detected, 20% of those sequenced",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(colour = "grey35", size = 10)))
}
