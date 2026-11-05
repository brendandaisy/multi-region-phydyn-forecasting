## =====================================================================
## 05-plots.R -- figures.
## =====================================================================
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})

REG_COL <- c(A = "#4C72B0", B = "#DD8452", C = "#55A868", D = "#C44E52")

theme_sim <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}
asof_layer <- function() list(
  geom_vline(xintercept = AS_OF, linetype = "dashed", colour = "grey30"),
  annotate("text", x = AS_OF - .4, y = Inf, label = "as-of", hjust = 1, vjust = 1.4,
           size = 3, colour = "grey30"))

p_prevalence <- function(prev) {
  ggplot(prev, aes(week, I, colour = region)) +
    geom_line(linewidth = .8) + asof_layer() +
    scale_colour_manual(values = REG_COL) +
    facet_wrap(~region, scales = "free_y", ncol = 4) +
    labs(title = "Latent prevalence (truth)", subtitle = "infectious individuals",
         x = "week", y = "I(t)") +
    theme_sim() + theme(legend.position = "none")
}

p_cases <- function(cases_cmp) {
  ggplot(cases_cmp, aes(week, n, colour = series, linetype = series)) +
    geom_line(linewidth = .75) + asof_layer() +
    scale_colour_manual(values = c(`true detections` = "grey25",
                                   `reported by as-of` = "#C44E52")) +
    scale_linetype_manual(values = c(`true detections` = 1, `reported by as-of` = 1)) +
    facet_wrap(~region, scales = "free_y", ncol = 4) +
    labs(title = "True vs reported cases",
         subtitle = "the gap at the right edge is the reporting triangle, not a real decline",
         x = "week", y = "cases") +
    theme_sim() + theme(legend.position = "bottom", legend.title = element_blank())
}

p_triangle <- function(tri_long) {
  ggplot(tri_long, aes(delay, week, fill = n)) +
    geom_tile(colour = "white", linewidth = .2) +
    scale_fill_viridis_c(option = "mako", direction = -1, na.value = "grey92") +
    scale_y_reverse(breaks = seq(0, AS_OF, 5)) +
    facet_wrap(~region, ncol = 4) +
    labs(title = "Reporting triangles", subtitle = "grey = not yet observed at the as-of date",
         x = "reporting delay (weeks)", y = "week of detection", fill = "cases") +
    theme_sim()
}

p_seq_dropoff <- function(seq_cmp) {
  ggplot(seq_cmp, aes(week, n, colour = series)) +
    geom_line(linewidth = .75) + asof_layer() +
    scale_colour_manual(values = c(`sampled` = "grey25", `available at as-of` = "#4C72B0")) +
    facet_wrap(~region, scales = "free_y", ncol = 4) +
    labs(title = "Sequencing drop-off near the present",
         subtitle = "sequences exist but have not come back from the lab",
         x = "week", y = "sequences") +
    theme_sim() + theme(legend.position = "bottom", legend.title = element_blank())
}

p_ltt <- function(ltt_long) {
  ggplot(ltt_long, aes(week, lineages, colour = region)) +
    geom_step(linewidth = .8) +
    scale_colour_manual(values = REG_COL) +
    labs(title = "Lineages through time by reconstructed region",
         subtitle = "from the sequenced tree as observed at the as-of date",
         x = "week", y = "lineages") +
    theme_sim()
}

p_matrix <- function(m, title, sub, lab_fmt = "%d") {
  d <- as.data.frame.table(m, responseName = "v")
  names(d)[1:2] <- c("from", "to")
  ggplot(d, aes(to, from, fill = v)) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_text(aes(label = sprintf(lab_fmt, v)), size = 3.2,
              colour = ifelse(d$v > max(d$v, na.rm = TRUE) * .55, "white", "grey15")) +
    scale_fill_viridis_c(option = "rocket", direction = -1) +
    scale_y_discrete(limits = rev(levels(d$from))) +
    labs(title = title, subtitle = sub, x = "to", y = "from", fill = NULL) +
    theme_sim() + theme(legend.position = "none")
}

p_xmrca <- function(xm) {
  ggplot(xm, aes(mrca, colour = pair)) +
    stat_ecdf(linewidth = .8) +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = "Cross-region MRCA times",
         subtitle = "reconstruction-free; a curve that rises late means recent exchange",
         x = "MRCA time (week)", y = "ECDF of tip pairs", colour = "pair") +
    theme_sim()
}
