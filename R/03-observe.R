## =====================================================================
## 03-observe.R -- the observation model.
##
## ReMASTER gives the latent truth; everything here is applied in R:
##   (1) aggregate events to weeks
##   (2) convolve detections with a reporting delay -> reporting triangle
##   (3) truncate at AS_OF, so recent weeks are incomplete exactly as in
##       real data. This truncation IS part of the simulator.
## =====================================================================

suppressMessages({library(dplyr); library(tidyr)})

read_traj <- function(path) {
  read.delim(path) |>
    mutate(region = REGIONS[index + 1L])
}

#' Event times, recovered by differencing the step-function trajectory.
#' Infections are read off S (S only ever decreases, by transmission);
#' detections are read off the cumulative sample populations.
extract_events <- function(traj) {
  traj |>
    group_by(population, region) |>
    arrange(t, .by_group = TRUE) |>
    mutate(d = value - lag(value)) |>
    ungroup() |>
    filter(!is.na(d), d != 0) |>
    mutate(n = abs(d),
           event = case_when(population == "S"    & d < 0 ~ "infection",
                             population == "case" & d > 0 ~ "case",
                             population == "seq"  & d > 0 ~ "seq",
                             TRUE ~ NA_character_)) |>
    filter(!is.na(event)) |>
    select(event, region, t, n) |>
    uncount(n)
}

#' Prevalence on a weekly grid (for plotting the latent truth).
weekly_prevalence <- function(traj, weeks = 0:WEEKS_TOTAL) {
  traj |>
    filter(population == "I") |>
    group_by(region) |>
    reframe(week = weeks,
            I = approx(t, value, xout = weeks, method = "constant", rule = 2)$y)
}

#' Weekly counts of an event type, on a complete region x week grid.
weekly_counts <- function(events, type, weeks = 0:WEEKS_TOTAL) {
  events |>
    filter(event == type) |>
    mutate(week = pmin(floor(t), max(weeks))) |>
    count(region, week, name = "n") |>
    complete(region = REGIONS, week = weeks, fill = list(n = 0))
}

## ---------------------------------------------------------------------
## Reporting delay
## ---------------------------------------------------------------------

#' Discretised gamma delay pmf over 0..max_delay weeks.
delay_pmf <- function(shape = REPORT_SHAPE, rate = REPORT_RATE,
                      max_delay = MAX_DELAY) {
  p <- pgamma(1:(max_delay + 1), shape, rate) - pgamma(0:max_delay, shape, rate)
  p / sum(p)
}

#' Assign each detected case a reporting delay (stochastic, not a
#' deterministic convolution -- this is a simulator, the noise is real).
assign_report_delay <- function(cases, pmf = delay_pmf()) {
  cases |>
    mutate(delay    = sample.int(length(pmf), n(), replace = TRUE, prob = pmf) - 1L,
           rep_week = week + delay)
}

#' Reporting triangle: rows = detection week, cols = delay.
#' Cells with week + delay > as_of are NA: not yet observed.
reporting_triangle <- function(cases_delayed, region_sel, as_of = AS_OF,
                               weeks = 0:as_of, max_delay = MAX_DELAY) {
  d <- cases_delayed |> filter(region == region_sel, week <= as_of)
  m <- matrix(0L, length(weeks), max_delay + 1,
              dimnames = list(week = weeks, delay = 0:max_delay))
  tb <- table(factor(d$week, weeks), factor(d$delay, 0:max_delay))
  m[] <- as.integer(tb)
  for (i in seq_along(weeks)) {
    unobs <- (weeks[i] + 0:max_delay) > as_of
    m[i, unobs] <- NA_integer_
  }
  m
}

#' What a user sees on the as-of date: cumulative reported-so-far by week.
reported_as_of <- function(cases_delayed, as_of = AS_OF, weeks = 0:AS_OF) {
  cases_delayed |>
    filter(rep_week <= as_of, week <= as_of) |>
    count(region, week, name = "reported") |>
    complete(region = REGIONS, week = weeks, fill = list(reported = 0))
}

## ---------------------------------------------------------------------
## True migration events, recovered from the trajectory.
## A migration is the only event that moves one unit of I between two
## regions at a single instant: I[i] -1 and I[j] +1, nothing else.
## Used ONLY for validation -- never as a summary statistic, since it is
## not computable from real data.
## ---------------------------------------------------------------------
extract_migrations <- function(traj) {
  d <- traj |>
    filter(population == "I") |>
    group_by(region) |> arrange(t, .by_group = TRUE) |>
    mutate(d = value - lag(value)) |> ungroup() |>
    filter(!is.na(d), d != 0)

  d |>
    group_by(t) |>
    filter(n() == 2, sum(d) == 0, any(d == 1), any(d == -1)) |>
    summarise(from = region[d == -1], to = region[d == 1], .groups = "drop")
}

true_migration_matrix <- function(mig_events) {
  m <- matrix(0L, NREG, NREG, dimnames = list(from = REGIONS, to = REGIONS))
  tb <- table(factor(mig_events$from, REGIONS), factor(mig_events$to, REGIONS))
  m[] <- as.integer(tb); m
}
