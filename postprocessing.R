library(dplyr)
library(purrr)
library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(jsonlite)
library(magrittr)
library(xml2)
library(lubridate)

library(timtamslamR)

set.seed(1)

## Define some hex codes to use as a colour scheme. These are from the
## ColorBrewer2 project: \code{https://colorbrewer2.org/}

my_colours <- list(
  truth = list(dark = "#d73027", light = "#fee090"),
  estimate = list(dark = "#4575b4", light = "#e0f3f8"),
  other = list(
    grey = list(
      light = "#e0e0e0",
      medium = "#999999",
      dark = "#4d4d4d"
    )
  )
)

## Add a helper function that is helpful for smoothing out the
## lubridate interface.

#' @title Convert numeric days to a lubridate period
#'
#' @param x A numeric vector of days.
#'
#' @return A lubridate period object.
#'
days_to_period <- function(x) {
  seconds(round(x * 24 * 60 * 60))
}

## Define a list of output files so they are all in one place.

output <- list(
  params_png = "tutorial-3-marginals.png",
  prev_png = "tutorial-3-prevalence.png"
)

## Define and extract relevant information for the input files after
## checking that they exist.

sim_xml <- "epidemic-simulation.xml"
stopifnot(file.exists(sim_xml))
sim_config <- read_xml(sim_xml)
timtam_xml <- "timtam-with-varying-obs.xml"
stopifnot(file.exists(timtam_xml))
mcmc_config <- read_xml(timtam_xml)

timtam_log <-
  mcmc_config |>
  xml_find_first("//logger[@id='tracelog']") |>
  xml_attr("fileName")
## If there is a \code{$(filebase)} in the log file name, replace it
## with the the base name of the XML file specifying the analysis.
timtam_log <- str_replace(timtam_log,
                          "\\$\\(filebase\\)",
                          str_remove(timtam_xml, "\\.xml"))
stopifnot(file.exists(timtam_log))
timtam_tree_id <-
  mcmc_config |>
  xml_find_first("//data[@id]/child::sequence/..") |>
  xml_attr("id")
timtam_tree_log <-
  mcmc_config |>
  xml_find_first("//logger[@mode='tree']") |>
  xml_attr("fileName")
timtam_tree_log <- str_replace(timtam_tree_log,
                               "\\$\\(filebase\\)",
                               str_remove(timtam_xml, "\\.xml"))
timtam_tree_log <- str_replace(timtam_tree_log,
                               "\\$\\(tree\\)",
                               timtam_tree_id)
stopifnot(file.exists(timtam_tree_log))

sim_traj <-
  sim_config |>
  xml_find_first("//logger[not(@mode)]") |>
  xml_attr("fileName")

## Because it is non-trivial to extract these parameters from the
## current simulation XML (due to the way the XML is structured), we
## will just hard-code them here.

true_inf_rate <- 0.0001
true_death_rate <- 0.1
true_occ_rate <- 0.08
true_sampling_rate <- 0.02
net_removal_rate <- true_death_rate + true_occ_rate + true_sampling_rate
r0_change_times <- c(14.0, 7.0)
num_r0_changes <- length(r0_change_times)
sim_duration <- 40.0

sim_surv_change_time <-
  sim_config |>
  xml_find_first(xpath = "//reaction[@changeTimes]") |>
  xml_attr("changeTimes") |>
  as.numeric()
sim_datetime_0 <- ymd_hms("2023-08-01 00:00:00")
sim_surv_start <- sim_datetime_0 + days_to_period(sim_surv_change_time)
sim_datetime_n <- sim_datetime_0 + days_to_period(sim_duration)

num_to_burn <- 200
timtam_post <- timtam_log %>%
  read_beast2_log(burn = num_to_burn) %>%
  mutate(r0_1 = TTR0.1, r0_2 = TTR0.2, r0_3 = TTR0.3)

## START R_EFF SUMMARY
##
## Since the effective reproduction number of this process depends on
## the realisation of the simulation, computing it is a bit fiddly, so
## feel free to skip this part of the code on a first reading. All you
## need to know is that we are building a dataframe that contains the
## effective reproduction number through time based on the simulation.

make_record <- function(x) {
  x |>
    as.numeric() |>
    as.list() |>
    set_names(c("t", "R", "S", "psi", "I", "omega")) |>
    as.data.frame()
}

sim_data <- sim_traj |>
  readLines() |>
  pluck(2) |>
  str_split("\t") |>
  unlist() |>
  pluck(2) |>
  str_split(";") |>
  unlist() |>
  str_replace_all("(t|R|S|psi|I|omega)=", "") |>
  str_split(":") |>
  map(make_record) |>
  bind_rows()

sim_r_eff <- data.frame(
  time = sim_data$t,
  r_eff = sim_data$S * true_inf_rate / net_removal_rate
)

## START STEP FUNCTION SUMMARY
##
## Summarising the step function is a bit fiddly, so feel free to skip
## this part of the code on a first reading. All you need to know is
## that we are building a dataframe that contains a summary of the
## estimated values of the basic reproduction number through time
## based on the parameters of the model.

#' Compute a posterior summary of the samples of a piece-wise constant
#' function.
#'
#' @param ts the times at which there is a step.
#' @param vss the posterior samples of the value of the function at
#'   the step.
#'
step_function_cri <- function(ts, vss) {
  .ts <- c(head(ts, 1), rep(tail(head(ts, -1), -1), each = 2), tail(ts, 1))
  n <- length(.ts)
  result <- data.frame(
    t = .ts,
    lower_bound = rep(NA, n),
    median = rep(NA, n),
    upper_bound = rep(NA, n)
  )
  for (ix in seq.int(length(vss))) {
    summ <- as.numeric(quantile(x = vss[[ix]], probs = c(0.025, 0.5, 0.975)))
    .ix <- 2 * ix - 1
    result[.ix, 2:4] <- summ
    result[.ix + 1, 2:4] <- summ
  }
  return(result)
}

sf_times <- sim_duration - c(sim_duration, r0_change_times, 0.0)
eval(parse(text = paste0(c("sf_samples <- list(", paste0(sprintf("timtam_post$r0_%d", 1:(num_r0_changes + 1)), collapse = ","), ")"), collapse = "")))
estimate_cri <- step_function_cri(sf_times, sf_samples)

#' END STEP FUNCTION SUMMARY

gg_r_eff <-
  ggplot() +
  geom_ribbon(
    data = estimate_cri,
    mapping = aes(x = t, ymin = lower_bound, ymax = upper_bound),
    fill = my_colours$estimate$light
  ) +
  geom_line(
    data = estimate_cri,
    mapping = aes(x = t, y = median),
    colour = my_colours$estimate$dark
  ) +
  geom_line(
    data = sim_r_eff,
    mapping = aes(x = time, y = r_eff),
    colour = my_colours$truth$dark
  ) +
  geom_vline(
    xintercept = c(sim_surv_change_time, sim_duration)
  ) +
  labs(
    x = "Days since outbreak began",
    y = "Effective reproduction number",
    title = "Estimated effective reproduction number",
    subtitle = "The values from the SIR model and piece-wise constant estimate"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13)
  )

ggsave(
  filename = output$params_png,
  plot = gg_r_eff,
  height = 10.5, width = 14.8,
  units = "cm"
)

my_present <- readRDS("last-sequence.Rds")
my_units <- "days"


my_beast2_log <- read_beast2_log(timtam_log, burn = num_to_burn)
my_beast2_trees <- read_beast2_trees(timtam_tree_log, burn = num_to_burn)
my_beast2_model <- read_beast2_xml(timtam_xml, my_present, my_units)

## Prevalence calculation can take a while so we will subsample the
## posterior samples first to keep things managable.

log_len <- nrow(my_beast2_log)
subsample_ixs <-
  round(seq(from = 1, to = nrow(my_beast2_log), length.out = 200))
prev_helper <- total_prevalence_log(
  my_beast2_log[subsample_ixs,],
  my_beast2_trees[subsample_ixs],
  my_beast2_model
)
prev_times  <- my_beast2_model$hist_times

## We can make a diagnostic plot of the prevalence through time
## estimates to double check that the MCMC has converged.

trace_plot <-
  ggplot() +
  geom_line(
    data = melt(prev_helper,
                id.vars = "Sample",
                value.name = "Prevalence",
                variable.name = "Time"),
    mapping = aes(x = Sample, y = Prevalence, colour = Time)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

## The following is just some boring munging to get the prevalence
## estimates into a sensible format for plotting.

bkwd_hist_times <- as.numeric(my_beast2_model$present$date_time - my_beast2_model$hist_times) / (60 * 60 * 24)

time_df <- data.frame(
  point = sprintf("size%d", seq_along(prev_times)),
  time = sim_duration - bkwd_hist_times
)

summary_vec <- function(prev_samples) {
  as.numeric(quantile(prev_samples, probs = c(0.025, 0.5, 0.975)))
}

## Use some meta-programming to get the prevalence estimates into a
## suitable data frame because I am lazy.

prev_df <- data.frame(type = c("lower", "mid", "upper"))
for (ix in seq_along(prev_times)) {
  eval(parse(text = paste0(c("prev_df$size", ix, " <- summary_vec(prev_helper$TTPrevalence.", ix, ")"), collapse = "")))
}

prev_df <-
  prev_df |>
  melt(id.vars = "type",
       variable.name = "point",
       value.name = "size") |>
  mutate(point = as.character(point)) |>
  full_join(time_df, by = "point") |>
  select(time, type, size) |>
  dcast(time ~ type, value.var = c("size"))
sim_prev_df <-
  data.frame(
    time = as.numeric(sim_data$t),
    num_infectious = as.numeric(sim_data$I)
  )

## Make a nice plot that demonstrates all of the prevalence results.

prev_fig <-
  ggplot() +
  geom_linerange(
    data = prev_df,
    mapping = aes(x = time, ymin = lower, ymax = upper),
    colour = my_colours$estimate$light,
    linewidth = 3
  ) +
  geom_point(
    data = prev_df,
    mapping = aes(x = time, y = mid),
    colour = my_colours$estimate$dark,
    size = 2
  ) +
  geom_line(
    data = sim_prev_df,
    mapping = aes(x = time, y = num_infectious),
    colour = my_colours$truth$dark
  ) +
  geom_vline(
    xintercept = c(sim_surv_change_time, sim_duration)
  ) +
  labs(
    x = "Days since outbreak began",
    y = "Prevalence of infection",
    title = "Prevalence of infection",
    subtitle = "Estimates and true values from simulation"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13)
  )

ggsave(
  filename = output$prev_png,
  plot = prev_fig,
  height = 10.5, width = 14.8,
  units = "cm"
)
