## =====================================================================
## run-simulation.R -- end-to-end driver.
##   ReMASTER truth -> observation model -> tree statistics -> figures
## =====================================================================
suppressMessages({library(glue); library(stringr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
source("R/01-config.R"); source("R/02-xml.R"); source("R/03-observe.R")
source("R/04-tree.R");   source("R/05-plots.R"); source("R/06-schematic.R")
set.seed(SEED)
.t <- Sys.time(); tick <- function(lab){ now <- Sys.time(); cat(sprintf("  [%-22s %5.1fs]\n", paste0(lab,"]"), as.numeric(difftime(now,.t,units="secs")))); .t <<- now }

## ---- 1. simulate -----------------------------------------------------
xml <- build_remaster_xml(BETA, MIG, POP,
  must_have = "case[0]>=20 &amp;&amp; case[1]>=20 &amp;&amp; case[2]>=3 &amp;&amp; case[3]>=20")
writeLines(xml, "output/mtbd-4region-sir.xml")
if (!file.exists(file.path(SCRATCH, "sim-4region.traj"))) {
  run_remaster(xml, "sim-4region", seed = SEED)
}
tick("simulate (BEAST2)")
traj <- read_traj(file.path(SCRATCH, "sim-4region.traj"))
tick("read trajectory")

## ---- 2. latent truth -------------------------------------------------
prev   <- weekly_prevalence(traj)
events <- extract_events(traj)
migs   <- extract_migrations(traj)
true_mig_counts <- true_migration_matrix(migs)

tick("latent truth")

## ---- 3. case observation model --------------------------------------
det <- events |>
  filter(event %in% c("case", "seq")) |>
  mutate(week = pmin(floor(t), WEEKS_TOTAL)) |>
  assign_report_delay()

true_cases <- det |> count(region, week, name = "n") |>
  complete(region = REGIONS, week = 0:WEEKS_TOTAL, fill = list(n = 0))
rep_cases  <- reported_as_of(det)

cases_cmp <- bind_rows(
  true_cases |> filter(week <= AS_OF) |> mutate(series = "true detections"),
  rep_cases  |> rename(n = reported)  |> mutate(series = "reported by as-of")) |>
  mutate(series = factor(series, c("true detections", "reported by as-of")))

tri <- lapply(REGIONS, function(r) {
  m <- reporting_triangle(det, r)
  as.data.frame.table(m, responseName = "n") |>
    mutate(region = r, week = as.integer(as.character(week)),
           delay = as.integer(as.character(delay)))
}) |> bind_rows()

tick("case observation")

## ---- 4. sequence observation model ----------------------------------
tr  <- read_remaster_tree(file.path(SCRATCH, "sim-4region.sampled.trees"))
obs <- apply_seq_observation(tr)
obs$phy <- round_tips_to_week(obs$phy, obs$tips$time)
obs$tips$time <- ceiling(obs$tips$time)

seq_cmp <- bind_rows(
  obs$all_tips |> mutate(week = floor(time)) |> count(region, week, name = "n") |>
    mutate(series = "sampled"),
  obs$all_tips |> filter(available) |> mutate(week = floor(time)) |>
    count(region, week, name = "n") |> mutate(series = "available at as-of")) |>
  complete(region = REGIONS, week = 0:WEEKS_TOTAL, series, fill = list(n = 0)) |>
  filter(week <= WEEKS_TOTAL)

tick("sequence observation")

## ---- 5. tree statistics ---------------------------------------------
S <- summarise_tree(obs$phy, obs$tips)

ltt_long <- as.data.frame(S$ltt) |> mutate(week = S$grid) |>
  pivot_longer(-week, names_to = "region", values_to = "lineages")

states <- factor(obs$tips$region, REGIONS)
nt  <- node_times(obs$phy, obs$tips$time)
co  <- cophenetic.phylo(obs$phy)
mrc <- outer(obs$tips$time, obs$tips$time, "+")/2 - co/2
prs <- t(combn(NREG, 2))
xm <- lapply(seq_len(nrow(prs)), function(r) {
  i <- which(as.integer(states) == prs[r,1]); j <- which(as.integer(states) == prs[r,2])
  if (!length(i) || !length(j)) return(NULL)
  data.frame(pair = paste0(REGIONS[prs[r,1]], "-", REGIONS[prs[r,2]]),
             mrca = as.numeric(mrc[i, j]))
}) |> bind_rows()

tick("tree statistics")

## ---- 6. outputs ------------------------------------------------------
write.csv(true_cases, "output/true-cases.csv", row.names = FALSE)
write.csv(rep_cases,  "output/reported-cases-as-of.csv", row.names = FALSE)
write.csv(tri,        "output/reporting-triangles.csv", row.names = FALSE)
write.csv(data.frame(statistic = names(S$vec), value = as.numeric(S$vec)),
          "output/tree-summary-statistics.csv", row.names = FALSE)
write.csv(as.data.frame.table(S$trans, responseName = "n"),
          "output/reconstructed-transitions.csv", row.names = FALSE)
write.csv(as.data.frame.table(true_mig_counts, responseName = "n"),
          "output/true-migration-counts.csv", row.names = FALSE)

pngs <- function(f, p, w, h) ggsave(file.path("figures", f), p, width = w, height = h, dpi = 150)
ggsave("figures/00-model.png", p_model(), width = 11.5, height = 5.2, dpi = 160, bg = "white")
pngs("01-prevalence.png",   p_prevalence(prev),      11, 2.9)
pngs("02-cases.png",        p_cases(cases_cmp),      11, 3.4)
pngs("03-triangles.png",    p_triangle(tri),         12, 3.6)
pngs("04-seq-dropoff.png",  p_seq_dropoff(seq_cmp),  11, 3.4)
pngs("05-ltt.png",          p_ltt(ltt_long),          7, 3.6)
pngs("07-xmrca.png",        p_xmrca(xm),              7, 3.6)

mig_rate_m <- MIG; diag(mig_rate_m) <- NA
pngs("06-migration.png",
  (p_matrix(mig_rate_m,       "True per-lineage rate", "model input (weeks^-1)", "%.3f") |
   p_matrix(true_mig_counts,  "True migration events", "from trajectory (unobservable)") |
   p_matrix(S$trans,          "Reconstructed", "max-parsimony on sequenced tree")) +
   plot_annotation(title = "Migration: truth vs what the tree recovers",
                   theme = theme(plot.title = element_text(face = "bold"))),
  13, 3.8)

png("figures/08-tree.png", width = 1500, height = 1500, res = 150)
bs <- branch_states(obs$phy, states)
plot(obs$phy, type = "fan", show.tip.label = FALSE, no.margin = TRUE,
     edge.color = REG_COL[REGIONS[bs$edge]], edge.width = 1.1)
legend("topleft", legend = REGIONS, col = REG_COL[REGIONS], lwd = 3, bty = "n", title = "region")
dev.off()

saveRDS(list(S = S, obs = obs, det = det, traj_summary = prev,
             true_mig = true_mig_counts), file.path(SCRATCH, "analysis.rds"))

cat("\n=== DONE ===\n")
tick("figures")
cat("tips analysed:", Ntip(obs$phy), " summary stats:", length(S$vec), "\n")
