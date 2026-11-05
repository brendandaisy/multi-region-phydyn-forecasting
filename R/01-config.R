## =====================================================================
## 01-config.R -- all simulation settings in one place.
##
## Time unit: WEEKS throughout. Simulation runs forward from week 0.
## =====================================================================

SCRATCH <- path.expand("~/phylo-forecast-sim")   # BEAST runs: local disk, NOT iCloud
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)

BEAST <- "/Applications/BEAST 2.7.7/bin/beast"

## ---- time axis ------------------------------------------------------
WEEKS_TOTAL <- 40    # simulate this much truth
AS_OF       <- 30    # observation cutoff ("today"); weeks 31-40 are held out

## ---- regions --------------------------------------------------------
REGIONS <- c("A", "B", "C", "D")
NREG    <- length(REGIONS)
POP     <- c(A = 10000, B = 20000, C = 800, D = 1500)
SEED_REGION <- "D"

## ---- epidemiology ---------------------------------------------------
## Flu-like transmission with a deliberately high removal rate.
## delta = mu + psi_case + psi_seq = total becoming-uninfectious rate.
## Mean infectious period = 1/delta = 1 week.
DELTA     <- 1.0
MU        <- 0.6     # unobserved recovery (constant across regions)
PSI_TOTAL <- 0.4     # detection-with-removal; p_detect = 0.4
P_SEQ     <- 0.20    # fraction of detected cases that get sequenced
PSI_CASE  <- PSI_TOTAL * (1 - P_SEQ)   # 0.32  detected, not sequenced
PSI_SEQ   <- PSI_TOTAL * P_SEQ         # 0.08  detected and sequenced

## R0 per region (= lambda/delta since delta = 1)
R0 <- c(A = 1.3, B = 1.5, C = 1.2, D = 1.6)

## ReMASTER's S+I->2I propensity is rate*S*I, so the XML rate is
## beta_i = R0_i * delta / N_i  (frequency-dependent transmission).
BETA <- R0 * DELTA / POP

## ---- migration ------------------------------------------------------
## Per-lineage rate per week, row = origin, col = destination.
## Taken as given ("migration data is available") rather than from gravity.
LEVELS <- c(none = 0, very_low = 0.002, low = 0.008,
            medium = 0.02, high = 0.04, very_high = 0.05)

MIG_LEVEL <- rbind(
  A = c(A = NA,          B = "very_high", C = "very_low", D = "very_low"),
  B = c(A = "high",      B = NA,          C = "very_low", D = "very_low"),
  C = c(A = "low",       B = "very_low",  C = NA,         D = "very_low"),
  D = c(A = "medium",    B = "medium",    C = "none",     D = NA)
)

MIG <- matrix(0, NREG, NREG, dimnames = list(REGIONS, REGIONS))
for (i in REGIONS) for (j in REGIONS)
  if (i != j) MIG[i, j] <- LEVELS[[ MIG_LEVEL[i, j] ]]

## ---- observation model ----------------------------------------------
## Case reporting delay: gamma, mean 3 weeks.
REPORT_SHAPE <- 3; REPORT_RATE <- 1        # mean 3, sd 1.73
MAX_DELAY    <- 12                          # weeks of triangle to track

## Sequencing delay: gamma, mean 2 weeks -> sequencing coverage falls off
## over roughly the 5 weeks before AS_OF.
SEQ_SHAPE <- 4; SEQ_RATE <- 2               # mean 2, sd 1

SEED <- 20260301
