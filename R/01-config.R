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

## ---- migration: three-parameter gravity model ----------------------
## Per-lineage rate per week from region i to region j:
##
##     m_ij = c * P_j^beta / d_ij^gamma
##
## P_j (destination population) and d_ij (distance) are DATA. Note there is
## no P_i term: ReMASTER's propensity for I[i] -> I[j] is already m_ij * I[i],
## so the origin's size enters through the lineage count, not the rate.
##
## Asymmetry is therefore controlled by the single scalar beta:
##     m_ij / m_ji = (P_j / P_i)^beta
##
## Distances between the four (fictitious) regions, km. Chosen so the model
## approximates the qualitative rates this project started from; see
## MIG_TARGET below and SIMULATION.md for how well it does.
DIST <- matrix(c(  0,  40, 116,  95,
                  40,   0, 207, 107,
                 116, 207,   0, 148,
                  95, 107, 148,   0),
               NREG, NREG, byrow = TRUE, dimnames = list(REGIONS, REGIONS))

GRAVITY <- list(c = 0.088, beta = 0.508, gamma = 1.5)

make_migration_gravity <- function(pop, dist, c, beta, gamma) {
  m <- c * outer(rep(1, length(pop)), pop^beta) / dist^gamma
  diag(m) <- 0
  dimnames(m) <- list(names(pop), names(pop))
  m
}

MIG <- make_migration_gravity(POP, DIST, GRAVITY$c, GRAVITY$beta, GRAVITY$gamma)

## The qualitative matrix this project used before gravity, kept for
## comparison only. It is NOT internally consistent with a single beta:
## the A-D pair implies beta = 1.21 while B-C implies beta = 0, so no
## three-parameter gravity model can reproduce all of it at once.
LEVELS <- c(none = 0, very_low = 0.002, low = 0.008,
            medium = 0.02, high = 0.04, very_high = 0.05)
MIG_TARGET <- matrix(c(    0, 0.050, 0.002, 0.002,
                       0.040,     0, 0.002, 0.002,
                       0.008, 0.002,     0, 0.002,
                       0.020, 0.020,     0,     0),
                     NREG, NREG, byrow = TRUE, dimnames = list(REGIONS, REGIONS))

## ---- observation model ----------------------------------------------
## Case reporting delay: gamma, mean 3 weeks.
REPORT_SHAPE <- 3; REPORT_RATE <- 1        # mean 3, sd 1.73
MAX_DELAY    <- 12                          # weeks of triangle to track

## Sequencing delay: gamma, mean 2 weeks -> sequencing coverage falls off
## over roughly the 5 weeks before AS_OF.
SEQ_SHAPE <- 4; SEQ_RATE <- 2               # mean 2, sd 1

SEED <- 20260301
