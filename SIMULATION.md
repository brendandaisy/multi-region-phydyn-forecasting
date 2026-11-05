# Multi-region phylodynamic simulation suite

Forward simulation of a 4-region SIR epidemic in ReMASTER, with a case-reporting
and sequencing observation model applied in R, producing case data, a reporting
triangle, and a typed phylogeny with summary statistics.

## Layout

| file | role |
|---|---|
| `R/01-config.R` | every parameter: regions, rates, migration matrix, delays |
| `R/02-xml.R` | parametric ReMASTER XML generator + BEAST2 runner |
| `R/03-observe.R` | event extraction, reporting delay, triangle, truncation |
| `R/04-tree.R` | NEXUS parser, sequencing observation, 4 blocks of tree statistics |
| `R/05-plots.R` | figures |
| `R/06-schematic.R` | model diagram (`figures/00-model.png`), driven by the config |
| `run-simulation.R` | end-to-end driver |

Run with `Rscript run-simulation.R`. Outputs go to `output/` and `figures/`;
raw BEAST2 files go to `~/phylo-forecast-sim` (local disk, deliberately **not**
the iCloud-synced project folder).

## Model

Time in **weeks**, forward from week 0. Simulated to week 40; observed as of
week 30, so weeks 31-40 are held-out truth for forecast scoring.

Per region *i*, with `delta = mu + psi_case + psi_seq = 1.0` (1-week infectious period):

```
S[i] + I[i] -> 2I[i]     beta_i = R0_i * delta / N_i
I[i] -> 0                mu       = 0.60     unobserved recovery
I[i] -> case[i]          psi_case = 0.32     detected, not sequenced
I[i] -> seq[i]           psi_seq  = 0.08     detected and sequenced
I[i] -> I[j]             m_ij                migration
```

Detection probability is 0.4; 20% of detected cases are sequenced. Both
detections remove the individual, so `case` and `seq` are a thinning of one
detection process — sequences are a genuine subset of cases, which is what makes
the two data streams mutually informative.

`R0 = (A 1.3, B 1.5, C 1.2, D 1.6)`, `N = (A 10000, B 20000, C 800, D 1500)`,
outbreak seeded with one infection in **D**.

### Migration: three-parameter gravity model

```
m_ij = c * P_j^beta / d_ij^gamma        theta_mig = (c, beta, gamma)
```

with `c = 0.088`, `beta = 0.508`, `gamma = 1.5`. Destination population `P_j`
and distance `d_ij` are **data**, fixed outside the simulation loop, so the
parameter count does not grow with the number of regions.

There is deliberately no `P_i` term: ReMASTER's propensity for `I[i] -> I[j]`
is already `m_ij * I[i]`, so putting origin population in the rate would count
it twice. Asymmetry therefore comes from the single scalar `beta`:

```
m_ij / m_ji = (P_j / P_i)^beta
```

Distances (km, fictitious) were solved so that the model approximates the
qualitative matrix the project started from (`MIG_TARGET` in the config):

|        | A | B | C | D |
|---|---|---|---|---|
| **A** | 0 | 40 | 116 | 95 |
| **B** | 40 | 0 | 207 | 107 |
| **C** | 116 | 207 | 0 | 148 |
| **D** | 95 | 107 | 148 | 0 |

**Two things gravity cannot reproduce**, both inherent rather than fixable by
tuning:

1. *No single beta fits every pair.* The A-D asymmetry implies `beta = 1.21`,
   A-B implies `0.32`, and B-C implies `0.00`. The fitted `0.508` is a
   least-squares compromise, so D's total outflow lands at 0.024 per lineage
   per week against a target of 0.040.
2. *No structural zeros.* `d_ij` is finite for every pair, so every rate is
   strictly positive. The old matrix had `D -> C = 0`; gravity gives 0.0015,
   which seeds C directly from the index region and moves C's peak from week 32
   to week 18. If a hard zero matters, multiply the gravity matrix by a fixed
   adjacency mask -- still three parameters, with the mask as data.

## Observation model (all in R, none of it in ReMASTER)

**Cases.** Detection times are binned to weeks, then each case is assigned a
reporting delay drawn from a discretised Gamma(3, 1) (mean 3 weeks). Cases with
`detection week + delay > 30` are **not yet observed**. This truncation is part
of the simulator: training on complete counts and applying to a truncated
triangle would teach the network that incidence always falls at the right edge.

**Sequences.** Each sequenced tip draws a Gamma(4, 2) lab delay (mean 2 weeks);
tips with `sample time + delay > 30` are pruned. Pruning is the correct
treatment rather than a time-varying psi: detection already removed the
individual from the infectious pool whether or not the sequence has come back.
Only the observation is missing, not the epidemiological event.

Tip dates are rounded **up** to the end of their week, which only lengthens
terminal branches and so can never produce a negative branch length.

## Tree statistics (104 values)

Hard rule: every statistic is computable from a tree whose regions are known
only at the **tips**. ReMASTER annotates internal branches with the truth; that
is discarded and reconstructed by max parsimony, because real data has no such
annotation. True migration counts are extracted from the trajectory for
**validation only**.

1. **Global shape** — `treestats`, 51 statistics. Blind to regions by
   construction: permuting tip labels leaves every one unchanged.
2. **Per-region marginals** — tip counts, per-region TMRCA, per-region LTT.
3. **Association** — Fitch parsimony score, association index, maximum
   monophyletic clade size. Measures clustering, i.e. total migration intensity.
4. **Directionality** — parsimony-reconstructed i->j transition counts, pairwise
   asymmetry index, the population-ratio slope that identifies gravity's beta,
   and reconstruction-free cross-region MRCA quantiles.

## Timing (1215-tip tree, this machine)

| stage | time |
|---|---|
| BEAST2 / ReMASTER | 2.0s (JVM startup dominates) |
| trajectory parse | 0.4s |
| observation model | <0.1s |
| tree statistics | 0.6s |

`stats_global(fast = FALSE)` switches to `treestats::calc_all_stats`, which adds
~20 spectral statistics but costs **4.6s** instead of 0.02s. At SBI scale that
one call would dominate the entire budget.
