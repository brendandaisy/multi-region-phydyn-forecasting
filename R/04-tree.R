## =====================================================================
## 04-tree.R -- read ReMASTER trees, apply the sequencing observation
## process, and compute summary statistics.
##
## HARD RULE: every statistic here must be computable from a tree whose
## regions are known only at the TIPS. ReMASTER annotates internal
## branches with the truth; we deliberately throw that away and
## reconstruct, because real data has no such annotation.
## =====================================================================

suppressMessages({library(ape); library(dplyr); library(tidyr)})

## ---------------------------------------------------------------------
## Parsing. ape::read.nexus silently DROPS the [&type=...] annotations,
## which is the whole signal, so we parse the NEXUS ourselves.
## ---------------------------------------------------------------------
read_remaster_tree <- function(path) {
  txt  <- readLines(path, warn = FALSE)
  line <- grep("^\\s*tree ", txt, value = TRUE)[1]
  if (is.na(line)) stop("no tree line found in ", path)

  ## tips carry a samp= attribute; internal nodes do not
  tip_re <- '([0-9]+)\\[&type="I\\{([0-9]+)\\}",samp="[A-Za-z]+\\{([0-9]+)\\}",time=([0-9.eE+-]+)\\]'
  m <- stringr::str_match_all(line, tip_re)[[1]]
  if (nrow(m) == 0) stop("no annotated tips matched in ", path)

  tips <- tibble(label  = m[, 2],
                 region = REGIONS[as.integer(m[, 3]) + 1L],
                 time   = as.numeric(m[, 5]))

  newick <- gsub("\\[&[^]]*\\]", "", sub("^[^(]*", "", line))
  phy <- read.tree(text = newick)
  if (is.null(phy)) stop("failed to parse newick from ", path)

  tips <- tips[match(phy$tip.label, tips$label), ]
  stopifnot(!anyNA(tips$region))
  list(phy = phy, tips = tips)
}

## ---------------------------------------------------------------------
## Sequencing observation process (applied post hoc, in R).
##
## Pruning tips is the CORRECT way to do this rather than lowering psi
## near the present: detection already removed the individual from the
## infectious pool, whether or not the sequence has come back yet. Only
## the observation is missing, not the epidemiological event.
## ---------------------------------------------------------------------
apply_seq_observation <- function(tr, as_of = AS_OF,
                                  shape = SEQ_SHAPE, rate = SEQ_RATE) {
  tips <- tr$tips |>
    mutate(seq_delay = rgamma(n(), shape, rate),
           available = time + seq_delay <= as_of)
  keep <- tips$label[tips$available]
  if (length(keep) < 3) stop("fewer than 3 sequences available at as_of")
  phy <- keep.tip(tr$phy, keep)
  list(phy = phy, tips = tips[match(phy$tip.label, tips$label), ], all_tips = tips)
}

#' Tip dates are only known to the week, so push each tip to its week end.
#' Rounding UP only lengthens terminal branches, so this can never produce
#' a negative branch length.
round_tips_to_week <- function(phy, times) {
  target <- ceiling(times)
  term <- which(phy$edge[, 2] <= Ntip(phy))
  tip_of_edge <- phy$edge[term, 2]
  phy$edge.length[term] <- phy$edge.length[term] + (target[tip_of_edge] - times[tip_of_edge])
  phy
}

## ---------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------

#' Per-node tip-state counts, via a postorder sweep.
node_state_counts <- function(phy, states) {
  n <- Ntip(phy); K <- nlevels(states)
  cnt <- matrix(0L, n + phy$Nnode, K, dimnames = list(NULL, levels(states)))
  cnt[cbind(seq_len(n), as.integer(states))] <- 1L
  ph <- reorder(phy, "postorder")
  for (e in seq_len(nrow(ph$edge)))
    cnt[ph$edge[e, 1], ] <- cnt[ph$edge[e, 1], ] + cnt[ph$edge[e, 2], ]
  cnt
}

#' Forward time of every node (tips and internal), given tip times.
node_times <- function(phy, tip_times) {
  depth <- node.depth.edgelength(phy)          # root = 0
  root_time <- mean(tip_times - depth[seq_along(tip_times)])
  depth + root_time
}

#' Max-parsimony ancestral states, then the state of every branch
#' (taken as the state of its child node).
branch_states <- function(phy, states) {
  K <- nlevels(states)
  asr <- castor::asr_max_parsimony(phy, as.integer(states), Nstates = K)
  node_best <- max.col(asr$ancestral_likelihoods, ties.method = "first")
  all_state <- c(as.integer(states), node_best)
  list(node = all_state,
       edge = all_state[phy$edge[, 2]],
       parent = all_state[phy$edge[, 1]])
}

## ---------------------------------------------------------------------
## BLOCK 1 -- global tree shape (treestats). Blind to regions by
## construction: permuting tip labels leaves every one of these unchanged.
## ---------------------------------------------------------------------
## `fast = TRUE` uses the topology + branching-time groups (51 stats, ~0.02s).
## `calc_all_stats` adds ~20 spectral/distance statistics but costs ~4.6s on a
## 1200-tip tree, which would dominate an SBI simulation budget.
stats_global <- function(phy, fast = TRUE) {
  s <- if (fast) c(treestats::calc_topology_stats(phy), treestats::calc_brts_stats(phy))
       else treestats::calc_all_stats(phy, normalize = FALSE)
  s <- s[!vapply(s, function(z) length(z) != 1 || !is.finite(z), logical(1))]
  setNames(as.numeric(s), paste0("global_", names(s)))
}

## ---------------------------------------------------------------------
## BLOCK 2 -- per-region marginals
## ---------------------------------------------------------------------
stats_region <- function(phy, states, tip_times, grid = seq(0, AS_OF, by = 1)) {
  bs <- branch_states(phy, states)
  nt <- node_times(phy, tip_times)
  e_start <- nt[phy$edge[, 1]]; e_end <- nt[phy$edge[, 2]]

  ltt <- vapply(seq_len(nlevels(states)), function(k) {
    vapply(grid, function(tt) sum(e_start <= tt & e_end > tt & bs$edge == k), numeric(1))
  }, numeric(length(grid)))
  colnames(ltt) <- levels(states)

  tipcount <- table(factor(states, levels(states)))
  tmrca <- vapply(levels(states), function(k) {
    w <- which(states == k)
    if (length(w) < 2) return(NA_real_)
    nt[getMRCA(phy, phy$tip.label[w])]
  }, numeric(1))

  list(ltt = ltt, grid = grid,
       vec = c(setNames(as.numeric(tipcount), paste0("ntip_", levels(states))),
               setNames(tmrca,                paste0("tmrca_", levels(states)))))
}

## ---------------------------------------------------------------------
## BLOCK 3 -- phylogeny/trait association (BaTS-style).
## Measures how strongly regions cluster: total migration intensity.
## ---------------------------------------------------------------------
stats_association <- function(phy, states) {
  cnt <- node_state_counts(phy, states)
  internal <- (Ntip(phy) + 1):(Ntip(phy) + phy$Nnode)
  ci <- cnt[internal, , drop = FALSE]
  tot <- rowSums(ci); fmax <- apply(ci, 1, max)

  ## Association index. 2^(n-1) overflows to Inf for big clades, which
  ## correctly sends those terms to 0.
  AI <- sum((1 - fmax / tot) / (2^(tot - 1) - 1), na.rm = TRUE)

  ## Maximum monophyletic clade size, overall and per state.
  mono <- fmax == tot
  MC <- max(tot[mono], 1)
  MC_k <- vapply(seq_len(ncol(ci)), function(k) {
    sel <- mono & ci[, k] == tot
    if (any(sel)) max(tot[sel]) else 1
  }, numeric(1))

  ## Fitch parsimony score for the region trait on this tree.
  pd <- phangorn::phyDat(matrix(as.character(states), ncol = 1,
                                dimnames = list(phy$tip.label, NULL)),
                         type = "USER", levels = levels(states))
  PS <- as.numeric(phangorn::parsimony(phy, pd))

  c(assoc_AI = AI, assoc_PS = PS, assoc_PS_norm = PS / Ntip(phy),
    assoc_MC = MC,
    setNames(MC_k, paste0("assoc_MC_", levels(states))))
}

## ---------------------------------------------------------------------
## BLOCK 4 -- directionality. The ONLY block that can see asymmetry.
## ---------------------------------------------------------------------
stats_directional <- function(phy, states, tip_times, pop = POP) {
  K <- nlevels(states); lv <- levels(states)
  bs <- branch_states(phy, states)

  ## (a) reconstructed transition counts
  trans <- matrix(0L, K, K, dimnames = list(from = lv, to = lv))
  chg <- which(bs$parent != bs$edge)
  if (length(chg))
    for (e in chg) trans[bs$parent[e], bs$edge[e]] <- trans[bs$parent[e], bs$edge[e]] + 1L

  ## (b) pairwise asymmetry index, and its relationship to population
  ##     ratio -- under gravity, m_ij/m_ji = (P_j/P_i)^beta, so this is
  ##     the slope that identifies beta.
  pairs <- t(combn(K, 2))
  asym <- apply(pairs, 1, function(p) {
    a <- trans[p[1], p[2]]; b <- trans[p[2], p[1]]
    if (a + b == 0) NA_real_ else (a - b) / (a + b)
  })
  logratio <- apply(pairs, 1, function(p) log(pop[lv[p[2]]] / pop[lv[p[1]]]))
  ok <- !is.na(asym)
  beta_slope <- if (sum(ok) >= 2) unname(coef(lm(asym[ok] ~ logratio[ok]))[2]) else NA_real_

  ## (c) cross-region MRCA times, reconstruction-free.
  ##     For tip times t_i, t_j and patristic distance d_ij,
  ##     MRCA time = (t_i + t_j - d_ij)/2.
  co <- cophenetic.phylo(phy)
  tt <- tip_times
  mr <- outer(tt, tt, "+") / 2 - co / 2
  gi <- as.integer(states)
  xr <- lapply(seq_len(nrow(pairs)), function(r) {
    i <- which(gi == pairs[r, 1]); j <- which(gi == pairs[r, 2])
    if (!length(i) || !length(j)) return(c(NA, NA, NA))
    v <- as.numeric(mr[i, j]); unname(quantile(v, c(.1, .5, .9), na.rm = TRUE))
  })
  xr <- do.call(rbind, xr)
  pnm <- paste0(lv[pairs[, 1]], lv[pairs[, 2]])

  list(trans = trans,
       vec = c(setNames(as.numeric(trans), paste0("trans_", outer(lv, lv, paste0))),
               setNames(asym, paste0("asym_", pnm)),
               beta_slope = beta_slope,
               setNames(as.numeric(xr), paste0("xmrca_", rep(c("q10","q50","q90"), each = nrow(pairs)), "_", pnm))))
}

## ---------------------------------------------------------------------
## Everything, assembled.
## ---------------------------------------------------------------------
summarise_tree <- function(phy, tips, fast = TRUE) {
  states <- factor(tips$region, levels = REGIONS)
  tt <- tips$time
  reg <- stats_region(phy, states, tt)
  dir <- stats_directional(phy, states, tt)
  list(
    vec = c(stats_global(phy, fast = fast), reg$vec, stats_association(phy, states), dir$vec),
    ltt = reg$ltt, grid = reg$grid, trans = dir$trans, states = states)
}
