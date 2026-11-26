### Benchmark script: when is parallel computation of the simulations worth it? (graft_missing_taxa(), `parallel` argument)

## Design (same number of taxa to graft as tips in the backbone: 100, 500 and 1000 tips; 10, 100 and 1000 simulations; both methods, no simmap)
##  Taxa to graft: the taxa are spread over `n_groups = ceiling(log(n_tips))` target clades (disjoint clades of the backbone, 5 tips or more (the default `min_tips` of the birth-death fits), at most half of the tips),
##          with a broken stick distribution of the number of taxa per clade (a few large groups and many small ones, as in a taxonomy; the largest groups go to the largest clades, which are more likely to be chosen). Taxa of a group share the same target clade,
##          so the grafting of "BD_clades" fits one birth-death model per clade (not one per taxon).
##  Part 1. Fixed costs. (a) Preparation (done once, in the main process, before the simulations) and time of one simulation: sequential runs with 1 and 2 simulations
##          (time = preparation + N x one simulation). (b) Parallel fixed costs with instant simulations: "cold" (plan created and workers started inside the timing,
##          package loaded by the workers) and "warm" (workers already running). Measured for each number of workers of `workers_measured`.
##  Part 2. Wall-clock times of N simulations: sequential (once), then parallel cold and parallel warm for each number of workers of `workers_measured`
##          (warm only for N = 10, since the fixed cost is then negligible and the runs are long for N >= 100).
##          Cells too long to run (more than `max_minutes` of sequential computation, e.g., 1000 tips x 1000 simulations) are not run: they come from the model of Part 3.
##  Part 3. Model: time_parallel(N, W) = preparation + fixed cost(W) + N x t1 x (1 + overhead) / W, with t1 = time of one simulation, W = workers, and the overhead
##          (transfer of the data and of the results, slower workers) estimated on the runs for each size of phylogeny and number of measured workers. It is checked on the runs,
##          then used for the number of workers of `workers_model` (workers that were not measured, e.g. 10, use the overhead of the largest number of measured workers and a linear
##          extrapolation of the fixed cost: ideal workers, no contention for memory or cores, so these are upper bounds of the speed-up) and for the break-even
##          (number of simulations from which parallel is faster).
##
## How to run: install the package first (the workers need it), then source this file. Results (csv) are written in `out_dir`.
##   `quick = TRUE`: only the 100 tips case, N = 10 and 100 (a few minutes).
##   Use as many physical cores as the largest value of `workers_measured`, and close heavy applications during the run.

library(phylograftR)
library(ape)

workers_measured <- c(2, 4)   # Workers of the plans used for the measures
workers_model <- c(2, 4, 10)  # Workers of the model
max_minutes <- 15             # Cells with more than this expected sequential time (in minutes) are not run (model)
quick <- FALSE
out_dir <- "."
reuse_csv <- FALSE            # TRUE: skip the measures and redo Part 3 (model and recommendations) from the csv files already written in `out_dir`

# --- Backbone and taxa to graft ---

# Random coalescent backbone (root age 50). The `n_missing` taxa are spread over `n_groups` disjoint target clades (nodes with `min_clade` tips or more, at most half of the tips)
# with a broken stick distribution (each group has at least one taxon; the most taxa go to the largest clades). Taxa of the same group share the same target clade (all its tips).
make_case <- function (n_tips, n_missing, n_groups = ceiling(log(n_tips)), min_clade = 5, seed = 1)
{
  set.seed(seed)
  phy <- rcoal(n_tips); phy$edge.length <- phy$edge.length / max(node.depth.edgelength(phy)) * 50
  clades <- ape::prop.part(phy)                      # Tip indices descending from each internal node
  sizes <- lengths(clades)
  candidates <- which(sizes >= min_clade & sizes <= n_tips / 2)
  best <- integer(0)
  for (attempt in 1:200)                             # Random disjoint clades (greedy), retried until `n_groups` are found
  {
    chosen <- integer(0); used <- integer(0)
    for (i in sample(candidates, prob = sizes[candidates]))   # Larger clades are more likely to be chosen
      if (!any(clades[[i]] %in% used)) { chosen <- c(chosen, i); used <- c(used, clades[[i]]); if (length(chosen) == n_groups) break }
    if (length(chosen) > length(best)) best <- chosen
    if (length(best) == n_groups) break
  }
  if (length(best) < n_groups) warning("Only ", length(best), " disjoint target clades found (", n_groups, " requested).")
  n_groups <- length(best)
  share <- sort(diff(c(0, sort(stats::runif(n_groups - 1)), 1)), decreasing = TRUE)   # Broken stick, the largest shares going to the largest clades
  best <- best[order(sizes[best], decreasing = TRUE)]
  n_per_group <- as.vector(stats::rmultinom(1, n_missing - n_groups, share)) + 1
  clade_tips <- lapply(best, function (i) phy$tip.label[clades[[i]]])
  nm <- paste0("X", seq_len(n_missing))
  tl <- stats::setNames(clade_tips[rep(seq_len(n_groups), n_per_group)], nm)
  structure(list(phy = phy, mtl = list(Taxon_name = nm, Target_clade = tl)), n_per_group = sort(n_per_group, decreasing = TRUE), clade_sizes = sizes[best])
}
graft <- function (cs, method, nb_simulations = 1, parallel = FALSE)
{
  graft_missing_taxa(cs$phy, cs$mtl, method = method, nb_simulations = nb_simulations, parallel = parallel, verbose = FALSE, seed = 1)
}
elapsed <- function (expr) unname(system.time(suppressWarnings(expr))[3])

# Package loaded and every namespace used by the grafting loaded in the workers (a short run)
warm_up <- function (cs, W) invisible(graft(cs, "uniform", nb_simulations = W * 2, parallel = TRUE))
with_plan <- function (expr, n_workers) { old <- future::plan(future::multisession, workers = n_workers); on.exit(future::plan(old)); expr }

sizes <- if (quick) 100 else c(100, 500, 1000)
methods <- c("uniform", "BD_clades")
n_sims_grid <- if (quick) c(10, 100) else c(10, 100, 1000)

if (reuse_csv)
{
  unit <- utils::read.csv(file.path(out_dir, "benchmark_parallel_part1a.csv"))
  fx <- utils::read.csv(file.path(out_dir, "benchmark_parallel_part1b.csv"))
  cells <- utils::read.csv(file.path(out_dir, "benchmark_parallel_part2.csv"))
} else {

  # --- Part 1a: preparation and time of one simulation (sequential) ---
  cat("\n== Part 1a: preparation and time of one simulation (sequential) ==\n")
  unit <- list()
  for (n in sizes) for (m in methods)
  {
    cs <- make_case(n, n)
    t1 <- elapsed(graft(cs, m, 1)); t2 <- elapsed(graft(cs, m, 2))
    one <- t2 - t1
    unit[[length(unit) + 1]] <- data.frame(n_tips = n, n_missing = n, n_groups = length(attr(cs, "n_per_group")), method = m, t_one_sim = one, t_prep = max(t1 - one, 0))
    cat(sprintf("tips %4d | %-9s | %d target clades (taxa per clade: %s) | one simulation %7.2f s | preparation %5.2f s\n", n, m, length(attr(cs, "n_per_group")),
                paste(attr(cs, "n_per_group"), collapse = "/"), one, max(t1 - one, 0)))
  }
  unit <- do.call(rbind, unit)

  # --- Part 1b: fixed costs of the parallel run, with instant simulations ---
  cat("\n== Part 1b: fixed costs of the parallel run, with instant simulations ==\n")
  cs_small <- make_case(30, 4, n_groups = 1)
  fast <- function (cs, n_sims) graft_missing_taxa(cs$phy, cs$mtl, method = "uniform", nb_simulations = n_sims, parallel = TRUE, verbose = FALSE, seed = 1)
  fx <- list()
  for (W in workers_measured)
  {
    fixed_cold <- vapply(1:3, function (i) elapsed(with_plan(fast(cs_small, W * 2), W)), numeric(1))   # a new plan each time
    old <- future::plan(future::multisession, workers = W); warm_up(cs_small, W)
    fixed_warm <- vapply(1:3, function (i) elapsed(fast(cs_small, W * 2)), numeric(1)); future::plan(old)
    seq_small <- elapsed(graft_missing_taxa(cs_small$phy, cs_small$mtl, method = "uniform", nb_simulations = W * 2, verbose = FALSE, seed = 1))
    fx[[length(fx) + 1]] <- data.frame(workers = W, fixed_cold = median(fixed_cold) - seq_small, fixed_warm = median(fixed_warm) - seq_small)
    cat(sprintf("%d workers | fixed cost, cold (plan creation, workers start, package loaded in the workers): %.2f s | warm (workers already running): %.2f s\n",
                W, fx[[length(fx)]]$fixed_cold, fx[[length(fx)]]$fixed_warm))
  }
  fx <- do.call(rbind, fx)

  # --- Part 2: wall-clock times ---
  cat("\n== Part 2: sequential vs parallel (", paste(workers_measured, collapse = " and "), " workers) ==\n", sep = "")
  cells <- list()
  for (n in sizes) for (m in methods)
  {
    u <- unit[unit$n_tips == n & unit$method == m, ]
    for (N in n_sims_grid)
    {
      # Budget: skip cells with more than `max_minutes` of sequential computation
      t_seq_expected <- u$t_prep + N * u$t_one_sim
      if (t_seq_expected > max_minutes * 60) { cat(sprintf("tips %4d | %-9s | N %4d | not run (expected sequential time %.0f s): model\n", n, m, N, t_seq_expected)); next }
      cs <- make_case(n, n)
      t_seq <- elapsed(graft(cs, m, N, FALSE))
      for (W in workers_measured)
      {
        t_cold <- elapsed(with_plan(graft(cs, m, N, TRUE), W))
        t_warm <- NA_real_
        if (N <= 10) { old <- future::plan(future::multisession, workers = W); warm_up(cs, W); t_warm <- elapsed(graft(cs, m, N, TRUE)); future::plan(old) }
        cells[[length(cells) + 1]] <- data.frame(n_tips = n, n_missing = n, method = m, n_sims = N, workers = W, t_seq = t_seq, t_par_cold = t_cold, t_par_warm = t_warm)
        cat(sprintf("tips %4d | %-9s | N %4d | %d workers | seq %8.1f s | cold %8.1f s (x%.2f) | warm %s\n", n, m, N, W, t_seq, t_cold, t_seq / t_cold,
                    if (is.na(t_warm)) "-" else sprintf("%8.1f s (x%.2f)", t_warm, t_seq / t_warm)))
        utils::write.csv(do.call(rbind, cells), file.path(out_dir, "benchmark_parallel_part2.csv"), row.names = FALSE)   # Saved after each run
      }
    }
  }
  cells <- do.call(rbind, cells)
  utils::write.csv(unit, file.path(out_dir, "benchmark_parallel_part1a.csv"), row.names = FALSE)
  utils::write.csv(fx, file.path(out_dir, "benchmark_parallel_part1b.csv"), row.names = FALSE)
}

# --- Part 3: model ---

cat("\n== Part 3: model ==\n")
cells <- merge(cells, unit, by = c("n_tips", "n_missing", "method"))
# Fixed cost of the parallel run for W workers: measured for `workers_measured`, linear extrapolation otherwise (workers start one after the other)
fixed_cost <- function (W, cold = TRUE)
{
  y <- if (cold) fx$fixed_cold else fx$fixed_warm
  vapply(W, function (w) if (w %in% fx$workers) y[fx$workers == w] else
    max(0, unname(stats::predict(stats::lm(y ~ workers, data.frame(y = y, workers = fx$workers)), data.frame(workers = w)))), numeric(1))
}
# Overhead of the parallel computation on the simulations: time of the simulations in parallel / ideal time (N x t1 / W), minus 1, estimated on the runs with N >= 100 (cold) and N = 10 (warm)
fc <- ifelse(is.na(cells$t_par_warm), fixed_cost(cells$workers, TRUE), fixed_cost(cells$workers, FALSE))
cells$overhead <- (ifelse(is.na(cells$t_par_warm), cells$t_par_cold, cells$t_par_warm) - fc - cells$t_prep) / (cells$n_sims * cells$t_one_sim / cells$workers) - 1
print(cells[, c("n_tips", "method", "n_sims", "workers", "t_seq", "t_par_cold", "t_par_warm", "overhead")], digits = 3, row.names = FALSE)
# The overhead depends on the size of the phylogeny (a larger share of the time is spent in the simulation itself) and on the number of workers: median per number of tips and workers
ov <- stats::aggregate(overhead ~ n_tips + workers, data = cells, FUN = stats::median)
cat("\nMedian overhead on the simulations, by number of tips and workers:\n"); print(ov, digits = 2, row.names = FALSE)
overhead_of <- function (n_tips, W)   # Not measured workers (e.g., 10) use the overhead of the largest number of measured workers
{
  k <- max(length(n_tips), length(W)); n_tips <- rep_len(n_tips, k); W <- rep_len(W, k)
  Wm <- ifelse(W %in% ov$workers, W, max(ov$workers))
  vapply(seq_len(k), function (i) ov$overhead[ov$n_tips == n_tips[i] & ov$workers == Wm[i]], numeric(1))
}
predict_par <- function (u, N, W, cold = TRUE) u$t_prep + fixed_cost(W, cold) + N * u$t_one_sim * (1 + overhead_of(u$n_tips, W)) / W
# Check of the model on the runs (cold)
chk <- cells[!is.na(cells$t_par_cold), ]
chk$predicted_cold <- predict_par(chk, chk$n_sims, chk$workers)
cat("\nModel check (cold), measured vs predicted time of the parallel runs:\n")
print(chk[, c("n_tips", "method", "n_sims", "workers", "t_par_cold", "predicted_cold")], digits = 3, row.names = FALSE)

# Predictions: speed-up (sequential / parallel, cold plan) for each size, method, N and number of workers
pred <- expand.grid(n_tips = sizes, method = methods, n_sims = n_sims_grid, workers = workers_model, stringsAsFactors = FALSE)
pred <- merge(pred, unit, by = c("n_tips", "method"))
pred$t_seq <- pred$t_prep + pred$n_sims * pred$t_one_sim
pred$t_par <- predict_par(pred, pred$n_sims, pred$workers)
pred$speedup <- pred$t_seq / pred$t_par
pred$measured <- paste(pred$n_tips, pred$method, pred$n_sims, pred$workers) %in% paste(cells$n_tips, cells$method, cells$n_sims, cells$workers)
pred <- pred[order(pred$method, pred$n_tips, pred$n_sims, pred$workers), c("n_tips", "method", "n_sims", "workers", "t_seq", "t_par", "speedup", "measured")]
print(pred, digits = 3, row.names = FALSE)
utils::write.csv(pred, file.path(out_dir, "benchmark_parallel_part3_predictions.csv"), row.names = FALSE)

# --- Recommendations: from which number of grafted phylogenies (P maps x N simulations) is parallel computation worth it? ---

# Worth it = (1) faster than sequential, cold plan (break-even), and (2) clearly worth it: speed-up of at least 1.5 and at least 1 minute saved.
cat("\nRecommendations: minimum number of grafted phylogenies (maps x simulations) from which parallel computation is worth it (cold plan):\n")
grid_n <- unique(round(exp(seq(log(2), log(5000), length.out = 400))))
rec <- list()
for (i in seq_len(nrow(unit))) for (W in workers_model)
{
  u <- unit[i, ]
  t_seq <- u$t_prep + grid_n * u$t_one_sim
  t_par <- predict_par(u, grid_n, W)
  first <- function (ok) if (any(ok)) grid_n[which(ok)[1]] else NA
  rec[[length(rec) + 1]] <- data.frame(n_tips = u$n_tips, method = u$method, workers = W, time_one_simulation_s = u$t_one_sim,
                                       min_n_break_even = first(t_seq > t_par),
                                       min_n_worth_it = first(t_seq / t_par >= 1.5 & t_seq - t_par >= 60),
                                       speedup_at_100 = (u$t_prep + 100 * u$t_one_sim) / predict_par(u, 100, W),
                                       speedup_at_1000 = (u$t_prep + 1000 * u$t_one_sim) / predict_par(u, 1000, W))
}
rec <- do.call(rbind, rec)
print(rec, digits = 3, row.names = FALSE)
utils::write.csv(rec, file.path(out_dir, "benchmark_parallel_recommendations.csv"), row.names = FALSE)
sessionInfo()
