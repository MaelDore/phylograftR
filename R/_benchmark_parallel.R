
### Benchmark script: when is parallel computation of the simulations worth it? (graft_missing_taxa(), `parallel` argument)

## Design (same number of taxa to graft as tips in the backbone: 100, 500 and 1000 tips; 10, 100 and 1000 simulations; both methods, no simmap)
##  Part 1. Fixed costs. (a) Preparation (done once, in the main process, before the simulations) and time of one simulation: sequential runs with 1 and 2 simulations
##          (time = preparation + N x one simulation). (b) Parallel fixed costs with instant simulations: "cold" (plan created and workers started inside the timing,
##          package loaded by the workers) and "warm" (workers already running).
##  Part 2. Wall-clock times of N simulations: sequential, parallel cold and parallel warm (cold only for N >= 100, since the fixed cost is then negligible and the runs are long).
##          Cells too long to run (e.g., 1000 tips x 1000 simulations) are not run: they come from the model of Part 3.
##  Part 3. Model: time_parallel(N, W) = preparation + fixed cost + N x t1 x (1 + overhead) / W, with t1 = time of one simulation, W = workers, and the overhead
##          (transfer of the data and of the results, slower workers) estimated on the runs for each size of phylogeny. It is checked on the runs, then used for 4 and 10 workers (ideal workers: no contention
##          for memory or cores, so these are upper bounds of the speed-up) and for the break-even (number of simulations from which parallel is faster).
##
## How to run: install the package first (the workers need it), then source this file. Results (csv) are written in `out_dir`.
##   `quick = TRUE`: only the 100 tips case, few simulations (a few minutes).

library(phylograftR)
library(ape)

workers <- 2           # Workers of the plan used for the measures
workers_model <- c(2, 4, 10)  # Workers of the model
quick <- FALSE
out_dir <- "."
reuse_csv <- FALSE     # TRUE: skip the measures and redo Part 3 (model and recommendations) from the csv files already written in `out_dir`

# --- Markov simulator of a discrete map on a tree (not used by default: no simmap) ---

make_case <- function (n_tips, n_missing, seed = 1)
{
  set.seed(seed)
  phy <- rcoal(n_tips); phy$edge.length <- phy$edge.length / max(node.depth.edgelength(phy)) * 50
  # Each missing taxon: target clade = 2 tips of a random internal node
  nodes <- sample((n_tips + 1):(n_tips + phy$Nnode), n_missing, replace = TRUE)
  tl <- lapply(nodes, function (nd) { d <- phytools::getDescendants(phy, nd); d <- d[d <= n_tips]; unique(phy$tip.label[sample(d, min(2, length(d)))]) })
  nm <- paste0("X", seq_len(n_missing))
  list(phy = phy, mtl = list(Taxon_name = nm, Target_clade = stats::setNames(tl, nm)))
}
graft <- function (cs, method, nb_simulations = 1, parallel = FALSE)
{
  graft_missing_taxa(cs$phy, cs$mtl, method = method, nb_simulations = nb_simulations, parallel = parallel, verbose = FALSE, seed = 1)
}
elapsed <- function (expr) unname(system.time(suppressWarnings(expr))[3])

# Package loaded and every namespace used by the grafting loaded in the workers (a short run)
warm_up <- function (cs) invisible(graft(cs, "uniform", nb_simulations = workers * 2, parallel = TRUE))
with_plan <- function (expr, n_workers = workers) { old <- future::plan(future::multisession, workers = n_workers); on.exit(future::plan(old)); expr }

sizes <- if (quick) 100 else c(100, 500, 1000)
methods <- c("uniform", "BD_clades")

if (reuse_csv)
{
  unit <- utils::read.csv(file.path(out_dir, "benchmark_parallel_part1a.csv"))
  fx <- utils::read.csv(file.path(out_dir, "benchmark_parallel_part1b.csv")); fixed_cold <- fx$fixed_cold; fixed_warm <- fx$fixed_warm; workers <- fx$workers
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
    unit[[length(unit) + 1]] <- data.frame(n_tips = n, n_missing = n, method = m, t_one_sim = one, t_prep = max(t1 - one, 0))
    cat(sprintf("tips %4d | %-9s | one simulation %7.2f s | preparation %5.2f s\n", n, m, one, max(t1 - one, 0)))
  }
  unit <- do.call(rbind, unit)

  # --- Part 1b: fixed costs of the parallel run, with instant simulations ---
  cat("\n== Part 1b: fixed costs of the parallel run (", workers, " workers) ==\n", sep = "")
  cs_small <- make_case(30, 4)
  fast <- function (cs, n_sims) graft_missing_taxa(cs$phy, cs$mtl, method = "uniform", nb_simulations = n_sims, parallel = TRUE, verbose = FALSE, seed = 1)
  fixed_cold <- vapply(1:3, function (i) elapsed(with_plan(fast(cs_small, workers * 2))), numeric(1))   # a new plan each time
  old <- future::plan(future::multisession, workers = workers); warm_up(cs_small)
  fixed_warm <- vapply(1:3, function (i) elapsed(fast(cs_small, workers * 2)), numeric(1)); future::plan(old)
  seq_small <- elapsed(graft_missing_taxa(cs_small$phy, cs_small$mtl, method = "uniform", nb_simulations = workers * 2, verbose = FALSE, seed = 1))
  fixed_cold <- median(fixed_cold) - seq_small; fixed_warm <- median(fixed_warm) - seq_small
  cat(sprintf("Fixed cost, cold (plan creation, workers start, package loaded in the workers): %.2f s\nFixed cost, warm (workers already running, package loaded):                      %.2f s\n", fixed_cold, fixed_warm))

  # --- Part 2: wall-clock times ---
  cat("\n== Part 2: sequential vs parallel, ", workers, " workers ==\n", sep = "")
  cells <- list()
  for (n in sizes) for (m in methods)
  {
    u <- unit[unit$n_tips == n & unit$method == m, ]
    for (N in c(10, 100, 1000))
    {
      if (quick && N > 100) next
      # Budget: skip cells with more than ~15 minutes of sequential computation
      t_seq_expected <- u$t_prep + N * u$t_one_sim
      if (t_seq_expected > 15 * 60) { cat(sprintf("tips %4d | %-9s | N %4d | not run (expected sequential time %.0f s): model\n", n, m, N, t_seq_expected)); next }
      cs <- make_case(n, n)
      t_seq <- elapsed(graft(cs, m, N, FALSE))
      t_cold <- elapsed(with_plan(graft(cs, m, N, TRUE)))
      t_warm <- NA_real_
      if (N <= 10) { old <- future::plan(future::multisession, workers = workers); warm_up(cs); t_warm <- elapsed(graft(cs, m, N, TRUE)); future::plan(old) }
      cells[[length(cells) + 1]] <- data.frame(n_tips = n, n_missing = n, method = m, n_sims = N, workers = workers, t_seq = t_seq, t_par_cold = t_cold, t_par_warm = t_warm)
      cat(sprintf("tips %4d | %-9s | N %4d | seq %8.1f s | cold %8.1f s (x%.2f) | warm %s\n", n, m, N, t_seq, t_cold, t_seq / t_cold,
                  if (is.na(t_warm)) "-" else sprintf("%8.1f s (x%.2f)", t_warm, t_seq / t_warm)))
      utils::write.csv(do.call(rbind, cells), file.path(out_dir, "benchmark_parallel_part2.csv"), row.names = FALSE)
    }
  }
  cells <- do.call(rbind, cells)
  utils::write.csv(unit, file.path(out_dir, "benchmark_parallel_part1a.csv"), row.names = FALSE)
  utils::write.csv(data.frame(fixed_cold = fixed_cold, fixed_warm = fixed_warm, workers = workers), file.path(out_dir, "benchmark_parallel_part1b.csv"), row.names = FALSE)
}

# --- Part 3: model ---

cat("\n== Part 3: model ==\n")
# Overhead of the parallel computation on the simulations: time of the simulations in parallel / ideal time (N x t1 / W), minus 1, estimated on the runs with N >= 100 (cold) and N = 10 (warm)
cells <- merge(cells, unit, by = c("n_tips", "n_missing", "method"))
ideal <- cells$n_sims * cells$t_one_sim / cells$workers
cells$overhead <- (ifelse(is.na(cells$t_par_warm), cells$t_par_cold - fixed_cold, cells$t_par_warm - fixed_warm) - cells$t_prep) / ideal - 1
print(cells[, c("n_tips", "method", "n_sims", "t_seq", "t_par_cold", "t_par_warm", "overhead")], digits = 3, row.names = FALSE)
# The overhead depends on the size of the phylogeny (a larger share of the time is spent in the simulation itself): median per number of tips
overhead_by_size <- tapply(cells$overhead, cells$n_tips, stats::median)
cat("\nMedian overhead on the simulations, by number of tips:\n"); print(round(overhead_by_size, 2))
predict_par <- function (u, N, W, cold = TRUE) u$t_prep + (if (cold) fixed_cold else fixed_warm) + N * u$t_one_sim * (1 + unname(overhead_by_size[as.character(u$n_tips)])) / W
# Check of the model on the runs (cold)
chk <- cells[!is.na(cells$t_par_cold), ]
chk$predicted_cold <- predict_par(chk, chk$n_sims, chk$workers)
cat("\nModel check (cold), measured vs predicted time of the parallel runs:\n")
print(chk[, c("n_tips", "method", "n_sims", "t_par_cold", "predicted_cold")], digits = 3, row.names = FALSE)

# Predictions: speed-up (sequential / parallel, cold plan) for each size, method, N and number of workers
pred <- expand.grid(n_tips = sizes, method = methods, n_sims = c(10, 100, 1000), workers = workers_model, stringsAsFactors = FALSE)
pred <- merge(pred, unit, by = c("n_tips", "method"))
pred$t_seq <- pred$t_prep + pred$n_sims * pred$t_one_sim
pred$t_par <- predict_par(pred, pred$n_sims, pred$workers)
pred$speedup <- pred$t_seq / pred$t_par
pred$measured <- paste(pred$n_tips, pred$method, pred$n_sims) %in% paste(cells$n_tips, cells$method, cells$n_sims) & pred$workers == workers
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
