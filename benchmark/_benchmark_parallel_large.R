### Benchmark script (large backbones and memory): complements `_benchmark_parallel.R` (graft_missing_taxa(), `parallel` argument)

## Why: in `_benchmark_parallel.R`, as many taxa are grafted as tips (a worst case). The cost of a simulation depends much more on the number of taxa to graft than on the number
##      of tips, so 5000 and 10,000 tips are benchmarked here with a realistic share of missing taxa (`frac_missing`, 20% by default: 1000 and 2000 taxa to graft).
##
## Design
##  Part A. Peak memory of the standard design (taxa to graft = tips; 100, 500 and 1000 tips): one sequential simulation, then one round of simulations with 2 and with 4 workers (cold plan).
##          Peak memory = sum of the memory (resident set) of the main R process and of its workers, sampled every 0.5 s by a separate R process (needs the package `ps`).
##  Part B. Large backbones (`sizes_large`, `frac_missing` of the tips as missing taxa, methods `uniform` and `BD_clades`):
##          (1) one sequential simulation (time t1 and peak memory); (2) `n_sims_par` simulations in parallel (cold plan) with 4 and then 2 workers (time, speed-up vs `n_sims_par` x t1,
##          peak memory). Sequential time is linear in the number of simulations, so it is not measured for several simulations: it is `n_sims_par` x t1.
##          Cells are run in order of priority (see `queue`), and skipped when they would exceed `max_hours` in total. Results (csv) are saved after each step.
##  Part C. Model: time_parallel(N, W) = fixed cost(W) + ceiling(N / W) x t1 x (1 + overhead), with the overhead estimated on the runs of Part B (largest backbone first).
##          Predictions for N = 10, 100, 1000 simulations and 2, 4, 10 workers (10 workers are not measured: upper bound), written in a csv.
##
## How to run: install the package first (the workers need it), plus `ps` for the memory, then source this file. Use at least 4 free physical cores and close heavy applications.
##   Duration: roughly 3 to 5 hours with the default settings (1 simulation of 10,000 tips takes several minutes to more than an hour, depending on the machine). Set `max_hours` to
##   the time you have: cells that would not fit are skipped. `quick = TRUE`: tiny test of the script (a few minutes).

library(phylograftR)
library(ape)

sizes_large <- c(1000, 5000, 10000)   # Number of tips of the backbones of Part B
frac_missing <- 0.2                   # Share of missing taxa (number of taxa to graft / number of tips)
workers_measured <- c(4, 2)           # Workers of the plans of Part B (in this order)
n_sims_par <- 8                       # Simulations of the parallel runs (2 rounds with 4 workers, 4 rounds with 2 workers)
max_hours <- 7                        # Total budget: a step is skipped if its expected end is beyond this time
fixed_cold <- c(`2` = 13.5, `4` = 18.2)   # Fixed costs (s) of a cold plan, measured on the standard benchmark (cold minus warm), used for the budget and the model
quick <- FALSE
out_dir <- "."

if (quick) { sizes_large <- c(100, 200); max_hours <- 0.5 }
t_start <- Sys.time()
hours_elapsed <- function () as.numeric(difftime(Sys.time(), t_start, units = "hours"))

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


# --- Peak memory of a process and its workers ---

have_ps <- requireNamespace("ps", quietly = TRUE)
if (!have_ps) message("Package 'ps' not installed: the peak memory will not be measured (NA).")
start_mem_monitor <- function ()
{
  if (!have_ps) return(NULL)
  stop_file <- tempfile("stop_"); out_file <- tempfile("peak_"); script <- tempfile(fileext = ".R")
  # The sampler sums the memory (rss) of: the main process, its children (recursively), and the processes launched by it as workers (their command line contains "parallelly.parent=<pid>")
  writeLines(c(
    "a <- commandArgs(TRUE); parent <- as.integer(a[1]); stop_file <- a[2]; out_file <- a[3]; me <- Sys.getpid(); peak <- 0",
    "total_mem <- function () {",
    "  pids <- parent",
    "  kids <- tryCatch(ps::ps_children(ps::ps_handle(parent), recursive = TRUE), error = function (e) list())",
    "  pids <- c(pids, vapply(kids, ps::ps_pid, numeric(1)))",
    "  procs <- tryCatch(ps::ps(), error = function (e) NULL)",
    "  if (!is.null(procs)) for (i in which(grepl('^[Rr]', procs$name) & !(procs$pid %in% pids)))",
    "    if (isTRUE(tryCatch(grepl(paste0('parallelly.parent=', parent, '[^0-9]'), paste(ps::ps_cmdline(procs$ps_handle[[i]]), collapse = ' ')), error = function (e) FALSE))) pids <- c(pids, procs$pid[i])",
    "  sum(vapply(setdiff(unique(pids), me), function (p) tryCatch(unname(ps::ps_memory_info(ps::ps_handle(p))[['rss']]), error = function (e) 0), numeric(1)))",
    "}",
    "repeat {",
    "  tot <- tryCatch(total_mem(), error = function (e) 0)",
    "  if (tot > peak) { peak <- tot; writeLines(as.character(peak), out_file) }",
    "  if (file.exists(stop_file)) break",
    "  Sys.sleep(0.5)",
    "}"), script)
  system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), Sys.getpid(), shQuote(stop_file), shQuote(out_file)), wait = FALSE, stdout = FALSE, stderr = FALSE)
  Sys.sleep(1)
  list(stop_file = stop_file, out_file = out_file)
}
stop_mem_monitor <- function (mon)
{
  if (is.null(mon)) return(NA_real_)
  file.create(mon$stop_file); Sys.sleep(1.5)
  if (file.exists(mon$out_file)) as.numeric(readLines(mon$out_file)) / 2^20 else NA_real_   # MB
}
# Time and peak memory (MB) of an expression (a plan created inside is a cold plan)
measure <- function (expr_fun) { mon <- start_mem_monitor(); t <- elapsed(expr_fun()); c(time_s = t, peak_mb = stop_mem_monitor(mon)) }

# --- Part A: peak memory of the standard design (taxa to graft = tips) ---

cat("\n== Part A: peak memory, standard design (taxa to graft = tips) ==\n")
mem_std <- list()
for (n in (if (quick) 100 else c(100, 500, 1000)))
{
  cs <- make_case(n, n)
  for (W in c(1, 2, 4))
  {
    N <- W   # One round: one simulation per worker
    r <- if (W == 1) measure(function () graft(cs, "uniform", 1, FALSE)) else measure(function () with_plan(graft(cs, "uniform", N, TRUE), W))
    mem_std[[length(mem_std) + 1]] <- data.frame(n_tips = n, n_missing = n, method = "uniform", workers = W, n_sims = N, time_s = r[["time_s"]], peak_mb = r[["peak_mb"]])
    cat(sprintf("tips %5d | %d worker(s) | %d simulation(s) | %8.1f s | peak memory %7.0f MB\n", n, W, N, r[["time_s"]], r[["peak_mb"]]))
    utils::write.csv(do.call(rbind, mem_std), file.path(out_dir, "benchmark_parallel_large_memory_standard.csv"), row.names = FALSE)
  }
}

# --- Part B: large backbones, realistic share of missing taxa ---

cat("\n== Part B: large backbones (", frac_missing * 100, "% of missing taxa) ==\n", sep = "")
queue <- rbind(expand.grid(n_tips = sizes_large, method = "uniform", stringsAsFactors = FALSE),    # Priority: `uniform` for all sizes, then `BD_clades`
               expand.grid(n_tips = sizes_large, method = "BD_clades", stringsAsFactors = FALSE))
large <- list()
save_large <- function () utils::write.csv(do.call(rbind, large), file.path(out_dir, "benchmark_parallel_large.csv"), row.names = FALSE)
add_row <- function (n, nm, m, W, N, r, t1 = NA) {
  large[[length(large) + 1]] <<- data.frame(n_tips = n, n_missing = nm, method = m, workers = W, n_sims = N, time_s = r[["time_s"]], peak_mb = r[["peak_mb"]],
                                            speedup_vs_sequential = if (W > 1) N * t1 / r[["time_s"]] else NA)
  save_large()
}
t1_est <- NULL   # Expected time of one simulation, from the previous size (scaling with the square of the number of tips) for the budget
for (i in seq_len(nrow(queue)))
{
  n <- queue$n_tips[i]; m <- queue$method[i]; nm <- round(n * frac_missing)
  prev <- Filter(function (x) x$method == m && x$workers == 1, large)
  if (length(prev) > 0) { p <- prev[[length(prev)]]; t1_est <- p$time_s * (n / p$n_tips)^2 } else t1_est <- NULL
  if (!is.null(t1_est) && hours_elapsed() + t1_est / 3600 > max_hours) { cat(sprintf("tips %5d | %-9s | not run (expected %.1f h for one simulation, %.1f h already used)\n", n, m, t1_est / 3600, hours_elapsed())); next }
  cs <- make_case(n, nm)
  r1 <- measure(function () graft(cs, m, 1, FALSE)); t1 <- r1[["time_s"]]
  add_row(n, nm, m, 1, 1, r1)
  cat(sprintf("tips %5d | %5d taxa | %-9s | sequential, one simulation: %8.1f s | peak memory %7.0f MB | %.2f h used\n", n, nm, m, t1, r1[["peak_mb"]], hours_elapsed()))
  for (W in workers_measured)
  {
    expected <- fixed_cold[[as.character(W)]] + ceiling(n_sims_par / W) * t1 * 1.1
    if (hours_elapsed() + expected / 3600 > max_hours) { cat(sprintf("  %d workers | not run (expected %.1f h, %.1f h already used)\n", W, expected / 3600, hours_elapsed())); next }
    r <- measure(function () with_plan(graft(cs, m, n_sims_par, TRUE), W))
    add_row(n, nm, m, W, n_sims_par, r, t1)
    cat(sprintf("  %d workers | %d simulations | %8.1f s | x%.2f vs sequential | peak memory %7.0f MB | %.2f h used\n", W, n_sims_par, r[["time_s"]], n_sims_par * t1 / r[["time_s"]], r[["peak_mb"]], hours_elapsed()))
  }
}
large <- do.call(rbind, large)

# --- Part C: model and predictions ---

cat("\n== Part C: model ==\n")
seq_rows <- large[large$workers == 1, ]
par_rows <- large[large$workers > 1, ]
par_rows <- merge(par_rows, seq_rows[, c("n_tips", "method", "time_s")], by = c("n_tips", "method"), suffixes = c("", "_t1"))
par_rows$overhead <- (par_rows$time_s - fixed_cold[as.character(par_rows$workers)]) / (ceiling(par_rows$n_sims / par_rows$workers) * par_rows$time_s_t1) - 1
cat("Overhead of the simulations in parallel (time / ideal time - 1):\n"); print(par_rows[, c("n_tips", "method", "workers", "n_sims", "time_s", "overhead")], digits = 3, row.names = FALSE)
# Overhead per number of workers: median over the runs lasting at least 5 times the fixed cost (shorter runs are dominated by it), floored at 0
reliable <- par_rows[par_rows$time_s >= 5 * fixed_cold[as.character(par_rows$workers)], ]
ov <- if (nrow(reliable) > 0) stats::aggregate(overhead ~ workers, data = reliable, FUN = function (x) max(stats::median(x), 0)) else data.frame(workers = integer(0), overhead = numeric(0))
ov_of <- function (W) { if (nrow(ov) == 0) return(0); if (W %in% ov$workers) ov$overhead[ov$workers == W] else max(ov$overhead[ov$workers == max(ov$workers)], 0) }
fixed_of <- function (W) if (as.character(W) %in% names(fixed_cold)) fixed_cold[[as.character(W)]] else fixed_cold[["4"]] + (W - 4) * (fixed_cold[["4"]] - fixed_cold[["2"]]) / 2
pred <- do.call(rbind, lapply(seq_len(nrow(seq_rows)), function (i) { u <- seq_rows[i, ]
  expand.grid(n_tips = u$n_tips, n_missing = u$n_missing, method = u$method, n_sims = c(10, 100, 1000), workers = c(2, 4, 10), t1 = u$time_s, stringsAsFactors = FALSE) }))
pred$t_seq <- pred$n_sims * pred$t1
pred$t_par <- mapply(function (N, W, t1) fixed_of(W) + ceiling(N / W) * t1 * (1 + ov_of(W)), pred$n_sims, pred$workers, pred$t1)
pred$speedup <- pred$t_seq / pred$t_par
pred$measured <- pred$workers %in% workers_measured & pred$n_sims == n_sims_par
cat("\nPredictions (10 workers: extrapolation, upper bound):\n"); print(pred[order(pred$method, pred$n_tips, pred$n_sims, pred$workers), c("n_tips", "n_missing", "method", "n_sims", "workers", "t_seq", "t_par", "speedup")], digits = 3, row.names = FALSE)
utils::write.csv(pred, file.path(out_dir, "benchmark_parallel_large_predictions.csv"), row.names = FALSE)
sessionInfo()
