
##### Test suite for the parallel computation of the simulations #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Helper functions to manage parallelize architecture

# Skip if the packages are missing, or if the workers cannot load the package phylograftR (case when it is only loaded with devtools::load_all())
skip_if_no_workers <- function ()
{
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  if (is.null(.parallel_cache$ok))
  {
    old <- future::plan(future::multisession, workers = 2)
    on.exit(future::plan(old), add = TRUE)
    .parallel_cache$ok <- tryCatch(isTRUE(future::value(future::future(requireNamespace("phylograftR", quietly = TRUE)))), error = function (e) FALSE)
  }
  if (!.parallel_cache$ok) { testthat::skip("the workers cannot load the package phylograftR (install it first)") }
}
.parallel_cache <- new.env()

# Run an expression with a plan with 2 workers, and reset the previous plan afterwards
with_workers <- function (expr, workers = 2)
{
  old <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old), add = TRUE)
  force(expr)
}


## Build toy lists of missing taxa / missing clades, and multiSimmap (same as in test-multiple_simmaps.R)

maps_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}
maps_missing_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}

# The fake simulations are defined in test-utils_simulations.R / test-multiple_simmaps.R and are not available in this file
fake_sim <- function (failing = integer(0), warn = FALSE)
{
  function (sim_seed, sim_index, verbose_sim)
  {
    if (warn) { warning("same warning in each simulation") }
    if (!is.null(sim_seed) && sim_seed %in% failing) { stop_simulation_failure("failed for seed ", sim_seed) }
    set_local_seed(sim_seed) # As the grafting functions: the random number generator is restored on exit
    phy <- ape::rtree(5)
    phy$info <- c(seed = sim_seed, index = sim_index)
    phy
  }
}


## ---- Arguments and rules: no worker needed ---- ##

test_that("check_parallel() accepts NULL, TRUE and FALSE only",
{
  expect_no_error(check_parallel(NULL)); expect_no_error(check_parallel(TRUE)); expect_no_error(check_parallel(FALSE))
  for (bad in list("yes", 1, NA, c(TRUE, FALSE), list(TRUE))) { expect_error(check_parallel(bad), "'parallel' must be NULL") }
})

test_that("resolve_parallel() follows the rules: FALSE and NULL never stop, TRUE stops without packages or plan, and warns with a single task",
{
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  old <- future::plan(future::sequential); on.exit(future::plan(old), add = TRUE)
  # No plan (sequential): NULL and FALSE are sequential, silently. TRUE stops
  expect_false(resolve_parallel(FALSE, 10)); expect_false(resolve_parallel(NULL, 10))
  expect_error(resolve_parallel(TRUE, 10), "needs a 'future' plan.*future::plan\\(future::multisession")
  # A plan with workers: NULL and TRUE use it, FALSE does not. A single task is run sequentially (and TRUE warns)
  future::plan(future::multisession, workers = 2)
  expect_true(resolve_parallel(NULL, 10)); expect_true(resolve_parallel(TRUE, 10)); expect_false(resolve_parallel(FALSE, 10))
  expect_no_warning(expect_false(resolve_parallel(NULL, 1)))
  expect_warning(expect_false(resolve_parallel(TRUE, 1)), "no effect with a single simulation")
  future::plan(future::sequential)
  # Missing packages: only TRUE stops, with the installation instruction
  local_mocked_bindings(have_parallel_packages = function () FALSE)
  expect_false(resolve_parallel(NULL, 10)); expect_false(resolve_parallel(FALSE, 10))
  expect_error(resolve_parallel(TRUE, 10), "requires the packages 'future' and 'future.apply'.*install.packages")
})

test_that("the grafting functions validate 'parallel' and stop with an informative error for parallel = TRUE without a plan",
{
  skip_if_not_installed("future.apply")
  toy <- toy_setup()
  old <- future::plan(future::sequential); on.exit(future::plan(old), add = TRUE)
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$phy$tip.label[1:2]))
  expect_error(graft_missing_taxa(toy$phy, mtl, parallel = "yes", verbose = FALSE), "'parallel' must be NULL")
  expect_error(graft_missing_clades(toy$phy, list(Taxon_name = "X1", Sister_clade = list(X1 = toy$phy$tip.label[1:2])), parallel = NA, verbose = FALSE), "'parallel' must be NULL")
  expect_error(graft_missing_taxa(toy$phy, mtl, nb_simulations = 3, parallel = TRUE, verbose = FALSE), "needs a 'future' plan")
  expect_no_error(graft_missing_taxa(toy$phy, mtl, nb_simulations = 3, parallel = NULL, verbose = FALSE))
})

test_that("a sequential run does not use future, even with a plan set (parallel = FALSE) and prints no line about the workers",
{
  skip_if_no_workers()
  with_workers({
    expect_silent(run_simulations(fake_sim(), 3, FALSE, seed = 7, verbose = FALSE, parallel = FALSE))
    expect_silent(run_simulations(fake_sim(), 3, FALSE, seed = 7, verbose = TRUE, parallel = FALSE))
  })
})


## ---- Same output with and without workers ---- ##

test_that("run_simulations() in parallel gives the same output as sequentially: seeds, names, failures, forced attempts and warnings",
{
  skip_if_no_workers()
  seeds <- derive_simulation_seeds(7, 40)
  for (case in list(list(n = 6, force = FALSE, failing = integer(0)),
                    list(n = 6, force = FALSE, failing = seeds[c(2, 5)]),
                    list(n = 6, force = TRUE, failing = seeds[c(1, 2, 4, 8)]),
                    list(n = 3, force = TRUE, failing = seeds[-(1:3)][c(1, 3)])))
  {
    run <- function (parallel) suppressWarnings(run_simulations(fake_sim(case$failing), case$n, case$force, seed = 7, parallel = parallel))
    sequential <- run(FALSE)
    with_workers(parallel <- run(TRUE))
    expect_identical(parallel, sequential)
    expect_identical(attr(parallel, "seeds"), attr(sequential, "seeds"))
    expect_identical(attr(parallel, "failed"), attr(sequential, "failed"))
  }
  # The warnings are gathered and sent once, as sequentially
  with_workers(expect_warning(run_simulations(fake_sim(warn = TRUE), 4, FALSE, seed = 7, parallel = TRUE), "same warning in each simulation"))
  messages <- character(0)
  with_workers(withCallingHandlers(run_simulations(fake_sim(warn = TRUE), 4, FALSE, seed = 7, parallel = TRUE),
                                   warning = function (w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }))
  expect_equal(messages, "same warning in each simulation")
})

test_that("an error that is not specific to a simulation stops the run, in parallel too, with the same message",
{
  skip_if_no_workers()
  generic <- function (sim_seed, sim_index, verbose_sim) stop("generic problem")
  with_workers(expect_error(run_simulations(generic, 4, FALSE, seed = 7, parallel = TRUE), "generic problem"))
  # No simulation obtained: same error as sequentially
  failing_all <- fake_sim(failing = derive_simulation_seeds(7, 3))
  with_workers(expect_error(suppressWarnings(run_simulations(failing_all, 3, FALSE, seed = 7, parallel = TRUE)), "None of the 3 simulation\\(s\\) could be obtained"))
})

test_that("without seed, the simulations in parallel use independent random number streams",
{
  skip_if_no_workers()
  random_sim <- function (sim_seed, sim_index, verbose_sim) { phy <- ape::rtree(5); phy$info <- stats::runif(1); phy }
  with_workers(expect_no_warning(out <- run_simulations(random_sim, 4, FALSE, seed = NULL, parallel = TRUE)))
  expect_length(out, 4)
  expect_null(attr(out, "seeds"))
  expect_equal(length(unique(vapply(out, function (x) x$info, numeric(1)))), 4)
})

test_that("with verbose = TRUE, the main session prints the number of workers and a line per simulation done, in parallel only",
{
  skip_if_no_workers()
  with_workers({
    out <- capture.output(run_simulations(fake_sim(), 3, FALSE, seed = 7, verbose = TRUE, parallel = TRUE))
    expect_match(out[1], "Running 3 simulation\\(s\\) in parallel on 2 worker\\(s\\) \\(plan: multisession\\)")
    expect_equal(sum(grepl("Simulation [0-9]/3 done", out)), 3)
    expect_silent(run_simulations(fake_sim(), 3, FALSE, seed = 7, verbose = FALSE, parallel = TRUE))
  })
})

test_that("run_maps() runs the maps x simulations as one set of tasks, with the same output as sequentially (names, seeds, maps, dropped maps)",
{
  skip_if_no_workers()
  prepared <- lapply(1:3, function (k) list(tracker = k, missing_taxa_list = list()))
  seeds <- derive_simulation_seeds(4, 4)
  sim <- function (tracker, missing_taxa_list, map_name, sim_seed, sim_index, verbose_sim)
  {
    if (map_name == "b" && sim_seed %in% seeds[2:3]) { stop_simulation_failure("no usable age for the map ", map_name) }
    if (map_name == "c") { stop_simulation_failure("no usable age for the map ", map_name) }
    set_local_seed(sim_seed); phy <- ape::rtree(5); phy$info <- c(tracker = tracker, index = sim_index); phy
  }
  run <- function (parallel, nb = 4, force = FALSE) suppressWarnings(run_maps(sim, prepared, c("a", "b", "c"), TRUE, nb, force, seed = 4, parallel = parallel))
  for (settings in list(list(nb = 4, force = FALSE), list(nb = 2, force = TRUE), list(nb = 1, force = FALSE)))
  {
    sequential <- run(FALSE, settings$nb, settings$force)
    with_workers(parallel <- run(TRUE, settings$nb, settings$force))
    expect_identical(parallel, sequential)
    expect_equal(names(parallel), c("a", "b"))
    expect_identical(attr(parallel, "failed_maps"), attr(sequential, "failed_maps"))
  }
  # The warning of a dropped map and the prefix of the map in the messages are kept
  messages <- character(0)
  with_workers(withCallingHandlers(run_maps(sim, prepared, c("a", "b", "c"), TRUE, 2, FALSE, seed = 4, parallel = TRUE),
                                   warning = function (w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }))
  expect_true(any(grepl("^Map 'c': none of the 2 simulation", messages)))
  expect_true(any(grepl("^Map 'b': Only 1 of the 2 simulations", messages)))
  # An error that is not a simulation failure stops the run and names its map. In parallel, the maps run at the same time and the first error that arrives is the one reported,
  # so it can be the one of any map (sequentially, it is always the first map)
  expect_error(run_maps(function (...) stop("generic problem"), prepared, c("a", "b", "c"), TRUE, 2, FALSE, seed = 4, parallel = FALSE), "Map 'a': generic problem")
  with_workers(expect_error(run_maps(function (...) stop("generic problem"), prepared, c("a", "b", "c"), TRUE, 2, FALSE, seed = 4, parallel = TRUE), "Map '[abc]': generic problem"))
})


## ---- Same output for the grafting functions ---- ##

for (name in c("taxa_uniform", "taxa_BD", "clades"))
{
  test_that(paste0("graft_missing_", sub("_.*", "", name), "() (", name, "): parallel and sequential give identical outputs for several simulations, and several maps x simulations"),
  {
    skip_if_no_workers()
    toy <- toy_setup()
    maps <- lapply(1:3, function (k) make_toy_simmap(toy$phy, seed = k))
    multi <- structure(maps, class = c("multiSimmap", "multiPhylo"))
    graft <- function (parallel, ..., simmap = NULL, allowed_states = NULL) suppressWarnings(
      if (name == "clades") graft_missing_clades(toy$phy, maps_missing_clades(toy), simmap = simmap, allowed_states = allowed_states, verbose = FALSE, seed = 11, parallel = parallel, ...)
      else graft_missing_taxa(toy$phy, maps_missing_taxa(toy), simmap = simmap, allowed_states = allowed_states, method = if (name == "taxa_BD") "BD_clades" else "uniform",
                              verbose = FALSE, seed = 11, parallel = parallel, ...))
    # Several simulations, no map
    sequential <- graft(FALSE, nb_simulations = 4)
    with_workers(parallel <- graft(TRUE, nb_simulations = 4))
    expect_identical(parallel, sequential)
    # The default (NULL) uses the plan
    with_workers(default <- graft(NULL, nb_simulations = 4))
    expect_identical(default, sequential)
    # Maps x simulations, with the states of the maps
    sequential <- graft(FALSE, nb_simulations = 3, simmap = multi, allowed_states = c("A", "B", "C"))
    with_workers(parallel <- graft(TRUE, nb_simulations = 3, simmap = multi, allowed_states = c("A", "B", "C")))
    expect_identical(parallel, sequential)
    expect_length(parallel, 3)
    expect_equal(names(parallel[[1]]), paste0("map_1_sim_", 1:3))
  })
}
