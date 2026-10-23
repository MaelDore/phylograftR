
##### Test suite for the helpers of the multiple simulations (utils_simulations.R) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Test the check of the arguments
test_that("check_nb_simulations() accepts positive integers and stops otherwise",
{
  expect_silent(check_nb_simulations(1, FALSE))
  expect_silent(check_nb_simulations(5L, TRUE))
  expect_silent(check_nb_simulations(3, FALSE)) # double with an integer value
  for (value in list(0, -1, 1.5, NA, NA_real_, "2", c(1, 2), NULL, TRUE))
  {
    expect_error(check_nb_simulations(value, FALSE), "'nb_simulations' must be a single positive integer")
  }
  for (value in list(NA, "TRUE", 1, c(TRUE, FALSE), NULL))
  {
    expect_error(check_nb_simulations(2, value), "'force_simulations' must be TRUE or FALSE")
  }
})

## Test the failure specific to a simulation
test_that("stop_simulation_failure() raises an error of class 'phylograftR_simulation_error' with the pasted message",
{
  err <- tryCatch(stop_simulation_failure("No usable age for ", 3, " taxa."), error = function (e) e)
  expect_s3_class(err, "phylograftR_simulation_error")
  expect_s3_class(err, "error")
  expect_equal(conditionMessage(err), "No usable age for 3 taxa.")
  # It is also caught as a usual error
  expect_error(stop_simulation_failure("failure"), "failure")
})

## Test the seeds of the simulations
test_that("derive_simulation_seeds() starts with the seed, gives distinct seeds, is a prefix-stable stream, and leaves the RNG of the session untouched",
{
  expect_null(derive_simulation_seeds(NULL, 5))
  expect_equal(derive_simulation_seeds(42, 1), 42)

  seeds <- derive_simulation_seeds(42, 10)
  expect_length(seeds, 10)
  expect_equal(seeds[1], 42)
  expect_equal(length(unique(seeds)), 10)
  expect_true(all(seeds[-1] >= 1 & seeds[-1] <= .Machine$integer.max))
  # Same seeds when called again, and the first seeds of a longer run are the seeds of a shorter run
  expect_equal(derive_simulation_seeds(42, 10), seeds)
  expect_equal(derive_simulation_seeds(42, 4), seeds[1:4])
  expect_equal(derive_simulation_seeds(42, 20)[1:10], seeds)
  # Another seed gives other seeds
  expect_false(identical(derive_simulation_seeds(43, 10)[-1], seeds[-1]))

  # The RNG of the session is not modified
  set.seed(99)
  before <- .Random.seed
  derive_simulation_seeds(42, 10)
  expect_identical(.Random.seed, before)
})

## Test the temporary file of each simulation
test_that("temp_file_path_simulation() adds the suffix of the simulation only with several simulations",
{
  expect_null(temp_file_path_simulation(NULL, 1, multiple = FALSE))
  expect_equal(temp_file_path_simulation("a/b.rds", 1, multiple = FALSE), "a/b.rds")
  expect_equal(temp_file_path_simulation("a/b.rds", 3, multiple = TRUE), "a/b_sim_3.rds")
  expect_equal(temp_file_path_simulation(NULL, 2, multiple = TRUE), "./temp_phy_from_grafting_sim_2.rds")
  expect_equal(temp_file_path_simulation("a.rds.rds", 2, multiple = TRUE), "a.rds_sim_2.rds")
})

## Helper: a simulation function returning a small phylogeny, failing for some seeds
fake_simulation <- function (failing_seeds = integer(0), warn = FALSE, generic_fail_seeds = integer(0))
{
  function (sim_seed, sim_index, verbose_sim)
  {
    if (warn) { warning("same warning in each simulation") }
    if (!is.null(sim_seed) && sim_seed %in% failing_seeds) { stop_simulation_failure("failed for seed ", sim_seed) }
    if (!is.null(sim_seed) && sim_seed %in% generic_fail_seeds) { stop("generic problem") }
    phy <- ape::rtree(4)
    phy$info <- c(seed = if (is.null(sim_seed)) NA else sim_seed, index = sim_index, verbose = verbose_sim)
    phy
  }
}

## Test the run of the simulations
test_that("run_simulations() returns the phylogeny for one simulation and a named multiPhylo with seeds for several",
{
  one <- run_simulations(fake_simulation(), 1, FALSE, seed = 7)
  expect_s3_class(one, "phylo")
  expect_false(inherits(one, "multiPhylo"))
  expect_equal(unname(one$info["seed"]), 7)

  several <- run_simulations(fake_simulation(), 3, FALSE, seed = 7)
  expect_s3_class(several, "multiPhylo")
  expect_equal(names(several), c("sim_1", "sim_2", "sim_3"))
  seeds <- derive_simulation_seeds(7, 3)
  expect_equal(attr(several, "seeds"), setNames(seeds, names(several)))
  expect_equal(unname(vapply(several, function (p) p$info[["seed"]], numeric(1))), seeds)
  expect_equal(nrow(attr(several, "failed")), 0)
  # The simulation 1 is the same as a single simulation
  expect_equal(several$sim_1$info[["seed"]], one$info[["seed"]])

  # A single multiPhylo element of length 1 stays a phylo (nb_simulations = 1)
  expect_equal(length(run_simulations(fake_simulation(), 2, FALSE, seed = 7)), 2)

  # Without seed: no seeds attribute, simulations follow the session
  without_seed <- run_simulations(fake_simulation(), 2, FALSE, seed = NULL)
  expect_null(attr(without_seed, "seeds"))
  expect_equal(names(without_seed), c("sim_1", "sim_2"))
})

test_that("run_simulations() gathers the warnings of the simulations and sends each once",
{
  expect_warning(run_simulations(fake_simulation(warn = TRUE), 3, FALSE, seed = 7), "same warning in each simulation")
  n_warnings <- 0
  withCallingHandlers(run_simulations(fake_simulation(warn = TRUE), 3, FALSE, seed = 7),
                      warning = function (w) { n_warnings <<- n_warnings + 1; invokeRestart("muffleWarning") })
  expect_equal(n_warnings, 1)
  # With a single simulation, the warning is sent as it is
  n_warnings <- 0
  withCallingHandlers(run_simulations(fake_simulation(warn = TRUE), 1, FALSE, seed = 7),
                      warning = function (w) { n_warnings <<- n_warnings + 1; invokeRestart("muffleWarning") })
  expect_equal(n_warnings, 1)
})

test_that("run_simulations() skips the simulations with a specific failure, keeps their index in the names, and warns with the number obtained",
{
  seeds <- derive_simulation_seeds(7, 4)
  expect_warning(out <- run_simulations(fake_simulation(failing_seeds = seeds[c(2, 3)]), 4, FALSE, seed = 7),
                 "Only 2 of the 4 simulations requested were obtained: 2 simulation\\(s\\) failed.*force_simulations = TRUE.*failed for seed")
  expect_equal(names(out), c("sim_1", "sim_4"))
  expect_equal(attr(out, "seeds"), setNames(seeds[c(1, 4)], c("sim_1", "sim_4")))
  expect_equal(attr(out, "failed")$seed, seeds[c(2, 3)])
  expect_match(attr(out, "failed")$reason[1], "failed for seed")

  # Warning of the failure message with the seeds of the failed simulations
  expect_warning(run_simulations(fake_simulation(failing_seeds = seeds[2]), 4, FALSE, seed = 7), paste0("Seeds of the failed simulations: ", seeds[2]))
})

test_that("run_simulations() with 'force_simulations' tries new seeds until the number requested is obtained, with consecutive names",
{
  seeds <- derive_simulation_seeds(7, 30)
  failing <- seeds[c(2, 3)]
  expect_no_warning(out <- run_simulations(fake_simulation(failing_seeds = failing), 4, TRUE, seed = 7))
  expect_equal(names(out), paste0("sim_", 1:4))
  # The seeds of the phylogenies obtained skip the failed seeds
  expect_equal(attr(out, "seeds"), setNames(seeds[c(1, 4, 5, 6)], names(out)))
  expect_equal(attr(out, "failed")$seed, failing)

  # Up to 10 x the number of simulations: here all the seeds after the first fail, so 1 of 3 simulations is obtained after 30 attempts
  expect_warning(out2 <- run_simulations(fake_simulation(failing_seeds = seeds[-1]), 3, TRUE, seed = 7),
                 "Only 1 of the 3 simulations requested were obtained: 29 simulation\\(s\\) failed after 30 attempts")
  expect_equal(names(out2), "sim_1")
  expect_equal(nrow(attr(out2, "failed")), 29)
})

test_that("run_simulations() stops when no simulation is obtained, and stops at once on a generic error",
{
  seeds <- derive_simulation_seeds(7, 3)
  expect_error(run_simulations(fake_simulation(failing_seeds = seeds), 3, FALSE, seed = 7), "None of the 3 simulation\\(s\\) could be obtained.*Reasons:.*failed for seed")
  # With a single simulation the error is raised as it is
  err <- tryCatch(run_simulations(fake_simulation(failing_seeds = 7), 1, FALSE, seed = 7), error = function (e) e)
  expect_s3_class(err, "phylograftR_simulation_error")
  expect_equal(conditionMessage(err), "failed for seed 7")

  # A generic error stops the whole run, even if other simulations could succeed
  calls <- 0
  counting <- function (sim_seed, sim_index, verbose_sim) { calls <<- calls + 1; fake_simulation(generic_fail_seeds = seeds[2])(sim_seed, sim_index, verbose_sim) }
  expect_error(run_simulations(counting, 3, FALSE, seed = 7), "generic problem")
  expect_equal(calls, 2)
  expect_error(run_simulations(counting, 3, TRUE, seed = 7), "generic problem")
})

test_that("run_simulations() passes 'verbose' to the simulations as it is, whatever the number of simulations, and prints nothing itself",
{
  seen <- NULL
  recorder <- function (sim_seed, sim_index, verbose_sim) { seen <<- c(seen, verbose_sim); fake_simulation()(sim_seed, sim_index, verbose_sim) }
  expect_silent(run_simulations(recorder, 1, FALSE, seed = 7, verbose = TRUE))
  expect_silent(run_simulations(recorder, 2, FALSE, seed = 7, verbose = TRUE))
  expect_silent(run_simulations(recorder, 2, FALSE, seed = 7, verbose = FALSE))
  expect_equal(seen, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})
