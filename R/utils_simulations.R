
##### Internal helpers to run several grafting simulations in the same run #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Helper functions to produce several simulations (grafted phylogenies) from the same backbone and the same list of missing taxa #####

# Rationale
#   The grafting is stochastic, so analyses should be run on a set of grafted phylogenies rather than on a single one.
#   The grafting functions validate their inputs once, then run 'nb_simulations' simulations (see `run_simulations()`).
#   * Seeds: the first simulation uses the 'seed' argument itself (it is the same phylogeny as a run with a single simulation), and the seeds of the other
#     simulations are drawn from a stream started by 'seed' (see `derive_simulation_seeds()`). Each simulation can thus be reproduced alone with its seed,
#     and the first simulations of a longer run are the same as the ones of a shorter run.
#   * Failures: a failure specific to a simulation (stochastic: raised with `stop_simulation_failure()`) skips this simulation with a warning,
#     whereas any other error (inputs, validation, a problem that would affect all the simulations) stops the whole run.


#' @title Check the arguments to run several simulations
#'
#' @description Stops if `nb_simulations` is not a single positive integer, or if `force_simulations` is not `TRUE` or `FALSE`.
#'
#' @param nb_simulations Number of simulations (grafted phylogenies) requested.
#' @param force_simulations Whether to try new seeds until `nb_simulations` simulations are obtained.
#'
#' @noRd
#'
#' @return `NULL` (invisibly).
#'

check_nb_simulations <- function (nb_simulations, force_simulations)
{
  if (!is.numeric(nb_simulations) || length(nb_simulations) != 1 || is.na(nb_simulations) || nb_simulations < 1 || nb_simulations != round(nb_simulations))
  {
    stop("'nb_simulations' must be a single positive integer.")
  }
  if (!is.logical(force_simulations) || length(force_simulations) != 1 || is.na(force_simulations))
  {
    stop("'force_simulations' must be TRUE or FALSE.")
  }
  invisible(NULL)
}


#' @title Stop with a failure that is specific to a simulation
#'
#' @description Stop with an error of class `"phylograftR_simulation_error"`. When several simulations are run, such an error skips the simulation (see `run_simulations()`),
#'   whereas any other error stops the whole run. Use it for failures due to the random draws of a simulation (e.g., no usable divergence age after the maximum number of redraws).
#'
#' @param ... Pieces of the message, pasted together (as in `stop()`).
#'
#' @noRd
#'
#' @return Does not return.
#'

stop_simulation_failure <- function (...)
{
  stop(structure(class = c("phylograftR_simulation_error", "error", "condition"),
                 list(message = paste0(...), call = NULL)))
}


#' @title Derive the seeds of the simulations from the seed of the run
#'
#' @description The first seed is `seed` itself (so the first simulation is the same as a run with a single simulation). The other seeds are drawn with `sample.int()`
#'   after `set.seed(seed)`, without altering the random number generator of the session. The seeds of a shorter run are the first seeds of a longer run.
#'
#' @param seed `NULL`, or integer. Seed of the run.
#' @param n Integer. Number of seeds.
#'
#' @noRd
#'
#' @return `NULL` if `seed` is `NULL` (the simulations then follow the random number generator of the session), otherwise a vector of `n` distinct seeds.
#'

derive_simulation_seeds <- function (seed, n)
{
  if (is.null(seed)) { return(NULL) }
  if (n == 1) { return(seed) }
  set_local_seed(seed)
  others <- sample.int(.Machine$integer.max, size = n)
  others <- others[others != seed]
  c(seed, others[seq_len(n - 1)])
}


#' @title Path of the temporary file of a simulation
#'
#' @description With several simulations, the temporary phylogenies are saved in one file per simulation, to avoid overwriting: `"_sim_k"` is added before the extension `".rds"`.
#'
#' @param temp_file_path `NULL` (default path `"./temp_phy_from_grafting.rds"`), or path ending with `".rds"`.
#' @param index Integer. Index of the simulation.
#' @param multiple Logical. Whether several simulations are run. If `FALSE`, `temp_file_path` is returned as it is.
#'
#' @noRd
#'
#' @return A path (or `NULL` if `multiple = FALSE` and `temp_file_path` is `NULL`).
#'

temp_file_path_simulation <- function (temp_file_path, index, multiple)
{
  if (!multiple) { return(temp_file_path) }
  if (is.null(temp_file_path)) { temp_file_path <- "./temp_phy_from_grafting.rds" }
  sub("\\.rds$", paste0("_sim_", index, ".rds"), temp_file_path)
}


#' @title Run the simulations of a grafting function
#'
#' @description Run a simulation function once per simulation, with the seeds of `derive_simulation_seeds()`, and gather the grafted phylogenies.
#'
#'   With one simulation (and without `force_simulations`), the phylogeny is returned as it is and every error is raised as it is.
#'   With several simulations, the warnings sent by the simulations are gathered and sent once at the end, and a failure specific to a simulation
#'   (error of class `"phylograftR_simulation_error"`, see `stop_simulation_failure()`) skips the simulation.
#'   With `force_simulations`, new seeds are tried until `nb_simulations` phylogenies are obtained, up to `10 * nb_simulations` attempts.
#'
#' @param simulate_one Function with 3 arguments: `sim_seed` (seed of the simulation, `NULL` to use the random number generator of the session), `sim_index` (index of the attempt)
#'   and `verbose_sim` (whether the simulation displays its progress). Returns one grafted phylogeny.
#' @param nb_simulations Integer. Number of simulations requested.
#' @param force_simulations Logical. Whether to try new seeds until `nb_simulations` simulations are obtained.
#' @param seed `NULL`, or integer. Seed of the run.
#' @param verbose Logical. Whether to display the progress.
#'
#' @noRd
#'
#' @return With `nb_simulations = 1`, a phylogeny. Otherwise a `"multiPhylo"` object with the phylogenies named `"sim_1"`, `"sim_2"`, ...
#'   (the index of the simulation, so some can be missing if simulations failed, or consecutive numbers with `force_simulations`), with the attributes:
#'  * `"seeds"`: named vector with the seeds of the phylogenies (`NULL` if `seed` is `NULL`).
#'  * `"failed"`: data.frame with the columns `seed` and `reason` of the simulations that failed (0 row if none).
#'

run_simulations <- function (simulate_one, nb_simulations, force_simulations, seed, verbose = FALSE)
{
  several <- nb_simulations > 1
  max_attempts <- if (force_simulations) 10 * nb_simulations else nb_simulations
  seeds <- derive_simulation_seeds(seed, max_attempts)
  verbose_sim <- verbose && !several

  results <- list()
  failures <- list()
  warnings_sent <- character(0)
  attempt <- 0
  while (length(results) < nb_simulations && attempt < max_attempts)
  {
    attempt <- attempt + 1
    sim_seed <- if (is.null(seeds)) NULL else seeds[attempt]
    if (verbose && several) { cat(paste0(Sys.time(), " - Simulation n\u00b0", attempt, " (", length(results), "/", nb_simulations, " obtained)\n")) }

    run <- function () simulate_one(sim_seed = sim_seed, sim_index = attempt, verbose_sim = verbose_sim)
    if (several || force_simulations)
    {
      # Warnings are gathered, to be sent once at the end
      collect_warnings <- function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") }
      out <- tryCatch(withCallingHandlers(run(), warning = collect_warnings), phylograftR_simulation_error = function (e) e)
    } else {
      out <- tryCatch(run(), phylograftR_simulation_error = function (e) e)
    }

    if (inherits(out, "phylograftR_simulation_error"))
    {
      failures[[length(failures) + 1]] <- list(seed = if (is.null(sim_seed)) NA else sim_seed, reason = conditionMessage(out), condition = out)
    } else {
      results[[length(results) + 1]] <- out
      names(results)[length(results)] <- paste0("sim_", if (force_simulations) length(results) else attempt)
    }
  }

  # Warnings of the simulations, once
  for (message in unique(warnings_sent)) { warning(message, call. = FALSE) }

  failed <- data.frame(seed = vapply(failures, function (f) as.numeric(f$seed), numeric(1)),
                       reason = vapply(failures, function (f) f$reason, character(1)), stringsAsFactors = FALSE)

  ## No simulation obtained
  if (length(results) == 0)
  {
    # The error of the only attempt is raised as it is
    if (attempt == 1) { stop(failures[[1]]$condition) }
    stop("None of the ", attempt, " simulation(s) could be obtained", if (force_simulations) " (up to 10 times the number of simulations requested)" else "", ".\n",
         "Reasons:\n", paste0("  - ", unique(failed$reason), collapse = "\n"))
  }

  ## Fewer simulations than requested
  if (length(results) < nb_simulations)
  {
    n_show <- min(nrow(failed), 10)
    reasons <- unique(failed$reason)
    warning("Only ", length(results), " of the ", nb_simulations, " simulations requested were obtained: ", nrow(failed), " simulation(s) failed",
            if (force_simulations) paste0(" after ", attempt, " attempts (10 times the number of simulations requested)") else "",
            if (!force_simulations) " (use 'force_simulations = TRUE' to try new seeds until the number requested is reached)" else "", ".\n",
            "Seeds of the failed simulations: ", paste(utils::head(failed$seed, n_show), collapse = ", "), if (nrow(failed) > n_show) ", ..." else "", ".\n",
            "Reasons:\n", paste0("  - ", utils::head(reasons, 5), collapse = "\n"), if (length(reasons) > 5) "\n  - ..." else "", call. = FALSE)
  }

  if (!several) { return(results[[1]]) }

  seeds_obtained <- NULL
  if (!is.null(seeds))
  {
    used <- if (force_simulations) setdiff(seq_len(attempt), match(failed$seed, seeds)) else as.integer(sub("sim_", "", names(results)))
    seeds_obtained <- stats::setNames(seeds[used], names(results))
  }
  structure(results, class = "multiPhylo", seeds = seeds_obtained, failed = failed)
}
