
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
#   * Parallel computation: the simulations (of all the maps) are independent tasks that can be run with the 'future' package (see `run_attempts()`).
#     The seed of a simulation only depends on its index, so the results are the same as when run sequentially.
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
#' @description With several simulations (or several maps), the temporary phylogenies are saved in one file per simulation (and map), to avoid overwriting:
#'   `"_map_<name>"` (name of the map, with the characters other than letters, digits, `_`, `.` and `-` replaced by `_`) and `"_sim_k"` are added before the extension `".rds"`.
#'
#' @param temp_file_path `NULL` (default path `"./temp_phy_from_grafting.rds"`), or path ending with `".rds"`.
#' @param index Integer. Index of the simulation.
#' @param multiple Logical. Whether several simulations are run. If `FALSE`, no suffix for the simulation is added.
#' @param map_name `NULL` (a single map), or name of the map. The suffix of the map is added if not `NULL`.
#' @param map_tag `NULL`, or tag of the map used in the suffix (see `temp_map_tags()`). By default, the name of the map with the characters other than letters, digits, `_`, `.` and `-` replaced by `_`.
#'
#' @noRd
#'
#' @return A path (or `NULL` if there is no suffix to add and `temp_file_path` is `NULL`).
#'

temp_file_path_simulation <- function (temp_file_path, index, multiple, map_name = NULL, map_tag = NULL)
{
  if (!multiple && is.null(map_name)) { return(temp_file_path) }
  if (is.null(temp_file_path)) { temp_file_path <- "./temp_phy_from_grafting.rds" }
  # The tag of the map is computed once for all the maps (see temp_map_tags()), so maps never share a file; otherwise from the name
  if (is.null(map_tag) && !is.null(map_name)) { map_tag <- gsub("[^A-Za-z0-9_.-]", "_", map_name) }
  suffix <- paste0(if (!is.null(map_name)) paste0("_map_", map_tag),
                   if (multiple) paste0("_sim_", index), collapse = "")
  sub("\\.rds$", paste0(suffix, ".rds"), temp_file_path)
}


#' @title Check the `parallel` argument
#'
#' @description Stops if `parallel` is not `NULL`, `TRUE` or `FALSE`.
#'
#' @param parallel `NULL`, `TRUE` or `FALSE`.
#'
#' @noRd
#'
#' @return `NULL` (invisibly).
#'

check_parallel <- function (parallel)
{
  if (!is.null(parallel) && (!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)))
  {
    stop("'parallel' must be NULL (use the 'future' plan if one is set), TRUE or FALSE.")
  }
  invisible(NULL)
}


#' @title Are the packages used for the parallel computation available?
#'
#' @description The packages `future` and `future.apply` are only suggested. A function of its own, so tests can mock it.
#'
#' @noRd
#'
#' @return Logical.
#'

have_parallel_packages <- function ()
{
  requireNamespace("future", quietly = TRUE) && requireNamespace("future.apply", quietly = TRUE)
}


#' @title Decide whether the tasks are run in parallel
#'
#' @description Rules (see the `parallel` argument of the grafting functions):
#'  * `FALSE`: sequential.
#'  * `NULL`: parallel if the packages `future` and `future.apply` are available, a `future` plan with several workers is set, and there are several tasks. Otherwise sequential, silently.
#'  * `TRUE`: same as `NULL`, but stops if the packages are missing or if no plan is set (nothing to run on), and warns if there is a single task (run sequentially).
#'
#' @param parallel `NULL`, `TRUE` or `FALSE`.
#' @param n_tasks Integer. Number of tasks (simulations) to run.
#'
#' @noRd
#'
#' @return Logical. Whether to run the tasks with `future`.
#'

resolve_parallel <- function (parallel, n_tasks)
{
  if (isFALSE(parallel)) { return(FALSE) }
  if (!have_parallel_packages())
  {
    if (isTRUE(parallel))
    {
      stop("'parallel = TRUE' requires the packages 'future' and 'future.apply', which are not installed: install.packages(c(\"future\", \"future.apply\")). ",
           "Use 'parallel = FALSE' to run sequentially.", call. = FALSE)
    }
    return(FALSE)
  }
  if (inherits(future::plan(), "sequential"))
  {
    if (isTRUE(parallel))
    {
      stop("'parallel = TRUE' needs a 'future' plan with several workers, set before the call, e.g., future::plan(future::multisession, workers = 4). ",
           "Use 'parallel = FALSE' to run sequentially.", call. = FALSE)
    }
    return(FALSE)
  }
  if (n_tasks < 2)
  {
    if (isTRUE(parallel)) { warning("'parallel = TRUE' has no effect with a single simulation to run: it is run sequentially.", call. = FALSE) }
    return(FALSE)
  }
  TRUE
}


#' @title Apply a function to a list of tasks, sequentially or with `future`
#'
#' @description `lapply()`, or `future.apply::future_lapply()` (with the plan set by the user), which returns the same list. The messages, the warnings and the printed output
#'   of the workers are relayed by `future` to the main session when the tasks are collected. Each task sets its own seed (with `set_local_seed()`, which restores the random number generator),
#'   so `future` must not manage the random numbers (`future.seed = FALSE`), the results are then the same as when run sequentially.
#'   Without seed (`seeded = FALSE`), `future` provides parallel-safe random number streams.
#'
#'   The globals are not scanned (`future.globals = FALSE`): the function is sent with its environment, and the internal functions of the package are found by the workers in the installed package
#'   (`future.packages`). Otherwise `future` would export each internal function used by the simulations as a separate object (about 50 messages per worker, about 1 s per worker).
#'
#' @param tasks List of tasks.
#' @param FUN Function applied to each task.
#' @param use_future Logical. Whether to use `future`.
#' @param seeded Logical. Whether each task sets its own seed.
#'
#' @noRd
#'
#' @return List with the result of each task, in the order of `tasks`.
#'

apply_tasks <- function (tasks, FUN, use_future, seeded)
{
  if (!use_future || length(tasks) < 2) { return(lapply(tasks, FUN)) }
  future.apply::future_lapply(tasks, FUN, future.seed = if (seeded) FALSE else TRUE, future.globals = FALSE, future.packages = "phylograftR")
}


#' @title Run the attempts (simulations) of one or several maps
#'
#' @description Runs the simulations of each map as ONE flat set of independent tasks (a task is a map x an attempt), sequentially or in parallel (see `apply_tasks()`).
#'   The seed of an attempt only depends on the index of the attempt, so the results are the same whatever the way the tasks are run.
#'   With `force_simulations`, the attempts are run in rounds: first `nb_simulations` attempts per map, then as many new attempts (with the next seeds) as simulations missing,
#'   and so on up to `10 * nb_simulations` attempts per map. A round only runs the attempts that the sequential rule (try new seeds until `nb_simulations` are obtained) would run, so the result is the same.
#'
#'   A task catches the failures specific to a simulation (error of class `"phylograftR_simulation_error"`, returned as a value), and gathers the warnings of the simulation if `collect` is `TRUE`
#'   (returned as a character vector, to be sent once by `summarise_simulations()`). The errors and the warnings that are not gathered carry the name of the map (see `with_map_context()`).
#'   In parallel, the output and the messages printed by the simulations are relayed by `future` when the tasks are collected, and each task prints one line when it is done.
#'
#' @param simulate_attempt Function with 3 arguments: `k` (index of the map), `sim_seed` and `sim_index`. Returns one grafted phylogeny.
#' @param n_maps Integer. Number of maps.
#' @param map_names `NULL` (a single map, or no map), or names of the maps (used in the messages).
#' @param nb_simulations,force_simulations See `run_simulations()`.
#' @param seeds `NULL`, or vector with the seed of each attempt (see `derive_simulation_seeds()`).
#' @param collect Logical. Whether to gather the warnings of the simulations.
#' @param parallel,verbose See `resolve_parallel()` and `run_simulations()`.
#'
#' @noRd
#'
#' @return A list with one element per map: a list with the elements `results` (named list of the phylogenies obtained), `failures` (list of failed attempts, with `seed`, `reason` and `condition`),
#'   `warnings` (warnings gathered) and `attempt` (index of the last attempt run).
#'

run_attempts <- function (simulate_attempt, n_maps, map_names, nb_simulations, force_simulations, seeds, collect, parallel, verbose)
{
  max_attempts <- if (force_simulations) 10 * nb_simulations else nb_simulations
  use_future <- resolve_parallel(parallel, n_tasks = n_maps * nb_simulations)
  if (use_future && verbose)
  {
    cat(paste0(Sys.time(), " - Running ", n_maps * nb_simulations, " simulation(s) in parallel on ", future::nbrOfWorkers(), " worker(s) (plan: ", setdiff(class(future::plan()), c("FutureStrategy", "tweaked"))[1], ")\n"))
  }

  run_task <- function (task)
  {
    k <- task$map
    attempt <- task$attempt
    map_name <- if (is.null(map_names)) NULL else map_names[k]
    sim_seed <- if (is.null(seeds)) NULL else seeds[attempt]
    if (verbose && !use_future && !is.null(map_name) && attempt == 1) { cat(paste0(Sys.time(), " - Map ", map_name, " (", k, "/", n_maps, ")\n")) }

    warnings_sent <- character(0)
    collect_warnings <- function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") }
    run <- function () simulate_attempt(k = k, sim_seed = sim_seed, sim_index = attempt)
    out <- with_map_context(
      if (collect) { tryCatch(withCallingHandlers(run(), warning = collect_warnings), phylograftR_simulation_error = function (e) e)
      } else { tryCatch(run(), phylograftR_simulation_error = function (e) e) },
      map_name)

    if (verbose && use_future)
    {
      cat(paste0(Sys.time(), " - Simulation ", task$id, "/", task$total, " done (", if (!is.null(map_name)) paste0("map ", map_name, ", "), "simulation n°", attempt, ")\n"))
    }
    list(out = out, warnings = warnings_sent, seed = sim_seed)
  }

  state <- lapply(seq_len(n_maps), function (k) list(results = list(), failures = list(), warnings = character(0), attempt = 0))
  repeat
  {
    # Attempts to run in this round: the ones still needed by each map
    tasks <- list()
    for (k in seq_len(n_maps))
    {
      n_run <- min(nb_simulations - length(state[[k]]$results), max_attempts - state[[k]]$attempt)
      if (n_run > 0) { tasks <- c(tasks, lapply(state[[k]]$attempt + seq_len(n_run), function (a) list(map = k, attempt = a))) }
    }
    if (length(tasks) == 0) { break }
    for (i in seq_along(tasks)) { tasks[[i]]$id <- i; tasks[[i]]$total <- length(tasks) }

    outcomes <- apply_tasks(tasks, run_task, use_future = use_future, seeded = !is.null(seeds))

    # Gather the outcomes, in the order of the attempts of each map
    for (i in seq_along(tasks))
    {
      k <- tasks[[i]]$map
      outcome <- outcomes[[i]]
      s <- state[[k]]
      s$attempt <- tasks[[i]]$attempt
      s$warnings <- c(s$warnings, outcome$warnings)
      if (inherits(outcome$out, "phylograftR_simulation_error"))
      {
        s$failures[[length(s$failures) + 1]] <- list(seed = if (is.null(outcome$seed)) NA else outcome$seed, reason = conditionMessage(outcome$out), condition = outcome$out)
      } else {
        phy <- outcome$out
        # The seed of the simulation is kept with the phylogeny
        if (!is.null(outcome$seed)) { attr(phy, "seed") <- as.integer(outcome$seed) }
        s$results[[length(s$results) + 1]] <- phy
        names(s$results)[length(s$results)] <- paste0("sim_", if (force_simulations) length(s$results) else s$attempt)
      }
      state[[k]] <- s
    }
  }
  state
}


#' @title Summarise the simulations obtained for a map
#'
#' @description Sends the warnings of the simulations (each distinct message once), and the warning or the error if fewer simulations than requested were obtained,
#'   and builds the output of a map (see `run_simulations()`).
#'
#' @param state Element of the output of `run_attempts()` for the map.
#' @param nb_simulations,force_simulations,seeds,allow_empty See `run_simulations()` and `run_attempts()`.
#'
#' @noRd
#'
#' @return See `run_simulations()`.
#'

summarise_simulations <- function (state, nb_simulations, force_simulations, seeds, allow_empty)
{
  several <- nb_simulations > 1
  results <- state$results
  failures <- state$failures
  attempt <- state$attempt

  # Warnings of the simulations, once
  for (message in unique(state$warnings)) { warning(message, call. = FALSE) }

  failed <- data.frame(seed = vapply(failures, function (f) as.numeric(f$seed), numeric(1)),
                       reason = vapply(failures, function (f) f$reason, character(1)), stringsAsFactors = FALSE)

  ## No simulation obtained
  if (length(results) == 0)
  {
    if (allow_empty) { return(structure(list(), class = "phylograftR_no_simulation", failed = failed)) }
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


#' @title Run the simulations of a grafting function
#'
#' @description Run a simulation function once per simulation, with the seeds of `derive_simulation_seeds()`, and gather the grafted phylogenies.
#'
#'   With one simulation (and without `force_simulations`), the phylogeny is returned as it is and every error is raised as it is.
#'   With several simulations, the warnings sent by the simulations are gathered and sent once at the end, and a failure specific to a simulation
#'   (error of class `"phylograftR_simulation_error"`, see `stop_simulation_failure()`) skips the simulation.
#'   With `force_simulations`, new seeds are tried until `nb_simulations` phylogenies are obtained, up to `10 * nb_simulations` attempts.
#'   The simulations can be run in parallel (see `run_attempts()`).
#'
#' @param simulate_one Function with 3 arguments: `sim_seed` (seed of the simulation, `NULL` to use the random number generator of the session), `sim_index` (index of the attempt)
#'   and `verbose_sim` (whether the simulation displays its progress). Returns one grafted phylogeny.
#' @param nb_simulations Integer. Number of simulations requested.
#' @param force_simulations Logical. Whether to try new seeds until `nb_simulations` simulations are obtained.
#' @param seed `NULL`, or integer. Seed of the run.
#' @param allow_empty Logical. If `TRUE`, no simulation obtained does not stop the function: an empty object of class `"phylograftR_no_simulation"` with the attribute `"failed"` is returned (used to drop a map, see `run_maps()`).
#' @param verbose Logical. Whether the simulations display their progress (passed to `simulate_one()` as `verbose_sim`).
#' @param parallel `NULL`, `TRUE` or `FALSE`. Whether to run the simulations in parallel with `future` (see `resolve_parallel()`). Default is `FALSE`.
#'
#' @noRd
#'
#' @return With `nb_simulations = 1`, a phylogeny (with the attribute `"seed"`, its seed, unless `seed` is `NULL`). Otherwise a `"multiPhylo"` object with the phylogenies named `"sim_1"`, `"sim_2"`, ...
#'   (the index of the simulation, so some can be missing if simulations failed, or consecutive numbers with `force_simulations`), with the attributes:
#'  * `"seeds"`: named vector with the seeds of the phylogenies (`NULL` if `seed` is `NULL`).
#'  * `"failed"`: data.frame with the columns `seed` and `reason` of the simulations that failed (0 row if none).
#'

run_simulations <- function (simulate_one, nb_simulations, force_simulations, seed, verbose = FALSE, allow_empty = FALSE, parallel = FALSE)
{
  max_attempts <- if (force_simulations) 10 * nb_simulations else nb_simulations
  seeds <- derive_simulation_seeds(seed, max_attempts)
  simulate_attempt <- function (k, sim_seed, sim_index) simulate_one(sim_seed = sim_seed, sim_index = sim_index, verbose_sim = verbose)

  state <- run_attempts(simulate_attempt = simulate_attempt, n_maps = 1, map_names = NULL, nb_simulations = nb_simulations, force_simulations = force_simulations,
                        seeds = seeds, collect = nb_simulations > 1 || force_simulations, parallel = parallel, verbose = verbose)
  summarise_simulations(state[[1]], nb_simulations = nb_simulations, force_simulations = force_simulations, seeds = seeds, allow_empty = allow_empty)
}


#' @title Check the `simmap` argument of the grafting functions
#'
#' @description Accepts `NULL`, a single stochastic character map (`"simmap"`), a set of maps of the same trait (`"multiSimmap"`), or a plain list of the maps of several traits
#'   (one element per trait, each a `"simmap"` or a `"multiSimmap"`; see `check_trait_maps()`). The maps are only listed here: each map is validated against the backbone
#'   when it is prepared (see `utils_simmap.R`).
#'
#' @param simmap `NULL`, a `"simmap"`, a `"multiSimmap"` or a list of them (several traits).
#'
#' @noRd
#'
#' @return A list with the elements `maps` (list of maps, `list(NULL)` without `simmap`; with several traits, each map is a list of the simmaps of the traits, class `"phylograftR_trait_maps"`),
#'   `names` (names of the maps: the names of the `"multiSimmap"` if it has some, otherwise `"map_1"`, `"map_2"`, ...; `NULL` for a single map or no map),
#'   `multiple` (`TRUE` for a `"multiSimmap"`, even with one map: the output is then a list with one element per map) and `traits` (names of the traits, `NULL` for a single trait).
#'

check_simmap_input <- function (simmap)
{
  if (is.null(simmap)) { return(list(maps = list(NULL), names = NULL, multiple = FALSE, traits = NULL)) }
  if (inherits(simmap, "multiSimmap"))
  {
    got <- get_multisimmap_maps(simmap)
    map_names <- if (is.null(got$names)) paste0("map_", seq_along(got$maps)) else got$names
    return(list(maps = got$maps, names = map_names, multiple = TRUE, traits = NULL))
  }
  if (inherits(simmap, "simmap")) { return(list(maps = list(simmap), names = NULL, multiple = FALSE, traits = NULL)) }
  if (is.list(simmap) && !inherits(simmap, "phylo")) { return(check_trait_maps(simmap)) }
  stop("'simmap' must be a stochastic character map of class 'simmap', as returned by phytools::make.simmap() or phytools::sim.history(), a set of maps of the same trait of class 'multiSimmap', ",
       "or a list of the maps of several traits.")
}


#' @title Add the name of the map to the messages of a step of the preparation
#'
#' @description Evaluates `expr`. If `map_name` is not `NULL`, the errors and the warnings of `expr` are sent again with the message starting with the name of the map,
#'   to know which map is concerned when several maps are used.
#'
#' @param expr Expression to evaluate.
#' @param map_name `NULL`, or name of the map.
#'
#' @noRd
#'
#' @return The value of `expr`.
#'

with_map_context <- function (expr, map_name)
{
  if (is.null(map_name)) { return(expr) }
  tryCatch(
    withCallingHandlers(expr, warning = function (w) { warning("Map '", map_name, "': ", conditionMessage(w), call. = FALSE); invokeRestart("muffleWarning") }),
    error = function (e) stop("Map '", map_name, "': ", conditionMessage(e), call. = FALSE))
}


#' @title Run the simulations of each map
#'
#' @description Runs the simulations of all the maps with `run_attempts()` (one flat set of tasks: the maps x the simulations, run in parallel if requested), then builds the output of each map.
#'
#'   Without `"multiSimmap"` (`multiple = FALSE`), this is the same as `run_simulations()` on the single map (or no map).
#'   With a `"multiSimmap"`, each map has the same seeds (the maps are compared with the same random draws) and the result is a list with one element per map, named after the map:
#'   a phylogeny if `nb_simulations = 1`, otherwise a `"multiPhylo"` (as in `run_simulations()`) named `"<map>_sim_1"`, `"<map>_sim_2"`, ...
#'   A map with no simulation obtained is dropped with a warning, and the function stops if no map gives a simulation.
#'
#' @param simulate_one Function with 6 arguments: `tracker` and `missing_taxa_list` (of the map), `map_name` (`NULL` without `"multiSimmap"`), `sim_seed`, `sim_index` and `verbose_sim`.
#' @param prepared List with one element per map, each a list with the elements `tracker` and `missing_taxa_list`.
#' @param map_names Names of the maps (`NULL` if `multiple = FALSE`).
#' @param multiple Logical. Whether `simmap` was a `"multiSimmap"`.
#' @param nb_simulations,force_simulations,seed,verbose,parallel See `run_simulations()`.
#'
#' @noRd
#'
#' @return Without `"multiSimmap"`, the output of `run_simulations()`. Otherwise a named list with one element per map kept, each a phylogeny or a `"multiPhylo"`, with the attribute `"failed_maps"`
#'   (data.frame with the columns `map` and `reason` of the maps dropped, 0 row if none). Each phylogeny has the attributes `"seed"` and `"map"`.
#'

run_maps <- function (simulate_one, prepared, map_names, multiple, nb_simulations, force_simulations, seed, verbose = FALSE, parallel = FALSE)
{
  n_maps <- length(prepared)
  if (multiple && n_maps * nb_simulations > 1000)
  {
    warning(n_maps, " maps x ", nb_simulations, " simulations = ", n_maps * nb_simulations, " grafted phylogenies will be produced (more than 1000): ",
            "this can take a long time and use a lot of memory.", call. = FALSE)
  }

  max_attempts <- if (force_simulations) 10 * nb_simulations else nb_simulations
  seeds <- derive_simulation_seeds(seed, max_attempts)
  simulate_attempt <- function (k, sim_seed, sim_index)
  {
    simulate_one(tracker = prepared[[k]]$tracker, missing_taxa_list = prepared[[k]]$missing_taxa_list, map_name = if (multiple) map_names[k] else NULL,
                 sim_seed = sim_seed, sim_index = sim_index, verbose_sim = verbose)
  }
  state <- run_attempts(simulate_attempt = simulate_attempt, n_maps = n_maps, map_names = if (multiple) map_names else NULL, nb_simulations = nb_simulations,
                        force_simulations = force_simulations, seeds = seeds, collect = nb_simulations > 1 || force_simulations, parallel = parallel, verbose = verbose)
  if (!multiple)
  {
    return(summarise_simulations(state[[1]], nb_simulations = nb_simulations, force_simulations = force_simulations, seeds = seeds, allow_empty = FALSE))
  }

  results <- list()
  failed_maps <- list()
  for (k in seq_len(n_maps))
  {
    out <- with_map_context(summarise_simulations(state[[k]], nb_simulations = nb_simulations, force_simulations = force_simulations, seeds = seeds, allow_empty = TRUE), map_names[k])
    if (inherits(out, "phylograftR_no_simulation"))
    {
      reasons <- unique(attr(out, "failed")$reason)
      failed_maps[[length(failed_maps) + 1]] <- data.frame(map = map_names[k], reason = paste(reasons, collapse = "\n"), stringsAsFactors = FALSE)
      warning("Map '", map_names[k], "': none of the ", nb_simulations, " simulation(s) could be obtained, the map is dropped from the output.\n",
              "Reasons:\n", paste0("  - ", reasons, collapse = "\n"), call. = FALSE)
      next
    }
    if (inherits(out, "multiPhylo"))
    {
      for (j in seq_along(out)) { attr(out[[j]], "map") <- map_names[k] }
      names(out) <- paste0(map_names[k], "_", names(out))
      if (!is.null(attr(out, "seeds"))) { names(attr(out, "seeds")) <- names(out) }
    } else {
      attr(out, "map") <- map_names[k]
    }
    results[[map_names[k]]] <- out
  }
  failed_maps <- if (length(failed_maps) > 0) do.call(rbind, failed_maps) else data.frame(map = character(0), reason = character(0), stringsAsFactors = FALSE)

  if (length(results) == 0)
  {
    stop("None of the ", n_maps, " maps gave a simulation.\nReasons:\n",
         paste0("  - ", failed_maps$map, ": ", sub("\n.*", "", failed_maps$reason), collapse = "\n"))
  }
  structure(results, failed_maps = failed_maps)
}


#' @title Evaluate an expression and send its warnings once
#'
#' @description Evaluates `expr`, gathers its warnings, and sends each distinct message once at the end. Used when the same step is applied to several phylogenies
#'   (e.g., the same warning about the inputs would otherwise be sent for each phylogeny).
#'
#' @param expr Expression to evaluate.
#'
#' @noRd
#'
#' @return The value of `expr`.
#'

with_unique_warnings <- function (expr)
{
  messages <- character(0)
  value <- withCallingHandlers(expr, warning = function (w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
  for (message in unique(messages)) { warning(message, call. = FALSE) }
  value
}


#' @title Describe the structure of grafted phylogenies given to a function
#'
#' @description Lists the phylogenies of the output of the grafting functions: a phylogeny, a `"multiPhylo"` (the simulations of one map), or a list of them (the maps of a `"multiSimmap"`).
#'
#' @param x A `"phylo"`, a `"multiPhylo"`, or a plain list (named after the maps) of `"phylo"` and `"multiPhylo"` objects.
#' @param arg Name of the argument, for the error messages.
#'
#' @noRd
#'
#' @return A list with the elements `phylos` (list of phylogenies, in the order of the input), `map` (name of the map of each phylogeny), `simulation` (integer index of the simulation of each phylogeny:
#'   the index in the name `"sim_k"` or `"<map>_sim_k"` of a `"multiPhylo"`, otherwise its position, 1 for a single phylogeny), `seed` (integer, `NA` if unknown) and `is_list` (whether `x` is a list of maps).
#'   The map is the name of the list of maps, or `"map_1"` if the phylogeny has a single `$simmap` without name, otherwise `NA` (no map was used).
#'

list_grafted_phylogenies <- function (x, arg = "phy_grafted")
{
  phylos <- list(); map <- character(0); simulation <- integer(0)
  index_of <- function (names_x, n) {
    if (is.null(names_x)) { return(seq_len(n)) }
    idx <- suppressWarnings(as.integer(sub("^(.*_)?sim_([0-9]+)$", "\\2", names_x)))
    ifelse(is.na(idx) | !grepl("sim_[0-9]+$", names_x), seq_len(n), idx)
  }
  add <- function (phy, default_map, forced_map = NULL, index = 1L)
  {
    phylos[[length(phylos) + 1]] <<- phy
    # Name of the map: the name in the list of maps, otherwise the attribute "map" of the phylogeny, otherwise "map_1" if the phylogeny has the simmap of a single map, otherwise NA (no map used)
    map <<- c(map, if (!is.null(forced_map)) forced_map else if (!is.null(attr(phy, "map"))) attr(phy, "map") else if (!is.null(phy$simmap)) default_map else NA_character_)
    simulation <<- c(simulation, index)
  }
  bad <- function () stop("'", arg, "' must be a phylogeny (class 'phylo'), a set of phylogenies (class 'multiPhylo'), or a list of them with one element per map (output of the grafting functions with a 'multiSimmap').")

  if (inherits(x, "phylo")) {
    add(x, "map_1")
    is_list <- FALSE
  } else if (inherits(x, "multiPhylo")) {
    idx <- index_of(names(x), length(x))
    for (j in seq_along(x)) { add(x[[j]], "map_1", index = idx[j]) }
    is_list <- FALSE
  } else if (is.list(x) && length(x) > 0) {
    map_names <- names(x)
    if (is.null(map_names) || anyNA(map_names) || any(map_names == "")) { map_names <- paste0("map_", seq_along(x)) }
    for (k in seq_along(x))
    {
      el <- x[[k]]
      if (inherits(el, "phylo")) { add(el, map_names[k], forced_map = map_names[k]) }
      else if (inherits(el, "multiPhylo")) {
        idx <- index_of(names(el), length(el))
        for (j in seq_along(el)) { add(el[[j]], map_names[k], forced_map = map_names[k], index = idx[j]) }
      } else { bad() }
    }
    is_list <- TRUE
  } else { bad() }
  if (length(phylos) == 0) { stop("'", arg, "' does not contain any phylogeny.") }
  seed <- vapply(phylos, function (p) { s <- attr(p, "seed"); if (is.null(s)) NA_integer_ else as.integer(s) }, integer(1))
  list(phylos = phylos, map = map, simulation = simulation, seed = seed, is_list = is_list)
}
