
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
#' @description With several simulations (or several maps), the temporary phylogenies are saved in one file per simulation (and map), to avoid overwriting:
#'   `"_map_<name>"` (name of the map, with the characters other than letters, digits, `_`, `.` and `-` replaced by `_`) and `"_sim_k"` are added before the extension `".rds"`.
#'
#' @param temp_file_path `NULL` (default path `"./temp_phy_from_grafting.rds"`), or path ending with `".rds"`.
#' @param index Integer. Index of the simulation.
#' @param multiple Logical. Whether several simulations are run. If `FALSE`, no suffix for the simulation is added.
#' @param map_name `NULL` (a single map), or name of the map. The suffix of the map is added if not `NULL`.
#'
#' @noRd
#'
#' @return A path (or `NULL` if there is no suffix to add and `temp_file_path` is `NULL`).
#'

temp_file_path_simulation <- function (temp_file_path, index, multiple, map_name = NULL)
{
  if (!multiple && is.null(map_name)) { return(temp_file_path) }
  if (is.null(temp_file_path)) { temp_file_path <- "./temp_phy_from_grafting.rds" }
  suffix <- paste0(if (!is.null(map_name)) paste0("_map_", gsub("[^A-Za-z0-9_.-]", "_", map_name)),
                   if (multiple) paste0("_sim_", index), collapse = "")
  sub("\\.rds$", paste0(suffix, ".rds"), temp_file_path)
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
#' @param allow_empty Logical. If `TRUE`, no simulation obtained does not stop the function: an empty object of class `"phylograftR_no_simulation"` with the attribute `"failed"` is returned (used to drop a map, see `run_maps()`).
#' @param verbose Logical. Whether the simulations display their progress (passed to `simulate_one()` as `verbose_sim`).
#'
#' @noRd
#'
#' @return With `nb_simulations = 1`, a phylogeny (with the attribute `"seed"`, its seed, unless `seed` is `NULL`). Otherwise a `"multiPhylo"` object with the phylogenies named `"sim_1"`, `"sim_2"`, ...
#'   (the index of the simulation, so some can be missing if simulations failed, or consecutive numbers with `force_simulations`), with the attributes:
#'  * `"seeds"`: named vector with the seeds of the phylogenies (`NULL` if `seed` is `NULL`).
#'  * `"failed"`: data.frame with the columns `seed` and `reason` of the simulations that failed (0 row if none).
#'

run_simulations <- function (simulate_one, nb_simulations, force_simulations, seed, verbose = FALSE, allow_empty = FALSE)
{
  several <- nb_simulations > 1
  max_attempts <- if (force_simulations) 10 * nb_simulations else nb_simulations
  seeds <- derive_simulation_seeds(seed, max_attempts)
  verbose_sim <- verbose

  results <- list()
  failures <- list()
  warnings_sent <- character(0)
  attempt <- 0
  while (length(results) < nb_simulations && attempt < max_attempts)
  {
    attempt <- attempt + 1
    sim_seed <- if (is.null(seeds)) NULL else seeds[attempt]

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
      # The seed of the simulation is kept with the phylogeny
      if (!is.null(sim_seed)) { attr(out, "seed") <- as.integer(sim_seed) }
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
#' @description Runs `run_simulations()` for each map prepared, and gathers the results.
#'
#'   Without `"multiSimmap"` (`multiple = FALSE`), this is `run_simulations()` on the single map (or no map).
#'   With a `"multiSimmap"`, each map has the same seeds (the maps are compared with the same random draws) and the result is a list with one element per map, named after the map:
#'   a phylogeny if `nb_simulations = 1`, otherwise a `"multiPhylo"` (as in `run_simulations()`) named `"<map>_sim_1"`, `"<map>_sim_2"`, ...
#'   A map with no simulation obtained is dropped with a warning, and the function stops if no map gives a simulation.
#'
#' @param simulate_one Function with 6 arguments: `tracker` and `missing_taxa_list` (of the map), `map_name` (`NULL` without `"multiSimmap"`), `sim_seed`, `sim_index` and `verbose_sim`.
#' @param prepared List with one element per map, each a list with the elements `tracker` and `missing_taxa_list`.
#' @param map_names Names of the maps (`NULL` if `multiple = FALSE`).
#' @param multiple Logical. Whether `simmap` was a `"multiSimmap"`.
#' @param nb_simulations,force_simulations,seed,verbose See `run_simulations()`.
#'
#' @noRd
#'
#' @return Without `"multiSimmap"`, the output of `run_simulations()`. Otherwise a named list with one element per map kept, each a phylogeny or a `"multiPhylo"`, with the attribute `"failed_maps"`
#'   (data.frame with the columns `map` and `reason` of the maps dropped, 0 row if none). Each phylogeny has the attributes `"seed"` and `"map"`.
#'

run_maps <- function (simulate_one, prepared, map_names, multiple, nb_simulations, force_simulations, seed, verbose = FALSE)
{
  run_map <- function (k)
  {
    name <- if (multiple) map_names[k] else NULL
    simulate_map <- function (sim_seed, sim_index, verbose_sim)
    {
      simulate_one(tracker = prepared[[k]]$tracker, missing_taxa_list = prepared[[k]]$missing_taxa_list, map_name = name,
                   sim_seed = sim_seed, sim_index = sim_index, verbose_sim = verbose_sim)
    }
    run_simulations(simulate_one = simulate_map, nb_simulations = nb_simulations, force_simulations = force_simulations, seed = seed, verbose = verbose, allow_empty = multiple)
  }
  if (!multiple) { return(run_map(1)) }

  n_maps <- length(prepared)
  if (n_maps * nb_simulations > 1000)
  {
    warning(n_maps, " maps x ", nb_simulations, " simulations = ", n_maps * nb_simulations, " grafted phylogenies will be produced (more than 1000): ",
            "this can take a long time and use a lot of memory.", call. = FALSE)
  }

  results <- list()
  failed_maps <- list()
  for (k in seq_len(n_maps))
  {
    if (verbose) { cat(paste0(Sys.time(), " - Map ", map_names[k], " (", k, "/", n_maps, ")\n")) }
    out <- with_map_context(run_map(k), map_names[k])
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
