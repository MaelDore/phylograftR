
##### Function to flatten the output of the grafting functions into a single set of phylogenies #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

## Rationale:
  # If the inputs of gratfing functions contain multiple stochastic maps (simmaps) for a single trait,
  # the output will be a list of phylogenies (or multiPhylo, if 'nb_simulations' > 1),
  # with each element associated to a unique simmap.
  # This function provides a tool for the user to 'flatten' this list of phylogenies into a unique multiPhylo object
  # were each phylogeny is named after is associated map and simulation index.

### 1/ Flatten the grafted phylogenies of several maps and simulations ####

#' @title Flatten grafted phylogenies into a single set of phylogenies
#'
#' @description The grafting functions return a list with one element per stochastic map when they are used with several maps (a `"multiSimmap"`),
#' each element being a phylogeny or a set of phylogenies (one per simulation). This function flattens such an output into a single object of class `"multiPhylo"`,
#' which can be used directly by the functions that expect a set of phylogenies under this class (e.g., from the packages `ape` or `phytools`).
#' The associated stochastic map ID, the simulation ID and the seed of each phylogeny are kept, so the phylogenies can still be traced back and selected.
#'
#' The arguments `maps` and `simulations` can be used to subset the object and retain only a specific sets of grafted phylogenies
#' based on their associated stochastic map, and/or their simulation index.
#'
#' Additionally, the function can trim out the updated stochastic maps stored in the `$simmap` element of each grafted function to save memory space
#' if those updated grafted stochastic maps are not required for downstream analyses.
#'
#' @param phy_grafted The output of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades]:
#'   a list with one element per stochastic map (each a `"phylo"` or a `"multiPhylo"`), a `"multiPhylo"`, or a `"phylo"`.
#'   A `"multiPhylo"` or a `"phylo"` object is already flat, but the function can still be used to select a subset of phylogenies and remove their associated simmaps if requested.
#'   A `"phylo"` is returned as a `"phylo"` (the only possible action is then to remove its `$simmap`).
#' @param maps Character vector. Names of the maps to keep (the names of the elements of the list, as in the column `map` of the attribute `"phylogenies"`).
#'   Default is `NULL`, which keeps all the maps. The function stops if some names are not found.
#' @param simulations Integer vector. Indices of the simulations to keep (the numbers in the names `"sim_k"` of the phylogenies, which are the same in each map:
#'   `1` for the first simulation of each map, etc.). Default is `NULL`, which keeps all the simulations.
#' @param keep_simmap Logical. Whether to keep the element `$simmap` (the updated simmap) of the phylogenies, if they have one. Set to `FALSE` to save memory
#'   when the stochastic maps are not needed in downstream analyses. Default is `TRUE`.
#'
#' @details The phylogenies are returned in the order of the input (map by map, then simulation by simulation), with their names as in the input:
#' the names of the maps for a list of phylogenies (one simulation per map), and `"<map>_sim_k"` for a list of sets of phylogenies (several simulations per map).
#' If two phylogenies have the same name, the names are made unique with a warning.
#'
#' Each phylogeny keeps its attributes `"map"` and `"seed"`. When `maps` and `simulations` are used, the selection is applied to all the phylogenies together:
#' the phylogenies kept are the ones from the maps in `maps` and with the simulation index in `simulations`.
#'
#' @return For a `"phylo"` as `phy_grafted`, the phylogeny (class `"phylo"`), with its `$simmap` removed if `keep_simmap = FALSE` (otherwise unchanged).
#'   Otherwise, an object of class `"multiPhylo"` with all the phylogenies kept, with the following attributes:
#'  * `"phylogenies"`: data.frame with one row per phylogeny, in the order of the object, with the columns `name`, `map` (`NA` if no map was used), `simulation` (index of the simulation) and `seed` (`NA` if unknown).
#'  * `"seeds"`: named integer vector with the seed of each phylogeny (`NULL` if the seeds are unknown, i.e., if the grafting was run with `seed = NULL`).
#'  * `"failed"`: data.frame with the simulations that failed in the maps kept, with the columns `map`, `seed` and `reason` (0 row if none). Not filtered by `simulations`.
#'  * `"failed_maps"`: data.frame with the maps dropped from the output of the grafting function because none of their simulations could be obtained, with the columns `map` and `reason` (0 row if none).
#'
#' @seealso Grafting functions that produce the output to flatten:
#'   * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade with constraints (target clades and excluded subclades),
#'   * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as clades sister to a clade of the backbone.
#'
#'   Functions that also accept several grafted phylogenies:
#'   * [`phylograftR::check_grafting_constraints()`][check_grafting_constraints] to check a posteriori that grafted taxa respect their constraints,
#'   * [`phylograftR::identify_grafted()`][identify_grafted] to identify the grafted tips and branches.
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Load backbone phylogeny, associated simmap, and missing taxa list
#' data(toy_backbone)
#' data(toy_simmap)
#' data(toy_missing_taxa_list_with_allowed_states)
#'
#' # Create a toy backbone phylogeny, with a set of 3 stochastic maps of a trait
#' tip_states <- phytools::getStates(toy_simmap, type = "tips")
#' had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
#' if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
#' set.seed(118)
#' toy_multiSimmap <- phytools::make.simmap(toy_backbone, tip_states, model = "ER", nsim = 3, message = FALSE)
#' if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
#'
#' # --- Graft missing taxa --- #
#'
#' # Graft with the 3 maps and 2 simulations per map as constraints:
#' # the output is a list with one element per map hosting 2 simulations each.
#' grafted_list <- graft_missing_taxa(phy = toy_backbone,
#'                               missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
#'                               simmap = toy_multiSimmap, allowed_states = c("blue", "green"),
#'                               method = "BD_clades", nb_simulations = 2, verbose = FALSE, seed = 1234)
#' names(grafted_list)
#' str(grafted_list, 2)
#'
#' # --- Flatten the output --- #
#'
#' # Flatten the output into a single multiPhylo with the 6 grafted phylogenies
#' grafted_phys <- flatten_grafted_phylogenies(grafted_list)
#' names(grafted_phys)
#' # The map, simulation and seed of each phylogeny are stored as an attribute
#' attr(grafted_phys, "phylogenies")
#'
#' # --- Flatten the output, subset, and remove simmaps all-at-once --- #
#'
#' # Keep the first simulation of two maps only, without the simmaps (to save memory)
#' grafted_subset <- flatten_grafted_phylogenies(grafted_list, maps = c("map_1", "map_3"),
#'                                               simulations = 1, keep_simmap = FALSE)
#' attr(grafted_subset, "phylogenies")
#'
#' # Use the set of phylogenies with functions that expect a multiPhylo, e.g., to compute the tree heights
#' sapply(grafted_phys, function (x) max(ape::node.depth.edgelength(x)))
#'
#'


flatten_grafted_phylogenies <- function (phy_grafted, # Output of the grafting functions: list of maps (phylo or multiPhylo), multiPhylo or phylo
                                         maps = NULL, # Names of the maps to keep. NULL: all
                                         simulations = NULL, # Indices of the simulations to keep. NULL: all
                                         keep_simmap = TRUE) # Keep the $simmap of the phylogenies?
{
  ### ---- Check input validity ---- ###
  {
    if (!is.null(maps) && (!is.character(maps) || length(maps) == 0 || anyNA(maps)))
    {
      stop("'maps' must be NULL or a character vector of map names.")
    }
    if (!is.null(simulations) && (!is.numeric(simulations) || length(simulations) == 0 || anyNA(simulations) || any(simulations < 1) || any(simulations != round(simulations))))
    {
      stop("'simulations' must be NULL or a vector of positive integers (indices of simulations).")
    }
    if (!is.logical(keep_simmap) || length(keep_simmap) != 1 || is.na(keep_simmap))
    {
      stop("'keep_simmap' must be TRUE or FALSE.")
    }
  }

  ### ---- List the phylogenies, with their names and the information on the failures ---- ###

  listed <- list_grafted_phylogenies(phy_grafted, arg = "phy_grafted")
  n_phylos <- length(listed$phylos)

  # Extract names of the phylogenies, as in the input (and names of the maps of a list)
  map_names <- NULL
  if (listed$is_list)
  {
    map_names <- names(phy_grafted)
    if (is.null(map_names) || anyNA(map_names) || any(map_names == "")) { map_names <- paste0("map_", seq_along(phy_grafted)) }
  }
  phylo_names <- if (inherits(phy_grafted, "phylo"))
  {
    if (is.na(listed$map[1])) "sim_1" else listed$map[1]
  } else if (inherits(phy_grafted, "multiPhylo")) {
    if (is.null(names(phy_grafted))) paste0("sim_", listed$simulation) else names(phy_grafted)
  } else {
    unlist(lapply(seq_along(phy_grafted), function (k)
    {
      el <- phy_grafted[[k]]
      if (inherits(el, "phylo")) { map_names[k] }
      else if (is.null(names(el))) { paste0(map_names[k], "_sim_", seq_along(el)) }
      else { names(el) }
    }), use.names = FALSE)
  }
  if (anyDuplicated(phylo_names) > 0)
  {
    warning("Some phylogenies have the same name: the names are made unique.", call. = FALSE)
    phylo_names <- make.unique(phylo_names)
  }

  # Record failed simulations and failed maps, with the map of each failed simulation
  empty_failed <- data.frame(map = character(0), seed = integer(0), reason = character(0), stringsAsFactors = FALSE)
  failed_of <- function (set, map)
  {
    f <- attr(set, "failed")
    if (is.null(f) || nrow(f) == 0) { return(empty_failed) }
    data.frame(map = rep(map, nrow(f)), seed = as.integer(f$seed), reason = f$reason, stringsAsFactors = FALSE)
  }
  if (listed$is_list) {
    failed <- do.call(rbind, c(list(empty_failed), lapply(seq_along(phy_grafted), function (k) failed_of(phy_grafted[[k]], map_names[k]))))
    failed_maps <- attr(phy_grafted, "failed_maps")
  } else {
    failed <- failed_of(phy_grafted, if (n_phylos > 0) listed$map[1] else NA_character_)
    failed_maps <- NULL
  }
  if (is.null(failed_maps)) { failed_maps <- data.frame(map = character(0), reason = character(0), stringsAsFactors = FALSE) }

  ### ---- Select the maps and the simulations ---- ###

  keep <- rep(TRUE, n_phylos)
  if (!is.null(maps))
  {
    available <- unique(listed$map[!is.na(listed$map)])
    if (length(available) == 0) { stop("'maps' cannot be used: the phylogenies were not grafted with a 'multiSimmap', so they do not have a map.") }
    unknown <- setdiff(maps, available)
    if (length(unknown) > 0)
    {
      stop("Some maps listed in 'maps' are not found: ", paste(unknown, collapse = ", "), ".\nAvailable maps: ", paste(available, collapse = ", "), ".")
    }
    keep <- keep & listed$map %in% maps
    failed <- failed[failed$map %in% maps, , drop = FALSE]
    failed_maps <- failed_maps[failed_maps$map %in% maps, , drop = FALSE]
  }
  if (!is.null(simulations))
  {
    keep <- keep & listed$simulation %in% simulations
  }
  if (!any(keep))
  {
    stop("No phylogeny left after the selection of the maps and the simulations.",
         if (!is.null(simulations)) paste0("\nSimulations available: ", paste(sort(unique(listed$simulation)), collapse = ", "), ".") else "")
  }

  ### ---- Build the multiPhylo output ---- ###

  phylos <- listed$phylos[keep]
  for (j in seq_along(phylos))
  {
    if (!is.na(listed$map[keep][j])) { attr(phylos[[j]], "map") <- listed$map[keep][j] }
    if (!keep_simmap) { phylos[[j]]$simmap <- NULL }
  }
  names(phylos) <- phylo_names[keep]
  table_phylos <- data.frame(name = phylo_names[keep], map = listed$map[keep], simulation = listed$simulation[keep], seed = listed$seed[keep], stringsAsFactors = FALSE)
  rownames(table_phylos) <- NULL
  seeds <- if (all(is.na(table_phylos$seed))) NULL else stats::setNames(table_phylos$seed, table_phylos$name)
  rownames(failed) <- NULL; rownames(failed_maps) <- NULL

  # A single phylogeny stays a single phylogeny
  if (inherits(phy_grafted, "phylo")) { return(phylos[[1]]) }

  structure(phylos, class = "multiPhylo", phylogenies = table_phylos, seeds = seeds, failed = failed, failed_maps = failed_maps)
}
