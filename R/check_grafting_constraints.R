##### Function to check grafting constraints #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Main function to check that grafted taxa respect their clade assignments ####

#' @title Check that grafted taxa respect their grafting constraints
#'
#' @description This function checks, a posteriori, that each taxon grafted with `graft_missing_taxa()` sits at a valid position in the grafted phylogeny:
#' within its target clade (or on its stem branch if allowed), and outside its excluded subclades.
#' It can also check that the backbone phylogeny is preserved by the grafting once the grafted taxa are pruned.
#'
#' @param phy_grafted An object of class `"phylo"`. The phylogeny after grafting. It can also be a set of grafted phylogenies (class `"multiPhylo"`, from `nb_simulations > 1`),
#'   or a list with one element per stochastic map (a `"phylo"` or a `"multiPhylo"`, from the grafting functions used with a `"multiSimmap"` to account for multiple alternative trait evolutionary histories).
#'   See 'Several grafted phylogenies' in Details section.
#' @param phy_backbone An object of class `"phylo"`. The backbone phylogeny used for grafting (i.e., before grafting).
#' @param missing_taxa_list A list with the information on the grafted taxa. Typically the same list as used for grafting.
#'   Items: `$Taxon_name`, `$Target_clade`, and optionally `$include_target_clade_stem`, `$Excluded_subclades` and `$exclude_subclades_stem`.
#'   Every item provided must be named after `$Taxon_name` (same names, same order), as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa], see this function for details.
#' @param include_target_clade_stem Logical. Same as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]: whether grafting on the stem branch of target clades is accepted.
#'   Taxa with a value in `missing_taxa_list$include_target_clade_stem` use it instead (`NA` or no item: this global value is used). Default is `TRUE`.
#' @param use_excluded_subclades Logical. Same as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]: whether to check the excluded subclades.
#'   If no `$Excluded_subclades` item is found in `missing_taxa_list`, a warning is issued and excluded subclades are not checked. Default is `TRUE`.
#' @param exclude_subclades_stem Logical. Same as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]: whether the stem of the excluded subclades is also excluded
#'   (a taxon sister to an excluded subclade is then a violation). Excluded subclades with a value in `missing_taxa_list$exclude_subclades_stem` use it instead
#'   (`NA` or no item: this global value is used). Default is `TRUE`.
#' @param allowed_states Character vector. Same as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]: the states in which grafting is allowed.
#'   This is the default for all taxa: taxa with a value in `missing_taxa_list$allowed_states` use this specific value instead.
#'   Only used if `phy_grafted` has a `$simmap` item (the updated simmap returned by [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] with a `simmap`).
#'   If the `$simmap` item is a list of the simmaps of several traits (grafting constrained with several traits), `allowed_states` is a list with one element per trait
#'   (in the order of the traits, or named after the traits; `NULL` = all the states of the trait), and the taxon-level values of `missing_taxa_list$allowed_states`
#'   are given per trait as in the grafting functions (see 'Several traits' in Details).
#'   The function stops if some states are not in the `simmap`. Default is `NULL` (all states allowed for the taxa without a taxon-level value).
#' @param treat_states_as_ranges Logical. Same as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]: whether the states of the `simmap` and the allowed states are treated as biogeographic ranges:
#'   an allowed range enables grafting on branch pieces in any range that contains it: `"A"` can be grafted on `"A"`, `"AB"` and `"ABC"`, whereas `"AB"` can be grafted on `"AB"` and `"ABC"` but not `"A"`.
#'   Only used if `phy_grafted` has a `$simmap` item, otherwise the function stops if `TRUE`. Not supported yet if the `$simmap` item covers several traits (the function stops).
#'   States / ranges are coded as in BioGeoBEARS (see the Details section "Biogeographic ranges as states").
#'   The function stops if the states of the `simmap` or the allowed states are not correctly coded, or if an area of the allowed states is not in the `simmap`.
#'   Default is `FALSE` (states are matched exactly).
#' @param check_backbone Logical. Whether to check that the topology and branch lengths of the backbone are preserved once the grafted taxa are pruned. Default is `TRUE`.
#' @param tolerance Numeric. Tolerance for the maximum deviation in backbone branch lengths, relative to the depth (root age) of the backbone.
#'   Default is `1e-6`.
#' @param stop_on_failure Logical. Whether to stop with an error if any constraint is violated. Default is `FALSE`.
#' @param verbose Logical. Whether to print a summary of the checks. Default is `TRUE`.
#'
#' @return The function returns a data.frame with one row per taxon in `missing_taxa_list` (in the original order) for each phylogeny checked, and the following columns:
#'
#' * `$map`: Name of the stochastic map used to constrain grafting based on species attributes. `NA` if none was used.
#' * `$simulation`: Integer. Index of the simulation (the index in the name `"sim_k"` of the phylogeny; 1 for a single phylogeny).
#' * `$seed`: Integer. Seed of the simulation (attribute `"seed"` of the phylogeny, `NA` if unknown).
#' * `$Taxon_name`: Name of the taxon.
#' * `$in_tree`: Logical. Whether the taxon is found in the grafted phylogeny.
#' * `$include_target_clade_stem`: Logical. Whether the stem of the target clade was accepted for this taxon (`missing_taxa_list$include_target_clade_stem` of the taxon if provided, otherwise the global `include_target_clade_stem` argument).
#' * `$target_position`: Character string. Position of the taxon relative to its target clade: `"crown"` (within the clade), `"stem"` (sister to the clade, or to other missing taxa grafted on the stem),
#'      `"tip_branch"` (on the branch of the tip, for single-tip target clades), or `"outside"` (neither within the clade nor on its stem = an error).
#' * `$target_ok`: Logical. Whether the position respects the target clade (`"crown"` or, if the stem is accepted for the taxon, `"stem"`; `"tip_branch"` for single-tip clades).
#' * `$excluded_ok`: Logical. Whether the taxon is outside all its excluded subclades (and not on their stem, for the subclades with an excluded stem).
#' * `$n_excluded_violated`: Integer. Number of excluded subclades violated.
#' * `$excluded_violated`: Character string. Names of the excluded subclades violated (or their indices if the subclades are not named), separated by commas (`NA` if none).
#' * `$states_ok`: Logical. Only checked if `phy_grafted` has a `$simmap` item: whether the taxon was grafted in one of its allowed states
#'      (`missing_taxa_list$allowed_states` of the taxon if provided, otherwise the global `allowed_states` argument; all the states if neither is provided).
#'      With `treat_states_as_ranges = TRUE`, the state of the terminal branch of the taxon (the range it inherited) must be one of its allowed ranges (any range if the taxon has no allowed states).
#'      `NA` if `phy_grafted` has no `$simmap`.
#'      With the maps of several traits, `TRUE` if the taxon was grafted in an allowed state for all the traits.
#' * `$states_ok_<trait>`: Logical. Only with the maps of several traits (one column per trait, named after the traits): whether the taxon was grafted in one of its allowed states for this trait.
#'      Used to find which trait is violated when `$states_ok` is `FALSE`.
#' * `$ranges_ok`: Logical. Only checked with `treat_states_as_ranges = TRUE` and a `$simmap` item: whether the range inherited by the taxon is contained in the state of the lineage at the graft point
#'      (see the Details section "Biogeographic ranges as states"). `NA` if not checked, or if the taxon is grafted at the root (no parent lineage).
#'      This differs from `$states_ok` that checks only the validity of the terminal tip state. Here, the validity of the inherited range is assessed against the range of the backbone lineage it was grafted on.
#' * `$ok`: Logical. Whether the taxon respects its target clade, its excluded subclades, and its allowed states and ranges (if checked).
#'
#' With a list of maps as `phy_grafted`, the function returns a list of such data.frames (one per map, named after the maps, with the attribute `"all_ok"` for all the maps),
#' since the maps have their own constraints: in each data.frame, the rows of all the simulations of the map are merged.
#'
#' The data.frame has the following attributes (for several phylogenies: `"topology_preserved"` and `"backbone_preserved"` are `TRUE` if they are `TRUE` for all the phylogenies,
#' the deviations are the maximum over the phylogenies, and `"all_ok"` is `TRUE` if it is `TRUE` for all of them):
#'
#' * `"topology_preserved"`: Logical. Whether the backbone topology is preserved once the grafted taxa are pruned (`NA` if `check_backbone = FALSE`).
#' * `"max_branch_length_deviation"`: Numeric. Maximum absolute deviation in backbone branch lengths (`NA` if the topology is not preserved).
#' * `"max_branch_length_deviation_relative"`: Numeric. Same, relative to the depth of the backbone.
#' * `"backbone_preserved"`: Logical. Whether both topology and branch lengths are preserved, within `tolerance` (`NA` if `check_backbone = FALSE`).
#' * `"all_ok"`: Logical. Whether all taxa respect their constraints and the backbone is preserved (if checked).
#' * `"trait_names"`: Character vector. Names of the traits checked (only with the maps of several traits, otherwise `NULL`).
#' * `"phylogenies"`: data.frame with one row per phylogeny checked (columns `map`, `simulation`, `seed`, `n_failed` (number of taxa that failed), `backbone_preserved` and `all_ok`).
#'
#' Attributes are accessed with `attr(x = output, which = "attribute_name")`.
#'
#'
#' @details
#' A taxon respects its target clade if it is found within the clade, or on its stem branch if the stem is accepted for this taxon:
#' global `include_target_clade_stem = TRUE`, or `missing_taxa_list$include_target_clade_stem` of this taxon
#' (the taxon is then sister to the target clade, or to other missing taxa sister to the target clade, with no other backbone tip in between).
#' For a single-tip target clade, the taxon must be on the branch of the tip, with no other backbone tip in between.
#'
#' The 'stem' position is not restricted to a taxon being directly sister to the target clade: taxa grafted on the stem can receive further grafts,
#' so several missing taxa can form together a group sister to the target clade (i.e, a missing-only subtree grafted on the the stem).
#' Only backbone tips are considered when checking the position of a taxon (missing taxa are ignored), so such configurations are valid.
#'
#' If `phy_grafted` has a `$simmap` item (the updated `simmap` returned by the grafting functions when a `simmap` is provided), the states are also checked:
#' a grafted taxon is in the state of the point where it was grafted, along its whole terminal branch, so the states of its terminal branch in the `simmap` must be allowed.
#' The `simmap` item is used as it is: it is not compared to the backbone.
#'
#' A taxon respects its excluded subclades if it is not nested within any of them.
#' If the stem of an excluded subclade is excluded (`exclude_subclades_stem = TRUE`, the default), the taxon must also not be on its stem, i.e., sister to the subclade
#' (or to other missing taxa sister to the subclade), with no other backbone tip in between. For a single-tip excluded subclade, the taxon must not be on the branch of the tip.
#' As in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa], if the stem of the subclade is not excluded, the taxon can be sister to the subclade, not within it.
#'
#' Branches of the backbone are matched between the pruned grafted phylogeny and the backbone by the set of tips they subtend.
#' Branch lengths are considered preserved if the maximum deviation, relative to the depth of the backbone, is lower than `tolerance`.
#'
#' ### Several traits
#'
#' If the grafting used several traits as constraints, the grafting functions return a `$simmap` item that is a list of simmaps, one per trait, named after the traits,
#' (or `"trait_1"`, `"trait_2"`, ... if the input was not named). Then this function checks the validity of states trait by trait:
#' the states of the terminal branch of the taxon in the simmap of each trait must be allowed for this trait.
#' The allowed states of each taxon and trait are resolved as in the grafting functions: the value of `missing_taxa_list$allowed_states` for this trait if provided,
#' otherwise the global `allowed_states` for this trait, otherwise all the states of the trait (see [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#' `$states_ok_<trait>` gives the result for each trait, and `$states_ok` is `TRUE` only if all traits are allowed.
#' With `verbose = TRUE`, the number of violations is also given for each trait.
#' `treat_states_as_ranges = TRUE` is not supported yet with several traits (the function stops).
#'
#' ### Several grafted phylogenies
#'
#' The phylogenies are checked one by one against the same `phy_backbone` and `missing_taxa_list`. When `allowed_states` or `treat_states_as_ranges` are used,
#' each phylogeny is checked with its own `$simmap` item (the maps of a `"multiSimmap"` can differ). The warnings about the inputs are sent once.
#' With several phylogenies, `verbose = TRUE` prints one summary for all of them (number of phylogenies respecting all the constraints, number of violations of each type,
#' and the first phylogenies with violations), and `stop_on_failure = TRUE` stops if any phylogeny violates a constraint.
#'
#' ### Biogeographic ranges as states
#'
#' With `treat_states_as_ranges = TRUE`, the states are biogeographic ranges (one capital letter per area, e.g., `"A"`, `"AB"`, `"BC"`), and a range is a set of areas.
#' A grafted lineage inherits one of its allowed ranges, which can be a subset of the state of the lineage it is grafted on.
#' Example: the range `"A"` can be inherited from a lineage in the state `"AB"`. This then describes a subset sympatric speciation event: `(AB) -> (AB), (A)`.
#'
#' Two checks are then performed:
#' * `$states_ok`: the state of the terminal branch of the taxon is exactly one of its allowed ranges (the range inherited).
#' * `$ranges_ok`: the range inherited is contained in the state of the lineage at the graft point (e.g., `"A"` is contained in `"AB"`, but `"BC"` is not).
#'   The check is skipped (`NA`) for taxa grafted at the root.
#'
#' @importFrom ape all.equal.phylo getMRCA keep.tip node.depth.edgelength prop.part
#'
#' @importFrom stats setNames
#'
#' @export
#'
#' @seealso Grafting functions:
#'   * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade with constraints (target clades and excluded subclades),
#'   * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as clades sister to a clade of the backbone.
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' # List the missing taxa to graft
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   # Items are named after Taxon_name (same names, same order)
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")),
#'   Excluded_subclades = list(X1 = list(c("t10", "t19")), X2 = NULL, X3 = NULL, X4 = NULL))
#'
#' # --- Perform grafting --- #
#'
#' # Graft the missing taxa
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   verbose = FALSE, seed = 1234)
#'
#' # --- Check one grafted phylogenies --- #
#'
#' # Check that the grafted taxa respect their constraints, and that the backbone is preserved
#' check_results <- check_grafting_constraints(phy_grafted = grafted_phy, phy_backbone = phy,
#'                                             missing_taxa_list = missing_taxa_list)
#' check_results
#' attr(check_results, "all_ok")
#'
#' # --- Check several grafted phylogenies --- #
#'
#' # Build 3 grafted phylognies (i.e., simulations) with the same backbone and constraints
#' grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                    nb_simulations = 3, verbose = FALSE, seed = 1234)
#'
#'# Check all resulting grafted phylogenies at once
#' check_results_multi <- check_grafting_constraints(phy_grafted = grafted_phys, phy_backbone = phy,
#'                                                   missing_taxa_list = missing_taxa_list)
#'
#' # With several simulations: the data.frame has the columns 'map', 'simulation' and 'seed',
#' # and the summary is printed once for all the phylogenies
#' head(check_results_multi[, c("map", "simulation", "seed", "Taxon_name", "ok")])
#'
#' # An additional attribute provide a summary per phylogeny:
#' # number of taxa that failed, backbone preserved, all constraints respected
#' attr(check_results_multi, "phylogenies")
#'
#'


check_grafting_constraints <- function (phy_grafted, # Phylogeny after grafting
                                        phy_backbone, # Backbone phylogeny used for the grafting (before grafting)
                                        missing_taxa_list, # Same list as used for grafting. Items: Taxon_name, Target_clade, include_target_clade_stem, Excluded_subclades, exclude_subclades_stem
                                        include_target_clade_stem = TRUE, # Same as in graft_missing_taxa(): accept grafting on the stem branch of the target clade?
                                        use_excluded_subclades = TRUE, # Same as in graft_missing_taxa(): check the excluded subclades?
                                        exclude_subclades_stem = TRUE, # Same as in graft_missing_taxa(): is the stem of the excluded subclades also excluded?
                                        allowed_states = NULL, # Same as in graft_missing_taxa(): states where grafting is allowed. Only checked if phy_grafted has a $simmap item
                                        treat_states_as_ranges = FALSE, # Same as in graft_missing_taxa(): states are biogeographic ranges (allowed ranges inherited, contained in the state of the graft point)
                                        check_backbone = TRUE, # Check that the backbone topology and branch lengths are preserved after pruning the grafted taxa
                                        tolerance = 1e-6, # Tolerance for the max deviation in backbone branch lengths, relative to the depth (root age) of the backbone
                                        stop_on_failure = FALSE, # Stop with an error if any constraint is violated
                                        verbose = TRUE) # Print a summary
{
  ## List the phylogenies to check: a phylogeny, a multiPhylo (simulations of a map), or a list of maps
  listed <- list_grafted_phylogenies(phy_grafted, arg = "phy_grafted")
  n_phylos <- length(listed$phylos)

  ## A single phylogeny: the summary and the error are the ones of the check of one phylogeny
  one <- function (phy, verbose, stop_on_failure, who = NULL)
  {
    with_context <- function (expr) {
      if (is.null(who)) return(expr)
      tryCatch(expr, error = function (e) stop(who, ": ", conditionMessage(e), call. = FALSE))
    }
    with_context(check_grafting_constraints_one(
      phy_grafted = phy, phy_backbone = phy_backbone, missing_taxa_list = missing_taxa_list, include_target_clade_stem = include_target_clade_stem,
      use_excluded_subclades = use_excluded_subclades, exclude_subclades_stem = exclude_subclades_stem, allowed_states = allowed_states,
      treat_states_as_ranges = treat_states_as_ranges, check_backbone = check_backbone, tolerance = tolerance,
      stop_on_failure = stop_on_failure, verbose = verbose))
  }
  add_columns <- function (res, k) { cbind(map = listed$map[k], simulation = listed$simulation[k], seed = listed$seed[k], res, stringsAsFactors = FALSE) }
  copy_attributes <- function (res, from) { for (a in c("topology_preserved", "max_branch_length_deviation", "max_branch_length_deviation_relative", "backbone_preserved", "all_ok", "trait_names")) { attr(res, a) <- attr(from, a) }; res }

  if (n_phylos == 1 && !listed$is_list)
  {
    res <- one(listed$phylos[[1]], verbose = verbose, stop_on_failure = stop_on_failure)
    out <- copy_attributes(add_columns(res, 1), res)
    attr(out, "phylogenies") <- data.frame(map = listed$map[1], simulation = listed$simulation[1], seed = listed$seed[1],
                                           n_failed = sum(!res$ok | !res$in_tree), backbone_preserved = attr(res, "backbone_preserved"), all_ok = attr(res, "all_ok"), stringsAsFactors = FALSE)
    return(out)
  }

  ## Several phylogenies: checked one by one, the warnings about the inputs are sent once
  results <- with_unique_warnings(lapply(seq_len(n_phylos), function (k)
  {
    one(listed$phylos[[k]], verbose = FALSE, stop_on_failure = FALSE, who = paste0(if (is.na(listed$map[k])) "Phylogeny" else paste0("Map '", listed$map[k], "'"), ", simulation ", listed$simulation[k]))
  }))

  # Summary table of the phylogenies
  table_phylos <- data.frame(map = listed$map, simulation = listed$simulation, seed = listed$seed,
                             n_failed = vapply(results, function (r) sum(!r$ok | !r$in_tree), integer(1)),
                             backbone_preserved = vapply(results, function (r) attr(r, "backbone_preserved"), logical(1)),
                             all_ok = vapply(results, function (r) attr(r, "all_ok"), logical(1)), stringsAsFactors = FALSE)

  # One data.frame per map, with the rows of all its simulations
  with_columns <- lapply(seq_len(n_phylos), function (k) add_columns(results[[k]], k))
  aggregate_attributes <- function (ks, res)
  {
    parts <- results[ks]
    all_of <- function (a) { v <- vapply(parts, function (r) attr(r, a), logical(1)); if (all(is.na(v))) NA else all(v, na.rm = TRUE) && !anyNA(v) }
    max_of <- function (a) { v <- vapply(parts, function (r) attr(r, a), numeric(1)); if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) }
    attr(res, "topology_preserved") <- if (check_backbone) all(vapply(parts, function (r) isTRUE(attr(r, "topology_preserved")), logical(1))) else NA
    attr(res, "max_branch_length_deviation") <- max_of("max_branch_length_deviation")
    attr(res, "max_branch_length_deviation_relative") <- max_of("max_branch_length_deviation_relative")
    attr(res, "backbone_preserved") <- if (check_backbone) all(vapply(parts, function (r) isTRUE(attr(r, "backbone_preserved")), logical(1))) else NA
    attr(res, "all_ok") <- all(vapply(parts, function (r) isTRUE(attr(r, "all_ok")), logical(1)))
    attr(res, "trait_names") <- attr(parts[[1]], "trait_names")
    attr(res, "phylogenies") <- table_phylos[ks, , drop = FALSE]
    rownames(attr(res, "phylogenies")) <- NULL
    res
  }
  merge_map <- function (ks) { res <- do.call(rbind, with_columns[ks]); rownames(res) <- NULL; aggregate_attributes(ks, res) }
  if (listed$is_list)
  {
    maps <- unique(listed$map)
    output <- stats::setNames(lapply(maps, function (m) merge_map(which(listed$map == m))), maps)
    attr(output, "all_ok") <- all(table_phylos$all_ok)
    attr(output, "phylogenies") <- table_phylos
  } else {
    output <- merge_map(seq_len(n_phylos))
  }
  all_ok <- all(table_phylos$all_ok)

  # Aggregated summary
  if (verbose)
  {
    all_rows <- do.call(rbind, results)
    n_taxa <- nrow(results[[1]])
    check_states <- !all(is.na(all_rows$states_ok))
    trait_names <- attr(results[[1]], "trait_names")
    check_ranges <- !all(is.na(all_rows$ranges_ok))
    n_maps <- length(unique(stats::na.omit(listed$map)))
    deviation <- suppressWarnings(max(vapply(results, function (r) attr(r, "max_branch_length_deviation_relative"), numeric(1)), na.rm = TRUE))
    backbone_msg <- if (!check_backbone) { "not checked" } else {
      paste0("preserved in ", sum(table_phylos$backbone_preserved %in% TRUE), "/", n_phylos, " phylogenies, branch lengths shifted by max ",
             format(deviation, digits = 3), " of backbone depth (tolerance = ", format(tolerance), ")")
    }
    violating <- which(!table_phylos$all_ok)
    n_show <- min(length(violating), 10)
    cat(paste0(n_phylos, " phylogenies checked (", if (n_maps == 0) "no map" else paste0(n_maps, " map(s)"), "), ", n_taxa, " taxa each\n",
               "  - Phylogenies respecting all the constraints: ", sum(table_phylos$all_ok), "/", n_phylos, "\n",
               "  Violations, over all the phylogenies (number of taxa x phylogenies):\n",
               "  - Missing from the grafted tree: ", sum(!all_rows$in_tree), "\n",
               "  - Target clade violated:          ", sum(all_rows$target_ok %in% FALSE), "\n",
               "  - Excluded subclade violated:     ", sum(all_rows$excluded_ok %in% FALSE), "\n",
               if (check_states) paste0("  - Allowed states violated:        ", sum(all_rows$states_ok %in% FALSE), trait_breakdown(all_rows, trait_names), "\n") else "",
               if (check_ranges) paste0("  - Range not in graft-point state: ", sum(all_rows$ranges_ok %in% FALSE), "\n") else "",
               "  - Backbone: ", backbone_msg, "\n",
               if (length(violating) > 0) paste0("  Phylogenies with violations:\n",
                                                 paste0("    - ", ifelse(is.na(table_phylos$map[violating[seq_len(n_show)]]), "no map", paste0("map ", table_phylos$map[violating[seq_len(n_show)]])), ", simulation ", table_phylos$simulation[violating[seq_len(n_show)]],
                                                        if (!all(is.na(table_phylos$seed))) paste0(" (seed ", table_phylos$seed[violating[seq_len(n_show)]], ")") else "",
                                                        ": ", table_phylos$n_failed[violating[seq_len(n_show)]], " taxa failed",
                                                        ifelse(table_phylos$backbone_preserved[violating[seq_len(n_show)]] %in% FALSE, ", backbone NOT preserved", ""), collapse = "\n"),
                                                 if (length(violating) > n_show) paste0("\n    ... and ", length(violating) - n_show, " more") else "", "\n") else "",
               "  => ", ifelse(all_ok, "All constraints respected\n", "Some constraints NOT respected (see returned data.frame(s) and attributes)\n")))
  }

  if (stop_on_failure && !all_ok)
  {
    violating <- which(!table_phylos$all_ok)
    n_show <- min(length(violating), 10)
    stop("Some grafting constraints are not respected in ", length(violating), " of ", n_phylos, " phylogenies:\n",
         paste0("  - ", ifelse(is.na(table_phylos$map[violating[seq_len(n_show)]]), "no map", paste0("map ", table_phylos$map[violating[seq_len(n_show)]])), ", simulation ", table_phylos$simulation[violating[seq_len(n_show)]], ": ",
                table_phylos$n_failed[violating[seq_len(n_show)]], " taxa failed; backbone preserved = ", table_phylos$backbone_preserved[violating[seq_len(n_show)]], collapse = "\n"),
         if (length(violating) > n_show) paste0("\n  ... and ", length(violating) - n_show, " more") else "")
  }

  output
}


## Per-trait breakdown of the violations of the allowed states for the summary (empty string for a single trait): " (color: 1, elevation: 0)"
trait_breakdown <- function (rows, trait_names)
{
  if (is.null(trait_names)) { return("") }
  paste0(" (", paste0(trait_names, ": ", vapply(trait_names, function (t) as.character(sum(rows[[paste0("states_ok_", t)]] %in% FALSE)), character(1)), collapse = ", "), ")")
}


## Check the grafting constraints for a single phylogeny (see check_grafting_constraints())
check_grafting_constraints_one <- function (phy_grafted, # Phylogeny after grafting
                                        phy_backbone, # Backbone phylogeny used for the grafting (before grafting)
                                        missing_taxa_list, # Same list as used for grafting. Items: Taxon_name, Target_clade, include_target_clade_stem, Excluded_subclades, exclude_subclades_stem
                                        include_target_clade_stem = TRUE, # Same as in graft_missing_taxa(): accept grafting on the stem branch of the target clade?
                                        use_excluded_subclades = TRUE, # Same as in graft_missing_taxa(): check the excluded subclades?
                                        exclude_subclades_stem = TRUE, # Same as in graft_missing_taxa(): is the stem of the excluded subclades also excluded?
                                        allowed_states = NULL, # Same as in graft_missing_taxa(): states where grafting is allowed. Only checked if phy_grafted has a $simmap item
                                        treat_states_as_ranges = FALSE, # Same as in graft_missing_taxa(): states are biogeographic ranges (allowed ranges inherited, contained in the state of the graft point)
                                        check_backbone = TRUE, # Check that the backbone topology and branch lengths are preserved after pruning the grafted taxa
                                        tolerance = 1e-6, # Tolerance for the max deviation in backbone branch lengths, relative to the depth (root age) of the backbone
                                        stop_on_failure = FALSE, # Stop with an error if any constraint is violated
                                        verbose = TRUE) # Print a summary
{
  ## 1.1/ Prepare inputs ####

  # States: arguments on the states require the updated simmap in the grafted phylogeny (the states are not checked otherwise)
  if (!is.logical(treat_states_as_ranges) || length(treat_states_as_ranges) != 1 || is.na(treat_states_as_ranges))
  {
    stop("'treat_states_as_ranges' must be TRUE or FALSE.")
  }
  if (is.null(phy_grafted$simmap) && (!is.null(allowed_states) || treat_states_as_ranges))
  {
    stop("'", paste(c(if (!is.null(allowed_states)) "allowed_states", if (treat_states_as_ranges) "treat_states_as_ranges"), collapse = "' and '"),
         "' can only be used if 'phy_grafted' has a '$simmap' item (the updated simmap returned by the grafting functions when a 'simmap' is provided), ",
         "otherwise the states cannot be checked. Grafting with a 'simmap' and 'return_simmap = TRUE' (the default) provides it.")
  }

  # Check the structure of missing_taxa_list (Taxon_name, and every item provided with the same length and names as Taxon_name), and extract taxa
  taxa_names <- validate_missing_taxa_list(missing_taxa_list, clade_item = "Target_clade")
  n_taxa <- length(taxa_names)

  # Stem availability of each taxon: 'missing_taxa_list$include_target_clade_stem' of the taxon if provided, otherwise the global 'include_target_clade_stem'
  include_target_clade_stem_taxa <- resolve_include_target_clade_stem(include_target_clade_stem_taxa = missing_taxa_list$include_target_clade_stem,
                                                      include_target_clade_stem = include_target_clade_stem, n_taxa = n_taxa,
                                                      taxon_names = taxa_names)

  # Check presence of excluded subclades
  excluded_list <- missing_taxa_list$Excluded_subclades
  if (use_excluded_subclades & is.null(excluded_list))
  {
    warning("No 'Excluded_subclades' item found in missing_taxa_list: excluded subclades are not checked")
    use_excluded_subclades <- FALSE
  }

  # Stem exclusion of each excluded subclade: 'missing_taxa_list$exclude_subclades_stem' of the subclade if provided, otherwise the global 'exclude_subclades_stem'
  exclude_stem_list <- resolve_exclude_subclades_stem(exclude_subclades_stem_taxa = missing_taxa_list$exclude_subclades_stem,
                                                      exclude_subclades_stem = exclude_subclades_stem,
                                                      excluded_subclades = if (use_excluded_subclades) excluded_list else NULL, n_taxa = n_taxa,
                                                      taxon_names = taxa_names)

  # States: only checked if the grafted phylogeny has the updated simmap. The states of each taxon are resolved as in the grafting functions
  simmap_grafted <- phy_grafted$simmap
  check_states <- !is.null(simmap_grafted)
  check_ranges <- check_states && treat_states_as_ranges # Check the ranges inherited
  allowed_states_taxa <- NULL
  trait_names <- NULL # Names of the traits if the simmap is a list of the simmaps of several traits
  allowed_sets_taxa <- NULL # Allowed states of each trait for each taxon (several traits)
  if (check_states)
  {
    # Maps of several traits: a list of simmaps (one per trait). The states are checked trait by trait
    if (is.list(simmap_grafted) && !inherits(simmap_grafted, "simmap") && length(simmap_grafted) > 0 && all(vapply(simmap_grafted, inherits, logical(1), what = "simmap")))
    {
      if (treat_states_as_ranges)
      {
        stop("'treat_states_as_ranges = TRUE' is not supported with the maps of several traits ('phy_grafted$simmap' is a list of simmaps): the ranges cannot be checked with several traits yet.")
      }
      trait_names <- names(simmap_grafted)
      if (is.null(trait_names) || anyNA(trait_names) || any(trait_names == "")) { trait_names <- paste0("trait_", seq_along(simmap_grafted)) }
      if (anyDuplicated(trait_names)) { stop("The names of the traits in 'phy_grafted$simmap' must be unique.") }
      names(simmap_grafted) <- trait_names
      for (t in trait_names)
      {
        if (is.null(simmap_grafted[[t]]$maps) || length(simmap_grafted[[t]]$maps) != nrow(phy_grafted$edge))
        {
          stop("The simmap of trait '", t, "' in 'phy_grafted$simmap' must have one map per branch of the grafted phylogeny.")
        }
      }
      trait_states <- lapply(simmap_grafted, function (sm)
      {
        st <- colnames(sm$mapped.edge)
        if (is.null(st)) sort(unique(unlist(lapply(sm$maps, names)))) else st
      })
      allowed_sets_taxa <- resolve_trait_state_sets(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                    traits = list(names = trait_names, states = trait_states), taxon_names = taxa_names)
    } else {
    if (!inherits(simmap_grafted, "simmap") || is.null(simmap_grafted$maps) || length(simmap_grafted$maps) != nrow(phy_grafted$edge))
    {
      stop("'phy_grafted$simmap' must be a stochastic character map (class 'simmap') of the grafted phylogeny, with one map per branch.")
    }
    }
    if (is.null(trait_names))
    {
    states_simmap <- colnames(simmap_grafted$mapped.edge)
    if (is.null(states_simmap)) { states_simmap <- sort(unique(unlist(lapply(simmap_grafted$maps, names)))) }
    if (check_ranges)
    {
      # Allowed ranges of each taxon: the range inherited must be one of them. NULL = no allowed ranges, so any range is accepted
      # (the ranges are not expanded into the states containing them here, as the states of the grafted simmap are the ranges inherited)
      validate_range_coding(states_simmap, who = "the simmap of 'phy_grafted'")
      allowed_states_taxa <- resolve_allowed_values(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states, taxon_names = taxa_names,
                                                    check_values = function (values, who)
                                                    {
                                                      validate_range_coding(values, who = who)
                                                      check_allowed_areas(values, states_simmap, who = who)
                                                    })
    } else {
      allowed_states_taxa <- resolve_allowed_states(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                    states = states_simmap, taxon_names = taxa_names)
    }
    }
  }

  # Check tips used to define the clades are in the backbone
  defining_tips <- unique(unlist(c(missing_taxa_list$Target_clade, if (use_excluded_subclades) excluded_list)))
  if (!all(defining_tips %in% phy_backbone$tip.label))
  {
    stop(paste0("Some tip labels used to define clades are missing from the backbone phylogeny\n",
                "Missing tips: ", paste(setdiff(defining_tips, phy_backbone$tip.label), collapse = ", ")))
  }

  # Check the backbone tips are all in the grafted phylogeny
  if (!all(phy_backbone$tip.label %in% phy_grafted$tip.label))
  {
    stop("Some backbone tips are missing from the grafted phylogeny")
  }

  # Precompute tips descending from each internal node in both trees
  n_tips_bb <- length(phy_backbone$tip.label)
  n_tips_g <- length(phy_grafted$tip.label)
  parts_bb <- ape::prop.part(phy_backbone, check.labels = FALSE)
  parts_g <- ape::prop.part(phy_grafted, check.labels = FALSE)

  # Functions to extract descending tips for a given node in both the backbone and grafted trees
  tips_of_node_bb <- function (node) { phy_backbone$tip.label[parts_bb[[node - n_tips_bb]]] }
  tips_of_node_g <- function (node) { parts_g[[node - n_tips_g]] }

  # Backbone tips found in the grafted tree
  is_backbone_tip_g <- phy_grafted$tip.label %in% phy_backbone$tip.label

  # Cache for the tips of the backbone clades
  clade_cache <- new.env()
  # Function to retrieve tips of a given node in the backbone trees, but using a cache system
  # to not have to recompute it for each taxa if it was already done for its associated clade
  get_backbone_clade_tips <- function (tips)
  {
    tips <- unique(tips)
    if (length(tips) == 1) { return(tips) }
    key <- paste(sort(tips), collapse = "|")
    if (is.null(clade_cache[[key]]))
    {
      clade_cache[[key]] <- tips_of_node_bb(ape::getMRCA(phy = phy_backbone, tip = tips))
    }
    return(clade_cache[[key]])
  }

  ## 1.2/ Check each taxon ####

  # Initiate resulting data.frame
  res <- data.frame(Taxon_name = taxa_names,
                    in_tree = taxa_names %in% phy_grafted$tip.label,
                    include_target_clade_stem = include_target_clade_stem_taxa,
                    target_position = NA, # "crown", "stem", "outside", or "tip_branch" (single-tip clade)
                    target_ok = NA,
                    excluded_ok = NA,
                    n_excluded_violated = NA,
                    excluded_violated = NA, # Indices of the violated excluded subclades
                    states_ok = NA, # Only checked with a simmap in phy_grafted
                    ranges_ok = NA, # Only checked with a simmap in phy_grafted and treat_states_as_ranges = TRUE
                    ok = FALSE,
                    stringsAsFactors = FALSE)
  if (!is.null(trait_names))
  {
    # One column per trait, after the joint 'states_ok'
    trait_cols <- paste0("states_ok_", trait_names)
    for (tc in rev(trait_cols)) { res[[tc]] <- NA }
    res <- res[, c(setdiff(names(res), c(trait_cols, "ranges_ok", "ok")), trait_cols, "ranges_ok", "ok")]
  }
  graft_point_state <- rep(NA_character_, n_taxa) # State of the lineage at the graft point (only with treat_states_as_ranges = TRUE), for the messages

  ## Loop per taxa
  for (i in seq_len(n_taxa))
  {
    # i <- 1
    if (!res$in_tree[i]) { next }

    taxon_idx <- match(taxa_names[i], phy_grafted$tip.label)
    target_clade_tips_i <- unique(unname(missing_taxa_list$Target_clade[[i]]))

    ## Target clade
    C_tips <- get_backbone_clade_tips(target_clade_tips_i) # Backbone tips of the target clade

    # Smallest clade of the grafted tree containing the target clade and the taxon
    Y_node <- ape::getMRCA(phy = phy_grafted, tip = c(target_clade_tips_i, taxa_names[i]))
    Y_tips_idx <- tips_of_node_g(Y_node)
    Y_backbone_tips <- phy_grafted$tip.label[Y_tips_idx[is_backbone_tip_g[Y_tips_idx]]]
    # Check that no foreign tip (outside of the backbone target clade) are present in the associated clade in the grafted tree
    no_foreign_backbone_tip <- all(Y_backbone_tips %in% C_tips)

    ## Deal with the different cases
    if (length(target_clade_tips_i) == 1)
    {
      ## Case with a single tip

      # Single-tip clade: taxon must be on the tip branch, with no other backbone tip in between (no crown possible)
      res$target_position[i] <- ifelse(no_foreign_backbone_tip, "tip_branch", "outside")
      res$target_ok[i] <- no_foreign_backbone_tip
    } else {
      ## Case with a clade (defined by two tips)

      # Get MRCA of the target clade in the grafted tree
      X_node <- ape::getMRCA(phy = phy_grafted, tip = target_clade_tips_i)
      # Check that the taxa is found within the target clade on the grafted phylogeny
      in_crown <- taxon_idx %in% tips_of_node_g(X_node)
      # If not, use 'no_foreign_backbone_tip' to distinguish between taxa grafted on the stem branch, and taxa outside of the target clade
      res$target_position[i] <- ifelse(in_crown, "crown", ifelse(no_foreign_backbone_tip, "stem", "outside"))
      # Target is okay if in_crown or on the stem when the stem is accepted for this taxon
      res$target_ok[i] <- in_crown | (include_target_clade_stem_taxa[i] & no_foreign_backbone_tip)
    }

    ## Excluded subclades
    if (use_excluded_subclades)
    {
      excluded_i <- excluded_list[[i]]
      stem_excluded_i <- exclude_stem_list[[i]] # Is the stem of each excluded subclade also excluded?
      violated <- c()
      for (j in seq_along(excluded_i))
      {
        E_tips <- unique(unname(excluded_i[[j]]))
        if (length(E_tips) == 0) { next }

        # Case with only a single tip as Excluded subclade: it has no descendants, only its terminal branch (its stem) can be excluded
        if (length(E_tips) == 1)
        {
          in_E_crown <- FALSE
        } else {
          # Case with a clade (pairs of tips) as Excluded subclade(s)
          E_node <- ape::getMRCA(phy = phy_grafted, tip = E_tips)
          # Check that the clade containing the target taxon on the grafted tree does not include tips from the excluded subclade(s)
          in_E_crown <- taxon_idx %in% tips_of_node_g(E_node)
        }
        if (in_E_crown) { violated <- c(violated, j); next }

        # If the stem of the excluded subclade is excluded: the taxon must not be on its stem, i.e., sister to the subclade
        # (or to other missing taxa sister to it), with no other backbone tip in between
        if (stem_excluded_i[j])
        {
          Z_node <- ape::getMRCA(phy = phy_grafted, tip = c(E_tips, taxa_names[i]))
          Z_tips_idx <- tips_of_node_g(Z_node)
          Z_backbone_tips <- phy_grafted$tip.label[Z_tips_idx[is_backbone_tip_g[Z_tips_idx]]]
          if (all(Z_backbone_tips %in% get_backbone_clade_tips(E_tips))) { violated <- c(violated, j) }
        }
      }
      # Record the number of excluded subclades violated
      res$n_excluded_violated[i] <- length(violated)
      # Record the excluded subclades violated: their names if they are named, otherwise their indices
      subclade_labels <- names(excluded_i)
      violated_labels <- if (is.null(subclade_labels)) as.character(violated) else ifelse(is.na(subclade_labels[violated]) | subclade_labels[violated] == "", as.character(violated), subclade_labels[violated])
      res$excluded_violated[i] <- ifelse(length(violated) == 0, NA, paste(violated_labels, collapse = ","))
      # Taxon is ok if no excluded subclades are violated
      res$excluded_ok[i] <- length(violated) == 0
    } else {
      res$excluded_ok[i] <- TRUE
      res$n_excluded_violated[i] <- 0L
    }

    ## States of the terminal branch of the taxon (the whole branch is in the state of the graft point)
    if (check_states)
    {
      tip_edge <- match(taxon_idx, phy_grafted$edge[, 2])
      # With ranges, a taxon without allowed states (NULL) inherits any range
      if (!is.null(trait_names))
      {
        # Each trait: the states of the terminal branch in the simmap of the trait must be allowed for this trait
        for (t in trait_names)
        {
          res[[paste0("states_ok_", t)]][i] <- all(names(simmap_grafted[[t]]$maps[[tip_edge]]) %in% allowed_sets_taxa[[i]][[t]])
        }
        res$states_ok[i] <- all(unlist(res[i, paste0("states_ok_", trait_names)]))
      } else {
        res$states_ok[i] <- is.null(allowed_states_taxa[[i]]) | all(names(simmap_grafted$maps[[tip_edge]]) %in% allowed_states_taxa[[i]])
      }

      ## Range inherited: it must be contained in the state of the lineage at the graft point
      if (check_ranges)
      {
        inherited_range <- names(simmap_grafted$maps[[tip_edge]])[1]
        parent_edge <- match(phy_grafted$edge[tip_edge, 1], phy_grafted$edge[, 2]) # Branch above the node where the taxon is attached (NA at the root)
        while (!is.na(parent_edge))
        {
          parent_map_states <- names(simmap_grafted$maps[[parent_edge]])
          graft_point_state[i] <- parent_map_states[length(parent_map_states)] # State at the tipward end of the branch
          # Go up through branches that only lead to missing taxa and are in the range inherited (the clade or lineage the taxon is grafted on)
          child_node <- phy_grafted$edge[parent_edge, 2]
          leads_to_backbone <- child_node <= n_tips_g || any(is_backbone_tip_g[tips_of_node_g(child_node)])
          if (leads_to_backbone || graft_point_state[i] != inherited_range) { break }
          parent_edge <- match(phy_grafted$edge[parent_edge, 1], phy_grafted$edge[, 2])
          if (is.na(parent_edge)) { graft_point_state[i] <- NA_character_ }
        }
        if (!is.na(graft_point_state[i])) { res$ranges_ok[i] <- range_in_states(inherited_range, graft_point_state[i]) }
      }
    }

    # Overall, a taxon is okay if the rules for target clade and excluded subclades are respected, and the states are allowed (if checked)
    res$ok[i] <- res$target_ok[i] & res$excluded_ok[i] & !isFALSE(res$states_ok[i]) & !isFALSE(res$ranges_ok[i])
  }

  ## 1.3/ Check the backbone is preserved ####

  # Two separate results:
  #   - topology_preserved: same topology (and tip labels) once the grafted taxa are pruned
  #   - max_branch_length_deviation: maximum absolute difference in branch length between matched backbone branches,
  #     also given relative to the depth (root age) of the backbone. NA if the topology is not preserved.
  # Branch lengths are considered preserved if the relative deviation <= tolerance.

  # Initiate variables
  topology_preserved <- NA
  max_dev_abs <- NA
  max_dev_rel <- NA
  lengths_preserved <- NA

  # Check backbone only if requested
  if (check_backbone)
  {
    # Prune the grafted tree from the grafted tips
    pruned <- ape::keep.tip(phy = phy_grafted, tip = phy_backbone$tip.label)
    # Compare the topology of the pruned tree to the initial backbone tree
    topology_preserved <- isTRUE(ape::all.equal.phylo(target = phy_backbone, current = pruned,
                                                      use.edge.length = FALSE, use.tip.label = TRUE))

    # if the topology match, check branch lengths
    if (topology_preserved & !is.null(phy_backbone$edge.length) & !is.null(pruned$edge.length))
    {
      # Function to identify each branch by the set of backbone tips it subtends, to match branches between the two trees
      branch_keys <- function (tree)
      {
        n <- length(tree$tip.label)
        parts <- ape::prop.part(tree, check.labels = FALSE)
        node_keys <- character(n + tree$Nnode)
        node_keys[1:n] <- as.character(match(tree$tip.label, phy_backbone$tip.label))
        for (k in seq_along(parts))
        {
          node_keys[n + k] <- paste(sort(match(tree$tip.label[parts[[k]]], phy_backbone$tip.label)), collapse = ",")
        }
        return(node_keys[tree$edge[, 2]])
      }
      # Extract matching branch lengths
      lengths_bb <- stats::setNames(phy_backbone$edge.length, branch_keys(phy_backbone))
      lengths_pruned <- stats::setNames(pruned$edge.length, branch_keys(pruned))
      # Compute maximum deviation in length across all matching branches
      max_dev_abs <- max(abs(lengths_pruned[names(lengths_bb)] - lengths_bb))
      # Adjust deviation to the depth of the tree = relative deviation
      backbone_depth <- max(ape::node.depth.edgelength(phy_backbone))
      max_dev_rel <- max_dev_abs / backbone_depth
      # Compare maximum relative deviation to tolerance as the check for preserved branch lengths
      lengths_preserved <- max_dev_rel <= tolerance
    }
  }
  # backbone is fully preserved if both topology and branch lengths are matching (within tolerance)
  backbone_preserved <- if (check_backbone) isTRUE(topology_preserved) & isTRUE(lengths_preserved) else NA

  ## 1.4/ Summary ####

  # Summarize check results as attributes
  all_ok <- all(res$ok) & !isFALSE(backbone_preserved)
  attr(res, "topology_preserved") <- topology_preserved
  attr(res, "max_branch_length_deviation") <- max_dev_abs
  attr(res, "max_branch_length_deviation_relative") <- max_dev_rel
  attr(res, "backbone_preserved") <- backbone_preserved
  attr(res, "all_ok") <- all_ok
  attr(res, "trait_names") <- trait_names

  # Describe the taxa that failed (max 10), for the summary and the error message
  describe_failures <- function ()
  {
    failed <- which(!res$ok | !res$in_tree)
    if (length(failed) == 0) { return("") }
    n_show <- min(length(failed), 10)
    reason <- vapply(failed[seq_len(n_show)], function (i)
    {
      if (!res$in_tree[i]) { return("missing from the grafted tree") }
      paste(c(if (isFALSE(res$target_ok[i])) paste0("outside its target clade (position: ", res$target_position[i], ")"),
              if (isFALSE(res$excluded_ok[i])) paste0("within or on the stem of excluded subclade(s): ", res$excluded_violated[i]),
              if (isFALSE(res$states_ok[i]) && !is.null(trait_names))
                {
                  tip_edge_i <- match(match(taxa_names[i], phy_grafted$tip.label), phy_grafted$edge[, 2])
                  bad_traits <- trait_names[!unlist(res[i, paste0("states_ok_", trait_names)])]
                  paste0("grafted in a state that is not allowed (", paste0(bad_traits, ": ", vapply(bad_traits, function (t) paste(names(simmap_grafted[[t]]$maps[[tip_edge_i]]), collapse = ", "), character(1)),
                                                                           " (allowed: ", vapply(bad_traits, function (t) paste(allowed_sets_taxa[[i]][[t]], collapse = ", "), character(1)), ")", collapse = "; "), ")")
                } else if (isFALSE(res$states_ok[i])) paste0(if (check_ranges) "inherited a range that is not allowed (" else "grafted in a state that is not allowed (",
                                                    paste(names(simmap_grafted$maps[[match(match(taxa_names[i], phy_grafted$tip.label), phy_grafted$edge[, 2])]]), collapse = ", "),
                                                    "; allowed: ", paste(allowed_states_taxa[[i]], collapse = ", "), ")"),
              if (isFALSE(res$ranges_ok[i])) paste0("inherited a range that is not contained in the state of the graft point (",
                                                   names(simmap_grafted$maps[[match(match(taxa_names[i], phy_grafted$tip.label), phy_grafted$edge[, 2])]])[1],
                                                   " not in ", graft_point_state[i], ")")), collapse = "; ")
    }, character(1))
    paste0("\n  Failed taxa:\n", paste0("    - ", taxa_names[failed[seq_len(n_show)]], ": ", reason, collapse = "\n"),
           if (length(failed) > n_show) paste0("\n    ... and ", length(failed) - n_show, " more") else "")
  }

  # Display check results if requested
  if (verbose)
  {
    backbone_msg <- if (!check_backbone) { "not checked" } else if (!topology_preserved) {
      "topology NOT preserved"
    } else {
      paste0("topology preserved, branch lengths shifted by max ", format(max_dev_abs, digits = 3),
             " (", format(max_dev_rel, digits = 3), " of backbone depth; tolerance = ", format(tolerance), ")")
    }
    cat(paste0(n_taxa, " taxa checked\n",
               "  - Missing from the grafted tree: ", sum(!res$in_tree), "\n",
               "  - Target clade violated:          ", sum(res$target_ok %in% FALSE), "\n",
               "  - Excluded subclade violated:     ", sum(res$excluded_ok %in% FALSE), "\n",
               if (check_states) paste0("  - Allowed states violated:        ", sum(res$states_ok %in% FALSE), trait_breakdown(res, trait_names), "\n") else "",
               if (check_ranges) paste0("  - Range not in graft-point state: ", sum(res$ranges_ok %in% FALSE), "\n") else "",
               "  - Backbone: ", backbone_msg, if (!all(res$ok)) paste0(describe_failures(), "\n") else "", "\n",
               "  => ", ifelse(all_ok, "All constraints respected\n", "Some constraints NOT respected (see returned data.frame and attributes)\n")))
  }

  # Display error message
  if (stop_on_failure & !all_ok)
  {
    stop("Some grafting constraints are not respected: ", sum(!res$ok), " taxa failed; backbone preserved = ", backbone_preserved, describe_failures())
  }

  return(res)
}
