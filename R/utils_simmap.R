
##### Internal helpers dealing with Stochastic character maps (simmap) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Validation, states of the branch segments, and tracking of the states through grafting #####

# Rationale
#   A simmap maps the evolution of a discrete trait (or of ranges) along the branches of the backbone phylogeny:
#   each branch is cut into pieces, each piece being in one 'state'.
#   Grafting can be restricted to the pieces of the branches in 'allowed states'. These restrictions are expressed as branch segments
#   (see utils_segments.R), so the draw functions do not need to know about states.
#
#   The phylogeny changes after each graft (branches are split, new branches are added, node IDs are renumbered),
#   but the input simmap describes the backbone phylogeny only.
#   A grafted branch is assumed to be in the state of the point where it was grafted, along its whole length. Therefore:
#     * Any branch of the grafted tree that belongs to a lineage of the backbone (an original branch, or one of the pieces created by grafting on it)
#       takes its states from the corresponding branch of the simmap.
#     * Any branch of a 'missing-only subtree' (a grafted taxon, and the taxa grafted on its branch, etc.) is entirely in one state:
#       the state of the point of the backbone where the first taxon of the subtree was grafted.
#   To follow these lineages after each graft, the internal nodes of the tree are labelled (in `phy$node.label`) with a key,
#   because labels are kept by `phytools::bind.tip()` while node IDs are not.
#   The backbone is labelled by the grafting functions (see `label_backbone_nodes()` and `restore_node_labels()` in utils_grafting.R), and the new nodes by `track_graft()`:
#     * "bb<k>": original node k of the backbone. The branch leading to it belongs to the backbone branch k (identified by its tipward node).
#     * "gr<k>": node created by a graft on the backbone branch k.
#     * "ms<j>": node created by a graft on a branch of a missing-only subtree in the state number j (index in the states of the simmap).
#   The tips are identified by their label. The original `node.label` of the backbone (e.g., support values) is restored at the end.


#' @title Check that a simmap matches with the backbone phylogeny
#'
#' @description Check that `simmap` is a single stochastic character map (class `"simmap"`) of the backbone phylogeny `phy`:
#'   same tips, same topology and same branch lengths (the order of the tips and of the branches can differ),
#'   and the pieces of each branch sum to the branch length.
#'
#' @param simmap Object of class `"simmap"` (e.g., from `phytools::make.simmap()`).
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param tol Numeric. Relative tolerance (to the longest branch) for the comparison of branch lengths. Default is `1e-6`.
#'
#' @importFrom ape Ntip prop.part
#'
#' @noRd
#'
#' @return A list with 2 elements:
#'  * `states` Character vector. States found in the simmap (sorted).
#'  * `edge_map` Integer vector. For each branch of `phy` (row of `phy$edge`), the index of the same branch in the simmap.
#'

validate_simmap <- function (simmap, phy, tol = 1e-6)
{
  if (inherits(simmap, "multiSimmap"))
  {
    stop("'simmap' must be a single stochastic character map (class 'simmap'), not a set of maps (class 'multiSimmap').\n",
         "Provide one map, e.g., 'simmap = simmaps[[1]]', or run the function once per map.")
  }
  if (!inherits(simmap, "simmap") || !inherits(simmap, "phylo") || is.null(simmap$maps))
  {
    stop("'simmap' must be a stochastic character map of class 'simmap', as returned by phytools::make.simmap() or phytools::sim.history().")
  }

  # Same tips
  if (anyDuplicated(simmap$tip.label) > 0 || !setequal(simmap$tip.label, phy$tip.label) || length(simmap$tip.label) != length(phy$tip.label))
  {
    stop("The tips of 'simmap' differ from the tips of 'phy'. The simmap must be built on the backbone phylogeny.")
  }

  # Match the branches of both trees (same topology)
  edge_map <- match_edges_between_trees(phy, simmap)
  if (anyNA(edge_map))
  {
    stop("The topology of 'simmap' differs from the topology of 'phy'. The simmap must be built on the backbone phylogeny.")
  }

  # Same branch lengths
  tol_abs <- tol * max(phy$edge.length)
  if (any(abs(simmap$edge.length[edge_map] - phy$edge.length) > tol_abs))
  {
    stop("The branch lengths of 'simmap' differ from the branch lengths of 'phy'. The simmap must be built on the backbone phylogeny.")
  }

  # Valid maps: named numeric vectors (state = name, duration = value) summing to the branch length
  maps <- simmap$maps
  valid_maps <- vapply(maps, function (x) is.numeric(x) && length(x) > 0 && !is.null(names(x)) && !anyNA(names(x)) && all(nzchar(names(x))) && !anyNA(x) && all(x >= 0),
                       logical(1))
  if (!all(valid_maps))
  {
    stop("'simmap$maps' must contain, for each branch, a named numeric vector of non-negative durations (names = states).")
  }
  if (any(abs(vapply(maps, sum, numeric(1))[edge_map] - phy$edge.length) > tol_abs))
  {
    stop("The mapped durations of 'simmap$maps' do not sum to the branch lengths.")
  }

  list(states = sort(unique(unlist(lapply(maps, names), use.names = FALSE))), edge_map = edge_map)
}


#' @title Match the branches of two trees with the same tips
#'
#' @description For each branch of `phy`, find the branch of `other` leading to the same set of tips. Returns `NA` for the branches without match (different topology).
#'   If both trees have the same tip order and the same edge matrix, the branches are matched by position (fast path).
#'
#' @param phy,other Phylogenies of class `"phylo"` with the same tip labels (possibly in a different order).
#'
#' @importFrom ape Ntip prop.part
#'
#' @noRd
#'
#' @return Integer vector with one value per branch of `phy`: row of `other$edge` of the same branch.
#'

match_edges_between_trees <- function (phy, other)
{
  if (identical(phy$tip.label, other$tip.label) && identical(phy$edge, other$edge)) { return(seq_len(nrow(phy$edge))) }

  # Key of each node: the sorted ranks (in the alphabetical order of the labels) of the tips it subtends
  ranks <- rank(phy$tip.label, ties.method = "first")
  names(ranks) <- phy$tip.label
  node_keys <- function (tree)
  {
    n <- ape::Ntip(tree)
    parts <- ape::prop.part(tree)
    internal <- vapply(parts, function (i) paste(sort(ranks[tree$tip.label[i]]), collapse = "-"), character(1))
    c(as.character(ranks[tree$tip.label]), internal)
  }
  keys_phy <- node_keys(phy)
  keys_other <- node_keys(other)
  match(keys_phy[phy$edge[, 2]], keys_other[other$edge[, 2]])
}


#' @title Check that the allowed states exist in the simmap
#'
#' @description Stops with an informative message if some states are not in the simmap (e.g., typo), listing the states available.
#'
#' @param allowed_states Character vector. States to check.
#' @param states Character vector. States of the simmap.
#' @param who Character string. Describes where the states come from, used in the message (e.g., "'allowed_states'", or "taxon 'X1'").
#'
#' @noRd
#'
#' @return `allowed_states` (invisibly).
#'

check_allowed_states <- function (allowed_states, states, who = "'allowed_states'")
{
  if (!is.character(allowed_states) || anyNA(allowed_states) || length(allowed_states) == 0)
  {
    stop("The states in ", who, " must be a non-empty character vector without NA.")
  }
  unknown <- setdiff(allowed_states, states)
  if (length(unknown) > 0)
  {
    stop("The states in ", who, " are not in the simmap: ", paste0("'", unknown, "'", collapse = ", "), ".\n",
         "States available in the simmap: ", paste0("'", states, "'", collapse = ", "), ".")
  }
  invisible(allowed_states)
}


#' @title Resolve the allowed values of each taxon
#'
#' @description Common part of `resolve_allowed_states()` and `resolve_allowed_ranges()`: get, for each taxon, the values of `missing_taxa_list$allowed_states` if provided,
#'   otherwise the global `allowed_states` argument. A taxon-level value replaces the global one (no intersection).
#'   Same logic as `resolve_include_target_clade_stem()`.
#'   `NULL`, an empty vector or `NA` for a taxon use the global argument.
#'
#' @param allowed_states_taxa `NULL`, or list (named after `Taxon_name`) of character vectors: `missing_taxa_list$allowed_states`.
#' @param allowed_states `NULL`, or character vector. The global `allowed_states` argument.
#' @param taxon_names Character vector. `missing_taxa_list$Taxon_name`.
#' @param check_values Function with 2 arguments (`values`, `who`) that stops if the values are not valid
#'   (`who` describes where the values come from, used in the message).
#'
#' @noRd
#'
#' @return A list named after `taxon_names`, with the values (character vector, no duplicate) of each taxon. `NULL` for the taxa without values (neither a taxon-level value nor a global one).
#'

resolve_allowed_values <- function (allowed_states_taxa, allowed_states, taxon_names, check_values)
{
  n_taxa <- length(taxon_names)
  if (!is.null(allowed_states))
  {
    check_values(allowed_states, "the global 'allowed_states' argument")
    allowed_states <- unique(allowed_states)
  }
  if (!is.null(allowed_states_taxa) && !is.list(allowed_states_taxa))
  {
    stop("'missing_taxa_list$allowed_states' must be a list (named after 'Taxon_name') of character vectors of states.")
  }
  if (!is.null(allowed_states_taxa) && length(allowed_states_taxa) != n_taxa)
  {
    stop("'missing_taxa_list$allowed_states' must have the same length as 'Taxon_name' (", length(allowed_states_taxa), " vs ", n_taxa, ").")
  }

  resolved <- vector("list", n_taxa)
  names(resolved) <- taxon_names
  for (i in seq_len(n_taxa))
  {
    value <- if (is.null(allowed_states_taxa)) NULL else allowed_states_taxa[[i]]
    use_global <- length(value) == 0 || (length(value) == 1 && is.na(value))
    if (use_global)
    {
      resolved[i] <- list(allowed_states)
    } else {
      check_values(value, paste0("'missing_taxa_list$allowed_states' of taxon '", taxon_names[i], "'"))
      resolved[[i]] <- unique(value)
    }
  }
  resolved
}


#' @title Resolve the allowed states of each taxon
#'
#' @description Get the states in which each taxon can be grafted: the states of `missing_taxa_list$allowed_states` for the taxon if provided,
#'   otherwise the global `allowed_states` argument. A taxon-level value replaces the global one (no intersection).
#'   Same logic as `resolve_include_target_clade_stem()`.
#'   `NULL`, an empty vector or `NA` for a taxon use the global argument. A taxon with neither a taxon-level value nor a global one has all the states
#'   of the simmap allowed (the states are then only followed, and grafting is not restricted). The function stops if some states are not in the simmap.
#'   States are matched exactly with the states of the simmap (see `resolve_allowed_ranges()` to treat the states as biogeographic ranges).
#'
#' @param allowed_states_taxa `NULL`, or list (named after `Taxon_name`) of character vectors: `missing_taxa_list$allowed_states`.
#' @param allowed_states `NULL`, or character vector. The global `allowed_states` argument.
#' @param states Character vector. States of the simmap.
#' @param taxon_names Character vector. `missing_taxa_list$Taxon_name`.
#'
#' @noRd
#'
#' @return A list named after `taxon_names`, with the allowed states (character vector, no duplicate) of each taxon.
#'

resolve_allowed_states <- function (allowed_states_taxa = NULL, allowed_states = NULL, states, taxon_names)
{
  resolved <- resolve_allowed_values(allowed_states_taxa = allowed_states_taxa, allowed_states = allowed_states, taxon_names = taxon_names,
                                     check_values = function (values, who) { check_allowed_states(values, states, who = who) })

  # No states provided (neither taxon-level nor global): all the states are allowed
  without_states <- which(vapply(resolved, is.null, logical(1)))
  for (i in without_states) { resolved[[i]] <- states }
  resolved
}


#' @title Resolve the allowed ranges of each taxon (states treated as biogeographic ranges)
#'
#' @description Same as `resolve_allowed_states()` with `treat_states_as_ranges = TRUE`: the states of the simmap and the allowed states
#'   are biogeographic ranges (see utils_ranges.R). An allowed range R matches a state S of the simmap if R is contained in S.
#'
#'   The states of the simmap and the allowed ranges must be correctly coded (stops otherwise, naming the incorrect states and why),
#'   and the areas of the allowed ranges must be in the simmap (stops otherwise, with the areas available).
#'
#'   The allowed ranges of each taxon are then resolved into the acceptable states (the states of the simmap that contain at least one of them) and,
#'   for each acceptable state, the allowed ranges that a lineage grafted in this state can inherit.
#'   A taxon without allowed states (neither taxon-level nor global) can be grafted in all the states and inherits the whole state.
#'   The function stops if no state of the simmap contains any of the allowed ranges of a taxon, and warns if some allowed ranges are contained in no state (they are ignored).
#'
#' @details
#' Ranges are coded as in BioGeoBEARS: one capital letter per area, and the letters of a multi-area range in alphabetical order (e.g., "A", "AB", "ABC").
#' An allowed range matches a state of the simmap if it is contained in it, i.e., if all its areas are in the state.
#' Acceptable states are the states of the simmap on which the taxa can be grafted because at least one of its allowed range is contained in them.
#' Candidate ranges are the allowed ranged that can be inherited when grafted on an edge in each of the acceptable states.
#'
#' The states of the simmap below are `c("A", "B", "AB", "BC")`:
#' * Allowed range "A" matches the states "A" and "AB". The acceptable states are `c("A", "AB")`, and a lineage grafted in either of them inherits "A".
#'   On a branch in "AB", this is a subset sympatric speciation event, (AB) -> (AB), (A).
#' * Allowed range "B" matches "B", "AB" and "BC". The acceptable states are `c("B", "AB", "BC")`, and the lineage inherits "B" in each of them.
#' * Allowed range "AB" alone matches only the states that contain both areas: "AB". The lineage inherits "AB", never "A" or "B" alone.
#' * Allowed ranges `c("A", "B")` match all the states. In "AB" both are candidates, and the range inherited is drawn among them when grafting,
#'   whereas "A" is the only candidate in "A", and "B" the only one in "B" and in "BC".
#' * Allowed ranges `c("A", "AC")`: "AC" is made of areas of the simmap but is contained in no state, so it is ignored with a warning, and "A" is used alone.
#'   If it was the only allowed range, the function would stop.
#' * No taxa-level allowed states provided: all the states are acceptable and each state is its own only candidate, so the lineage inherits the whole state of the branch it is grafted on.
#'
#' With single-area states ("A", "B", "C"), an allowed range only matches the identical state: the result is the same as `resolve_allowed_states()`
#' (the acceptable states are the allowed states and each one is its own candidate).
#'
#' A taxon-level value of `missing_taxa_list$allowed_states` replaces the global `allowed_states` for this taxon (no intersection), as in `resolve_allowed_states()`.
#' Incorrect codings (e.g., "ab", "BA", "A1", "AA") stop the function, naming the incorrect states and why, and so does an allowed area absent from the simmap
#' (with the areas available). For example, with the states above, `allowed_states = "D"` stops and lists the areas available: "A", "B" and "C".
#'
#' Examples of outputs, with the states `c("A", "B", "AB", "BC")` present in the `simmap` and one taxon "X1":
#'
#' * `allowed_states = "B"` gives `$allowed_states$X1` = `c("B", "AB", "BC")` and `$range_candidates$X1` = `list(B = "B", AB = "B", BC = "B")`.
#' * `allowed_states = c("A", "B")` gives `$allowed_states$X1` = `c("A", "B", "AB", "BC")` and `$range_candidates$X1` = `list(A = "A", B = "B", AB = c("A", "B"), BC = "B")`.
#' * `allowed_states = "AB"` gives `$allowed_states$X1` = `"AB"` and `$range_candidates$X1` = `list(AB = "AB")`.
#' * No allowed states gives `$allowed_states$X1` = `c("A", "B", "AB", "BC")` and `$range_candidates$X1` = `list(A = "A", B = "B", AB = "AB", BC = "BC")`.
#'
#' A lineage can inherit a range that is not a state of the simmap. For example, with the states `c("AB", "BC")` and the allowed range "B", a lineage grafted in "AB" or "BC" inherits "B"
#' (an allowed range of the taxon), which is not a state of the simmap. This range is then a state of the lineages grafted, so taxa grafted later can be grafted on their branches
#' if it contains one of their allowed ranges (or whatever the range if they have no allowed states). These ranges are listed in `$inherited_states`, and are added to the acceptable states and to the candidates.
#' Here, with the states `c("A", "B", "AB", "BC")` present in the `simmap` and one taxon "X1":
#' * `allowed_states = "B"` gives `$allowed_states$X1` = `c("AB", "BC", "B")`, `$range_candidates$X1` = `list(AB = "B", BC = "B", B = "B")` and `$inherited_states` = `"B"`.
#'
#' An allowed range that is not a state of the simmap but is contained in some of them is accepted (it does not stop the function), with a warning explaining the mismatch:
#' the grafting is run assuming that the range can be inherited by lineages grafted in the states that contain it, which results in subset sympatric speciation events.
#' Here, the warning concerns "B", contained in "AB" and "BC". The function only stops if an area is not in the simmap, or if a range is contained in no state (or warns and ignores it if other allowed ranges can be used).
#'
#' @param allowed_states_taxa `NULL`, or list (named after `Taxon_name`) of character vectors: `missing_taxa_list$allowed_states`.
#' @param allowed_states `NULL`, or character vector. The global `allowed_states` argument.
#' @param states Character vector. States of the simmap.
#' @param taxon_names Character vector. `missing_taxa_list$Taxon_name`.
#'
#' @noRd
#'
#' @return A list with 3 items:
#'  * `allowed_states` List named after `taxon_names` of character vectors. The acceptable states of each taxon (states of the simmap, then inherited ranges, in this order),
#'    to be used as the allowed states of `resolve_allowed_states()`.
#'  * `range_candidates` List named after `taxon_names` of named lists. For each taxon, a list named after its acceptable states, giving the character vector of the allowed ranges
#'    contained in each state (the ranges a lineage grafted in this state can inherit). Without allowed states: each state is its own (only) candidate.
#'  * `inherited_states` Character vector. The ranges that can be inherited by grafted lineages and are not states of the simmap
#'    (to add to the states of the tracker with `extend_tracker_states()`), sorted. The acceptable states and candidates already include them.
#'

resolve_allowed_ranges <- function (allowed_states_taxa = NULL, allowed_states = NULL, states, taxon_names)
{
  validate_range_coding(states, who = "the simmap")
  check_values <- function (values, who)
  {
    validate_range_coding(values, who = who)
    check_allowed_areas(values, states, who = who)
  }
  resolved <- resolve_allowed_values(allowed_states_taxa = allowed_states_taxa, allowed_states = allowed_states, taxon_names = taxon_names, check_values = check_values)

  # Acceptable states and candidate ranges of each taxon, among a set of states
  expand_taxa <- function (state_set)
  {
    acceptable <- candidates <- unmatched <- vector("list", length(taxon_names))
    names(acceptable) <- names(candidates) <- names(unmatched) <- taxon_names
    for (i in seq_along(taxon_names))
    {
      if (is.null(resolved[[i]]))
      {
        # No allowed states: all the states are acceptable, and a lineage inherits the whole state of the branch
        acceptable[[i]] <- state_set
        candidates[[i]] <- as.list(state_set)
        names(candidates[[i]]) <- state_set
      } else {
        expanded <- expand_allowed_ranges(allowed_ranges = resolved[[i]], states = state_set)
        acceptable[[i]] <- expanded$acceptable_states
        candidates[[i]] <- expanded$candidates
        unmatched[[i]] <- expanded$unmatched_ranges
      }
    }
    list(acceptable = acceptable, candidates = candidates, unmatched = unmatched)
  }

  ## States of the simmap
  expanded <- expand_taxa(states)
  no_state <- which(lengths(expanded$acceptable) == 0)
  if (length(no_state) > 0)
  {
    stop("No state of the simmap contains the allowed range(s) of taxon '", taxon_names[no_state[1]], "': ", paste0("'", resolved[[no_state[1]]], "'", collapse = ", "), ".\n",
         "States available in the simmap: ", paste0("'", states, "'", collapse = ", "), ".")
  }

  # Allowed ranges contained in no state of the simmap are never used
  unmatched <- list()
  for (i in seq_along(taxon_names)) { for (range in expanded$unmatched[[i]]) { unmatched[[range]] <- c(unmatched[[range]], taxon_names[i]) } }
  if (length(unmatched) > 0)
  {
    describe <- function (range) { paste0("'", range, "' (taxa: ", paste(utils::head(unmatched[[range]], 5), collapse = ", "), if (length(unmatched[[range]]) > 5) ", ..." else "", ")") }
    warning("With 'treat_states_as_ranges = TRUE', some allowed ranges are not contained in any state of the simmap, so they are ignored: ",
            paste(vapply(names(unmatched), describe, character(1)), collapse = "; "), ".\n",
            "States available in the simmap: ", paste0("'", states, "'", collapse = ", "), ".", call. = FALSE)
  }

  # Allowed ranges that are not states of the simmap, but are contained in some of them: accepted, with a warning explaining the mismatch
  not_states <- list()
  for (i in seq_along(taxon_names))
  {
    for (range in setdiff(resolved[[i]], c(states, expanded$unmatched[[i]]))) { not_states[[range]] <- c(not_states[[range]], taxon_names[i]) }
  }
  if (length(not_states) > 0)
  {
    describe <- function (range)
    {
      paste0("'", range, "' (contained in the states ", paste0("'", states[range_in_states(range, states)], "'", collapse = ", "),
             "; taxa: ", paste(utils::head(not_states[[range]], 5), collapse = ", "), if (length(not_states[[range]]) > 5) ", ..." else "", ")")
    }
    # Example of event given in the message: the first range with the first state of the simmap containing it
    example_range <- names(not_states)[1]
    example_state <- states[range_in_states(example_range, states)][1]
    warning("With 'treat_states_as_ranges = TRUE', some allowed ranges are not states of the simmap: ", paste(vapply(names(not_states), describe, character(1)), collapse = "; "), ".\n",
            "The grafting is run assuming that these ranges can be inherited by lineages grafted in the states of the simmap that contain them, ",
            "which results in subset sympatric speciation events (e.g., (", example_state, ") -> (", example_state, "), (", example_range, ") for the range '", example_range,
            "' inherited from the state '", example_state, "'). ",
            "The ranges inherited are added to the states of the updated simmap.", call. = FALSE)
  }

  ## Ranges inherited by the grafted lineages that are not states of the simmap (e.g., "A" inherited on a branch in "AB")
  # They are states of the lineages grafted: a taxon grafted later can be grafted on their branches if they contain one of its allowed ranges
  # (or whatever the range if it has no allowed states), so the acceptable states and candidates are resolved again among the states of the simmap and these ranges
  inherited_states <- sort(setdiff(unique(unlist(expanded$candidates, use.names = FALSE)), states), method = "radix")
  if (length(inherited_states) > 0) { expanded <- expand_taxa(c(states, inherited_states)) }

  list(allowed_states = expanded$acceptable, range_candidates = expanded$candidates, inherited_states = inherited_states)
}


#' @title Warn if a simmap is provided but has no effect
#'
#' @description A simmap has an effect if it restricts grafting (some allowed states, global or taxon-level) or if the updated simmap is returned.
#'   Without any allowed states and with `return_simmap = FALSE`, the simmap is useless: this is likely a mistake of the user, so a warning is sent.
#'
#' @param allowed_states_taxa `NULL`, or list (named after `Taxon_name`) of character vectors: `missing_taxa_list$allowed_states` (before resolution).
#' @param allowed_states `NULL`, or character vector. The global `allowed_states` argument.
#' @param return_simmap Logical. The `return_simmap` argument.
#'
#' @noRd
#'
#' @return `NULL` (invisibly).
#'

warn_if_simmap_unused <- function (allowed_states_taxa = NULL, allowed_states = NULL, return_simmap = TRUE)
{
  given_in_list <- any(vapply(allowed_states_taxa, function (x) length(x) > 0 && !all(is.na(x)), logical(1)))
  if (!return_simmap && is.null(allowed_states) && !given_in_list)
  {
    warning("A 'simmap' is provided without any allowed states (neither the global 'allowed_states' argument nor 'missing_taxa_list$allowed_states') and 'return_simmap = FALSE': ",
            "the simmap does not restrict grafting and is not returned, so it has no effect.\n",
            "Provide 'allowed_states' to restrict grafting, or set 'return_simmap = TRUE' to get the states of the grafted taxa.", call. = FALSE)
  }
  invisible(NULL)
}


#' @title Table of the states along the branches of a simmap
#'
#' @description Convert the maps of a simmap in a table with one row per piece of branch (segment) in one state.
#'   The segments are described as in `get_branch_segments()`: distances from the tipward node of the branch.
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param simmap Object of class `"simmap"` matching `phy`.
#' @param edge_map Integer vector. Output of `validate_simmap()`.
#'
#' @noRd
#'
#' @return A data.frame with one row per piece of branch (zero-length pieces are dropped), ordered by node and distance, and 4 columns:
#'  * `node` Integer. Tipward node ID of the branch in `phy`.
#'  * `from` Numeric. Distance from the tipward node to the tipward end of the piece.
#'  * `to` Numeric. Distance from the tipward node to the rootward end of the piece. The tipward piece starts at 0 and the rootward piece ends at the branch length (exactly).
#'  * `state` Character. State of the piece.
#'

build_state_table <- function (phy, simmap, edge_map)
{
  maps <- simmap$maps[edge_map]
  n_pieces <- lengths(maps)
  edge_row <- rep(seq_along(maps), n_pieces)
  duration <- unlist(maps, use.names = FALSE)
  state <- unlist(lapply(maps, names), use.names = FALSE)
  branch_length <- phy$edge.length[edge_row]

  # Maps list the pieces from the rootward node to the tipward node: distance from the rootward node to the rootward end (start) and tipward end (end) of each piece
  end_from_root <- stats::ave(duration, edge_row, FUN = cumsum)
  start_from_root <- end_from_root - duration
  first_piece <- c(TRUE, edge_row[-1] != edge_row[-length(edge_row)])
  last_piece <- c(edge_row[-1] != edge_row[-length(edge_row)], TRUE)

  # Distances from the tipward node. The ends of the branch are set exactly (the sum of the durations can differ from the length by rounding errors)
  to <- branch_length - start_from_root
  from <- branch_length - end_from_root
  to[first_piece] <- branch_length[first_piece]
  from[last_piece] <- 0

  table_states <- data.frame(node = phy$edge[edge_row, 2], from = from, to = to, state = state, stringsAsFactors = FALSE)
  table_states <- table_states[table_states$to > table_states$from, , drop = FALSE]
  table_states <- table_states[order(table_states$node, table_states$from), , drop = FALSE]
  rownames(table_states) <- NULL
  table_states
}


#' @title Start tracking the states along the branches during grafting
#'
#' @description Prepare the tracking of the states of the simmap along the lineages of the backbone while taxa are grafted (see the header of this file).
#'
#' @param phy Backbone phylogeny of class `"phylo"`, ultrametric, with its internal nodes labelled by `label_backbone_nodes()`.
#' @param simmap Object of class `"simmap"` matching `phy`.
#' @param tol Numeric. Tolerance used to validate the simmap. Default is `1e-6`.
#'
#' @importFrom phytools nodeHeights
#'
#' @noRd
#'
#' @return A list used by the functions below to follow the states (the tracker): the table of the states of the backbone, its states,
#'   and the information on the backbone nodes.
#'

start_state_tracking <- function (phy, simmap, tol = 1e-6)
{
  # The internal nodes of the backbone must be labelled to be followed during grafting
  if (!identical(phy$node.label, paste0("bb", ape::Ntip(phy) + seq_len(phy$Nnode))))
  {
    stop("Internal error: the nodes of the backbone must be labelled with label_backbone_nodes() to track the states during grafting.")
  }
  checked <- validate_simmap(simmap = simmap, phy = phy, tol = tol)
  state_table <- build_state_table(phy = phy, simmap = simmap, edge_map = checked$edge_map)
  n_tips <- ape::Ntip(phy)
  n_nodes <- n_tips + phy$Nnode

  # Age of each node of the backbone (time before present)
  node_heights <- phytools::nodeHeights(phy)
  root_age <- max(node_heights[, 2])
  node_age <- numeric(n_nodes)
  node_age[phy$edge[, 2]] <- root_age - node_heights[, 2]
  node_age[n_tips + 1] <- root_age

  tracker <- list(
    states = checked$states,
    state_table = state_table,
    rows_per_node = split(seq_len(nrow(state_table)), factor(state_table$node, levels = seq_len(n_nodes))),
    tip_labels = phy$tip.label,
    n_tips = n_tips,
    node_age = node_age,
    extra_states = character(0), # States added by extend_tracker_states() (ranges that can be inherited and are not states of the simmap)
    missing_state = character(0)) # State of each grafted taxon (named by its label)

  tracker
}


#' @title Add states to the tracker
#'
#' @description With `treat_states_as_ranges = TRUE`, a lineage can inherit a range that is not a state of the simmap (e.g., "A" grafted on a branch in "AB",
#'   when the simmap has no state "A"). These ranges are added to the states of the tracker (after the states of the simmap, so their indices do not change),
#'   so they can be recorded for the grafted lineages and be the state of the branches a taxon grafted later is grafted on.
#'
#' @param tracker Output of `start_state_tracking()`.
#' @param extra_states Character vector. States to add (the states already in the tracker are ignored).
#'
#' @noRd
#'
#' @return The updated `tracker`. The added states are listed in `tracker$extra_states`.
#'

extend_tracker_states <- function (tracker, extra_states)
{
  new_states <- setdiff(extra_states, tracker$states)
  tracker$states <- c(tracker$states, new_states)
  tracker$extra_states <- c(tracker$extra_states, new_states)
  tracker
}


#' @title Identify the lineage of branches of the current tree
#'
#' @description For branches of the current tree (identified by their tipward node), find the lineage they belong to
#'   (a branch of the backbone, or a missing-only subtree in a given state), see the header of this file.
#'
#' @param tracker Output of `start_state_tracking()`, updated by `track_graft()`.
#' @param phy Current phylogeny of class `"phylo"`, with the nodes labelled by the tracking functions.
#' @param nodes Integer vector. Tipward node IDs of the branches in `phy`.
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy` in the order of `phy$edge`: `root_age - phytools::nodeHeights(phy)`.
#'
#' @noRd
#'
#' @return A data.frame with one row per node, and 3 columns:
#'  * `backbone_node` Integer. For the branches of a backbone lineage, the tipward node ID of the backbone branch. `NA` for the missing-only subtrees.
#'  * `distance` Numeric. For the branches of a backbone lineage, distance from the tipward node of the backbone branch to the tipward node of the branch in the current tree
#'     (0 for the original branches).
#'  * `state` Character. For the branches of a missing-only subtree, its state. `NA` for the branches of the backbone lineages.
#'

get_lineage_info <- function (tracker, phy, nodes, edge_ages)
{
  n_tips <- tracker$n_tips
  n_current_tips <- ape::Ntip(phy)
  is_tip <- nodes <= n_current_tips
  label <- character(length(nodes))
  label[is_tip] <- phy$tip.label[nodes[is_tip]]
  label[!is_tip] <- phy$node.label[nodes[!is_tip] - n_current_tips]

  backbone_node <- rep(NA_integer_, length(nodes))
  distance <- rep(NA_real_, length(nodes))
  state <- rep(NA_character_, length(nodes))

  # Tips: backbone tips, or grafted taxa
  tip_in_backbone <- match(label[is_tip], tracker$tip_labels)
  is_backbone_tip <- !is.na(tip_in_backbone)
  backbone_node[is_tip][is_backbone_tip] <- tip_in_backbone[is_backbone_tip]
  distance[is_tip][is_backbone_tip] <- 0
  grafted_state <- unname(tracker$missing_state[label[is_tip][!is_backbone_tip]])
  state[is_tip][!is_backbone_tip] <- grafted_state

  # Internal nodes: key in the label
  internal_label <- label[!is_tip]
  kind <- substr(internal_label, 1, 2)
  value <- suppressWarnings(as.integer(substring(internal_label, 3)))
  internal_backbone <- !is.na(kind) & kind %in% c("bb", "gr") & !is.na(value)
  internal_missing <- !is.na(kind) & kind == "ms" & !is.na(value)

  internal_backbone_node <- rep(NA_integer_, length(internal_label))
  internal_distance <- rep(NA_real_, length(internal_label))
  internal_state <- rep(NA_character_, length(internal_label))
  internal_backbone_node[internal_backbone] <- value[internal_backbone]
  internal_distance[internal_backbone & kind == "bb"] <- 0
  graft_node <- internal_backbone & kind == "gr"
  if (any(graft_node))
  {
    age_node <- edge_ages[match(nodes[!is_tip][graft_node], phy$edge[, 2]), 2]
    internal_distance[graft_node] <- age_node - tracker$node_age[value[graft_node]]
  }
  internal_state[internal_missing] <- tracker$states[value[internal_missing]]
  backbone_node[!is_tip] <- internal_backbone_node
  distance[!is_tip] <- internal_distance
  state[!is_tip] <- internal_state

  # Every branch must belong to a lineage: otherwise the tree was modified outside of the tracking functions
  untracked <- is.na(backbone_node) & is.na(state)
  if (any(untracked))
  {
    stop("Internal error: ", sum(untracked), " branch(es) of the tree are not tracked (tipward nodes: ", paste(utils::head(nodes[untracked], 5), collapse = ", "), ").\n",
         "The states along the branches can only be followed for trees modified with the tracking functions.")
  }
  data.frame(backbone_node = backbone_node, distance = distance, state = state, stringsAsFactors = FALSE)
}


#' @title Get the segments of the branches of the current tree in allowed states
#'
#' @description For branches of the current tree, get the segments (pieces of the branches) in one of the allowed states.
#'   The segments in contiguous allowed states are merged, so that if all the states are allowed each branch is a single segment covering the whole branch.
#'
#' @param tracker Output of `start_state_tracking()`, updated by `track_graft()`.
#' @param phy Current phylogeny of class `"phylo"`, with the nodes labelled by the tracking functions.
#' @param nodes Integer vector. Tipward node IDs of the branches in `phy`.
#' @param allowed_states Character vector. States in which grafting is allowed.
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy` in the order of `phy$edge`: `root_age - phytools::nodeHeights(phy)`.
#'
#' @noRd
#'
#' @return Segments in the format of `get_branch_segments()` (columns `node`, `from`, `to`; `from` and `to` are distances from the tipward node of the
#'   branch in the current tree), ordered as `nodes`. Branches without any segment in an allowed state are absent.
#'

get_allowed_segments_tracked <- function (tracker, phy, nodes, allowed_states, edge_ages)
{
  info <- get_lineage_info(tracker = tracker, phy = phy, nodes = nodes, edge_ages = edge_ages)
  branch_length <- phy$edge.length[match(nodes, phy$edge[, 2])]
  cand <- integer(0); from <- numeric(0); to <- numeric(0)

  ## Branches of the backbone lineages: pieces of the backbone branch in an allowed state, cut to the extent of the branch
  is_backbone <- !is.na(info$backbone_node)
  if (any(is_backbone))
  {
    idx <- which(is_backbone)
    rows <- tracker$rows_per_node[info$backbone_node[idx]]
    n_rows <- lengths(rows)
    cand_rep <- rep(idx, n_rows)
    r <- unlist(rows, use.names = FALSE)
    keep <- tracker$state_table$state[r] %in% allowed_states
    cand_rep <- cand_rep[keep]; r <- r[keep]
    lo <- info$distance[cand_rep]
    hi <- lo + branch_length[cand_rep]
    # A piece that reaches an end of the branch ends exactly there (0 or the branch length), without rounding errors from the subtractions,
    # so a branch entirely in allowed states is a single segment covering exactly the whole branch (as without simmap).
    # The tolerance is relative: phytools::bind.tip() rounds all the branch lengths to 10 significant digits at each graft,
    # so the lengths of the current tree differ slightly from the lengths of the simmap (relative differences of the order of 1e-10)
    eps <- 1e-8 * hi
    f <- ifelse(tracker$state_table$from[r] <= lo + eps, 0, tracker$state_table$from[r] - lo)
    t <- ifelse(tracker$state_table$to[r] >= hi - eps, branch_length[cand_rep], tracker$state_table$to[r] - lo)
    keep <- t > f
    cand <- c(cand, cand_rep[keep]); from <- c(from, f[keep]); to <- c(to, t[keep])
  }

  ## Branches of the missing-only subtrees: the whole branch is in one state
  is_missing <- !is_backbone & info$state %in% allowed_states
  if (any(is_missing))
  {
    cand <- c(cand, which(is_missing)); from <- c(from, rep(0, sum(is_missing))); to <- c(to, branch_length[is_missing])
  }

  # Order as the input nodes, then along the branch, and merge the contiguous segments
  ord <- order(cand, from)
  cand <- cand[ord]; from <- from[ord]; to <- to[ord]
  n <- length(cand)
  if (n > 1)
  {
    new_segment <- c(TRUE, cand[-1] != cand[-n] | from[-1] != to[-n])
    last_of_segment <- c(new_segment[-1], TRUE)
    cand <- cand[new_segment]; from <- from[new_segment]; to <- to[last_of_segment]
  }
  data.frame(node = nodes[cand], from = from, to = to)
}


#' @title State of a position along a branch of the current tree
#'
#' @description Get the state at a given position along a branch of the current tree. If the position is exactly at the limit of two pieces in different states,
#'   the state among the allowed ones is returned (if possible).
#'
#' @param tracker Output of `start_state_tracking()`, updated by `track_graft()`.
#' @param phy Current phylogeny of class `"phylo"`, with the nodes labelled by the tracking functions.
#' @param node Integer. Tipward node ID of the branch in `phy`.
#' @param position Numeric. Distance from the tipward node of the branch.
#' @param allowed_states Character vector. States in which grafting is allowed. Used to solve ties at the limit of two pieces. Default is `NULL` (all states).
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy`.
#'
#' @noRd
#'
#' @return Character. The state.
#'

get_state_at_position <- function (tracker, phy, node, position, allowed_states = NULL, edge_ages)
{
  info <- get_lineage_info(tracker = tracker, phy = phy, nodes = node, edge_ages = edge_ages)
  if (is.na(info$backbone_node)) { return(info$state) }

  d <- info$distance + position
  r <- tracker$rows_per_node[[info$backbone_node]]
  tab <- tracker$state_table[r, , drop = FALSE]
  tol <- 1e-9 * max(1, abs(d))
  hit <- tab$from - tol <= d & d <= tab$to + tol
  candidates <- tab$state[hit]
  if (length(candidates) == 0) { stop("Internal error: no state found at position ", signif(d, 6), " of the backbone branch ", info$backbone_node) }
  if (!is.null(allowed_states) && any(candidates %in% allowed_states)) { candidates <- candidates[candidates %in% allowed_states] }
  candidates[1]
}


#' @title Choose the state inherited by a grafted lineage
#'
#' @description Without candidate ranges, a grafted lineage takes the state of the point where it was grafted. With candidate ranges
#'   (`treat_states_as_ranges = TRUE`), it inherits one of the allowed ranges contained in this state, drawn uniformly. No random number is drawn
#'   if there is a single candidate, so single-area states give the same trees as without ranges.
#'
#' @param state Character string. State of the point of the branch where the lineage is grafted.
#' @param range_candidates `NULL`, or named list (named after the states) of character vectors: the ranges that can be inherited in each state.
#'
#' @noRd
#'
#' @return Character string. The state of the grafted lineage.
#'

choose_inherited_state <- function (state, range_candidates = NULL)
{
  if (is.null(range_candidates)) { return(state) }
  candidates <- range_candidates[[state]]
  if (is.null(candidates))
  {
    stop("Internal error: no range can be inherited when grafting on a branch in state '", state, "'.")
  }
  sample_one(candidates)
}


#' @title Update the tracking after a graft
#'
#' @description After a taxon was grafted with `phytools::bind.tip()`, record its state (the state of the point where it was grafted) and label the new node,
#'   so the new branches keep being tracked.
#'
#' @param tracker Output of `start_state_tracking()`, updated by `track_graft()`.
#' @param phy_before Phylogeny before the graft, with the nodes labelled by the tracking functions.
#' @param phy_after Phylogeny after the graft (output of `phytools::bind.tip()`).
#' @param node Integer. Tipward node ID, in `phy_before`, of the branch the taxon was grafted on.
#' @param position Numeric. Position of the graft on this branch: distance from its tipward node.
#' @param tip_label Character string. Label of the grafted taxon.
#' @param allowed_states Character vector. States in which grafting was allowed (used to solve ties at the limit of two pieces). Default is `NULL` (all states).
#' @param range_candidates `NULL`, or named list (named after the states) of character vectors: the ranges the new taxon can inherit when grafted on a branch in each state
#'   (`treat_states_as_ranges = TRUE`, see `resolve_allowed_ranges()`). One of the candidates of the state of the point of the graft is drawn (without drawing a random number if there is a single one).
#'   Default is `NULL`: the taxon takes the state of the point where it was grafted.
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy_before`.
#'
#' @noRd
#'
#' @return A list with 2 elements: the updated `phy` (`phy_after` with its new node labelled) and the updated `tracker`.
#'

track_graft <- function (tracker, phy_before, phy_after, node, position, tip_label, allowed_states = NULL, range_candidates = NULL, edge_ages)
{
  info <- get_lineage_info(tracker = tracker, phy = phy_before, nodes = node, edge_ages = edge_ages)
  state <- get_state_at_position(tracker = tracker, phy = phy_before, node = node, position = position, allowed_states = allowed_states, edge_ages = edge_ages)

  # State of the new taxon: the state of the point of the graft, or one of the ranges it can inherit there
  inherited_state <- choose_inherited_state(state = state, range_candidates = range_candidates)

  # Label of the new node: backbone lineage or missing-only subtree (the node is on the branch where the taxon was grafted, so it has the state of this branch, not the inherited one)
  new_label <- if (!is.na(info$backbone_node)) paste0("gr", info$backbone_node) else paste0("ms", match(state, tracker$states))

  # The new node is the parent of the new tip. It is unlabelled unless the tip was attached directly to an existing node
  # (phytools::bind.tip() gives the character string "NA" to the new nodes of a tree with node labels)
  new_tip <- match(tip_label, phy_after$tip.label)
  new_node <- phy_after$edge[match(new_tip, phy_after$edge[, 2]), 1]
  node_index <- new_node - ape::Ntip(phy_after)
  current_label <- phy_after$node.label[node_index]
  if (is.na(current_label) || !nzchar(current_label) || current_label == "NA")
  {
    phy_after$node.label[node_index] <- new_label
  }

  tracker$missing_state[tip_label] <- inherited_state
  list(phy = phy_after, tracker = tracker)
}


#' @title Update the tracking after the graft of a subtree
#'
#' @description Same as `track_graft()`, for a clade of several taxa grafted with `ape::bind.tree()` (e.g., a missing clade grafted on the stem branch of its sister clade).
#'   The state of the whole subtree is the state of the point of the backbone where it was grafted.
#'   The nodes created by the graft are labelled: the node where the subtree is attached to the branch, and the nodes within the subtree.
#'
#' @param tracker Output of `start_state_tracking()`, updated by `track_graft()` and `track_graft_subtree()`.
#' @param phy_before Phylogeny before the graft, with the nodes labelled by the tracking functions.
#' @param phy_after Phylogeny after the graft (output of `ape::bind.tree()`).
#' @param node Integer. Tipward node ID, in `phy_before`, of the branch the subtree was grafted on.
#' @param position Numeric. Position of the graft on this branch: distance from its tipward node.
#' @param tip_labels Character vector. Labels of the taxa of the grafted subtree (at least 2).
#' @param allowed_states Character vector. States in which grafting was allowed (used to solve ties at the limit of two pieces). Default is `NULL` (all states).
#' @param range_candidates `NULL`, or named list of character vectors: the ranges the subtree can inherit when grafted on a branch in each state (see `track_graft()`).
#'   One range is drawn for the whole subtree. Default is `NULL`: the subtree takes the state of the point where it was grafted.
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy_before`.
#'
#' @importFrom ape getMRCA Ntip
#'
#' @noRd
#'
#' @return A list with 2 elements: the updated `phy` (`phy_after` with the new nodes labelled) and the updated `tracker`.
#'

track_graft_subtree <- function (tracker, phy_before, phy_after, node, position, tip_labels, allowed_states = NULL, range_candidates = NULL, edge_ages)
{
  info <- get_lineage_info(tracker = tracker, phy = phy_before, nodes = node, edge_ages = edge_ages)
  state <- get_state_at_position(tracker = tracker, phy = phy_before, node = node, position = position, allowed_states = allowed_states, edge_ages = edge_ages)
  inherited_state <- choose_inherited_state(state = state, range_candidates = range_candidates)
  edge_label <- paste0("ms", match(state, tracker$states)) # State of the branch where the subtree is grafted
  missing_label <- paste0("ms", match(inherited_state, tracker$states)) # State of the subtree

  # Root of the subtree and node where it is attached to the branch (its parent)
  n_tips <- ape::Ntip(phy_after)
  subtree_root <- ape::getMRCA(phy_after, tip_labels)
  attach_node <- phy_after$edge[match(subtree_root, phy_after$edge[, 2]), 1]
  subtree_nodes <- phytools::getDescendants(phy_after, subtree_root)
  subtree_nodes <- c(subtree_root, subtree_nodes[subtree_nodes > n_tips])

  # The attachment node belongs to the lineage of the branch (backbone lineage or missing-only subtree), and the nodes of the subtree to a missing-only subtree
  attach_label <- if (!is.na(info$backbone_node)) paste0("gr", info$backbone_node) else edge_label
  phy_after$node.label[attach_node - n_tips] <- attach_label
  phy_after$node.label[subtree_nodes - n_tips] <- missing_label

  tracker$missing_state[tip_labels] <- inherited_state
  list(phy = phy_after, tracker = tracker)
}


#' @title Build the updated simmap of the grafted tree
#'
#' @description At the end of grafting, convert the states followed by the tracker into a stochastic character map (class `"simmap"`) of the grafted tree:
#'   the pieces of the backbone branches (cut by the grafts) keep their states, and every grafted branch is entirely in the state of the point where it was grafted.
#'
#' @param tracker Output of `start_state_tracking()`, updated by `track_graft()` and `track_graft_subtree()` during grafting.
#' @param phy_tracked Final phylogeny of class `"phylo"`, with the nodes labelled by the tracking functions (before the original node labels are restored).
#' @param phy Same phylogeny (same branches in the same order), with the final node labels. The simmap is built on this tree.
#'
#' @noRd
#'
#' @return An object of class `c("simmap", "phylo")`: `phy` with the elements `maps` (list with one named vector per branch, in the order of `phy$edge`,
#'   giving the lengths of the pieces in each state from the rootward to the tipward node) and `mapped.edge` (matrix of the total time in each state per branch).
#'

build_updated_simmap <- function (tracker, phy_tracked, phy)
{
  n_edges <- nrow(phy_tracked$edge)
  edge_length <- phy_tracked$edge.length
  root_age <- max(phytools::nodeHeights(phy_tracked)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy_tracked)
  info <- get_lineage_info(tracker = tracker, phy = phy_tracked, nodes = phy_tracked$edge[, 2], edge_ages = edge_ages)

  ## Pieces of each branch: edge index, state, distance from the tipward node of the end of the piece (to order them), and duration
  # Branches of a missing-only subtree: a single piece in the state of the subtree
  is_missing <- is.na(info$backbone_node)
  idx_missing <- which(is_missing)
  piece_edge <- idx_missing
  piece_state <- info$state[idx_missing]
  piece_from <- rep(0, length(idx_missing))
  piece_duration <- edge_length[idx_missing]

  # Branches of the backbone lineages: pieces of the backbone branch that overlap the extent of the branch
  idx <- which(!is_missing)
  if (length(idx) > 0)
  {
    rows <- tracker$rows_per_node[info$backbone_node[idx]]
    cand <- rep(idx, lengths(rows))
    r <- unlist(rows, use.names = FALSE)
    lo <- info$distance[cand]
    hi <- lo + edge_length[cand]
    # Ends of the pieces are snapped to the extent of the branch, with a relative tolerance (bind.tip() rounds the branch lengths to 10 significant digits)
    eps <- 1e-8 * hi
    f <- ifelse(tracker$state_table$from[r] <= lo + eps, lo, tracker$state_table$from[r])
    t <- ifelse(tracker$state_table$to[r] >= hi - eps, hi, tracker$state_table$to[r])
    f <- pmax(f, lo); t <- pmin(t, hi)
    keep <- (t - f) > eps
    piece_edge <- c(piece_edge, cand[keep])
    piece_state <- c(piece_state, tracker$state_table$state[r][keep])
    piece_from <- c(piece_from, f[keep] - lo[keep])
    piece_duration <- c(piece_duration, (t - f)[keep])
  }

  ## Maps: pieces from the rootward to the tipward node of each branch, rescaled to sum exactly to the branch length
  ord <- order(piece_edge, -piece_from)
  piece_edge <- piece_edge[ord]; piece_state <- piece_state[ord]; piece_duration <- piece_duration[ord]
  if (!identical(sort(unique(piece_edge)), seq_len(n_edges)))
  {
    stop("Internal error: some branches of the grafted tree have no state.")
  }
  durations <- split(stats::setNames(piece_duration, piece_state), factor(piece_edge, levels = seq_len(n_edges)))
  maps <- lapply(seq_len(n_edges), function (e) { d <- durations[[e]]; d * edge_length[e] / sum(d) })
  names(maps) <- NULL

  ## Total time in each state per branch
  # The states of the simmap are all kept, but the ranges added to the states of the tracker are only kept if a branch is in this state
  states <- tracker$states[!(tracker$states %in% tracker$extra_states) | tracker$states %in% unlist(lapply(maps, names), use.names = FALSE)]
  mapped_edge <- matrix(0, nrow = n_edges, ncol = length(states), dimnames = list(paste(phy$edge[, 1], phy$edge[, 2], sep = ","), states))
  cell <- (match(unlist(lapply(maps, names), use.names = FALSE), states) - 1) * n_edges + rep(seq_len(n_edges), lengths(maps))
  total <- tapply(unlist(maps, use.names = FALSE), cell, sum)
  mapped_edge[as.integer(names(total))] <- total

  simmap <- phy
  simmap$maps <- maps
  simmap$mapped.edge <- mapped_edge
  class(simmap) <- c("simmap", "phylo")
  simmap
}
