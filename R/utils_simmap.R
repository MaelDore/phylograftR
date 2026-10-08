
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


#' @title Resolve the allowed states of each taxon
#'
#' @description Get the states in which each taxon can be grafted: the states of `missing_taxa_list$allowed_states` for the taxon if provided,
#'   otherwise the global `allowed_states` argument. A taxon-level value replaces the global one (no intersection).
#'   Same logic as `resolve_include_target_clade_stem()`.
#'   `NULL`, an empty vector or `NA` for a taxon use the global argument. The function stops if some states are not in the simmap,
#'   or if some taxa have no allowed states at all (neither a taxon-level value, nor a global one).
#'
#' @param allowed_states_taxa `NULL`, or list (named after `Taxon_name`) of character vectors: `missing_taxa_list$allowed_states`.
#' @param allowed_states `NULL`, or character vector. The global `allowed_states` argument.
#' @param states Character vector. States of the simmap.
#' @param taxon_names Character vector. `missing_taxa_list$Taxon_name`.
#' @param max_show Integer. Maximum number of taxa listed in the error message.
#'
#' @noRd
#'
#' @return A list named after `taxon_names`, with the allowed states (character vector, no duplicate) of each taxon.
#'

resolve_allowed_states <- function (allowed_states_taxa = NULL, allowed_states = NULL, states, taxon_names, max_show = 10)
{
  n_taxa <- length(taxon_names)
  if (!is.null(allowed_states))
  {
    check_allowed_states(allowed_states, states, who = "the global 'allowed_states' argument")
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
      check_allowed_states(value, states, who = paste0("'missing_taxa_list$allowed_states' of taxon '", taxon_names[i], "'"))
      resolved[[i]] <- unique(value)
    }
  }

  without_states <- which(vapply(resolved, is.null, logical(1)))
  if (length(without_states) > 0)
  {
    n_show <- min(length(without_states), max_show)
    stop("A 'simmap' was provided, so the states in which grafting is allowed are required for every taxon: either the global 'allowed_states' argument, ",
         "or 'missing_taxa_list$allowed_states' for the taxon.\n",
         length(without_states), " taxa have none:\n",
         paste0("  - ", taxon_names[without_states[seq_len(n_show)]], collapse = "\n"),
         if (length(without_states) > n_show) paste0("\n  ... and ", length(without_states) - n_show, " more") else "",
         "\nStates available in the simmap: ", paste0("'", states, "'", collapse = ", "), ".", call. = FALSE)
  }
  resolved
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
    missing_state = character(0)) # State of each grafted taxon (named by its label)

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
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy_before`.
#'
#' @noRd
#'
#' @return A list with 2 elements: the updated `phy` (`phy_after` with its new node labelled) and the updated `tracker`.
#'

track_graft <- function (tracker, phy_before, phy_after, node, position, tip_label, allowed_states = NULL, edge_ages)
{
  info <- get_lineage_info(tracker = tracker, phy = phy_before, nodes = node, edge_ages = edge_ages)
  state <- get_state_at_position(tracker = tracker, phy = phy_before, node = node, position = position, allowed_states = allowed_states, edge_ages = edge_ages)

  # Label of the new node: backbone lineage or missing-only subtree
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

  tracker$missing_state[tip_label] <- state
  list(phy = phy_after, tracker = tracker)
}
