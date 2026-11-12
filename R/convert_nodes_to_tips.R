
##### Convert MRCA node IDs into tips #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - convert_nodes_to_tips(): converts the node IDs of a 'missing_taxa_list' (and of 'total_richness_per_clades') into tips
#  - Internal helpers shared with build_missing_taxa_list(): node_to_tips(), resolve_clade_tokens(), resolve_clade_value(), report_node_conversions()


##### 1/ Exported function #####

#' @title Convert the MRCA node IDs of a `missing_taxa_list` into tips
#'
#' @description Clades (target clades, sister clades, and excluded subclades) of a `missing_taxa_list` can be defined by tips
#'   or by the ID of the node of the backbone phylogeny they correspond to.
#'   This function **converts all the node IDs into tips**, so that the list has the tips-only format required by
#'   [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'   **Tips and node IDs can be mixed, within the list and within the items.** The conversions are printed for transparency:
#'   the warnings and errors of the grafting functions then refer to tips that can be traced back to the node IDs provided.
#'
#'   The same conversion is done internally by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]
#'   and [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] when a table includes node IDs.
#'
#' @param missing_taxa_list A list with the information on the taxa to graft, as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#'   or [`phylograftR::graft_missing_clades()`][graft_missing_clades], except that clades can be defined by node IDs. **Every item provided must be named after `$Taxon_name`.**
#'   Items:
#'   * `$Taxon_name`: Character vector of the names of the taxa to graft.
#'   * `$Target_clade` (or `$Sister_clade`): List (named after `$Taxon_name`) of clades, one per taxon. Each clade is defined by a
#'     character vector of tips, or by a single node ID (see Details).
#'   * `$Excluded_subclades`: (Optional) List (named after `$Taxon_name`) of lists of clades. For each taxon, the subclades of its target clade where the taxon
#'     is known not to belong, each defined by a character vector of tips or by a single node ID. Use `NULL` for taxa with no excluded subclades.
#'     The inner names of the subclades are kept.
#'   * Other items (e.g., `$include_target_clade_stem`, `$exclude_subclades_stem`) are returned unchanged.
#'
#'   It can also be **the list of both types of missing lineages** as `list(taxa = ..., clades = ...)`, as returned by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]
#'   when both types are provided, to use with [`phylograftR::graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades]. Both parts are converted in the same call (see Details).
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which the taxa will be grafted. The node IDs refer to the nodes of this phylogeny.
#' @param total_richness_per_clades (Optional) List with the total richness of clades, as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#'   or [`phylograftR::graft_missing_clades()`][graft_missing_clades]. Its `$Target_clade` (or `$Sister_clade`) clades can also be defined by node IDs, and are converted the same way.
#'   Default is `NULL`. **It must be provided in the same call as its `missing_taxa_list`** (see Details).
#'   Its optional `$Target_clade_name` (or `$Sister_clade_name`) item is returned unchanged.
#'   With a `list(taxa, clades)` as `missing_taxa_list`, it must be a `list(taxa = ..., clades = ...)` too, and each part must match a type present in `missing_taxa_list`.
#' @param verbose Logical. Whether to print the conversions of the node IDs into tips. Default is `TRUE`.
#'
#' @details
#' ### Node IDs and tips
#'
#' Node IDs are integers: `1` to N, with N being the number of tips (`ape::Ntip(phy)`), correspond to terminal nodes, and thus define single-tip clades.
#' Following IDs are the internal nodes, as shown by `ape::nodelabels()`, and correspond to clades with multiple descending species.
#' A node ID is recognized as being exclusively made of digits: either a number (e.g., `34`), or a text with only digits (e.g., `"34"`, as read from a table).
#' **Tip labels must therefore include at least one other character than digits** to be distinguished.
#' The function stops if a node ID is used while some tip labels of `phy` are only made of digits. In this case, you need to rename these tips first.
#'
#' **A clade is defined either by tips, or by one node ID. The two cannot be mixed within the same clade** (e.g., `c("t8", "34")` is an error),
#' and a node ID alone defines one clade (several node IDs, e.g., `c(34, 35)`, are an error). However, different clades can be defined differently.
#' For example, `Excluded_subclades = list(X1 = list(50, c("t5", "t12"), "t30"))` mixes an internal node ID, a pair of tips and a single tip.
#'
#' Single tips and clades with mutiple descending species are distinguished by the type of node: a tip node (`<= ape::Ntip(phy)`) is equivalent to the single tip,
#' i.e., the terminal branch of the tip, and an internal node is equivalent to the clade it subtends.
#'
#' ### Conversion
#'
#' A **tip node** is converted into the label of the tip. An **internal node** is converted into one pair of tips whose most recent common ancestor (MRCA) is the node.
#' The pair is chosen in deterministic fashion: for each of the first two descendant lineages of the node (in the order of `phy$edge`), the first tip of each lineage is selected.
#' This ensures that **the MRCA node of this pair always defines the same clade**, even when other species are grafted inside.
#' The conversions are printed, grouped by node ID, with the items and taxa where each node ID was used.
#'
#' ### Same clade, same tips
#'
#' If a clade is already defined by tips in `missing_taxa_list` (written by hand, or converted before), a node ID with the same MRCA node is converted
#' into these tips instead of the default pair, so the same clade is always defined by the same tips. This holds for `missing_taxa_list` and `total_richness_per_clades`.
#' `total_richness_per_clades` **must be converted at the same time** as its associated `missing_taxa_list` (same call to the function), so that the tips defining a clade match between the two objects.
#' It can be converted together with the initial `missing_taxa_list` (both converted at the same time), or together with the `missing_taxa_list` that was already converted,
#' used as a template (the tips it already uses for the clades are matched to the node IDs).
#'
#' ### Both types of taxa
#'
#' With a `list(taxa, clades)` provided as `missing_taxa_list`, the two parts (and their total richness) are converted in one call with a shared record of the clades: a clade is defined by the same tips
#' wherever it appears (target clades of the taxa, sister clades of the missing clades, excluded subclades, total richness), so that **all the objects remain consistent when used together**
#' with [`phylograftR::graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades]. The tips already used for a clade are those of the first part where it appears (`taxa` before `clades`).
#' The conversions of both parts are printed together.
#'
#' The function also checks that the structure of the list is valid (items named after `$Taxon_name`), that node IDs are in `phy`, and that the tips are in `phy`.
#' It stops with an informative message otherwise.
#'
#' @return A list with the same items as `missing_taxa_list`, with all the clades defined by tips.
#' The items, the names of the items and of the subclades are unchanged.
#' If `total_richness_per_clades` is provided, a list with two items:
#'   * `$missing_taxa_list` The converted `missing_taxa_list`.
#'   * `$total_richness_per_clades` The converted `total_richness_per_clades`.
#'
#' **With a `list(taxa, clades)`, the output has the same structure as the input**: a list with the elements `taxa` and `clades`,
#' each converted (and, if `total_richness_per_clades` is provided, the two items above, each being a `list(taxa, clades)`).
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] to build a `missing_taxa_list` from an internal data.frame, with the same conversion of node IDs.
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read a table and build the list directly from an external file, with the same conversion of node IDs.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa based on target clades and excluded subclades using the list.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades using the list.
#'
#' @importFrom ape Ntip
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Load the backbone phylogeny
#' data(toy_backbone)
#'
#' # Plot it with the node IDs
#' ape::plot.phylo(toy_backbone, cex = 0.7, no.margin = TRUE)
#' ape::nodelabels(cex = 0.6, frame = "none", col = "blue")
#'
#' # A 'missing_taxa_list' with clades defined by node IDs, tips, or both
#' my_missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3"),
#'   Target_clade = list(X1 = 34, X2 = c("t5", "t12"), X3 = "t26"),
#'   Excluded_subclades = list(X1 = list(group_1 = 50), X2 = NULL, X3 = NULL))
#'
#' # Convert the node IDs into tips
#' tips_missing_taxa_list <- convert_nodes_to_tips(my_missing_taxa_list, phy = toy_backbone)
#' str(tips_missing_taxa_list, 2)
#'
#' # The total richness of the clades can be converted too
#' total_richness_per_clades <- list(Target_clade = list(32, c("t5", "t12")),
#'                                   total_richness = c(25, 8))
#' converted <- convert_nodes_to_tips(my_missing_taxa_list, phy = toy_backbone,
#'                                    total_richness_per_clades = total_richness_per_clades)
#' converted$total_richness_per_clades
#'
#' # Node 34 is the same clade in both objects, and is defined by the same tips
#' converted$missing_taxa_list$Target_clade$X1
#' # An already converted list can be used:
#' # the node IDs of the total richness are matched to the tips in the list
#' converted_again <- convert_nodes_to_tips(converted$missing_taxa_list, phy = toy_backbone,
#'                                          total_richness_per_clades = list(Target_clade = list(34),
#'                                                                           total_richness = 14))
#' converted_again$total_richness_per_clades$Target_clade
#'
#' # --- Convert nodes ID for both types of missing lineages --- #
#'
#' # With list(taxa, clades), both parts are converted at once, with a shared record of the clades
#' both_list <- list(
#'   taxa = my_missing_taxa_list,
#'   clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = 34, M2 = "t26")))
#'
#' converted_both <- convert_nodes_to_tips(both_list, phy = toy_backbone)
#' names(converted_both)
#'
#' # Node 34 is defined by the same tips in the target clade of X1 and the sister clade of M1
#' converted_both$taxa$Target_clade$X1
#' converted_both$clades$Sister_clade$M1
#'

convert_nodes_to_tips <- function (missing_taxa_list, # List of missing taxa with clades defined by tips and/or node IDs, or list(taxa, clades)
                                   phy, # Backbone phylogeny
                                   total_richness_per_clades = NULL, # Optional list with the total richness of clades (list(taxa, clades) with both types)
                                   verbose = TRUE) # Print the conversions
{
  ## 1/ Check the inputs ####

  if (!inherits(phy, "phylo")) { stop("'phy' must be a phylogeny of class 'phylo'") }
  if (!is.list(missing_taxa_list)) { stop("'missing_taxa_list' must be a list") }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) { stop("'verbose' must be TRUE or FALSE") }

  ## 2/ Identify the types of taxa ####

  # Type ("taxa" or "clades") of a single-type list, from its clade item (an invalid list is left to the checks of the conversion)
  type_of <- function (x, items) { if (!is.null(x[[items[1]]]) && is.null(x[[items[2]]])) "taxa" else if (is.null(x[[items[1]]]) && !is.null(x[[items[2]]])) "clades" else "taxa" }
  typed <- split_typed_list(missing_taxa_list, arg = "missing_taxa_list")
  as_typed <- function (x) { out <- list(x); names(out) <- type_of(x, c("Target_clade", "Sister_clade")); out }
  lists <- if (typed$typed) typed$lists else as_typed(missing_taxa_list)

  richness <- NULL
  if (!is.null(total_richness_per_clades))
  {
    if (!is.list(total_richness_per_clades)) { stop("'total_richness_per_clades' must be a list") }
    tr_typed <- split_typed_list(total_richness_per_clades, arg = "total_richness_per_clades")
    richness <- if (tr_typed$typed) tr_typed$lists else as_typed(total_richness_per_clades)
    if (!all(names(richness) %in% names(lists)))
    {
      stop("'total_richness_per_clades' includes the type '", paste(setdiff(names(richness), names(lists)), collapse = "', '"), "' that is not in 'missing_taxa_list'.", call. = FALSE)
    }
    if (length(lists) == 2 && !tr_typed$typed)
    {
      stop("'missing_taxa_list' includes both types of taxa: 'total_richness_per_clades' must be a list(taxa = , clades = ) too.", call. = FALSE)
    }
  }

  ## 3/ Convert each type, with a shared record of the clades ####

  # Record of the conversions, and of the clades already defined by tips (written by hand, or converted earlier) in 'missing_taxa_list':
  # a node ID is converted into the tips already used for this node, so the same clade is always defined by the same tips, in all the objects
  recorder <- new.env(parent = emptyenv())
  recorder$rows <- list()
  known <- unlist(lapply(lists, function (x)
  {
    k <- list(if (!is.null(x$Target_clade)) x$Target_clade else x$Sister_clade)
    if (is.list(x$Excluded_subclades)) { k <- c(k, unlist(lapply(x$Excluded_subclades, function (y) if (is.list(y)) list(y) else NULL), recursive = FALSE)) }
    unlist(k, recursive = FALSE)
  }), recursive = FALSE)
  lookup <- build_tips_lookup(phy, known)

  converted <- lapply(names(lists), function (type)
  {
    convert_nodes_to_tips_typed(missing_taxa_list = lists[[type]], phy = phy, total_richness_per_clades = richness[[type]], recorder = recorder, lookup = lookup)
  })
  names(converted) <- names(lists)

  ## 4/ Report the conversions ####

  n_converted <- report_node_conversions(recorder$rows, phy = phy, verbose = verbose)
  if (verbose && n_converted == 0)
  {
    message("No node ID found: all the clades are already defined by tips.")
  }

  ## 5/ Output with the structure of the input ####

  if (length(converted) == 1)
  {
    if (is.null(richness)) { return(converted[[1]]$missing_taxa_list) }
    return(list(missing_taxa_list = converted[[1]]$missing_taxa_list, total_richness_per_clades = converted[[1]]$total_richness_per_clades))
  }
  converted_lists <- lapply(converted, function (x) x$missing_taxa_list)
  if (is.null(richness)) { return(converted_lists) }
  converted_richness <- lapply(converted[names(richness)], function (x) x$total_richness_per_clades)
  return(list(missing_taxa_list = converted_lists, total_richness_per_clades = converted_richness))
}


##### 2/ Internal helpers #####

## Convert the node IDs of one type of taxa (see convert_nodes_to_tips())
# 'recorder': environment recording the conversions (shared by the types); 'lookup': tips already used for the nodes (shared by the types)
convert_nodes_to_tips_typed <- function (missing_taxa_list, phy, total_richness_per_clades, recorder, lookup)
{
  ## 1/ Check the inputs ####

  has_target <- !is.null(missing_taxa_list$Target_clade)
  has_sister <- !is.null(missing_taxa_list$Sister_clade)
  if (has_target && has_sister) { stop("'missing_taxa_list' has both 'Target_clade' and 'Sister_clade': use only one, depending on the grafting function (graft_missing_taxa() or graft_missing_clades())") }
  if (!has_target && !has_sister) { stop("'missing_taxa_list' must have an item 'Target_clade' (for graft_missing_taxa()) or 'Sister_clade' (for graft_missing_clades())") }
  clade_item <- if (has_target) "Target_clade" else "Sister_clade"

  # Check structure: items must be named after Taxon_name
  validate_missing_taxa_list(missing_taxa_list, clade_item = clade_item, allow_node_ids = TRUE)
  taxon_name <- missing_taxa_list$Taxon_name

  ## 2/ Convert the clades of the taxa ####

  clades <- missing_taxa_list[[clade_item]]
  if (!is.list(clades)) { stop("'missing_taxa_list$", clade_item, "' must be a list with one clade (vector of tips, or a node ID) per taxon") }
  missing_clades <- which(vapply(clades, is.null, logical(1)))
  if (length(missing_clades) > 0) { stop("'missing_taxa_list$", clade_item, "' is empty for: ", paste(taxon_name[missing_clades], collapse = ", ")) }
  new_clades <- lapply(seq_along(clades), function (i)
  {
    resolve_clade_value(clades[[i]], phy = phy, item = clade_item, who = taxon_name[i], recorder = recorder, lookup = lookup)
  })
  names(new_clades) <- names(clades)
  missing_taxa_list[[clade_item]] <- new_clades
  check_tips_in_phy(new_clades, phy, taxon_name, clade_item)

  ## 3/ Convert the excluded subclades ####

  if (!is.null(missing_taxa_list$Excluded_subclades))
  {
    excluded <- missing_taxa_list$Excluded_subclades
    if (!is.list(excluded)) { stop("'missing_taxa_list$Excluded_subclades' must be a list with one list of subclades (or NULL) per taxon") }
    new_excluded <- vector("list", length(excluded))
    for (i in seq_along(excluded))
    {
      subclades <- excluded[[i]]
      if (is.null(subclades)) { next }
      if (!is.list(subclades))
      {
        stop("'missing_taxa_list$Excluded_subclades' of '", taxon_name[i], "' must be a list of subclades (each a vector of tips, or a node ID), or NULL.\n",
             "Use e.g. list(c(\"t1\", \"t2\"), 34) for two subclades.")
      }
      new_excluded[[i]] <- lapply(subclades, function (x)
      {
        resolve_clade_value(x, phy = phy, item = "Excluded_subclades", who = taxon_name[i], recorder = recorder, lookup = lookup)
      })
    }
    names(new_excluded) <- names(excluded)
    missing_taxa_list["Excluded_subclades"] <- list(new_excluded)
    check_tips_in_phy(lapply(new_excluded, function (x) unique(unlist(x))), phy, taxon_name, "Excluded_subclades")
  }

  ## 4/ Convert the clades of the total richness ####

  if (!is.null(total_richness_per_clades))
  {
    if (!is.list(total_richness_per_clades)) { stop("'total_richness_per_clades' must be a list") }
    has_target_tr <- !is.null(total_richness_per_clades[["Target_clade"]])
    has_sister_tr <- !is.null(total_richness_per_clades[["Sister_clade"]])
    if (has_target_tr == has_sister_tr)
    {
      stop("'total_richness_per_clades' must have one item 'Target_clade' (or 'Sister_clade'), and not both")
    }
    tr_item <- if (has_target_tr) "Target_clade" else "Sister_clade"
    tr_clades <- total_richness_per_clades[[tr_item]]
    if (!is.list(tr_clades)) { stop("'total_richness_per_clades$", tr_item, "' must be a list of clades (vectors of tips, or node IDs)") }
    new_tr_clades <- lapply(seq_along(tr_clades), function (k)
    {
      resolve_clade_value(tr_clades[[k]], phy = phy, item = paste0("total_richness_per_clades$", tr_item),
                          who = paste0("clade ", k), recorder = recorder, lookup = lookup)
    })
    names(new_tr_clades) <- names(tr_clades)
    check_tips_in_phy(new_tr_clades, phy, paste0("clade ", seq_along(new_tr_clades)), paste0("total_richness_per_clades$", tr_item))
    total_richness_per_clades[[tr_item]] <- new_tr_clades
  }

  return(list(missing_taxa_list = missing_taxa_list, total_richness_per_clades = total_richness_per_clades))
}


#' @title Convert a node ID into tips
#'
#' @description A tip node gives its tip label. An internal node gives one pair of tips whose MRCA is the node:
#'   the first tip of each of the first two descendant lineages of the node, in the order of `phy$edge` (deterministic).
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param node Integer. Node ID, between 1 and `Ntip + Nnode`.
#'
#' @noRd
#'
#' @return Character vector of one tip (tip node) or two tips (internal node).
#'

node_to_tips <- function (phy, node)
{
  n_tip <- ape::Ntip(phy)
  if (node <= n_tip) { return(phy$tip.label[node]) }

  kids <- phy$edge[phy$edge[, 1] == node, 2]
  if (length(kids) < 2) { stop("The node ", node, " has fewer than two descendant lineages: it cannot be defined by tips") }
  # First tip of a lineage (follow the first descendant until a tip)
  first_tip <- function (nd)
  {
    while (nd > n_tip) { nd <- phy$edge[match(nd, phy$edge[, 1]), 2] }
    nd
  }
  tips <- phy$tip.label[vapply(kids[1:2], first_tip, numeric(1))]

  # Safety: the pair must define the node
  if (!identical(as.numeric(ape::getMRCA(phy, tips)), as.numeric(node)))
  {
    stop("Internal error: the tips ", paste(tips, collapse = ", "), " do not define the node ", node)
  }
  return(tips)
}


#' @title Resolve the tokens of a clade into tips
#'
#' @description Tokens (tips, or node ID written as digits) of ONE clade are converted into tips.
#'   Tokens are either all tips (returned as is), or one node ID (converted).
#'   A mix of tips and node IDs, or several node IDs, stops with an informative message.
#'
#' @param tokens Character vector (or numeric). Tokens of the clade.
#' @param phy Phylogeny of class `"phylo"`.
#' @param item Character string. Item or column, for messages and records.
#' @param who Character string. Taxon (or clade), for messages and records.
#' @param recorder Environment with a list `rows` where conversions are recorded.
#' @param lookup Named list (names are node IDs) of the tips already used for a node: used instead of the default pair of tips. Default is `NULL`.
#' @param cell Character string. Optional text of the cell, for messages.
#'
#' @noRd
#'
#' @return Character vector of tips. Unchanged tokens if there is no node ID.
#'

resolve_clade_tokens <- function (tokens, phy, item, who, recorder, cell = NULL, lookup = NULL)
{
  shown <- if (is.null(cell)) paste(tokens, collapse = ", ") else cell
  where <- paste0("'", item, "' of '", who, "'")

  if (length(tokens) == 0 || anyNA(tokens)) { stop(where, " has a missing value in a clade: ", shown) }

  if (is.numeric(tokens))
  {
    is_node <- rep(TRUE, length(tokens))
  } else {
    tokens <- trimws(as.character(tokens))
    is_node <- grepl("^[0-9]+$", tokens)
  }

  # Only tips
  if (!any(is_node)) { return(tokens) }

  # Mixed or multiple node IDs
  if (!all(is_node))
  {
    stop(where, " mixes node IDs and tips in the same clade: ", shown, ".\n",
         "A clade is defined either by tips (labels with at least one other character than digits), or by one node ID, not both.")
  }
  if (length(tokens) > 1)
  {
    stop(where, " has several node IDs in the same clade: ", shown, ".\n",
         "A clade defined by a node ID needs exactly one node ID. To define a clade by several tips, use their labels.")
  }

  # One node ID
  node <- suppressWarnings(as.numeric(tokens))
  n_node_max <- ape::Ntip(phy) + phy$Nnode
  if (is.na(node) || node != round(node) || node < 1 || node > n_node_max)
  {
    stop(where, " has an invalid node ID: ", shown, ". Node IDs of 'phy' are integers between 1 and ", n_node_max,
         " (1 to ", ape::Ntip(phy), " for the tips).")
  }
  numeric_tips <- grep("^[0-9]+$", phy$tip.label, value = TRUE)
  if (length(numeric_tips) > 0)
  {
    stop(where, " uses a node ID (", shown, "), but some tip labels of 'phy' are only made of digits (e.g., ",
         paste(utils::head(numeric_tips, 3), collapse = ", "), "): node IDs and these tips cannot be distinguished.\n",
         "Rename these tips with at least one other character than digits (e.g., 'sp_", numeric_tips[1], "').")
  }

  # Tips already used for this node in the objects converted together, otherwise the default tips of the node
  reused <- !is.null(lookup[[as.character(node)]])
  tips <- if (reused) lookup[[as.character(node)]] else node_to_tips(phy, node)
  recorder$rows[[length(recorder$rows) + 1]] <- data.frame(node = node, tips = paste(tips, collapse = ", "), reused = reused,
                                                           item = item, who = who, stringsAsFactors = FALSE)
  return(tips)
}


#' @title Resolve a clade of a list into tips
#'
#' @param x Character vector of tips, or node ID (number, or text of digits), defining one clade.
#' @inheritParams resolve_clade_tokens
#'
#' @noRd
#'
#' @return Character vector of tips, keeping no names.
#'

resolve_clade_value <- function (x, phy, item, who, recorder, lookup = NULL)
{
  if (is.null(x) || length(x) == 0) { stop("'", item, "' of '", who, "' has an empty clade") }
  if (is.factor(x)) { x <- as.character(x) }
  if (is.list(x)) { stop("'", item, "' of '", who, "' has a clade defined by a list: use a character vector of tips, or one node ID") }
  if (!is.numeric(x) && !is.character(x))
  {
    stop("'", item, "' of '", who, "' must define its clades with tips (character) or node IDs (numbers), not '", class(x)[1], "'")
  }
  resolve_clade_tokens(unname(x), phy = phy, item = item, who = who, recorder = recorder, lookup = lookup)
}


#' @title Print the conversions of node IDs into tips
#'
#' @param rows List of one-row data frames (`node`, `tips`, `reused`, `item`, `who`).
#' @param phy Phylogeny of class `"phylo"`.
#' @param verbose Logical. Whether to print.
#' @param max_who Integer. Maximum number of taxa listed per item.
#'
#' @noRd
#'
#' @return Integer. Number of distinct node IDs converted (invisibly printed if `verbose`).
#'

report_node_conversions <- function (rows, phy, verbose = TRUE, max_who = 4)
{
  if (length(rows) == 0) { return(0L) }
  df <- do.call(rbind, rows)
  nodes <- unique(df$node)
  if (!verbose) { return(length(nodes)) }

  n_tip <- ape::Ntip(phy)
  lines <- vapply(nodes, function (nd)
  {
    sub <- df[df$node == nd, , drop = FALSE]
    if (nd <= n_tip)
    {
      type <- "tip node"
    } else {
      type <- paste0("internal node, ", sum(get_descendant_nodes(phy, nd) <= n_tip), " tips")
    }
    used <- vapply(unique(sub$item), function (it)
    {
      w <- unique(sub$who[sub$item == it])
      shown <- paste(utils::head(w, max_who), collapse = ", ")
      if (length(w) > max_who) shown <- paste0(shown, ", ... (", length(w) - max_who, " more)")
      paste0(it, ": ", shown)
    }, character(1))
    note <- if (any(sub$reused)) " (tips already used in the list)" else ""
    paste0("  - node ", nd, " (", type, ") -> ", sub$tips[1], note, "   [", paste(used, collapse = " | "), "]")
  }, character(1))

  message(length(nodes), " node ID(s) converted into tips (a tip node gives its tip, an internal node gives one pair of tips with this MRCA node):\n",
          paste(lines, collapse = "\n"))
  return(length(nodes))
}


#' @title Tips already used for the nodes of clades
#'
#' @description Associate the node of each clade defined by tips to the tips used. The first clade found is kept for a node.
#'   Clades defined by a node ID, or by tips absent from `phy`, are skipped (they are checked and converted elsewhere).
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param clades List of clades (character vectors of tips, or node IDs).
#'
#' @noRd
#'
#' @return Named list (names are node IDs) of character vectors of tips. Empty list if none.
#'

build_tips_lookup <- function (phy, clades)
{
  lookup <- list()
  for (x in clades)
  {
    if (!is.character(x) || length(x) == 0 || anyNA(x) || is_node_id_clade(x) || !all(x %in% phy$tip.label)) { next }
    node <- as.character(get_clade_node(phy, unname(x)))
    if (is.null(lookup[[node]])) { lookup[[node]] <- unname(x) }
  }
  lookup
}
