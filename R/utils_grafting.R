
##### Internal helpers shared by graft_missing_taxa() and graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

## Utility functions used internally by grafting functions


#' @title Get all descendant nodes of a node
#'
#' @description Fast, non-recursive search for descendants nodes.
#'   Include both internal and terminal nodes.
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param node Integer. Node ID.
#'
#' @noRd
#'
#' @return Integer vector of descendant node IDs including both internal and terminal nodes.
#'

get_descendant_nodes <- function (phy, node)
{
  if (length(node) == 0) return(integer(0))
  edge <- phy$edge
  out <- integer(0)
  frontier <- node
  while (length(frontier) > 0)
  {
    kids <- edge[edge[, 1] %in% frontier, 2]
    out <- c(out, kids)
    frontier <- kids
  }
  out
}


#' @title Get the node ID of a clade defined by tips
#'
#' @description Identify node ID associated to a clade based on tips:
#'  * If defined by one tip, returns the tip node
#'  * If defined by two or more tips, returns the MRCA node
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param tips Character vector. One tip (tip node) or two or more tips (MRCA node).
#'
#' @importFrom ape getMRCA
#'
#' @noRd
#'
#' @return Integer. Node ID.

get_clade_node <- function (phy, tips)
{
  tips <- unname(tips)

  if (length(tips) == 1)
  {
    node_ID <- match(tips, phy$tip.label)
  } else {
    node_ID <- ape::getMRCA(phy = phy, tip = tips)
  }

  return(node_ID)
}


#' @title Get the parent node of a node
#'
#' @description Identify node ID of the parent node of a node.
#'   A parent node is defined as the rootward node of the branch leading to the focal node (i.e., stem branch).
#'   Return `NA` for the root as it does not have a stem branch.
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param node Integer. Node ID.
#'
#' @noRd
#'
#' @return Integer node ID, or `NA` for the root.
#'

get_parent_node <- function (phy, node)
{
  parent_node <- phy$edge[phy$edge[, 2] == node, 1]
  if (length(parent_node ) == 0)
  {
    parent_node <- NA
  }
  return(parent_node)
}


#' @title Identify the branches of the stem region of a clade
#'
#' @description Identify the "stem region" of a clade: the set of branches where a taxon can be grafted
#'   to become sister to the clade (i.e., grafted on its stem).
#'   Because grafting splits branches and adds new branches, the stem region is not just the original stem branch.
#'   It is made of:
#'   * The stem branch of the clade, in as many pieces as it has been split by previous grafts.
#'   * The branches of 'missing-only subtrees' attached to the stem, i.e., subtrees hanging from the stem
#'     that contain only missing taxa (taxa listed for grafting, never backbone tips).
#'     Those subtrees are made of taxa previously grafted on the stem, and of taxa later grafted on them.
#'
#'   The region is explored by walking rootward from the clade node. At each step, the parent node is a "stem-insertion node"
#'   if all the other children of that parent are missing-only subtrees. The walk stops at the first parent node
#'   with a side subtree containing at least one backbone tip: this is the parent node of the clade in the backbone phylogeny.
#'
#' @param phy Phylogeny of class `"phylo"`, possibly including taxa already grafted.
#' @param node Integer. Node ID of the clade (tip node for a single-tip clade).
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not).
#'   Tips with these labels are considered as missing taxa, all the others as backbone tips.
#'
#' @importFrom ape Ntip
#'
#' @noRd
#'
#' @return A list with 2 elements:
#'  * `nodes` Integer vector. Tipward node IDs of all the branches of the stem region. Empty for the root, which has no stem.
#'  * `top_parent` Integer. Node ID of the parent of the stem in the backbone phylogeny, i.e., the rootward end of the stem region
#'    (`NA` for the root).
#'

get_stem_region_nodes <- function (phy, node, missing_tip_labels = character(0))
{
  # The root has no stem
  if (is.na(get_parent_node(phy, node)))
  {
    return(list(nodes = integer(0), top_parent = NA))
  }

  n_tips <- ape::Ntip(phy)

  # A subtree is 'missing-only' if all its tips are missing taxa (i.e., it contains no backbone tips).
  is_missing_only <- function (subtree_node)
  {
    subtree_nodes <- c(subtree_node, get_descendant_nodes(phy, subtree_node))
    subtree_tip_nodes <- subtree_nodes[subtree_nodes <= n_tips]
    all(phy$tip.label[subtree_tip_nodes] %in% missing_tip_labels)
  }

  # The branch leading to the clade node is always part of the stem region
  # (if the stem has been split by a graft, this is its most tipward piece)
  stem_nodes <- node
  current_node <- node

  repeat
  {
    parent_node <- get_parent_node(phy, current_node)

    # Reaching the root = no backbone tip found on the other side (not expected with a valid backbone)
    if (is.na(parent_node))
    {
      top_parent <- NA
      break
    }

    # Subtrees hanging from the stem at this parent node (the children of the parent, but the current node)
    side_nodes <- setdiff(phy$edge[phy$edge[, 1] == parent_node, 2], current_node)
    side_is_missing_only <- length(side_nodes) > 0 && all(vapply(side_nodes, is_missing_only, logical(1)))

    # If a backbone tip is found in a side subtree, the parent node is a node of the backbone phylogeny:
    # the walk is over, and this node is the true parent of the stem
    if (!side_is_missing_only)
    {
      top_parent <- parent_node
      break
    }

    # Otherwise, the parent node was created by a graft on the stem.
    # Its side subtrees (made of missing taxa only) are part of the stem region: add all their branches
    stem_nodes <- c(stem_nodes, side_nodes, unlist(lapply(side_nodes, function (n) get_descendant_nodes(phy, n))))
    # The branch above the parent node is the next piece of the split stem (if the parent is not the root)
    if (!is.na(get_parent_node(phy, parent_node))) { stem_nodes <- c(stem_nodes, parent_node) }
    # Keep walking rootward
    current_node <- parent_node
  }

  return(list(nodes = stem_nodes, top_parent = top_parent))
}


#' @title Identify the branches available for grafting within the target clade
#'
#' @description Identify the tipward node IDs of all the branches available for grafting a taxon within a target clade,
#'   before removing the branches of the excluded subclades:
#'   * All branches within the clade (descending from its MRCA node).
#'   * The stem region of the clade (see `get_stem_region_nodes()`), if `include_target_clade_stem = TRUE`.
#'     A taxon grafted on the stem region becomes sister to the clade, or sister to other missing taxa already grafted on the stem.
#'
#'   Single-tip clades always use their stem region (their tip branch is their only branch), whatever the value of `include_target_clade_stem`.
#'   The root has no stem, thus no stem region.
#'
#' @param phy Phylogeny of class `"phylo"`, possibly including taxa already grafted.
#' @param MRCA_node Integer. Node ID of the clade (tip node for a single-tip clade).
#' @param include_target_clade_stem Logical. Whether the stem region of the clade is available. Default is `TRUE`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not). Default is `character(0)`.
#'
#' @importFrom ape Ntip
#'
#' @noRd
#'
#' @return Integer vector of tipward node IDs.
#'

get_available_branch_nodes <- function (phy, MRCA_node, include_target_clade_stem = TRUE, missing_tip_labels = character(0))
{
  # All branches within the clade
  available_nodes <- get_descendant_nodes(phy = phy, node = MRCA_node)

  # Single-tip clades: the tip branch is the only branch available (it is the stem of the clade), so include_target_clade_stem is ignored (set to TRUE)
  if (MRCA_node <= ape::Ntip(phy)) { include_target_clade_stem <- TRUE }

  # Add the stem region if requested. It is empty for the root, which has no stem
  if (include_target_clade_stem)
  {
    available_nodes <- c(available_nodes, get_stem_region_nodes(phy = phy, node = MRCA_node, missing_tip_labels = missing_tip_labels)$nodes)
  }

  return(available_nodes)
}


#' @title Sample randomly one element from a vector
#'
#' @description Sample randomly one element from a vector in a safe fashion that
#'   handles vector of length 0 or 1. Return an error if the vector is empty.
#'
#' @param x Vector to sample from.
#' @param prob Numeric vector. Optional sampling probabilities. Default is `NULL`.
#'
#' @noRd
#'
#' @return One element of `x`. Stops if `x` is empty.
#'

sample_one <- function (x, prob = NULL)
{
  if (length(x) == 0) stop("No valid option to sample from")
  if (length(x) == 1) return(x)
  x[sample.int(length(x), size = 1, prob = prob)]
}


#' @title Fit a birth-death model on a clade enforcing a minimum number of tips
#'
#' @description Fit a birth-death model on a clade enforcing a minimum number of tips.
#'   If the clade has fewer tips than `min_tips`, climb towards the root until it has enough (stops at the root in any case).
#'   This ensures the enclosing clade chosen has enough divergence times to be able to fit a relevant BD model,
#'   while retaining close relationship with the initial focal clade.
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param node Integer. Node ID of the clade.
#' @param min_tips Integer. Minimum number of tips for the fit. Default is `5`.
#' @param richness_info `NULL`, or the output of `build_richness_info()`. If provided, the sampling fraction of the clade used for the fit is
#'   `n_present / (n_present + n_pending + n_unsampled)`, with:
#'   * `n_present`: number of tips of the clade in the current tree, i.e., backbone tips plus taxa already grafted (they are not 'missing' anymore).
#'   * `n_pending`: number of taxa still to graft that belong to the clade.
#'   * `n_unsampled`: number of unsampled species of the clade, derived from `total_richness_per_clades`.
#'   Default is `NULL`: no taxa to graft and no unsampled species are accounted for, so the sampling fraction is 1.
#'   The birth-death model is always fitted on the crown group of the selected clade, with a sampling fraction describing this same crown group.
#'
#' @details The goal is to use the fitted rates to simulate divergence ages of missing taxa using `TreeSim::corsim()`.
#'   Because no new divergence times can be drawn from a negative net diversification process,
#'   the function ensures that the net diversification rate is always positive by capping the death rate to 99% of the birth rate.
#'   A death rate that is not finite or negative is set to zero.
#'
#' @importFrom ape Ntip extract.clade
#'
#' @importFrom phytools fit.bd
#'
#' @noRd
#'
#' @return A list with 7 elements:
#'  * `b` Numeric. Birth rate = speciation rate in events per time per lineage.
#'  * `d` Numeric. Death rate = extinction rate in events per time per lineage.
#'  * `n_tips` Integer. Number of tips in the enclosing clade used for fitting.
#'  * `node` Integer. MRCA node ID of the selected enclosing clade.
#'  * `n_pending` Numeric. Number of taxa still to graft in the selected clade, used for the sampling fraction (0 if `richness_info = NULL`).
#'  * `n_unsampled` Numeric. Number of unsampled species used for the sampling fraction.
#'  * `rho` Numeric. Sampling fraction used for the fit.
#'

fit_bd_on_clade <- function (phy, node, min_tips = 5, richness_info = NULL)
{
  clade_node <- node

  ## Case with node as a tip
  # extract.clade() does not accept tips: start from the parent node instead
  if (clade_node <= ape::Ntip(phy))
  {
    clade_node <- get_parent_node(phy, clade_node)
  }

  ## Search rootward for an enclosing clade with at least 'min_tips' tips
  # Ensure that the BD model is fit on a clade with at least 'min_tips' tips
  repeat
  {
    # Extract current target clade
    clade_tree <- ape::extract.clade(phy = phy, node = clade_node)

    # Stop once the current clade has enough tips
    if (ape::Ntip(clade_tree) >= min_tips) break
    # If not, move rootward to the next enclosing parent clade
    parent_node <- get_parent_node(phy, clade_node)
    # If reaching the root (i.e., parent = NA), stop.
    if (is.na(parent_node)) break
    # Replace current target clade by its next enclosing clade / parent clade
    clade_node <- parent_node
  }

  # Extract final nb of tips in the selected enclosing clade
  n_tips <- ape::Ntip(clade_tree)

  # Stop if no clade can cover 'min_tips'
  if (n_tips < min_tips)
  {
    stop("Backbone phylogeny is too small to fit a birth-death model (< ",min_tips," tips)")
  }

  ## Fit BD with sampling fraction accounting for missing tips (the tips we are aiming to graft) and unsampled species (listed neither in the backbone nor in missing species)
  # Numbers of taxa still to graft and of unsampled species in the selected clade (none by default)
  n_pending <- 0
  n_unsampled <- 0
  if (!is.null(richness_info))
  {
    # Locate the selected clade in the backbone: its backbone tips (grafted taxa removed) define the same clade in the original backbone
    backbone_tips <- setdiff(clade_tree$tip.label, richness_info$taxa_names)
    node_in_backbone <- get_clade_node(richness_info$phy_backbone, backbone_tips)
    # Missing taxa still to graft that belong to this clade. Taxa already grafted are in the tree: they are in 'n_tips', not missing anymore
    still_to_graft <- !(richness_info$taxa_names %in% phy$tip.label)
    in_clade <- vapply(richness_info$taxon_sets, function (s) node_in_backbone %in% s, logical(1))
    n_pending <- sum(still_to_graft & in_clade)
    # Unsampled species of the clade, derived from the total richness
    n_unsampled <- richness_info$unsampled_by_node[node_in_backbone]
  }
  # Sampling fraction of the crown group of the selected clade (1 if nothing is pending nor unsampled in it)
  rho <- n_tips / (n_tips + n_pending + n_unsampled)
  clade_tree <- ensure_ultrametry(clade_tree) # Ensure the selected clade is ultrametric
  fit <- phytools::fit.bd(tree = clade_tree, rho = rho)

  ## Extract rates
  b <- unname(fit$b); d <- unname(fit$d)

  ## Keep rates usable by TreeSim (mu < lambda)
  # If death rate is infinite or negative, fix it to zero
  if (!is.finite(d) || d < 0) d <- 0
  d <- min(d, 0.99 * b) # Limit death rate to 99% of birth rate to ensure the net diversification is never negative
  # No divergence times can be drawn from a negative net diversification process

  # Record rates, nb of tips in selected enclosing clade, and MRCA node of the selected clade
  BD_output <- list(b = b, d = d, n_tips = n_tips, node = clade_node, n_pending = n_pending, n_unsampled = n_unsampled, rho = rho)

  # Export output
  return(BD_output)
}


#' @title Validate the structure of `missing_taxa_list`
#'
#' @description Check that `missing_taxa_list` has a valid `Taxon_name` item, the item defining the clades (`Target_clade` or `Sister_clade`),
#'   and that every other item provided is consistent with `Taxon_name`: same length and names identical to `Taxon_name`
#'   (same names, same order). Items that are `NULL` or absent are not checked, except the item defining the target/sister clades, which is required.
#'
#' @param missing_taxa_list List with the items `Taxon_name`, `Target_clade` (or `Sister_clade`), and optionally
#'   `include_target_clade_stem`, `Excluded_subclades` and `exclude_subclades_stem`.
#' @param clade_item Character string. `"Target_clade"` (for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::check_grafting_constraints()`][check_grafting_constraints]) or `"Sister_clade"` (for [`phylograftR::graft_missing_clades()`][graft_missing_clades]).
#' @param optional_items Character vector. Other items validated if they are provided.
#' @param allow_node_ids Logical. Whether clades defined by node IDs are accepted. If `FALSE` (default), the function stops if a clade
#'   (target or sister clade, or excluded subclade) is defined by a node ID, and points to [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips].
#'   `TRUE` is used by [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] itself.
#'
#' @noRd
#'
#' @return Invisibly, the character vector of the taxon names. The function stops if the list is not valid.
#'

validate_missing_taxa_list <- function (missing_taxa_list,
                                        clade_item = c("Target_clade", "Sister_clade"),
                                        optional_items = c("include_target_clade_stem", "Excluded_subclades", "exclude_subclades_stem"),
                                        allow_node_ids = FALSE)
{
  clade_item <- match.arg(clade_item)
  if (!is.list(missing_taxa_list)) { stop("'missing_taxa_list' must be a list") }

  ## Taxon_name: a non-empty character vector of unique and non-missing names
  taxon_name <- missing_taxa_list$Taxon_name
  if (is.null(taxon_name) || length(taxon_name) == 0) { stop("'missing_taxa_list$Taxon_name' is missing or empty") }
  if (!is.character(taxon_name) || anyNA(taxon_name) || any(taxon_name == ""))
  {
    stop("'missing_taxa_list$Taxon_name' must be a character vector without missing values or empty strings")
  }
  if (anyDuplicated(taxon_name))
  {
    stop("'missing_taxa_list$Taxon_name' must contain unique names. Duplicated: ", paste(unique(taxon_name[duplicated(taxon_name)]), collapse = ", "))
  }
  n_taxa <- length(taxon_name)

  ## Required item (clades) and optional items
  if (is.null(missing_taxa_list[[clade_item]])) { stop("'missing_taxa_list$", clade_item, "' is missing") }
  for (item in c(clade_item, optional_items))
  {
    value <- missing_taxa_list[[item]]
    if (is.null(value)) { next }

    if (length(value) != n_taxa)
    {
      stop("'missing_taxa_list$", item, "' must have the same length as 'Taxon_name' (", length(value), " vs ", n_taxa, ")")
    }
    item_names <- names(value)
    if (is.null(item_names))
    {
      stop("'missing_taxa_list$", item, "' must be named after 'Taxon_name' (same names, same order), but it has no names.\n",
           "Use e.g. names(missing_taxa_list$", item, ") <- missing_taxa_list$Taxon_name")
    }
    if (!identical(unname(item_names), taxon_name))
    {
      differ <- which(is.na(item_names) | item_names != taxon_name)
      n_show <- min(length(differ), 5)
      reason <- if (setequal(item_names, taxon_name)) "the names are the same but not in the same order" else "the names do not match"
      stop("The names of 'missing_taxa_list$", item, "' must match 'Taxon_name' (same names, same order): ", reason, ".\n",
           paste0("  - position ", differ[seq_len(n_show)], ": '", item_names[differ[seq_len(n_show)]], "' in '", item, "' vs '", taxon_name[differ[seq_len(n_show)]], "' in 'Taxon_name'", collapse = "\n"),
           if (length(differ) > n_show) paste0("\n  ... and ", length(differ) - n_show, " more") else "")
    }
  }

  ## Clades must be defined by tips: node IDs have to be converted first
  if (!allow_node_ids) { check_no_node_ids(missing_taxa_list, clade_item = clade_item) }

  invisible(taxon_name)
}


#' @title Detect a clade defined by a node ID
#'
#' @description A clade defined by a node ID is a number, or a text made only of digits (tip labels must include at least one other character).
#'
#' @param x A clade: character vector of tips, or node ID.
#'
#' @noRd
#'
#' @return Logical. `TRUE` if the clade seems defined by a node ID.
#'

is_node_id_clade <- function (x)
{
  if (is.factor(x)) { x <- as.character(x) }
  if (is.numeric(x)) { return(length(x) > 0) }
  if (is.character(x)) { return(length(x) > 0 && any(!is.na(x) & grepl("^[[:space:]]*[0-9]+[[:space:]]*$", x))) }
  FALSE
}


#' @title Stop if clades are defined by node IDs
#'
#' @description Check the clades of `missing_taxa_list` (`Target_clade` or `Sister_clade`, and `Excluded_subclades`), and stop with a message
#'   that lists some of the node IDs and points to [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] if some are found.
#'
#' @param missing_taxa_list List of missing taxa.
#' @param clade_item Character string. `"Target_clade"` or `"Sister_clade"`.
#' @param max_show Integer. Maximum number of node IDs listed in the message.
#'
#' @noRd
#'
#' @return Nothing. Stops if a clade is defined by a node ID.
#'

check_no_node_ids <- function (missing_taxa_list, clade_item, max_show = 5)
{
  found <- character(0)
  taxon_name <- missing_taxa_list$Taxon_name
  show_value <- function (x) paste0(if (length(x) > 1) "c(" else "", paste(if (is.character(x)) paste0("\"", x, "\"") else x, collapse = ", "), if (length(x) > 1) ")" else "")

  clades <- missing_taxa_list[[clade_item]]
  if (is.list(clades))
  {
    for (i in seq_along(clades))
    {
      if (is_node_id_clade(clades[[i]])) { found <- c(found, paste0("'", clade_item, "' of '", taxon_name[i], "': ", show_value(clades[[i]]))) }
    }
  }
  excluded <- missing_taxa_list$Excluded_subclades
  if (is.list(excluded))
  {
    for (i in seq_along(excluded))
    {
      subclades <- excluded[[i]]
      if (!is.list(subclades)) { next }
      for (k in seq_along(subclades))
      {
        if (is_node_id_clade(subclades[[k]]))
        {
          label <- if (!is.null(names(subclades)) && !is.na(names(subclades)[k]) && names(subclades)[k] != "") paste0("'", names(subclades)[k], "'") else paste0("#", k)
          found <- c(found, paste0("'Excluded_subclades' of '", taxon_name[i], "' (subclade ", label, "): ", show_value(subclades[[k]])))
        }
      }
    }
  }

  if (length(found) > 0)
  {
    n_show <- min(length(found), max_show)
    stop(length(found), " clade(s) of 'missing_taxa_list' are defined by node IDs, but tips are required:\n",
         paste0("  - ", found[seq_len(n_show)], collapse = "\n"),
         if (length(found) > n_show) paste0("\n  ... and ", length(found) - n_show, " more") else "",
         "\nUse phylograftR::convert_nodes_to_tips(missing_taxa_list, phy) to convert the node IDs into pairs of tips (the conversions are printed), ",
         "and pass the converted list.", call. = FALSE)
  }
  invisible(NULL)
}


#' @title Format labels of excluded subclades of a taxon
#'
#' @description Format the excluded subclades of a taxon for messages: the name of the subclade followed by its tips
#'   if the subclades are named (names inside `missing_taxa_list$Excluded_subclades[[i]]`), otherwise only the tips.
#'
#' @param excluded_subclades_i List of character vectors. The excluded subclades of one taxon (`missing_taxa_list$Excluded_subclades[[i]]`).
#' @param max_tips Integer. Maximum number of tips printed per subclade.
#'
#' @noRd
#'
#' @return A character vector with one label per excluded subclade.
#'

label_excluded_subclades <- function (excluded_subclades_i, max_tips = 5)
{
  if (length(excluded_subclades_i) == 0) { return(character(0)) }
  tips_label <- vapply(excluded_subclades_i, function (tips)
  {
    paste0(paste(utils::head(tips, max_tips), collapse = ", "), if (length(tips) > max_tips) paste0(", ... (", length(tips), " tips)") else "")
  }, character(1))
  nm <- names(excluded_subclades_i)
  has_name <- if (is.null(nm)) rep(FALSE, length(tips_label)) else (!is.na(nm) & nm != "")
  unname(ifelse(has_name, paste0("'", nm, "' (", tips_label, ")"), tips_label))
}


#' @title Resolve whether the stem of the target clade is available for each taxon
#'
#' @description Combine the taxon-level values `$include_target_clade_stem` (from `missing_taxa_list`) with the global `include_target_clade_stem` argument:
#'   taxa with `NA` (or all taxa if it is `NULL`) use the global value.
#'
#' @param include_target_clade_stem_taxa `NULL`, or a logical vector (`TRUE`, `FALSE` or `NA`) with one value per taxon.
#' @param include_target_clade_stem Logical. Global value (single `TRUE` or `FALSE`).
#' @param n_taxa Integer. Number of taxa.
#' @param taxon_names `NULL`, or the character vector of the taxon names, used to name the output.
#'
#' @noRd
#'
#' @return A logical vector without `NA`, with one value per taxon (named after the taxa if `taxon_names` is provided).
#'

resolve_include_target_clade_stem <- function (include_target_clade_stem_taxa = NULL, include_target_clade_stem = TRUE, n_taxa, taxon_names = NULL)
{
  if (!is.logical(include_target_clade_stem) || length(include_target_clade_stem) != 1 || is.na(include_target_clade_stem))
  {
    stop("the global 'include_target_clade_stem' must be a single logical value (TRUE or FALSE)")
  }
  # Default: the global argument applies to all taxa
  if (is.null(include_target_clade_stem_taxa)) { return(stats::setNames(rep(include_target_clade_stem, n_taxa), taxon_names)) }

  if (length(include_target_clade_stem_taxa) != n_taxa)
  {
    stop("'include_target_clade_stem' must have the same length as 'Taxon_name'")
  }
  # A vector with only NA of another type (e.g., NA_integer_) is accepted as logical
  if (!is.logical(include_target_clade_stem_taxa))
  {
    if (!all(is.na(include_target_clade_stem_taxa)))
    {
      stop("'include_target_clade_stem' must be a logical vector (TRUE, FALSE or NA)")
    }
    include_target_clade_stem_taxa <- as.logical(include_target_clade_stem_taxa)
  }
  # NA: use the global argument
  include_target_clade_stem_taxa[is.na(include_target_clade_stem_taxa)] <- include_target_clade_stem

  # Return the output as a logical vector named after the taxa (without names if the taxon names are not provided)
  return(stats::setNames(unname(include_target_clade_stem_taxa), taxon_names))
}


#' @title Resolve whether the stem of each excluded subclade is excluded
#'
#' @description Combine the subclade-level values `exclude_subclades_stem` (from `missing_taxa_list`) with the global `exclude_subclades_stem` argument:
#'   subclades with `NA` (or all subclades if the item is `NULL`) use the global value.
#'
#' @param exclude_subclades_stem_taxa `NULL`, a logical vector with one value per taxon (recycled over the excluded subclades of the taxon),
#'   or a list with one element per taxon. Each element is `NULL` (global value), a single value (recycled), or a logical vector with one value per excluded subclade of the taxon.
#'   `NA` values use the global value.
#' @param exclude_subclades_stem Logical. Global value (single `TRUE` or `FALSE`).
#' @param excluded_subclades `NULL`, or the list with one element per taxon of lists of tips defining the excluded subclades (`missing_taxa_list$Excluded_subclades`).
#' @param n_taxa Integer. Number of taxa.
#' @param taxon_names `NULL`, or the character vector of the taxon names, used to name the output and in error messages.
#'
#' @noRd
#'
#' @return A list with one element per taxon (named after the taxa if `taxon_names` is provided): a logical vector without `NA`, with one value per excluded subclade of the taxon (`logical(0)` if none),
#'   named after the excluded subclades if they are named.
#'

resolve_exclude_subclades_stem <- function (exclude_subclades_stem_taxa = NULL, exclude_subclades_stem = TRUE, excluded_subclades = NULL, n_taxa, taxon_names = NULL)
{
  if (!is.logical(exclude_subclades_stem) || length(exclude_subclades_stem) != 1 || is.na(exclude_subclades_stem))
  {
    stop("'exclude_subclades_stem' must be a single logical value (TRUE or FALSE)")
  }
  # Number of excluded subclades of each taxon
  n_subclades <- if (is.null(excluded_subclades)) rep(0L, n_taxa) else vapply(excluded_subclades, length, integer(1))

  # Default: the global argument applies to all subclades
  if (is.null(exclude_subclades_stem_taxa))
  {
    resolved <- lapply(seq_len(n_taxa), function (i) stats::setNames(rep(exclude_subclades_stem, n_subclades[i]), names(excluded_subclades[[i]])))
    return(stats::setNames(resolved, taxon_names))
  }

  # A vector with one value per taxon is recycled over the excluded subclades of each taxon
  if (!is.list(exclude_subclades_stem_taxa)) { exclude_subclades_stem_taxa <- as.list(exclude_subclades_stem_taxa) }
  if (length(exclude_subclades_stem_taxa) != n_taxa)
  {
    stop("'missing_taxa_list$exclude_subclades_stem' must have the same length as 'Taxon_name'")
  }
  # Label of the taxon in messages
  taxon_label <- function (i) if (is.null(taxon_names)) paste0("element ", i) else paste0("taxon '", taxon_names[i], "'")

  resolved <- vector("list", n_taxa)
  for (i in seq_len(n_taxa))
  {
    value <- exclude_subclades_stem_taxa[[i]]
    # No excluded subclade for this taxon: nothing to set
    if (n_subclades[i] == 0) { resolved[[i]] <- logical(0); next }
    # No value for this taxon: global value
    if (is.null(value)) { resolved[[i]] <- stats::setNames(rep(exclude_subclades_stem, n_subclades[i]), names(excluded_subclades[[i]])); next }

    # A vector with only NA of another type (e.g., NA_integer_) is accepted as logical
    if (!is.logical(value))
    {
      if (!all(is.na(value))) { stop("'missing_taxa_list$exclude_subclades_stem' must contain logical values (TRUE, FALSE or NA), see ", taxon_label(i)) }
      value <- as.logical(value)
    }
    if (!(length(value) %in% c(1L, n_subclades[i])))
    {
      stop("'missing_taxa_list$exclude_subclades_stem': ", taxon_label(i), " has ", length(value), " values, but this taxon has ", n_subclades[i],
           " excluded subclades.\nProvide one value per excluded subclade, or a single value for all of them.")
    }
    # If both the values and the excluded subclades are named, the names must match (same names, same order)
    subclade_names <- names(excluded_subclades[[i]])
    if (length(value) > 1 && !is.null(names(value)) && !is.null(subclade_names) && !identical(names(value), subclade_names))
    {
      stop("'missing_taxa_list$exclude_subclades_stem': the names of the values of ", taxon_label(i), " (", paste(names(value), collapse = ", "),
           ") must match the names of its excluded subclades (", paste(subclade_names, collapse = ", "), "), in the same order.")
    }
    value <- rep_len(unname(value), n_subclades[i])
    value[is.na(value)] <- exclude_subclades_stem
    resolved[[i]] <- stats::setNames(value, subclade_names)
  }
  return(stats::setNames(resolved, taxon_names))
}


#' @title Get the values of exclude_subclades_stem of one taxon
#'
#' @description Extract the stem exclusion values of the excluded subclades of one taxon, from the list resolved by `resolve_exclude_subclades_stem()`,
#'   or from a single logical value recycled over the excluded subclades.
#'
#' @param exclude_subclades_stem List with one element per taxon, or a single logical value.
#' @param taxon_index Integer. Index of the taxon.
#' @param n_subclades Integer. Number of excluded subclades of the taxon.
#'
#' @noRd
#'
#' @return A logical vector of length `n_subclades`.
#'

get_exclude_stem_flags <- function (exclude_subclades_stem, taxon_index, n_subclades)
{
  value <- if (is.list(exclude_subclades_stem)) exclude_subclades_stem[[taxon_index]] else exclude_subclades_stem
  if (n_subclades == 0) { return(logical(0)) }
  return(rep_len(value, n_subclades))
}


#' @title Identify the branches of the excluded subclades
#'
#' @description Identify the tipward node IDs of the branches to remove from the pool of branches available for grafting a taxon:
#'   * The branches within each excluded subclade (all the descendant nodes of its MRCA node).
#'   * The stem region of each excluded subclade (see `get_stem_region_nodes()`), if its stem is excluded.
#'     A taxon grafted on the stem region would be sister to the excluded subclade, or sister to other missing taxa grafted on its stem.
#'
#'   For a single-tip subclade, the terminal branch is the only branch (it is its stem): it is excluded only if its stem is excluded,
#'   otherwise the subclade excludes nothing (a warning would be sent by the grafting function).
#'
#' @param phy Phylogeny of class `"phylo"`, possibly including taxa already grafted.
#' @param excluded_subclades_tips_list List of character vectors of tips, defining the excluded subclades of one taxon. Can be `NULL`.
#' @param exclude_subclades_stem Logical vector with one value per excluded subclade, or a single value recycled. Default is `TRUE`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not). Default is `character(0)`.
#'
#' @importFrom ape Ntip
#'
#' @noRd
#'
#' @return Integer vector of tipward node IDs (without duplicates). Empty if nothing is excluded.
#'

get_excluded_branch_nodes <- function (phy, excluded_subclades_tips_list, exclude_subclades_stem = TRUE, missing_tip_labels = character(0))
{
  n_subclades <- length(excluded_subclades_tips_list)
  if (n_subclades == 0) { return(integer(0)) }
  stem_flags <- rep_len(exclude_subclades_stem, n_subclades)

  excluded_nodes <- integer(0)
  for (j in seq_len(n_subclades))
  {
    tips <- excluded_subclades_tips_list[[j]]
    if (length(tips) == 0) { next }
    node <- get_clade_node(phy, tips)
    # A single tip excludes its terminal branch (its stem) only if its stem is excluded
    if (node <= ape::Ntip(phy) && !stem_flags[j]) { next }
    # Branches within the subclade
    excluded_nodes <- c(excluded_nodes, get_descendant_nodes(phy, node))
    # Stem region of the subclade
    if (stem_flags[j])
    {
      excluded_nodes <- c(excluded_nodes, get_stem_region_nodes(phy, node, missing_tip_labels = missing_tip_labels)$nodes)
    }
  }
  return(unique(excluded_nodes))
}


#' @title Group missing taxa by clade and order the clades by nestedness
#'
#' @description Group missing taxa by the node of their MRCA and by the availability of the stem of their target clade,
#'   then order groups by size so nested ones come first
#'   (a clade always has strictly more descendants than any clade nested within it,
#'   thus sorting on clade size is enough to list nested clades first).
#'   The rationale is to graft taxa in nested clades first, so they are accounted for when grafting missing taxa in higher clades
#'   comprising those nested clades.
#'
#'   Taxa with the same target clade but different `include_target_clade_stem` values form different groups, because the time window available
#'   for their divergence times differs (stem age or crown age). The groups of a same clade are consecutive, with `include_target_clade_stem = FALSE` first
#'   (taxa restricted to the crown are grafted first, so they are present when drawing the taxa that can also use the stem).
#'
#'   Use tips as inputs to define the associated clades:
#'     * One tip defines a terminal branch.
#'     * Two or more tips define a clade based on their MRCA.
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param clade_tips List of character vectors defining the clade of each taxon.
#' @param include_target_clade_stem Logical vector with one value per taxon (no `NA`), or `NULL` (all `TRUE`).
#'   Whether the stem of the target clade is available for each taxon. It is always `TRUE` for single-tip clades (their only branch is their stem).
#'
#' @noRd
#'
#' @return A list with 6 elements.
#'
#'  The first two are taxon-level elements (one value per taxon):
#'  * `group_ID` Integer vector. ID of the group (i.e., unique clade and `include_target_clade_stem` value) the taxon belongs to.
#'  * `taxa_order` Integer vector. Indices of the taxa ordered by group, so the taxa of nested clades come first.
#'     Used to graft taxa one by one while processing nested clades first.
#'
#'  The last four elements are group-level elements (one value per group):
#'  * `group_clade_tips` List with one element per group. Character vector of the tips used to identify the associated clade.
#'  * `group_clade_ID` Integer vector. ID of the unique clade of the group. Groups of the same clade share it.
#'  * `group_include_target_clade_stem` Logical vector. Whether the stem of the target clade is available for the taxa of the group.
#'  * `group_order` Integer vector with one value per group. Group IDs ranked by clade size, so nested clades come first.
#'

group_and_order_clades <- function (phy, clade_tips, include_target_clade_stem = NULL)
{
  n_taxa <- length(clade_tips)
  # Extract (MRCA) nodes associated to clades
  nodes <- vapply(clade_tips, function (tips) get_clade_node(phy, tips), numeric(1))
  # ID of the unique clades, based on similarity of (MRCA) nodes
  clade_ID <- match(nodes, unique(nodes))

  # Availability of the stem. Single-tip clades always use their stem
  if (is.null(include_target_clade_stem)) { include_target_clade_stem <- rep(TRUE, n_taxa) }
  include_target_clade_stem <- include_target_clade_stem | (nodes <= ape::Ntip(phy))

  # Group taxa by clade and stem availability
  key <- paste(clade_ID, include_target_clade_stem, sep = "_")
  group <- match(key, unique(key))
  # Record first appearance of each group in the list
  first <- match(seq_along(unique(key)), group)
  group_clade_ID <- clade_ID[first]
  group_include_target_clade_stem <- include_target_clade_stem[first]

  # Compute number of descendants of each clade, and record it per group
  n_desc_clade <- vapply(unique(nodes), function (n) length(get_descendant_nodes(phy, n)), numeric(1))
  n_desc <- n_desc_clade[group_clade_ID]

  # Order of the groups: smallest clades first. Groups of the same clade stay together, with 'include_target_clade_stem = FALSE' first
  group_order <- order(n_desc, group_clade_ID, group_include_target_clade_stem)
  # Rank of each group: 1 for the first group to graft
  rank_of_group <- match(seq_along(group_order), group_order)

  # Record information: taxon-level elements (one value per taxon), then group-level elements (one value per group)
  groups_per_taxa <- list(group_ID = group, # Group ID of each taxon
                          taxa_order = order(rank_of_group[group]), # Order of the taxa by group, nested clades first.
                          group_clade_tips = clade_tips[first], # Tips defining the clade of each group
                          group_clade_ID = group_clade_ID, # ID of the unique clade of each group
                          group_include_target_clade_stem = group_include_target_clade_stem, # Stem availability of each group
                          group_order = group_order) # Order of the groups to use in the loop

  # Export output
  return(groups_per_taxa)
}


#' @title Extract the new divergence times returned by `TreeSim::corsim()`
#'
#' @description Extract the new divergence times returned by `TreeSim::corsim()`.
#'   `TreeSim::corsim()` function returns the input times `x` plus `missing` new times associated to the missing taxa, sorted in decreasing order.
#'   This function ensure a safe comparison between the input times and output times to extract the new divergence times associated with missing taxa,
#'   by removing ONE instance of each input time from the output.
#'   Contrary to `setdiff()`, this is robust to a new time being numerically identical to an existing one, and to duplicated input times.
#'
#' @param all_times Numeric vector. Output of `TreeSim::corsim()`.
#' @param old_times Numeric vector. Input times given to `TreeSim::corsim()` as `x`.
#'
#' @noRd
#'
#' @return Numeric vector of the new times.
#'

extract_new_divergence_times <- function (all_times, # Output of TreeSim::corsim()
                                          old_times) # Input times 'x' given to TreeSim::corsim()
{
  all_times <- unname(all_times)
  old_times <- unname(old_times)
  remaining <- all_times

  # Remove one instance of each old time
  for (t in old_times)
  {
    idx <- match(t, remaining)
    remaining <- remaining[-idx]
  }

  new_times <- sort(remaining, decreasing = TRUE)

  return(new_times)
}


#' @title Set a seed for the duration of the calling function only
#'
#' @description Set a seed used internally by the function, without altering the global `.Random.seed` of the R session.
#'   Restore the global `.Random.seed` when the calling function exits or yields an error.
#'   Does nothing if seed is `NULL`.
#'
#' @param seed Integer. Seed. Nothing is done if `NULL`.
#' @param envir Environment of the calling function. Default is `parent.frame()`.
#'
#' @noRd
#'
#' @return Invisible `NULL`.
#'

set_local_seed <- function (seed, envir = parent.frame())
{
  if (is.null(seed)) return(invisible(NULL))

  global_env <- globalenv()
  had_seed <- exists(".Random.seed", envir = global_env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = global_env, inherits = FALSE) else NULL

  restore_seed <- function ()
  {
    if (is.null(old_seed))
    {
      if (exists(".Random.seed", envir = global_env, inherits = FALSE)) rm(".Random.seed", envir = global_env)
    } else {
      assign(".Random.seed", old_seed, envir = global_env)
    }
  }
  # Register the restoration on exit of the calling function
  do.call(base::on.exit, list(as.call(list(restore_seed)), TRUE), envir = envir)

  set.seed(seed = seed)
  invisible(NULL)
}


#' @title Check the position of the excluded subclades relative to their target clade
#'
#' @description Classify, for each taxon, the position of each of its excluded subclades relative to its target clade in the backbone phylogeny.
#'   Excluded subclades are expected to be nested within the target clade. Because clades of a tree are either nested or disjoint,
#'   there are four possible relationships (single-tip subclades with a stem that is not excluded exclude nothing and are ignored):
#'   * `"nested"`: The excluded subclade is within the target clade (expected).
#'   * `"identical"`: The excluded subclade is the target clade itself. The whole crown is excluded, so no branch (or only the stem) is left.
#'     This is flagged as an ERROR as `graft_missing_clades()` should be used to graft specifically on the stem as a sister clade.
#'   * `"disjoint"`: The excluded subclade does not overlap with the target clade. The exclusion has no effect.
#'     This is flagged as a WARNING as it is an erroneous setting from the user.
#'   * `"contains_target"`: The excluded subclade contains the target clade. All branches of the target clade (stem included) are excluded, so no branch is available.
#'     This is flagged as an ERROR as grafting will be impossible.
#'
#'   Only the problematic relationships are returned.
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param missing_taxa_list List with `Taxon_name`, `Target_clade` and `Excluded_subclades`, and optionally `exclude_subclades_stem`
#'   (list resolved by `resolve_exclude_subclades_stem()`; all stems are assumed excluded if absent).
#'
#' @importFrom ape getMRCA
#'
#' @noRd
#'
#' @return A data.frame with one row per problematic excluded subclade (0 rows if none), and the columns:
#'  * `Taxon_name` Character string. Name of the missing taxon.
#'  * `Target_clade` Character string. Tips defining the target clade, separated by commas.
#'  * `Excluded_subclade` Character string. Tips defining the excluded subclade, separated by commas (preceded by the name of the subclade if it is named).
#'  * `relation` Character string. One of `"identical"`, `"disjoint"` or `"contains_target"`.
#'

check_excluded_subclades_position <- function (phy, missing_taxa_list)
{
  n_missing <- length(missing_taxa_list$Taxon_name)

  # Cache the MRCA nodes and descendant nodes of clades, as many taxa share the same target clades / excluded subclades
  # (same approach as in find_taxa_without_allowed_branch() used to compute MRCA nodes and descending nodes once per similar target clades / excluded subclades)
  cache <- new.env()
  get_node <- function (tips)
  {
    key <- paste0("mrca_", paste(sort(tips), collapse = "|"))
    if (is.null(cache[[key]])) { cache[[key]] <- get_clade_node(phy, tips) }
    cache[[key]]
  }
  get_desc <- function (node)
  {
    key <- paste0("desc_", node)
    if (is.null(cache[[key]])) { cache[[key]] <- get_descendant_nodes(phy, node) }
    cache[[key]]
  }

  # Initiate the output
  out <- list()

  # Loop per missing taxa
  for (i in seq_len(n_missing))
  {
    all_excluded_i <- missing_taxa_list$Excluded_subclades[[i]]
    if (length(all_excluded_i) == 0) next
    # Stem exclusion of each excluded subclade (all excluded by default)
    stem_flags_i <- get_exclude_stem_flags(if (is.null(missing_taxa_list$exclude_subclades_stem)) TRUE else missing_taxa_list$exclude_subclades_stem,
                                           i, length(all_excluded_i))

    target_clade_tips_i <- missing_taxa_list$Target_clade[[i]]
    target_node_i <- get_node(target_clade_tips_i)

    for (j in seq_along(all_excluded_i))
    {
      tips <- all_excluded_i[[j]]
      # Empty subclades exclude nothing. Neither do single tips with a stem that is not excluded
      if (length(tips) == 0 || (length(tips) == 1 && !stem_flags_i[j])) next
      excluded_node <- get_node(tips)

      # Compare the position of the excluded subclade to the target clade
      relation <- if (excluded_node == target_node_i) {
        "identical" # The excluded subclade is the target clade. ERROR: no branch remains, or only the stem branch if it is available and not excluded ; advise to use graft_missing_clades() with the target clade as Sister_clade instead.
      } else if (excluded_node %in% get_desc(target_node_i)) {
        "nested" # Expected case: the excluded subclade is nested in the target clade
      } else if (target_node_i %in% get_desc(excluded_node)) {
        "contains_target" # The target clade is nested in the excluded subclade. ERROR: No valid branches remain.
      } else {
        "disjoint" # The two clades are polyphyletic and share no descending tips. WARNING: The exclusion is useless.
      }

      # Format output df wiht only the problematic cases
      if (relation != "nested")
      {
        out[[length(out) + 1]] <- data.frame(
          Taxon_name = missing_taxa_list$Taxon_name[i],
          Target_clade = paste(target_clade_tips_i, collapse = ", "),
          Excluded_subclade = label_excluded_subclades(all_excluded_i[j], max_tips = Inf),
          relation = relation, stringsAsFactors = FALSE)
      }
    }
  }

  # Return a df, with 0 rows if no problem was found
  if (length(out) == 0)
  {
    return(data.frame(Taxon_name = character(0), Target_clade = character(0),
                      Excluded_subclade = character(0), relation = character(0), stringsAsFactors = FALSE))
  }
  return(do.call(rbind, out))
}


#' @title Find taxa with no branch available for grafting according to constraints
#'
#' @description Identify the taxa with no branch available for grafting, given their target clade and excluded subclades.
#'   This function is called as an early check, before any grafting, to save time instead of discovering potential problems in the middle of a long grafting run.
#'
#'   # How it works, for each taxon:
#'   1. Start from the pool of branches of its target clade.
#'      A branch is identified by its tipward node (the node it leads to), so the pool is made of
#'      all the descendant nodes of the MRCA node of the clade, plus the MRCA node itself if the stem is allowed for grafting (`include_target_clade_stem = TRUE` for this taxon).
#'   2. Remove from the pool the branches inside each excluded subclade (= all descendant nodes of the subclade's MRCA node),
#'      and the stem branch of each excluded subclade whose stem is excluded (`exclude_subclades_stem = TRUE`).
#'      A single-tip subclade has no branch but its stem: it excludes nothing if its stem is not excluded.
#'   3. If the pool is empty, the taxon is flagged.
#'
#' It mirrors the rules of the grafting functions:
#'   - the stem of the target clade is available if `include_target_clade_stem = TRUE` (never for the root, which has no stem)
#'   - the stem of a single-tip clade is always available (it is its only branch), ignoring `include_target_clade_stem`
#'   - branches inside excluded subclades are unavailable, and so is their own stem if it is excluded (`exclude_subclades_stem = TRUE`)
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param missing_taxa_list List with `Taxon_name`, `Target_clade` and `Excluded_subclades`.
#' @param include_target_clade_stem Logical. Single value used for all taxa, or one value per taxon (resolved from `include_target_clade_stem`, see `resolve_include_target_clade_stem()`). Default is `TRUE`.
#' @param use_excluded_subclades Logical. Default is `TRUE`.
#' @param exclude_subclades_stem Logical, or list with one element per taxon of logical vectors (one value per excluded subclade), as resolved by
#'   `resolve_exclude_subclades_stem()`. Whether the stem of the excluded subclades is excluded. Default is `TRUE`.
#'
#' @importFrom ape getMRCA
#'
#' @noRd
#'
#' @return A data.frame with two columns summarizing the flagged taxa. It has 0 rows if no taxa were flagged.
#'   * `Taxon_name` Character string. Name of the missing taxon to graft.
#'   * `Target_clade` Character string. Tips used to describe the target clade associated with the missing taxon; separated by commas if multiple.
#'

find_taxa_without_allowed_branch <- function (phy,
                                              missing_taxa_list,
                                              include_target_clade_stem = TRUE,
                                              use_excluded_subclades = TRUE,
                                              exclude_subclades_stem = TRUE)
{
  n_missing <- length(missing_taxa_list$Taxon_name)
  root_node <- setdiff(phy$edge[, 1], phy$edge[, 2])
  # One value per taxon
  include_target_clade_stem_taxa <- rep_len(include_target_clade_stem, n_missing)

  # Use temporary internal environment as cache to store associated MRCA nodes and descendant nodes information.
  # Rationale: Many taxa share the same target/excluded clades. Instead of recomputing descendant nodes for each target/excluded clades,
  # the function 'cached' is used to access the environment as a library, and check if its results have already been computed / stored.
  # For each clade, a key is generated by concatenating the tips defining a clade (target clade or excluded subclade) or MRCA nodes, and its associated MRCA nodes or descending nodes (respectively) are stored under each key once.
  # If a key is already present, the results can be extracted directly instead of computed again.
  cache <- new.env()
  cached <- function (key, expr) # Return cache[[key]] if stored; otherwise evaluate 'expr', store it in the internal environment and return it
  {
    if (is.null(cache[[key]])) { cache[[key]] <- expr }
    cache[[key]]
  }

  # Initiate vector of flagged
  flagged <- logical(n_missing)

  ## Loop per taxa
  for (i in seq_len(n_missing))
  {
    # Extract the tips defining the target clade
    target_clade_tips_i <- missing_taxa_list$Target_clade[[i]]
    # Retrieve associated MRCA node (from cache env. or compute for the first time)
    MRCA_node_i <- cached(paste0("mrca_", paste(sort(target_clade_tips_i), collapse = "|")), get_clade_node(phy, target_clade_tips_i))

    # Detect case where the target clade is a unique tip
    is_single_tip <- length(target_clade_tips_i) == 1

    # Is the stem available for grafting?
    # Only if requested or the 'clade' is a tip; but the clade cannot be the whole tree.
    stem_available <- (include_target_clade_stem_taxa[i] || is_single_tip) && (MRCA_node_i != root_node)

    # Identify the tipward nodes associated with the valid branches (from cache env. or compute for the first time)
    allowed_nodes <- c(cached(paste0("desc_", MRCA_node_i), get_descendant_nodes(phy, MRCA_node_i)),
                       if (stem_available) MRCA_node_i)

    # Remove branches from excluded subclades
    if (use_excluded_subclades)
    {
      excluded_i <- missing_taxa_list$Excluded_subclades[[i]]
      stem_flags_i <- get_exclude_stem_flags(exclude_subclades_stem, i, length(excluded_i))
      for (j in seq_along(excluded_i))
      {
        tips <- excluded_i[[j]]
        # Empty subclades exclude nothing. Neither do single tips with a stem that is not excluded (their terminal branch is their only branch)
        if (length(tips) == 0 || (length(tips) == 1 && !stem_flags_i[j])) { next }
        # Identify the tipward nodes associated with the INvalid branches (from cache env. or compute for the first time)
        excluded_node <- cached(paste0("mrca_", paste(sort(tips), collapse = "|")), get_clade_node(phy, tips))
        # Remove them from available nodes: branches within the subclade, and its stem if excluded
        allowed_nodes <- setdiff(allowed_nodes, cached(paste0("desc_", excluded_node), get_descendant_nodes(phy, excluded_node)))
        if (stem_flags_i[j]) { allowed_nodes <- setdiff(allowed_nodes, excluded_node) }
      }
    }
    # The taxa if flagged only if no valid branches (tipward nodes) remains in the pool
    flagged[i] <- length(allowed_nodes) == 0
  }

  # Build data.frame summarizing all taxa that did not pass the check
  # If no taxon is flagged, return an empty df with 0 rows.
  flagged_taxa_df <- data.frame(
    Taxon_name = missing_taxa_list$Taxon_name[flagged],
    Target_clade = vapply(missing_taxa_list$Target_clade[flagged], paste, character(1), collapse = ", "), stringsAsFactors = FALSE)

  # Return df
  return(flagged_taxa_df)
}


#' @title Check the path used for temporary saves
#'
#' @description Check that `temp_file_path` is valid, and stops with an informative error if not.
#'   Nothing is checked if `temp_file_path` is `NULL`, as the default path is in the working directory.
#'   The path must be a single character string that:
#'   * Ends with `.rds`.
#'   * Is located in an existing directory.
#'   * Is located in a writable directory.
#'
#' @param temp_file_path Character string or `NULL`. Path of the `.rds` file used for temporary saves.
#'
#' @noRd
#'
#' @return Invisible `TRUE`. Stops if the path is not valid.
#'

check_temp_file_path <- function (temp_file_path)
{
  if (is.null(temp_file_path)) return(invisible(TRUE))

  # Must be a single path
  if (!is.character(temp_file_path) || length(temp_file_path) != 1 || is.na(temp_file_path))
  {
    stop("'temp_file_path' must be a single character string")
  }
  # Must end with ".rds"
  if (!grepl(pattern = "\\.rds$", x = temp_file_path))
  {
    stop("'temp_file_path' must end with '.rds'")
  }
  # Directory must exist
  temp_dir <- dirname(temp_file_path)
  if (!dir.exists(temp_dir))
  {
    stop("The directory of 'temp_file_path' does not exist: ", temp_dir)
  }
  # Directory must be writable
  if (file.access(names = temp_dir, mode = 2) != 0)
  {
    stop("The directory of 'temp_file_path' is not writable: ", temp_dir)
  }

  invisible(TRUE)
}


#' @title Build the richness information used for the sampling fractions
#'
#' @description Validate the `total_richness_per_clades` input, and derive from it, for every node of the backbone,
#'   the number of unsampled species, i.e., species of the clade that are neither in the backbone nor in the list of taxa to graft.
#'   The output is used:
#'    * (1) to compute the sampling fraction `rho` of the birth-death models fitted during the grafting process;
#'    * (2) to build the `sampling_fraction_df` table added to the grafted phylogeny that describes the distribution of species among
#'      `nb_backbone_species`, `nb_missing_species`, and `nb_unsampled_species`, for each clade listed in `total_richness_per_clades`.
#'
#' @details
#' ### Counts
#' For a clade k with a total richness provided (not `NA`):
#'   `nb_unsampled = total_richness - nb_backbone - nb_missing`.
#'   `nb_backbone` is the number of tips of the clade in the backbone.
#'   `nb_missing` is the number of taxa to graft that belong to the clade:
#'   * `mode = "target"` (`graft_missing_taxa()`): taxa whose target clade is the clade itself, or is nested within it.
#'   * `mode = "sister"` (`graft_missing_clades()`): taxa whose sister clade is strictly nested within the clade.
#'     The taxa grafted on the stem branch of the clade itself are sister to the clade: they are not part of it.
#'
#' For a clade without a total richness (`NA`, or not listed, e.g., an enclosing clade), `nb_unsampled` is the sum of `nb_unsampled`
#' over the largest clades with a total richness that are nested within it.
#' The rationale is that even if `nb_unsampled` cannot be derived from its `total_richness` (as it is not provided for this clade),
#' it has to at least account for the `nb_unsampled` species that have being identified in its own nested clades.
#' It is `0` if there is none.
#'
#' ### Checks
#' The function stops if:
#'  * The structure of `total_richness_per_clades` is not valid.
#'  * Tips defining a clade are missing from the backbone, or a clade is listed twice.
#'  * A target clade (or sister clade) used in `missing_taxa_list` is not listed in `total_richness_per_clades`.
#'  * A total richness is lower than the number of backbone species plus missing species it contains (`nb_unsampled < 0`).
#'  * A clade has fewer unsampled species than the sum of those of the clades with a total richness nested within it.
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param missing_taxa_names Character vector. Names of the taxa to graft.
#' @param missing_clades List of character vectors, same length as `missing_taxa_names`. Tips defining the target clade
#'   (or the sister clade) of each taxon.
#' @param mode Character string. `"target"` or `"sister"`. See details.
#' @param total_richness_per_clades `NULL`, or a list with two elements: the clades (`Target_clade` or `Sister_clade`, depending on `mode`),
#'   a list of character vectors of tips defining each clade, and `total_richness`, a vector of postive integers of the same length.
#'   `NA` is allowed when no information for a given clade is available.
#'
#' @importFrom ape node.depth
#'
#' @importFrom stats setNames
#'
#' @noRd
#'
#' @return A list with 6 elements:
#'  * `phy_backbone` The backbone phylogeny.
#'  * `mode` Character string. The mode used (`"target"` or `"sister"`).
#'  * `taxa_names` Vector of character strings. Names of the taxa to graft.
#'  * `taxon_sets` List (one element per taxon) of the node IDs of the backbone whose clade is considered to contain the taxon.
#'  * `unsampled_by_node` Numeric vector (one value per node of the backbone, tips included). Number of unsampled species per node.
#'  * `table` `NULL` if `total_richness_per_clades` is `NULL`. Otherwise, a data frame with one row per clade listed in `total_richness_per_clades`.
#'

build_richness_info <- function (phy, missing_taxa_names, missing_clades, mode = c("target", "sister"), total_richness_per_clades = NULL)
{
  mode <- match.arg(mode)
  clade_element <- if (mode == "target") "Target_clade" else "Sister_clade"
  n_nodes <- length(phy$tip.label) + phy$Nnode

  ## Identify parent of every node (NA for the root) and ancestors of a node
  parent <- rep(NA_integer_, n_nodes)
  parent[phy$edge[, 2]] <- phy$edge[, 1]
  get_ancestors <- function (node)
  {
    out <- integer(0)
    a <- parent[node]
    while (!is.na(a)) { out <- c(out, a); a <- parent[a] }
    out
  }

  ## Identify nodes of the clades of the missing taxa, and the backbone nodes whose clade contains each taxon
  missing_nodes <- vapply(missing_clades, function (tips) get_clade_node(phy, tips), numeric(1))
  taxon_sets <- lapply(missing_nodes, function (nd)
  {
    anc <- get_ancestors(nd)
    if (mode == "target") c(nd, anc) else anc
  })

  unsampled_by_node <- numeric(n_nodes)
  richness_table <- NULL

  if (!is.null(total_richness_per_clades))
  {
    ## ---- Check the structure of 'total_richness_per_clades' ---- ##
    if (!is.list(total_richness_per_clades) || is.null(total_richness_per_clades[[clade_element]]) || is.null(total_richness_per_clades$total_richness))
    {
      stop("'total_richness_per_clades' must be a list with two elements: '", clade_element, "' (list of character vectors of tips defining the clades) ",
           "and 'total_richness' (vector of the total number of species of each clade)")
    }
    clades <- total_richness_per_clades[[clade_element]]
    values <- total_richness_per_clades$total_richness
    node_clades <- which(vapply(clades, is_node_id_clade, logical(1)))
    if (is.list(clades) && length(node_clades) > 0)
    {
      n_show <- min(length(node_clades), 5)
      stop(length(node_clades), " clade(s) of 'total_richness_per_clades$", clade_element, "' are defined by node IDs, but tips are required:\n",
           paste0("  - clade ", node_clades[seq_len(n_show)], ": ", vapply(clades[node_clades[seq_len(n_show)]], paste, character(1), collapse = ", "), collapse = "\n"),
           if (length(node_clades) > n_show) paste0("\n  ... and ", length(node_clades) - n_show, " more") else "",
           "\n'total_richness_per_clades' must be converted AT THE SAME TIME as its associated 'missing_taxa_list', so that the tips defining a clade match between both objects:\n",
           "  converted <- phylograftR::convert_nodes_to_tips(missing_taxa_list, phy, total_richness_per_clades = total_richness_per_clades)\n",
           "and pass converted$missing_taxa_list and converted$total_richness_per_clades.\n",
           "If 'missing_taxa_list' was already converted, use it in this call as a template: the tips it already uses for the clades are matched to the node IDs of 'total_richness_per_clades'.", call. = FALSE)
    }
    if (!is.list(clades) || length(clades) == 0 || !all(vapply(clades, function (x) is.character(x) && length(x) >= 1, logical(1))))
    {
      stop("'total_richness_per_clades$", clade_element, "' must be a non-empty list of character vectors of tips")
    }
    if (length(values) != length(clades))
    {
      stop("'total_richness_per_clades$total_richness' must have the same length as 'total_richness_per_clades$", clade_element, "'")
    }
    # A vector with only NA is logical: accept it
    if (is.logical(values) && all(is.na(values))) values <- as.numeric(values)
    if (!is.numeric(values) || any(!is.na(values) & (values < 0 | values != round(values))))
    {
      stop("'total_richness_per_clades$total_richness' must be a numeric vector of non-negative whole numbers (NA allowed)")
    }

    ## ---- Locate the clades ---- ##
    all_tips <- unique(unlist(clades))
    if (!all(all_tips %in% phy$tip.label))
    {
      stop("Some tip labels provided in 'total_richness_per_clades$", clade_element, "' are missing from the backbone phylogeny\n",
           "Missing tips: ", paste(setdiff(all_tips, phy$tip.label), collapse = ", "))
    }
    clade_nodes <- vapply(clades, function (tips) get_clade_node(phy, tips), numeric(1))
    clade_labels <- vapply(clades, function (tips) paste(tips, collapse = ", "), character(1))

    if (anyDuplicated(clade_nodes))
    {
      dup <- unique(clade_nodes[duplicated(clade_nodes)])
      stop("Some clades are listed more than once in 'total_richness_per_clades$", clade_element, "' (clades are matched by their MRCA node):\n",
           paste0("  - ", vapply(dup, function (nd) paste(clade_labels[clade_nodes == nd], collapse = "  =  "), character(1)), collapse = "\n"))
    }

    ## Every clade used in 'missing_taxa_list' must be listed
    used_nodes <- unique(missing_nodes)
    not_covered <- setdiff(used_nodes, clade_nodes)
    if (length(not_covered) > 0)
    {
      tips_not_covered <- vapply(not_covered, function (nd) paste(missing_clades[[match(nd, missing_nodes)]], collapse = ", "), character(1))
      n_show <- min(length(tips_not_covered), 10)
      stop(length(not_covered), " ", if (mode == "target") "target" else "sister", " clade(s) of 'missing_taxa_list' are not listed in 'total_richness_per_clades$", clade_element, "'. ",
           "When 'total_richness_per_clades' is provided, it must cover all the unique clades (use NA in 'total_richness' if the richness is unknown).\n",
           paste0("  - ", tips_not_covered[seq_len(n_show)], collapse = "\n"),
           if (length(not_covered) > n_show) paste0("\n  ... and ", length(not_covered) - n_show, " more") else "")
    }

    ## ---- Counts ---- ##
    n_backbone_all <- ape::node.depth(phy, method = 1) # Number of tips below each node (1 for tips)
    nb_backbone <- n_backbone_all[clade_nodes]
    nb_missing <- vapply(clade_nodes, function (k) sum(vapply(taxon_sets, function (s) k %in% s, logical(1))), numeric(1))

    provided <- !is.na(values)
    nb_unsampled <- rep(NA_real_, length(clades))
    nb_unsampled[provided] <- values[provided] - nb_backbone[provided] - nb_missing[provided]

    ## A total richness cannot be lower than the species known in the clade
    too_low <- which(provided & nb_unsampled < 0)
    if (length(too_low) > 0)
    {
      n_show <- min(length(too_low), 10)
      stop(length(too_low), " clade(s) have a 'total_richness' lower than the number of backbone species plus missing species they contain.\n",
           paste0("  - ", clade_labels[too_low[seq_len(n_show)]], ": total_richness = ", values[too_low[seq_len(n_show)]],
                  ", backbone species = ", nb_backbone[too_low[seq_len(n_show)]], ", missing species = ", nb_missing[too_low[seq_len(n_show)]], collapse = "\n"),
           if (length(too_low) > n_show) paste0("\n  ... and ", length(too_low) - n_show, " more") else "")
    }

    ## ---- Nesting between clades with a total richness ---- ##
    prov_idx <- which(provided)
    prov_nodes <- clade_nodes[prov_idx]
    children_sum <- stats::setNames(numeric(length(prov_idx)), prov_nodes)
    for (j in seq_along(prov_idx))
    {
      # Climb rootward until reaching a clade with a total richness: it is the parent of this clade
      a <- parent[prov_nodes[j]]
      while (!is.na(a) && !(a %in% prov_nodes))
      {
        # Nodes without total richness (enclosing clades) get the unsampled species of the largest clades with a total richness below them
        unsampled_by_node[a] <- unsampled_by_node[a] + nb_unsampled[prov_idx[j]]
        a <- parent[a]
      }
      if (!is.na(a)) children_sum[as.character(a)] <- children_sum[as.character(a)] + nb_unsampled[prov_idx[j]]
    }
    # Unsampled species of a clade are never lower than the sum of those of the clades with a total richness nested within it
    inconsistent <- which(children_sum > nb_unsampled[prov_idx])
    if (length(inconsistent) > 0)
    {
      n_show <- min(length(inconsistent), 10)
      stop(length(inconsistent), " clade(s) have a 'total_richness' too low to contain the clades with a 'total_richness' nested within them.\n",
           paste0("  - ", clade_labels[prov_idx[inconsistent[seq_len(n_show)]]], ": unsampled species = ", nb_unsampled[prov_idx[inconsistent[seq_len(n_show)]]],
                  ", but nested clades have ", children_sum[inconsistent[seq_len(n_show)]], " unsampled species in total", collapse = "\n"),
           if (length(inconsistent) > n_show) paste0("\n  ... and ", length(inconsistent) - n_show, " more") else "")
    }

    # Nodes with a total richness use their own value
    unsampled_by_node[prov_nodes] <- nb_unsampled[prov_idx]
    # Clades without a total richness: derived from nested clades
    nb_unsampled[!provided] <- unsampled_by_node[clade_nodes[!provided]]

    richness_table <- data.frame(
      Clade_tips = clade_labels,
      total_richness = nb_backbone + nb_missing + nb_unsampled,
      total_richness_provided = provided,
      nb_backbone_species = nb_backbone,
      nb_missing_species = nb_missing,
      nb_unsampled_species = nb_unsampled,
      sampling_fraction = (nb_backbone + nb_missing) / (nb_backbone + nb_missing + nb_unsampled),
      stringsAsFactors = FALSE)
    rownames(richness_table) <- NULL
  }

  list(phy_backbone = phy, mode = mode, taxa_names = missing_taxa_names, taxon_sets = taxon_sets,
       unsampled_by_node = unsampled_by_node, table = richness_table)
}


#' @title Label the internal nodes of the backbone to follow them during grafting
#'
#' @description Node IDs are renumbered each time a tip is grafted, but node labels are kept by `phytools::bind.tip()` and `ape::bind.tree()`.
#'   The internal nodes of the backbone are therefore labelled `"bb<node ID>"`, so they can be identified in the grafted tree
#'   (and the original node labels restored with `restore_node_labels()`).
#'   The labelling is done in all cases, so the labels of the final tree do not depend on the options used to graft (e.g., with or without a simmap).
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#'
#' @noRd
#'
#' @return A list with 3 elements:
#'  * `phy` The backbone phylogeny with its internal nodes labelled `"bb<node ID>"`.
#'  * `original_node_labels` The original `node.label` of the backbone (`NULL` if it had none), to restore at the end.
#'  * `n_tips` Integer. Number of tips of the backbone, used to find the original label of a node.
#'

label_backbone_nodes <- function (phy)
{
  original_node_labels <- phy$node.label
  n_tips <- length(phy$tip.label)
  phy$node.label <- paste0("bb", n_tips + seq_len(phy$Nnode))
  list(phy = phy, original_node_labels = original_node_labels, n_tips = n_tips)
}


#' @title Restore the original node labels of the backbone after grafting
#'
#' @description Replace the labels used to follow the nodes of the backbone (see `label_backbone_nodes()`) by their original labels.
#'   * If the backbone had no node labels, the grafted tree has no node labels either.
#'   * If the backbone had node labels (e.g., support values), its nodes keep them, and the nodes created by grafting are labelled with an empty string (`""`).
#'
#' @param phy Phylogeny of class `"phylo"`, grafted from a backbone labelled with `label_backbone_nodes()`.
#' @param labelled Output of `label_backbone_nodes()` on the backbone.
#'
#' @noRd
#'
#' @return The phylogeny with its original node labels.
#'

restore_node_labels <- function (phy, labelled)
{
  if (is.null(labelled$original_node_labels))
  {
    phy$node.label <- NULL
    return(phy)
  }
  labels <- phy$node.label
  backbone_node <- suppressWarnings(as.integer(substring(labels, 3)))
  is_backbone_node <- !is.na(labels) & substr(labels, 1, 2) == "bb" & !is.na(backbone_node)
  restored <- rep("", length(labels))
  restored[is_backbone_node] <- labelled$original_node_labels[backbone_node[is_backbone_node] - labelled$n_tips]
  phy$node.label <- restored
  phy
}
