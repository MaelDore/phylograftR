##### Internal helpers shared by graft_missing_taxa() and graft_missing_clades() #####

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
#'   * The stem region of the clade (see `get_stem_region_nodes()`), if `include_stem = TRUE`.
#'     A taxon grafted on the stem region becomes sister to the clade, or sister to other missing taxa already grafted on the stem.
#'
#'   Single-tip clades always use their stem region (their tip branch is their only branch), whatever the value of `include_stem`.
#'   The root has no stem, thus no stem region.
#'
#' @param phy Phylogeny of class `"phylo"`, possibly including taxa already grafted.
#' @param MRCA_node Integer. Node ID of the clade (tip node for a single-tip clade).
#' @param include_stem Logical. Whether the stem region of the clade is available. Default is `TRUE`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not). Default is `character(0)`.
#'
#' @importFrom ape Ntip
#'
#' @noRd
#'
#' @return Integer vector of tipward node IDs.
#'

get_available_branch_nodes <- function (phy, MRCA_node, include_stem = TRUE, missing_tip_labels = character(0))
{
  # All branches within the clade
  available_nodes <- get_descendant_nodes(phy = phy, node = MRCA_node)

  # Single-tip clades: the tip branch is the only branch available (it is the stem of the clade), so include_stem is ignored (set to TRUE)
  if (MRCA_node <= ape::Ntip(phy)) { include_stem <- TRUE }

  # Add the stem region if requested. It is empty for the root, which has no stem
  if (include_stem)
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
#' @param n_missing Integer. Number of missing taxa, used for the sampling fraction. Default is `1`.
#' @param min_tips Integer. Minimum number of tips for the fit. Default is `5`.
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
#' @return A list with 4 elements:
#'  * `b` Numeric. Birth rate = speciation rate in events per time per lineage.
#'  * `d` Numeric. Death rate = extinction rate in events per time per lineage.
#'  * `n_tips` Integer. Number of tips in the enclosing clade used for fitting.
#'  * `node` Integer. MRCA node ID of the selected enclosing clade.
#'

fit_bd_on_clade <- function (phy, node, n_missing = 1, min_tips = 5)
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

  ## Fit BD with sampling fraction accounting for missing tips (the tips we are aiming to graft)
  rho <- n_tips / (n_tips + n_missing) # Compute sampling fraction
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
  BD_output <- list(b = b, d = d, n_tips = n_tips, node = clade_node)

  # Export output
  return(BD_output)
}


#' @title Group missing taxa by clade and order the clades by nestedness
#'
#' @description Group missing taxa by the node of their MRCA, then order groups by size so nested ones come first
#'   (a clade always has strictly more descendants than any clade nested within it,
#'   thus sorting on clade size is enough to list nested clades first).
#'   The rationale is to graft taxa in nested clades first, so they are accounted for when grafting missing taxa in higher clades
#'   comprising those nested clades.
#'
#'   Use tips as inputs to define the associated clades:
#'     * One tip defines a terminal branch.
#'     * Two or more tips define a clade based on their MRCA.
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param clade_tips List of character vectors defining the clade of each taxon.
#'
#' @noRd
#'
#' @return A list with 3 elements:
#'  * `group` Integer vector with one value per taxon. ID of the group (i.e., unique clade) the taxon belongs to.
#'  * `clade_tips` List with one element per group. Character vector of the tips used to identify the associated clade.
#'  * `order` Integer vector with one value per group. Group IDs ranked by clade size, so nested clades come first.
#'

group_and_order_clades <- function (phy, clade_tips)
{
  # Extract (MRCA) nodes associated to clades
  nodes <- vapply(clade_tips, function (tips) get_clade_node(phy, tips), numeric(1))
  # Group them based on similarity of (MRCA) nodes
  group <- match(nodes, unique(nodes))
  # Record first appearance of each group in the list
  first <- match(seq_along(unique(nodes)), group)

  # Compute number of descendants in each group
  n_desc <- vapply(unique(nodes), function (n) length(get_descendant_nodes(phy, n)), numeric(1))
  # Record information per taxa
  groups_per_taxa <- list(group = group, # Group ID
                          clade_tips = clade_tips[first], # Tips
                          order = order(n_desc)) # Order of the group to use in the loop

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
#'   there are four possible relationships (single-tip subclades exclude nothing and are ignored):
#'   * `"nested"`: The excluded subclade is within the target clade (expected).
#'   * `"identical"`: The excluded subclade is the target clade itself. The whole crown is excluded, so no branch (or only the stem) is left.
#'   * `"disjoint"`: The excluded subclade does not overlap with the target clade. The exclusion has no effect.
#'   * `"contains_target"`: The excluded subclade contains the target clade. All branches of the target clade (stem included) are excluded, so no branch is available.
#'
#'   Only the problematic relationships are returned.
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param missing_taxa_list List with `Taxa_name`, `Target_clade` and `Excluded_subclades`.
#'
#' @importFrom ape getMRCA
#'
#' @noRd
#'
#' @return A data.frame with one row per problematic excluded subclade (0 rows if none), and the columns:
#'  * `Taxa_name` Character string. Name of the missing taxon.
#'  * `Target_clade` Character string. Tips defining the target clade, separated by commas.
#'  * `Excluded_subclade` Character string. Tips defining the excluded subclade, separated by commas.
#'  * `relation` Character string. One of `"identical"`, `"disjoint"` or `"contains_target"`.
#'

check_excluded_subclades_position <- function (phy, missing_taxa_list)
{
  n_missing <- length(missing_taxa_list$Taxa_name)

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
    # Excluded subclades with at least 2 tips (single-tip or empty subclades exclude nothing)
    excluded_i <- Filter(function (x) length(x) >= 2, missing_taxa_list$Excluded_subclades[[i]])
    if (length(excluded_i) == 0) next

    target_clade_tips_i <- missing_taxa_list$Target_clade[[i]]
    target_node_i <- get_node(target_clade_tips_i)

    for (tips in excluded_i)
    {
      excluded_node <- get_node(tips)

      # Compare the position of the excluded subclade to the target clade
      relation <- if (excluded_node == target_node_i) {
        "identical" # The excluded subclade is the target clade. ERROR: only the stem branch remain (if include_stem = TRUE) ; advise to use graft_missing_clades with the target calde as Sister_clade instead.
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
          Taxa_name = missing_taxa_list$Taxa_name[i],
          Target_clade = paste(target_clade_tips_i, collapse = ", "),
          Excluded_subclade = paste(tips, collapse = ", "),
          relation = relation, stringsAsFactors = FALSE)
      }
    }
  }

  # Return a df, with 0 rows if no problem was found
  if (length(out) == 0)
  {
    return(data.frame(Taxa_name = character(0), Target_clade = character(0),
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
#'      all the descendant nodes of the MRCA node of the clade, plus the MRCA node itself if the stem is allowed for grafting (`include_stem = TRUE`).
#'   2. Remove from the pool the branches inside each excluded subclade (= all descendant nodes of the subclade's MRCA node).
#'   3. If the pool is empty, the taxon is flagged.
#'
#' It mirrors the rules of the grafting functions:
#'   - the stem of the target clade is available if `include_stem = TRUE` (never for the root, which has no stem)
#'   - the stem of a single-tip clade is always available (it is its only branch), ignoring `include_stem`
#'   - branches inside excluded subclades are unavailable, but their own stem stays available
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param missing_taxa_list List with `Taxa_name`, `Target_clade` and `Excluded_subclades`.
#' @param include_stem Logical. Default is `TRUE`.
#' @param use_excluded_subclades Logical. Default is `TRUE`.
#'
#' @importFrom ape getMRCA
#'
#' @noRd
#'
#' @return A data.frame with two columns summarizing the flagged taxa. It has 0 rows if no taxa were flagged.
#'   * `Taxa_name` Character string. Name of the missing taxon to graft.
#'   * `Target_clade` Character string. Tips used to describe the target clade associated with the missing taxon; separated by commas if multiple.
#'


find_taxa_without_allowed_branch <- function (phy,
                                              missing_taxa_list,
                                              include_stem = TRUE,
                                              use_excluded_subclades = TRUE)
{
  n_missing <- length(missing_taxa_list$Taxa_name)
  root_node <- setdiff(phy$edge[, 1], phy$edge[, 2])

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
    stem_available <- (include_stem || is_single_tip) && (MRCA_node_i != root_node)

    # Identify the tipward nodes associated with the valid branches (from cache env. or compute for the first time)
    allowed_nodes <- c(cached(paste0("desc_", MRCA_node_i), get_descendant_nodes(phy, MRCA_node_i)),
                       if (stem_available) MRCA_node_i)

    # Remove branches from excluded subclades
    if (use_excluded_subclades)
    {
      # Excluded subclades with at least 2 tips (single-tip or empty subclades exclude nothing)
      for (tips in Filter(function (x) length(x) >= 2, missing_taxa_list$Excluded_subclades[[i]]))
      {
        # Identify the tipward nodes associated with the INvalid branches (from cache env. or compute for the first time)
        excluded_node <- cached(paste0("mrca_", paste(sort(tips), collapse = "|")), ape::getMRCA(phy, tips))
        # Remove them from available nodes
        allowed_nodes <- setdiff(allowed_nodes, cached(paste0("desc_", excluded_node), get_descendant_nodes(phy, excluded_node)))
      }
    }
    # The taxa if flagged only if no valid branches (tipward nodes) remains in the pool
    flagged[i] <- length(allowed_nodes) == 0
  }

  # Build data.frame summarizing all taxa that did not pass the check
  # If no taxon is flagged, return an empty df with 0 rows.
  flagged_taxa_df <- data.frame(
    Taxa_name = missing_taxa_list$Taxa_name[flagged],
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
