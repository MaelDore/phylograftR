##### Internal helpers shared by graft_missing_taxa() and graft_missing_clades() #####

## Identify all descendant nodes from a given node in a tree
# Fast, non-recursive search for descendants (all internal + terminal nodes below `node`)
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


## Identify node ID associated to a clade
# Defined by 1 tip -> tip node
# Defined by >= 2 tips -> MRCA node
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


## Identify the parent node of a node (NA for the root)
# Parent node = tipward node of the branch leading to the focal node
get_parent_node <- function (phy, node)
{
  parent_node <- phy$edge[phy$edge[, 2] == node, 1]
  if (length(parent_node ) == 0)
  {
    parent_node <- NA
  }
  return(parent_node)
}


## Safe sample of one element handling edge cases with one or no elements in x
sample_one <- function (x, prob = NULL)
{
  if (length(x) == 0) stop("No valid option to sample from")
  if (length(x) == 1) return(x)
  x[sample.int(length(x), size = 1, prob = prob)]
}


## Fit a birth-death model on the clade at `node`.
# If it has less tips than min_tips, climb towards the root until it has enough (stops at the root in any case).
# Rationale = have enough divergence times to be able to fit a relevant BD model
# Returns list(b, d, n_tips, node).
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

    # Stop once the current clade have enough tips
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
  d <- min(d, 0.99 * b) # Limit death rate to 99% of birth rate to ensure no net diversification is negative
  # No divergence times can be drawn from a neagtive net diversification process

  # Record rates, nb of tips in selected enclosing clade, and MRCA node of the selected clade
  BD_output <- list(b = b, d = d, n_tips = n_tips, node = clade_node)

  # Export output
  return(BD_output)
}


## Identify similar target/excluded clades and order them by nestedness

# Group missing taxa by the node of their MRCA,
# then order groups by size so nested ones come first
# (a clade always has strictly more descendants than any clade nested within it => sorting on clade size is enough).
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
  groups_per_taxa <- list(group = group, # Order
                          clade_tips = clade_tips[first], # Tips
                          order = order(n_desc)) # Order to use in the loop

  # Export output
  return(groups_per_taxa)
}


## Extract the NEW divergence times returned by TreeSim::corsim()
# corsim() returns the input times 'x' plus 'missing' new times, sorted in decreasing order.
# The new times are recovered by removing ONE instance of each input time from the output
# Contrary to setdiff(), this is robust to a new time being numerically identical to an existing one, and to duplicated input times.
extract_new_divergence_times <- function (all_times, # Output of TreeSim::corsim()
                                          old_times, # Input times 'x' given to TreeSim::corsim()
                                          tol = 1e-8) # Tolerance for the bound checks
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


## Seed the RNG for the duration of the calling function only
# The global RNG state (.Random.seed) is saved and restored when the calling function exits (also on error),
# so setting a seed inside the grafting functions does not alter the random stream of the user's session.
# Does nothing if seed is NULL.
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


## Identify the taxa with no branch available for grafting, given their target clade and excluded subclades
# Static check done on the backbone, before any grafting to save time instead of discovering the problem in the middle of a long grafting run.
#
# How it works, for each taxon:
#   1. Start from the pool of branches of its target clade.
#      A branch is identified by its tipward node (the node it leads to), so the pool is made of
#      all the descendant nodes of the MRCA node of the clade, plus the MRCA node itself if the stem is available.
#   2. Remove from the pool the branches inside each excluded subclade (= all descendant nodes of the subclade's MRCA node).
#   3. If the pool is empty, the taxon is flagged.
#
# It mirrors the rules of the grafting functions:
#   - the stem of the target clade is available if include_stem = TRUE (never for the root, which has no stem)
#   - the stem of a single-tip clade is always available (it is its only branch), ignoring include_stem
#   - branches inside excluded subclades are unavailable, but their own stem stays available
#
# Returns a data.frame with the flagged taxa (0 rows if none).

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
