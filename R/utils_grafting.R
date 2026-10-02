##### Internal helpers shared by grafting_missing_taxa() and grafting_missing_clades() #####

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


## Force ultrametry only when needed
# (nnls method is very costly: dense ntips x nedges matrix)
ensure_ultrametric <- function (phy, method = "extend", tol = 1e-8)
{
  if (ape::is.ultrametric(phy, tol = tol)) return(phy)
  phytools::force.ultrametric(tree = phy, method = method, message = FALSE)
}


## Fit a birth-death model on the clade at `node`.
# If it has less tips than min_tips, climb towards the root until it has enough (stops at the root in any case).
# Rationale = have enough divergence times to be able to fit a relevant BD model
# Returns list(b, d, n_tips, node).
fit_bd_on_clade <- function (phy, node, n_missing = 1, min_tips = 3)
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


## Identify similar inclusive/exclusive clades and order them by nestedness

# Group missing taxa by the node of their MRCA,
# then order groups by size so nested ones come first
# (a clade always has strictly more descendants than any clade nested within it => sorting on clade size is enough).
group_and_order_clades <- function (phy, MRCA_tips)
{
  # Extract (MRCA) nodes associated to clades
  nodes <- vapply(MRCA_tips, function (tips) get_clade_node(phy, tips), numeric(1))
  # Group them based on similarity of (MRCA) nodes
  group <- match(nodes, unique(nodes))
  # Record first appearance of each group in the list
  first <- match(seq_along(unique(nodes)), group)

  # Compute number of descendants in each group
  n_desc <- vapply(unique(nodes), function (n) length(get_descendant_nodes(phy, n)), numeric(1))
  # Record information per taxa
  groups_per_taxa <- list(group = group, # Order
                          MRCA_tips = MRCA_tips[first], # Tips
                          order = order(n_desc)) # Order to use in the loop

  # Export output
  return(groups_per_taxa)
}

