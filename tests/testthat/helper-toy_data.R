
##### Helpers shared by all the test files #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Loaded automatically by testthat
# Produce small toy phylogenies and clades: kept tiny so tests stay fast.


## Toy backbone: random coalescent tree, ultrametric, with a given root age.
# The RNG state of the session is restored on exit.
make_toy_phy <- function (n_tips = 20, seed = 1234, root_age = 50)
{
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  phy <- ape::rcoal(n_tips)
  phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * root_age
  phy
}

## Build a backbone with polytomies by collapsing internal branches (the given number of them, starting from the highest node IDs).
# The length of a collapsed branch is added to the branches below it, so the backbone stays ultrametric.
make_toy_polytomies <- function (phy, n_collapsed = 1)
{
  internal_edges <- which(phy$edge[, 2] > ape::Ntip(phy))
  collapsed <- internal_edges[order(phy$edge[internal_edges, 2], decreasing = TRUE)][seq_len(n_collapsed)]
  for (e in collapsed)
  {
    below <- which(phy$edge[, 1] == phy$edge[e, 2])
    phy$edge.length[below] <- phy$edge.length[below] + phy$edge.length[e]
    phy$edge.length[e] <- 0
  }
  ape::di2multi(phy, tol = 1e-8)
}

# Two tips whose MRCA is the given node (a clade is defined by a pair of tips)
pair_for_node <- function (phy, node)
{
  tips <- phy$tip.label[ape::prop.part(phy)[[node - ape::Ntip(phy)]]]
  pairs <- utils::combn(tips, 2)
  for (i in seq_len(ncol(pairs)))
  {
    if (ape::getMRCA(phy, pairs[, i]) == node) return(pairs[, i])
  }
}

## Toy setup: backbone + a clade A that contains two disjoint nested clades B and C,
# + a clade D disjoint from A, and a single tip 's' in A (neither in B nor C) and a single tip 'o' in D.
# Clades are given as tips (pairs of tips), as expected by the grafting functions.
toy_setup <- function (n_tips = 20, seed = 1234)
{
  phy <- make_toy_phy(n_tips = n_tips, seed = seed)
  n <- ape::Ntip(phy)
  parts <- ape::prop.part(phy)
  size <- lengths(parts)
  desc_nodes <- function (node) { d <- phytools::getDescendants(phy, node); d[d > n] }
  for (a in n + which(size >= 8 & size <= 12))
  {
    inside <- desc_nodes(a)
    cand_B <- inside[size[inside - n] %in% 3:6]
    for (b in cand_B)
    {
      cand_C <- setdiff(inside[size[inside - n] %in% 3:6], c(b, desc_nodes(b)))
      cand_C <- cand_C[!vapply(cand_C, function (x) b %in% desc_nodes(x), logical(1))]
      if (length(cand_C) > 0)
      {
        cc <- cand_C[1]
        tipsA <- phy$tip.label[parts[[a - n]]]
        tipsB <- phy$tip.label[parts[[b - n]]]
        tipsC <- phy$tip.label[parts[[cc - n]]]
        other <- setdiff(phy$tip.label, tipsA)

        # A clade D disjoint from A: not nested in A, and not containing A
        anc_A <- a; while (!is.na(par <- phy$edge[match(anc_A[length(anc_A)], phy$edge[, 2]), 1])) anc_A <- c(anc_A, par)
        cand_D <- setdiff(n + which(size >= 2), c(anc_A, inside))
        d <- cand_D[1]
        return(list(phy = phy, root_age = max(ape::node.depth.edgelength(phy)),
                    A = pair_for_node(phy, a), B = pair_for_node(phy, b), C = pair_for_node(phy, cc),
                    nA = length(tipsA), nB = length(tipsB), nC = length(tipsC),
                    tipsA = tipsA, tipsB = tipsB, tipsC = tipsC,
                    D = pair_for_node(phy, d), s = setdiff(tipsA, c(tipsB, tipsC))[1], o = other[1],
                    node_A = a, node_B = b, node_C = cc, node_D = d))
      }
    }
  }
  stop("toy_setup(): no suitable clades found, try another seed")
}

## Graft a tip on a branch (the branch leading to 'node'), at a fraction of the branch length from its rootward end measured from the tipward node.
# The new tip reaches the present (ultrametric)
graft_tip_on_branch <- function (phy, label, node, frac = 0.5)
{
  branch_length <- phy$edge.length[match(node, phy$edge[, 2])]
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  phytools::bind.tip(phy, tip.label = label, where = node, position = branch_length * frac,
                     edge.length = branch_length * frac + (root_age - phytools::nodeheight(phy, node)))
}

## Run a grafting call and return the result with its warnings silenced
quiet_graft <- function (fun, ...) { suppressWarnings(fun(..., verbose = FALSE)) }

## Expect a grafted tree to be valid: ultrametric, backbone preserved and constraints respected
expect_valid_graft <- function (grafted, phy, missing_taxa_list, ...)
{
  expect_true(ape::is.ultrametric(grafted))
  expect_true(all(missing_taxa_list$Taxon_name %in% grafted$tip.label))
  checks <- suppressWarnings(check_grafting_constraints(grafted, phy, missing_taxa_list, verbose = FALSE, ...))
  expect_true(isTRUE(attr(checks, "all_ok")))
}
