
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

## Remove the updated simmap from a grafted tree, to compare it with a tree grafted without simmap
drop_simmap <- function (phy) { phy$simmap <- NULL; phy }

## Expect a grafted tree to be valid: ultrametric, backbone preserved and constraints respected
expect_valid_graft <- function (grafted, phy, missing_taxa_list, ...)
{
  expect_true(ape::is.ultrametric(grafted))
  expect_true(all(missing_taxa_list$Taxon_name %in% grafted$tip.label))
  checks <- suppressWarnings(check_grafting_constraints(grafted, phy, missing_taxa_list, verbose = FALSE, ...))
  expect_true(isTRUE(attr(checks, "all_ok")))
}


## --- Toy simmaps --- ##

# Used to test the grafting restricted to states

## Function to  produce a toy simmap on a backbone.
# Includes three states with several changes along the branches (some branches have several pieces).
# The RNG state of the session is restored on exit.
make_toy_simmap <- function (phy, seed = 1, rate = 0.04)
{
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  Q <- matrix(rate, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  diag(Q) <- -2 * rate
  phytools::sim.history(phy, Q, anc = "A", message = FALSE)
}

## Relabel the states of a simmap (named vector: old state = new state), e.g., to turn the toy states into biogeographic ranges
relabel_simmap <- function (simmap, mapping)
{
  simmap$maps <- lapply(simmap$maps, function (x) { names(x) <- unname(mapping[names(x)]); x })
  colnames(simmap$mapped.edge) <- unname(mapping[colnames(simmap$mapped.edge)])
  simmap
}

## Toy simmap with biogeographic ranges: states A, AB and BC (the toy states A, B and C relabelled)
make_toy_range_simmap <- function (phy, seed = 1, rate = 0.04)
{
  relabel_simmap(make_toy_simmap(phy, seed = seed, rate = rate), c(A = "A", B = "AB", C = "BC"))
}

## Is each range contained in the state at the same position? (vectorised version of range_in_states())
range_in_states_pairwise <- function (ranges, states) { mapply(function (range, state) range_in_states(range, state), ranges, states, USE.NAMES = FALSE) }

## State of each grafted taxon in the returned simmap, and state of its parent branch at the graft point (end of the branch above the node where the taxon is attached)
# Returns a data.frame with one row per taxon: taxon, state (the range inherited), parent_state
graft_states_in_simmap <- function (grafted, taxa)
{
  edge <- grafted$edge
  one_taxon <- function (taxon)
  {
    tip_edge <- match(match(taxon, grafted$tip.label), edge[, 2])
    parent_edge <- match(edge[tip_edge, 1], edge[, 2])
    parent_map <- if (is.na(parent_edge)) NULL else grafted$simmap$maps[[parent_edge]]
    data.frame(taxon = taxon, state = names(grafted$simmap$maps[[tip_edge]])[1],
               parent_state = if (is.null(parent_map)) NA_character_ else names(parent_map)[length(parent_map)], stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(taxa, one_taxon))
}

## Function to create a backbone phylogeny and its associated simmaps
toy_simmap_setup <- function ()
{
  phy <- make_toy_phy(n_tips = 20)
  list(phy = phy, simmap = make_toy_simmap(phy, seed = 1))
}


## Independent computation of the state of the simmap at a given age along the lineage leading to a set of backbone tips
# (the lineage is the backbone branch leading to the MRCA of the tips, or to the tip if only one)
# This is used to check that grafted taxa inherited the state of the branch they were grafted on
state_of_lineage_at_age <- function (simmap, backbone_tips, age)
{
  node <- if (length(backbone_tips) == 1) match(backbone_tips, simmap$tip.label) else ape::getMRCA(simmap, backbone_tips)
  edge_row <- match(node, simmap$edge[, 2])
  root_age <- max(phytools::nodeHeights(simmap)[, 2])
  tip_age <- root_age - phytools::nodeHeights(simmap)[edge_row, 2]
  distance_from_tip <- age - tip_age
  branch_length <- simmap$edge.length[edge_row]
  map <- simmap$maps[[edge_row]] # from the rootward node to the tipward node
  distance_from_root <- branch_length - distance_from_tip
  ends <- cumsum(map)
  names(map)[which(distance_from_root <= ends + 1e-9)[1]]
}


## Independent computation of the state of the simmap at the point where a taxon was grafted:
# the lineage is identified by the backbone tips below the sister branch of the taxon in the grafted tree, at the age of the node of the graft.
# NA if the taxon was grafted on a branch of other grafted taxa (no backbone tip below it)
state_at_graft <- function (grafted, simmap, taxon, missing_taxa)
{
  n_tips <- ape::Ntip(grafted)
  tip <- match(taxon, grafted$tip.label)
  parent <- grafted$edge[match(tip, grafted$edge[, 2]), 1]
  sister <- setdiff(grafted$edge[grafted$edge[, 1] == parent, 2], tip)
  below <- if (sister <= n_tips) sister else phytools::getDescendants(grafted, sister)
  backbone_tips <- setdiff(grafted$tip.label[below[below <= n_tips]], missing_taxa)
  if (length(backbone_tips) == 0) { return(NA_character_) }
  age <- max(phytools::nodeHeights(grafted)[, 2]) - phytools::nodeheight(grafted, parent)
  state_of_lineage_at_age(simmap, backbone_tips, age)
}

## Age of the node where a grafted taxon (or the root of a grafted clade) is attached to the rest of the tree
graft_age <- function (grafted, taxa)
{
  node <- if (length(taxa) == 1) match(taxa, grafted$tip.label) else ape::getMRCA(grafted, taxa)
  parent <- grafted$edge[match(node, grafted$edge[, 2]), 1]
  max(phytools::nodeHeights(grafted)[, 2]) - phytools::nodeheight(grafted, parent)
}

## Change the states of a branch of a toy simmap: the branch leading to the MRCA of the tips (or to the tip),
# cut in pieces of the given proportions of the branch length, from the rootward end to the tipward end
set_branch_states <- function (simmap, tips, states, proportions)
{
  node <- if (length(tips) == 1) match(tips, simmap$tip.label) else ape::getMRCA(simmap, tips)
  edge_row <- match(node, simmap$edge[, 2])
  branch_length <- simmap$edge.length[edge_row]
  simmap$maps[[edge_row]] <- stats::setNames(branch_length * proportions / sum(proportions), states)
  simmap
}

## Build a toy simmap with the state "B" in the youngest part of the tree only (ages below 'age'), and the state "A" before
# All the branches crossing the age are cut in two pieces. Used to test narrow and young windows of allowed states.
make_young_state_simmap <- function (phy, age)
{
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  ages <- root_age - phytools::nodeHeights(phy) # rootward, tipward
  maps <- lapply(seq_len(nrow(phy$edge)), function (e)
  {
    if (ages[e, 2] >= age) { return(stats::setNames(phy$edge.length[e], "A")) }
    if (ages[e, 1] <= age) { return(stats::setNames(phy$edge.length[e], "B")) }
    stats::setNames(c(ages[e, 1] - age, age - ages[e, 2]), c("A", "B"))
  })
  simmap <- phy
  simmap$maps <- maps
  class(simmap) <- c("simmap", "phylo")
  simmap
}
