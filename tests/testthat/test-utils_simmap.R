
##### Test suite for functions handling simmaps #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


##### 1/ Helper functions #####

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

## Function to create abackbone phylogeny and its associated simmaps
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


##### 2/ Tests simmap-related functions #####

test_that("validate_simmap() accepts a matching simmap and returns its states and the matching of branches",
{
  s <- toy_simmap_setup()
  checked <- validate_simmap(s$simmap, s$phy)
  expect_equal(checked$states, c("A", "B", "C"))
  expect_equal(checked$edge_map, seq_len(nrow(s$phy$edge)))

  # A backbone with its tips and branches in another order is matched through the topology
  # (rows of the edge matrix in reverse order, and tips renumbered in reverse order)
  shuffled <- s$simmap
  n <- ape::Ntip(shuffled)
  renumber <- function (x) ifelse(x <= n, n + 1 - x, x)
  perm <- rev(seq_len(nrow(shuffled$edge)))
  shuffled$edge <- matrix(renumber(shuffled$edge[perm, ]), ncol = 2)
  shuffled$edge.length <- shuffled$edge.length[perm]
  shuffled$maps <- shuffled$maps[perm]
  shuffled$tip.label <- rev(shuffled$tip.label)
  expect_false(identical(shuffled$edge, s$simmap$edge))
  checked_shuffled <- validate_simmap(shuffled, s$phy)
  expect_equal(shuffled$edge.length[checked_shuffled$edge_map], s$phy$edge.length)
  tips_below <- function (tree, node) sort(tree$tip.label[phytools::getDescendants(tree, node)[phytools::getDescendants(tree, node) <= ape::Ntip(tree)]])
  for (i in c(1, 7, 20, 30))
  {
    expect_equal(tips_below(shuffled, shuffled$edge[checked_shuffled$edge_map[i], 2]), tips_below(s$phy, s$phy$edge[i, 2]))
  }
})

test_that("validate_simmap() rejects objects that are not a single simmap of the backbone",
{
  s <- toy_simmap_setup()
  expect_error(validate_simmap(s$phy, s$phy), "stochastic character map")
  expect_error(validate_simmap("A", s$phy), "stochastic character map")
  multi <- structure(list(s$simmap, s$simmap), class = c("multiSimmap", "multiPhylo"))
  expect_error(validate_simmap(multi, s$phy), "single stochastic character map")

  # Other tips
  other_tips <- s$simmap; other_tips$tip.label[1] <- "not_a_tip"
  expect_error(validate_simmap(other_tips, s$phy), "tips of 'simmap' differ")

  # Other topology
  other_topology <- make_toy_phy(n_tips = 20, seed = 99)
  other_topology$tip.label <- s$phy$tip.label
  expect_error(validate_simmap(make_toy_simmap(other_topology), s$phy), "topology|branch lengths")

  # Other branch lengths
  longer <- s$simmap; longer$edge.length <- longer$edge.length * 1.1
  expect_error(validate_simmap(longer, s$phy), "branch lengths")

  # Maps not summing to the branch lengths, or invalid
  bad_sum <- s$simmap; bad_sum$maps[[3]] <- bad_sum$maps[[3]] * 0.5
  expect_error(validate_simmap(bad_sum, s$phy), "do not sum")
  unnamed <- s$simmap; unnamed$maps[[3]] <- unname(unnamed$maps[[3]])
  expect_error(validate_simmap(unnamed, s$phy), "named numeric")
})

test_that("check_allowed_states() lists the unknown states and the states available",
{
  expect_identical(check_allowed_states(c("A", "B"), c("A", "B", "C")), c("A", "B"))
  expect_error(check_allowed_states(c("A", "D", "E"), c("A", "B", "C"), who = "taxon 'X1'"), "taxon 'X1'.*'D', 'E'.*'A', 'B', 'C'")
  expect_error(check_allowed_states(character(0), c("A")), "non-empty")
  expect_error(check_allowed_states(NA_character_, c("A")), "non-empty|NA")
})

test_that("build_state_table() describes the pieces of each branch from its tipward node, and matches the simmap",
{
  s <- toy_simmap_setup()
  checked <- validate_simmap(s$simmap, s$phy)
  tab <- build_state_table(s$phy, s$simmap, checked$edge_map)

  # Total duration per state is the one of the simmap
  per_state <- tapply(tab$to - tab$from, tab$state, sum)
  expect_equal(as.numeric(per_state[colnames(s$simmap$mapped.edge)]), unname(colSums(s$simmap$mapped.edge)))

  # Every branch is fully covered, from 0 to its length, without gap or overlap
  for (node in s$phy$edge[, 2])
  {
    piece <- tab[tab$node == node, ]
    len <- s$phy$edge.length[match(node, s$phy$edge[, 2])]
    expect_identical(min(piece$from), 0)
    expect_identical(max(piece$to), len)
    expect_equal(piece$from[-1], piece$to[-nrow(piece)])
  }

  # The tipward piece is the last of the map (maps go from the rootward node to the tipward node)
  multi <- which(lengths(s$simmap$maps) > 1)[1]
  node <- s$phy$edge[multi, 2]
  piece <- tab[tab$node == node, ]
  expect_equal(piece$state, rev(names(s$simmap$maps[[multi]])))
})

test_that("get_allowed_segments_tracked() returns whole branches when all states are allowed, and only the allowed pieces otherwise",
{
  s <- toy_simmap_setup()
  started <- start_state_tracking(s$phy, s$simmap)
  phy <- started$phy
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy)
  nodes <- phy$edge[, 2]

  all_states <- get_allowed_segments_tracked(started$tracker, phy, nodes, c("A", "B", "C"), edge_ages)
  expect_identical(all_states, get_branch_segments(phy, nodes))

  only_A <- get_allowed_segments_tracked(started$tracker, phy, nodes, "A", edge_ages)
  expect_equal(sum(only_A$to - only_A$from), unname(colSums(s$simmap$mapped.edge)["A"]))
  expect_true(all(only_A$to > only_A$from))
  # The order of the input nodes is kept
  reversed <- get_allowed_segments_tracked(started$tracker, phy, rev(nodes), c("A", "B"), edge_ages)
  expect_equal(sum(reversed$to - reversed$from), sum(colSums(s$simmap$mapped.edge)[c("A", "B")]), ignore_attr = TRUE)
  expect_false(is.unsorted(match(reversed$node, rev(nodes))))
})

test_that("node labels are labelled for tracking and restored at the end, with or without original labels",
{
  s <- toy_simmap_setup()
  expect_null(s$phy$node.label)
  started <- start_state_tracking(s$phy, s$simmap)
  expect_equal(started$phy$node.label, paste0("bb", ape::Ntip(s$phy) + seq_len(s$phy$Nnode)))
  expect_null(stop_state_tracking(started$tracker, started$phy)$node.label)

  # With support values
  phy_support <- s$phy; phy_support$node.label <- paste0("support", seq_len(phy_support$Nnode))
  started <- start_state_tracking(phy_support, s$simmap)
  expect_equal(started$tracker$original_node_labels, phy_support$node.label)
  restored <- stop_state_tracking(started$tracker, started$phy)
  expect_identical(restored$node.label, phy_support$node.label)
})

## Graft taxa at random positions among the allowed segments, following the tracking, and check the tracking against independent computations
run_random_grafts <- function (n_grafts, allowed_per_graft, seed = 1, support = FALSE)
{
  s <- toy_simmap_setup()
  backbone <- s$phy
  if (support) backbone$node.label <- paste0("support", seq_len(backbone$Nnode))
  started <- start_state_tracking(backbone, s$simmap)
  phy <- started$phy; tracker <- started$tracker
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  expected_grafted_length <- c(A = 0, B = 0, C = 0) # Length of the grafted branches in each state, computed independently
  set.seed(seed)
  log <- list()
  for (i in seq_len(n_grafts))
  {
    allowed <- allowed_per_graft[[(i - 1) %% length(allowed_per_graft) + 1]]
    node_heights <- phytools::nodeHeights(phy)
    edge_ages <- root_age - node_heights
    nodes <- phy$edge[, 2]
    segments <- get_allowed_segments_tracked(tracker, phy, nodes, allowed, edge_ages)
    drawn <- draw_point_in_segments(segments)
    node <- drawn$node; position <- drawn$position

    # Independent state at the graft point: backbone tips below the node in the current tree, and the age of the point
    graft_age <- root_age - phytools::nodeheight(phy, node) + position
    below <- phytools::getDescendants(phy, node); below <- below[below <= ape::Ntip(phy)]
    below_backbone <- intersect(phy$tip.label[below], s$phy$tip.label)
    expected_state <- if (length(below_backbone) > 0) state_of_lineage_at_age(s$simmap, below_backbone, graft_age) else NA

    new_length <- position + (root_age - phytools::nodeheight(phy, node))
    label <- paste0("X", i)
    phy_after <- phytools::bind.tip(phy, tip.label = label, edge.length = new_length, where = node, position = position)
    tracked <- track_graft(tracker, phy, phy_after, node, position, label, allowed_states = allowed, edge_ages = edge_ages)
    phy <- tracked$phy; tracker <- tracked$tracker

    state <- unname(tracker$missing_state[label])
    if (!is.na(expected_state)) expect_identical(state, expected_state) # grafted on a backbone lineage
    expect_true(state %in% allowed) # grafting happened in an allowed state
    expected_grafted_length[state] <- expected_grafted_length[state] + new_length
    log[[i]] <- list(on_missing_only = is.na(expected_state))
  }

  # Final consistency check on the tracking: the total length in each state over the whole tree is the length of the backbone plus the grafted branches
  node_heights <- phytools::nodeHeights(phy)
  edge_ages <- root_age - node_heights
  for (st in c("A", "B", "C"))
  {
    seg <- get_allowed_segments_tracked(tracker, phy, phy$edge[, 2], st, edge_ages)
    expect_equal(sum(seg$to - seg$from), unname(colSums(s$simmap$mapped.edge)[st] + expected_grafted_length[st]), tolerance = 1e-8)
  }
  list(phy = phy, tracker = tracker, simmap = s$simmap, backbone = backbone, log = log)
}

test_that("tracking follows the states through successive grafts, including grafts on previously grafted branches",
{
  res <- run_random_grafts(60, list(c("A", "B", "C")), seed = 1)
  expect_true(any(vapply(res$log, function (x) x$on_missing_only, logical(1)))) # some grafts were on missing-only subtrees
  expect_equal(ape::Ntip(res$phy), 20 + 60)
  expect_true(ape::is.ultrametric(res$phy))
})

test_that("tracking follows the states when the allowed states change from one graft to the next",
{
  allowed <- list("A", c("B", "C"), "B", c("A", "C"), "C", c("A", "B"))
  res <- run_random_grafts(60, allowed, seed = 2)
  # The grafted taxa are in the state of their graft point, so they are in one of the allowed states of their graft
  states <- res$tracker$missing_state
  for (i in seq_along(states)) expect_true(states[[paste0("X", i)]] %in% allowed[[(i - 1) %% length(allowed) + 1]])
})

test_that("restoring the labels after grafting gives the original node labels, NA for the new nodes, and no labels if there were none",
{
  res <- run_random_grafts(15, list(c("A", "B", "C")), seed = 3)
  expect_null(stop_state_tracking(res$tracker, res$phy)$node.label)

  res <- run_random_grafts(15, list(c("A", "B", "C")), seed = 3, support = TRUE)
  restored <- stop_state_tracking(res$tracker, res$phy)
  expect_length(restored$node.label, restored$Nnode)
  expect_equal(sum(!is.na(restored$node.label)), res$backbone$Nnode)
  expect_setequal(restored$node.label[!is.na(restored$node.label)], res$backbone$node.label)
  # Each support value is still on the node subtending the same backbone tips
  for (lab in res$backbone$node.label[1:5])
  {
    node_before <- ape::Ntip(res$backbone) + match(lab, res$backbone$node.label)
    node_after <- ape::Ntip(restored) + match(lab, restored$node.label)
    tips <- function (tree, node) { d <- phytools::getDescendants(tree, node); sort(setdiff(tree$tip.label[d[d <= ape::Ntip(tree)]], paste0("X", 1:15))) }
    expect_equal(tips(restored, node_after), tips(res$backbone, node_before))
  }
})

test_that("get_state_at_position() solves ties at the limit of two pieces with the allowed states",
{
  s <- toy_simmap_setup()
  started <- start_state_tracking(s$phy, s$simmap)
  phy <- started$phy
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy)
  multi <- which(lengths(s$simmap$maps) > 1)[1]
  node <- phy$edge[multi, 2]
  map <- s$simmap$maps[[multi]]
  limit <- phy$edge.length[multi] - map[[1]] # distance from the tipward node to the limit between the first two pieces (from the root)
  state_root_side <- names(map)[1]; state_tip_side <- names(map)[2]
  expect_identical(get_state_at_position(started$tracker, phy, node, limit, allowed_states = state_root_side, edge_ages = edge_ages), state_root_side)
  expect_identical(get_state_at_position(started$tracker, phy, node, limit, allowed_states = state_tip_side, edge_ages = edge_ages), state_tip_side)
})
