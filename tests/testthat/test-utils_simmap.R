
##### Test suite for functions handling simmaps #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


##### 1/ Tests simmap-related functions #####

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
  phy <- label_backbone_nodes(s$phy)$phy
  tracker <- start_state_tracking(phy, s$simmap)
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy)
  nodes <- phy$edge[, 2]

  all_states <- get_allowed_segments_tracked(tracker, phy, nodes, c("A", "B", "C"), edge_ages)
  expect_identical(all_states, get_branch_segments(phy, nodes))

  only_A <- get_allowed_segments_tracked(tracker, phy, nodes, "A", edge_ages)
  expect_equal(sum(only_A$to - only_A$from), unname(colSums(s$simmap$mapped.edge)["A"]))
  expect_true(all(only_A$to > only_A$from))
  # The order of the input nodes is kept
  reversed <- get_allowed_segments_tracked(tracker, phy, rev(nodes), c("A", "B"), edge_ages)
  expect_equal(sum(reversed$to - reversed$from), sum(colSums(s$simmap$mapped.edge)[c("A", "B")]), ignore_attr = TRUE)
  expect_false(is.unsorted(match(reversed$node, rev(nodes))))
})

## Test that the tracking requires the nodes of the backbone to be labelled with label_backbone_nodes()
test_that("start_state_tracking() stops if the nodes of the backbone were not labelled for tracking",
{
  s <- toy_simmap_setup()
  expect_error(start_state_tracking(s$phy, s$simmap), "label_backbone_nodes")
  phy_support <- s$phy
  phy_support$node.label <- paste0("support", seq_len(phy_support$Nnode))
  expect_error(start_state_tracking(phy_support, s$simmap), "label_backbone_nodes")
  expect_type(start_state_tracking(label_backbone_nodes(phy_support)$phy, s$simmap), "list")
})

## Test that a branch entirely in allowed states stays a single segment covering exactly the whole branch after grafting (as without simmap)
test_that("get_allowed_segments_tracked() gives the same whole-branch segments as get_branch_segments() all along the grafting when all the states are allowed",
{
  s <- toy_simmap_setup()
  phy <- label_backbone_nodes(s$phy)$phy
  tracker <- start_state_tracking(phy, s$simmap)
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  set.seed(1)
  for (i in 1:15)
  {
    edge_ages <- root_age - phytools::nodeHeights(phy)
    nodes <- phy$edge[, 2]
    segments <- get_allowed_segments_tracked(tracker, phy, nodes, c("A", "B", "C"), edge_ages)
    expect_identical(segments, get_branch_segments(phy, nodes))

    # Graft on a random point, whatever its state (phytools::bind.tip() rounds the branch lengths, so they differ slightly from the ones of the simmap)
    drawn <- draw_point_in_segments(segments)
    new_length <- drawn$position + (root_age - phytools::nodeheight(phy, drawn$node))
    phy_after <- phytools::bind.tip(phy, tip.label = paste0("X", i), edge.length = new_length, where = drawn$node, position = drawn$position)
    tracked <- track_graft(tracker, phy, phy_after, drawn$node, drawn$position, paste0("X", i), edge_ages = edge_ages)
    phy <- tracked$phy
    tracker <- tracked$tracker
  }
})

## Test the resolution of the allowed states of each taxon: taxon-level values replace the global value, and NULL or NA use it
test_that("resolve_allowed_states() uses the taxon-level states, and the global states for NULL or NA, and validates its inputs",
{
  states <- c("A", "B", "C")
  taxa <- c("X1", "X2", "X3", "X4")

  # Global only
  resolved <- resolve_allowed_states(allowed_states = c("A", "B"), states = states, taxon_names = taxa)
  expect_named(resolved, taxa)
  expect_true(all(vapply(resolved, identical, logical(1), c("A", "B"))))

  # Taxon-level values replace the global ones (no intersection); NULL, NA and empty vectors use the global ones
  resolved <- resolve_allowed_states(allowed_states_taxa = list(X1 = "C", X2 = NULL, X3 = NA, X4 = character(0)), allowed_states = c("A", "B"),
                                     states = states, taxon_names = taxa)
  expect_identical(resolved, list(X1 = "C", X2 = c("A", "B"), X3 = c("A", "B"), X4 = c("A", "B")))

  # Duplicated states are removed
  expect_identical(resolve_allowed_states(allowed_states = c("A", "A", "B"), states = states, taxon_names = "X1")$X1, c("A", "B"))

  # Only taxon-level values: fine if all the taxa have one
  resolved <- resolve_allowed_states(allowed_states_taxa = list(X1 = "A", X2 = "B"), states = states, taxon_names = c("X1", "X2"))
  expect_identical(resolved, list(X1 = "A", X2 = "B"))

  # Taxa without any allowed state (neither taxon-level nor global): all the states of the simmap are allowed
  expect_identical(resolve_allowed_states(states = states, taxon_names = taxa), stats::setNames(rep(list(states), 4), taxa))
  resolved <- resolve_allowed_states(allowed_states_taxa = list(X1 = "A", X2 = NULL, X3 = "B", X4 = NA), states = states, taxon_names = taxa)
  expect_identical(resolved, list(X1 = "A", X2 = states, X3 = "B", X4 = states))

  # Unknown states: stop, naming the origin of the value and listing the states of the simmap
  expect_error(resolve_allowed_states(allowed_states = c("A", "D"), states = states, taxon_names = taxa), "global 'allowed_states' argument.*'D'.*'A', 'B', 'C'")
  expect_error(resolve_allowed_states(allowed_states_taxa = list(X1 = "red", X2 = NULL, X3 = NULL, X4 = NULL), allowed_states = "A", states = states, taxon_names = taxa),
               "taxon 'X1'.*'red'")
  # Exact matching
  expect_error(resolve_allowed_states(allowed_states = "a", states = states, taxon_names = taxa), "'a'")

  # Invalid types
  expect_error(resolve_allowed_states(allowed_states = 1, states = states, taxon_names = taxa), "character vector")
  expect_error(resolve_allowed_states(allowed_states_taxa = c(X1 = "A"), allowed_states = "A", states = states, taxon_names = "X1"), "must be a list")
  expect_error(resolve_allowed_states(allowed_states_taxa = list(X1 = "A"), allowed_states = "A", states = states, taxon_names = taxa), "same length")
})

## Graft taxa at random positions among the allowed segments, following the tracking, and check the tracking against independent computations
run_random_grafts <- function (n_grafts, allowed_per_graft, seed = 1, support = FALSE)
{
  s <- toy_simmap_setup()
  backbone <- s$phy
  if (support) backbone$node.label <- paste0("support", seq_len(backbone$Nnode))
  labelled <- label_backbone_nodes(backbone)
  phy <- labelled$phy
  tracker <- start_state_tracking(phy, s$simmap)
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
  list(phy = phy, tracker = tracker, simmap = s$simmap, backbone = backbone, labelled = labelled, log = log)
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

## Test that the node labels used for tracking are replaced by the original ones after grafting with tracking
test_that("restoring the labels after tracked grafts gives the original node labels, and an empty string for the new nodes",
{
  # Backbone without node labels: none after grafting
  res <- run_random_grafts(15, list(c("A", "B", "C")), seed = 3)
  expect_null(restore_node_labels(res$phy, res$labelled)$node.label)

  # Backbone with node labels
  res <- run_random_grafts(15, list(c("A", "B", "C")), seed = 3, support = TRUE)
  restored <- restore_node_labels(res$phy, res$labelled)
  expect_length(restored$node.label, restored$Nnode)
  expect_equal(sum(restored$node.label != ""), res$backbone$Nnode)
  expect_setequal(restored$node.label[restored$node.label != ""], res$backbone$node.label)
  # Each node label is still on the node subtending the same backbone tips
  tips_below <- function (tree, node)
  {
    d <- phytools::getDescendants(tree, node)
    sort(setdiff(tree$tip.label[d[d <= ape::Ntip(tree)]], paste0("X", 1:15)))
  }
  for (lab in res$backbone$node.label[1:5])
  {
    node_before <- ape::Ntip(res$backbone) + match(lab, res$backbone$node.label)
    node_after <- ape::Ntip(restored) + match(lab, restored$node.label)
    expect_equal(tips_below(restored, node_after), tips_below(res$backbone, node_before))
  }
})

test_that("get_state_at_position() solves ties at the limit of two pieces with the allowed states",
{
  s <- toy_simmap_setup()
  phy <- label_backbone_nodes(s$phy)$phy
  tracker <- start_state_tracking(phy, s$simmap)
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy)
  multi <- which(lengths(s$simmap$maps) > 1)[1]
  node <- phy$edge[multi, 2]
  map <- s$simmap$maps[[multi]]
  limit <- phy$edge.length[multi] - map[[1]] # distance from the tipward node to the limit between the first two pieces (from the root)
  state_root_side <- names(map)[1]; state_tip_side <- names(map)[2]
  expect_identical(get_state_at_position(tracker, phy, node, limit, allowed_states = state_root_side, edge_ages = edge_ages), state_root_side)
  expect_identical(get_state_at_position(tracker, phy, node, limit, allowed_states = state_tip_side, edge_ages = edge_ages), state_tip_side)
})


##### 2/ Tests of the toy simmap #####

## Test that the toy simmap can be used with the toy backbone, and is balanced as described in its documentation
test_that("toy_simmap is a valid simmap of toy_backbone with three states of similar importance and a moderate number of transitions",
{
  toy_backbone <- phylograftR::toy_backbone
  toy_simmap <- phylograftR::toy_simmap

  checked <- validate_simmap(toy_simmap, toy_backbone)
  expect_equal(checked$states, c("blue", "green", "red"))
  expect_equal(checked$edge_map, seq_len(nrow(toy_backbone$edge)))

  # Each state covers between a quarter and 40% of the total branch length
  share_of_length <- colSums(toy_simmap$mapped.edge) / sum(toy_simmap$mapped.edge)
  expect_true(all(share_of_length > 0.25 & share_of_length < 0.40))

  # Moderate number of transitions: some but not all the branches include a change of state
  share_of_branches_with_change <- mean(lengths(toy_simmap$maps) > 1)
  expect_true(share_of_branches_with_change >= 0.15 && share_of_branches_with_change <= 0.30)
  expect_lte(max(lengths(toy_simmap$maps)), 4)

  # Each state is found in several tips
  expect_true(all(table(toy_simmap$states) >= 7))
})

test_that("track_graft_subtree() gives the state of the graft point to all the taxa of a grafted clade, and labels the new nodes",
{
  s <- toy_simmap_setup()
  labelled <- label_backbone_nodes(s$phy)
  phy <- labelled$phy
  tracker <- start_state_tracking(phy, s$simmap)
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy)

  # Graft a clade of three taxa on a branch with several states, in the middle of its second piece (from the tipward end)
  table_states <- tracker$state_table
  node <- as.integer(names(which(table(table_states$node) >= 2))[1])
  pieces <- table_states[table_states$node == node, ]
  position <- mean(c(pieces$from[2], pieces$to[2]))
  expected_state <- pieces$state[2]
  subtree <- ape::rcoal(3, tip.label = c("M1", "M2", "M3"))
  subtree$edge.length <- subtree$edge.length / max(ape::node.depth.edgelength(subtree)) * 2
  subtree$root.edge <- 1
  phy_after <- ape::bind.tree(phy, subtree, where = node, position = position)

  tracked <- track_graft_subtree(tracker, phy, phy_after, node, position, c("M1", "M2", "M3"), allowed_states = unique(table_states$state), edge_ages = edge_ages)
  expect_identical(unname(tracked$tracker$missing_state[c("M1", "M2", "M3")]), rep(expected_state, 3))

  # The new nodes are labelled: the attachment node on the backbone branch, and the nodes within the clade (all the labels exist)
  new_labels <- tracked$phy$node.label[!tracked$phy$node.label %in% phy$node.label]
  expect_false(anyNA(tracked$phy$node.label))
  expect_equal(sum(grepl("^gr", new_labels)), 1)
  expect_equal(sum(grepl("^ms", new_labels)), 2) # the two internal nodes of a clade of three taxa
  expect_identical(unique(new_labels[grepl("^ms", new_labels)]), paste0("ms", match(expected_state, tracker$states)))
  expect_identical(new_labels[grepl("^gr", new_labels)], paste0("gr", node))

  # All the branches of the tree are tracked, and the branches of the clade are entirely in the state of the graft point
  edge_ages_after <- root_age - phytools::nodeHeights(tracked$phy)
  clade_nodes <- c(match(c("M1", "M2", "M3"), tracked$phy$tip.label), setdiff(phytools::getDescendants(tracked$phy, ape::getMRCA(tracked$phy, c("M1", "M2", "M3"))), 1:ape::Ntip(tracked$phy)))
  clade_nodes <- c(clade_nodes, ape::getMRCA(tracked$phy, c("M1", "M2", "M3")))
  segments <- get_allowed_segments_tracked(tracked$tracker, tracked$phy, clade_nodes, expected_state, edge_ages_after)
  expect_equal(nrow(segments), length(unique(clade_nodes)))
  expect_equal(segments$to - segments$from, tracked$phy$edge.length[match(segments$node, tracked$phy$edge[, 2])])
  expect_equal(nrow(get_allowed_segments_tracked(tracked$tracker, tracked$phy, clade_nodes, setdiff(tracker$states, expected_state), edge_ages_after)), 0)

  # The two pieces of the branch of the backbone keep their states
  info <- get_lineage_info(tracked$tracker, tracked$phy, tracked$phy$edge[, 2], edge_ages_after)
  expect_false(any(is.na(info$backbone_node) & is.na(info$state)))
})
