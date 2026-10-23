
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


# ---- Allowed ranges (states treated as biogeographic ranges) ----

## Test that resolve_allowed_ranges() resolves the allowed ranges into acceptable states and candidate ranges, with the global or taxon-level values
test_that("resolve_allowed_ranges() gives the acceptable states and the candidate ranges of each taxon, taxon-level values replacing the global one",
{
  states <- c("A", "AB", "BC")
  taxa <- c("X1", "X2", "X3", "X4")
  # Global "A" matches A and AB; the taxon-level values replace it (NULL and NA use the global value)
  resolved <- suppressWarnings(resolve_allowed_ranges(allowed_states_taxa = list(X1 = NULL, X2 = "BC", X3 = c("A", "B"), X4 = NA), allowed_states = "A",
                                                      states = states, taxon_names = taxa))
  expect_identical(names(resolved), c("allowed_states", "range_candidates", "inherited_states"))
  expect_identical(names(resolved$allowed_states), taxa)
  expect_identical(names(resolved$range_candidates), taxa)
  expect_identical(resolved$allowed_states$X1, c("A", "AB"))
  expect_identical(resolved$range_candidates$X1, list(A = "A", AB = "A"))
  expect_identical(resolved$allowed_states$X2, "BC")
  expect_identical(resolved$range_candidates$X2, list(BC = "BC"))
  # "B" is inherited in AB and BC but is not a state of the simmap: it is added to the acceptable states
  expect_identical(resolved$inherited_states, "B")
  expect_identical(resolved$allowed_states$X3, c("A", "AB", "BC", "B"))
  expect_identical(resolved$range_candidates$X3, list(A = "A", AB = c("A", "B"), BC = "B", B = "B"))
  expect_identical(resolved$allowed_states$X4, c("A", "AB"))
})

## Test that without allowed states all the states are acceptable and each state is its own candidate (the whole state is inherited), without warning
test_that("resolve_allowed_ranges() without allowed states allows all the states, each inherited as such, without warning",
{
  states <- c("A", "AB", "BC")
  expect_no_warning(resolved <- resolve_allowed_ranges(states = states, taxon_names = c("X1", "X2")))
  expect_identical(resolved$allowed_states, list(X1 = states, X2 = states))
  expect_identical(resolved$range_candidates$X1, list(A = "A", AB = "AB", BC = "BC"))
  expect_identical(resolved$range_candidates$X2, resolved$range_candidates$X1)
})

## Test that single-area states give the same acceptable states as the exact matching, without warning
test_that("resolve_allowed_ranges() with single-area states is the same as the exact matching of resolve_allowed_states(), without warning",
{
  states <- c("A", "B", "C")
  taxa <- c("X1", "X2", "X3")
  allowed_taxa <- list(X1 = c("C", "A"), X2 = NULL, X3 = "B")
  expect_no_warning(ranges <- resolve_allowed_ranges(allowed_states_taxa = allowed_taxa, allowed_states = "A", states = states, taxon_names = taxa))
  exact <- resolve_allowed_states(allowed_states_taxa = allowed_taxa, allowed_states = "A", states = states, taxon_names = taxa)
  # The acceptable states are in the order of the states of the simmap (the exact matching keeps the order of the allowed states)
  expect_identical(lapply(ranges$allowed_states, sort), lapply(exact, sort))
  expect_true(all(lengths(unlist(ranges$range_candidates, recursive = FALSE)) == 1))
  expect_identical(ranges$range_candidates$X1, list(A = "A", C = "C"))
})

## Test that no warning is sent when an allowed range is a state of the simmap and is smaller than other states containing it (e.g., "A" in "AB")
test_that("resolve_allowed_ranges() does not warn when a candidate range is smaller than the state",
{
  states <- c("A", "AB", "BC")
  expect_no_warning(resolved <- resolve_allowed_ranges(allowed_states = "A", states = states, taxon_names = "X1"))
  expect_identical(resolved$range_candidates$X1, list(A = "A", AB = "A"))
  expect_no_warning(resolve_allowed_ranges(allowed_states = c("AB", "BC"), states = states, taxon_names = "X1"))
  expect_no_warning(resolved <- resolve_allowed_ranges(allowed_states = "BC", states = c("A", "AB", "BC", "ABC"), taxon_names = "X1"))
  expect_identical(resolved$range_candidates$X1, list(BC = "BC", ABC = "BC"))
})

## Test the examples given in the documentation of resolve_allowed_ranges()
test_that("resolve_allowed_ranges() gives the outputs described in its documentation",
{
  states <- c("A", "B", "AB", "BC")
  resolve <- function (...) resolve_allowed_ranges(states = states, taxon_names = "X1", ...)
  expect_identical(resolve(allowed_states = "B")$allowed_states$X1, c("B", "AB", "BC"))
  expect_identical(resolve(allowed_states = "B")$range_candidates$X1, list(B = "B", AB = "B", BC = "B"))
  expect_identical(resolve(allowed_states = c("A", "B"))$allowed_states$X1, c("A", "B", "AB", "BC"))
  expect_identical(resolve(allowed_states = c("A", "B"))$range_candidates$X1, list(A = "A", B = "B", AB = c("A", "B"), BC = "B"))
  expect_identical(resolve(allowed_states = "AB")$allowed_states$X1, "AB")
  expect_identical(resolve(allowed_states = "AB")$range_candidates$X1, list(AB = "AB"))
  expect_identical(resolve()$allowed_states$X1, states)
  expect_identical(resolve()$range_candidates$X1, list(A = "A", B = "B", AB = "AB", BC = "BC"))
  expect_error(resolve(allowed_states = "D"), "Areas available in the simmap: 'A', 'B', 'C'")
})

## Test that the ranges inherited that are not states of the simmap are added to the acceptable states of all the taxa
test_that("resolve_allowed_ranges() adds the inherited ranges that are not states of the simmap, also for the taxa without allowed states",
{
  # Example of the documentation: allowed range "B" with the states AB and BC
  expect_warning(resolved <- resolve_allowed_ranges(allowed_states_taxa = list(X1 = "B", X2 = NULL, X3 = "BC"), states = c("AB", "BC"), taxon_names = c("X1", "X2", "X3")),
                 "not states of the simmap: 'B'")
  expect_identical(resolved$inherited_states, "B")
  expect_identical(resolved$allowed_states$X1, c("AB", "BC", "B"))
  expect_identical(resolved$range_candidates$X1, list(AB = "B", BC = "B", B = "B"))
  # A taxon without allowed states can be grafted on the branches of the lineages that inherited "B", and inherits the whole state
  expect_identical(resolved$allowed_states$X2, c("AB", "BC", "B"))
  expect_identical(resolved$range_candidates$X2, list(AB = "AB", BC = "BC", B = "B"))
  # "BC" is a state of the simmap: it matches BC only here, and not "B" (not contained in "BC")
  expect_identical(resolved$allowed_states$X3, "BC")
  expect_identical(resolved$range_candidates$X3, list(BC = "BC"))

  # No inherited range outside of the states of the simmap
  expect_identical(resolve_allowed_ranges(allowed_states = "BC", states = c("AB", "BC"), taxon_names = "X1")$inherited_states, character(0))
  expect_identical(resolve_allowed_ranges(states = c("AB", "BC"), taxon_names = "X1")$inherited_states, character(0))
})

## Test that an allowed range that is not a state of the simmap but is contained in some of them is accepted, with a warning explaining the mismatch
test_that("resolve_allowed_ranges() warns about the allowed ranges that are not states of the simmap but are contained in some of them, and still resolves them",
{
  states <- c("A", "AB", "BC")
  # "B" is contained in AB and BC but is not a state of the simmap: no stop, and a warning explaining the mismatch
  expect_warning(resolved <- resolve_allowed_ranges(allowed_states = c("A", "B"), states = states, taxon_names = c("X1", "X2")),
                 "some allowed ranges are not states of the simmap: 'B' \\(contained in the states 'AB', 'BC'; taxa: X1, X2\\)\\.\nThe grafting is run assuming that these ranges can be inherited.*subset sympatric speciation events \\(e\\.g\\., \\(AB\\) -> \\(AB\\), \\(B\\) for the range 'B' inherited from the state 'AB'\\).*added to the states of the updated simmap")
  expect_identical(resolved$inherited_states, "B") # inherited, and added to the acceptable states
  expect_identical(resolved$allowed_states$X1, c("A", "AB", "BC", "B"))
  expect_identical(resolved$range_candidates$X1, list(A = "A", AB = c("A", "B"), BC = "B", B = "B"))

  # Taxon-level values: only the taxa concerned are listed. Several ranges are described in one warning
  warnings <- testthat::capture_warnings(resolve_allowed_ranges(allowed_states_taxa = list(X1 = "B", X2 = "A", X3 = c("B", "C")), states = c("A", "AB", "BC", "C"), taxon_names = c("X1", "X2", "X3")))
  expect_length(warnings, 1)
  expect_match(warnings, "'B' \\(contained in the states 'AB', 'BC'; taxa: X1, X3\\)")
  expect_no_match(warnings, "'A' \\(contained")
  expect_no_match(warnings, "'C' \\(contained")

  # No warning when the allowed ranges are states of the simmap, or without allowed states
  expect_no_warning(resolve_allowed_ranges(allowed_states = c("A", "AB", "BC"), states = states, taxon_names = "X1"))
  expect_no_warning(resolve_allowed_ranges(states = states, taxon_names = "X1"))
  # A range that is not contained in any state is not described in this warning (it is ignored with its own warning)
  warnings <- testthat::capture_warnings(resolve_allowed_ranges(allowed_states = c("A", "AC"), states = c("A", "AB", "BC"), taxon_names = "X1"))
  expect_length(warnings, 1)
  expect_match(warnings, "not contained in any state of the simmap, so they are ignored: 'AC'")
})

## Test that allowed ranges contained in no state are ignored with a warning, and that the function stops if none is contained in a state
test_that("resolve_allowed_ranges() warns about allowed ranges contained in no state, and stops if a taxon has no state containing its allowed ranges",
{
  states <- c("A", "B", "AB", "BC")
  # "AC" is made of areas of the simmap but no state contains it: X2 has no acceptable state
  message <- tryCatch(suppressWarnings(resolve_allowed_ranges(allowed_states_taxa = list(X1 = c("A", "AC"), X2 = "AC", X3 = NULL), allowed_states = "B",
                                                              states = states, taxon_names = c("X1", "X2", "X3"))),
                      error = function (e) conditionMessage(e))
  expect_match(message, "No state of the simmap contains the allowed range\\(s\\) of taxon 'X2': 'AC'.*States available in the simmap: 'A', 'B', 'AB', 'BC'")

  # Only some ranges are contained in no state: they are ignored
  warnings <- testthat::capture_warnings(resolved <- resolve_allowed_ranges(allowed_states_taxa = list(X1 = c("A", "AC")), states = states, taxon_names = "X1"))
  expect_true(any(grepl("not contained in any state of the simmap, so they are ignored: 'AC' \\(taxa: X1\\)", warnings)))
  expect_identical(resolved$allowed_states$X1, c("A", "AB"))
  expect_identical(resolved$range_candidates$X1, list(A = "A", AB = "A"))
})

## Test that the states of the simmap and the allowed ranges are strictly checked, naming the incorrect states
test_that("resolve_allowed_ranges() stops if the states or the allowed ranges are not coded as in BioGeoBEARS, or if an allowed area is not in the simmap",
{
  states <- c("A", "AB", "BC")
  # States of the simmap
  expect_error(resolve_allowed_ranges(states = c("A", "red", "BA"), taxon_names = "X1"),
               "the states in the simmap must be biogeographic ranges.*'red' \\(lowercase letters \\(r, e, d\\)\\); 'BA' \\(areas not in alphabetical order \\(expected 'AB'\\)\\)")
  # Global and taxon-level allowed ranges
  expect_error(resolve_allowed_ranges(allowed_states = c("A", "ba"), states = states, taxon_names = "X1"),
               "the global 'allowed_states' argument.*'ba' \\(lowercase letters")
  expect_error(resolve_allowed_ranges(allowed_states_taxa = list(X1 = "BA"), states = states, taxon_names = "X1"),
               "'missing_taxa_list\\$allowed_states' of taxon 'X1'.*'BA' \\(areas not in alphabetical order")
  # Area not in the simmap
  expect_error(resolve_allowed_ranges(allowed_states = c("A", "D"), states = states, taxon_names = "X1"),
               "the global 'allowed_states' argument are not in the simmap: 'D'.*Areas available in the simmap: 'A', 'B', 'C'")
  # Same checks of the list as with the exact matching
  expect_error(resolve_allowed_ranges(allowed_states_taxa = list(X1 = "A"), states = states, taxon_names = c("X1", "X2")), "same length as 'Taxon_name'")
  expect_error(resolve_allowed_ranges(allowed_states_taxa = c(X1 = "A"), states = states, taxon_names = "X1"), "must be a list")
})


# ---- Inheritance of ranges at graft ----

## Test that the inherited state is the state of the point of the graft without candidates, and a candidate range otherwise, without random draw for a single candidate
test_that("choose_inherited_state() returns the state without candidates, the only candidate without drawing, and draws uniformly among several candidates",
{
  expect_identical(choose_inherited_state("AB"), "AB")
  candidates <- list(A = "A", AB = c("A", "B"), BC = "B")
  expect_identical(choose_inherited_state("A", candidates), "A")

  # No random number is drawn with a single candidate
  set.seed(1)
  seed_before <- .Random.seed
  choose_inherited_state("BC", candidates)
  expect_identical(.Random.seed, seed_before)

  # Both candidates are drawn, with a similar frequency
  set.seed(2)
  draws <- vapply(1:400, function (i) choose_inherited_state("AB", candidates), character(1))
  expect_setequal(unique(draws), c("A", "B"))
  expect_true(all(table(draws) > 140))
  expect_error(choose_inherited_state("C", candidates), "no range can be inherited when grafting on a branch in state 'C'")
})

## Test that tracker states are extended after the states of the simmap, without changing their index
test_that("extend_tracker_states() adds the new states after the states of the simmap",
{
  s <- toy_simmap_setup()
  tracker <- start_state_tracking(label_backbone_nodes(s$phy)$phy, s$simmap)
  extended <- extend_tracker_states(tracker, c("AB", "A", "B"))
  expect_identical(extended$states, c("A", "B", "C", "AB"))
  expect_identical(extended$extra_states, "AB")
  expect_identical(tracker$extra_states, character(0))
  expect_identical(extend_tracker_states(extended, c("AB", "A"))$states, extended$states)
})

## Run two successive grafts on the toy simmap with ranges (A, AB, BC): the first one on a piece in AB, the second one on the branch of the first one
graft_with_ranges <- function (candidates, extra_states)
{
  s <- toy_simmap_setup()
  simmap <- relabel_simmap(s$simmap, c(A = "A", B = "AB", C = "BC"))
  phy <- label_backbone_nodes(s$phy)$phy
  tracker <- extend_tracker_states(start_state_tracking(phy, simmap), extra_states)
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  acceptable <- names(candidates)

  graft <- function (phy, tracker, node, position, label)
  {
    edge_ages <- root_age - phytools::nodeHeights(phy)
    new_length <- position + (root_age - phytools::nodeheight(phy, node))
    phy_after <- phytools::bind.tip(phy, tip.label = label, edge.length = new_length, where = node, position = position)
    track_graft(tracker, phy, phy_after, node, position, label, allowed_states = acceptable, range_candidates = candidates, edge_ages = edge_ages)
  }

  # First graft: in the middle of a piece in AB
  pieces <- tracker$state_table[tracker$state_table$state == "AB", ]
  first <- graft(phy, tracker, pieces$node[1], mean(c(pieces$from[1], pieces$to[1])), "T1")
  # Second graft: in the middle of the branch of the first taxon
  node_2 <- match("T1", first$phy$tip.label)
  second <- graft(first$phy, first$tracker, node_2, first$phy$edge.length[match(node_2, first$phy$edge[, 2])] / 2, "T2")
  list(first = first, second = second, simmap = simmap, root_age = root_age)
}

## Test that a grafted taxon inherits a candidate range, and the lineages of the grafted taxa keep this state
test_that("track_graft() records the inherited range, and the new node on the branch of a grafted taxon keeps the state of this branch",
{
  # In AB, the taxon inherits "A" (it is a state of the simmap)
  res <- graft_with_ranges(list(A = "A", AB = "A", BC = "B", B = "B"), character(0))
  expect_identical(unname(res$first$tracker$missing_state["T1"]), "A")
  # The second taxon is grafted on the branch of the first one, in state A: it inherits the candidate of this state, and the new node has the state of the branch
  expect_identical(unname(res$second$tracker$missing_state["T2"]), "A")
  states <- res$second$tracker$states
  node_labels <- res$second$phy$node.label
  expect_true(paste0("ms", match("A", states)) %in% node_labels[!node_labels %in% res$first$phy$node.label])

  # The returned simmap: the grafted branches are entirely in the inherited range, the branch of T1 is cut in two pieces in this range
  simmap <- build_updated_simmap(res$second$tracker, res$second$phy, res$second$phy)
  edge <- simmap$edge
  tip_states <- vapply(c("T1", "T2"), function (taxon) paste(names(simmap$maps[[match(match(taxon, simmap$tip.label), edge[, 2])]]), collapse = "+"), character(1))
  expect_identical(unname(tip_states), c("A", "A"))
  node_T1 <- res$second$phy$edge[match(match("T1", res$second$phy$tip.label), res$second$phy$edge[, 2]), 1]
  expect_identical(names(simmap$maps[[match(node_T1, edge[, 2])]]), "A") # branch above the node where T2 is attached: still the branch of T1
  expect_identical(colnames(simmap$mapped.edge), c("A", "AB", "BC")) # states of the simmap
  expect_equal(unname(rowSums(simmap$mapped.edge)), simmap$edge.length, tolerance = 1e-8)
})

## Test that a range inherited that is not a state of the simmap is a state of the tracker and of the returned simmap, and that taxa can be grafted on its branches
test_that("track_graft() records an inherited range that is not a state of the simmap, which becomes a state of the grafted lineages and of the returned simmap",
{
  res <- graft_with_ranges(list(A = "A", AB = "B", BC = "B", B = "B"), "B")
  expect_identical(res$second$tracker$states, c("A", "AB", "BC", "B"))
  expect_identical(unname(res$first$tracker$missing_state["T1"]), "B")
  expect_identical(unname(res$second$tracker$missing_state["T2"]), "B")

  # The branches of the grafted taxa are in the allowed state "B": the segments in this state are the branches of the grafted taxa (and the new node on the branch of T1)
  root_age <- res$root_age
  phy <- res$second$phy
  segments <- get_allowed_segments_tracked(res$second$tracker, phy, phy$edge[, 2], "B", root_age - phytools::nodeHeights(phy))
  expect_true(all(match(c("T1", "T2"), phy$tip.label) %in% segments$node))
  expect_equal(sum(segments$to - segments$from), sum(phy$edge.length[match(segments$node, phy$edge[, 2])]))
  expect_equal(nrow(get_allowed_segments_tracked(res$first$tracker, res$first$phy, res$first$phy$edge[, 2], "B", root_age - phytools::nodeHeights(res$first$phy))), 1)

  # The returned simmap has the column "B", with the total length of the branches of T1 and T2 and of the pieces of the branch of T1
  simmap <- build_updated_simmap(res$second$tracker, phy, phy)
  expect_identical(colnames(simmap$mapped.edge), c("A", "AB", "BC", "B"))
  expect_true(simmap$mapped.edge[match(match("T1", phy$tip.label), phy$edge[, 2]), "B"] > 0)

  # An added state that no lineage inherited is not in the returned simmap
  tracker <- extend_tracker_states(res$second$tracker, "D")
  expect_false("D" %in% colnames(build_updated_simmap(tracker, phy, phy)$mapped.edge))
  expect_true(all(c("A", "AB", "BC") %in% colnames(build_updated_simmap(res$second$tracker, phy, phy)$mapped.edge)))
})


# ---- Simmap without effect ----

## Test that the warning is only sent if no allowed states are provided anywhere and the simmap is not returned
test_that("warn_if_simmap_unused() warns only without any allowed states and with return_simmap = FALSE",
{
  expect_warning(warn_if_simmap_unused(return_simmap = FALSE), "no effect.*'return_simmap = TRUE'")
  expect_warning(warn_if_simmap_unused(allowed_states_taxa = list(X1 = NULL, X2 = NA, X3 = character(0)), return_simmap = FALSE), "no effect")
  # The simmap is returned: it only follows the states, which is a documented use
  expect_no_warning(warn_if_simmap_unused(return_simmap = TRUE))
  expect_no_warning(warn_if_simmap_unused(allowed_states_taxa = list(X1 = NULL, X2 = NA), return_simmap = TRUE))
  # Allowed states, global or taxon-level
  expect_no_warning(warn_if_simmap_unused(allowed_states = "A", return_simmap = FALSE))
  expect_no_warning(warn_if_simmap_unused(allowed_states_taxa = list(X1 = NULL, X2 = "B"), return_simmap = FALSE))
  expect_null(suppressWarnings(warn_if_simmap_unused(return_simmap = FALSE)))
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
