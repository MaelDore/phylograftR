
##### Test suite for check_grafting_plan() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Name the items of a list of missing taxa after the taxa, as required by the functions
plan_mtl <- function (l)
{
  if (!is.list(l) || is.null(l$Taxon_name)) { return(l) }
  for (item in setdiff(names(l), "Taxon_name")) { if (is.null(names(l[[item]])) && length(l[[item]]) == length(l$Taxon_name)) { names(l[[item]]) <- l$Taxon_name } }
  l
}
plan_check <- function (phy, missing_taxa_list, ...) { check_grafting_plan(phy, plan_mtl(missing_taxa_list), ...) }

plan_multi_simmap <- function (maps) { structure(maps, class = c("multiSimmap", "multiPhylo")) }

## Age of the rootward end (stem age) and tipward end (crown age) of the branch leading to a node
plan_stem_age <- function (phy, node) { max(phytools::nodeHeights(phy)[, 2]) - phytools::nodeHeights(phy)[match(node, phy$edge[, 2]), 1] }
plan_crown_age <- function (phy, node) { max(phytools::nodeHeights(phy)[, 2]) - phytools::nodeheight(phy, node) }


test_that("check_grafting_plan() reports the branches, the time window and the enclosing clade of each taxon, without simmap",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("xA", "xB", "xs"), Target_clade = list(toy$A, toy$B, toy$s))
  plan <- plan_check(toy$phy, mtl, min_tips = 3, verbose = FALSE)

  expect_s3_class(plan, "data.frame")
  expect_equal(plan$Taxon_name, mtl$Taxon_name)
  expect_equal(plan$target_n_tips, c(toy$nA, toy$nB, 1))
  # Pool of branches: the descendant nodes of the target clade, and its stem
  expect_equal(plan$n_branches, c(length(get_descendant_nodes(toy$phy, toy$node_A)) + 1, length(get_descendant_nodes(toy$phy, toy$node_B)) + 1, 1))
  # Time window: up to the stem age of the target clade
  s_node <- match(toy$s, toy$phy$tip.label)
  expect_equal(plan$window_age, c(plan_stem_age(toy$phy, toy$node_A), plan_stem_age(toy$phy, toy$node_B), plan_stem_age(toy$phy, s_node)))
  # With min_tips = 3, the clades of 3 tips or more are their own enclosing clade, and the tip uses its parent clade
  expect_equal(plan$enclosing_n_tips[1:2], c(toy$nA, toy$nB))
  expect_equal(plan$enclosing_age[1:2], c(plan_crown_age(toy$phy, toy$node_A), plan_crown_age(toy$phy, toy$node_B)))
  expect_gte(plan$enclosing_n_tips[3], 3)
  expect_true(all(plan$feasible))
  expect_equal(plan$cause, rep("none", 3))
  # No state columns without simmap
  expect_false(any(c("valid_fraction", "time_coverage", "allowed_states", "states_available") %in% names(plan)))
  expect_equal(attr(plan, "type"), "taxa")
  expect_equal(attr(plan, "n_maps"), 0)

  # The enclosing clade is the smallest one with at least 'min_tips' tips
  big <- plan_check(toy$phy, mtl, min_tips = toy$nA + 1, verbose = FALSE)
  expect_true(all(big$enclosing_n_tips >= toy$nA + 1))
  expect_true(all(big$enclosing_n_tips >= big$target_n_tips))
  parent_A <- get_parent_node(toy$phy, toy$node_A)
  expect_equal(big$enclosing_n_tips[1], ape::Ntip(ape::extract.clade(toy$phy, parent_A)))
})

test_that("check_grafting_plan() accounts for the stem availability and the excluded subclades, and flags the ones that stop the grafting functions",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("with_stem", "no_stem", "excl_B", "excl_contains", "excl_same"),
              Target_clade = list(toy$A, toy$A, toy$A, toy$B, toy$B),
              include_target_clade_stem = c(TRUE, FALSE, TRUE, TRUE, TRUE),
              Excluded_subclades = list(NULL, NULL, list(toy$B), list(toy$A), list(toy$B)))
  names(mtl$include_target_clade_stem) <- mtl$Taxon_name
  plan <- plan_check(toy$phy, mtl, verbose = FALSE)
  n_A <- length(get_descendant_nodes(toy$phy, toy$node_A))
  n_B <- length(get_descendant_nodes(toy$phy, toy$node_B))
  expect_equal(plan$n_branches[1:3], c(n_A + 1, n_A, n_A + 1 - n_B - 1))
  # Without the stem, the window ends at the crown age
  expect_equal(plan$window_age[2], plan_crown_age(toy$phy, toy$node_A))
  expect_equal(plan$cause, c("none", "none", "none", "excluded_contains_target", "excluded_identical"))
  expect_equal(plan$feasible, c(TRUE, TRUE, TRUE, FALSE, FALSE))

  # The grafting function stops on the same taxa
  expect_error(graft_missing_taxa(toy$phy, plan_mtl(mtl), verbose = FALSE), "excluded subclade")
  # Without the excluded subclades, everything can be grafted
  expect_true(all(plan_check(toy$phy, mtl, use_excluded_subclades = FALSE, verbose = FALSE)$feasible))

  # Global stem availability
  plan_no_stem <- plan_check(toy$phy, list(Taxon_name = "x", Target_clade = list(toy$A)), include_target_clade_stem = FALSE, verbose = FALSE)
  expect_equal(plan_no_stem$n_branches, n_A)
})

test_that("check_grafting_plan() reports the share of branches and of time in allowed states, and the cause of the taxa that cannot be grafted, with a simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- list(Taxon_name = c("in_B", "not_in_B", "tip_A", "tip_B"), Target_clade = list(toy$B, toy$B, toy$s, toy$s),
              allowed_states = list(in_B = "B", not_in_B = c("A", "C"), tip_A = "A", tip_B = "B"))
  plan <- plan_check(toy$phy, mtl, simmap = simmap, verbose = FALSE)

  # Same values as the upfront check of the grafting functions
  phy <- label_backbone_nodes(toy$phy)$phy
  tracker <- start_state_tracking(phy, simmap)
  flagged <- find_taxa_without_allowed_branch(phy, mtl, tracker = tracker, allowed_states = mtl$allowed_states)
  expect_equal(plan$valid_fraction, unname(attr(flagged, "valid_fraction")))
  expect_equal(plan$time_coverage, unname(attr(flagged, "time_coverage")))
  expect_equal(plan$feasible, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(plan$cause, c("none", "allowed_states", "none", "allowed_states"))
  expect_equal(plan$allowed_states, c("B", "A, C", "A", "B"))
  # The states of the pool: the clade B is entirely in "B", the tip in "A"
  expect_equal(plan$states_available, c("B", "B", "A", "A"))

  # The grafting function stops with the same taxa
  expect_error(graft_missing_taxa(toy$phy, plan_mtl(mtl), simmap = simmap, verbose = FALSE), "no branch available for grafting in their allowed states")
  ok <- list(Taxon_name = mtl$Taxon_name[c(1, 3)], Target_clade = mtl$Target_clade[c(1, 3)], allowed_states = mtl$allowed_states[c(1, 3)])
  expect_no_error(suppressWarnings(graft_missing_taxa(toy$phy, plan_mtl(ok), simmap = simmap, verbose = FALSE)))

  # Global allowed states: used for the taxa without their own
  plan_global <- plan_check(toy$phy, list(Taxon_name = c("a", "b"), Target_clade = list(toy$B, toy$B)), simmap = simmap, allowed_states = "C", verbose = FALSE)
  expect_equal(plan_global$allowed_states, c("C", "C"))
  expect_false(any(plan_global$feasible))
})

test_that("check_grafting_plan() flags the taxa with a low share of branches in allowed states, and the time window covered",
{
  toy <- toy_setup()
  # The state "B" is found below the age 4 only
  simmap <- make_young_state_simmap(toy$phy, 4)
  mtl <- list(Taxon_name = c("with_stem", "crown_only"), Target_clade = list(toy$A, toy$A), include_target_clade_stem = c(TRUE, FALSE))
  names(mtl$include_target_clade_stem) <- mtl$Taxon_name
  plan <- plan_check(toy$phy, mtl, simmap = simmap, allowed_states = "B", min_valid_fraction = 0.5, verbose = FALSE)
  expect_equal(plan$time_coverage, c(4 / plan_stem_age(toy$phy, toy$node_A), 4 / plan_crown_age(toy$phy, toy$node_A)))
  expect_true(all(plan$valid_fraction < 0.5))
  expect_equal(plan$cause, c("low_valid_fraction", "low_valid_fraction"))
  expect_true(all(plan$feasible)) # A low share is a warning, not an error
  # Below the threshold only if it is high enough
  expect_equal(plan_check(toy$phy, mtl, simmap = simmap, allowed_states = "B", min_valid_fraction = 0, verbose = FALSE)$cause, c("none", "none"))
})

test_that("check_grafting_plan() summarizes the maps of a multiSimmap, with the table of each map in an attribute",
{
  toy <- toy_setup()
  maps <- lapply(1:3, function (i) make_toy_simmap(toy$phy, seed = i))
  mtl <- list(Taxon_name = c("xA", "xB", "xs"), Target_clade = list(toy$A, toy$B, toy$s))
  plan <- plan_check(toy$phy, mtl, simmap = plan_multi_simmap(maps), allowed_states = "A", verbose = FALSE)
  per_map <- attr(plan, "per_map")

  expect_equal(attr(plan, "n_maps"), 3)
  expect_equal(nrow(plan), 3)
  expect_equal(nrow(per_map), 9)
  expect_equal(unique(per_map$map), c("map_1", "map_2", "map_3"))
  # Each map gives the same table as the check of the map alone
  single <- lapply(maps, function (m) plan_check(toy$phy, mtl, simmap = m, allowed_states = "A", verbose = FALSE))
  expect_equal(per_map$valid_fraction, unlist(lapply(single, `[[`, "valid_fraction")))
  expect_equal(per_map$feasible, unlist(lapply(single, `[[`, "feasible")))
  # Summary over the maps
  vf <- matrix(per_map$valid_fraction, nrow = 3)
  expect_equal(plan$valid_fraction_min, apply(vf, 1, min))
  expect_equal(plan$valid_fraction_mean, rowMeans(vf))
  expect_equal(plan$valid_fraction_max, apply(vf, 1, max))
  feasible <- matrix(per_map$feasible, nrow = 3)
  expect_equal(plan$n_maps_infeasible, 3L - as.integer(rowSums(feasible)))
  expect_equal(plan$feasible, plan$n_maps_infeasible == 0)
  # Named maps
  named <- plan_multi_simmap(stats::setNames(maps, c("a", "b", "c")))
  expect_equal(unique(attr(plan_check(toy$phy, mtl, simmap = named, verbose = FALSE), "per_map")$map), c("a", "b", "c"))
})

test_that("check_grafting_plan() with a 'Sister_clade' list checks the missing clades: stem, allowed states and consistency",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- list(Taxon_name = c("M1", "M2", "M3"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$s))
  plan <- plan_check(toy$phy, mtl, verbose = FALSE)
  # One row per missing clade (M1 and M2 are grafted together)
  expect_equal(nrow(plan), 2)
  expect_equal(plan$Taxa, c("M1, M2", "M3"))
  expect_equal(plan$sister_n_tips, c(toy$nB, 1))
  expect_equal(plan$stem_length, toy$phy$edge.length[match(c(toy$node_B, match(toy$s, toy$phy$tip.label)), toy$phy$edge[, 2])])
  expect_equal(plan$stem_age, c(plan_stem_age(toy$phy, toy$node_B), plan_stem_age(toy$phy, match(toy$s, toy$phy$tip.label))))
  expect_equal(attr(plan, "type"), "clades")
  expect_true(all(plan$feasible))

  # With a simmap: the clade B is entirely in "B", and the tip 's' in "A"
  mtl$allowed_states <- list(M1 = "B", M2 = "B", M3 = "B")
  plan_states <- plan_check(toy$phy, mtl, simmap = simmap, verbose = FALSE)
  expect_equal(plan_states$valid_fraction, c(1, 0))
  expect_equal(plan_states$time_coverage, plan_states$valid_fraction)
  expect_equal(plan_states$cause, c("none", "allowed_states"))
  expect_equal(plan_states$feasible, c(TRUE, FALSE))
  expect_error(graft_missing_clades(toy$phy, plan_mtl(mtl), simmap = simmap, verbose = FALSE), "no piece of the stem branch")

  # The taxa of a missing clade must share their allowed states
  mtl$allowed_states <- list(M1 = "B", M2 = "A", M3 = "A")
  plan_diff <- plan_check(toy$phy, mtl, simmap = simmap, verbose = FALSE)
  expect_equal(plan_diff$cause, c("inconsistent_allowed_states", "none"))
  expect_equal(plan_diff$feasible, c(FALSE, TRUE))
  expect_error(graft_missing_clades(toy$phy, plan_mtl(mtl), simmap = simmap, verbose = FALSE), "different allowed states")

  # A sister clade that is the whole tree has no stem
  whole <- plan_check(toy$phy, list(Taxon_name = "W", Sister_clade = list(W = toy$phy$tip.label)), verbose = FALSE)
  expect_equal(whole$cause, "whole_tree")
  expect_false(whole$feasible)

  # Several maps
  maps <- lapply(1:3, function (i) make_toy_simmap(toy$phy, seed = i))
  plan_multi <- plan_check(toy$phy, list(Taxon_name = c("M1", "M3"), Sister_clade = list(M1 = toy$B, M3 = toy$s)), simmap = plan_multi_simmap(maps), allowed_states = "A", verbose = FALSE)
  expect_equal(nrow(attr(plan_multi, "per_map")), 6)
  expect_true(all(c("valid_fraction_min", "n_maps_infeasible") %in% names(plan_multi)))
})

test_that("check_grafting_plan() works with ranges and several traits",
{
  toy <- toy_setup()
  ranges <- make_toy_range_simmap(toy$phy)
  mtl <- list(Taxon_name = c("x1", "x2"), Target_clade = list(toy$A, toy$B), allowed_states = list(x1 = "A", x2 = "B"))
  plan <- suppressWarnings(plan_check(toy$phy, mtl, simmap = ranges, treat_states_as_ranges = TRUE, verbose = FALSE))
  # The allowed ranges are resolved into the states that contain them (the states of the simmap, and the ranges that can be inherited from them)
  expect_true(all(c("A", "AB") %in% strsplit(plan$allowed_states[1], ", ")[[1]]))
  expect_true(all(c("AB", "BC") %in% strsplit(plan$allowed_states[2], ", ")[[1]]))
  expect_false("A" %in% strsplit(plan$allowed_states[2], ", ")[[1]])
  # Same as the upfront check of the grafting function
  expect_equal(plan$feasible, c(TRUE, TRUE) & plan$valid_fraction > 0)

  # Several traits: one simmap per trait, composite states
  traits <- list(biogeo = ranges, elev = make_toy_simmap(toy$phy, seed = 2))
  mtl_traits <- list(Taxon_name = "x", Target_clade = list(toy$A))
  plan_traits <- suppressWarnings(plan_check(toy$phy, mtl_traits, simmap = traits, allowed_states = list(c("A", "AB"), NULL), treat_states_as_ranges = c(biogeo = TRUE, elev = FALSE), verbose = FALSE))
  expect_equal(nrow(plan_traits), 1)
  expect_true(is.numeric(plan_traits$valid_fraction))
})

test_that("check_grafting_plan() stops on invalid inputs",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = "x", Target_clade = list(toy$A))
  expect_error(plan_check(toy$phy, list(Taxon_name = "x"), verbose = FALSE), "'Target_clade' item")
  expect_error(plan_check(toy$phy, list(Taxon_name = "x", Target_clade = list(toy$A), Sister_clade = list(toy$A)), verbose = FALSE), "both 'Target_clade' and 'Sister_clade'")
  expect_error(plan_check(toy$phy, mtl, min_tips = 2, verbose = FALSE), "'min_tips'")
  expect_error(plan_check(toy$phy, mtl, min_valid_fraction = 2, verbose = FALSE), "'min_valid_fraction'")
  expect_error(plan_check(toy$phy, mtl, allowed_states = "A", verbose = FALSE), "requires a 'simmap'")
  expect_error(plan_check(toy$phy, mtl, treat_states_as_ranges = TRUE, verbose = FALSE), "requires a 'simmap'")
  expect_error(plan_check(toy$phy, list(Taxon_name = toy$phy$tip.label[1], Target_clade = list(toy$A)), verbose = FALSE), "absent from the backbone")
  expect_error(plan_check(toy$phy, list(Taxon_name = "x", Target_clade = list(c("nope", "t1"))), verbose = FALSE), "missing from the backbone")
  expect_error(check_grafting_plan(make_toy_polytomies(toy$phy), mtl, verbose = FALSE), "binary")
  expect_error(plan_check(toy$phy, "a", verbose = FALSE), "must be a list")
  # Map names appear in the errors
  maps <- plan_multi_simmap(list(make_toy_simmap(toy$phy), make_toy_simmap(make_toy_phy(n_tips = 20, seed = 5))))
  expect_error(plan_check(toy$phy, mtl, simmap = maps, verbose = FALSE), "Map 'map_2'")
})

test_that("check_grafting_plan() prints a summary, and nothing with verbose = FALSE",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- list(Taxon_name = c("in_B", "tip_B"), Target_clade = list(toy$B, toy$s), allowed_states = list(in_B = "B", tip_B = "B"))
  expect_output(plan_check(toy$phy, mtl, simmap = simmap), "Grafting plan: 2 taxa")
  expect_output(plan_check(toy$phy, mtl, simmap = simmap), "Cannot be grafted \\(1\\)")
  expect_output(plan_check(toy$phy, mtl, simmap = simmap), "tip_B")
  expect_silent(plan_check(toy$phy, mtl, simmap = simmap, verbose = FALSE))
  expect_output(plan_check(toy$phy, list(Taxon_name = "M", Sister_clade = list(M = toy$B))), "1 missing clades")
})
