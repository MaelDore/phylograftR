
##### Test suite for check_grafting_constraints() with missing clades (output of graft_missing_clades(), 'Sister_clade' lists) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

## Toy list of missing clades: M1 and M2 share a sister clade (one missing clade of 2 taxa), M3 and M4 have their own sister clade (single taxa, M4 on a tip branch)
cc_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}
cc_graft <- function (toy, mcl = cc_clades(toy), ...) suppressWarnings(graft_missing_clades(toy$phy, mcl, verbose = FALSE, seed = 4, ...))
cc_check <- function (grafted, toy, mcl = cc_clades(toy), ...) suppressWarnings(check_grafting_constraints(grafted, toy$phy, mcl, verbose = FALSE, ...))

## State names of the taxa of a clade in a simmap: replace the state of a branch by another one (to build violations)
cc_set_edge_state <- function (grafted, taxon, state, trait = NULL)
{
  edge <- match(match(taxon, grafted$tip.label), grafted$edge[, 2])
  sm <- if (is.null(trait)) grafted$simmap else grafted$simmap[[trait]]
  sm$maps[[edge]] <- stats::setNames(sum(sm$maps[[edge]]), state)
  if (is.null(trait)) { grafted$simmap <- sm } else { grafted$simmap[[trait]] <- sm }
  grafted
}


test_that("A missing clades list is checked, with one row per missing clade (the taxa sharing a sister clade), in the order of first appearance", {
  toy <- toy_setup()
  grafted <- cc_graft(toy)
  res <- cc_check(grafted, toy)
  expect_equal(res$missing_clade, 1:3)
  expect_equal(res$Taxon_names, c("M1, M2", "M3", "M4"))
  expect_equal(res$n_taxa, c(2, 1, 1))
  expect_equal(res$n_in_tree, c(2, 1, 1))
  expect_equal(res$Sister_clade_tips, c(paste(toy$B, collapse = ", "), paste(toy$C, collapse = ", "), toy$s))
  expect_equal(names(res), c("map", "simulation", "seed", "missing_clade", "Taxon_names", "n_taxa", "Sister_clade_tips", "in_tree", "n_in_tree", "monophyletic",
                             "sister_position", "sister_ok", "states_ok", "ranges_ok", "ok"))
  expect_true(all(res$in_tree) && all(res$monophyletic) && all(res$sister_ok) && all(res$ok))
  expect_equal(res$sister_position, c("stem", "stem", "tip_branch"))
  expect_true(all(is.na(res$states_ok)) && all(is.na(res$ranges_ok)))
  expect_true(attr(res, "all_ok"))
  expect_true(attr(res, "backbone_preserved"))
  expect_equal(attr(res, "phylogenies")$n_failed, 0)
})

test_that("The arguments of the taxa (stem, excluded subclades) are ignored with missing clades, without warning", {
  toy <- toy_setup()
  grafted <- cc_graft(toy)
  expect_no_warning(res <- check_grafting_constraints(grafted, toy$phy, cc_clades(toy), include_target_clade_stem = FALSE, use_excluded_subclades = TRUE, exclude_subclades_stem = FALSE, verbose = FALSE))
  expect_true(attr(res, "all_ok"))
})

test_that("A wrong sister clade is detected: position 'outside', and 'crown' if the clade is nested in the sister clade given", {
  toy <- toy_setup()
  grafted <- cc_graft(toy)
  # M3 is on the stem of C, checked against the sister clade D (disjoint from C): outside
  wrong <- cc_clades(toy); wrong$Sister_clade$M3 <- toy$D
  res <- cc_check(grafted, toy, wrong)
  expect_equal(res$sister_position[res$Taxon_names == "M3"], "outside")
  expect_false(res$sister_ok[res$Taxon_names == "M3"])
  expect_false(res$ok[res$Taxon_names == "M3"])
  expect_true(all(res$ok[res$Taxon_names != "M3"]))
  expect_false(attr(res, "all_ok"))
  # Checked against the whole clade A (that contains B and C): the missing clades are nested in it (crown), not sister to it
  nested <- cc_clades(toy); nested$Sister_clade$M3 <- toy$A
  res <- cc_check(grafted, toy, nested)
  expect_equal(res$sister_position[res$Taxon_names == "M3"], "crown")
  expect_false(res$sister_ok[res$Taxon_names == "M3"])
})

test_that("Taxa that do not form one clade are detected (monophyletic = FALSE), also when they are all on the stem", {
  toy <- toy_setup()
  # M1 and M2 grafted as two separate missing clades (different sister clades), checked as one missing clade of the same sister clade
  separate <- cc_clades(toy); separate$Sister_clade$M2 <- toy$C
  grafted <- cc_graft(toy, separate)
  res <- cc_check(grafted, toy, cc_clades(toy))
  expect_false(res$monophyletic[res$Taxon_names == "M1, M2"])
  expect_false(res$ok[res$Taxon_names == "M1, M2"])
  # Two taxa grafted separately on the same stem: still not one clade (the checker uses the sister clade of the first taxon)
  same_stem <- cc_clades(toy); same_stem$Sister_clade$M2 <- toy$B
  grafted_tips <- suppressWarnings(graft_missing_taxa(toy$phy, list(Taxon_name = c("M1", "M2"),
                                                                    Target_clade = list(M1 = toy$B, M2 = toy$B), Excluded_subclades = list(M1 = NULL, M2 = NULL)),
                                                      verbose = FALSE, seed = 2, include_target_clade_stem = TRUE))
  res <- cc_check(grafted_tips, toy, list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B)))
  expect_equal(nrow(res), 1)
  expect_true(is.logical(res$monophyletic))
})

test_that("A taxon missing from the grafted tree is reported, and the summary counts it", {
  toy <- toy_setup()
  grafted <- cc_graft(toy)
  extra <- cc_clades(toy)
  extra$Taxon_name <- c(extra$Taxon_name, "M5"); extra$Sister_clade$M5 <- toy$B
  res <- cc_check(grafted, toy, extra)
  expect_equal(res$Taxon_names[1], "M1, M2, M5")
  expect_equal(res$n_in_tree[1], 2)
  expect_false(res$in_tree[1])
  expect_false(res$ok[1])
  expect_true(is.na(res$sister_position[1]))
  expect_output(check_grafting_constraints(grafted, toy$phy, extra, verbose = TRUE), "Missing from the grafted tree:\\s+1")
  expect_error(check_grafting_constraints(grafted, toy$phy, extra, verbose = FALSE, stop_on_failure = TRUE), "1 missing clades failed")
})

test_that("The backbone check is the one of the taxa: a modified backbone is detected", {
  toy <- toy_setup()
  grafted <- cc_graft(toy)
  changed <- toy$phy
  changed$edge.length[1] <- changed$edge.length[1] * 0.5
  res <- suppressWarnings(check_grafting_constraints(grafted, changed, cc_clades(toy), verbose = FALSE))
  expect_false(attr(res, "backbone_preserved"))
  expect_false(attr(res, "all_ok"))
  res <- cc_check(grafted, toy, check_backbone = FALSE)
  expect_true(is.na(attr(res, "backbone_preserved")))
  expect_true(attr(res, "all_ok"))
})

test_that("The list must have a 'Target_clade' or a 'Sister_clade' item, not both", {
  toy <- toy_setup()
  grafted <- cc_graft(toy)
  both <- cc_clades(toy); both$Target_clade <- both$Sister_clade
  expect_error(cc_check(grafted, toy, both), "both 'Target_clade' and 'Sister_clade'")
  none <- cc_clades(toy); none$Sister_clade <- NULL
  expect_error(cc_check(grafted, toy, none), "must have an item 'Target_clade'")
  bad_tips <- cc_clades(toy); bad_tips$Sister_clade$M4 <- "not_a_tip"
  expect_error(cc_check(grafted, toy, bad_tips), "missing from the backbone")
})

test_that("States: the whole clade must be in one single allowed state", {
  toy <- toy_setup()
  sm <- make_toy_simmap(toy$phy, seed = 2)
  states <- colnames(sm$mapped.edge)
  grafted <- cc_graft(toy, simmap = sm)
  res <- cc_check(grafted, toy)
  expect_true(all(res$states_ok))
  expect_true(all(res$ok))
  # Allowed states of the grafting are respected; the global allowed_states of the check and the taxon-level ones are used
  used <- unique(unlist(lapply(grafted$simmap$maps[match(match(c("M1", "M2", "M3", "M4"), grafted$tip.label), grafted$edge[, 2])], names)))
  other <- setdiff(states, used)
  skip_if(length(other) == 0, "all the states are used by the toy clades")
  res <- cc_check(grafted, toy, allowed_states = other)
  expect_false(any(res$states_ok))
  expect_false(attr(res, "all_ok"))
  mcl_states <- cc_clades(toy); mcl_states$allowed_states <- stats::setNames(list(other, other, NULL, NULL), mcl_states$Taxon_name)
  res <- cc_check(grafted, toy, mcl_states, allowed_states = states)
  expect_equal(res$states_ok, c(FALSE, TRUE, TRUE))
  # Taxa of the same missing clade with different allowed states stop, as in the grafting function
  mcl_bad <- cc_clades(toy); mcl_bad$allowed_states <- stats::setNames(list(other, states, NULL, NULL), mcl_bad$Taxon_name)
  expect_error(cc_check(grafted, toy, mcl_bad), "different allowed states")
  # A taxon of the clade changed to another state: not one single state
  changed <- cc_set_edge_state(grafted, "M2", other[1])
  res <- cc_check(changed, toy)
  expect_false(res$states_ok[res$Taxon_names == "M1, M2"])
  expect_false(res$ok[res$Taxon_names == "M1, M2"])
  expect_output(check_grafting_constraints(changed, toy$phy, cc_clades(toy), verbose = TRUE), "not in a single state")
  # Without the updated simmap the states are not checked, and the arguments on the states stop
  expect_true(all(is.na(cc_check(drop_simmap(grafted), toy)$states_ok)))
  expect_error(cc_check(drop_simmap(grafted), toy, allowed_states = states), "can only be used if")
})

test_that("Ranges: the range inherited must be allowed and contained in the state of the lineage at the graft point", {
  toy <- toy_setup()
  sm <- make_toy_range_simmap(toy$phy, seed = 2)
  states <- colnames(sm$mapped.edge)
  grafted <- suppressWarnings(try(cc_graft(toy, simmap = sm, treat_states_as_ranges = TRUE), silent = TRUE))
  skip_if(inherits(grafted, "try-error"), "no valid grafting on the toy ranges")
  res <- cc_check(grafted, toy, treat_states_as_ranges = TRUE)
  expect_true(all(res$ranges_ok %in% TRUE))
  expect_true(all(res$ok))
  expect_false("states_ok_biogeo" %in% names(res))
  # A range that is not contained in the state of the lineage at the graft point is detected
  node <- ape::getMRCA(grafted, c("M1", "M2"))
  parent_edge <- match(grafted$edge[match(node, grafted$edge[, 2]), 1], grafted$edge[, 2])
  parent_state <- utils::tail(names(grafted$simmap$maps[[parent_edge]]), 1)
  not_contained <- states[!vapply(states, range_in_states, logical(1), state = parent_state)]
  skip_if(length(not_contained) == 0, "no range outside the graft state")
  # The whole clade (its stem branch and its two terminal branches) inherits a range that is not in the graft point state
  changed <- cc_set_edge_state(grafted, "M1", not_contained[1])
  changed <- cc_set_edge_state(changed, "M2", not_contained[1])
  stem_edge <- match(node, grafted$edge[, 2])
  changed$simmap$maps[[stem_edge]] <- stats::setNames(sum(changed$simmap$maps[[stem_edge]]), not_contained[1])
  res <- cc_check(changed, toy, treat_states_as_ranges = TRUE)
  expect_false(res$ranges_ok[res$Taxon_names == "M1, M2"])
  expect_false(res$ok[res$Taxon_names == "M1, M2"])
})

test_that("Several traits: the states are checked trait by trait for the whole clade, with the columns states_ok_<trait>", {
  toy <- toy_setup()
  sm1 <- make_toy_simmap(toy$phy, seed = 2)
  sm2 <- make_toy_simmap(toy$phy, seed = 3)
  s1 <- colnames(sm1$mapped.edge)
  maps <- list(color = sm1, shape = sm2)
  grafted <- suppressWarnings(try(cc_graft(toy, simmap = maps), silent = TRUE))
  skip_if(inherits(grafted, "try-error"), "no valid grafting with the toy traits")
  res <- cc_check(grafted, toy)
  expect_true(all(c("states_ok_color", "states_ok_shape") %in% names(res)))
  expect_true(all(res$states_ok) && all(res$ok))
  expect_equal(attr(res, "trait_names"), c("color", "shape"))
  # One trait changed for one taxon of a clade: only this trait is violated
  other <- setdiff(s1, names(grafted$simmap$color$maps[[match(match("M1", grafted$tip.label), grafted$edge[, 2])]]))
  skip_if(length(other) == 0)
  changed <- cc_set_edge_state(grafted, "M1", other[1], trait = "color")
  res <- cc_check(changed, toy)
  expect_false(res$states_ok_color[res$Taxon_names == "M1, M2"])
  expect_true(res$states_ok_shape[res$Taxon_names == "M1, M2"])
  expect_false(res$states_ok[res$Taxon_names == "M1, M2"])
  expect_output(check_grafting_constraints(changed, toy$phy, cc_clades(toy), verbose = TRUE), "color: 1, shape: 0")
})

test_that("Several simulations and several maps: the missing clades are checked in each phylogeny, with the aggregated summary", {
  toy <- toy_setup()
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  class(maps) <- c("multiSimmap", "multiPhylo"); names(maps) <- c("a", "b")
  grafted <- cc_graft(toy, simmap = maps, nb_simulations = 2)
  res <- cc_check(grafted, toy)
  expect_named(res, c("a", "b"))
  expect_equal(nrow(res$a), 6)
  expect_equal(res$a$simulation, rep(1:2, each = 3))
  expect_true(attr(res, "all_ok"))
  expect_output(check_grafting_constraints(grafted, toy$phy, cc_clades(toy), verbose = TRUE), "missing clades each")
  expect_output(check_grafting_constraints(grafted, toy$phy, cc_clades(toy), verbose = TRUE), "Not sister to their sister clade:\\s+0")
  # A phylogeny with a violation is listed
  wrong <- cc_clades(toy); wrong$Sister_clade$M3 <- toy$D
  expect_output(check_grafting_constraints(grafted, toy$phy, wrong, verbose = TRUE), "missing clades failed")
  expect_error(check_grafting_constraints(grafted, toy$phy, wrong, verbose = FALSE, stop_on_failure = TRUE), "missing clades failed")
  # The output of graft_missing_clades() with several simulations and no map
  plain <- cc_graft(toy, nb_simulations = 3)
  res <- cc_check(plain, toy)
  expect_equal(nrow(res), 9)
  expect_true(attr(res, "all_ok"))
})

test_that("The check of the taxa is unchanged by the detection of the missing clades", {
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$B), Excluded_subclades = list(X1 = list(toy$B), X2 = NULL))
  grafted <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, verbose = FALSE, seed = 1))
  res <- suppressWarnings(check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE))
  expect_true("Taxon_name" %in% names(res) && "target_ok" %in% names(res))
  expect_false("sister_ok" %in% names(res))
  expect_true(attr(res, "all_ok"))
  expect_output(check_grafting_constraints(grafted, toy$phy, mtl, verbose = TRUE), "Target clade violated")
})

test_that("The missing clades are labelled with the labels of their sister clades if the list has them, and numbered otherwise", {
  toy <- toy_setup()
  mcl <- cc_clades(toy)
  grafted <- cc_graft(toy, mcl)
  expect_identical(cc_check(grafted, toy, mcl)$missing_clade, 1:3)
  # Labels (one per taxon, shared by the taxa of a clade): used as the names of the missing clades
  mcl$Sister_clade_name <- c("cladeB", "cladeB", "cladeC", "tipS")
  res <- cc_check(grafted, toy, mcl)
  expect_identical(res$missing_clade, c("cladeB", "cladeC", "tipS"))
  expect_equal(res$Taxon_names, c("M1, M2", "M3", "M4"))
  # A clade without label keeps its ID as text
  mcl$Sister_clade_name <- c(NA, "cladeB", NA, "tipS")
  expect_identical(cc_check(grafted, toy, mcl)$missing_clade, c("cladeB", "2", "tipS"))
  # Several simulations: the same labels in each phylogeny
  sims <- cc_graft(toy, mcl, nb_simulations = 2)
  expect_identical(cc_check(sims, toy, mcl)$missing_clade, rep(c("cladeB", "2", "tipS"), 2))
})

