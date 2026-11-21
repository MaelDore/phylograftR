##### Test suite for convert_values_to_states() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

states3 <- c("low", "medium", "high")

test_that("convert_values_to_states() assigns the states from the cut points, open-ended at both ends, and a cut point belongs to the state above", {
  breaks <- c(-1, 1)
  values <- c(-5, -1.0001, -1, 0, 0.9999, 1, 7)
  out <- convert_values_to_states(values, breaks = breaks, states = states3)
  expect_equal(unlist(out), c("low", "low", "medium", "medium", "medium", "high", "high"))
  # One cut point, two states
  expect_equal(unlist(convert_values_to_states(c(0, 1, 2), breaks = 1, states = c("a", "b"))), c("a", "b", "b"))
})

test_that("convert_values_to_states() always returns a list with one element per value, with the names and NA of the values", {
  values <- c(Taxon_1 = -2, Taxon_2 = NA, Taxon_3 = 3)
  out <- convert_values_to_states(values, breaks = c(-1, 1), states = states3)
  expect_type(out, "list")
  expect_equal(out, list(Taxon_1 = "low", Taxon_2 = NA_character_, Taxon_3 = "high"))
  expect_true(all(vapply(out, is.character, logical(1))))
  # Same structure with a tolerance of 0 (explicit or per value)
  expect_equal(convert_values_to_states(values, breaks = c(-1, 1), states = states3, tolerance = 0), out)
  expect_equal(convert_values_to_states(values, breaks = c(-1, 1), states = states3, tolerance = c(0, 0, 0)), out)
  # No names in the input: no names in the output
  expect_null(names(convert_values_to_states(c(0, 5), breaks = 1, states = c("a", "b"))))
  expect_equal(convert_values_to_states(c(0, 5), breaks = 1, states = c("a", "b")), list("a", "b"))
  # No value
  expect_equal(convert_values_to_states(numeric(0), breaks = 1, states = c("a", "b")), list())
  # The output can be unlisted into a named character vector when each value has a single state
  expect_equal(unlist(out), c(Taxon_1 = "low", Taxon_2 = NA, Taxon_3 = "high"))
})

test_that("convert_values_to_states() gives the same states as the simmap produced by convert_contMap_to_simmap(), single or list of contMaps", {
  toy <- toy_setup()
  set.seed(5)
  tips <- stats::setNames(phytools::fastBM(toy$phy), toy$phy$tip.label)
  cm <- suppressWarnings(phytools::contMap(toy$phy, tips, plot = FALSE))
  sm <- convert_contMap_to_simmap(cm, probs = c(1/3, 2/3), states = states3, verbose = FALSE)
  expect_equal(attr(sm, "states"), states3)
  # The states of the tips of the simmap are the states of the tip values
  tip_states <- phytools::getStates(sm, type = "tips")
  expect_equal(unlist(convert_values_to_states(tips, simmap = sm))[names(tip_states)], tip_states)
  # Same result with the cut points and states given directly
  expect_equal(convert_values_to_states(tips, simmap = sm), convert_values_to_states(tips, breaks = attr(sm, "breaks"), states = states3))
  # List of contMaps: multiSimmap with the same attributes
  sms <- convert_contMap_to_simmap(list(cm, cm), breaks = c(-1, 1), states = states3, verbose = FALSE)
  expect_equal(attr(sms, "breaks"), c(-1, 1))
  expect_equal(attr(sms, "states"), states3)
  expect_equal(unlist(convert_values_to_states(c(-3, 0, 3), simmap = sms)), c("low", "medium", "high"))
})

test_that("convert_values_to_states() with a tolerance returns all the states covered by [value - tolerance, value + tolerance]", {
  breaks <- c(-1, 1)
  values <- c(Far_low = -5, Near_low_cut = -0.9, Wide = 0, Near_high_cut = 0.8, Far_high = 5)
  out <- convert_values_to_states(values, breaks = breaks, states = states3, tolerance = 0.3)
  expect_type(out, "list")
  expect_named(out, names(values))
  expect_equal(out, list(Far_low = "low", Near_low_cut = c("low", "medium"), Wide = "medium",
                         Near_high_cut = c("medium", "high"), Far_high = "high"))
  # A large tolerance covers all the states, from the lowest to the highest
  expect_equal(convert_values_to_states(0, breaks = breaks, states = states3, tolerance = 5), list(states3))
  # The values far from the cut points keep the state they have without tolerance
  single <- convert_values_to_states(values, breaks = breaks, states = states3)
  expect_equal(out[c("Far_low", "Wide", "Far_high")], single[c("Far_low", "Wide", "Far_high")])
})

test_that("convert_values_to_states() with a tolerance follows the rule of the cut points at the edges of the interval", {
  # A cut point belongs to the state above it
  expect_equal(convert_values_to_states(0.5, breaks = 1, states = c("a", "b"), tolerance = 0.5), list(c("a", "b"))) # value + tolerance = cut point
  expect_equal(convert_values_to_states(1.5, breaks = 1, states = c("a", "b"), tolerance = 0.5), list("b")) # value - tolerance = cut point
  expect_equal(convert_values_to_states(0.5, breaks = 1, states = c("a", "b"), tolerance = 0.4999), list("a"))
})

test_that("convert_values_to_states() accepts one tolerance per value, and keeps names and NA with a tolerance", {
  values <- c(Taxon_1 = 0.9, Taxon_2 = 0.9, Taxon_3 = NA)
  out <- convert_values_to_states(values, breaks = 1, states = c("a", "b"), tolerance = c(0, 0.2, 0.2))
  expect_equal(out, list(Taxon_1 = "a", Taxon_2 = c("a", "b"), Taxon_3 = NA_character_))
  expect_null(names(convert_values_to_states(c(0.9, 5), breaks = 1, states = c("a", "b"), tolerance = 0.2)))
  expect_equal(convert_values_to_states(numeric(0), breaks = 1, states = c("a", "b"), tolerance = 0.2), list())
})

test_that("convert_values_to_states() works with a simmap, and gives a valid 'allowed_states' for the grafting plan, with or without tolerance", {
  data(toy_backbone, package = "phylograftR")
  data(toy_contMap, package = "phylograftR")
  data(toy_tip_values, package = "phylograftR")
  data(toy_missing_taxa_list, package = "phylograftR")
  sm <- convert_contMap_to_simmap(toy_contMap, probs = c(1/3, 2/3), states = states3, tip_values = toy_tip_values, verbose = FALSE)
  values <- c(New_sp1 = 0.5, New_sp2 = 2.5, New_sp3 = 1.4, New_sp4 = -3, New_sp5 = -1, New_sp6 = 3.5, New_sp7 = 1.6)
  allowed_single <- convert_values_to_states(values, simmap = sm)
  allowed <- convert_values_to_states(values, simmap = sm, tolerance = 0.3)
  expect_equal(allowed$New_sp3, c("medium", "high"))
  expect_equal(allowed$New_sp7, c("medium", "high"))
  expect_equal(lengths(allowed[c("New_sp1", "New_sp2", "New_sp4", "New_sp5", "New_sp6")]), c(New_sp1 = 1, New_sp2 = 1, New_sp4 = 1, New_sp5 = 1, New_sp6 = 1))
  expect_equal(lengths(allowed_single), c(New_sp1 = 1, New_sp2 = 1, New_sp3 = 1, New_sp4 = 1, New_sp5 = 1, New_sp6 = 1, New_sp7 = 1))
  # The output is accepted as the allowed states of the missing taxa: more states means a larger share of valid branches
  taxa_single <- taxa_tolerance <- toy_missing_taxa_list
  taxa_single$allowed_states <- allowed_single
  taxa_tolerance$allowed_states <- allowed
  plan_single <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_single, simmap = sm, verbose = FALSE)
  plan_tolerance <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_tolerance, simmap = sm, verbose = FALSE)
  expect_true(all(plan_tolerance$valid_fraction >= plan_single$valid_fraction - 1e-12))
  expect_true(all(plan_tolerance$valid_fraction[c(3, 7)] > plan_single$valid_fraction[c(3, 7)]))
})

test_that("convert_values_to_states() checks its arguments", {
  expect_error(convert_values_to_states("a", breaks = 1, states = c("a", "b")), "numeric vector")
  expect_error(convert_values_to_states(c(1, Inf), breaks = 1, states = c("a", "b")), "infinite")
  expect_error(convert_values_to_states(1), "either 'simmap', or both")
  expect_error(convert_values_to_states(1, breaks = 1), "either 'simmap', or both")
  expect_error(convert_values_to_states(1, states = c("a", "b")), "either 'simmap', or both")
  expect_error(convert_values_to_states(1, breaks = c(2, 1), states = states3), "increasing order")
  expect_error(convert_values_to_states(1, breaks = c(1, 1), states = states3), "increasing order")
  expect_error(convert_values_to_states(1, breaks = c(0, NA), states = states3), "finite")
  expect_error(convert_values_to_states(1, breaks = c(0, 1), states = c("a", "b")), "one more element than 'breaks'")
  expect_error(convert_values_to_states(1, breaks = c(0, 1), states = c("a", "b", "b")), "unique")
  expect_error(convert_values_to_states(1, breaks = c(0, 1), states = 1:3), "character")
  toy <- toy_setup()
  simmap_without <- make_toy_simmap(toy$phy, seed = 1)
  expect_error(convert_values_to_states(1, simmap = simmap_without), "no attributes")
  expect_error(convert_values_to_states(1, simmap = "x"), "'simmap' must be")
  expect_error(convert_values_to_states(1, simmap = simmap_without, breaks = 0, states = c("a", "b")), "not both")
  # Tolerance
  expect_error(convert_values_to_states(1, breaks = 1, states = c("a", "b"), tolerance = -0.1), "zero or positive")
  expect_error(convert_values_to_states(1, breaks = 1, states = c("a", "b"), tolerance = NA), "zero or positive")
  expect_error(convert_values_to_states(1, breaks = 1, states = c("a", "b"), tolerance = Inf), "zero or positive")
  expect_error(convert_values_to_states(1, breaks = 1, states = c("a", "b"), tolerance = "a"), "zero or positive")
  expect_error(convert_values_to_states(c(1, 2, 3), breaks = 1, states = c("a", "b"), tolerance = c(0.1, 0.2)), "one value per element")
})
