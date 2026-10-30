##### Test suite for convert_values_to_states() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

states3 <- c("low", "medium", "high")

test_that("convert_values_to_states() assigns the states from the cut points, open-ended at both ends, and a cut point belongs to the state above", {
  breaks <- c(-1, 1)
  values <- c(-5, -1.0001, -1, 0, 0.9999, 1, 7)
  expect_equal(convert_values_to_states(values, breaks = breaks, states = states3), c("low", "low", "medium", "medium", "medium", "high", "high"))
  # One cut point, two states
  expect_equal(convert_values_to_states(c(0, 1, 2), breaks = 1, states = c("a", "b")), c("a", "b", "b"))
})

test_that("convert_values_to_states() keeps the names and the NA, and returns a character vector", {
  values <- c(Taxon_1 = -2, Taxon_2 = NA, Taxon_3 = 3)
  out <- convert_values_to_states(values, breaks = c(-1, 1), states = states3)
  expect_type(out, "character")
  expect_equal(out, c(Taxon_1 = "low", Taxon_2 = NA, Taxon_3 = "high"))
  expect_null(names(convert_values_to_states(c(0, 5), breaks = 1, states = c("a", "b"))))
  expect_equal(convert_values_to_states(numeric(0), breaks = 1, states = c("a", "b")), character(0))
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
  expect_equal(convert_values_to_states(tips, simmap = sm)[names(tip_states)], tip_states)
  # Same result with the cut points and states given directly
  expect_equal(convert_values_to_states(tips, simmap = sm), convert_values_to_states(tips, breaks = attr(sm, "breaks"), states = states3))
  # List of contMaps: multiSimmap with the same attributes
  sms <- convert_contMap_to_simmap(list(cm, cm), breaks = c(-1, 1), states = states3, verbose = FALSE)
  expect_equal(attr(sms, "breaks"), c(-1, 1))
  expect_equal(attr(sms, "states"), states3)
  expect_equal(convert_values_to_states(c(-3, 0, 3), simmap = sms), c("low", "medium", "high"))
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
})
