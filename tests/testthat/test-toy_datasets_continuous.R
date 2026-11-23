##### Test suite for the toy datasets of a continuous trait #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Test that the tip values match the backbone, and the documented numbers
test_that("toy_tip_values is a named vector of the tip values of toy_backbone",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_tip_values, package = "phylograftR")
  expect_type(toy_tip_values, "double")
  expect_identical(names(toy_tip_values), toy_backbone$tip.label)
  expect_false(anyNA(toy_tip_values))
  # Documented range
  expect_equal(round(range(toy_tip_values), 1), c(-7.5, 4.2))
})

## Test that the contMap is the one of the tip values on the backbone
test_that("toy_contMap is the contMap of toy_tip_values on toy_backbone",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_tip_values, package = "phylograftR")
  data(toy_contMap, package = "phylograftR")
  expect_s3_class(toy_contMap, "contMap")
  expect_identical(toy_contMap$tree$tip.label, toy_backbone$tip.label)
  expect_identical(toy_contMap$tree$edge, toy_backbone$edge)
  expect_equal(toy_contMap$tree$edge.length, toy_backbone$edge.length)
  expect_equal(toy_contMap$lims, range(toy_tip_values))
  # The discretization is a valid simmap: three states of similar shares of branch length (documented numbers of toy_tip_values)
  my_simmap <- convert_contMap_to_simmap(toy_contMap, probs = c(1/3, 2/3), states = c("low", "medium", "high"), tip_values = toy_tip_values, verbose = FALSE)
  expect_s3_class(my_simmap, "simmap")
  expect_no_error(validate_simmap(my_simmap, toy_backbone))
  share <- colSums(my_simmap$mapped.edge) / sum(my_simmap$mapped.edge)
  expect_equal(round(unname(share[c("low", "medium", "high")]), 2), c(0.35, 0.43, 0.22))
})

## Test that the list of contMaps is a set of contMaps of the backbone, that can be discretized into a multiSimmap
test_that("toy_contMaps is a list of 10 contMaps of toy_backbone that can be discretized into a multiSimmap with three states",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_tip_values, package = "phylograftR")
  data(toy_contMaps, package = "phylograftR")
  expect_type(toy_contMaps, "list")
  expect_length(toy_contMaps, 10)
  for (i in c(1, 5, 10))
  {
    expect_s3_class(toy_contMaps[[i]], "contMap")
    expect_identical(toy_contMaps[[i]]$tree$tip.label, toy_backbone$tip.label)
    expect_identical(toy_contMaps[[i]]$tree$edge, toy_backbone$edge)
    expect_equal(toy_contMaps[[i]]$tree$edge.length, toy_backbone$edge.length)
    # The observed tip values are within the range of the values of the map
    tolerance <- 1e-6 * diff(toy_contMaps[[i]]$lims)
    expect_true(all(toy_tip_values >= toy_contMaps[[i]]$lims[1] - tolerance & toy_tip_values <= toy_contMaps[[i]]$lims[2] + tolerance))
  }
  # Discretization: the first step of the use of the list
  my_multiSimmap <- convert_contMap_to_simmap(toy_contMaps, probs = c(1/3, 2/3), states = c("low", "medium", "high"),
                                              tip_values = toy_tip_values, verbose = FALSE)
  expect_s3_class(my_multiSimmap, "multiSimmap")
  expect_length(my_multiSimmap, 10)
  expect_identical(attr(my_multiSimmap, "states"), c("low", "medium", "high"))
  breaks <- attr(my_multiSimmap, "breaks")
  expect_length(breaks, 2)
  expect_true(!is.unsorted(breaks, strictly = TRUE))
  # The cut points are the terciles of the observed tip values, and each tip is in the state of its observed value
  expect_equal(unname(breaks), unname(stats::quantile(toy_tip_values, probs = c(1/3, 2/3))))
  observed_states <- unlist(convert_values_to_states(toy_tip_values, simmap = my_multiSimmap))
  for (i in c(1, 5, 10))
  {
    expect_no_error(validate_simmap(my_multiSimmap[[i]], toy_backbone))
    expect_setequal(colnames(my_multiSimmap[[i]]$mapped.edge), c("low", "medium", "high"))
    tip_states <- phytools::getStates(my_multiSimmap[[i]], type = "tips")
    expect_identical(tip_states[names(observed_states)], observed_states)
  }
})
