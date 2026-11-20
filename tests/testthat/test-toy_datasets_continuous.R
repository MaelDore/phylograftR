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

## Test that the multiSimmap is a valid set of simmaps of the backbone, with the attributes to convert values
test_that("toy_multiSimmap_continuous is a multiSimmap of 100 valid simmaps of toy_backbone with three states",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_tip_values, package = "phylograftR")
  data(toy_multiSimmap_continuous, package = "phylograftR")
  expect_s3_class(toy_multiSimmap_continuous, "multiSimmap")
  expect_length(toy_multiSimmap_continuous, 100)
  expect_identical(attr(toy_multiSimmap_continuous, "states"), c("low", "medium", "high"))
  breaks <- attr(toy_multiSimmap_continuous, "breaks")
  expect_length(breaks, 2)
  expect_true(!is.unsorted(breaks, strictly = TRUE))
  for (i in c(1, 50, 100))
  {
    expect_no_error(validate_simmap(toy_multiSimmap_continuous[[i]], toy_backbone))
    expect_identical(toy_multiSimmap_continuous[[i]]$tip.label, toy_backbone$tip.label)
    expect_setequal(colnames(toy_multiSimmap_continuous[[i]]$mapped.edge), c("low", "medium", "high"))
  }
  # The cut points are the terciles of the observed tip values, and each tip is in the state of its observed value
  expect_equal(unname(breaks), unname(stats::quantile(toy_tip_values, probs = c(1/3, 2/3))))
  observed_states <- convert_values_to_states(toy_tip_values, simmap = toy_multiSimmap_continuous)
  for (i in c(1, 50, 100))
  {
    tip_states <- phytools::getStates(toy_multiSimmap_continuous[[i]], type = "tips")
    expect_identical(tip_states[names(observed_states)], observed_states)
  }
})
