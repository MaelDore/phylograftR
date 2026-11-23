##### Test suite for the toy datasets used to illustrate the preparation of stochastic maps #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Test that the tip data are the states at the tips of the toy simmaps
test_that("toy_tip_data holds the states of the tips of the toy simmaps, in the order of the tips of toy_backbone",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_tip_data, package = "phylograftR")
  data(toy_simmap, package = "phylograftR")
  data(toy_simmap_elevation, package = "phylograftR")
  data(toy_simmap_ranges, package = "phylograftR")
  expect_s3_class(toy_tip_data, "data.frame")
  expect_identical(row.names(toy_tip_data), toy_backbone$tip.label)
  expect_identical(names(toy_tip_data), c("color", "elevation", "range"))
  expect_equal(dim(toy_tip_data), c(30, 3))
  expect_false(anyNA(toy_tip_data))
  expect_true(all(vapply(toy_tip_data, is.character, logical(1))))
  expect_identical(toy_tip_data$color, unname(phytools::getStates(toy_simmap, type = "tips")[toy_backbone$tip.label]))
  expect_identical(toy_tip_data$elevation, unname(phytools::getStates(toy_simmap_elevation, type = "tips")[toy_backbone$tip.label]))
  expect_identical(toy_tip_data$range, unname(phytools::getStates(toy_simmap_ranges, type = "tips")[toy_backbone$tip.label]))
  expect_setequal(toy_tip_data$color, c("blue", "green", "red"))
  expect_setequal(toy_tip_data$elevation, c("high", "low"))
  expect_setequal(toy_tip_data$range, c("A", "B", "AB", "BC"))
})

## Test that the three multiSimmaps are valid sets of maps of the backbone, conditioned on the tip data
test_that("toy_multiSimmap, toy_multiSimmap_elevation and toy_multiSimmap_ranges are multiSimmaps of 10 maps of toy_backbone consistent with toy_tip_data",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_tip_data, package = "phylograftR")
  data(toy_multiSimmap, package = "phylograftR")
  data(toy_multiSimmap_elevation, package = "phylograftR")
  data(toy_multiSimmap_ranges, package = "phylograftR")
  all_maps <- list(color = toy_multiSimmap, elevation = toy_multiSimmap_elevation, range = toy_multiSimmap_ranges)
  tip_data <- list(color = toy_tip_data$color, elevation = toy_tip_data$elevation, range = toy_tip_data$range)
  for (trait in names(all_maps))
  {
    maps <- all_maps[[trait]]
    # The class multiSimmap is what identifies the maps as alternative maps of one trait
    expect_s3_class(maps, "multiSimmap")
    expect_s3_class(maps, "multiPhylo")
    expect_length(maps, 10)
    for (i in seq_along(maps))
    {
      expect_s3_class(maps[[i]], "simmap")
      expect_no_error(validate_simmap(maps[[i]], toy_backbone))
      expect_identical(maps[[i]]$edge, toy_backbone$edge)
      expect_equal(maps[[i]]$edge.length, toy_backbone$edge.length)
      # The maps are conditioned on the tip data
      tip_states <- phytools::getStates(maps[[i]], type = "tips")[toy_backbone$tip.label]
      expect_identical(unname(tip_states), tip_data[[trait]])
    }
    # The maps are not all the same: different histories on the internal branches
    n_changes <- vapply(maps, function (x) sum(lengths(x$maps) - 1), numeric(1))
    expect_gt(length(unique(n_changes)), 1)
    # One map extracted with [[ is a simmap, a subset with [ is still a multiSimmap
    expect_s3_class(maps[[1]], "simmap")
    expect_s3_class(maps[1:3], "multiSimmap")
  }
  # The states of the maps include the states of the tips
  expect_setequal(colnames(toy_multiSimmap[[1]]$mapped.edge), c("blue", "green", "red"))
  expect_setequal(colnames(toy_multiSimmap_elevation[[1]]$mapped.edge), c("high", "low"))
  all_ranges <- unique(unlist(lapply(toy_multiSimmap_ranges, function (x) colnames(x$mapped.edge))))
  expect_true(all(c("A", "B", "AB", "BC") %in% all_ranges))
  # The maps of colors and elevation have the same number of maps: they can be combined as two traits
  expect_equal(length(toy_multiSimmap), length(toy_multiSimmap_elevation))
})

## Test that the maps can be used as constraints, alone or as the maps of several traits
test_that("the toy multiSimmaps are accepted by the grafting plan, as one trait or as several traits",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_missing_taxa_list, package = "phylograftR")
  data(toy_multiSimmap, package = "phylograftR")
  data(toy_multiSimmap_elevation, package = "phylograftR")
  plan_one <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list, simmap = toy_multiSimmap, verbose = FALSE)
  expect_equal(attr(plan_one, "n_maps"), 10)
  plan_two <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                  simmap = list(color = toy_multiSimmap, elevation = toy_multiSimmap_elevation), verbose = FALSE)
  expect_equal(attr(plan_two, "n_maps"), 10)
  # A plain list of the same maps (class lost) is read as the maps of 10 traits, not as 10 maps of one trait
  plan_list <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list, simmap = lapply(toy_multiSimmap, identity), verbose = FALSE)
  expect_equal(attr(plan_list, "n_maps"), 1)
  expect_error(check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list, simmap = lapply(toy_multiSimmap, identity),
                                   allowed_states = c("blue", "green"), verbose = FALSE), "10 traits")
})
