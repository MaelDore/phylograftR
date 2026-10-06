
##### Test suite for check_grafting_constraints() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Test correct validation with a correctly grafted phylogeny
test_that("check_grafting_constraints() validates a correct grafting",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxa_name = c("X1", "X2", "X3"),
                            Target_clade = list(toy$A, toy$B, toy$s),
                            Excluded_subclades = list(list(toy$B), NULL, NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, missing_taxa_list, method = "uniform", seed = 1)

  checks <- check_grafting_constraints(grafted, toy$phy, missing_taxa_list, verbose = FALSE)
  expect_s3_class(checks, "data.frame")
  expect_equal(nrow(checks), 3)
  expect_true(all(checks$target_ok))
  expect_true(all(checks$excluded_ok))
  expect_true(attr(checks, "topology_preserved"))
  expect_true(attr(checks, "backbone_preserved"))
  expect_true(attr(checks, "all_ok"))
})

## Test flagging for a taxon grafted inside its excluded subclade
test_that("check_grafting_constraints() flags a taxon grafted inside an excluded subclade",
{
  toy <- toy_setup()
  # X1 is grafted on a tip of B, but B is excluded for X1 (target clade A)
  inside_B <- graft_tip_on_branch(toy$phy, "X1", match(toy$tipsB[1], toy$phy$tip.label), 0.5)
  checks <- check_grafting_constraints(inside_B, toy$phy,
                                       list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(list(toy$B))),
                                       verbose = FALSE)
  expect_false(checks$excluded_ok)
  expect_false(attr(checks, "all_ok"))
})

## Test flagging for a taxon grafted outside its target clade
test_that("check_grafting_constraints() flags a taxon grafted outside its target clade",
{
  toy <- toy_setup()
  outside <- graft_tip_on_branch(toy$phy, "X1", match(toy$o, toy$phy$tip.label), 0.5)
  checks <- check_grafting_constraints(outside, toy$phy,
                                       list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(NULL)),
                                       verbose = FALSE)
  expect_false(checks$target_ok)
  expect_false(attr(checks, "all_ok"))
})

## Test flagging of a shift in backbone branch length, but not topology
test_that("check_grafting_constraints() reports a shift of the backbone branch lengths",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, missing_taxa_list, method = "uniform", seed = 1)
  shifted <- grafted
  shifted$edge.length[3] <- shifted$edge.length[3] * 1.01

  checks <- check_grafting_constraints(shifted, toy$phy, missing_taxa_list, verbose = FALSE)
  expect_true(attr(checks, "topology_preserved"))
  expect_false(attr(checks, "backbone_preserved"))
  expect_false(attr(checks, "all_ok"))
})

## Test flagging of a change in backbone topology
test_that("check_grafting_constraints() reports a change in the backbone topology",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, missing_taxa_list, method = "uniform", seed = 1)

  # Swap the labels of a backbone tip of clade A and of a tip outside it: same branch lengths, but a different backbone topology
  swapped <- grafted
  i <- match(c(toy$tipsA[1], toy$o), swapped$tip.label)
  swapped$tip.label[i] <- swapped$tip.label[rev(i)]

  checks <- check_grafting_constraints(swapped, toy$phy, missing_taxa_list, verbose = FALSE)
  expect_false(attr(checks, "topology_preserved"))
  expect_false(attr(checks, "backbone_preserved"))
  expect_false(attr(checks, "all_ok"))
  expect_error(check_grafting_constraints(swapped, toy$phy, missing_taxa_list, stop_on_failure = TRUE, verbose = FALSE))

  # Control: the unchanged grafted phylogeny is validated
  checks_ok <- check_grafting_constraints(grafted, toy$phy, missing_taxa_list, verbose = FALSE)
  expect_true(attr(checks_ok, "topology_preserved"))
  expect_true(attr(checks_ok, "all_ok"))
})

## Test the stem availability of each taxon
test_that("check_grafting_constraints() uses 'missing_species_list$include_target_clade_stem' of each taxon, and the global 'include_target_clade_stem' for NA or no item",
{
  toy <- toy_setup()
  # X1 and X2 are grafted on the stem branch of A
  on_stem <- graft_tip_on_branch(toy$phy, "X1", toy$node_A, 0.5)
  on_stem <- graft_tip_on_branch(on_stem, "X2", ape::getMRCA(on_stem, toy$A), 0.5)
  mtl <- list(Taxa_name = c("X1", "X2"), Target_clade = list(toy$A, toy$A), Excluded_subclades = list(NULL, NULL))
  check <- function (stem_item, global = TRUE)
  {
    mtl$include_target_clade_stem <- stem_item
    check_grafting_constraints(on_stem, toy$phy, mtl, include_target_clade_stem = global, verbose = FALSE)
  }

  # Stem accepted for X1 only
  checks <- check(c(TRUE, FALSE))
  expect_equal(checks$target_position, c("stem", "stem"))
  expect_equal(checks$include_target_clade_stem, c(TRUE, FALSE))
  expect_equal(checks$target_ok, c(TRUE, FALSE))
  expect_false(attr(checks, "all_ok"))

  # NA uses the global value
  expect_equal(check(c(NA, FALSE), global = TRUE)$target_ok, c(TRUE, FALSE))
  expect_equal(check(c(NA, TRUE), global = FALSE)$target_ok, c(FALSE, TRUE))
  # No item: the global value applies to all taxa
  expect_true(all(check(NULL, global = TRUE)$target_ok))
  expect_false(any(check(NULL, global = FALSE)$target_ok))

  # Invalid length
  expect_error(check(TRUE), "same length as 'Taxa_name'")
})

## Test that the function detect issues with taxa grafted on the stem of its excluded subclades
test_that("check_grafting_constraints() flags a taxon sister to an excluded subclade only if its stem is excluded",
{
  toy <- toy_setup()
  # X1 is grafted on the stem of B, which is excluded; X2 is grafted on the terminal branch of the excluded tip 's'
  on_stem_B <- graft_tip_on_branch(toy$phy, "X1", toy$node_B, 0.5)
  on_tip_s <- graft_tip_on_branch(toy$phy, "X2", match(toy$s, toy$phy$tip.label), 0.5)
  mtl <- list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(list(toy$B)))
  mtl_tip <- list(Taxa_name = "X2", Target_clade = list(toy$A), Excluded_subclades = list(list(toy$s)))

  check <- function (phy_grafted, mtl, ...) check_grafting_constraints(phy_grafted, toy$phy, mtl, verbose = FALSE, ...)

  # Stem excluded (default): violation
  checks <- check(on_stem_B, mtl)
  expect_false(checks$excluded_ok)
  expect_equal(checks$excluded_violated, "1")
  expect_false(check(on_tip_s, mtl_tip)$excluded_ok)
  # Stem not excluded (global argument): the taxon is sister to the subclade, not within it
  expect_true(check(on_stem_B, mtl, exclude_subclades_stem = FALSE)$excluded_ok)
  expect_true(check(on_tip_s, mtl_tip, exclude_subclades_stem = FALSE)$excluded_ok)
  # Setting of the subclade in the list overrides the global argument
  expect_true(check(on_stem_B, c(mtl, list(exclude_subclades_stem = list(FALSE))))$excluded_ok)
  expect_false(check(on_stem_B, c(mtl, list(exclude_subclades_stem = list(NA))), exclude_subclades_stem = TRUE)$excluded_ok)
  expect_false(check(on_tip_s, c(mtl_tip, list(exclude_subclades_stem = list(TRUE))), exclude_subclades_stem = FALSE)$excluded_ok)

  # A taxon within the excluded subclade is always a violation
  inside_B <- graft_tip_on_branch(toy$phy, "X1", match(toy$tipsB[1], toy$phy$tip.label), 0.5)
  expect_false(check(inside_B, mtl, exclude_subclades_stem = FALSE)$excluded_ok)

  # Several excluded subclades with their own setting: only the violated one is reported
  mtl_two <- list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(list(toy$C, toy$B)),
                  exclude_subclades_stem = list(c(TRUE, FALSE)))
  expect_true(check(on_stem_B, mtl_two)$excluded_ok) # B's stem is not excluded
  mtl_two$exclude_subclades_stem <- list(c(FALSE, TRUE))
  checks <- check(on_stem_B, mtl_two)
  expect_false(checks$excluded_ok)
  expect_equal(checks$excluded_violated, "2")

  # A taxon sister to a clade, but with another backbone tip in between, is not on its stem
  elsewhere <- graft_tip_on_branch(toy$phy, "X1", match(toy$tipsC[1], toy$phy$tip.label), 0.5)
  expect_true(check(elsewhere, mtl)$excluded_ok)

  # Invalid values
  expect_error(check(on_stem_B, c(mtl, list(exclude_subclades_stem = list(c(TRUE, FALSE))))), "excluded subclades")
  expect_error(check(on_stem_B, mtl, exclude_subclades_stem = NA), "single logical value")
})

## Test invalid input stops
test_that("check_grafting_constraints() stops on invalid inputs, and can stop on failure",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, missing_taxa_list, method = "uniform", seed = 1)

  # No excluded subclades in the list: warning, and not checked
  expect_warning(check_grafting_constraints(grafted, toy$phy, list(Taxa_name = "X1", Target_clade = list(toy$A)), verbose = FALSE),
                 "excluded subclades are not checked")
  # Inconsistent lengths
  expect_error(check_grafting_constraints(grafted, toy$phy, list(Taxa_name = c("X1", "X2"), Target_clade = list(toy$A),
                                                                 Excluded_subclades = list(NULL)), verbose = FALSE),
               "same length")
  # Tips missing from the backbone
  expect_error(check_grafting_constraints(grafted, toy$phy, list(Taxa_name = "X1", Target_clade = list(c("zz", "yy")),
                                                                 Excluded_subclades = list(NULL)), verbose = FALSE),
               "missing from the backbone")

  # stop_on_failure
  outside <- graft_tip_on_branch(toy$phy, "X1", match(toy$o, toy$phy$tip.label), 0.5)
  expect_error(check_grafting_constraints(outside, toy$phy, missing_taxa_list, stop_on_failure = TRUE, verbose = FALSE))
})
