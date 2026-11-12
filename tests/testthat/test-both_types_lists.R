##### Test suite for the functions handling list(taxa, clades) #####
# check_grafting_plan(), check_grafting_constraints(), extract_grafting_info(), extract_sampling_info()

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## ---- Helpers ---- ##

both_lists <- function (toy)
{
  list(taxa = list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$s)),
       clades = list(Taxon_name = c("M1", "M2", "M3"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C)))
}


## ---- Argument handling ---- ##

test_that("split_typed_list() recognises single lists, typed lists, one-element lists and invalid elements",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  s <- phylograftR:::split_typed_list(lists, "missing_taxa_list")
  expect_true(s$typed); expect_true(s$both); expect_equal(names(s$lists), c("taxa", "clades"))
  s1 <- phylograftR:::split_typed_list(lists$taxa, "missing_taxa_list")
  expect_false(s1$typed); expect_false(s1$both)
  s2 <- phylograftR:::split_typed_list(list(taxa = lists$taxa, clades = NULL), "missing_taxa_list")
  expect_true(s2$typed); expect_false(s2$both); expect_equal(names(s2$lists), "taxa")
  expect_error(phylograftR:::split_typed_list(list(taxa = NULL, clades = NULL), "missing_taxa_list"), "both are NULL")
  expect_error(phylograftR:::split_typed_list(list(taxa = lists$clades, clades = lists$clades), "missing_taxa_list"), "Target_clade")
})


## ---- check_grafting_plan() ---- ##

test_that("check_grafting_plan() checks both types and returns a named list",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  plans <- suppressWarnings(check_grafting_plan(toy$phy, lists, verbose = FALSE))
  expect_equal(names(plans), c("taxa", "clades"))
  expect_equal(plans$taxa, suppressWarnings(check_grafting_plan(toy$phy, lists$taxa, verbose = FALSE)))
  expect_equal(plans$clades, suppressWarnings(check_grafting_plan(toy$phy, lists$clades, verbose = FALSE)))
  expect_true(is.logical(attr(plans, "all_feasible")))
  # A single element is handled as that type
  one <- suppressWarnings(check_grafting_plan(toy$phy, list(taxa = lists$taxa, clades = NULL), verbose = FALSE))
  expect_equal(one, plans$taxa)
})


## ---- check_grafting_constraints() ---- ##

test_that("check_grafting_constraints() checks both types after grafting, in either order",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  for (ord in list(c("taxa", "clades"), c("clades", "taxa")))
  {
    g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, grafting_order = ord, method = "uniform", seed = 21)
    ch <- suppressWarnings(check_grafting_constraints(g, toy$phy, lists, verbose = FALSE))
    expect_equal(names(ch), c("taxa", "clades"))
    expect_true(isTRUE(attr(ch, "all_ok")))
    expect_true(isTRUE(attr(ch$taxa, "all_ok"))); expect_true(isTRUE(attr(ch$clades, "all_ok")))
    expect_true(all(ch$taxa$ok)); expect_true(all(ch$clades$ok))
  }
})

test_that("check_grafting_constraints() tolerates taxa grafted within a missing clade, and works for several simulations",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  # Taxon nested in a missing clade: X3 targets s; M3 (sister C) -- grafting clades first then the taxa
  g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, grafting_order = c("clades", "taxa"), method = "uniform", nb_simulations = 3, seed = 4)
  ch <- suppressWarnings(check_grafting_constraints(g, toy$phy, lists, verbose = FALSE))
  expect_true(isTRUE(attr(ch, "all_ok")))
  expect_equal(nrow(attr(ch$taxa, "phylogenies")), 3)
})

test_that("check_grafting_constraints() reports a failure and stops only at the end when stop_on_failure = TRUE",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", seed = 8)
  bad <- lists; bad$taxa$Target_clade$X2 <- toy$C  # wrong target clade for X2 (it was grafted in B)
  ch <- suppressWarnings(check_grafting_constraints(g, toy$phy, bad, verbose = FALSE))
  expect_false(isTRUE(attr(ch, "all_ok")))
  expect_false(isTRUE(attr(ch$taxa, "all_ok")))
  expect_true(isTRUE(attr(ch$clades, "all_ok")))
  expect_error(suppressWarnings(check_grafting_constraints(g, toy$phy, bad, verbose = FALSE, stop_on_failure = TRUE)), "taxa")
  expect_output(suppressWarnings(check_grafting_constraints(g, toy$phy, lists, verbose = TRUE)), "Missing clades")
})

test_that("a list with one non-NULL element gives the same output as that single type",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  g <- quiet_graft(graft_missing_taxa, toy$phy, lists$taxa, method = "uniform", seed = 2)
  a <- suppressWarnings(check_grafting_constraints(g, toy$phy, lists$taxa, verbose = FALSE))
  b <- suppressWarnings(check_grafting_constraints(g, toy$phy, list(taxa = lists$taxa, clades = NULL), verbose = FALSE))
  expect_equal(a, b)
})


## ---- extract_grafting_info() ---- ##

test_that("extract_grafting_info() returns the information of both types, and single-type calls work on a two-step tree",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  for (ord in list(c("taxa", "clades"), c("clades", "taxa")))
  {
    g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, grafting_order = ord, method = "uniform", seed = 6)
    info <- suppressWarnings(extract_grafting_info(g, lists, phy_backbone = toy$phy))
    expect_equal(names(info), c("taxa", "clades"))
    expect_equal(info$taxa$Taxon_name, lists$taxa$Taxon_name)
    expect_equal(nrow(info$clades), 2) # One row per missing clade: sister B (M1, M2) and sister C (M3)
    expect_false(anyNA(info$taxa$graft_age)); expect_false(anyNA(info$clades$graft_age))
    # Without the backbone: rebuilt by dropping all the flagged tips
    info2 <- suppressWarnings(extract_grafting_info(g, lists))
    expect_equal(names(info2), c("taxa", "clades"))
    # Single type on the two-step tree: the other type is dropped using the flags
    single <- suppressWarnings(extract_grafting_info(g, lists$taxa, phy_backbone = toy$phy))
    expect_equal(single$Taxon_name, lists$taxa$Taxon_name)
  }
})

test_that("extract_grafting_info() works for several simulations and with a one-element list",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", nb_simulations = 2, seed = 9)
  info <- suppressWarnings(extract_grafting_info(g, lists, phy_backbone = toy$phy))
  expect_equal(names(info), c("taxa", "clades"))
  expect_equal(length(unique(info$taxa$simulation)), 2)
  one <- suppressWarnings(extract_grafting_info(g, list(taxa = lists$taxa, clades = NULL), phy_backbone = toy$phy))
  expect_equal(one, suppressWarnings(extract_grafting_info(g, lists$taxa, phy_backbone = toy$phy)))
})


## ---- extract_sampling_info() ---- ##

test_that("extract_sampling_info() returns the sampling of both types as a named list",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  richness <- list(taxa = list(Target_clade = list(toy$A, toy$B, toy$s), total_richness = c(toy$nA + 15, toy$nB + 5, NA)),
                   clades = list(Sister_clade = list(toy$B, toy$C), total_richness = c(toy$nB + 12, NA)))
  g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, total_richness_per_clades = richness, method = "uniform", seed = 2)
  info <- suppressWarnings(extract_sampling_info(g, richness))
  expect_equal(names(info), c("taxa", "clades"))
  expect_equal(info$taxa, suppressWarnings(extract_sampling_info(g, richness$taxa)))
  expect_equal(info$clades, suppressWarnings(extract_sampling_info(g, richness$clades)))
  one <- suppressWarnings(extract_sampling_info(g, list(taxa = richness$taxa, clades = NULL)))
  expect_equal(one, info$taxa)
})
