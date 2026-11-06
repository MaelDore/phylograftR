
##### Test suite for graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Build toy 'missing_taxa_list' to test grafting of missing caldes
# Missing taxa: M1 and M2 share the sister clade B (grafted together as one clade), M3 is sister to C, M4 is sister to a single tip
toy_missing_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}

## Check if removing the grafted taxa gives back the backbone (topology and branch lengths)
expect_backbone_preserved <- function (grafted, phy, taxa)
{
  pruned <- ape::drop.tip(grafted, taxa)
  expect_true(ape::all.equal.phylo(pruned, phy, use.edge.length = TRUE, tolerance = 1e-6))
}

## Test invalid input stops
test_that("graft_missing_clades() stops on invalid inputs, before any grafting",
{
  toy <- toy_setup()
  phy <- toy$phy
  mcl <- toy_missing_clades(toy)

  not_ultra <- phy
  not_ultra$edge.length[1] <- not_ultra$edge.length[1] * 1.5
  expect_error(graft_missing_clades(not_ultra, mcl, verbose = FALSE), "not ultrametric")

  poly <- make_toy_polytomies(phy, n_collapsed = 1)
  expect_error(graft_missing_clades(poly, mcl, verbose = FALSE), "not binary")
  expect_error(graft_missing_clades(poly, mcl, verbose = FALSE), "ape::multi2di")

  expect_error(graft_missing_clades(phy, mcl, min_tips = 2, verbose = FALSE), "min_tips")
  expect_error(graft_missing_clades(make_toy_phy(n_tips = 4), list(Taxon_name = "M1", Sister_clade = list(M1 = c("t1", "t2"))), verbose = FALSE),
               "too small to fit a birth-death model")
  expect_error(graft_missing_clades(phy, list(Taxon_name = character(0), Sister_clade = list()), verbose = FALSE), "empty")
  expect_error(graft_missing_clades(phy, list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B)), verbose = FALSE), "same length as 'Taxon_name'")
  expect_error(graft_missing_clades(phy, list(Taxon_name = c("M1", "M1"), Sister_clade = list(M1 = toy$B, M1 = toy$B)), verbose = FALSE), "unique names")
  expect_error(graft_missing_clades(phy, list(Taxon_name = "M1", Sister_clade = list(M1 = c("zz", "yy"))), verbose = FALSE), "missing from the backbone")
  expect_error(graft_missing_clades(phy, mcl, temp_file_path = "temp.txt", verbose = FALSE), "must end with '.rds'")

  # A sister clade that is the whole tree has no stem
  whole_tree <- pair_for_node(phy, ape::Ntip(phy) + 1)
  expect_error(graft_missing_clades(phy, list(Taxon_name = "M1", Sister_clade = list(M1 = whole_tree)), verbose = FALSE), "whole tree")
})

## Test valid location of grafted taxa on the stem of their sister clade
test_that("graft_missing_clades() grafts missing clades on the stem of their sister clade",
{
  toy <- toy_setup()
  mcl <- toy_missing_clades(toy)

  for (seed in 1:2)
  {
    grafted <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = seed)
    expect_true(ape::is.ultrametric(grafted))
    expect_equal(ape::Ntip(grafted), ape::Ntip(toy$phy) + 4)
    expect_backbone_preserved(grafted, toy$phy, mcl$Taxon_name)

    # M1 and M2 form a monophyletic missing clade, sister to B
    node_M <- ape::getMRCA(grafted, c("M1", "M2"))
    expect_setequal(grafted$tip.label[ape::prop.part(grafted)[[node_M - ape::Ntip(grafted)]]], c("M1", "M2"))
    sister_of_M <- setdiff(grafted$edge[grafted$edge[, 1] == grafted$edge[match(node_M, grafted$edge[, 2]), 1], 2], node_M)
    expect_equal(sister_of_M, ape::getMRCA(grafted, toy$B))
  }
})

## Test handling of nested sister clades
test_that("graft_missing_clades() handles nested sister clades, and is reproducible",
{
  toy <- toy_setup()
  # A and B are nested sister clades, each with a missing taxon
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$A))
  a <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 4)
  b <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 4)
  expect_equal(a, b)
  expect_true(ape::is.ultrametric(a))
  expect_backbone_preserved(a, toy$phy, mcl$Taxon_name)
  # M2 is sister to A, so M1 (sister to B, in A) is not sister to M2
  expect_false(ape::getMRCA(a, c("M1", "M2")) == ape::getMRCA(a, c("M1", toy$B[1])))
})

## Test the temporary saves
test_that("graft_missing_clades() returns the grafted tips and branches, and saves temporary phylogenies",
{
  toy <- toy_setup()
  mcl <- toy_missing_clades(toy)
  temp_file <- file.path(tempdir(), "temp_graft_clades.rds")
  unlink(temp_file)

  grafted <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1, save_temp = TRUE, save_temp_every = 1, clean_temp_files = FALSE, temp_file_path = temp_file)
  expect_setequal(grafted$tip.label[grafted$grafted_tips], mcl$Taxon_name)
  expect_identical(names(grafted$grafted_edges), as.character(grafted$edge[, 2]))
  expect_null(grafted$sampling_fraction_df)
  expect_s3_class(readRDS(temp_file), "phylo")
  unlink(temp_file)
})

## Test the sampling fraction of the fits without total richness information
test_that("graft_missing_clades() fits the crown group, accounting for the taxa still to graft nested in the fitted clade only",
{
  toy <- toy_setup()
  # M2, M3: sister to B (small clade: the fit uses an enclosing clade). M1: sister to A (outside of A)
  mcl <- list(Taxon_name = c("M1", "M2", "M3"), Sister_clade = list(M1 = toy$A, M2 = toy$B, M3 = toy$B))

  fits <- list()
  original <- fit_bd_on_clade
  local_mocked_bindings(fit_bd_on_clade = function (...) { out <- original(...); fits[[length(fits) + 1]] <<- out; out })
  quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1)

  # Fit for B (first): enclosing clade, in which M2 and M3 are nested and pending. M1 is outside (sister to A, which is not strictly nested)
  expect_gt(fits[[1]]$n_tips, toy$nB)
  expect_equal(fits[[1]]$n_pending, 2)
  expect_equal(fits[[1]]$n_unsampled, 0)
  expect_equal(fits[[1]]$rho, fits[[1]]$n_tips / (fits[[1]]$n_tips + 2))

  # Fit for A (second): A itself, M2 and M3 are grafted already (present), M1 is sister to A, outside of its crown group: sampling fraction of 1
  expect_equal(fits[[2]]$n_tips, toy$nA + 2)
  expect_equal(fits[[2]]$n_pending, 0)
  expect_equal(fits[[2]]$n_unsampled, 0)
  expect_equal(fits[[2]]$rho, 1)
})

## Test that NULL and uninformative total richness give the same fits
test_that("graft_missing_clades() gives the same result with NULL and with only NA in total_richness_per_clades",
{
  toy <- toy_setup()
  mcl <- toy_missing_clades(toy)
  totals_NA <- list(Sister_clade = list(toy$B, toy$C, toy$s), total_richness = rep(NA_real_, 3))

  grafted_null <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1)
  grafted_NA <- quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = totals_NA, seed = 1)

  expect_equal(grafted_NA$edge.length, grafted_null$edge.length)
  expect_equal(grafted_NA$tip.label, grafted_null$tip.label)

  expect_null(grafted_null$sampling_fraction_df)
  expect_null(grafted_NA$sampling_fraction_df)

  # Nothing is unsampled: the sampling fraction extracted from the tree is 1
  info <- extract_sampling_info(grafted_NA, totals_NA)
  expect_false(any(info$total_richness_provided))
  expect_equal(info$nb_unsampled_species, rep(0, 3))
  expect_equal(info$sampling_fraction, rep(1, 3))
})

## Test the sampling fraction of the fits with total richness information
test_that("total_richness_per_clades is used in the sampling fraction of the fits, with missing taxa counted only if nested",
{
  toy <- toy_setup()
  # M1: sister to B (on the stem of B: in A, not in B). M2: sister to A (outside A).
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$A))
  total_B <- toy$nB + 7
  total_A <- toy$nA + 1 + 11 # M1 is nested in A, M2 is not
  totals <- list(Sister_clade = list(toy$B, toy$A), total_richness = c(total_B, total_A))

  fits <- list()
  original <- fit_bd_on_clade
  local_mocked_bindings(fit_bd_on_clade = function (...) { out <- original(...); fits[[length(fits) + 1]] <<- out; out })
  # min_tips = 3 so that B (3 tips) is fitted on itself and not on an enclosing clade
  grafted <- quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = totals, min_tips = 3, seed = 1)

  # Fit on B (first): nothing missing inside B (M1 is on its stem), 7 unsampled
  expect_equal(fits[[1]]$n_pending, 0)
  expect_equal(fits[[1]]$n_unsampled, 7)
  expect_equal(fits[[1]]$rho, fits[[1]]$n_tips / (fits[[1]]$n_tips + 7))

  # Fit on A (second): M1 is already grafted (present), M2 (sister to A) is outside A, so nothing pending. 11 unsampled
  expect_equal(fits[[2]]$n_tips, toy$nA + 1)
  expect_equal(fits[[2]]$n_pending, 0)
  expect_equal(fits[[2]]$n_unsampled, 11)

  # Realized counts from the tree: M1 (sister to B) is in A but not in B, M2 (sister to A) is outside A
  expect_null(grafted$sampling_fraction_df)
  table <- extract_sampling_info(grafted, totals)
  expect_equal(table$nb_grafted_species, c(0, 1))
  expect_equal(table$nb_unsampled_species, c(7, 11))
})

## Test check stops for 'total_richness_per_clades'
test_that("graft_missing_clades() stops before grafting if total_richness_per_clades is invalid",
{
  toy <- toy_setup()
  mcl <- toy_missing_clades(toy)
  expect_error(graft_missing_clades(toy$phy, mcl, total_richness_per_clades = list(Sister_clade = list(toy$B), total_richness = toy$nB + 1), verbose = FALSE),
               "not listed")
})


## Test that Sister_clade must be named after Taxon_name
test_that("graft_missing_clades() requires Sister_clade to be named after Taxon_name (same names, same order)",
{
  toy <- toy_setup()
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$A))
  expect_no_error(quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1))

  unnamed <- mcl; names(unnamed$Sister_clade) <- NULL
  expect_error(graft_missing_clades(toy$phy, unnamed, verbose = FALSE), "Sister_clade.*must be named after 'Taxon_name'")
  wrong <- mcl; names(wrong$Sister_clade) <- c("M1", "M3")
  expect_error(graft_missing_clades(toy$phy, wrong, verbose = FALSE), "Sister_clade.*do not match.*'M3'")
  swapped <- mcl; swapped$Sister_clade <- rev(swapped$Sister_clade)
  expect_error(graft_missing_clades(toy$phy, swapped, verbose = FALSE), "not in the same order")
})

## Test that node labels of the backbone are kept, with an empty string for the nodes created by grafting, and no labels if the backbone had none
test_that("graft_missing_clades() keeps the node labels of the backbone and gives an empty string to the new nodes",
{
  toy <- toy_setup()
  mcl <- toy_missing_clades(toy)
  phy_labelled <- toy$phy
  phy_labelled$node.label <- paste0("support", seq_len(phy_labelled$Nnode))

  # Backbone without labels
  grafted <- suppressWarnings(graft_missing_clades(toy$phy, mcl, verbose = FALSE, seed = 1))
  expect_null(grafted$node.label)

  # Backbone with labels: same tree as without labels, original labels kept, new nodes labelled with an empty string
  grafted_labelled <- suppressWarnings(graft_missing_clades(phy_labelled, mcl, verbose = FALSE, seed = 1))
  expect_equal(grafted_labelled$edge, grafted$edge)
  expect_equal(grafted_labelled$edge.length, grafted$edge.length)
  expect_length(grafted_labelled$node.label, grafted_labelled$Nnode)
  expect_setequal(grafted_labelled$node.label[grafted_labelled$node.label != ""], phy_labelled$node.label)
  expect_equal(sum(grafted_labelled$node.label != ""), phy_labelled$Nnode)
  expect_false(anyNA(grafted_labelled$node.label))
})


# ---- Species attributes: allowed states of a simmap ----

## Test that all the states allowed gives exactly the same phylogeny as without simmap, with global or taxon-level states
test_that("graft_missing_clades() with all the states allowed gives the same phylogeny as without simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mcl <- toy_missing_clades(toy)

  for (seed in 1:3)
  {
    without <- graft_missing_clades(toy$phy, mcl, verbose = FALSE, seed = seed)

    # Global states
    with_all <- graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B", "C"), verbose = FALSE, seed = seed)
    expect_identical(drop_simmap(with_all), without)

    # Taxon-level states, in a different order, with the global states for some taxa (the taxa of a missing clade share the same states)
    mcl_taxa <- mcl
    mcl_taxa$allowed_states <- list(M1 = c("C", "A", "B"), M2 = c("B", "C", "A"), M3 = NULL, M4 = NA)
    with_taxa <- graft_missing_clades(toy$phy, mcl_taxa, simmap = simmap, allowed_states = c("A", "B", "C"), verbose = FALSE, seed = seed)
    expect_identical(drop_simmap(with_taxa), without)
  }
})

## Test that the missing clades are only grafted on the pieces of the stem branch in their allowed states
test_that("graft_missing_clades() grafts the missing clades on the parts of the stem branches in their allowed states, and redraws the divergence times otherwise",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # The stem branch of B: state "A" in its rootward half, state "C" in its tipward half. The terminal branch of the tip 's': "B" then "A" (from the root to the tip)
  simmap <- set_branch_states(simmap, toy$B, c("A", "C"), c(1, 1))
  simmap <- set_branch_states(simmap, toy$s, c("B", "A"), c(1, 1))
  mcl <- toy_missing_clades(toy)
  mcl$allowed_states <- list(M1 = "C", M2 = "C", M3 = NULL, M4 = "A")

  crown_age_B <- toy$root_age - phytools::nodeheight(toy$phy, toy$node_B)
  stem_B <- toy$phy$edge.length[match(toy$node_B, toy$phy$edge[, 2])]
  tip_s <- toy$phy$edge.length[match(match(toy$s, toy$phy$tip.label), toy$phy$edge[, 2])]

  ages_B <- ages_s <- numeric(0)
  for (seed in 1:15)
  {
    grafted <- quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B", "C"), seed = seed)
    expect_true(ape::is.ultrametric(grafted))
    expect_backbone_preserved(grafted, toy$phy, mcl$Taxon_name)

    # M1 and M2 (state "C") are grafted on the tipward half of the stem of B, and M4 (state "A") on the tipward half of the terminal branch of 's'
    ages_B <- c(ages_B, graft_age(grafted, c("M1", "M2")))
    ages_s <- c(ages_s, graft_age(grafted, "M4"))
    expect_identical(state_of_lineage_at_age(simmap, toy$B, tail(ages_B, 1)), "C")
    expect_identical(state_of_lineage_at_age(simmap, toy$s, tail(ages_s, 1)), "A")
  }
  expect_true(all(ages_B >= crown_age_B - 1e-6 & ages_B <= crown_age_B + stem_B / 2 + 1e-6))
  expect_true(all(ages_s >= 0 & ages_s <= tip_s / 2 + 1e-6))

  # The ages are not all the same: they are spread over the allowed half
  expect_gt(length(unique(round(ages_B, 6))), 5)
})

## Test that the divergence times are redrawn (a truncated draw), by counting the draws of TreeSim::corsim()
test_that("graft_missing_clades() draws the divergence time again until it is in an allowed state",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  simmap <- set_branch_states(simmap, toy$B, c("A", "C"), c(1, 1))
  mcl <- list(Taxon_name = "M1", Sister_clade = list(M1 = toy$B), allowed_states = list(M1 = "C"))

  n_draws <- 0
  original_corsim <- TreeSim::corsim
  local_mocked_bindings(corsim = function (...) { n_draws <<- n_draws + 1; original_corsim(...) }, .package = "TreeSim")
  for (seed in 1:20) { quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, seed = seed) }
  # Without restriction, one draw per run. Here a part of the draws falls in the rootward half and is drawn again
  expect_gt(n_draws, 20)
})

## Test the stops and warning before grafting
test_that("graft_missing_clades() stops or warns before grafting on invalid or problematic 'simmap' and 'allowed_states'",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  simmap <- set_branch_states(simmap, toy$B, c("A", "C"), c(1, 1))
  mcl <- toy_missing_clades(toy)
  graft <- function (...) graft_missing_clades(toy$phy, mcl, verbose = FALSE, ...)

  # Allowed states must be in the simmap (a simmap without allowed states is valid: all the states are allowed)
  expect_no_error(graft(simmap = simmap))
  expect_error(graft(simmap = simmap, allowed_states = c("A", "red")), "global 'allowed_states' argument.*'red'.*'A', 'B', 'C'")
  # Allowed states without simmap
  expect_error(graft(allowed_states = "A"), "requires a 'simmap'")
  mcl_states <- mcl
  mcl_states$allowed_states <- list(M1 = "A", M2 = "A", M3 = "A", M4 = "A")
  expect_error(graft_missing_clades(toy$phy, mcl_states, verbose = FALSE), "requires a 'simmap'")
  # A plain list of maps
  expect_error(graft(simmap = list(simmap, simmap), allowed_states = "A"), "must be given for each trait") # A plain list is the maps of several traits
  # Other arguments
  for (value in list(-0.1, 1.5, NA, "0.1", c(0.1, 0.2)))
  {
    expect_error(graft(simmap = simmap, allowed_states = "A", min_valid_fraction = value), "'min_valid_fraction' must be a single number between 0 and 1")
  }
  for (value in list(0, -1, 1.5, NA, "10", c(10, 20)))
  {
    expect_error(graft(simmap = simmap, allowed_states = "A", max_age_redraws = value), "'max_age_redraws' must be a single positive integer")
  }

  # Taxa of the same missing clade (M1 and M2, same sister clade) with different allowed states
  mcl_inconsistent <- mcl
  mcl_inconsistent$allowed_states <- list(M1 = "A", M2 = c("A", "C"), M3 = NULL, M4 = NULL)
  expect_error(graft_missing_clades(toy$phy, mcl_inconsistent, simmap = simmap, allowed_states = "A", verbose = FALSE),
               "1 missing clade\\(s\\) have taxa with different allowed states.*sister clade tips: .*M1 \\(A\\), M2 \\(A, C\\)")
  # The same set of states in a different order is not a difference
  mcl_inconsistent$allowed_states <- list(M1 = c("A", "C"), M2 = c("C", "A"), M3 = NULL, M4 = NULL)
  expect_no_error(graft_missing_clades(toy$phy, mcl_inconsistent, simmap = simmap, allowed_states = c("A", "B", "C"), verbose = FALSE))

  # A stem branch without any piece in the allowed states: the entire stem of the tip 's' is in state "A"
  mcl_none <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$s, M2 = toy$B), allowed_states = list(M1 = "B", M2 = c("A", "C")))
  expect_error(graft_missing_clades(toy$phy, mcl_none, simmap = simmap, verbose = FALSE),
               "1 missing clade\\(s\\) have no piece of the stem branch of their sister clade in their allowed states.*sister clade tips: .*allowed states: B")

  # A small share of the stem branch in the allowed states: warning listing the clades. The stem of B is half in state "C"
  mcl_B <- list(Taxon_name = "M1", Sister_clade = list(M1 = toy$B))
  expect_warning(graft_missing_clades(toy$phy, mcl_B, simmap = simmap, allowed_states = "C", min_valid_fraction = 0.6, verbose = FALSE),
                 "1 missing clade\\(s\\) have less than 60% of the stem branch of their sister clade in their allowed states.*50.00% in allowed states")
  expect_no_warning(graft_missing_clades(toy$phy, mcl_B, simmap = simmap, allowed_states = "C", min_valid_fraction = 0.4, verbose = FALSE))
})

## Test the stop when no usable divergence time is found
test_that("graft_missing_clades() stops after 'max_age_redraws' redraws of the divergence time, with an informative message",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # Only a tiny piece of the stem branch of B (at its tipward end) is in state "C"
  simmap <- set_branch_states(simmap, toy$B, c("A", "C"), c(1, 1e-7))
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B))
  expect_error(suppressWarnings(graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = "C", max_age_redraws = 3, verbose = FALSE, seed = 1)),
               "No usable divergence time.*M1, M2.*after 3 redraws.*'max_age_redraws'.*allowed states: C")

  # The message gives the time window of the divergence time (the stem branch of the sister clade) and the share of it in the allowed states
  crown_age_B <- toy$root_age - phytools::nodeheight(toy$phy, toy$node_B)
  stem_age_B <- toy$root_age - phytools::nodeheight(toy$phy, toy$phy$edge[match(toy$node_B, toy$phy$edge[, 2]), 1])
  expect_error(suppressWarnings(graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = "C", max_age_redraws = 3, verbose = FALSE, seed = 1)),
               paste0("Time window of the divergence time \\(stem branch of the sister clade\\): from age ", signif(crown_age_B, 4), " to age ", signif(stem_age_B, 4),
                      "\\. Share of this window in the allowed states: 0.00%"))
})

## Test the grafting of missing clades after the grafting of nested ones with a simmap (tracking of the states through bind.tree), and the node labels
test_that("graft_missing_clades() with a simmap handles nested sister clades, several taxa per missing clade and node labels",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # B and A are nested sister clades, each with a missing clade of several taxa
  mcl <- list(Taxon_name = c("M1", "M2", "M3", "M4", "M5"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$B, M4 = toy$A, M5 = toy$A))

  for (seed in 1:3)
  {
    without <- graft_missing_clades(toy$phy, mcl, verbose = FALSE, seed = seed)
    with_all <- graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B", "C"), verbose = FALSE, seed = seed)
    expect_identical(drop_simmap(with_all), without)
  }

  phy_support <- toy$phy
  phy_support$node.label <- paste0("support", seq_len(phy_support$Nnode))
  grafted <- graft_missing_clades(phy_support, mcl, simmap = simmap, allowed_states = c("A", "B"), verbose = FALSE, seed = 1)
  expect_length(grafted$node.label, grafted$Nnode)
  expect_setequal(grafted$node.label[grafted$node.label != ""], phy_support$node.label)
  expect_false(anyNA(grafted$node.label))
  expect_equal(sum(grafted$node.label == ""), grafted$Nnode - phy_support$Nnode)
  grafted_no_labels <- graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B"), verbose = FALSE, seed = 1)
  expect_null(grafted_no_labels$node.label)
})

## Test that a simmap without allowed states does not restrict grafting, and the updated simmap returned: each missing clade is in the state of its graft point
test_that("graft_missing_clades() with a simmap returns the updated simmap, also without allowed states, where each missing clade is in the state of its graft point",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # B and A are nested sister clades, each with a missing clade of several taxa
  mcl <- list(Taxon_name = c("M1", "M2", "M3", "M4", "M5"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$B, M4 = toy$A, M5 = toy$A))

  for (seed in 1:3)
  {
    without <- graft_missing_clades(toy$phy, mcl, verbose = FALSE, seed = seed)
    # No allowed states: all the states are allowed, same phylogeny as without simmap
    followed <- graft_missing_clades(toy$phy, mcl, simmap = simmap, verbose = FALSE, seed = seed)
    expect_identical(drop_simmap(followed), without)
    sm <- followed$simmap
    expect_s3_class(sm, "simmap")
    expect_identical(sm$edge, followed$edge)
    expect_no_error(validate_simmap(sm, followed))
    expect_equal(vapply(sm$maps, sum, numeric(1)), followed$edge.length)
    # All the grafted branches (of the whole missing clades) are in one state, and the taxa of a missing clade share it
    expect_true(all(lengths(sm$maps[followed$grafted_edges]) == 1))
    tip_state <- function (taxon) names(sm$maps[[match(match(taxon, followed$tip.label), followed$edge[, 2])]])
    expect_identical(tip_state("M1"), tip_state("M2"))
    expect_identical(tip_state("M1"), tip_state("M3"))
    expect_identical(tip_state("M4"), tip_state("M5"))
    # Pruning the grafted taxa gives back the time of the backbone in each state
    pruned <- phytools::drop.tip.simmap(sm, mcl$Taxon_name)
    expect_equal(colSums(pruned$mapped.edge)[colnames(simmap$mapped.edge)], colSums(simmap$mapped.edge))

    # With allowed states, the missing clades are in allowed states
    restricted <- graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B"), verbose = FALSE, seed = seed)
    tip_state_r <- function (taxon) names(restricted$simmap$maps[[match(match(taxon, restricted$tip.label), restricted$edge[, 2])]])
    expect_true(all(vapply(mcl$Taxon_name, tip_state_r, character(1)) %in% c("A", "B")))
  }

  # return_simmap
  # (without allowed states and without the returned simmap, the simmap has no effect: warning, silenced here)
  expect_null(quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, return_simmap = FALSE, seed = 1)$simmap)
  expect_identical(quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, return_simmap = FALSE, seed = 1),
                   drop_simmap(graft_missing_clades(toy$phy, mcl, simmap = simmap, verbose = FALSE, seed = 1)))
  expect_no_warning(graft_missing_clades(toy$phy, mcl, return_simmap = TRUE, verbose = FALSE, seed = 1))
  expect_null(graft_missing_clades(toy$phy, mcl, return_simmap = TRUE, verbose = FALSE, seed = 1)$simmap)
  expect_error(graft_missing_clades(toy$phy, mcl, simmap = simmap, return_simmap = NA, verbose = FALSE), "'return_simmap' must be TRUE or FALSE")
})


# ---- Biogeographic ranges as states (treat_states_as_ranges) ----

## Test that the argument is checked, and requires a simmap
test_that("graft_missing_clades() checks treat_states_as_ranges and requires a simmap with TRUE",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mcl <- toy_missing_clades(toy)
  graft <- function (...) graft_missing_clades(toy$phy, mcl, verbose = FALSE, ...)

  expect_error(graft(simmap = simmap, treat_states_as_ranges = "yes"), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(graft(simmap = simmap, treat_states_as_ranges = NA), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(graft(simmap = simmap, treat_states_as_ranges = c(TRUE, FALSE)), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(graft(treat_states_as_ranges = TRUE), "'treat_states_as_ranges = TRUE' requires a 'simmap'")
  expect_error(graft(treat_states_as_ranges = TRUE, allowed_states = "A"), "requires a 'simmap'")
  expect_no_error(graft(treat_states_as_ranges = FALSE))
})

## Test that single-area states give the same phylogeny and simmap as the exact matching of the states
test_that("graft_missing_clades() with treat_states_as_ranges = TRUE gives the same result as FALSE for single-area states",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mcl <- toy_missing_clades(toy)
  {
    for (seed in 1:2)
    {
      exact <- quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B"), seed = seed)
      areas <- quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, seed = seed)
      expect_identical(areas, exact)
    }
  }
  # Without allowed states
  expect_identical(quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, treat_states_as_ranges = TRUE, seed = 4),
                   quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = simmap, seed = 4))
})

## Test that allowed ranges restrict grafting to the pieces of branches in the states containing them
test_that("graft_missing_clades() with treat_states_as_ranges = TRUE grafts on the states containing the allowed ranges, without warning",
{
  toy <- toy_setup()
  # States of the simmap: A, AB and BC (the toy states A, B and C)
  range_simmap <- make_toy_range_simmap(toy$phy)
  mcl <- toy_missing_clades(toy)
  tip_state <- function (grafted, taxon) names(grafted$simmap$maps[[match(match(taxon, grafted$tip.label), grafted$edge[, 2])]])

  # Range "A" matches the states A and AB
  expect_no_warning(grafted <- graft_missing_clades(toy$phy, mcl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE, seed = 3))
  expect_true(all(vapply(mcl$Taxon_name, function (taxon) tip_state(grafted, taxon), character(1)) %in% c("A", "AB")))

  # A multi-area range alone only matches the states containing it: "AB" matches AB only here
  # (missing clades on the stem branch of the clade B, which is entirely in the state AB)
  mcl_b <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B))
  grafted_ab <- quiet_graft(graft_missing_clades, toy$phy, mcl_b, simmap = range_simmap, allowed_states = "AB", treat_states_as_ranges = TRUE, seed = 3)
  expect_true(all(vapply(mcl_b$Taxon_name, function (taxon) tip_state(grafted_ab, taxon), character(1)) == "AB"))

  # The exact matching of the same simmap: "A" only matches the state A, where some taxa have no branch
  expect_error(graft_missing_clades(toy$phy, mcl, simmap = range_simmap, allowed_states = "A", verbose = FALSE, seed = 3), "no piece of the stem branch of their sister clade in their allowed states")
})

## Test that a missing clade inherits one allowed range for all its taxa and branches
test_that("graft_missing_clades() with treat_states_as_ranges = TRUE gives the whole missing clade one of its allowed ranges, contained in the state of its stem branch",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy) # states A, AB and BC
  mcl <- toy_missing_clades(toy)

  # Allowed range "A": acceptable states A and AB. Same positions as the exact matching of A and AB (a single candidate: no additional random draw)
  areas <- quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, seed = 3)
  exact <- quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = range_simmap, allowed_states = c("A", "AB"), seed = 3)
  expect_identical(drop_simmap(areas), drop_simmap(exact))
  states_areas <- graft_states_in_simmap(areas, mcl$Taxon_name)
  expect_identical(states_areas$state, rep("A", 4))
  expect_true(all(range_in_states_pairwise(states_areas$state, states_areas$parent_state)))
  expect_true(any(graft_states_in_simmap(exact, mcl$Taxon_name)$state == "AB"))

  # M1 and M2 are grafted together as one clade (same sister clade): all the branches of the clade are in the inherited range
  clade_edges <- which(areas$edge[, 2] %in% c(match(c("M1", "M2"), areas$tip.label), ape::getMRCA(areas, c("M1", "M2"))))
  expect_length(clade_edges, 3)
  expect_identical(unique(unlist(lapply(areas$simmap$maps[clade_edges], names))), "A")

  # Several allowed ranges: one range for the whole clade, among the ranges contained in the state of the stem branch
  inherited <- lapply(1:20, function (seed)
  {
    grafted <- quiet_graft(graft_missing_clades, toy$phy, mcl, simmap = range_simmap, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, seed = seed)
    states <- graft_states_in_simmap(grafted, mcl$Taxon_name)
    expect_true(all(range_in_states_pairwise(states$state, states$parent_state)))
    expect_identical(states$state[1], states$state[2]) # the two taxa of the clade
    states$state[1]
  })
  expect_setequal(unique(unlist(inherited)), c("A", "B"))
})

## Test that taxa of the same missing clade must also share the ranges they inherit
test_that("graft_missing_clades() with treat_states_as_ranges = TRUE stops if the taxa of a missing clade have the same acceptable states but inherit different ranges",
{
  toy <- toy_setup()
  # A single state, AB: allowed range "A" gives the acceptable states AB and A (inherited), and so do the allowed ranges "A" and "AB" but with another range to inherit in AB
  one_state <- relabel_simmap(make_toy_simmap(toy$phy), c(A = "AB", B = "AB", C = "AB"))
  one_state$mapped.edge <- NULL
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B), allowed_states = list(M1 = "A", M2 = c("A", "AB")))
  expect_warning(expect_error(graft_missing_clades(toy$phy, mcl, simmap = one_state, treat_states_as_ranges = TRUE, verbose = FALSE),
               "have taxa with different allowed states.*M1 \\(A, AB; ranges inherited: A -> A, AB -> A\\).*M2 \\(A, AB; ranges inherited: A -> A, AB -> A\\|AB\\)"),
                 "not states of the simmap")
  # The same ranges give no error
  mcl$allowed_states <- list(M1 = c("AB", "A"), M2 = c("A", "AB"))
  suppressWarnings(expect_no_error(graft_missing_clades(toy$phy, mcl, simmap = one_state, treat_states_as_ranges = TRUE, verbose = FALSE)))
})

## Test that an allowed range that is not a state of the simmap but is contained in some of them does not stop the grafting, and sends a warning explaining the mismatch
test_that("graft_missing_clades() with treat_states_as_ranges = TRUE grafts with a warning when an allowed range is not a state of the simmap but is contained in some",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy) # states A, AB and BC: "B" is not a state
  mcl <- toy_missing_clades(toy)
  mcl_b <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B)) # the clade B is entirely in the state AB
  expect_warning(grafted <- graft_missing_clades(toy$phy, mcl_b, simmap = range_simmap, allowed_states = "B", treat_states_as_ranges = TRUE, verbose = FALSE, seed = 3),
                 "some allowed ranges are not states of the simmap: 'B' \\(contained in the states 'AB', 'BC'.*The grafting is run assuming that these ranges can be inherited.*subset sympatric speciation")
  states <- graft_states_in_simmap(grafted, mcl_b$Taxon_name)
  expect_identical(states$state, c("B", "B"))
  expect_true(all(range_in_states_pairwise(states$state, states$parent_state)))
  expect_true("B" %in% colnames(grafted$simmap$mapped.edge))

  # With the exact matching, the same allowed state stops with the list of the states available
  expect_error(graft_missing_clades(toy$phy, mcl, simmap = range_simmap, allowed_states = "B", verbose = FALSE), "not in the simmap: 'B'")
  expect_error(graft_missing_clades(toy$phy, mcl, simmap = range_simmap, allowed_states = "D", treat_states_as_ranges = TRUE, verbose = FALSE), "Areas available in the simmap")
})

## Test that the states and the allowed ranges are strictly checked with TRUE only
test_that("graft_missing_clades() with treat_states_as_ranges = TRUE stops if the states or the allowed ranges are not coded as in BioGeoBEARS",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  range_simmap <- make_toy_range_simmap(toy$phy)
  mcl <- toy_missing_clades(toy)
  graft <- function (...) graft_missing_clades(toy$phy, mcl, verbose = FALSE, ...)

  # Incorrect states in the simmap: named in the message, and accepted with the exact matching
  bad_simmap <- relabel_simmap(simmap, c(A = "A", B = "red", C = "CB"))
  expect_error(graft(simmap = bad_simmap, allowed_states = "A", treat_states_as_ranges = TRUE),
               "the states in the simmap must be biogeographic ranges.*'CB' \\(areas not in alphabetical order \\(expected 'BC'\\)\\).*'red' \\(lowercase letters")
  expect_no_error(graft(simmap = bad_simmap, allowed_states = c("A", "red", "CB"), treat_states_as_ranges = FALSE))

  # Incorrect allowed ranges, global and taxon-level
  expect_error(graft(simmap = range_simmap, allowed_states = "ba", treat_states_as_ranges = TRUE), "the global 'allowed_states' argument.*'ba' \\(lowercase letters")
  mcl_bad <- mcl
  mcl_bad$allowed_states <- list(M1 = "AB", M2 = "BA", M3 = NULL, M4 = NULL)
  expect_error(graft_missing_clades(toy$phy, mcl_bad, simmap = range_simmap, treat_states_as_ranges = TRUE, verbose = FALSE),
               "taxon 'M2'.*'BA' \\(areas not in alphabetical order")

  # An allowed area absent from the simmap: the areas available are listed
  expect_error(graft(simmap = range_simmap, allowed_states = c("A", "D"), treat_states_as_ranges = TRUE), "not in the simmap: 'D'.*Areas available in the simmap: 'A', 'B', 'C'")
  # A range contained in no state of the simmap
  expect_error(graft(simmap = range_simmap, allowed_states = "AC", treat_states_as_ranges = TRUE), "No state of the simmap contains the allowed range")
})


# ---- Species attributes: simmap without effect ----

## Test that a simmap that is neither restricting grafting nor returned sends a warning
test_that("graft_missing_clades() warns if the simmap has no effect: no allowed states and return_simmap = FALSE",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mcl <- toy_missing_clades(toy)
  graft <- function (...) graft_missing_clades(toy$phy, mcl, verbose = FALSE, ...)
  has_warning <- function (expr) any(grepl("no effect", testthat::capture_warnings(expr)))

  expect_warning(graft(simmap = simmap, return_simmap = FALSE), "without any allowed states.*'return_simmap = FALSE'.*no effect")
  expect_warning(graft(simmap = simmap, return_simmap = FALSE, treat_states_as_ranges = TRUE), "no effect")
  # The simmap is returned (tracking only), or restricts grafting (global or taxon-level allowed states)
  expect_false(has_warning(graft(simmap = simmap, return_simmap = TRUE)))
  expect_false(has_warning(graft(simmap = simmap, allowed_states = c("A", "B"), return_simmap = FALSE)))
  mcl_states <- mcl
  mcl_states$allowed_states <- list(M1 = c("A", "B"), M2 = c("A", "B"), M3 = NULL, M4 = NULL)
  expect_false(has_warning(graft_missing_clades(toy$phy, mcl_states, simmap = simmap, return_simmap = FALSE, verbose = FALSE)))
  # Without simmap, return_simmap is silently ignored
  expect_false(has_warning(graft(return_simmap = FALSE)))
})
