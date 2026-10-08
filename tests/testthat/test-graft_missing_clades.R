
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

  grafted <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1, save_temp = TRUE, save_temp_every = 1, temp_file_path = temp_file)
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

## Test that NULL and uninformative total richness give the same fits, but only the latter returns the table
test_that("graft_missing_clades() gives the same result with NULL and with only NA in total_richness_per_clades, but returns the table only with the latter",
{
  toy <- toy_setup()
  mcl <- toy_missing_clades(toy)
  totals_NA <- list(Sister_clade = list(toy$B, toy$C, toy$s), total_richness = rep(NA_real_, 3))

  grafted_null <- quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1)
  grafted_NA <- quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = totals_NA, seed = 1)

  expect_equal(grafted_NA$edge.length, grafted_null$edge.length)
  expect_equal(grafted_NA$tip.label, grafted_null$tip.label)

  expect_null(grafted_null$sampling_fraction_df)
  expect_s3_class(grafted_NA$sampling_fraction_df, "data.frame")
  expect_false(any(grafted_NA$sampling_fraction_df$total_richness_provided))
  expect_equal(grafted_NA$sampling_fraction_df$nb_unsampled_species, rep(0, 3))
  expect_equal(grafted_NA$sampling_fraction_df$sampling_fraction, rep(1, 3))
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

  # Table
  table <- grafted$sampling_fraction_df
  expect_equal(table$nb_missing_species, c(0, 1))
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
