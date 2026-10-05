
##### Test suite for graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Build toy 'missing_taxa_list' to test graft_missing_clades(): 2 taxa in A (one excluding B), 1 in B, 1 in C, 1 on a single tip
toy_missing_taxa <- function (toy)
{
  list(Taxa_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(toy$A, toy$A, toy$B, toy$C, toy$s),
       Excluded_subclades = list(list(toy$B), NULL, NULL, NULL, NULL))
}


# ---- Input validation ----

## Test validity stops for 'missing_taxa_list' input
test_that("graft_missing_taxa() stops on invalid inputs, before any grafting",
{
  toy <- toy_setup()
  phy <- toy$phy
  mtl <- toy_missing_taxa(toy)

  # Method
  expect_error(graft_missing_taxa(phy, mtl, method = "BD_taxa", verbose = FALSE), "must be either 'uniform' or 'BD_clades'")
  expect_error(graft_missing_taxa(phy, mtl, method = "BD", verbose = FALSE), "must be either 'uniform' or 'BD_clades'")

  # Backbone not ultrametric
  not_ultra <- phy
  not_ultra$edge.length[1] <- not_ultra$edge.length[1] * 1.5
  expect_error(graft_missing_taxa(not_ultra, mtl, verbose = FALSE), "not ultrametric")

  # min_tips, and backbone too small for the birth-death model
  expect_error(graft_missing_taxa(phy, mtl, min_tips = 2, verbose = FALSE), "min_tips")
  expect_error(graft_missing_taxa(make_toy_phy(n_tips = 4), list(Taxa_name = "X1", Target_clade = list(c("t1", "t2"))), verbose = FALSE),
               "too small to fit a birth-death model")

  # Lists
  expect_error(graft_missing_taxa(phy, list(Taxa_name = character(0), Target_clade = list()), verbose = FALSE), "empty")
  expect_error(graft_missing_taxa(phy, list(Taxa_name = c("X1", "X2"), Target_clade = list(toy$A)), verbose = FALSE), "same length as Taxa_name")
  expect_error(graft_missing_taxa(phy, list(Taxa_name = c("X1", "X1"), Target_clade = list(toy$A, toy$A)), verbose = FALSE), "unique and absent")
  expect_error(graft_missing_taxa(phy, list(Taxa_name = toy$s, Target_clade = list(toy$A)), verbose = FALSE), "unique and absent")
  expect_error(graft_missing_taxa(phy, list(Taxa_name = "X1", Target_clade = list(c("zz", "yy"))), verbose = FALSE), "missing from the backbone")
  expect_error(graft_missing_taxa(phy, list(Taxa_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(c("zz", "yy"))), verbose = FALSE),
               "missing from the backbone")

  # Temporary file path
  expect_error(graft_missing_taxa(phy, mtl, temp_file_path = "temp.txt", verbose = FALSE), "must end with '.rds'")
  expect_error(graft_missing_taxa(phy, mtl, temp_file_path = file.path(tempdir(), "no_dir", "t.rds"), verbose = FALSE), "does not exist")
})

## Test that exclusive subclades are valid relative to their target clade
test_that("graft_missing_taxa() checks the position of the excluded subclades relative to the target clade",
{
  toy <- toy_setup()
  one <- function (target, excluded) list(Taxa_name = "X1", Target_clade = list(target), Excluded_subclades = list(list(excluded)))

  # Contains the target clade, or is identical to it: stop, with a pointer to graft_missing_clades()
  expect_error(graft_missing_taxa(toy$phy, one(toy$B, toy$A), verbose = FALSE), "contain the whole target clade")
  for (method in c("uniform", "BD_clades")) for (stem in c(TRUE, FALSE))
  {
    expect_error(graft_missing_taxa(toy$phy, one(toy$B, toy$B), method = method, include_stem = stem, verbose = FALSE),
                 "graft_missing_clades\\(\\)")
  }

  # Disjoint from the target clade: harmless, a warning
  expect_warning(out <- graft_missing_taxa(toy$phy, one(toy$B, toy$D), method = "uniform", verbose = FALSE), "outside the target clade")
  expect_true("X1" %in% out$tip.label)

  # Nested: no warning
  expect_no_warning(graft_missing_taxa(toy$phy, one(toy$A, toy$B), method = "uniform", verbose = FALSE))

  # Not checked if the excluded subclades are not used
  expect_no_warning(graft_missing_taxa(toy$phy, one(toy$B, toy$B), use_excluded_subclades = FALSE, method = "uniform", verbose = FALSE))
})

## Test that exclusive subclades are valid relative to their target clade
test_that("graft_missing_taxa() warns if the target clade is the whole tree and the stem is requested",
{
  toy <- toy_setup()
  whole_tree <- pair_for_node(toy$phy, ape::Ntip(toy$phy) + 1)
  mtl <- list(Taxa_name = "X1", Target_clade = list(whole_tree))

  expect_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", include_stem = TRUE, verbose = FALSE), "whole tree")
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", include_stem = FALSE, verbose = FALSE))
})


# ---- Valid grafting ----

## Test that excluded subclades are respected when grafting
test_that("excluded subclades are respected by both methods",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = paste0("X", 1:4), Target_clade = rep(list(toy$A), 4),
              Excluded_subclades = rep(list(list(toy$B)), 4))

  for (method in c("uniform", "BD_clades"))
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, seed = 3)
    checks <- check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)
    expect_true(all(checks$excluded_ok))
    expect_true(attr(checks, "all_ok"))
  }
})

## Test that grafting constraints are met, using check_grafting_constraints()
test_that("graft_missing_taxa() grafts taxa respecting the constraints with both methods",
{
  toy <- toy_setup()
  mtl <- toy_missing_taxa(toy)

  for (method in c("uniform", "BD_clades"))
  {
    for (seed in 1:2)
    {
      grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, seed = seed)
      expect_valid_graft(grafted, toy$phy, mtl)
      expect_equal(ape::Ntip(grafted), ape::Ntip(toy$phy) + 5)
    }
  }
})

## Test reproducibility with internal seed
test_that("graft_missing_taxa() is reproducible with a seed, and does not alter the global RNG stream",
{
  toy <- toy_setup()
  mtl <- toy_missing_taxa(toy)

  for (method in c("uniform", "BD_clades"))
  {
    a <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, seed = 7)
    b <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, seed = 7)
    d <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, seed = 8)
    expect_equal(a, b)
    expect_false(isTRUE(all.equal(a, d)))
  }

  # The global stream is untouched when a seed is provided
  set.seed(5)
  expected <- stats::runif(1)
  set.seed(5)
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 7)
  expect_identical(stats::runif(1), expected)

  # seed = NULL uses (and advances) the global stream
  set.seed(5)
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = NULL)
  expect_false(identical(stats::runif(1), expected))
})

## Test saving of temporary phylogenies
test_that("graft_missing_taxa() saves temporary phylogenies if requested",
{
  toy <- toy_setup()
  mtl <- toy_missing_taxa(toy)
  temp_file <- file.path(tempdir(), "temp_graft_taxa.rds")
  unlink(temp_file)

  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", save_temp = TRUE, save_temp_every = 1, temp_file_path = temp_file, seed = 1)
  expect_true(file.exists(temp_file))
  expect_s3_class(readRDS(temp_file), "phylo")
  unlink(temp_file)
})

## Test that the output includes grafted tips and branches ID
test_that("graft_missing_taxa() returns the grafted tips and branches, and no sampling fraction table without totals",
{
  toy <- toy_setup()
  mtl <- toy_missing_taxa(toy)
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1)

  expect_setequal(grafted$tip.label[grafted$grafted_tips], mtl$Taxa_name)
  expect_identical(names(grafted$grafted_edges), as.character(grafted$edge[, 2]))
  expect_identical(grafted$grafted_tips, identify_grafted(grafted, mtl$Taxa_name)$grafted_tips)
  expect_null(grafted$sampling_fraction_df)
})


# ---- Stem region ----

## Test that include_stem = TRUE allows multiple missing taxa to be grafting in the stem region
test_that("with include_stem = TRUE, taxa grafted on the stem can receive further grafts",
{
  toy <- toy_setup()
  # Many taxa on a single tip: only its terminal branch is available, so the taxa must be grafted on each other
  mtl <- list(Taxa_name = paste0("X", 1:4), Target_clade = rep(list(toy$s), 4))
  sister_to_missing <- 0

  for (seed in 1:8)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = seed)
    expect_valid_graft(grafted, toy$phy, mtl)
    sister_to_missing <- sister_to_missing + sum(vapply(mtl$Taxa_name, function (taxon) {
      node <- match(taxon, grafted$tip.label)
      parent <- grafted$edge[match(node, grafted$edge[, 2]), 1]
      siblings <- setdiff(grafted$edge[grafted$edge[, 1] == parent, 2], node)
      length(siblings) == 1 && siblings <= ape::Ntip(grafted) && grafted$tip.label[siblings] %in% mtl$Taxa_name
    }, logical(1)))
  }
  expect_gt(sister_to_missing, 0)

  # Same with the birth-death method (one run)
  grafted_bd <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", seed = 1)
  expect_valid_graft(grafted_bd, toy$phy, mtl)
})

## Test that include_stem = FALSE does NOT allow missing taxa to be grafting in the stem region
test_that("with include_stem = FALSE, all taxa are grafted within the crown of their target clade",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = paste0("X", 1:4), Target_clade = rep(list(toy$A), 4), Excluded_subclades = rep(list(NULL), 4))

  for (method in c("uniform", "BD_clades"))
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, include_stem = FALSE, seed = 2)
    checks <- check_grafting_constraints(grafted, toy$phy, mtl, include_stem = FALSE, verbose = FALSE)
    expect_true(attr(checks, "all_ok"))
    expect_false(any(checks$target_position == "stem"))
  }
})


# ---- Nested target clades ----

## Test the order of clades during grafting
test_that("with method = 'uniform', taxa are grafted per target clade from the smallest to the largest",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = paste0("X", 1:7),
              Target_clade = list(toy$A, toy$B, toy$s, toy$A, toy$B, toy$C, toy$s))

  # Record the size of the target clade of each taxon, in the order of grafting
  clade_sizes <- numeric()
  original <- draw_random_branch_weighted_per_length
  local_mocked_bindings(draw_random_branch_weighted_per_length = function (phy, MRCA_node, ...) {
    clade_sizes <<- c(clade_sizes, length(get_descendant_nodes(phy, MRCA_node)))
    original(phy = phy, MRCA_node = MRCA_node, ...)
  })

  for (seed in 1:2)
  {
    clade_sizes <- numeric()
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = seed)
    expect_length(clade_sizes, 7)
    expect_false(is.unsorted(clade_sizes))
    expect_valid_graft(grafted, toy$phy, mtl)
  }
})

## Test that already grafted taxa are accounted for when grafting in higher enclosing clades
test_that("with method = 'BD_clades', nested clades are processed first, and grafted taxa count as present in larger clades",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = c("b1", "b2", "c1", "a1", "a2", "a3"),
              Target_clade = list(toy$B, toy$B, toy$C, toy$A, toy$A, toy$A))

  fits <- list()
  original <- fit_bd_on_clade
  local_mocked_bindings(fit_bd_on_clade = function (...) {
    out <- original(...)
    fits[[length(fits) + 1]] <<- out
    out
  })
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", seed = 1)
  expect_valid_graft(grafted, toy$phy, mtl)

  # One fit per target clade, smaller clades first: B and C (fitted on enclosing clades if small), then A
  expect_length(fits, 3)
  fit_A <- fits[[3]]
  # Fit on A: taxa of B and C (3) are already grafted and counted as present; only the 3 taxa of A are pending
  expect_equal(fit_A$n_tips, toy$nA + 3)
  expect_equal(fit_A$n_pending, 3)
})


# ---- Total richness of the clades ----

## Test that missing taxa already grafted are accounted in sampling fractions when fitting BD models
test_that("total_richness_per_clades is used in the sampling fraction of the birth-death fits, and returned as a table",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = c("b1", "b2", "b3", "c1", "c2", "a1", "a2", "a3", "a4"),
              Target_clade = c(rep(list(toy$B), 3), rep(list(toy$C), 2), rep(list(toy$A), 4)))
  total_B <- toy$nB + 3 + 12
  total_A <- toy$nA + 9 + 31
  totals <- list(Target_clade = list(toy$B, toy$C, toy$A), total_richness = c(total_B, NA, total_A))

  fits <- list()
  original <- fit_bd_on_clade
  local_mocked_bindings(fit_bd_on_clade = function (...) { out <- original(...); fits[[length(fits) + 1]] <<- out; out })

  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, total_richness_per_clades = totals, method = "BD_clades", seed = 2)
  expect_valid_graft(grafted, toy$phy, mtl)

  # Sampling fraction of the fit on A (last): backbone + grafted taxa of B and C are present, 4 taxa of A are pending, 31 unsampled
  fit_A <- fits[[3]]
  expect_equal(fit_A$n_tips, toy$nA + 5)
  expect_equal(fit_A$n_pending, 4)
  expect_equal(fit_A$n_unsampled, 31)
  expect_equal(fit_A$rho, (toy$nA + 5) / total_A)

  # Sampling fraction table, as an element of the phylogeny
  table <- grafted$sampling_fraction_df
  expect_s3_class(table, "data.frame")
  expect_equal(table$nb_unsampled_species, c(12, 0, 31))
  expect_equal(table$total_richness_provided, c(TRUE, FALSE, TRUE))
  expect_equal(table$sampling_fraction, (table$nb_backbone_species + table$nb_missing_species) / table$total_richness)
})

## Test warnings when using 'total_richness_per_clade' with 'uniform' method
test_that("with method = 'uniform', total_richness_per_clades only produces the table, with a warning",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = c("X1", "X2"), Target_clade = list(toy$B, toy$A))
  totals <- list(Target_clade = list(toy$B, toy$A), total_richness = c(toy$nB + 1 + 2, toy$nA + 2 + 5))

  expect_warning(grafted <- graft_missing_taxa(toy$phy, mtl, total_richness_per_clades = totals, method = "uniform", verbose = FALSE),
                 "only used to derive the 'sampling_fraction_df' table")
  expect_s3_class(grafted$sampling_fraction_df, "data.frame")

  # No warning with the birth-death method
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, total_richness_per_clades = totals, method = "BD_clades", verbose = FALSE))
})

## Test validity of 'total_richness_per_clades' input
test_that("graft_missing_taxa() stops before grafting if total_richness_per_clades is invalid",
{
  toy <- toy_setup()
  mtl <- list(Taxa_name = c("X1", "X2"), Target_clade = list(toy$B, toy$A))

  # A target clade is not listed
  expect_error(graft_missing_taxa(toy$phy, mtl, total_richness_per_clades = list(Target_clade = list(toy$B), total_richness = toy$nB + 5), verbose = FALSE),
               "not listed")
  # A total is too low
  expect_error(graft_missing_taxa(toy$phy, mtl, total_richness_per_clades = list(Target_clade = list(toy$B, toy$A), total_richness = c(toy$nB, NA)), verbose = FALSE),
               "lower than the number of backbone species")
})
