
##### Test suite for graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Build toy 'missing_taxa_list' to test graft_missing_clades(): 2 taxa in A (one excluding B), 1 in B, 1 in C, 1 on a single tip
toy_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
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

  # Backbone with polytomies
  poly <- make_toy_polytomies(phy, n_collapsed = 1)
  expect_error(graft_missing_taxa(poly, mtl, verbose = FALSE), "not binary")
  expect_error(graft_missing_taxa(poly, mtl, verbose = FALSE), "ape::multi2di")

  # min_tips, and backbone too small for the birth-death model
  expect_error(graft_missing_taxa(phy, mtl, min_tips = 2, verbose = FALSE), "min_tips")
  expect_error(graft_missing_taxa(make_toy_phy(n_tips = 4), list(Taxon_name = "X1", Target_clade = list(X1 = c("t1", "t2"))), verbose = FALSE),
               "too small to fit a birth-death model")

  # Lists
  expect_error(graft_missing_taxa(phy, list(Taxon_name = character(0), Target_clade = list()), verbose = FALSE), "empty")
  expect_error(graft_missing_taxa(phy, list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A)), verbose = FALSE), "same length as 'Taxon_name'")
  expect_error(graft_missing_taxa(phy, list(Taxon_name = c("X1", "X1"), Target_clade = list(X1 = toy$A, X1 = toy$A)), verbose = FALSE), "unique names")
  expect_error(graft_missing_taxa(phy, list(Taxon_name = toy$s, Target_clade = stats::setNames(list(toy$A), toy$s)), verbose = FALSE), "absent from the backbone")
  expect_error(graft_missing_taxa(phy, list(Taxon_name = "X1", Target_clade = list(X1 = c("zz", "yy"))), verbose = FALSE), "missing from the backbone")
  expect_error(graft_missing_taxa(phy, list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = c("zz", "yy"))), verbose = FALSE),
               "missing from the backbone")

  # Temporary file path
  expect_error(graft_missing_taxa(phy, mtl, temp_file_path = "temp.txt", verbose = FALSE), "must end with '.rds'")
  expect_error(graft_missing_taxa(phy, mtl, temp_file_path = file.path(tempdir(), "no_dir", "t.rds"), verbose = FALSE), "does not exist")
})

## Test that exclusive subclades are valid relative to their target clade
test_that("graft_missing_taxa() checks the position of the excluded subclades relative to the target clade",
{
  toy <- toy_setup()
  one <- function (target, excluded) list(Taxon_name = "X1", Target_clade = list(X1 = target), Excluded_subclades = list(X1 = list(excluded)))

  # Contains the target clade, or is identical to it: stop, with a pointer to graft_missing_clades()
  expect_error(graft_missing_taxa(toy$phy, one(toy$B, toy$A), verbose = FALSE), "contain the whole target clade")
  for (method in c("uniform", "BD_clades")) for (stem in c(TRUE, FALSE))
  {
    expect_error(graft_missing_taxa(toy$phy, one(toy$B, toy$B), method = method, include_target_clade_stem = stem, verbose = FALSE),
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
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = whole_tree))

  expect_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", include_target_clade_stem = TRUE, verbose = FALSE), "whole tree")
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", include_target_clade_stem = FALSE, verbose = FALSE))
})


# ---- Valid grafting ----

## Test that excluded subclades are respected when grafting
test_that("excluded subclades are respected by both methods",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A, X4 = toy$A),
              Excluded_subclades = list(X1 = list(toy$B), X2 = list(toy$B), X3 = list(toy$B), X4 = list(toy$B)))

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

  expect_setequal(grafted$tip.label[grafted$grafted_tips], mtl$Taxon_name)
  expect_identical(names(grafted$grafted_edges), as.character(grafted$edge[, 2]))
  expect_identical(grafted$grafted_tips, identify_grafted(grafted, mtl$Taxon_name)$grafted_tips)
  expect_null(grafted$sampling_fraction_df)
})


# ---- Stem availability per taxon ----

## Test validity checks for the taxon-level stem availability
test_that("graft_missing_taxa() validates 'include_target_clade_stem'",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A))

  expect_error(graft_missing_taxa(toy$phy, c(mtl, list(include_target_clade_stem = TRUE)), method = "uniform", verbose = FALSE),
                         "same length as 'Taxon_name'")
  expect_error(graft_missing_taxa(toy$phy, c(mtl, list(include_target_clade_stem = c(X1 = 1, X2 = 0))), method = "uniform", verbose = FALSE),
                         "logical vector")
  expect_error(graft_missing_taxa(toy$phy, mtl, include_target_clade_stem = NA, method = "uniform", verbose = FALSE), "single logical value")
})

## Test the warning for single-tip target clades with include_target_clade_stem = FALSE
test_that("graft_missing_taxa() warns once, listing the taxa with a single-tip target clade and include_target_clade_stem = FALSE",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(X1 = toy$s, X2 = toy$s, X3 = toy$A, X4 = toy$A),
              include_target_clade_stem = c(X1 = FALSE, X2 = TRUE, X3 = FALSE, X4 = NA), Excluded_subclades = list(X1 = NULL, X2 = NULL, X3 = NULL, X4 = NULL))
  warning_messages <- function (expr)
  {
    messages <- character()
    withCallingHandlers(expr, warning = function (w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
    messages
  }

  for (method in c("uniform", "BD_clades"))
  {
    messages <- warning_messages(graft_missing_taxa(toy$phy, mtl, method = method, verbose = FALSE))
    single_tip_messages <- messages[grepl("single-tip target clade", messages)]
    expect_length(single_tip_messages, 1) # once, at the beginning
    expect_match(single_tip_messages, "1 taxa")
    expect_match(single_tip_messages, "X1")
    expect_no_match(single_tip_messages, "X2|X3|X4") # only the single-tip taxon with FALSE
    expect_match(single_tip_messages, "ignored")
  }

  # No warning for NA / TRUE, nor for a global include_target_clade_stem = FALSE (documented as ignored for single tips)
  mtl$include_target_clade_stem <- c(X1 = NA, X2 = TRUE, X3 = FALSE, X4 = NA)
  expect_false(any(grepl("single-tip target clade", warning_messages(
    graft_missing_taxa(toy$phy, mtl, method = "uniform", include_target_clade_stem = FALSE, verbose = FALSE)))))

  # The taxon is grafted on the branch of its tip, whatever the value
  mtl$include_target_clade_stem <- c(X1 = FALSE, X2 = TRUE, X3 = NA, X4 = NA)
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1)
  checks <- check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)
  expect_equal(checks$target_position[1:2], c("tip_branch", "tip_branch"))
})

## Test that the taxon-level value overrides the global argument, with the 'uniform' method
test_that("with method = 'uniform', include_target_clade_stem overrides the global include_target_clade_stem per taxon, and NA uses it",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A, X4 = toy$A),
              include_target_clade_stem = c(X1 = TRUE, X2 = FALSE, X3 = NA, X4 = FALSE))

  stem_used <- NULL
  original <- draw_random_grafting_position_along_available_branch_segments
  local_mocked_bindings(draw_random_grafting_position_along_available_branch_segments = function (...) {
      stem_used <<- c(stem_used, list(...)$include_target_clade_stem); original(...) })

  # Global TRUE: NA means TRUE. Taxa without access to the stem come first within their clade
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", include_target_clade_stem = TRUE, seed = 1)
  expect_equal(stem_used[paste0("X", 1:4)], c(X1 = TRUE, X2 = FALSE, X3 = TRUE, X4 = FALSE)) # Names follow the taxa when they are reordered
  expect_equal(unname(stem_used), c(FALSE, FALSE, TRUE, TRUE)) # Taxa without stem first

  # Global FALSE: NA means FALSE
  stem_used <- NULL
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", include_target_clade_stem = FALSE, seed = 1)
  expect_equal(stem_used[paste0("X", 1:4)], c(X1 = TRUE, X2 = FALSE, X3 = FALSE, X4 = FALSE))
  expect_equal(unname(stem_used), c(FALSE, FALSE, FALSE, TRUE))

  # Without the item, the global value applies to all taxa
  stem_used <- NULL
  quiet_graft(graft_missing_taxa, toy$phy, mtl[c("Taxon_name", "Target_clade")], method = "uniform", include_target_clade_stem = FALSE, seed = 1)
  expect_equal(unname(stem_used), rep(FALSE, 4))
})

## Test the groups for the 'BD_clades' method: one fit per clade, taxa without stem first, and the time window of these taxa
test_that("with method = 'BD_clades', taxa of a clade are grouped by stem availability, with one BD fit per clade",
{
  toy <- toy_setup()
  # X1-X4 target A, with a different stem availability; X5 targets B, nested in A
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A, X4 = toy$A, X5 = toy$B),
              include_target_clade_stem = c(X1 = TRUE, X2 = FALSE, X3 = NA, X4 = FALSE, X5 = TRUE), Excluded_subclades = list(X1 = NULL, X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))

  n_fits <- 0
  grafts <- list()
  original_fit <- fit_bd_on_clade
  original_graft <- graft_from_divergence_age
  local_mocked_bindings(
    fit_bd_on_clade = function (...) { n_fits <<- n_fits + 1; original_fit(...) },
    graft_from_divergence_age = function (...) {
      args <- list(...)
      grafts[[length(grafts) + 1]] <<- data.frame(taxon = args$tip_label, include_target_clade_stem = args$include_target_clade_stem, age = args$divergence_age)
      original_graft(...) })

  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", seed = 1)
  grafts <- do.call(rbind, grafts)

  # One fit for B, one for A (although the taxa of A form two groups)
  expect_equal(n_fits, 2)
  # Order: B (nested) first, then the taxa of A without access to the stem, then the taxa of A with the stem
  expect_equal(grafts$taxon[1], "X5")
  expect_setequal(grafts$taxon[2:3], c("X2", "X4"))
  expect_setequal(grafts$taxon[4:5], c("X1", "X3"))
  expect_equal(grafts$include_target_clade_stem, c(TRUE, FALSE, FALSE, TRUE, TRUE))
  # Taxa without the stem have divergence ages within the crown of A
  crown_age_A <- toy$root_age - phytools::nodeheight(toy$phy, toy$node_A)
  expect_true(all(grafts$age[grafts$taxon %in% c("X2", "X4")] <= crown_age_A))

  # The result respects the constraints of each taxon (X2 and X4 within the crown of A), and the item is taken into account by the checks
  expect_valid_graft(grafted, toy$phy, mtl)
  checks <- check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)
  expect_false(any(checks$target_position[checks$Taxon_name %in% c("X2", "X4")] == "stem"))
  expect_equal(checks$include_target_clade_stem, c(TRUE, FALSE, TRUE, FALSE, TRUE))
})

## Test that the richness information does not depend on the stem availability for a shared clade
test_that("clades with different stem availability share the same row of sampling_fraction_df",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A), include_target_clade_stem = c(X1 = TRUE, X2 = FALSE, X3 = NA))
  totals <- list(Target_clade = list(toy$A), total_richness = toy$nA + 3 + 4)
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", total_richness_per_clades = totals, seed = 1)

  expect_equal(nrow(grafted$sampling_fraction_df), 1)
  expect_equal(grafted$sampling_fraction_df$nb_missing_species, 3)
  expect_equal(grafted$sampling_fraction_df$nb_unsampled_species, 4)
})

## Test the whole tree warning with a taxon-level stem availability
test_that("graft_missing_taxa() warns about the stem of the whole tree only for the taxa requesting it",
{
  toy <- toy_setup()
  whole_tree <- pair_for_node(toy$phy, ape::Ntip(toy$phy) + 1)
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = whole_tree, X2 = whole_tree), include_target_clade_stem = c(X1 = FALSE, X2 = TRUE))

  # Messages of all the warnings (the fit of a BD model on the whole tree can issue unrelated optimizer warnings)
  warning_messages <- function (expr)
  {
    messages <- character()
    withCallingHandlers(expr, warning = function (w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
    messages
  }
  mtl_no_stem <- mtl
  mtl_no_stem$include_target_clade_stem <- c(X1 = FALSE, X2 = NA)

  for (method in c("uniform", "BD_clades"))
  {
    expect_true(any(grepl("whole tree", warning_messages(graft_missing_taxa(toy$phy, mtl, method = method, verbose = FALSE)))))
    expect_false(any(grepl("whole tree", warning_messages(
      graft_missing_taxa(toy$phy, mtl_no_stem, method = method, include_target_clade_stem = FALSE, verbose = FALSE)))))
  }
})


# ---- Stem region ----

## Test that include_target_clade_stem = TRUE allows multiple missing taxa to be grafting in the stem region
test_that("with include_target_clade_stem = TRUE, taxa grafted on the stem can receive further grafts",
{
  toy <- toy_setup()
  # Many taxa on a single tip: only its terminal branch is available, so the taxa must be grafted on each other
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(X1 = toy$s, X2 = toy$s, X3 = toy$s, X4 = toy$s))
  sister_to_missing <- 0

  for (seed in 1:8)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = seed)
    expect_valid_graft(grafted, toy$phy, mtl)
    sister_to_missing <- sister_to_missing + sum(vapply(mtl$Taxon_name, function (taxon) {
      node <- match(taxon, grafted$tip.label)
      parent <- grafted$edge[match(node, grafted$edge[, 2]), 1]
      siblings <- setdiff(grafted$edge[grafted$edge[, 1] == parent, 2], node)
      length(siblings) == 1 && siblings <= ape::Ntip(grafted) && grafted$tip.label[siblings] %in% mtl$Taxon_name
    }, logical(1)))
  }
  expect_gt(sister_to_missing, 0)

  # Same with the birth-death method (one run)
  grafted_bd <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", seed = 1)
  expect_valid_graft(grafted_bd, toy$phy, mtl)
})

## Test that include_target_clade_stem = FALSE does NOT allow missing taxa to be grafting in the stem region
test_that("with include_target_clade_stem = FALSE, all taxa are grafted within the crown of their target clade",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A, X4 = toy$A), Excluded_subclades = list(X1 = NULL, X2 = NULL, X3 = NULL, X4 = NULL))

  for (method in c("uniform", "BD_clades"))
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, include_target_clade_stem = FALSE, seed = 2)
    checks <- check_grafting_constraints(grafted, toy$phy, mtl, include_target_clade_stem = FALSE, verbose = FALSE)
    expect_true(attr(checks, "all_ok"))
    expect_false(any(checks$target_position == "stem"))
  }
})


# ---- Stem of the excluded subclades ----

## Identify "cherries" on the toy backbone: a pair of sister tips, used as a target clade with only two terminal branches
# Used to test exclusion from one of its terminal branches => the grafted taxa must land on the remaining sister-branch.
toy_cherry <- function (phy)
{
  n <- ape::Ntip(phy)
  parents <- phy$edge[phy$edge[, 2] <= n, 1]
  node <- unique(parents[duplicated(parents)])[1]
  phy$tip.label[phy$edge[phy$edge[, 1] == node, 2]]
}

## Concatenate as a message all warnings raised when evaluating an expression
warning_messages <- function (expr)
{
  messages <- character()
  withCallingHandlers(expr, warning = function (w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
  messages
}

## Test the validation of the stem exclusion of the excluded subclades
test_that("graft_missing_taxa() validates 'exclude_subclades_stem'",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A), Excluded_subclades = list(X1 = list(toy$B, toy$C), X2 = NULL))

  expect_error(graft_missing_taxa(toy$phy, mtl, exclude_subclades_stem = NA, method = "uniform", verbose = FALSE), "single logical value")
  expect_error(graft_missing_taxa(toy$phy, c(mtl, list(exclude_subclades_stem = list(TRUE))), method = "uniform", verbose = FALSE),
               "same length as 'Taxon_name'")
  expect_error(graft_missing_taxa(toy$phy, c(mtl, list(exclude_subclades_stem = list(X1 = c(TRUE, FALSE, TRUE), X2 = NULL))), method = "uniform", verbose = FALSE),
               "excluded subclades")
  expect_error(graft_missing_taxa(toy$phy, c(mtl, list(exclude_subclades_stem = list(X1 = c(1, 0), X2 = NULL))), method = "uniform", verbose = FALSE),
               "logical values")
})

## Test all combinations of branch exclusions in the uniform method
test_that("with method = 'uniform', the stem of an excluded subclade is excluded by default, with the setting of each subclade",
{
  toy <- toy_setup()
  desc <- function (node) get_descendant_nodes(toy$phy, node)
  excluded_used <- list()
  original <- draw_random_grafting_position_along_available_branch_segments
  local_mocked_bindings(draw_random_grafting_position_along_available_branch_segments = function (...) {
    excluded_used[[length(excluded_used) + 1]] <<- list(...)$excluded_nodes; original(...) })

  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = list(toy$B, toy$C)))

  # Default: the stem is excluded for both subclades
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1)
  expect_setequal(excluded_used[[1]], c(desc(toy$node_B), toy$node_B, desc(toy$node_C), toy$node_C))

  # Global FALSE: only the branches within the subclades
  excluded_used <- list()
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", exclude_subclades_stem = FALSE, seed = 1)
  expect_setequal(excluded_used[[1]], c(desc(toy$node_B), desc(toy$node_C)))

  # Setting specific to each subclade, with NA using the global value
  excluded_used <- list()
  quiet_graft(graft_missing_taxa, toy$phy, c(mtl, list(exclude_subclades_stem = list(X1 = c(TRUE, NA)))), method = "uniform", exclude_subclades_stem = FALSE, seed = 1)
  expect_setequal(excluded_used[[1]], c(desc(toy$node_B), toy$node_B, desc(toy$node_C)))

  # A single value for all the subclades of the taxon, or no exclusion at all
  excluded_used <- list()
  quiet_graft(graft_missing_taxa, toy$phy, c(mtl, list(exclude_subclades_stem = c(X1 = FALSE))), method = "uniform", seed = 1)
  expect_setequal(excluded_used[[1]], c(desc(toy$node_B), desc(toy$node_C)))
  excluded_used <- list()
  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", use_excluded_subclades = FALSE, seed = 1)
  expect_length(excluded_used, 0) # nothing excluded: NULL passed to the draw
})

## Test the stem exclusion passed to each taxon in the 'BD_clades' method
test_that("with method = 'BD_clades', the stem exclusion of each excluded subclade is used for each taxon",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A), Excluded_subclades = list(X1 = list(toy$B, toy$C), X2 = list(toy$B)),
              exclude_subclades_stem = list(X1 = c(TRUE, FALSE), X2 = NA))
  used <- list()
  original <- graft_from_divergence_age
  local_mocked_bindings(graft_from_divergence_age = function (...) {
    args <- list(...); used[[args$tip_label]] <<- args$exclude_subclades_stem; original(...) })

  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", seed = 1)
  expect_equal(used$X1, c(TRUE, FALSE))
  expect_equal(used$X2, TRUE) # NA uses the global value
  expect_valid_graft(grafted, toy$phy, mtl)
})

## Test that a single tip can be excluded: its terminal branch is not available
test_that("a single-tip excluded subclade excludes its terminal branch by default, with both methods",
{
  toy <- toy_setup()
  cherry <- toy_cherry(toy$phy)
  # Target clade with only two terminal branches (no stem), one of them excluded: the only branch left is the other one
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = cherry), Excluded_subclades = list(X1 = list(cherry[1])))

  for (method in c("uniform", "BD_clades"))
  {
    for (seed in 1:3)
    {
      grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, include_target_clade_stem = FALSE, seed = seed)
      sister_tips <- grafted$tip.label[phytools::getDescendants(grafted, ape::getMRCA(grafted, c("X1", cherry[2])))]
      sister_tips <- sister_tips[!is.na(sister_tips)]
      expect_setequal(sister_tips, c("X1", cherry[2]))
      checks <- check_grafting_constraints(grafted, toy$phy, mtl, include_target_clade_stem = FALSE, verbose = FALSE)
      expect_true(checks$excluded_ok)
    }
  }

  # The same exclusion is an error when the terminal branch is the only one left to exclude: target and excluded tip are identical
  same_tip <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$s), Excluded_subclades = list(X1 = list(toy$s)))
  expect_error(graft_missing_taxa(toy$phy, same_tip, method = "uniform", verbose = FALSE), "graft_missing_clades\\(\\)")
  # ... but not if the exclusion has no effect
  expect_true("X1" %in% quiet_graft(graft_missing_taxa, toy$phy, same_tip, method = "uniform", exclude_subclades_stem = FALSE, seed = 1)$tip.label)
})

## Test the warning for single-tip excluded subclades that exclude nothing
test_that("graft_missing_taxa() warns once, listing the single-tip excluded subclades with exclude_subclades_stem = FALSE",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A),
              Excluded_subclades = list(X1 = list(toy$s, toy$B), X2 = list(toy$B), X3 = NULL),
              exclude_subclades_stem = list(X1 = c(FALSE, TRUE), X2 = FALSE, X3 = NULL))

  for (method in c("uniform", "BD_clades"))
  {
    messages <- warning_messages(graft_missing_taxa(toy$phy, mtl, method = method, verbose = FALSE))
    noop_messages <- messages[grepl("single-tip excluded subclade", messages)]
    expect_length(noop_messages, 1) # once, at the beginning
    expect_match(noop_messages, "1 single-tip")
    expect_match(noop_messages, "X1")
    expect_match(noop_messages, toy$s, fixed = TRUE)
    expect_no_match(noop_messages, "X2|X3") # X2 has a true clade with FALSE: the exclusion has an effect
  }

  # A global FALSE also warns (single tips then exclude nothing), and the default or TRUE values do not
  mtl$exclude_subclades_stem <- NULL
  expect_true(any(grepl("single-tip excluded subclade", warning_messages(
    graft_missing_taxa(toy$phy, mtl, method = "uniform", exclude_subclades_stem = FALSE, verbose = FALSE)))))
  expect_false(any(grepl("single-tip excluded subclade", warning_messages(
    graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE)))))
  # No warning if the excluded subclades are not used
  expect_false(any(grepl("single-tip excluded subclade", warning_messages(
    graft_missing_taxa(toy$phy, mtl, method = "uniform", exclude_subclades_stem = FALSE, use_excluded_subclades = FALSE, verbose = FALSE)))))
})


# ---- Nested target clades ----

## Test the order of clades during grafting
test_that("with method = 'uniform', taxa are grafted per target clade from the smallest to the largest",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4", "X5", "X6", "X7"),
              Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$s, X4 = toy$A, X5 = toy$B, X6 = toy$C, X7 = toy$s))

  # Record the size of the target clade of each taxon, in the order of grafting
  clade_sizes <- numeric()
  original <- draw_random_grafting_position_along_available_branch_segments
  local_mocked_bindings(draw_random_grafting_position_along_available_branch_segments = function (phy, MRCA_node, ...) {
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
  mtl <- list(Taxon_name = c("b1", "b2", "c1", "a1", "a2", "a3"),
              Target_clade = list(b1 = toy$B, b2 = toy$B, c1 = toy$C, a1 = toy$A, a2 = toy$A, a3 = toy$A))

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
  mtl <- list(Taxon_name = c("b1", "b2", "b3", "c1", "c2", "a1", "a2", "a3", "a4"),
              Target_clade = list(b1 = toy$B, b2 = toy$B, b3 = toy$B, c1 = toy$C, c2 = toy$C,
                                  a1 = toy$A, a2 = toy$A, a3 = toy$A, a4 = toy$A))
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
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$B, X2 = toy$A))
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
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$B, X2 = toy$A))

  # A target clade is not listed
  expect_error(graft_missing_taxa(toy$phy, mtl, total_richness_per_clades = list(Target_clade = list(toy$B), total_richness = toy$nB + 5), verbose = FALSE),
               "not listed")
  # A total is too low
  expect_error(graft_missing_taxa(toy$phy, mtl, total_richness_per_clades = list(Target_clade = list(toy$B, toy$A), total_richness = c(toy$nB, NA)), verbose = FALSE),
               "lower than the number of backbone species")
})


## Test that every item of missing_taxa_list must be named after Taxon_name
test_that("graft_missing_taxa() requires the items of missing_taxa_list to be named after Taxon_name (same names, same order)",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$B), include_target_clade_stem = c(X1 = TRUE, X2 = NA),
                         Excluded_subclades = list(X1 = list(toy$B), X2 = NULL), exclude_subclades_stem = list(X1 = FALSE, X2 = NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1)
  expect_valid_graft(grafted, toy$phy, mtl)

  for (item in c("Target_clade", "include_target_clade_stem", "Excluded_subclades", "exclude_subclades_stem"))
  {
    unnamed <- mtl; names(unnamed[[item]]) <- NULL
    expect_error(graft_missing_taxa(toy$phy, unnamed, method = "uniform", verbose = FALSE), paste0(item, ".*must be named after 'Taxon_name'"))

    wrong <- mtl; names(wrong[[item]]) <- c("X1", "X3")
    expect_error(graft_missing_taxa(toy$phy, wrong, method = "uniform", verbose = FALSE), paste0(item, ".*do not match.*'X3'"))

    # Same names in a different order: items matched by position would be assigned to the wrong taxa
    swapped <- mtl; swapped[[item]] <- rev(swapped[[item]])
    expect_error(graft_missing_taxa(toy$phy, swapped, method = "uniform", verbose = FALSE), paste0(item, ".*not in the same order"))
  }

  # Optional items that are NULL (or absent) do not need names
  expect_no_error(quiet_graft(graft_missing_taxa, toy$phy, mtl[c("Taxon_name", "Target_clade")], method = "uniform", seed = 1))
  mtl$include_target_clade_stem <- NULL
  expect_no_error(quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1))

  # The names are not required to be used by the algorithm: only the order matters once validated, so reordering all the items together gives a valid graft
  reordered <- list(Taxon_name = c("X2", "X1"), Target_clade = list(X2 = toy$B, X1 = toy$A), Excluded_subclades = list(X2 = NULL, X1 = list(toy$B)))
  expect_valid_graft(quiet_graft(graft_missing_taxa, toy$phy, reordered, method = "uniform", seed = 1), toy$phy, reordered)
})

## Test that the optional names of the excluded subclades are used in the warnings
test_that("graft_missing_taxa() uses the names of the excluded subclades in its warnings",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A),
                         Excluded_subclades = list(X1 = list(tip_s = toy$s, clade_B = toy$B), X2 = list(toy$D)),
                         exclude_subclades_stem = list(X1 = c(FALSE, TRUE), X2 = NULL))

  # Single-tip excluded subclade that excludes nothing
  messages <- warning_messages(graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE))
  noop <- messages[grepl("single-tip excluded subclade", messages)]
  expect_length(noop, 1)
  expect_match(noop, paste0("X1 (excluded subclade: 'tip_s' (", toy$s, "))"), fixed = TRUE)

  # Excluded subclade outside of the target clade (named or not)
  mtl$Excluded_subclades$X1 <- list(clade_D = toy$D)
  mtl$exclude_subclades_stem <- NULL
  messages <- warning_messages(graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE))
  outside <- messages[grepl("outside the target clade", messages)]
  expect_length(outside, 1)
  expect_match(outside, "X1 (target clade tips:", fixed = TRUE)
  expect_match(outside, "excluded subclade: 'clade_D' (", fixed = TRUE)
  expect_match(outside, paste0("X2 (target clade tips:"), fixed = TRUE)
  expect_match(outside, paste0("excluded subclade: ", paste(toy$D, collapse = ", ")), fixed = TRUE) # unnamed: tips only
})

# ---- Node labels ----

## Test that node labels of the backbone are kept, with an empty string for the nodes created by grafting, and no labels if the backbone had none
test_that("graft_missing_taxa() keeps the node labels of the backbone and gives an empty string to the new nodes, for both methods",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$A))
  phy_labelled <- toy$phy
  phy_labelled$node.label <- paste0("support", seq_len(phy_labelled$Nnode))

  for (method in c("uniform", "BD_clades"))
  {
    # Backbone without labels
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = method, seed = 1)
    expect_null(grafted$node.label)

    # Backbone with labels: same tree as without labels, original labels kept, new nodes labelled with an empty string
    grafted_labelled <- quiet_graft(graft_missing_taxa, phy_labelled, mtl, method = method, seed = 1)
    expect_equal(grafted_labelled$edge, grafted$edge)
    expect_equal(grafted_labelled$edge.length, grafted$edge.length)
    expect_length(grafted_labelled$node.label, grafted_labelled$Nnode)
    expect_equal(sum(grafted_labelled$node.label == ""), 3)
    expect_setequal(grafted_labelled$node.label[grafted_labelled$node.label != ""], phy_labelled$node.label)
    expect_false(anyNA(grafted_labelled$node.label))
  }
})

## Test that the temporary phylogenies saved during grafting do not contain the labels used to follow the nodes
test_that("graft_missing_taxa() saves temporary phylogenies with the original node labels",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A))
  phy_labelled <- toy$phy
  phy_labelled$node.label <- paste0("support", seq_len(phy_labelled$Nnode))
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  for (method in c("uniform", "BD_clades"))
  {
    quiet_graft(graft_missing_taxa, phy_labelled, mtl, method = method, seed = 1, save_temp = TRUE, save_temp_every = 1, temp_file_path = temp_file)
    saved <- readRDS(temp_file)
    expect_false(any(grepl("^bb", saved$node.label)))
    expect_setequal(saved$node.label[saved$node.label != ""], phy_labelled$node.label)
    unlink(temp_file)
  }
})


# ---- Species attributes: allowed states of a simmap ----

## Test that all the states allowed gives exactly the same phylogeny as without simmap, with global or taxon-level states
test_that("graft_missing_taxa() with all the states allowed gives the same phylogeny as without simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  without <- graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, seed = 5)

  # Global states
  with_all <- graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 5)
  expect_identical(drop_simmap(with_all), without)

  # Taxon-level states, in a different order, with the global states for some taxa
  mtl$allowed_states <- list(X1 = c("C", "A", "B"), X2 = NULL, X3 = c("B", "A", "C"), X4 = NA, X5 = NULL)
  with_taxa <- graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 5)
  expect_identical(drop_simmap(with_taxa), without)

  # The simmap does not need to have the same order of tips and branches as the backbone, and does not change the result
  expect_identical(drop_simmap(graft_missing_taxa(toy$phy, toy_missing_taxa(toy), simmap = simmap, allowed_states = c("A", "B", "C"), method = "uniform",
                                                  verbose = FALSE, seed = 5)), without)
})

## Test that taxa are only grafted on pieces of branches in their allowed states, with global and taxon-level states, and that the constraints are still respected
test_that("graft_missing_taxa() grafts the taxa in their allowed states only, global or taxon-level, and respects the other constraints",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  taxa <- c(paste0("X", 1:8), "Y1", "Y2", "Y3")
  mtl <- list(Taxon_name = taxa,
              Target_clade = c(rep(list(toy$A), 8), rep(list(toy$C), 3)),
              Excluded_subclades = c(rep(list(list(toy$B)), 8), rep(list(NULL), 3)),
              allowed_states = c(rep(list(NULL), 8), rep(list("C"), 3)))
  names(mtl$Target_clade) <- names(mtl$Excluded_subclades) <- names(mtl$allowed_states) <- taxa

  n_checked <- 0
  for (seed in 1:4)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = seed)
    expect_valid_graft(grafted, toy$phy, mtl)
    for (taxon in taxa)
    {
      state <- state_at_graft(grafted, simmap, taxon, taxa)
      if (is.na(state)) { next }
      n_checked <- n_checked + 1
      expect_true(state %in% (if (taxon %in% c("Y1", "Y2", "Y3")) "C" else c("A", "B")), info = paste(taxon, "seed", seed))
    }
  }
  # The states of many taxa could be checked independently (the others were grafted on the branches of other grafted taxa)
  expect_gt(n_checked, 10)

  # Restricting the states changes the result
  unrestricted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B", "C"), method = "uniform", seed = 1)
  expect_false(identical(unrestricted$edge.length, grafted$edge.length))

  # Same seed, same result
  again <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = 4)
  expect_identical(again, grafted)
})

## Test that grafted taxa have the state of the point where they were grafted, so taxa with other allowed states cannot be grafted on their branches
test_that("graft_missing_taxa() does not graft taxa on the branches of grafted taxa in other states",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # The terminal branch of the tip 's' is entirely in state "A": the taxa grafted there are in state "A"
  expect_identical(names(simmap$maps[[match(match(toy$s, toy$phy$tip.label), toy$phy$edge[, 2])]]), "A")
  tip_taxa <- paste0("S", 1:5)
  clade_taxa <- paste0("Y", 1:10)
  taxa <- c(tip_taxa, clade_taxa)
  mtl <- list(Taxon_name = taxa,
              Target_clade = c(rep(list(toy$s), 5), rep(list(toy$A), 10)),
              allowed_states = c(rep(list("A"), 5), rep(list("B"), 10)))
  names(mtl$Target_clade) <- names(mtl$allowed_states) <- taxa

  for (seed in 1:5)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, method = "uniform", seed = seed)
    # The S taxa are all grafted on the terminal branch of the tip 's': together with 's', they form a subtree (its stem is the branch of 's' above the first graft)
    s_subtree <- phytools::getDescendants(grafted, ape::getMRCA(grafted, c(toy$s, tip_taxa)))
    s_subtree_tips <- grafted$tip.label[s_subtree[s_subtree <= ape::Ntip(grafted)]]
    expect_true(all(tip_taxa %in% s_subtree_tips))
    # The Y taxa (state "B" only) are never grafted in this subtree, which is entirely in state "A"
    expect_false(any(clade_taxa %in% s_subtree_tips), info = paste("seed", seed))
    # When they are grafted on a branch of the backbone, it is in state "B"
    states <- vapply(clade_taxa, function (taxon) state_at_graft(grafted, simmap, taxon, taxa), character(1))
    expect_true(all(states[!is.na(states)] == "B"), info = paste("seed", seed))
  }
})

## Test the validity checks of the species attributes arguments
test_that("graft_missing_taxa() stops or warns on invalid 'simmap', 'allowed_states' and 'min_valid_fraction' inputs, before any grafting",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  graft <- function (...) graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, ...)

  # A simmap does not require allowed states (all the states are allowed), also for the taxa without a taxon-level value
  expect_no_error(graft(simmap = simmap))
  mtl_some <- mtl
  mtl_some$allowed_states <- list(X1 = "A", X2 = "A", X3 = NULL, X4 = "B", X5 = NA)
  expect_no_error(graft_missing_taxa(toy$phy, mtl_some, simmap = simmap, method = "uniform", verbose = FALSE))

  # Unknown states: global or taxon-level, listing the states of the simmap
  expect_error(graft(simmap = simmap, allowed_states = c("A", "red")), "global 'allowed_states' argument.*'red'.*'A', 'B', 'C'")
  mtl_unknown <- mtl
  mtl_unknown$allowed_states <- list(X1 = "red", X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL)
  expect_error(graft_missing_taxa(toy$phy, mtl_unknown, simmap = simmap, allowed_states = "A", method = "uniform", verbose = FALSE), "taxon 'X1'.*'red'")

  # Allowed states without simmap
  expect_error(graft(allowed_states = "A"), "requires a 'simmap'")
  expect_error(graft_missing_taxa(toy$phy, mtl_some, method = "uniform", verbose = FALSE), "requires a 'simmap'")
  # ... but a list without any state is ignored
  mtl_empty <- mtl
  mtl_empty$allowed_states <- list(X1 = NULL, X2 = NULL, X3 = NA, X4 = NULL, X5 = NULL)
  expect_no_error(graft_missing_taxa(toy$phy, mtl_empty, method = "uniform", verbose = FALSE))

  # The taxon-level item must be named after Taxon_name, like the other items
  mtl_unnamed <- mtl_some
  names(mtl_unnamed$allowed_states) <- NULL
  expect_error(graft_missing_taxa(toy$phy, mtl_unnamed, simmap = simmap, method = "uniform", verbose = FALSE), "must be named after 'Taxon_name'")

  # A plain list of maps is the maps of several traits (see test-multiple_traits.R), and a simmap of another tree
  expect_error(graft(simmap = list(simmap, simmap), allowed_states = "A"), "must be given for each trait")
  other <- make_toy_simmap(make_toy_phy(n_tips = 20, seed = 99))
  expect_error(graft(simmap = other, allowed_states = "A"), "tips|topology|branch lengths")

  # max_age_redraws
  for (value in list(0, -1, 1.5, NA, "10", c(10, 20)))
  {
    expect_error(graft(simmap = simmap, allowed_states = "A", max_age_redraws = value), "'max_age_redraws' must be a single positive integer")
  }

  # min_valid_fraction
  for (value in list(-0.1, 1.5, NA, "0.1", c(0.1, 0.2)))
  {
    expect_error(graft(simmap = simmap, allowed_states = "A", min_valid_fraction = value), "'min_valid_fraction' must be a single number between 0 and 1")
  }
})

## Test the stop before grafting if some taxa have no piece of branch in their allowed states, and the warning for the taxa with very few
test_that("graft_missing_taxa() stops if a taxon has no branch in its allowed states, and warns if the fraction in allowed states is below 'min_valid_fraction'",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)

  # The clade B is entirely in state "B": a taxon that can only be grafted in other states has no branch
  mtl <- list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$B, X2 = toy$B, X3 = toy$s),
              allowed_states = list(X1 = c("A", "C"), X2 = "B", X3 = "B"))
  expect_error(graft_missing_taxa(toy$phy, mtl, simmap = simmap, method = "uniform", verbose = FALSE),
               "2 taxa have no branch available for grafting in their allowed states.*X1 .*allowed states: A, C.*X3 .*allowed states: B")
  # Not flagged with other allowed states
  mtl$allowed_states <- list(X1 = "B", X2 = "B", X3 = "A")
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, simmap = simmap, method = "uniform", verbose = FALSE))

  # Small fraction of the available branches in the allowed states: warning, listing the taxa with their fraction
  mtl_low <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$s), allowed_states = list(X1 = "C", X2 = "A"))
  expect_warning(out <- graft_missing_taxa(toy$phy, mtl_low, simmap = simmap, method = "uniform", verbose = FALSE),
                 "1 taxa have less than 10% of the length of their available branches in their allowed states.*X1 .*% in allowed states")
  expect_true(all(c("X1", "X2") %in% out$tip.label))
  # A lower threshold, or all states allowed: no warning
  expect_no_warning(graft_missing_taxa(toy$phy, mtl_low, simmap = simmap, min_valid_fraction = 0.01, method = "uniform", verbose = FALSE))
  mtl_global <- mtl_low
  mtl_global$allowed_states <- NULL
  expect_no_warning(graft_missing_taxa(toy$phy, mtl_global, simmap = simmap, allowed_states = c("A", "B", "C"), min_valid_fraction = 0.99, method = "uniform", verbose = FALSE))
  # A threshold of 1 warns for every taxon with some branches outside the allowed states
  expect_warning(graft_missing_taxa(toy$phy, mtl_low, simmap = simmap, min_valid_fraction = 1, method = "uniform", verbose = FALSE), "1 taxa have less than 100%")
})

## Test the node labels of the grafted phylogeny with a simmap: no labels if the backbone had none, an empty string for the new nodes otherwise
test_that("graft_missing_taxa() with a simmap gives the same node labels as without",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)

  grafted <- graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", verbose = FALSE, seed = 1)
  expect_null(grafted$node.label)

  phy_support <- toy$phy
  phy_support$node.label <- paste0("support", seq_len(phy_support$Nnode))
  grafted <- graft_missing_taxa(phy_support, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", verbose = FALSE, seed = 1)
  expect_length(grafted$node.label, grafted$Nnode)
  expect_equal(sum(grafted$node.label == ""), length(mtl$Taxon_name))
  expect_setequal(grafted$node.label[grafted$node.label != ""], phy_support$node.label)
})


# ---- Biogeographic ranges as states (treat_states_as_ranges) ----

## Test that the argument is checked, and requires a simmap
test_that("graft_missing_taxa() checks treat_states_as_ranges and requires a simmap with TRUE",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  graft <- function (...) graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, ...)

  expect_error(graft(simmap = simmap, treat_states_as_ranges = "yes"), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(graft(simmap = simmap, treat_states_as_ranges = NA), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(graft(simmap = simmap, treat_states_as_ranges = c(TRUE, FALSE)), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(graft(treat_states_as_ranges = TRUE), "'treat_states_as_ranges = TRUE' requires a 'simmap'")
  expect_error(graft(treat_states_as_ranges = TRUE, allowed_states = "A"), "requires a 'simmap'")
  expect_no_error(graft(treat_states_as_ranges = FALSE))
})

## Test that single-area states give the same phylogeny and simmap as the exact matching of the states
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE gives the same result as FALSE for single-area states",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  for (method in c("uniform", "BD_clades"))
  {
    for (seed in 1:2)
    {
      exact <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = method, seed = seed)
      areas <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, method = method, seed = seed)
      expect_identical(areas, exact)
    }
  }
  # Without allowed states
  expect_identical(quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, treat_states_as_ranges = TRUE, method = "uniform", seed = 4),
                   quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, method = "uniform", seed = 4))
})

## Test that allowed ranges restrict grafting to the pieces of branches in the states containing them
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE grafts on the states containing the allowed ranges, without warning",
{
  toy <- toy_setup()
  # States of the simmap: A, AB and BC (the toy states A, B and C)
  range_simmap <- make_toy_range_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  tip_state <- function (grafted, taxon) names(grafted$simmap$maps[[match(match(taxon, grafted$tip.label), grafted$edge[, 2])]])

  # Range "A" matches the states A and AB
  expect_no_warning(grafted <- graft_missing_taxa(toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, method = "uniform", verbose = FALSE, seed = 3))
  expect_true(all(vapply(mtl$Taxon_name, function (taxon) tip_state(grafted, taxon), character(1)) %in% c("A", "AB")))

  # A multi-area range alone only matches the states containing it: "AB" matches AB only here
  # (taxa of the clade B, which is entirely in the state AB)
  mtl_b <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$B, X2 = toy$B))
  grafted_ab <- quiet_graft(graft_missing_taxa, toy$phy, mtl_b, simmap = range_simmap, allowed_states = "AB", treat_states_as_ranges = TRUE, method = "uniform", seed = 3)
  expect_true(all(vapply(mtl_b$Taxon_name, function (taxon) tip_state(grafted_ab, taxon), character(1)) == "AB"))

  # The exact matching of the same simmap: "A" only matches the state A, where some taxa have no branch
  expect_error(graft_missing_taxa(toy$phy, mtl, simmap = range_simmap, allowed_states = "A", method = "uniform", verbose = FALSE, seed = 3), "no branch available for grafting in their allowed states")
})

## Test that the taxa inherit an allowed range, contained in the state of the branch they were grafted on, and that the grafting positions only depend on the acceptable states
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE gives the taxa one of their allowed ranges, contained in the state of their parent branch",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy) # states A, AB and BC
  mtl <- toy_missing_taxa(toy)

  # Allowed range "A": acceptable states A and AB. Same positions as the exact matching of A and AB (a single candidate: no additional random draw)
  areas <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, method = "uniform", seed = 3)
  exact <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = c("A", "AB"), method = "uniform", seed = 3)
  expect_identical(drop_simmap(areas), drop_simmap(exact))
  states_areas <- graft_states_in_simmap(areas, mtl$Taxon_name)
  states_exact <- graft_states_in_simmap(exact, mtl$Taxon_name)
  expect_identical(states_areas$state, rep("A", length(mtl$Taxon_name)))
  # The exact matching gives the state of the branch: AB for the taxa grafted on a branch in AB
  expect_identical(states_exact$state, states_exact$parent_state)
  expect_true(any(states_exact$state == "AB"))
  expect_true(all(range_in_states_pairwise(states_areas$state, states_areas$parent_state)))

  # Allowed ranges "A" and "B": the range inherited is drawn among the ranges contained in the state of the branch
  for (seed in 1:3)
  {
    both <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, method = "uniform", seed = seed)
    states_both <- graft_states_in_simmap(both, mtl$Taxon_name)
    expect_true(all(states_both$state %in% c("A", "B")))
    expect_true(all(range_in_states_pairwise(states_both$state, states_both$parent_state)))
  }

  # Allowed range "AB" alone: inherited as such, from the states AB only here
  mtl_b <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$B, X2 = toy$B))
  grafted_ab <- quiet_graft(graft_missing_taxa, toy$phy, mtl_b, simmap = range_simmap, allowed_states = "AB", treat_states_as_ranges = TRUE, method = "uniform", seed = 3)
  expect_identical(graft_states_in_simmap(grafted_ab, mtl_b$Taxon_name)$state, c("AB", "AB"))

  # Both methods, and a taxon-level value
  mtl_taxon <- mtl
  mtl_taxon$allowed_states <- list(X1 = "A", X2 = c("A", "B"), X3 = NULL, X4 = NULL, X5 = NULL)
  for (method in c("uniform", "BD_clades"))
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl_taxon, simmap = range_simmap, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, method = method, seed = 2)
    states_taxon <- graft_states_in_simmap(grafted, mtl$Taxon_name)
    expect_true(all(range_in_states_pairwise(states_taxon$state, states_taxon$parent_state)))
    expect_identical(states_taxon$state[1], "A")
    expect_true(states_taxon$state[2] %in% c("A", "B"))
    expect_true(ape::is.ultrametric(grafted))
    expect_equal(unname(rowSums(grafted$simmap$mapped.edge)), grafted$edge.length, tolerance = 1e-6)
  }
})

## Test that both candidate ranges are drawn when several allowed ranges are contained in the state of the branch
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE draws among the ranges contained in the state of the branch",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy)
  mtl_b <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$B)) # the clade B is entirely in the state AB
  inherited <- vapply(1:30, function (seed)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl_b, simmap = range_simmap, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, method = "uniform", seed = seed)
    graft_states_in_simmap(grafted, "X1")$state
  }, character(1))
  expect_setequal(unique(inherited), c("A", "B"))
  expect_true(all(table(inherited) >= 5))
})

## Test that an allowed range that is not a state of the simmap but is contained in some of them does not stop the grafting, and sends a warning explaining the mismatch
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE grafts with a warning when an allowed range is not a state of the simmap but is contained in some",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy) # states A, AB and BC: "B" is not a state
  mtl <- toy_missing_taxa(toy)
  mtl_b <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$B, X2 = toy$B)) # the clade B is entirely in the state AB
  expect_warning(grafted <- graft_missing_taxa(toy$phy, mtl_b, simmap = range_simmap, allowed_states = "B", treat_states_as_ranges = TRUE, method = "uniform", verbose = FALSE, seed = 3),
                 "some allowed ranges are not states of the simmap: 'B' \\(contained in the states 'AB', 'BC'.*The grafting is run assuming that these ranges can be inherited.*subset sympatric speciation")
  states <- graft_states_in_simmap(grafted, mtl_b$Taxon_name)
  expect_identical(states$state, c("B", "B"))
  expect_true(all(range_in_states_pairwise(states$state, states$parent_state)))
  expect_true("B" %in% colnames(grafted$simmap$mapped.edge))

  # With the exact matching, the same allowed state stops with the list of the states available
  expect_error(graft_missing_taxa(toy$phy, mtl, simmap = range_simmap, allowed_states = "B", method = "uniform", verbose = FALSE), "not in the simmap: 'B'")
  # An area that is not in the simmap, or a range contained in no state, still stops
  expect_error(graft_missing_taxa(toy$phy, mtl, simmap = range_simmap, allowed_states = "D", treat_states_as_ranges = TRUE, method = "uniform", verbose = FALSE), "Areas available in the simmap")
  expect_error(graft_missing_taxa(toy$phy, mtl, simmap = range_simmap, allowed_states = "AC", treat_states_as_ranges = TRUE, method = "uniform", verbose = FALSE), "No state of the simmap contains")
})

## Test that a range inherited that is not a state of the simmap is in the returned simmap
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE adds the ranges inherited that are not states of the simmap to the returned simmap",
{
  toy <- toy_setup()
  # States AB and BC only (the toy state A is also AB)
  two_states <- relabel_simmap(make_toy_simmap(toy$phy), c(A = "AB", B = "BC", C = "AB"))
  two_states$mapped.edge <- NULL
  mtl <- toy_missing_taxa(toy)
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = two_states, allowed_states = "B", treat_states_as_ranges = TRUE, method = "uniform", seed = 3)
  expect_identical(graft_states_in_simmap(grafted, mtl$Taxon_name)$state, rep("B", 5))
  expect_true("B" %in% colnames(grafted$simmap$mapped.edge))
  expect_setequal(colnames(grafted$simmap$mapped.edge), c("AB", "BC", "B"))
  expect_equal(unname(rowSums(grafted$simmap$mapped.edge)), grafted$edge.length, tolerance = 1e-6)

  # A taxon without allowed states can be grafted on the branches of the taxa that inherited "B"
  chain_names <- paste0("X", 1:12)
  mtl_chain <- list(Taxon_name = chain_names, Target_clade = stats::setNames(rep(list(toy$B), 12), chain_names),
                    allowed_states = stats::setNames(c(as.list(rep("B", 6)), vector("list", 6)), chain_names))
  chain <- quiet_graft(graft_missing_taxa, toy$phy, mtl_chain, simmap = two_states, treat_states_as_ranges = TRUE, method = "uniform", seed = 5)
  expect_equal(unname(rowSums(chain$simmap$mapped.edge)), chain$edge.length, tolerance = 1e-6)
  expect_true(all(unlist(lapply(chain$simmap$maps, names)) %in% c("AB", "BC", "B")))
  chain_states <- graft_states_in_simmap(chain, chain_names)
  expect_true(any(chain_states$parent_state[7:12] == "B")) # taxa without allowed states grafted on the branch of a taxon that inherited "B"
  expect_identical(chain_states$state[chain_states$parent_state == "B"], rep("B", sum(chain_states$parent_state == "B")))
})

## Test that the states and the allowed ranges are strictly checked with TRUE only
test_that("graft_missing_taxa() with treat_states_as_ranges = TRUE stops if the states or the allowed ranges are not coded as in BioGeoBEARS",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  range_simmap <- make_toy_range_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  graft <- function (...) graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, ...)

  # Incorrect states in the simmap: named in the message, and accepted with the exact matching
  bad_simmap <- relabel_simmap(simmap, c(A = "A", B = "red", C = "CB"))
  expect_error(graft(simmap = bad_simmap, allowed_states = "A", treat_states_as_ranges = TRUE),
               "the states in the simmap must be biogeographic ranges.*'CB' \\(areas not in alphabetical order \\(expected 'BC'\\)\\).*'red' \\(lowercase letters")
  expect_no_error(graft(simmap = bad_simmap, allowed_states = c("A", "red", "CB"), treat_states_as_ranges = FALSE))

  # Incorrect allowed ranges, global and taxon-level
  expect_error(graft(simmap = range_simmap, allowed_states = "ba", treat_states_as_ranges = TRUE), "the global 'allowed_states' argument.*'ba' \\(lowercase letters")
  mtl_bad <- mtl
  mtl_bad$allowed_states <- list(X1 = "AB", X2 = "BA", X3 = NULL, X4 = NULL, X5 = NULL)
  expect_error(graft_missing_taxa(toy$phy, mtl_bad, simmap = range_simmap, treat_states_as_ranges = TRUE, method = "uniform", verbose = FALSE),
               "taxon 'X2'.*'BA' \\(areas not in alphabetical order")

  # An allowed area absent from the simmap: the areas available are listed
  expect_error(graft(simmap = range_simmap, allowed_states = c("A", "D"), treat_states_as_ranges = TRUE), "not in the simmap: 'D'.*Areas available in the simmap: 'A', 'B', 'C'")
  # A range contained in no state of the simmap
  expect_error(graft(simmap = range_simmap, allowed_states = "AC", treat_states_as_ranges = TRUE), "No state of the simmap contains the allowed range")
})


# ---- Species attributes: simmap without effect ----

## Test that a simmap that is neither restricting grafting nor returned sends a warning
test_that("graft_missing_taxa() warns if the simmap has no effect: no allowed states and return_simmap = FALSE",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  graft <- function (...) graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, ...)
  has_warning <- function (expr) any(grepl("no effect", testthat::capture_warnings(expr)))

  expect_warning(graft(simmap = simmap, return_simmap = FALSE), "without any allowed states.*'return_simmap = FALSE'.*no effect")
  expect_warning(graft(simmap = simmap, return_simmap = FALSE, treat_states_as_ranges = TRUE), "no effect")
  # The simmap is returned (tracking only), or restricts grafting (global or taxon-level allowed states)
  expect_false(has_warning(graft(simmap = simmap, return_simmap = TRUE)))
  expect_false(has_warning(graft(simmap = simmap, allowed_states = c("A", "B"), return_simmap = FALSE)))
  mtl_states <- mtl
  mtl_states$allowed_states <- list(X1 = c("A", "B"), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL)
  expect_false(has_warning(graft_missing_taxa(toy$phy, mtl_states, simmap = simmap, return_simmap = FALSE, method = "uniform", verbose = FALSE)))
  # Without simmap, return_simmap is silently ignored
  expect_false(has_warning(graft(return_simmap = FALSE)))
})


# ---- Species attributes with method = "BD_clades" ----

## Test that all the states allowed gives exactly the same phylogeny as without simmap, with the BD method
test_that("graft_missing_taxa() with method = 'BD_clades' and all the states allowed gives the same phylogeny as without simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)
  # Taxa without access to the stem, and several groups in the same clade
  mtl$include_target_clade_stem <- c(X1 = TRUE, X2 = FALSE, X3 = NA, X4 = NA, X5 = TRUE)

  for (seed in 1:3)
  {
    without <- graft_missing_taxa(toy$phy, mtl, method = "BD_clades", verbose = FALSE, seed = seed)

    # Global states
    with_all <- graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B", "C"), method = "BD_clades", verbose = FALSE, seed = seed)
    expect_identical(drop_simmap(with_all), without)

    # Taxon-level states, in a different order, with the global states for some taxa: the taxa are in the same groups
    mtl_taxa <- mtl
    mtl_taxa$allowed_states <- list(X1 = c("C", "A", "B"), X2 = NULL, X3 = c("B", "A", "C"), X4 = NA, X5 = NULL)
    with_taxa <- graft_missing_taxa(toy$phy, mtl_taxa, simmap = simmap, allowed_states = c("A", "B", "C"), method = "BD_clades", verbose = FALSE, seed = seed)
    expect_identical(drop_simmap(with_taxa), without)
  }
})

## Test that taxa are only grafted on pieces of branches in their allowed states with the BD method, and that the constraints are still respected
test_that("graft_missing_taxa() with method = 'BD_clades' grafts the taxa in their allowed states only, global or taxon-level, and respects the other constraints",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  taxa <- c(paste0("X", 1:8), "Y1", "Y2", "Y3")
  mtl <- list(Taxon_name = taxa,
              Target_clade = c(rep(list(toy$A), 8), rep(list(toy$C), 3)),
              Excluded_subclades = c(rep(list(list(toy$B)), 8), rep(list(NULL), 3)),
              allowed_states = c(rep(list(NULL), 8), rep(list("C"), 3)))
  names(mtl$Target_clade) <- names(mtl$Excluded_subclades) <- names(mtl$allowed_states) <- taxa

  n_checked <- 0
  for (seed in 1:4)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "BD_clades", seed = seed)
    expect_valid_graft(grafted, toy$phy, mtl)
    for (taxon in taxa)
    {
      state <- state_at_graft(grafted, simmap, taxon, taxa)
      if (is.na(state)) { next }
      n_checked <- n_checked + 1
      expect_true(state %in% (if (taxon %in% c("Y1", "Y2", "Y3")) "C" else c("A", "B")), info = paste(taxon, "seed", seed))
    }
  }
  # The states of many taxa could be checked independently (the others were grafted on the branches of other grafted taxa)
  expect_gt(n_checked, 10)

  # Restricting the states changes the result, and the same seed gives the same result
  unrestricted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B", "C"), method = "BD_clades", seed = 1)
  expect_false(identical(unrestricted$edge.length, grafted$edge.length))
  again <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "BD_clades", seed = 4)
  expect_identical(again, grafted)
})

## Test that grafted taxa have the state of the point where they were grafted, so taxa with other allowed states cannot be grafted on their branches (BD method)
test_that("graft_missing_taxa() with method = 'BD_clades' does not graft taxa on the branches of grafted taxa in other states",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # The terminal branch of the tip 's' is entirely in state "A": the taxa grafted there are in state "A"
  expect_identical(names(simmap$maps[[match(match(toy$s, toy$phy$tip.label), toy$phy$edge[, 2])]]), "A")
  tip_taxa <- paste0("S", 1:5)
  clade_taxa <- paste0("Y", 1:10)
  taxa <- c(tip_taxa, clade_taxa)
  mtl <- list(Taxon_name = taxa,
              Target_clade = c(rep(list(toy$s), 5), rep(list(toy$A), 10)),
              allowed_states = c(rep(list("A"), 5), rep(list("B"), 10)))
  names(mtl$Target_clade) <- names(mtl$allowed_states) <- taxa

  for (seed in 1:5)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, method = "BD_clades", seed = seed)
    s_subtree <- phytools::getDescendants(grafted, ape::getMRCA(grafted, c(toy$s, tip_taxa)))
    s_subtree_tips <- grafted$tip.label[s_subtree[s_subtree <= ape::Ntip(grafted)]]
    expect_true(all(tip_taxa %in% s_subtree_tips))
    # The Y taxa (state "B" only) are never grafted in this subtree, which is entirely in state "A"
    expect_false(any(clade_taxa %in% s_subtree_tips), info = paste("seed", seed))
    # When they are grafted on a branch of the backbone, it is in state "B"
    states <- vapply(clade_taxa, function (taxon) state_at_graft(grafted, simmap, taxon, taxa), character(1))
    expect_true(all(states[!is.na(states)] == "B"), info = paste("seed", seed))
  }
})

## Test the groups of taxa with the BD method: taxa of the same clade with different allowed states are drawn in separate groups, with one fit per clade
test_that("with method = 'BD_clades' and a simmap, taxa of a clade are grouped by allowed states, with one BD fit per clade",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # X1-X4 target A: two sets of allowed states, in alternating order (and the same set in a different order); X5 targets B, nested in A
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A, X4 = toy$A, X5 = toy$B),
              allowed_states = list(X1 = c("A", "B"), X2 = "C", X3 = c("B", "A"), X4 = "C", X5 = c("A", "B", "C")))

  n_fits <- 0
  grafts <- list()
  original_fit <- fit_bd_on_clade
  original_graft <- graft_from_divergence_age
  local_mocked_bindings(
    fit_bd_on_clade = function (...) { n_fits <<- n_fits + 1; original_fit(...) },
    graft_from_divergence_age = function (...) {
      args <- list(...)
      grafts[[length(grafts) + 1]] <<- data.frame(taxon = args$tip_label, states = paste(sort(args$allowed_states), collapse = "|"))
      original_graft(...) })

  quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "BD_clades", simmap = simmap, seed = 1)
  grafts <- do.call(rbind, grafts)

  # One fit for B, one for A (although the taxa of A form two groups)
  expect_equal(n_fits, 2)
  # B first (nested), then the taxa of A grouped by states: each set of states is in one block
  expect_equal(grafts$taxon[1], "X5")
  states_A <- grafts$states[grafts$taxon != "X5"]
  expect_equal(length(rle(states_A)$lengths), 2)
  expect_setequal(grafts$taxon[grafts$states == "A|B" & grafts$taxon != "X5"], c("X1", "X3"))
  expect_setequal(grafts$taxon[grafts$states == "C"], c("X2", "X4"))
})

## Test the grafting in a narrow and young window of allowed states: the ages drawn that no branch crosses are redrawn, and the taxa cluster in the window
test_that("graft_missing_taxa() with method = 'BD_clades' redraws the ages that no branch in an allowed state crosses, and grafts all the taxa in a narrow young window",
{
  toy <- toy_setup()
  window_age <- 4 # The state "B" is only found below this age (the root age is 50)
  simmap <- make_young_state_simmap(toy$phy, window_age)
  taxa <- paste0("Y", 1:8)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 8), taxa))

  n_unusable <- 0
  original_graft <- graft_from_divergence_age
  local_mocked_bindings(graft_from_divergence_age = function (...) {
    out <- original_graft(...)
    if (is.null(out)) { n_unusable <<- n_unusable + 1 }
    out })

  for (seed in 1:4)
  {
    grafted <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", verbose = FALSE, seed = seed))
    expect_valid_graft(grafted, toy$phy, mtl)
    # All the taxa are in state "B": the nodes created by the grafts are all younger than the window, the first ones on the backbone and the others on the lineages grafted before
    node_ages <- vapply(taxa, function (taxon) graft_age(grafted, taxon), numeric(1))
    expect_true(all(node_ages <= window_age + 1e-6), info = paste("seed", seed))
  }
  # Many ages drawn could not be used, and were redrawn (the clade A has a crown age of about 13 and the window is below 4)
  expect_gt(n_unusable, 0)

  # Without redraw, with the same simmap but all the states allowed, nothing is redrawn
  n_unusable <- 0
  suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "BD_clades", verbose = FALSE, seed = 1))
  expect_equal(n_unusable, 0)
})

## Test the stop when no usable age is found
test_that("graft_missing_taxa() with method = 'BD_clades' stops after 'max_age_redraws' redraws of the ages, with an informative message",
{
  toy <- toy_setup()
  # The state "B" only exists in a tiny window: no age drawn can be used, except by (very unlikely) chance
  simmap <- make_young_state_simmap(toy$phy, 1e-4)
  taxa <- paste0("Y", 1:3)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 3), taxa))
  expect_error(suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 3, verbose = FALSE, seed = 1)),
               "No usable divergence age.*after 3 redraws.*'max_age_redraws'.*Y1.*allowed states: B")

  # The message gives the time window of the divergence ages (here, up to the stem age of A) and the share of it crossed by available branches in allowed states
  stem_age_A <- toy$root_age - phytools::nodeheight(toy$phy, toy$phy$edge[match(toy$node_A, toy$phy$edge[, 2]), 1])
  expect_error(suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 3, verbose = FALSE, seed = 1)),
               paste0("Time window of the divergence ages: from age 0 to age ", signif(stem_age_A, 4), " \\(stem age of the target clade\\)\\.\\n",
                      "Share of this window crossed by the available branches in the allowed states \\(current tree\\):\\n  - Y1: 0.00% \\(allowed states: B\\)"))
})

## Test that the ages with no usable branch are also redrawn without simmap, e.g., ages where all the branches are excluded
test_that("graft_missing_taxa() with method = 'BD_clades' redraws the ages with no available branch without simmap, and stops with the time window and its share covered",
{
  toy <- toy_setup()
  # The two child clades of A are excluded with their stem: only the stem region of A remains, so the ages below the crown age of A cannot be used
  children <- toy$phy$edge[toy$phy$edge[, 1] == toy$node_A, 2]
  child_tips <- lapply(children, function (node) if (node <= ape::Ntip(toy$phy)) toy$phy$tip.label[node] else pair_for_node(toy$phy, node))
  crown_age_A <- toy$root_age - phytools::nodeheight(toy$phy, toy$node_A)
  stem_age_A <- toy$root_age - phytools::nodeheight(toy$phy, toy$phy$edge[match(toy$node_A, toy$phy$edge[, 2]), 1])
  taxa <- paste0("Y", 1:3)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 3), taxa),
              Excluded_subclades = setNames(rep(list(child_tips), 3), taxa))

  n_unusable <- 0
  original_graft <- graft_from_divergence_age
  local_mocked_bindings(graft_from_divergence_age = function (...) {
    out <- original_graft(...)
    if (is.null(out)) { n_unusable <<- n_unusable + 1 }
    out })

  for (seed in 1:5)
  {
    grafted <- graft_missing_taxa(toy$phy, mtl, method = "BD_clades", verbose = FALSE, seed = seed)
    expect_valid_graft(grafted, toy$phy, mtl)
  }
  expect_gt(n_unusable, 0)

  # Without any usable age: stop after 'max_age_redraws', with the window and the share of it crossed by the stem of A, now in the tree
  # (mocked draws of TreeSim::corsim() with ages that are too young)
  original_corsim <- TreeSim::corsim
  local_mocked_bindings(corsim = function (x, ..., missing) c(x, rep(0.001, missing)), .package = "TreeSim")
  expect_error(graft_missing_taxa(toy$phy, mtl, method = "BD_clades", max_age_redraws = 2, verbose = FALSE, seed = 1),
               paste0("No usable divergence age.*after 2 redraws.*Time window of the divergence ages: from age 0 to age ", signif(stem_age_A, 4), " \\(stem age of the target clade\\)\\.\\n",
                      "Share of this window crossed by the available branches \\(current tree\\):\\n  - Y1: ", sprintf("%.2f%%", 100 * (stem_age_A - crown_age_A) / stem_age_A),
                      ".*adjust the excluded subclades or the stem availability"))
})

## Test the warning on the share of the time window in allowed states, with the BD method
test_that("graft_missing_taxa() with method = 'BD_clades' warns about the taxa with a small share of the time window of their divergence ages in allowed states",
{
  toy <- toy_setup()
  simmap <- make_young_state_simmap(toy$phy, 4)
  taxa <- paste0("Y", 1:3)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 3), taxa))

  # The window of the state "B" covers 4 / 50 = 8% of the time window of the taxa (from the present to the stem age of A, which is the root)
  expect_warning(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", verbose = FALSE),
                 "3 taxa have less than 10% of the time window of their divergence ages covered by available branches in their allowed states.*Y1 .*8.00% in allowed states")
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", min_valid_fraction = 0.05, method = "BD_clades", verbose = FALSE))
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "BD_clades", verbose = FALSE))
  # Without the stem, the time window is shorter (up to the crown age of A, 12.8): a larger share is covered
  mtl_crown <- c(mtl, list(include_target_clade_stem = setNames(rep(FALSE, 3), taxa)))
  expect_no_warning(graft_missing_taxa(toy$phy, mtl_crown, simmap = simmap, allowed_states = "B", min_valid_fraction = 0.2, method = "BD_clades", verbose = FALSE))
  # The warning with the uniform method is based on the length of the branches
  expect_warning(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", min_valid_fraction = 0.5, method = "uniform", verbose = FALSE),
                 "of the length of their available branches in their allowed states")
})

## Test the node labels of the grafted phylogeny with a simmap and the BD method
test_that("graft_missing_taxa() with a simmap and method = 'BD_clades' gives the same node labels as without",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)

  grafted <- graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "BD_clades", verbose = FALSE, seed = 1)
  expect_null(grafted$node.label)

  phy_support <- toy$phy
  phy_support$node.label <- paste0("support", seq_len(phy_support$Nnode))
  grafted <- graft_missing_taxa(phy_support, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "BD_clades", verbose = FALSE, seed = 1)
  expect_length(grafted$node.label, grafted$Nnode)
  expect_equal(sum(grafted$node.label == ""), length(mtl$Taxon_name))
  expect_setequal(grafted$node.label[grafted$node.label != ""], phy_support$node.label)
})

## Test that a simmap without allowed states does not restrict grafting (same phylogeny as without simmap) but the states are followed and returned
test_that("graft_missing_taxa() with a simmap and no allowed states does not restrict the grafting, and returns the updated simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)

  for (method in c("uniform", "BD_clades"))
  {
    for (seed in 1:2)
    {
      without <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = method, verbose = FALSE, seed = seed))
      with_simmap <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, method = method, verbose = FALSE, seed = seed))
      expect_identical(drop_simmap(with_simmap), without)
      expect_s3_class(with_simmap$simmap, "simmap")
    }
  }

  # Taxa without a taxon-level value use all the states when no global states are given, the others are restricted
  mtl$allowed_states <- list(X1 = "A", X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL)
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, method = "uniform", seed = 1)
  expect_s3_class(grafted$simmap, "simmap")
})

## Test the updated simmap returned with the grafted phylogeny
test_that("the updated simmap matches the grafted phylogeny, keeps the states of the backbone, and gives the state of the graft point to the grafted taxa",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  taxa <- c(paste0("X", 1:8), "Y1", "Y2", "Y3")
  mtl <- list(Taxon_name = taxa,
              Target_clade = c(rep(list(toy$A), 8), rep(list(toy$C), 3)),
              Excluded_subclades = c(rep(list(list(toy$B)), 8), rep(list(NULL), 3)),
              allowed_states = c(rep(list(NULL), 8), rep(list("C"), 3)))
  names(mtl$Target_clade) <- names(mtl$Excluded_subclades) <- names(mtl$allowed_states) <- taxa

  n_checked <- 0
  for (method in c("uniform", "BD_clades"))
  {
    for (seed in 1:3)
    {
      grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = method, seed = seed)
      sm <- grafted$simmap
      info <- paste(method, "seed", seed)

      # The updated simmap is a simmap of the grafted phylogeny: same tips, branches and branch lengths, in the same order
      expect_s3_class(sm, "simmap")
      expect_s3_class(sm, "phylo")
      expect_identical(sm$edge, grafted$edge)
      expect_identical(sm$tip.label, grafted$tip.label)
      expect_identical(sm$edge.length, grafted$edge.length)
      expect_no_error(validate_simmap(sm, grafted))
      expect_equal(vapply(sm$maps, sum, numeric(1)), grafted$edge.length, info = info)
      expect_equal(unname(rowSums(sm$mapped.edge)), grafted$edge.length, info = info)
      expect_identical(colnames(sm$mapped.edge), c("A", "B", "C"))
      expect_identical(rownames(sm$mapped.edge), paste(grafted$edge[, 1], grafted$edge[, 2], sep = ","))

      # Every grafted branch is in a single state, and the terminal branch of each taxon is in the state of its graft point
      expect_true(all(lengths(sm$maps[grafted$grafted_edges]) == 1), info = info)
      for (taxon in taxa)
      {
        tip_edge <- match(match(taxon, grafted$tip.label), grafted$edge[, 2])
        expect_length(sm$maps[[tip_edge]], 1)
        expect_true(names(sm$maps[[tip_edge]]) %in% (if (taxon %in% c("Y1", "Y2", "Y3")) "C" else c("A", "B")), info = paste(taxon, info))
        state <- state_at_graft(grafted, simmap, taxon, taxa)
        if (is.na(state)) { next }
        n_checked <- n_checked + 1
        expect_identical(names(sm$maps[[tip_edge]]), state)
      }

      # Pruning the grafted taxa gives back the time of the backbone in each state
      pruned <- phytools::drop.tip.simmap(sm, taxa)
      expect_equal(colSums(pruned$mapped.edge)[colnames(simmap$mapped.edge)], colSums(simmap$mapped.edge))
      expect_false("simmap" %in% names(simmap))
    }
  }
  expect_gt(n_checked, 20)
})

## Test the argument return_simmap
test_that("graft_missing_taxa() returns the updated simmap only with a simmap and return_simmap = TRUE, which does not change the phylogeny", 
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- toy_missing_taxa(toy)

  with_simmap <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = 3)
  expect_s3_class(with_simmap$simmap, "simmap")
  # The default is TRUE
  expect_identical(quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = 3, return_simmap = TRUE), with_simmap)
  # FALSE removes the item only
  no_simmap_item <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = 3, return_simmap = FALSE)
  expect_null(no_simmap_item$simmap)
  expect_identical(no_simmap_item, drop_simmap(with_simmap))
  # Without a simmap, return_simmap is silently ignored
  without <- graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, seed = 3)
  expect_no_warning(ignored <- graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, seed = 3, return_simmap = TRUE))
  expect_identical(ignored, without)
  expect_null(without$simmap)
  # Invalid values
  expect_error(graft_missing_taxa(toy$phy, mtl, simmap = simmap, return_simmap = "yes", verbose = FALSE), "'return_simmap' must be TRUE or FALSE")
  expect_error(graft_missing_taxa(toy$phy, mtl, return_simmap = NA, verbose = FALSE), "'return_simmap' must be TRUE or FALSE")
  expect_error(graft_missing_taxa(toy$phy, mtl, return_simmap = c(TRUE, FALSE), verbose = FALSE), "'return_simmap' must be TRUE or FALSE")
})

## Test the updated simmap when the backbone has node labels (they must be restored in the simmap too)
test_that("the updated simmap has the same node labels as the grafted phylogeny",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  phy_support <- toy$phy
  phy_support$node.label <- paste0("support", seq_len(phy_support$Nnode))
  simmap_support <- simmap
  simmap_support$node.label <- phy_support$node.label
  mtl <- toy_missing_taxa(toy)
  grafted <- quiet_graft(graft_missing_taxa, phy_support, mtl, simmap = simmap_support, allowed_states = c("A", "B"), method = "uniform", seed = 1)
  expect_identical(grafted$simmap$node.label, grafted$node.label)
  expect_false(any(grepl("^(bb|gr|ms)[0-9]", grafted$simmap$node.label)))
})
