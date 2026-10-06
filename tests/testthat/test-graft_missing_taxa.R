
##### Test suite for graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Build toy 'missing_taxa_list' to test graft_missing_clades(): 2 taxa in A (one excluding B), 1 in B, 1 in C, 1 on a single tip
toy_missing_taxa <- function (toy)
{
  name_items(list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(toy$A, toy$A, toy$B, toy$C, toy$s),
       Excluded_subclades = list(list(toy$B), NULL, NULL, NULL, NULL)))
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
  expect_error(graft_missing_taxa(make_toy_phy(n_tips = 4), name_items(list(Taxon_name = "X1", Target_clade = list(c("t1", "t2")))), verbose = FALSE),
               "too small to fit a birth-death model")

  # Lists
  expect_error(graft_missing_taxa(phy, name_items(list(Taxon_name = character(0), Target_clade = list())), verbose = FALSE), "empty")
  expect_error(graft_missing_taxa(phy, name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A))), verbose = FALSE), "same length as 'Taxon_name'")
  expect_error(graft_missing_taxa(phy, name_items(list(Taxon_name = c("X1", "X1"), Target_clade = list(toy$A, toy$A))), verbose = FALSE), "unique names")
  expect_error(graft_missing_taxa(phy, name_items(list(Taxon_name = toy$s, Target_clade = list(toy$A))), verbose = FALSE), "absent from the backbone")
  expect_error(graft_missing_taxa(phy, name_items(list(Taxon_name = "X1", Target_clade = list(c("zz", "yy")))), verbose = FALSE), "missing from the backbone")
  expect_error(graft_missing_taxa(phy, name_items(list(Taxon_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(c("zz", "yy")))), verbose = FALSE),
               "missing from the backbone")

  # Temporary file path
  expect_error(graft_missing_taxa(phy, mtl, temp_file_path = "temp.txt", verbose = FALSE), "must end with '.rds'")
  expect_error(graft_missing_taxa(phy, mtl, temp_file_path = file.path(tempdir(), "no_dir", "t.rds"), verbose = FALSE), "does not exist")
})

## Test that exclusive subclades are valid relative to their target clade
test_that("graft_missing_taxa() checks the position of the excluded subclades relative to the target clade",
{
  toy <- toy_setup()
  one <- function (target, excluded) name_items(list(Taxon_name = "X1", Target_clade = list(target), Excluded_subclades = list(list(excluded))))

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
  mtl <- name_items(list(Taxon_name = "X1", Target_clade = list(whole_tree)))

  expect_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", include_target_clade_stem = TRUE, verbose = FALSE), "whole tree")
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", include_target_clade_stem = FALSE, verbose = FALSE))
})


# ---- Valid grafting ----

## Test that excluded subclades are respected when grafting
test_that("excluded subclades are respected by both methods",
{
  toy <- toy_setup()
  mtl <- name_items(list(Taxon_name = paste0("X", 1:4), Target_clade = rep(list(toy$A), 4),
              Excluded_subclades = rep(list(list(toy$B)), 4)))

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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A, toy$A)))

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
  mtl <- name_items(list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(toy$s, toy$s, toy$A, toy$A),
              include_target_clade_stem = c(FALSE, TRUE, FALSE, NA), Excluded_subclades = rep(list(NULL), 4)))
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
  mtl$include_target_clade_stem <- c(NA, TRUE, FALSE, NA)
  mtl <- name_items(mtl)
  expect_false(any(grepl("single-tip target clade", warning_messages(
    graft_missing_taxa(toy$phy, mtl, method = "uniform", include_target_clade_stem = FALSE, verbose = FALSE)))))

  # The taxon is grafted on the branch of its tip, whatever the value
  mtl$include_target_clade_stem <- c(FALSE, TRUE, NA, NA)
  mtl <- name_items(mtl)
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1)
  checks <- check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)
  expect_equal(checks$target_position[1:2], c("tip_branch", "tip_branch"))
})

## Test that the taxon-level value overrides the global argument, with the 'uniform' method
test_that("with method = 'uniform', include_target_clade_stem overrides the global include_target_clade_stem per taxon, and NA uses it",
{
  toy <- toy_setup()
  mtl <- name_items(list(Taxon_name = paste0("X", 1:4), Target_clade = rep(list(toy$A), 4),
              include_target_clade_stem = c(TRUE, FALSE, NA, FALSE)))

  stem_used <- NULL
  original <- draw_random_branch_weighted_per_length
  local_mocked_bindings(draw_random_branch_weighted_per_length = function (...) {
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
  mtl <- name_items(list(Taxon_name = paste0("X", 1:5), Target_clade = list(toy$A, toy$A, toy$A, toy$A, toy$B),
              include_target_clade_stem = c(TRUE, FALSE, NA, FALSE, TRUE), Excluded_subclades = rep(list(NULL), 5)))

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
  mtl <- name_items(list(Taxon_name = paste0("X", 1:3), Target_clade = rep(list(toy$A), 3), include_target_clade_stem = c(TRUE, FALSE, NA)))
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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(whole_tree, whole_tree), include_target_clade_stem = c(FALSE, TRUE)))

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
  mtl <- name_items(list(Taxon_name = paste0("X", 1:4), Target_clade = rep(list(toy$s), 4)))
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
  mtl <- name_items(list(Taxon_name = paste0("X", 1:4), Target_clade = rep(list(toy$A), 4), Excluded_subclades = rep(list(NULL), 4)))

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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A, toy$A), Excluded_subclades = list(list(toy$B, toy$C), NULL)))

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
  original <- draw_random_branch_weighted_per_length
  local_mocked_bindings(draw_random_branch_weighted_per_length = function (...) {
    excluded_used[[length(excluded_used) + 1]] <<- list(...)$excluded_nodes; original(...) })

  mtl <- name_items(list(Taxon_name = "X1", Target_clade = list(toy$A), Excluded_subclades = list(list(toy$B, toy$C))))

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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A, toy$A), Excluded_subclades = list(list(toy$B, toy$C), list(toy$B)),
              exclude_subclades_stem = list(c(TRUE, FALSE), NA)))
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
  mtl <- name_items(list(Taxon_name = "X1", Target_clade = list(cherry), Excluded_subclades = list(list(cherry[1]))))

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
  same_tip <- name_items(list(Taxon_name = "X1", Target_clade = list(toy$s), Excluded_subclades = list(list(toy$s))))
  expect_error(graft_missing_taxa(toy$phy, same_tip, method = "uniform", verbose = FALSE), "graft_missing_clades\\(\\)")
  # ... but not if the exclusion has no effect
  expect_true("X1" %in% quiet_graft(graft_missing_taxa, toy$phy, same_tip, method = "uniform", exclude_subclades_stem = FALSE, seed = 1)$tip.label)
})

## Test the warning for single-tip excluded subclades that exclude nothing
test_that("graft_missing_taxa() warns once, listing the single-tip excluded subclades with exclude_subclades_stem = FALSE",
{
  toy <- toy_setup()
  mtl <- name_items(list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(toy$A, toy$A, toy$A),
              Excluded_subclades = list(list(toy$s, toy$B), list(toy$B), NULL),
              exclude_subclades_stem = list(c(FALSE, TRUE), FALSE, NULL)))

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
  mtl <- name_items(list(Taxon_name = paste0("X", 1:7),
              Target_clade = list(toy$A, toy$B, toy$s, toy$A, toy$B, toy$C, toy$s)))

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
  mtl <- name_items(list(Taxon_name = c("b1", "b2", "c1", "a1", "a2", "a3"),
              Target_clade = list(toy$B, toy$B, toy$C, toy$A, toy$A, toy$A)))

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
  mtl <- name_items(list(Taxon_name = c("b1", "b2", "b3", "c1", "c2", "a1", "a2", "a3", "a4"),
              Target_clade = c(rep(list(toy$B), 3), rep(list(toy$C), 2), rep(list(toy$A), 4))))
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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$B, toy$A)))
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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$B, toy$A)))

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
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A, toy$B), include_target_clade_stem = c(TRUE, NA),
                         Excluded_subclades = list(list(toy$B), NULL), exclude_subclades_stem = list(FALSE, NULL)))
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
  reordered <- name_items(list(Taxon_name = c("X2", "X1"), Target_clade = list(toy$B, toy$A), Excluded_subclades = list(NULL, list(toy$B))))
  expect_valid_graft(quiet_graft(graft_missing_taxa, toy$phy, reordered, method = "uniform", seed = 1), toy$phy, reordered)
})

## Test that the optional names of the excluded subclades are used in the warnings
test_that("graft_missing_taxa() uses the names of the excluded subclades in its warnings",
{
  toy <- toy_setup()
  mtl <- name_items(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A, toy$A),
                         Excluded_subclades = list(list(tip_s = toy$s, clade_B = toy$B), list(toy$D)),
                         exclude_subclades_stem = list(c(FALSE, TRUE), NULL)))

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

## Test that the table of sampling fractions reports the names of the missing taxa of each clade
test_that("sampling_fraction_df lists the missing taxa counted in each clade",
{
  toy <- toy_setup()
  mtl <- name_items(list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(toy$A, toy$A, toy$B)))
  totals <- list(Target_clade = list(toy$A, toy$B), total_richness = c(NA, NA))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", total_richness_per_clades = totals, seed = 1)
  df <- grafted$sampling_fraction_df
  expect_equal(df$Missing_taxa, c("X1, X2, X3", "X3")) # B is nested in A
  expect_equal(df$nb_missing_species, c(3, 1))
})
