
##### Test suite for check_grafting_constraints() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Test correct validation with a correctly grafted phylogeny
test_that("check_grafting_constraints() validates a correct grafting",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxon_name = c("X1", "X2", "X3"),
                            Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$s),
                            Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL))
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
                                       list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = list(toy$B))),
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
                                       list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = NULL)),
                                       verbose = FALSE)
  expect_false(checks$target_ok)
  expect_false(attr(checks, "all_ok"))
})

## Test flagging of a shift in backbone branch length, but not topology
test_that("check_grafting_constraints() reports a shift of the backbone branch lengths",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = NULL))
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
  missing_taxa_list <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = NULL))
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
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A), Excluded_subclades = list(X1 = NULL, X2 = NULL))
  check <- function (stem_item, global = TRUE)
  {
    mtl$include_target_clade_stem <- stem_item
    check_grafting_constraints(on_stem, toy$phy, mtl, include_target_clade_stem = global, verbose = FALSE)
  }

  # Stem accepted for X1 only
  checks <- check(c(X1 = TRUE, X2 = FALSE))
  expect_equal(checks$target_position, c("stem", "stem"))
  expect_equal(checks$include_target_clade_stem, c(TRUE, FALSE))
  expect_equal(checks$target_ok, c(TRUE, FALSE))
  expect_false(attr(checks, "all_ok"))

  # NA uses the global value
  expect_equal(check(c(X1 = NA, X2 = FALSE), global = TRUE)$target_ok, c(TRUE, FALSE))
  expect_equal(check(c(X1 = NA, X2 = TRUE), global = FALSE)$target_ok, c(FALSE, TRUE))
  # No item: the global value applies to all taxa
  expect_true(all(check(NULL, global = TRUE)$target_ok))
  expect_false(any(check(NULL, global = FALSE)$target_ok))

  # Invalid length
  expect_error(check(TRUE), "same length as 'Taxon_name'")
})

## Test that the function detect issues with taxa grafted on the stem of its excluded subclades
test_that("check_grafting_constraints() flags a taxon sister to an excluded subclade only if its stem is excluded",
{
  toy <- toy_setup()
  # X1 is grafted on the stem of B, which is excluded; X2 is grafted on the terminal branch of the excluded tip 's'
  on_stem_B <- graft_tip_on_branch(toy$phy, "X1", toy$node_B, 0.5)
  on_tip_s <- graft_tip_on_branch(toy$phy, "X2", match(toy$s, toy$phy$tip.label), 0.5)
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = list(toy$B)))
  mtl_tip <- list(Taxon_name = "X2", Target_clade = list(X2 = toy$A), Excluded_subclades = list(X2 = list(toy$s)))

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
  expect_true(check(on_stem_B, c(mtl, list(exclude_subclades_stem = list(X1 = FALSE))))$excluded_ok)
  expect_false(check(on_stem_B, c(mtl, list(exclude_subclades_stem = list(X1 = NA))), exclude_subclades_stem = TRUE)$excluded_ok)
  expect_false(check(on_tip_s, c(mtl_tip, list(exclude_subclades_stem = list(X2 = TRUE))), exclude_subclades_stem = FALSE)$excluded_ok)

  # A taxon within the excluded subclade is always a violation
  inside_B <- graft_tip_on_branch(toy$phy, "X1", match(toy$tipsB[1], toy$phy$tip.label), 0.5)
  expect_false(check(inside_B, mtl, exclude_subclades_stem = FALSE)$excluded_ok)

  # Several excluded subclades with their own setting: only the violated one is reported
  mtl_two <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = list(toy$C, toy$B)),
                  exclude_subclades_stem = list(X1 = c(TRUE, FALSE)))
  expect_true(check(on_stem_B, mtl_two)$excluded_ok) # B's stem is not excluded
  mtl_two$exclude_subclades_stem <- list(X1 = c(FALSE, TRUE))
  checks <- check(on_stem_B, mtl_two)
  expect_false(checks$excluded_ok)
  expect_equal(checks$excluded_violated, "2")

  # A taxon sister to a clade, but with another backbone tip in between, is not on its stem
  elsewhere <- graft_tip_on_branch(toy$phy, "X1", match(toy$tipsC[1], toy$phy$tip.label), 0.5)
  expect_true(check(elsewhere, mtl)$excluded_ok)

  # Invalid values
  expect_error(check(on_stem_B, c(mtl, list(exclude_subclades_stem = list(X1 = c(TRUE, FALSE))))), "excluded subclades")
  expect_error(check(on_stem_B, mtl, exclude_subclades_stem = NA), "single logical value")
})

## Test invalid input stops
test_that("check_grafting_constraints() stops on invalid inputs, and can stop on failure",
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A), Excluded_subclades = list(X1 = NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, missing_taxa_list, method = "uniform", seed = 1)

  # No excluded subclades in the list: warning, and not checked
  expect_warning(check_grafting_constraints(grafted, toy$phy, list(Taxon_name = "X1", Target_clade = list(X1 = toy$A)), verbose = FALSE),
                 "excluded subclades are not checked")
  # Inconsistent lengths
  expect_error(check_grafting_constraints(grafted, toy$phy, list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A),
                                                                 Excluded_subclades = list(X1 = NULL)), verbose = FALSE),
               "same length")
  # Tips missing from the backbone
  expect_error(check_grafting_constraints(grafted, toy$phy, list(Taxon_name = "X1", Target_clade = list(X1 = c("zz", "yy")),
                                                                 Excluded_subclades = list(X1 = NULL)), verbose = FALSE),
               "missing from the backbone")

  # stop_on_failure
  outside <- graft_tip_on_branch(toy$phy, "X1", match(toy$o, toy$phy$tip.label), 0.5)
  expect_error(check_grafting_constraints(outside, toy$phy, missing_taxa_list, stop_on_failure = TRUE, verbose = FALSE))
})


## Test that the items must be named after Taxon_name, and that the taxa are reported by name
test_that("check_grafting_constraints() requires named items, and reports the failed taxa and violated excluded subclades by name",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A),
                         Excluded_subclades = list(X1 = list(sub_B = toy$B), X2 = NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", seed = 1)
  checks <- check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)
  expect_equal(checks$Taxon_name, c("X1", "X2"))

  # Names are required for Target_clade, and for the other items if they are provided
  unnamed <- mtl; names(unnamed$Target_clade) <- NULL
  expect_error(check_grafting_constraints(grafted, toy$phy, unnamed, verbose = FALSE), "Target_clade.*must be named after 'Taxon_name'")
  wrong <- mtl; names(wrong$Excluded_subclades) <- c("X2", "X1")
  expect_error(check_grafting_constraints(grafted, toy$phy, wrong, verbose = FALSE), "Excluded_subclades.*not in the same order")
  expect_no_error(check_grafting_constraints(grafted, toy$phy, mtl[c("Taxon_name", "Target_clade")], use_excluded_subclades = FALSE, verbose = FALSE))

  # X1 grafted within its excluded subclade (named 'sub_B'), X2 outside of its target clade
  bad <- graft_tip_on_branch(toy$phy, "X1", match(toy$tipsB[1], toy$phy$tip.label), 0.5)
  bad <- graft_tip_on_branch(bad, "X2", match(toy$o, bad$tip.label), 0.5)
  checks <- check_grafting_constraints(bad, toy$phy, mtl, verbose = FALSE)
  expect_equal(checks$excluded_violated[checks$Taxon_name == "X1"], "sub_B") # names of the violated subclades (indices if unnamed)
  expect_false(checks$target_ok[checks$Taxon_name == "X2"])

  # The failed taxa are listed in the summary and in the error message
  expect_output(check_grafting_constraints(bad, toy$phy, mtl, verbose = TRUE), "Failed taxa:.*X1: .*sub_B.*X2: outside its target clade", )
  expect_error(check_grafting_constraints(bad, toy$phy, mtl, stop_on_failure = TRUE, verbose = FALSE), "X1: .*sub_B")
  # Nothing is listed if all the taxa respect their constraints
  expect_no_match(capture.output(check_grafting_constraints(grafted, toy$phy, mtl, verbose = TRUE)), "Failed taxa")
})

## Test the check of the allowed states with the updated simmap of the grafted phylogeny
test_that("check_grafting_constraints() checks the allowed states if the grafted phylogeny has the updated simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  taxa <- paste0("X", 1:6)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(toy$A), 6), taxa), Excluded_subclades = stats::setNames(rep(list(NULL), 6), taxa))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = 1)
  taxon_states <- vapply(taxa, function (taxon) names(grafted$simmap$maps[[match(match(taxon, grafted$tip.label), grafted$edge[, 2])]]), character(1))

  # The allowed states used for grafting are respected
  checks <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = c("A", "B"), verbose = FALSE)
  expect_true(all(checks$states_ok))
  expect_true(all(checks$ok))
  expect_true(attr(checks, "all_ok"))
  # Without allowed states, all the states are allowed
  checks_all <- check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)
  expect_true(all(checks_all$states_ok))

  # Other states: the taxa grafted in a state that is not allowed are flagged
  other <- setdiff(c("A", "B", "C"), unique(taxon_states))
  forbidden <- unique(taxon_states)[1]
  allowed_other <- setdiff(c("A", "B", "C"), forbidden)
  checks_other <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = allowed_other, verbose = FALSE)
  expect_identical(unname(checks_other$states_ok), unname(taxon_states != forbidden))
  expect_identical(unname(checks_other$ok), unname(taxon_states != forbidden))
  expect_false(attr(checks_other, "all_ok"))
  expect_match(paste(capture.output(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = allowed_other, verbose = TRUE)), collapse = "\n"),
               "Allowed states violated: +[1-9].*not allowed \\(")
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = allowed_other, stop_on_failure = TRUE, verbose = FALSE), "not allowed")

  # The taxon-level states replace the global ones
  mtl_taxon <- mtl
  mtl_taxon$allowed_states <- stats::setNames(rep(list(NULL), 6), taxa)
  mtl_taxon$allowed_states[[taxa[1]]] <- setdiff(c("A", "B", "C"), taxon_states[1])
  checks_taxon <- check_grafting_constraints(grafted, toy$phy, mtl_taxon, allowed_states = c("A", "B", "C"), verbose = FALSE)
  expect_false(checks_taxon$states_ok[1])
  expect_true(all(checks_taxon$states_ok[-1]))

  # Unknown states are listed with the states of the simmap
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "red", verbose = FALSE), "'red'.*'A', 'B', 'C'")

  # Without the updated simmap, the states are not checked: NA, and the function stops if allowed states are provided
  without_simmap <- drop_simmap(grafted)
  expect_error(check_grafting_constraints(without_simmap, toy$phy, mtl, allowed_states = "A", verbose = FALSE), "'allowed_states' can only be used if 'phy_grafted' has a '\\$simmap' item")
  checks_ns <- check_grafting_constraints(without_simmap, toy$phy, mtl, verbose = FALSE)
  expect_true(all(is.na(checks_ns$states_ok)))
  expect_true(all(checks_ns$ok))
  expect_no_warning(check_grafting_constraints(without_simmap, toy$phy, mtl, verbose = FALSE))
  expect_false(grepl("Allowed states violated", paste(capture.output(check_grafting_constraints(without_simmap, toy$phy, mtl, verbose = TRUE)), collapse = "\n")))

  # An invalid simmap item
  bad <- grafted
  bad$simmap <- list(1)
  expect_error(check_grafting_constraints(bad, toy$phy, mtl, verbose = FALSE), "must be a stochastic character map")
})


# ---- Biogeographic ranges as states (treat_states_as_ranges) ----

## Test that the argument is checked, and that single-area states give the same checks as the exact matching
test_that("check_grafting_constraints() checks treat_states_as_ranges and gives the same results as the exact matching for single-area states",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  taxa <- paste0("X", 1:6)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(toy$A), 6), taxa), Excluded_subclades = stats::setNames(rep(list(NULL), 6), taxa))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = "uniform", seed = 1)

  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, treat_states_as_ranges = "yes", verbose = FALSE), "'treat_states_as_ranges' must be TRUE or FALSE")
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, treat_states_as_ranges = NA, verbose = FALSE), "'treat_states_as_ranges' must be TRUE or FALSE")

  # Default: the column ranges_ok is not checked
  checks_exact <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = c("A", "B"), verbose = FALSE)
  expect_true(all(is.na(checks_exact$ranges_ok)))
  # With ranges, a single-area range inherited from a state is the state itself, and is contained in the state of the graft point
  checks_ranges <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_true(all(checks_ranges$ranges_ok))
  expect_identical(checks_ranges$states_ok, checks_exact$states_ok)
  expect_identical(checks_ranges$ok, checks_exact$ok)
  expect_true(attr(checks_ranges, "all_ok"))
  # Other allowed states flag the same taxa
  other_exact <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "C", verbose = FALSE)
  other_ranges <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "C", treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_identical(other_ranges$states_ok, other_exact$states_ok)
  expect_false(attr(other_ranges, "all_ok"))
})

## Test that a phylogeny grafted with ranges passes both checks, with and without the allowed ranges
test_that("check_grafting_constraints() validates taxa grafted with ranges, including ranges that are subsets of the states",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy) # States: A, AB, BC
  taxa <- paste0("X", 1:6)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(toy$B), 6), taxa), Excluded_subclades = stats::setNames(rep(list(NULL), 6), taxa))
  # Clade B is entirely in the state AB: range 'A' is a subset of it (subset sympatric speciation), range 'AB' is the state itself
  for (allowed in list("A", "AB", c("A", "B"), c("A", "AB"), NULL))
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = allowed, treat_states_as_ranges = TRUE, method = "uniform", seed = 2)
    checks <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = allowed, treat_states_as_ranges = TRUE, verbose = FALSE)
    expect_true(all(checks$states_ok))
    expect_true(all(checks$ranges_ok))
    expect_true(all(checks$ok))
    expect_true(attr(checks, "all_ok"))
  }

  # The ranges inherited are the allowed ranges ('A' here), and not the state of the branch
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, method = "uniform", seed = 2)
  expect_true(all(graft_states_in_simmap(grafted, taxa)$state == "A"))
  # The allowed ranges are not the ranges inherited: the taxa inherited 'A', so 'AB' is not accepted
  checks_ab <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "AB", treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_false(any(checks_ab$states_ok))
  expect_true(all(checks_ab$ranges_ok))
  # Any range is accepted without allowed states
  expect_true(all(check_grafting_constraints(grafted, toy$phy, mtl, treat_states_as_ranges = TRUE, verbose = FALSE)$states_ok))
  # The taxon-level allowed ranges replace the global ones
  mtl_taxon <- mtl
  mtl_taxon$allowed_states <- stats::setNames(rep(list(NULL), 6), taxa)
  mtl_taxon$allowed_states[["X1"]] <- "AB"
  checks_taxon <- check_grafting_constraints(grafted, toy$phy, mtl_taxon, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_identical(unname(checks_taxon$states_ok), c(FALSE, rep(TRUE, 5)))
})

## Test that the range inherited must be contained in the state of the lineage at the graft point
test_that("check_grafting_constraints() flags the taxa that inherited a range not contained in the state of the graft point",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy)
  taxa <- paste0("X", 1:6)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(toy$B), 6), taxa), Excluded_subclades = stats::setNames(rep(list(NULL), 6), taxa))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, method = "uniform", seed = 2)
  tip_edge <- match(match("X1", grafted$tip.label), grafted$edge[, 2])

  # The terminal branch of X1 is changed to 'BC', which is allowed but not contained in the state 'AB' at the graft point
  bad <- grafted
  names(bad$simmap$maps[[tip_edge]]) <- "BC"
  checks <- check_grafting_constraints(bad, toy$phy, mtl, allowed_states = c("A", "BC"), treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_true(checks$states_ok[1])
  expect_false(checks$ranges_ok[1])
  expect_false(checks$ok[1])
  expect_true(all(checks$ranges_ok[-1]))
  expect_true(all(checks$ok[-1]))
  expect_false(attr(checks, "all_ok"))
  # Both checks fail if 'BC' is not allowed
  checks_both <- check_grafting_constraints(bad, toy$phy, mtl, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_false(checks_both$states_ok[1])
  expect_false(checks_both$ranges_ok[1])
  # Matching the states exactly, the range is not checked
  expect_true(is.na(check_grafting_constraints(bad, toy$phy, mtl, allowed_states = c("A", "BC"), verbose = FALSE)$ranges_ok[1]))

  # The failure is described in the summary and in the error message
  summary_text <- paste(capture.output(check_grafting_constraints(bad, toy$phy, mtl, allowed_states = c("A", "BC"), treat_states_as_ranges = TRUE, verbose = TRUE)), collapse = "\n")
  expect_match(summary_text, "Range not in graft-point state: +1")
  expect_match(summary_text, "X1: inherited a range that is not contained in the state of the graft point \\(BC not in [A-Z]+\\)")
  expect_error(check_grafting_constraints(bad, toy$phy, mtl, allowed_states = c("A", "BC"), treat_states_as_ranges = TRUE, stop_on_failure = TRUE, verbose = FALSE),
               "X1: inherited a range that is not contained in the state of the graft point")
  expect_match(paste(capture.output(check_grafting_constraints(bad, toy$phy, mtl, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = TRUE)), collapse = "\n"),
               "X1: inherited a range that is not allowed \\(BC; allowed: A\\); inherited a range that is not contained")
  # The line is not printed without ranges
  expect_no_match(paste(capture.output(check_grafting_constraints(bad, toy$phy, mtl, verbose = TRUE)), collapse = "\n"), "Range not in graft-point state")
})

## Test that the state at the graft point is searched above the missing-only branches that are in the range inherited
test_that("check_grafting_constraints() goes up the branches leading to missing taxa only to find the state of the graft point",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy)
  taxa <- paste0("X", 1:8)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(toy$B), 8), taxa), Excluded_subclades = stats::setNames(rep(list(NULL), 8), taxa))
  n_tips <- ape::Ntip(toy$phy)
  # Find a grafting where a taxon is attached on a branch leading to missing taxa only (grafted on the branch of another missing taxon)
  found <- NULL
  for (seed in 1:40)
  {
    grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, method = "uniform", seed = seed)
    tip_edges <- match(match(taxa, grafted$tip.label), grafted$edge[, 2])
    parent_edges <- match(grafted$edge[tip_edges, 1], grafted$edge[, 2])
    leads_to_missing <- vapply(grafted$edge[parent_edges, 2], function (node)
    {
      tips <- grafted$tip.label[phytools::getDescendants(grafted, node)]
      !any(tips %in% toy$phy$tip.label)
    }, logical(1))
    if (any(leads_to_missing)) { found <- grafted; break }
  }
  expect_false(is.null(found))
  grafted <- found
  expect_true(all(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE)$ranges_ok))

  # The state of the graft point is changed to 'BC' on the backbone branches above the graft nodes, so the range 'A' is not contained
  # in it for the taxa grafted directly on them, and also for the taxa grafted on the branches of these taxa (only found by going up)
  bad <- grafted
  for (e in seq_len(nrow(bad$edge)))
  {
    tips <- bad$tip.label[phytools::getDescendants(bad, bad$edge[e, 2])]
    if (bad$edge[e, 2] <= n_tips) { tips <- bad$tip.label[bad$edge[e, 2]] }
    if (any(tips %in% toy$phy$tip.label) && any(tips %in% taxa))
    {
      map <- bad$simmap$maps[[e]]
      if (names(map)[length(map)] == "AB") { names(bad$simmap$maps[[e]])[length(map)] <- "BC" }
    }
  }
  checks <- check_grafting_constraints(bad, toy$phy, mtl, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_true(all(checks$states_ok))
  expect_false(any(checks$ranges_ok))
  expect_false(attr(checks, "all_ok"))
})

## Test that a taxon grafted at the root is not checked
test_that("check_grafting_constraints() does not check the range of a taxon grafted at the root",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy)
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$phy$tip.label), Excluded_subclades = list(X1 = NULL))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, seed = 1)
  tip_edge <- match(match("X1", grafted$tip.label), grafted$edge[, 2])
  at_root <- grafted$edge[tip_edge, 1] == ape::Ntip(grafted) + 1
  checks <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE)
  expect_true(checks$states_ok)
  if (at_root) { expect_true(is.na(checks$ranges_ok)) } else { expect_true(checks$ranges_ok) }
  expect_true(checks$ok)
})

## Test the validation of the states and the allowed ranges, and the case without simmap
test_that("check_grafting_constraints() validates the coding of the ranges, and ignores treat_states_as_ranges without simmap",
{
  toy <- toy_setup()
  range_simmap <- make_toy_range_simmap(toy$phy)
  taxa <- paste0("X", 1:3)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(toy$B), 3), taxa), Excluded_subclades = stats::setNames(rep(list(NULL), 3), taxa))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = range_simmap, allowed_states = "A", treat_states_as_ranges = TRUE, seed = 1)

  # Allowed ranges that are incorrectly coded, or with areas that are not in the simmap
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "ba", treat_states_as_ranges = TRUE, verbose = FALSE), "the global 'allowed_states' argument.*'ba'")
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "AD", treat_states_as_ranges = TRUE, verbose = FALSE), "'D'.*Areas available in the simmap: 'A', 'B', 'C'")
  mtl_taxon <- mtl
  mtl_taxon$allowed_states <- stats::setNames(list("A", "red", NULL), taxa)
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl_taxon, treat_states_as_ranges = TRUE, verbose = FALSE), "taxon 'X2'.*'red'")
  # States of the simmap that are not ranges
  bad_states <- grafted
  names(bad_states$simmap$maps[[1]])[1] <- "ab"
  colnames(bad_states$simmap$mapped.edge)[1] <- "ab"
  expect_error(check_grafting_constraints(bad_states, toy$phy, mtl, treat_states_as_ranges = TRUE, verbose = FALSE), "the simmap of 'phy_grafted'.*'ab'")
  # Matching the states exactly, ranges are just states
  expect_no_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "A", verbose = FALSE))

  # Without simmap, the function stops early if treat_states_as_ranges = TRUE (and/or allowed states are provided), and the columns are NA otherwise
  without_simmap <- drop_simmap(grafted)
  expect_error(check_grafting_constraints(without_simmap, toy$phy, mtl, treat_states_as_ranges = TRUE, verbose = FALSE), "'treat_states_as_ranges' can only be used if 'phy_grafted' has a '\\$simmap' item")
  expect_error(check_grafting_constraints(without_simmap, toy$phy, mtl, allowed_states = "A", treat_states_as_ranges = TRUE, verbose = FALSE), "'allowed_states' and 'treat_states_as_ranges' can only be used")
  checks <- check_grafting_constraints(without_simmap, toy$phy, mtl, treat_states_as_ranges = FALSE, verbose = FALSE)
  expect_true(all(is.na(checks$states_ok)))
  expect_true(all(is.na(checks$ranges_ok)))
  expect_true(all(checks$ok))
  expect_no_warning(check_grafting_constraints(without_simmap, toy$phy, mtl, verbose = FALSE))
  expect_no_match(paste(capture.output(check_grafting_constraints(without_simmap, toy$phy, mtl, verbose = TRUE)), collapse = "\n"), "Range not in graft-point state")
})
