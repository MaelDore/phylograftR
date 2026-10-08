
##### Test suite for identify_grafted() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


# ---- Tree navigation helpers ----

## Test that get_clade_node() / get_parent_node() / get_descendant_nodes() navigate the tree properly
test_that("get_clade_node() / get_parent_node() / get_descendant_nodes() navigate the tree",
{
  toy <- toy_setup()
  phy <- toy$phy
  n <- ape::Ntip(phy)

  # One tip = tip node; a pair of tips = MRCA node
  expect_identical(get_clade_node(phy, toy$s), match(toy$s, phy$tip.label))
  expect_equal(get_clade_node(phy, toy$A), toy$node_A)

  # Parent: NA for the root, and the parent of a tip is its rootward node
  root <- n + 1
  expect_true(is.na(get_parent_node(phy, root)))
  expect_equal(get_parent_node(phy, toy$node_B), toy$node_A)

  # Descendants of a clade: all its nodes (internal nodes and tips), and nothing else
  desc <- get_descendant_nodes(phy, toy$node_B)
  expect_equal(sum(desc <= n), toy$nB)
  expect_false(toy$node_B %in% desc)
  expect_length(get_descendant_nodes(phy, match(toy$s, phy$tip.label)), 0)
})

# ---- Stem region ----

## Test that get_stem_region_nodes() identify the proper backbone parent node
test_that("get_stem_region_nodes() follows the stem pieces and the missing-only subtrees, up to the backbone parent",
{
  toy <- toy_setup()
  phy <- toy$phy
  node_B <- toy$node_B
  parent_B <- get_parent_node(phy, node_B)
  missing_labels <- c("X1", "X2", "X3")

  # No graft: the stem region is the stem branch of the clade, and the top parent is its backbone parent
  r0 <- get_stem_region_nodes(phy, node_B, missing_labels)
  expect_equal(as.integer(r0$nodes), as.integer(node_B))
  expect_equal(r0$top_parent, parent_B)

  # X1 grafted on the stem: 2 stem pieces + the tip branch of X1
  g1 <- graft_tip_on_branch(phy, "X1", node_B, 0.5)
  r1 <- get_stem_region_nodes(g1, ape::getMRCA(g1, toy$B), missing_labels)
  expect_length(r1$nodes, 3)
  expect_true(match("X1", g1$tip.label) %in% r1$nodes)
  # The top parent is the backbone parent (found through its tips, because node IDs change after grafting)
  parent_tips <- g1$tip.label[ape::prop.part(g1)[[r1$top_parent - ape::Ntip(g1)]]]
  expect_setequal(parent_tips, c(phy$tip.label[ape::prop.part(phy)[[parent_B - ape::Ntip(phy)]]], "X1"))

  # X2 grafted on X1: the missing-only cherry (X1, X2) is part of the stem region
  g2 <- graft_tip_on_branch(g1, "X2", match("X1", g1$tip.label), 0.5)
  r2 <- get_stem_region_nodes(g2, ape::getMRCA(g2, toy$B), missing_labels)
  expect_true(all(match(c("X1", "X2"), g2$tip.label) %in% r2$nodes))
  expect_length(r2$nodes, 5) # 2 stem pieces + branch to the cherry + 2 tip branches

  # The root has no stem
  expect_length(get_stem_region_nodes(phy, ape::Ntip(phy) + 1)$nodes, 0)
})

## Test that get_available_branch_nodes() deals with stem inclusion/exclusion and single tips properly
test_that("get_available_branch_nodes() adds the stem region only if requested, and always for single tips",
{
  toy <- toy_setup()
  phy <- toy$phy
  with_stem <- get_available_branch_nodes(phy, toy$node_B, include_target_clade_stem = TRUE)
  without_stem <- get_available_branch_nodes(phy, toy$node_B, include_target_clade_stem = FALSE)
  expect_true(toy$node_B %in% with_stem)
  expect_false(toy$node_B %in% without_stem)
  expect_setequal(setdiff(with_stem, without_stem), toy$node_B)

  # A single tip: its terminal branch is its stem, always available
  tip_node <- match(toy$s, phy$tip.label)
  expect_identical(get_available_branch_nodes(phy, tip_node, include_target_clade_stem = FALSE), tip_node)
})

# ---- Sampling and seeds ----

## Test that sample_one() handle edge cases properly
test_that("sample_one() samples a single element and stops on an empty vector",
{
  set.seed(1)
  expect_true(sample_one(c(5, 8, 13)) %in% c(5, 8, 13))
  expect_identical(sample_one(7), 7) # no confusion with sample(1:7, 1)
  expect_error(sample_one(numeric(0)))
})

## Test that set_local_seed() ensure reproducibility wihtout altering the global .Random.seed
test_that("set_local_seed() makes a function reproducible without altering the global RNG stream",
{
  with_seed <- function (seed) { set_local_seed(seed); stats::runif(1) }

  set.seed(99)
  expected_next <- stats::runif(1) # first draw of the stream
  set.seed(99)
  a <- with_seed(7)
  b <- with_seed(7)
  expect_identical(a, b)
  # The global stream is untouched: the next draw is the one expected after set.seed(99)
  expect_identical(stats::runif(1), expected_next)

  # Restored even after an error
  failing <- function () { set_local_seed(3); stop("boom") }
  set.seed(99)
  try(failing(), silent = TRUE)
  expect_identical(stats::runif(1), expected_next)

  # seed = NULL does nothing
  expect_null(set_local_seed(NULL))
})

# ---- Divergence times from TreeSim::corsim() ----

## Test that new divergence times produced by TreeSim::corsim() are properly extracted
test_that("extract_new_divergence_times() removes the input times from the corsim output",
{
  old_times <- c(5, 3, 1)
  expect_equal(extract_new_divergence_times(c(5, 4, 3, 2.5, 1), old_times), c(4, 2.5))
  # A new time equal to an old time is kept once
  expect_equal(extract_new_divergence_times(c(5, 3, 3, 1), old_times), 3)
  # Single old time
  expect_equal(extract_new_divergence_times(c(4, 2), 4), 2)
  # On real corsim output
  set.seed(1)
  all_times <- TreeSim::corsim(x = c(10, 6, 3), lambda = 0.3, mu = 0.05, missing = 2, told = 15)
  new_times <- extract_new_divergence_times(all_times, c(10, 6, 3))
  expect_length(new_times, 2)
  expect_true(all(new_times <= 15))
})

# ---- Birth-death fit ----

## Test fit_bd_on_clade() behavior regarding enclosing clades and stops if reaching the root
test_that("fit_bd_on_clade() climbs to an enclosing clade with enough tips, and stops if the backbone is too small",
{
  toy <- toy_setup()
  # B is a small clade: the fit uses a larger enclosing clade
  # Without richness information, nothing is pending nor unsampled: sampling fraction of 1
  fit <- fit_bd_on_clade(toy$phy, ape::getMRCA(toy$phy, toy$B), min_tips = 6)
  expect_gte(fit$n_tips, 6)
  expect_true(fit$b > 0 && fit$d >= 0 && fit$d < fit$b)
  expect_equal(fit$n_pending, 0)
  expect_equal(fit$rho, 1)

  # The 2 pending taxa nested in B are counted in the enclosing clade
  richness_B <- build_richness_info(toy$phy, c("b1", "b2"), list(toy$B, toy$B), "target")
  fit_pending <- fit_bd_on_clade(toy$phy, ape::getMRCA(toy$phy, toy$B), min_tips = 6, richness_info = richness_B)
  expect_equal(fit_pending$n_pending, 2)
  expect_equal(fit_pending$rho, fit_pending$n_tips / (fit_pending$n_tips + 2))

  # A clade with enough tips is used as it is
  fit_A <- fit_bd_on_clade(toy$phy, ape::getMRCA(toy$phy, toy$A), min_tips = 5)
  expect_equal(fit_A$n_tips, toy$nA)

  # A tip: starts from its parent
  fit_tip <- fit_bd_on_clade(toy$phy, match(toy$s, toy$phy$tip.label), min_tips = 5)
  expect_gte(fit_tip$n_tips, 5)

  # Backbone too small
  tiny <- make_toy_phy(n_tips = 4)
  expect_error(fit_bd_on_clade(tiny, ape::getMRCA(tiny, tiny$tip.label[1:2]), min_tips = 5), "too small")
})

## Test that fit_bd_on_clade() accounts for the pending taxa and the unsampled species when richness information is given
test_that("fit_bd_on_clade() accounts for the pending taxa and the unsampled species when richness information is given",
{
  toy <- toy_setup()
  taxa <- c("b1", "b2", "a1")
  clades <- list(toy$B, toy$B, toy$A)
  richness <- build_richness_info(toy$phy, taxa, clades, "target",
                                  list(Target_clade = list(toy$B, toy$A),
                                       total_richness = c(toy$nB + 2 + 4, toy$nA + 3 + 10)))
  node_A <- ape::getMRCA(toy$phy, toy$A)

  # Fit on A: 3 pending taxa; unsampled species of A = 10 + ... (A includes B)
  fit <- fit_bd_on_clade(toy$phy, node_A, min_tips = 5, richness_info = richness)
  expect_equal(fit$n_pending, 3)
  expect_equal(fit$n_unsampled, 10)
  expect_equal(fit$rho, toy$nA / (toy$nA + 3 + 10))

  # Taxa already grafted are in the tree (present): they are not pending anymore
  grafted <- graft_tip_on_branch(toy$phy, "b1", ape::getMRCA(toy$phy, toy$B) , 0.5)
  fit_after <- fit_bd_on_clade(grafted, ape::getMRCA(grafted, toy$A), min_tips = 5, richness_info = richness)
  expect_equal(fit_after$n_pending, 2)
  expect_equal(fit_after$n_tips, toy$nA + 1)
  # The denominator always equals the total richness of the clade
  expect_equal(fit_after$n_tips + fit_after$n_pending + fit_after$n_unsampled, toy$nA + 3 + 10)
})


# ---- Richness information ----

## Build example for richness_info. Missing taxa: 3 in B, 2 in C, 4 in A only
richness_example <- function (toy, total = NULL, mode = "target")
{
  taxa <- c(paste0("b", 1:3), paste0("c", 1:2), paste0("a", 1:4))
  clades <- c(rep(list(toy$B), 3), rep(list(toy$C), 2), rep(list(toy$A), 4))
  build_richness_info(toy$phy, taxa, clades, mode, total)
}

## Test that build_richness_info() derives the unsampled species of each clade, without counting nested clades twice
test_that("build_richness_info() derives the unsampled species of each clade, without counting nested clades twice",
{
  toy <- toy_setup()
  total_B <- toy$nB + 3 + 12
  total_A <- toy$nA + 9 + 31
  info <- richness_example(toy, list(Target_clade = list(toy$B, toy$C, toy$A),
                                     total_richness = c(total_B, NA, total_A)))
  tb <- info$table

  expect_equal(tb$nb_backbone_species, c(toy$nB, toy$nC, toy$nA))
  expect_equal(tb$nb_missing_species, c(3, 2, 9)) # missing taxa nested in the clade, nested clades included
  expect_equal(tb$nb_unsampled_species, c(12, 0, 31))
  expect_equal(tb$total_richness_provided, c(TRUE, FALSE, TRUE))
  # NA: no unsampled species, total derived
  expect_equal(tb$total_richness[2], toy$nC + 2)
  expect_equal(tb$sampling_fraction, (tb$nb_backbone_species + tb$nb_missing_species) / tb$total_richness)

  # unsampled species per node
  expect_equal(info$unsampled_by_node[c(toy$node_B, toy$node_C, toy$node_A)], c(12, 0, 31))
})

## Test that build_richness_info() derives the unsampled species of clades based on information from its nested clades
test_that("build_richness_info() derives the unsampled species of clades without a total from the clades nested within them",
{
  toy <- toy_setup()
  clades <- list(Target_clade = list(toy$B, toy$C, toy$A))

  # NA for A: sum of the nested clades with a total (here B only)
  info <- richness_example(toy, c(clades, list(total_richness = c(toy$nB + 3 + 12, NA, NA))))
  expect_equal(info$unsampled_by_node[toy$node_A], 12)
  expect_equal(info$table$nb_unsampled_species[3], 12)

  # Sum across disjoint nested clades (B and C)
  info2 <- richness_example(toy, c(clades, list(total_richness = c(toy$nB + 3 + 12, toy$nC + 2 + 5, NA))))
  expect_equal(info2$unsampled_by_node[toy$node_A], 17)

  # A node above A that is not listed uses the largest informed clades below it
  info3 <- richness_example(toy, c(clades, list(total_richness = c(toy$nB + 3 + 12, NA, toy$nA + 9 + 31))))
  parent_A <- get_parent_node(toy$phy, toy$node_A)
  expect_equal(info3$unsampled_by_node[parent_A], 31)
})

## Test that build_richness_info() only count grafted taxa from sister clades that are nested in the focal clade
test_that("build_richness_info() in 'sister' mode only counts the taxa whose sister clade is strictly nested",
{
  toy <- toy_setup()
  info <- richness_example(toy, list(Sister_clade = list(toy$B, toy$C, toy$A),
                                     total_richness = c(toy$nB + 10, toy$nC, toy$nA + 5 + 20)), mode = "sister")
  # Taxa sister to B or C are on their stems: not in B or C, but in A
  expect_equal(info$table$nb_missing_species, c(0, 0, 5))
})

## Test that uninformative total_richness_per_clades input (i.e., vector of NA) returns a table with no unsampled species; and no total_richness_per_clades input returns no table
test_that("build_richness_info() without totals ('NULL') returns no table; with only 'NA' returns a table with no unsampled species",
{
  toy <- toy_setup()
  info <- richness_example(toy)
  expect_null(info$table)
  expect_true(all(info$unsampled_by_node == 0))
  # A vector with only NA (logical) is accepted
  info_na <- richness_example(toy, list(Target_clade = list(toy$B, toy$C, toy$A), total_richness = c(NA, NA, NA)))
  expect_true(all(info_na$table$nb_unsampled_species == 0))
  expect_false(any(info_na$table$total_richness_provided))
})

## Test validity checks for build_richness_info() when totals are not valid
test_that("build_richness_info() stops with informative messages for invalid totals",
{
  toy <- toy_setup()
  ok_clades <- list(toy$B, toy$C, toy$A)
  total_B <- toy$nB + 3 + 12
  total_A <- toy$nA + 9 + 31

  # A clade used by the missing taxa is not listed
  expect_error(richness_example(toy, list(Target_clade = list(toy$B, toy$A), total_richness = c(total_B, total_A))),
               "not listed in 'total_richness_per_clades")
  # Total lower than the backbone + missing species
  expect_error(richness_example(toy, list(Target_clade = ok_clades, total_richness = c(toy$nB + 2, NA, total_A))),
               "lower than the number of backbone species")
  # Parent total too low to contain a nested clade with a total
  expect_error(richness_example(toy, list(Target_clade = ok_clades, total_richness = c(total_B, NA, toy$nA + 9 + 5))),
               "too low to contain")
  # Same clade listed twice
  expect_error(richness_example(toy, list(Target_clade = c(ok_clades, list(rev(toy$B))), total_richness = c(total_B, NA, total_A, total_B))),
               "more than once")
  # Structure and values
  expect_error(richness_example(toy, list(Target_clade = ok_clades, Total = 1)), "must be a list with two elements")
  expect_error(richness_example(toy, list(Target_clade = ok_clades, total_richness = c(1, 2))), "same length")
  expect_error(richness_example(toy, list(Target_clade = ok_clades, total_richness = c(total_B + 0.5, NA, total_A))), "whole numbers")
  expect_error(richness_example(toy, list(Target_clade = list(toy$B, toy$C, c("zz", "yy")), total_richness = c(total_B, NA, total_A))),
               "missing from the backbone")
})


# ---- Clade grouping and order ----

## Test the ordering of clades
test_that("group_and_order_clades() groups taxa by clade and orders the taxa so nested clades come first",
{
  toy <- toy_setup()
  clades <- list(toy$A, toy$B, toy$A, toy$C, toy$B)
  groups <- group_and_order_clades(toy$phy, clades)

  expect_length(groups$group_ID, 5)
  expect_length(groups$group_clade_tips, 3)
  expect_equal(groups$group_ID, c(1, 2, 1, 3, 2))
  expect_setequal(groups$taxa_order, 1:5)

  # Clade size never decreases along the order, and the input order is kept within a clade (stable)
  size <- vapply(clades[groups$taxa_order], function (tips) {
    length(get_descendant_nodes(toy$phy, get_clade_node(toy$phy, tips))) }, numeric(1))
  expect_false(is.unsorted(size))
  expect_identical(groups$taxa_order[groups$group_ID[groups$taxa_order] == 1], c(1L, 3L)) # the two taxa of A stay in input order
})


## Test grouping and ordering of clades with different stem availability
test_that("group_and_order_clades() splits the taxa of a clade by stem availability, 'FALSE' first, and ignores it for single tips",
{
  toy <- toy_setup()
  clades <- list(toy$A, toy$B, toy$A, toy$A, toy$s)
  include_target_clade_stem <- c(TRUE, TRUE, FALSE, TRUE, FALSE) # FALSE is ignored for the single tip 's'
  groups <- group_and_order_clades(toy$phy, clades, include_target_clade_stem = include_target_clade_stem)

  expect_equal(groups$group_ID, c(1, 2, 3, 1, 4))
  expect_equal(groups$group_clade_ID, c(1, 2, 1, 3)) # the two groups of A share the ID of their clade
  expect_equal(groups$group_include_target_clade_stem, c(TRUE, TRUE, FALSE, TRUE))
  # Single tip, then B (nested in A), then the groups of A with 'FALSE' first
  expect_equal(groups$group_order, c(4, 2, 3, 1))
  expect_equal(groups$taxa_order, c(5, 2, 3, 1, 4))

  # Without 'include_target_clade_stem': stem available for all, so one group per clade
  expect_equal(group_and_order_clades(toy$phy, clades)$group_ID, c(1, 2, 1, 1, 3))
})

## Test the resolution of the stem availability of each taxon
test_that("resolve_include_target_clade_stem() uses the taxon values, and the global value for NA or NULL, and validates its inputs",
{
  expect_equal(resolve_include_target_clade_stem(NULL, TRUE, 3), rep(TRUE, 3))
  expect_equal(resolve_include_target_clade_stem(NULL, FALSE, 3), rep(FALSE, 3))
  expect_equal(resolve_include_target_clade_stem(c(NA, FALSE, TRUE), TRUE, 3), c(TRUE, FALSE, TRUE))
  expect_equal(resolve_include_target_clade_stem(c(NA, FALSE, TRUE), FALSE, 3), c(FALSE, FALSE, TRUE))
  expect_equal(resolve_include_target_clade_stem(rep(NA_integer_, 3), FALSE, 3), rep(FALSE, 3)) # only NA, whatever the type

  expect_error(resolve_include_target_clade_stem(c(TRUE, FALSE), TRUE, 3), "same length as 'Taxon_name'")
  expect_error(resolve_include_target_clade_stem(c(1, 0, 1), TRUE, 3), "logical vector")
  expect_error(resolve_include_target_clade_stem(NULL, NA, 3), "single logical value")
  expect_error(resolve_include_target_clade_stem(NULL, c(TRUE, FALSE), 3), "single logical value")
})

## Test that the taxa without branch are identified with the stem availability of each taxon
test_that("find_taxa_without_allowed_branch() uses the stem availability of each taxon, and the stem exclusion of each excluded subclade",
{
  toy <- toy_setup()
  # The excluded subclade is identical to the target clade
  mtl <- list(Taxon_name = c("no_stem", "stem_open", "stem_closed"), Target_clade = list(toy$B, toy$B, toy$B),
              Excluded_subclades = list(list(toy$B), list(toy$B), list(toy$B)))
  flagged <- find_taxa_without_allowed_branch(toy$phy, mtl, include_target_clade_stem = c(FALSE, TRUE, TRUE),
                                              exclude_subclades_stem = list(TRUE, FALSE, TRUE))
  # No branch without the stem of the target clade; only the stem if the stem of the excluded subclade is not excluded; no branch if it is excluded
  expect_equal(flagged$Taxon_name, c("no_stem", "stem_closed"))

  # A single tip: its terminal branch is the only branch, excluded only if its stem is excluded
  tip <- list(Taxon_name = "x", Target_clade = list(toy$s), Excluded_subclades = list(list(toy$s)))
  expect_equal(find_taxa_without_allowed_branch(toy$phy, tip, exclude_subclades_stem = TRUE)$Taxon_name, "x")
  expect_equal(nrow(find_taxa_without_allowed_branch(toy$phy, tip, exclude_subclades_stem = FALSE)), 0)
})

## Test that the taxa without any segment in their allowed states are identified, with the fraction of the length of their branches in allowed states
test_that("find_taxa_without_allowed_branch() flags the taxa without segment in their allowed states, and returns the valid fraction of the others",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  phy <- label_backbone_nodes(toy$phy)$phy
  tracker <- start_state_tracking(phy, simmap)
  # The clade B is entirely in state "B", and the tip 's' entirely in state "A"
  mtl <- list(Taxon_name = c("in_B", "not_in_B", "tip_A", "tip_B", "whole_A"), Target_clade = list(toy$B, toy$B, toy$s, toy$s, toy$A))
  allowed <- list(in_B = "B", not_in_B = c("A", "C"), tip_A = "A", tip_B = "B", whole_A = c("A", "B", "C"))
  flagged <- find_taxa_without_allowed_branch(phy, mtl, tracker = tracker, allowed_states = allowed)

  expect_identical(flagged$Taxon_name, c("not_in_B", "tip_B"))
  expect_identical(flagged$Cause, c("allowed_states", "allowed_states"))
  valid_fraction <- attr(flagged, "valid_fraction")
  expect_named(valid_fraction, mtl$Taxon_name)
  expect_equal(unname(valid_fraction[c("in_B", "tip_A", "whole_A")]), c(1, 1, 1))
  expect_equal(unname(valid_fraction[c("not_in_B", "tip_B")]), c(0, 0))

  # Fraction of the pool in the allowed states, including the stem of the target clade: independent computation with the table of the states
  A_pool <- c(get_descendant_nodes(phy, toy$node_A), toy$node_A)
  state_table <- build_state_table(phy, simmap, validate_simmap(simmap, phy)$edge_map)
  in_pool <- state_table[state_table$node %in% A_pool, ]
  mtl_A <- list(Taxon_name = "x", Target_clade = list(toy$A))
  flagged_A <- find_taxa_without_allowed_branch(phy, mtl_A, tracker = tracker, allowed_states = list(x = "C"))
  expect_equal(unname(attr(flagged_A, "valid_fraction")["x"]), sum((in_pool$to - in_pool$from)[in_pool$state == "C"]) / sum(in_pool$to - in_pool$from))

  # The branches excluded by the constraints are flagged with their own cause, and have no fraction
  mtl_excl <- list(Taxon_name = "bad", Target_clade = list(toy$B), Excluded_subclades = list(list(toy$B)))
  flagged_excl <- find_taxa_without_allowed_branch(phy, mtl_excl, include_target_clade_stem = FALSE, tracker = tracker, allowed_states = list(bad = "B"))
  expect_identical(flagged_excl$Cause, "excluded")
  expect_true(is.na(attr(flagged_excl, "valid_fraction")["bad"]))

  # Without tracker: no state, no attribute
  flagged_none <- find_taxa_without_allowed_branch(toy$phy, mtl)
  expect_equal(nrow(flagged_none), 0)
  expect_null(attr(flagged_none, "valid_fraction"))
})

## Test the resolution of the stem exclusion of each excluded subclade
test_that("resolve_exclude_subclades_stem() uses the subclade values, and the global value for NA or NULL, and validates its inputs",
{
  excluded <- list(list("a", c("b", "c")), NULL, list("d"))
  resolve <- function (item, global = TRUE) resolve_exclude_subclades_stem(item, global, excluded, 3)

  # No item: global value for each subclade
  expect_equal(resolve(NULL), list(c(TRUE, TRUE), logical(0), TRUE))
  expect_equal(resolve(NULL, FALSE), list(c(FALSE, FALSE), logical(0), FALSE))
  # One value per subclade, one value recycled, NA and NULL use the global value
  expect_equal(resolve(list(c(FALSE, NA), NULL, FALSE)), list(c(FALSE, TRUE), logical(0), FALSE))
  expect_equal(resolve(list(FALSE, NULL, NA)), list(c(FALSE, FALSE), logical(0), TRUE))
  expect_equal(resolve(list(NULL, NULL, NULL), FALSE), list(c(FALSE, FALSE), logical(0), FALSE))
  # A vector with one value per taxon is recycled over the subclades of each taxon
  expect_equal(resolve(c(FALSE, NA, TRUE)), list(c(FALSE, FALSE), logical(0), TRUE))
  # Values for a taxon without excluded subclades are ignored
  expect_equal(resolve(list(NULL, c(TRUE, FALSE), NULL)), list(c(TRUE, TRUE), logical(0), TRUE))

  expect_error(resolve(list(TRUE, NULL)), "same length as 'Taxon_name'")
  expect_error(resolve(list(c(TRUE, FALSE, TRUE), NULL, NULL)), "element 1 has 3 values")
  expect_error(resolve(list(c(1, 0), NULL, NULL)), "logical values")
  expect_error(resolve(NULL, NA), "single logical value")
  expect_error(resolve(NULL, c(TRUE, FALSE)), "single logical value")
  # Without any excluded subclades
  expect_equal(resolve_exclude_subclades_stem(NULL, TRUE, NULL, 2), list(logical(0), logical(0)))
})

## Test the branches of the excluded subclades
test_that("get_excluded_branch_nodes() returns the branches within the subclades, plus their stem region if their stem is excluded",
{
  toy <- toy_setup()
  phy <- toy$phy
  desc_B <- get_descendant_nodes(phy, toy$node_B)
  tip_s <- match(toy$s, phy$tip.label)

  # A clade: crown, plus the stem if excluded
  expect_setequal(get_excluded_branch_nodes(phy, list(toy$B), TRUE), c(desc_B, toy$node_B))
  expect_setequal(get_excluded_branch_nodes(phy, list(toy$B), FALSE), desc_B)

  # A single tip: its terminal branch if its stem is excluded, otherwise nothing
  expect_equal(get_excluded_branch_nodes(phy, list(toy$s), TRUE), tip_s)
  expect_equal(get_excluded_branch_nodes(phy, list(toy$s), FALSE), integer(0))

  # Several subclades with their own setting; empty subclades exclude nothing
  mixed <- get_excluded_branch_nodes(phy, list(toy$B, toy$s, character(0), toy$C), c(TRUE, FALSE, TRUE, FALSE))
  expect_setequal(mixed, c(desc_B, toy$node_B, get_descendant_nodes(phy, toy$node_C)))
  expect_equal(get_excluded_branch_nodes(phy, NULL), integer(0))

  # The stem region follows the grafts on the stem of the excluded subclade: pieces of the stem and branches of the missing-only subtrees
  grafted <- graft_tip_on_branch(phy, "X1", toy$node_B, 0.5)
  stem_region <- get_excluded_branch_nodes(grafted, list(toy$B), TRUE, missing_tip_labels = "X1")
  expect_true(match("X1", grafted$tip.label) %in% stem_region)
  expect_gt(length(stem_region), length(desc_B) + 1)
  # Without declaring the taxon as missing, its subtree is a backbone subtree: the stem region stops there
  expect_false(match("X1", grafted$tip.label) %in% get_excluded_branch_nodes(grafted, list(toy$B), TRUE))
})


# ---- Excluded subclades and available branches ----

## Test flagging of excluded subclades position relative to their target clade in the backbone
test_that("check_excluded_subclades_position() flags excluded subclades that are not nested in the target clade",
{
  toy <- toy_setup()
  missing_taxa_list <- list(
    Taxon_name = c("nested", "disjoint", "identical", "contains"),
    Target_clade = list(toy$A, toy$A, toy$B, toy$B),
    Excluded_subclades = list(list(toy$B), list(toy$D), list(toy$B), list(toy$A)))
  flagged <- check_excluded_subclades_position(toy$phy, missing_taxa_list)

  # The nested exclusion is fine and not reported
  expect_setequal(flagged$Taxon_name, c("disjoint", "identical", "contains"))
  expect_identical(flagged$relation[flagged$Taxon_name == "disjoint"], "disjoint")
  expect_identical(flagged$relation[flagged$Taxon_name == "identical"], "identical")
  expect_identical(flagged$relation[flagged$Taxon_name == "contains"], "contains_target")

  # Nothing to flag without excluded subclades, or with a nested single-tip subclade
  none <- list(Taxon_name = c("x", "y"), Target_clade = list(toy$A, toy$A), Excluded_subclades = list(NULL, list(toy$s)))
  expect_equal(nrow(check_excluded_subclades_position(toy$phy, none)), 0)
})

## Test the position of single-tip excluded subclades
test_that("check_excluded_subclades_position() checks single-tip subclades, unless their stem is not excluded (they exclude nothing)",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("identical", "disjoint", "contains"), Target_clade = list(toy$s, toy$A, toy$s),
              Excluded_subclades = list(list(toy$s), list(toy$o), list(toy$A)))
  flagged <- check_excluded_subclades_position(toy$phy, mtl)
  expect_equal(flagged$relation[match(c("identical", "disjoint", "contains"), flagged$Taxon_name)], c("identical", "disjoint", "contains_target"))

  # Single-tip subclades with a stem that is not excluded exclude nothing: not checked
  mtl$exclude_subclades_stem <- list(FALSE, FALSE, TRUE)
  flagged <- check_excluded_subclades_position(toy$phy, mtl)
  expect_equal(flagged$Taxon_name, "contains")
})

## Test flagging of taxa whose constraints left no valid branch for grafting
test_that("find_taxa_without_allowed_branch() lists the taxa whose branches are all excluded",
{
  toy <- toy_setup()
  missing_taxa_list <- list(
    Taxon_name = c("bad", "tip", "fine"),
    Target_clade = list(toy$B, toy$s, toy$A),
    Excluded_subclades = list(list(toy$B), NULL, list(toy$B)))

  # Without the stem, the whole crown of B is excluded: no branch for 'bad'
  flagged <- find_taxa_without_allowed_branch(toy$phy, missing_taxa_list, include_target_clade_stem = FALSE)
  expect_identical(flagged$Taxon_name, "bad")

  # A single tip always has its terminal branch, and A has branches outside B
  expect_false(any(c("tip", "fine") %in% flagged$Taxon_name))

  # Ignoring the excluded subclades: nothing is flagged
  expect_equal(nrow(find_taxa_without_allowed_branch(toy$phy, missing_taxa_list, include_target_clade_stem = FALSE, use_excluded_subclades = FALSE)), 0)
})


# ---- Temporary file path ----

## Test validity check for 'temp_file_path'
test_that("check_temp_file_path() validates the path of the temporary file",
{
  expect_true(check_temp_file_path(NULL))

  good <- file.path(tempdir(), "temp_phy.rds")
  expect_true(check_temp_file_path(good))

  expect_error(check_temp_file_path(file.path(tempdir(), "temp_phy.txt")), "must end with '.rds'")
  expect_error(check_temp_file_path(file.path(tempdir(), "missing_dir", "temp_phy.rds")), "does not exist")
  expect_error(check_temp_file_path(1), "single character string")
  expect_error(check_temp_file_path(c(good, good)), "single character string")
})


# ---- Names of the items of missing_taxa_list ----

## Test the validation of the structure of missing_taxa_list: items named after Taxon_name
test_that("validate_missing_taxa_list() requires items named after Taxon_name, in the same order, only for the items provided",
{
  valid <- list(Taxon_name = c("X1", "X2"),
                Target_clade = list(X1 = c("t1", "t2"), X2 = "t3"),
                include_target_clade_stem = c(X1 = TRUE, X2 = NA),
                Excluded_subclades = list(X1 = list(c("t1", "t2")), X2 = NULL),
                exclude_subclades_stem = list(X1 = FALSE, X2 = NULL))
  expect_equal(validate_missing_taxa_list(valid), c("X1", "X2"))
  # Optional items that are NULL (or absent) do not need names
  expect_silent(validate_missing_taxa_list(valid[c("Taxon_name", "Target_clade")]))
  expect_silent(validate_missing_taxa_list(c(valid[c("Taxon_name", "Target_clade")], list(include_target_clade_stem = NULL))))
  # Names inside the items (e.g., named excluded subclades) are not required
  expect_silent(validate_missing_taxa_list(valid))

  # Taxon_name
  expect_error(validate_missing_taxa_list("a"), "must be a list")
  expect_error(validate_missing_taxa_list(valid[-1]), "missing or empty")
  expect_error(validate_missing_taxa_list(modifyList(valid, list(Taxon_name = character(0)))), "missing or empty")
  expect_error(validate_missing_taxa_list(modifyList(valid, list(Taxon_name = c("X1", NA)))), "character vector without missing values")
  expect_error(validate_missing_taxa_list(modifyList(valid, list(Taxon_name = 1:2))), "character vector")
  expect_error(validate_missing_taxa_list(modifyList(valid, list(Taxon_name = c("X1", "X1")))), "unique names.*X1")

  # Item defining the clades
  expect_error(validate_missing_taxa_list(valid[-2]), "'missing_taxa_list\\$Target_clade' is missing")
  expect_error(validate_missing_taxa_list(valid, clade_item = "Sister_clade"), "Sister_clade' is missing")

  # Every item provided: same length, names, same names in the same order
  for (item in c("Target_clade", "include_target_clade_stem", "Excluded_subclades", "exclude_subclades_stem"))
  {
    too_long <- valid; too_long[[item]] <- c(valid[[item]], list(X3 = NULL))
    expect_error(validate_missing_taxa_list(too_long), paste0(item, ".*same length as 'Taxon_name' \\(3 vs 2\\)"))

    unnamed <- valid; names(unnamed[[item]]) <- NULL
    expect_error(validate_missing_taxa_list(unnamed), paste0(item, ".*must be named after 'Taxon_name'.*no names"))

    wrong <- valid; names(wrong[[item]]) <- c("X1", "Y2")
    expect_error(validate_missing_taxa_list(wrong), paste0(item, ".*do not match.*position 2: 'Y2' in '", item, "' vs 'X2' in 'Taxon_name'"))

    swapped <- valid; names(swapped[[item]]) <- c("X2", "X1")
    expect_error(validate_missing_taxa_list(swapped), paste0(item, ".*not in the same order"))

    na_name <- valid; names(na_name[[item]]) <- c("X1", NA)
    expect_error(validate_missing_taxa_list(na_name), "must match 'Taxon_name'")
  }

  # Only the items listed are validated
  unnamed_excluded <- valid; names(unnamed_excluded$Excluded_subclades) <- NULL
  expect_silent(validate_missing_taxa_list(unnamed_excluded, optional_items = character(0)))
})

## Test the labels of the excluded subclades in messages
test_that("label_excluded_subclades() uses the optional names of the excluded subclades", 
{
  expect_equal(label_excluded_subclades(list(c("a", "b"), "c")), c("a, b", "c"))
  expect_equal(label_excluded_subclades(list(grp = c("a", "b"), "c")), c("'grp' (a, b)", "c")) # partly named
  expect_equal(label_excluded_subclades(list(g1 = "a", g2 = "b")), c("'g1' (a)", "'g2' (b)"))
  expect_equal(label_excluded_subclades(list(g1 = letters[1:8]), max_tips = 3), "'g1' (a, b, c, ... (8 tips))")
  expect_equal(label_excluded_subclades(list(g1 = letters[1:8]), max_tips = Inf), paste0("'g1' (", paste(letters[1:8], collapse = ", "), ")"))
  expect_equal(label_excluded_subclades(NULL), character(0))
})

## Test that the resolvers name their outputs after the taxa (and the excluded subclades), and use the names in messages
test_that("the resolvers of the stem settings return outputs named after the taxa, and report taxa by name", 
{
  expect_equal(resolve_include_target_clade_stem(c(NA, FALSE), TRUE, 2, taxon_names = c("X1", "X2")), c(X1 = TRUE, X2 = FALSE))
  expect_equal(resolve_include_target_clade_stem(NULL, FALSE, 2, taxon_names = c("X1", "X2")), c(X1 = FALSE, X2 = FALSE))
  expect_null(names(resolve_include_target_clade_stem(NULL, FALSE, 2))) # no names without taxon names

  excluded <- list(X1 = list(grp_B = c("a", "b"), "c"), X2 = NULL)
  out <- resolve_exclude_subclades_stem(list(X1 = c(FALSE, NA), X2 = NULL), TRUE, excluded, 2, taxon_names = c("X1", "X2"))
  expect_named(out, c("X1", "X2"))
  expect_equal(out$X1, stats::setNames(c(FALSE, TRUE), c("grp_B", ""))) # names of the excluded subclades are carried
  expect_equal(out$X2, logical(0))
  expect_named(resolve_exclude_subclades_stem(NULL, TRUE, excluded, 2, taxon_names = c("X1", "X2")), c("X1", "X2"))

  # Messages name the taxon
  expect_error(resolve_exclude_subclades_stem(list(X1 = c(TRUE, FALSE, TRUE), X2 = NULL), TRUE, excluded, 2, taxon_names = c("X1", "X2")), "taxon 'X1' has 3 values")
  expect_error(resolve_exclude_subclades_stem(list(X1 = c(1, 0), X2 = NULL), TRUE, excluded, 2, taxon_names = c("X1", "X2")), "logical values.*taxon 'X1'")

  # If both the values and the excluded subclades are named, the names must match
  expect_equal(unname(resolve_exclude_subclades_stem(list(X1 = stats::setNames(c(FALSE, TRUE), c("grp_B", "")), X2 = NULL), TRUE, excluded, 2)[[1]]), c(FALSE, TRUE))
  expect_error(resolve_exclude_subclades_stem(list(X1 = c(grp_C = FALSE, other = TRUE), X2 = NULL), TRUE, excluded, 2), "must match the names of its excluded subclades")
})

# ---- Node labels of the backbone ----

## Test that the internal nodes of the backbone are labelled with their ID, keeping the original labels aside
test_that("label_backbone_nodes() labels the internal nodes with their ID and keeps the original labels",
{
  phy <- make_toy_phy()
  labelled <- label_backbone_nodes(phy)
  expect_equal(labelled$phy$node.label, paste0("bb", ape::Ntip(phy) + seq_len(phy$Nnode)))
  expect_null(labelled$original_node_labels)
  expect_equal(labelled$n_tips, ape::Ntip(phy))

  phy$node.label <- paste0("support", seq_len(phy$Nnode))
  labelled <- label_backbone_nodes(phy)
  expect_equal(labelled$original_node_labels, phy$node.label)
  expect_equal(labelled$phy$node.label, paste0("bb", ape::Ntip(phy) + seq_len(phy$Nnode)))
})

## Test that no label is returned if the backbone had none, and an empty string is given to the new nodes otherwise
test_that("restore_node_labels() gives no labels if the backbone had none, and an empty string to the new nodes otherwise",
{
  phy <- make_toy_phy()
  toy_node <- phy$edge[5, 2]

  # No original labels: no labels after grafting
  labelled <- label_backbone_nodes(phy)
  grafted <- graft_tip_on_branch(labelled$phy, "X1", toy_node)
  expect_null(restore_node_labels(grafted, labelled)$node.label)
  expect_null(restore_node_labels(labelled$phy, labelled)$node.label)

  # Original labels: kept for the backbone nodes, empty string for the nodes created by grafting (whatever the label given by the binding function)
  phy$node.label <- paste0("support", seq_len(phy$Nnode))
  labelled <- label_backbone_nodes(phy)
  grafted <- graft_tip_on_branch(labelled$phy, "X1", toy_node)
  grafted <- graft_tip_on_branch(grafted, "X2", match("X1", grafted$tip.label))
  restored <- restore_node_labels(grafted, labelled)
  expect_length(restored$node.label, restored$Nnode)
  expect_equal(sum(restored$node.label == ""), 2)
  expect_setequal(restored$node.label[restored$node.label != ""], phy$node.label)
  expect_false(anyNA(restored$node.label))

  # Original labels that are NA or empty are kept as they are
  phy$node.label[c(2, 4)] <- c(NA, "")
  labelled <- label_backbone_nodes(phy)
  restored <- restore_node_labels(graft_tip_on_branch(labelled$phy, "X1", toy_node), labelled)
  expect_identical(sum(is.na(restored$node.label)), 1L)
  expect_identical(sum(restored$node.label == "", na.rm = TRUE), 2L)
})


# ---- Check of the backbone ----

## Test that a binary backbone is accepted, and that polytomies are named in the error with the suggested fix
test_that("check_backbone_binary() accepts a binary backbone and stops on polytomies, naming the nodes and suggesting ape::multi2di()",
{
  phy <- make_toy_phy()
  expect_true(check_backbone_binary(phy))

  # One polytomy: the node, its number of children and its tips are given
  poly <- make_toy_polytomies(phy, n_collapsed = 1)
  expect_false(ape::is.binary(poly))
  node <- which(tabulate(poly$edge[, 1]) > 2)
  expect_length(node, 1)
  expect_error(check_backbone_binary(poly), "not binary: 1 internal node\\(s\\)")
  expect_error(check_backbone_binary(poly), paste0("node ", node, ": 3 children"))
  expect_error(check_backbone_binary(poly), "ape::multi2di")
  expect_error(check_backbone_binary(poly), "own risk")
  descendants <- get_descendant_nodes(poly, node)
  expect_error(check_backbone_binary(poly), poly$tip.label[min(descendants[descendants <= ape::Ntip(poly)])], fixed = TRUE)

  # The suggested fix gives a backbone that is accepted
  expect_true(check_backbone_binary(ape::multi2di(poly)))
})

## Test that the list of nodes is truncated, that node labels are shown and that nodes with a single child are reported
test_that("check_backbone_binary() truncates long lists of nodes, shows node labels and reports nodes with a single child",
{
  phy <- make_toy_phy()
  poly <- make_toy_polytomies(phy, n_collapsed = 8)
  n_poly <- sum(tabulate(poly$edge[, 1]) > 2)
  expect_gt(n_poly, 2)
  expect_error(check_backbone_binary(poly, max_nodes = 2), paste0("not binary: ", n_poly, " internal node\\(s\\)"))
  expect_error(check_backbone_binary(poly, max_nodes = 2), paste0("and ", n_poly - 2, " more"))
  expect_error(check_backbone_binary(poly, max_tips = 2), "\\.\\.\\.\\)")

  # Node labels
  poly <- make_toy_polytomies(phy, n_collapsed = 1)
  poly$node.label <- paste0("support", seq_len(poly$Nnode))
  node <- which(tabulate(poly$edge[, 1]) > 2)
  expect_error(check_backbone_binary(poly), paste0("node ", node, " 'support", node - ape::Ntip(poly), "'"))

  # Node with a single child
  single <- phy
  single$edge <- single$edge[-which(single$edge[, 1] == phy$edge[phy$edge[, 2] > ape::Ntip(phy), 1][1])[1], ]
  single$edge.length <- single$edge.length[seq_len(nrow(single$edge))]
  expect_error(check_backbone_binary(single), "1 child")
})
