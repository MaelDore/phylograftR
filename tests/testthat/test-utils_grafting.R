
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
  with_stem <- get_available_branch_nodes(phy, toy$node_B, include_stem = TRUE)
  without_stem <- get_available_branch_nodes(phy, toy$node_B, include_stem = FALSE)
  expect_true(toy$node_B %in% with_stem)
  expect_false(toy$node_B %in% without_stem)
  expect_setequal(setdiff(with_stem, without_stem), toy$node_B)

  # A single tip: its terminal branch is its stem, always available
  tip_node <- match(toy$s, phy$tip.label)
  expect_identical(get_available_branch_nodes(phy, tip_node, include_stem = FALSE), tip_node)
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


# ---- Excluded subclades and available branches ----

## Test flagging of excluded subclades position relative to their target clade in the backbone
test_that("check_excluded_subclades_position() flags excluded subclades that are not nested in the target clade",
{
  toy <- toy_setup()
  missing_taxa_list <- list(
    Taxa_name = c("nested", "disjoint", "identical", "contains"),
    Target_clade = list(toy$A, toy$A, toy$B, toy$B),
    Excluded_subclades = list(list(toy$B), list(toy$D), list(toy$B), list(toy$A)))
  flagged <- check_excluded_subclades_position(toy$phy, missing_taxa_list)

  # The nested exclusion is fine and not reported
  expect_setequal(flagged$Taxa_name, c("disjoint", "identical", "contains"))
  expect_identical(flagged$relation[flagged$Taxa_name == "disjoint"], "disjoint")
  expect_identical(flagged$relation[flagged$Taxa_name == "identical"], "identical")
  expect_identical(flagged$relation[flagged$Taxa_name == "contains"], "contains_target")

  # Nothing to flag without excluded subclades, or with single-tip subclades
  none <- list(Taxa_name = c("x", "y"), Target_clade = list(toy$A, toy$A), Excluded_subclades = list(NULL, list(toy$s)))
  expect_equal(nrow(check_excluded_subclades_position(toy$phy, none)), 0)
})

## Test flagging of taxa whose constraints left no valid branch for grafting
test_that("find_taxa_without_allowed_branch() lists the taxa whose branches are all excluded",
{
  toy <- toy_setup()
  missing_taxa_list <- list(
    Taxa_name = c("bad", "tip", "fine"),
    Target_clade = list(toy$B, toy$s, toy$A),
    Excluded_subclades = list(list(toy$B), NULL, list(toy$B)))

  # Without the stem, the whole crown of B is excluded: no branch for 'bad'
  flagged <- find_taxa_without_allowed_branch(toy$phy, missing_taxa_list, include_stem = FALSE)
  expect_identical(flagged$Taxa_name, "bad")

  # A single tip always has its terminal branch, and A has branches outside B
  expect_false(any(c("tip", "fine") %in% flagged$Taxa_name))

  # Ignoring the excluded subclades: nothing is flagged
  expect_equal(nrow(find_taxa_without_allowed_branch(toy$phy, missing_taxa_list, include_stem = FALSE, use_excluded_subclades = FALSE)), 0)
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
