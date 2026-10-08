
##### Test suite to ensure branch segments provide similar resutls than whole branches #####

# Previously, the branches available for a taxon were handled as whole branches (a set of tipward node IDs).
# To be able to restrict grafting to a part of a branch (e.g., the section of the branch mapped in an allowed state),
# the available branches are now handled as SEGMENTS: portions of branches described by their tipward node and by two distances
# measured from this tipward node, up the branch (0 = tipward node, branch length = rootward node).

# The following regression tests compare results obtained with the previous implementation, with results obtained with the new one, for a given seed.

# The reference values were computed with the code that handled available branches as whole branches (before the introduction of branch segments).
# For a given seed, grafting must keep returning exactly the same trees (same branches and same positions), otherwise results already published
# with a seed would no longer be reproducible. If one of these tests fails, the change altered the random draws.


## Build toy missing_taxa lists for the tests
regression_taxa <- function (toy)
{
  list(mtl1 = list(Taxon_name = paste0("X", 1:7),
                   Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$s, X4 = toy$A, X5 = toy$B, X6 = toy$C, X7 = toy$s),
                   include_target_clade_stem = c(X1 = TRUE, X2 = FALSE, X3 = NA, X4 = FALSE, X5 = TRUE, X6 = NA, X7 = NA),
                   Excluded_subclades = list(X1 = list(toy$B), X2 = list(), X3 = list(), X4 = list(toy$B, toy$C), X5 = list(), X6 = list(), X7 = list()),
                   exclude_subclades_stem = list(X1 = FALSE, X2 = NA, X3 = NA, X4 = c(TRUE, FALSE), X5 = NA, X6 = NA, X7 = NA)),
       mtl2 = list(Taxon_name = paste0("Y", 1:12), Target_clade = stats::setNames(rep(list(toy$A, toy$D, toy$s), 4), paste0("Y", 1:12))))
}

## Function to compare grafted trees based on total edge length
expect_same_tree_signature <- function (phy, sum_edge, sum_edge_squared, n_tips)
{
  expect_equal(sum(phy$edge.length), sum_edge, tolerance = 1e-7)
  expect_equal(sum(phy$edge.length^2), sum_edge_squared, tolerance = 1e-7)
  expect_equal(ape::Ntip(phy), n_tips)
}


test_that("method = 'uniform' keeps returning the same trees for a given seed",
{
  toy <- toy_setup()
  taxa <- regression_taxa(toy)
  expect_same_tree_signature(quiet_graft(graft_missing_taxa, toy$phy, taxa$mtl1, method = "uniform", seed = 1), 192.641517, 3663.388406, 27)
  expect_same_tree_signature(quiet_graft(graft_missing_taxa, toy$phy, taxa$mtl1, method = "uniform", seed = 2), 176.138253, 3556.124778, 27)
  expect_same_tree_signature(quiet_graft(graft_missing_taxa, toy$phy, taxa$mtl2, method = "uniform", seed = 1), 331.344696, 5359.291248, 32)
})

test_that("method = 'BD_clades' keeps returning the same trees for a given seed",
{
  toy <- toy_setup()
  taxa <- regression_taxa(toy)
  expect_same_tree_signature(quiet_graft(graft_missing_taxa, toy$phy, taxa$mtl1, method = "BD_clades", seed = 1), 178.205603, 3547.972578, 27)
  expect_same_tree_signature(quiet_graft(graft_missing_taxa, toy$phy, taxa$mtl1, method = "BD_clades", seed = 2), 174.763661, 3538.858991, 27)
  expect_same_tree_signature(quiet_graft(graft_missing_taxa, toy$phy, taxa$mtl2, method = "BD_clades", seed = 1), 175.791329, 3545.542197, 32)
})
