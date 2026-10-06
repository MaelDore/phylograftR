
##### Test suite for identify_grafted() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Detect grafted edges for reference: a branch is grafted if all the tips below it are grafted taxa
reference_grafted_edges <- function (tree, grafted)
{
  grafted_ids <- which(tree$tip.label %in% grafted)
  vapply(tree$edge[, 2], function (node)
  {
    desc <- phytools::getDescendants(tree, node)
    desc <- desc[desc <= ape::Ntip(tree)]
    all(desc %in% grafted_ids)
  }, logical(1))
}

## Build toy phylogeny with taxa grafted in different location: taxa on a clade, on its stem, and on a single tip
toy_grafted <- function ()
{
  toy <- toy_setup()
  missing_taxa_list <- list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
                            Target_clade = list(X1 = toy$B, X2 = toy$B, X3 = toy$C, X4 = toy$s, X5 = toy$s))
  grafted <- quiet_graft(graft_missing_taxa, toy$phy, missing_taxa_list, method = "uniform", seed = 3)
  list(toy = toy, grafted = grafted, taxa = missing_taxa_list$Taxon_name)
}

## Test that the output phylogeny contains valid grafted elements
test_that("identify_grafted() adds the grafted tips and branches to the phylogeny",
{
  g <- toy_grafted()
  out <- identify_grafted(g$grafted, g$taxa)

  expect_s3_class(out, "phylo")
  # Tips: named by label, in the order of the tip labels
  expect_identical(names(out$grafted_tips), out$tip.label)
  expect_setequal(out$tip.label[out$grafted_tips], g$taxa)
  # Branches: named by tipward node, in the order of the edges, and equal to the reference
  expect_identical(names(out$grafted_edges), as.character(out$edge[, 2]))
  expect_identical(unname(out$grafted_edges), reference_grafted_edges(out, g$taxa))
  # The terminal branch of each grafted tip is a grafted branch
  expect_true(all(out$grafted_edges[as.character(match(g$taxa, out$tip.label))]))
  # The other elements of the phylogeny are untouched, and the input is not modified
  expect_identical(out$edge, g$grafted$edge)
  expect_identical(out$edge.length, g$grafted$edge.length)
  expect_identical(out$tip.label, g$grafted$tip.label)
})

## Test that branches tagged as grafted subtend missing-only subtrees
test_that("grafted branches include the branches of missing-only subtrees (grafted on grafted taxa)",
{
  toy <- toy_setup()
  single <- match(toy$s, toy$phy$tip.label)
  g1 <- graft_tip_on_branch(toy$phy, "X1", single, 0.5)
  g2 <- graft_tip_on_branch(g1, "X2", match("X1", g1$tip.label), 0.5)
  out <- identify_grafted(g2, c("X1", "X2"))
  expect_identical(unname(out$grafted_edges), reference_grafted_edges(g2, c("X1", "X2")))
  # The cherry (X1, X2) and its stem branch are grafted: 2 tips + 1 internal branch
  expect_equal(sum(out$grafted_edges), 3)
})

## Test update behavior with ladderize() and drop.tip()
test_that("identify_grafted() elements can be realigned after ladderize(), and are recomputed after drop.tip()",
{
  g <- toy_grafted()
  out <- identify_grafted(g$grafted, g$taxa)

  # ladderize() reorders the edges but node IDs are kept: names allow realigning the edge vector
  lad <- ape::ladderize(out)
  expect_false(identical(lad$edge, out$edge))
  expect_identical(unname(lad$grafted_edges[as.character(lad$edge[, 2])]), reference_grafted_edges(lad, g$taxa))
  expect_identical(unname(identify_grafted(lad, g$taxa)$grafted_edges), reference_grafted_edges(lad, g$taxa))

  # drop.tip() does not update the elements: they are stale, and identify_grafted() updates them
  pruned <- ape::drop.tip(out, c(g$taxa[1], g$toy$tipsA[1]))
  expect_false(length(pruned$grafted_tips) == ape::Ntip(pruned))
  remaining <- intersect(g$taxa, pruned$tip.label)
  pruned <- identify_grafted(pruned, remaining)
  expect_length(pruned$grafted_tips, ape::Ntip(pruned))
  expect_length(pruned$grafted_edges, nrow(pruned$edge))
  expect_identical(unname(pruned$grafted_edges), reference_grafted_edges(pruned, remaining))
  # One element each, replaced and not duplicated
  expect_equal(sum(names(pruned) == "grafted_tips"), 1)
})

## Test validity checks for the input
test_that("identify_grafted() handles edge cases and invalid inputs",
{
  phy <- make_toy_phy(n_tips = 10)

  none <- identify_grafted(phy, character(0))
  expect_false(any(none$grafted_tips))
  expect_false(any(none$grafted_edges))

  all_grafted <- identify_grafted(phy, phy$tip.label)
  expect_true(all(all_grafted$grafted_tips))
  expect_true(all(all_grafted$grafted_edges))

  # Taxa that are not tips are ignored with a warning
  expect_warning(out <- identify_grafted(phy, c("t1", "absent")), "not tips of 'phy'")
  expect_equal(sum(out$grafted_tips), 1)

  expect_error(identify_grafted(1, "t1"), "class 'phylo'")
  expect_error(identify_grafted(phy, 1), "character vector")
})
