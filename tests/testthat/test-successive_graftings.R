##### Test suite for successive graftings (taxa then clades) and accumulation of the grafted flags #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Reference: a branch is grafted if all the tips below it are grafted taxa
successive_reference_edges <- function (tree, grafted)
{
  grafted_ids <- which(tree$tip.label %in% grafted)
  vapply(tree$edge[, 2], function (node)
  {
    desc <- phytools::getDescendants(tree, node)
    desc <- desc[desc <= ape::Ntip(tree)]
    all(desc %in% grafted_ids)
  }, logical(1))
}

## Lists of taxa and of missing clades for the toy backbone
successive_lists <- function (toy)
{
  list(taxa = list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$B, X2 = toy$B, X3 = toy$C)),
       clades = list(Taxon_name = c("M1", "M2", "M3"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$s)))
}


test_that("the grafted flags accumulate when missing clades are grafted on a phylogeny with grafted taxa",
{
  toy <- toy_setup()
  lists <- successive_lists(toy)
  g1 <- quiet_graft(graft_missing_taxa, toy$phy, lists$taxa, method = "uniform", seed = 1)
  g2 <- quiet_graft(graft_missing_clades, g1, lists$clades, seed = 1)
  all_grafted <- c(lists$taxa$Taxon_name, lists$clades$Taxon_name)

  # Step 1 flags only the taxa, step 2 flags all the grafted taxa
  expect_equal(sort(names(g1$grafted_tips)[g1$grafted_tips]), sort(lists$taxa$Taxon_name))
  expect_equal(sort(names(g2$grafted_tips)[g2$grafted_tips]), sort(all_grafted))
  expect_equal(sum(g2$grafted_tips), 6)
  expect_equal(ape::Ntip(g2), ape::Ntip(toy$phy) + 6)
  expect_equal(names(g2$grafted_tips), g2$tip.label)
  # The backbone tips are never flagged
  expect_false(any(g2$grafted_tips[toy$phy$tip.label]))
  # The grafted edges are consistent with the grafted tips
  expect_equal(unname(g2$grafted_edges), unname(successive_reference_edges(g2, all_grafted)))
  expect_equal(names(g2$grafted_edges), as.character(g2$edge[, 2]))
  # Same result as recomputing the flags from the full list of grafted taxa
  expect_equal(identify_grafted(g2, all_grafted)$grafted_edges, g2$grafted_edges)
  # Backbone branch lengths untouched, and the final phylogeny is ultrametric
  expect_true(ape::is.ultrametric(g2))
})

test_that("a grafting on an ungrafted phylogeny only flags its own taxa (no carry-over)",
{
  toy <- toy_setup()
  lists <- successive_lists(toy)
  g <- quiet_graft(graft_missing_clades, toy$phy, lists$clades, seed = 1)
  expect_equal(sort(names(g$grafted_tips)[g$grafted_tips]), sort(lists$clades$Taxon_name))
})

test_that("the order of the successive graftings (clades first, then taxa) also accumulates the flags",
{
  toy <- toy_setup()
  lists <- successive_lists(toy)
  g1 <- quiet_graft(graft_missing_clades, toy$phy, lists$clades, seed = 2)
  g2 <- quiet_graft(graft_missing_taxa, g1, lists$taxa, method = "uniform", seed = 2)
  expect_equal(sort(names(g2$grafted_tips)[g2$grafted_tips]), sort(c(lists$taxa$Taxon_name, lists$clades$Taxon_name)))
})

test_that("successive graftings can be chained per phylogeny with lapply() for several simulations",
{
  toy <- toy_setup()
  lists <- successive_lists(toy)
  g1 <- quiet_graft(graft_missing_taxa, toy$phy, lists$taxa, method = "uniform", nb_simulations = 3, seed = 5)
  g2 <- lapply(g1, function (x) quiet_graft(graft_missing_clades, x, lists$clades, seed = 5))
  expect_length(g2, 3)
  expect_equal(names(g2), names(g1))
  for (x in g2) { expect_equal(sum(x$grafted_tips), 6) }
  # The taxa of the first step are the same in the second step
  for (i in seq_along(g1)) { expect_true(all(lists$taxa$Taxon_name %in% g2[[i]]$tip.label)) }
})

test_that("a phylogeny modified after grafting (grafted_tips out of date) warns and keeps only the tips still present",
{
  toy <- toy_setup()
  lists <- successive_lists(toy)
  g1 <- quiet_graft(graft_missing_taxa, toy$phy, lists$taxa, method = "uniform", seed = 1)
  dropped <- ape::drop.tip(g1, "X1")
  # drop.tip() keeps the elements as they are (out of date)
  expect_warning(g2 <- graft_missing_clades(phy = dropped, missing_taxa_list = lists$clades, seed = 1, verbose = FALSE),
                 "does not match the tips")
  expect_equal(sort(names(g2$grafted_tips)[g2$grafted_tips]), sort(c("X2", "X3", lists$clades$Taxon_name)))
  # No warning once the flags are recomputed
  fixed <- identify_grafted(dropped, c("X2", "X3"))
  expect_no_warning(graft_missing_clades(phy = fixed, missing_taxa_list = lists$clades, seed = 1, verbose = FALSE))
})

test_that("extract_sampling_info() and plot_grafted() use the accumulated flags of the final phylogeny",
{
  toy <- toy_setup()
  lists <- successive_lists(toy)
  g1 <- quiet_graft(graft_missing_taxa, toy$phy, lists$taxa, method = "uniform", seed = 1)
  g2 <- quiet_graft(graft_missing_clades, g1, lists$clades, seed = 1)
  totals <- list(Target_clade = list(toy$A), total_richness = toy$nA + 20)
  info <- extract_sampling_info(g2, totals)
  # Backbone species of the clade A (grafted taxa excluded), and unsampled species: total - backbone - grafted
  expect_equal(info$nb_backbone_species, toy$nA)
  desc <- phytools::getDescendants(g2, ape::getMRCA(g2, toy$A)); desc <- desc[desc <= ape::Ntip(g2)]
  expect_equal(info$nb_grafted_species, sum(g2$grafted_tips[desc]))
  expect_equal(info$nb_unsampled_species, info$total_richness - info$nb_backbone_species - info$nb_grafted_species)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_grafted(g2))
})
