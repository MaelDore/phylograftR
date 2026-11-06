
##### Test suite for extract_sampling_info() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## First child of the MRCA of a clade (the branch leads to a clade nested within it)
sampling_child <- function (phy, tips) { phy$edge[phy$edge[, 1] == ape::getMRCA(phy, tips), 2][1] }

## Toy tree with taxa grafted by hand at known positions:
# X1 inside B, X2 on the stem of B (sister to B, in A), X3 on the terminal branch of the tip s, X4 on the terminal branch of X3
sampling_toy_graft <- function (toy)
{
  g <- graft_tip_on_branch(toy$phy, "X1", sampling_child(toy$phy, toy$B), frac = 0.5)
  g <- graft_tip_on_branch(g, "X2", ape::getMRCA(g, toy$B), frac = 0.5)
  g <- graft_tip_on_branch(g, "X3", match(toy$s, g$tip.label), frac = 0.5)
  g <- graft_tip_on_branch(g, "X4", match("X3", g$tip.label), frac = 0.5)
  identify_grafted(g, c("X1", "X2", "X3", "X4"))
}


test_that("extract_sampling_info() counts the backbone and grafted species of the clades in the grafted phylogeny",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)
  totals <- list(Target_clade = list(toy$B, toy$A, toy$C), total_richness = c(toy$nB + 10, toy$nA + 12, NA))
  info <- extract_sampling_info(g, totals)

  expect_s3_class(info, "data.frame")
  expect_equal(names(info), c("map", "simulation", "seed", "Clade_tips", "total_richness", "total_richness_provided", "nb_backbone_species",
                              "nb_grafted_species", "nb_unsampled_species", "sampling_fraction"))
  expect_equal(nrow(info), 3)
  expect_equal(attr(info, "type"), "taxa")
  expect_true(all(is.na(info$map)))
  expect_equal(info$simulation, rep(1, 3))
  expect_equal(info$Clade_tips, vapply(totals$Target_clade, paste, character(1), collapse = ", "))
  expect_equal(info$total_richness_provided, c(TRUE, TRUE, FALSE))
  expect_equal(info$nb_backbone_species, c(toy$nB, toy$nA, toy$nC))

  # B: X1 inside, X2 on the stem (outside B). A: X1, X2 (sister to B, nested in A), X3 and X4 (grafted on the tip s of A). C: nothing
  expect_equal(info$nb_grafted_species, c(1, 4, 0))
  expect_equal(info$nb_unsampled_species, c(9, 8, 0))
  expect_equal(info$total_richness, c(toy$nB + 10, toy$nA + 12, toy$nC))
  expect_equal(info$sampling_fraction, c((toy$nB + 1) / (toy$nB + 10), (toy$nA + 4) / (toy$nA + 12), 1))
})

test_that("extract_sampling_info() derives the unsampled species of the clades without total richness from the nested clades",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)
  # A unknown: sum of the unsampled species of the largest clades with a total richness nested within it (B and C)
  totals <- list(Target_clade = list(toy$A, toy$B, toy$C), total_richness = c(NA, toy$nB + 10, toy$nC + 3))
  info <- extract_sampling_info(g, totals)
  expect_equal(info$total_richness_provided, c(FALSE, TRUE, TRUE))
  expect_equal(info$nb_unsampled_species, c(9 + 3, 9, 3))
  expect_equal(info$total_richness[1], toy$nA + 4 + 12)
  expect_equal(info$sampling_fraction[1], (toy$nA + 4) / (toy$nA + 4 + 12))

  # Only NA: nothing is unsampled
  info_na <- extract_sampling_info(g, list(Target_clade = list(toy$A, toy$B), total_richness = c(NA, NA)))
  expect_equal(info_na$nb_unsampled_species, c(0, 0))
  expect_equal(info_na$sampling_fraction, c(1, 1))
  # Unknown clade nested in unknown clade, the unsampled species of the nested one are reported up
  info_nest <- extract_sampling_info(g, list(Target_clade = list(toy$A, toy$B, toy$C), total_richness = c(NA, NA, toy$nC + 3)))
  expect_equal(info_nest$nb_unsampled_species, c(3, 0, 3))
})

test_that("extract_sampling_info() handles the clades of single tips and the sister clades",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)

  # Target clade of a single tip: its terminal branch is its only branch, so the taxa grafted on it (X3, then X4 on X3) are in the clade
  info <- extract_sampling_info(g, list(Target_clade = list(toy$s), total_richness = 4))
  expect_equal(info$nb_backbone_species, 1)
  expect_equal(info$nb_grafted_species, 2)
  expect_equal(info$nb_unsampled_species, 1)

  # Sister clade of a single tip: the missing clades are on its stem, so they are sister to the clade and never in it
  info_sister <- extract_sampling_info(g, list(Sister_clade = list(toy$s), total_richness = 4))
  expect_equal(attr(info_sister, "type"), "clades")
  expect_equal(info_sister$nb_grafted_species, 0)
  expect_equal(info_sister$nb_unsampled_species, 3)
  # Sister clade with two tips: only the taxa within are counted (X1 in B, not X2 on the stem)
  info_sister_B <- extract_sampling_info(g, list(Sister_clade = list(toy$B), total_richness = toy$nB + 5))
  expect_equal(info_sister_B$nb_grafted_species, 1)
})

test_that("extract_sampling_info() reports negative unsampled species and warns",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)
  # B has one more species (X1) than its total richness
  totals <- list(Target_clade = list(toy$B, toy$A, toy$C), total_richness = c(toy$nB, toy$nA + 12, NA))
  expect_warning(info <- extract_sampling_info(g, totals), "negative number of unsampled species")
  expect_equal(info$nb_unsampled_species[1], -1)
  expect_equal(info$sampling_fraction[1], (toy$nB + 1) / toy$nB)
  expect_gt(info$sampling_fraction[1], 1)
  expect_equal(info$nb_unsampled_species[2:3], c(8, 0))

  # Negative values of nested clades count as 0 for the clades derived from them
  totals_na <- list(Target_clade = list(toy$A, toy$B, toy$C), total_richness = c(NA, toy$nB, toy$nC + 3))
  expect_warning(info_na <- extract_sampling_info(g, totals_na), "negative")
  expect_equal(info_na$nb_unsampled_species, c(0 + 3, -1, 3))
})

test_that("extract_sampling_info() gives one table over simulations, maps and lists of maps",
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$A))
  totals <- list(Target_clade = list(toy$A, toy$B), total_richness = c(toy$nA + 20, toy$nB + 10))

  # Several simulations
  g_sim <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", nb_simulations = 3, seed = 5)
  info <- extract_sampling_info(g_sim, totals)
  expect_equal(nrow(info), 6)
  expect_equal(info$simulation, rep(1:3, each = 2))
  expect_true(all(is.na(info$map)))
  expect_false(anyNA(info$seed))
  # Each row is consistent with the tree of its simulation
  for (k in 1:3)
  {
    sub <- info[info$simulation == k, ]
    expect_equal(sub$nb_backbone_species + sub$nb_grafted_species, c(length(ape::extract.clade(g_sim[[k]], ape::getMRCA(g_sim[[k]], toy$A))$tip.label),
                                                                      length(ape::extract.clade(g_sim[[k]], ape::getMRCA(g_sim[[k]], toy$B))$tip.label)))
  }

  # List of maps (multiSimmap)
  maps <- structure(list(m1 = make_toy_simmap(toy$phy, seed = 1), m2 = make_toy_simmap(toy$phy, seed = 2)), class = c("multiSimmap", "multiPhylo"))
  g_maps <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", simmap = maps, nb_simulations = 2, seed = 5)
  info_maps <- extract_sampling_info(g_maps, totals)
  expect_equal(nrow(info_maps), 2 * 2 * 2)
  expect_equal(unique(info_maps$map), c("m1", "m2"))
  expect_equal(info_maps$simulation[info_maps$map == "m2"], rep(1:2, each = 2))
})

test_that("extract_sampling_info() checks its inputs",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)
  ok <- list(Target_clade = list(toy$B, toy$A), total_richness = c(toy$nB + 10, toy$nA + 12))

  expect_error(extract_sampling_info(g, 1), "must be a list")
  expect_error(extract_sampling_info(g, list(total_richness = 1)), "'Target_clade' item")
  expect_error(extract_sampling_info(g, list(Target_clade = list(toy$B), Sister_clade = list(toy$B), total_richness = 10)), "both")
  expect_error(extract_sampling_info(g, list(Target_clade = list(toy$B))), "'total_richness'")
  expect_error(extract_sampling_info(g, list(Target_clade = list(toy$B, toy$A), total_richness = 5)), "same length")
  expect_error(extract_sampling_info(g, list(Target_clade = list(toy$B), total_richness = 5.5)), "whole numbers")
  expect_error(extract_sampling_info(g, list(Target_clade = list(toy$B), total_richness = -1)), "whole numbers")
  expect_error(extract_sampling_info(g, list(Target_clade = list(32, toy$A), total_richness = c(10, 10))), "node IDs")
  expect_error(extract_sampling_info(g, list(Target_clade = list(c("zz", "yy")), total_richness = 10)), "not tips of the backbone")
  # Grafted taxa are not tips of the backbone
  expect_error(extract_sampling_info(g, list(Target_clade = list(c("X1", toy$B[1])), total_richness = 10)), "not tips of the backbone")
  # The same clade twice
  expect_error(extract_sampling_info(g, list(Target_clade = list(toy$B, rev(toy$B)), total_richness = c(10, 10))), "more than once")

  # Phylogenies without the grafted elements
  no_info <- g; no_info$grafted_tips <- NULL
  expect_error(extract_sampling_info(no_info, ok), "identify_grafted")
  expect_error(extract_sampling_info("not a phylogeny", ok))
  expect_no_error(extract_sampling_info(g, ok))
})

test_that("extract_sampling_info() reports the names of the clades when provided",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)
  totals <- list(Target_clade_name = c("clade_B", NA, "clade_C"), Target_clade = list(toy$B, toy$A, toy$C), total_richness = c(toy$nB + 10, toy$nA + 12, NA))
  info <- extract_sampling_info(g, totals)
  expect_equal(names(info), c("map", "simulation", "seed", "Clade_name", "Clade_tips", "total_richness", "total_richness_provided", "nb_backbone_species",
                              "nb_grafted_species", "nb_unsampled_species", "sampling_fraction"))
  expect_equal(info$Clade_name, c("clade_B", NA, "clade_C"))
  # Same counts as without names
  totals$Target_clade_name <- NULL
  expect_equal(info[, names(info) != "Clade_name"], extract_sampling_info(g, totals), ignore_attr = TRUE)

  # Several phylogenies: the names are repeated for each phylogeny
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$B))
  g_sim <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", nb_simulations = 3, seed = 5)
  totals_sim <- list(Target_clade_name = c("A", "B"), Target_clade = list(toy$A, toy$B), total_richness = c(toy$nA + 5, toy$nB + 5))
  info_sim <- extract_sampling_info(g_sim, totals_sim)
  expect_equal(info_sim$Clade_name, rep(c("A", "B"), 3))

  # Sister clades
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$C))
  g_clades <- quiet_graft(graft_missing_clades, toy$phy, mcl, nb_simulations = 3, seed = 5)
  totals_sis <- list(Sister_clade_name = c("sis_B", "sis_C"), Sister_clade = list(toy$B, toy$C), total_richness = c(toy$nB + 3, NA))
  info_sis <- extract_sampling_info(g_clades, totals_sis)
  expect_equal(info_sis$Clade_name, rep(c("sis_B", "sis_C"), 3))
  # Missing clades are always grafted on the stem of their sister clade: the counts are the same in all the phylogenies
  cols <- c("nb_backbone_species", "nb_grafted_species", "nb_unsampled_species", "sampling_fraction")
  expect_equal(info_sis[info_sis$simulation == 2, cols], info_sis[info_sis$simulation == 1, cols], ignore_attr = TRUE)
  expect_equal(info_sis[info_sis$simulation == 3, cols], info_sis[info_sis$simulation == 1, cols], ignore_attr = TRUE)

  # Negative unsampled species: the warning names the clade
  totals_neg <- list(Target_clade_name = c("clade_B", "clade_A", "clade_C"), Target_clade = list(toy$B, toy$A, toy$C), total_richness = c(toy$nB, toy$nA + 12, NA))
  expect_warning(extract_sampling_info(g, totals_neg), "clade_B")
})

test_that("extract_sampling_info() checks the names of the clades",
{
  toy <- toy_setup()
  g <- sampling_toy_graft(toy)
  base <- list(Target_clade = list(toy$B, toy$A), total_richness = c(toy$nB + 10, toy$nA + 12))
  expect_error(extract_sampling_info(g, c(list(Target_clade_name = "B"), base)), "one name per clade")
  expect_error(extract_sampling_info(g, c(list(Target_clade_name = c(1, 2)), base)), "one name per clade")
  expect_error(extract_sampling_info(g, c(list(Target_clade_name = c("same", "same")), base)), "must be unique")
  # NA and empty names are unnamed clades, and are not duplicates
  info <- extract_sampling_info(g, c(list(Target_clade_name = c(NA, "")), base))
  expect_equal(info$Clade_name, c(NA_character_, NA_character_))
  expect_no_error(extract_sampling_info(g, c(list(Target_clade_name = c(NA, NA)), base)))
})

test_that("extract_sampling_info() shows the variation of the richness of nested clades among simulations of graft_missing_taxa()",
{
  toy <- toy_setup()
  # Taxa assigned to A can be grafted within B (nested in A), in A outside B, or on the stem of A
  mtl <- list(Taxon_name = paste0("X", 1:4), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$A, X4 = toy$A))
  totals <- list(Target_clade_name = c("A", "B"), Target_clade = list(toy$A, toy$B), total_richness = c(toy$nA + 10, toy$nB + 6))
  g <- quiet_graft(graft_missing_taxa, toy$phy, mtl, method = "uniform", nb_simulations = 20, seed = 3)
  info <- extract_sampling_info(g, totals)
  expect_gt(length(unique(info$nb_grafted_species[info$Clade_name == "B"])), 1)
  expect_gt(length(unique(info$sampling_fraction[info$Clade_name == "B"])), 1)
  # The backbone species do not vary
  expect_equal(length(unique(info$nb_backbone_species[info$Clade_name == "A"])), 1)
})

