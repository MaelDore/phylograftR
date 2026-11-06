
##### Test suite for the names of the clades of total_richness_per_clades #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Backbone with 3 taxa: X1 and X2 in the clade B, X3 in the clade A (B is nested in A)
## Insert the names of the clades after the clades, before the total richness (as documented)
with_names <- function (totals, names, item = "Target_clade")
{
  out <- totals[item]
  out[[paste0(item, "_name")]] <- names
  c(out, totals[names(totals) != item])
}

names_setup <- function ()
{
  toy <- toy_setup()
  list(toy = toy,
       taxa = c("X1", "X2", "X3"),
       clades = list(X1 = toy$B, X2 = toy$B, X3 = toy$A),
       totals = list(Target_clade = list(toy$B, toy$A), total_richness = c(toy$nB + 5, toy$nA + 8)))
}

test_that("check_clade_names() validates the names and their consistency with the missing taxa list", 
{
  s <- names_setup()
  toy <- s$toy
  node_B <- ape::getMRCA(toy$phy, toy$B)
  node_A <- ape::getMRCA(toy$phy, toy$A)
  check <- function (names, missing_names = NULL, missing_nodes = c(node_B, node_B, node_A))
  {
    check_clade_names(names, clade_nodes = c(node_B, node_A), clade_labels = c("B", "A"), item = "tr_name", n_clades = 2,
                      missing_nodes = missing_nodes, missing_clade_names = missing_names, missing_item = "list_name")
  }
  expect_null(check(NULL))
  expect_no_error(check(c("B", "A")))
  expect_no_error(check(c("B", NA), c("B", "B", "A")))
  # Names of the list only (not in total_richness_per_clades): nothing to compare
  expect_no_error(check(c(NA, NA), c("B", "B", "A")))
  # Same name in both objects
  expect_no_error(check(c("B", "A"), c("B", "B", "A")))
  # Taxa without name in the list
  expect_no_error(check(c("B", "A"), c(NA, "B", NA)))

  expect_error(check(c("B", "A", "C")), "one name per clade")
  expect_error(check(1:2), "one name per clade")
  expect_error(check(c("B", "B")), "must be unique")
  # Different names for the same clade in both objects
  expect_error(check(c("B", "A"), c("B", "B", "other")), "different name")
  # Taxa sharing the same clade with different names in the list
  expect_error(check(c("B", "A"), c("B", "other", "A")), "share the same clade")
  # A name identifying different clades in both objects
  expect_error(check(c("B", NA), c("A", "A", "A"), missing_nodes = c(node_B, node_B, node_A)), "different name|identify different clades")
  expect_error(check(c("B", NA), c(NA, NA, "B")), "identify different clades")
})

test_that("graft_missing_taxa() checks the names of total_richness_per_clades against Target_clade_name", 
{
  s <- names_setup()
  mtl <- list(Taxon_name = s$taxa, Target_clade = s$clades, Target_clade_name = c("clade_B", "clade_B", "clade_A"))
  tr <- with_names(s$totals, c("clade_B", "clade_A"))
  g <- quiet_graft(graft_missing_taxa, s$toy$phy, mtl, total_richness_per_clades = tr, seed = 1)
  expect_s3_class(g, "phylo")
  # Names are not required in both objects
  expect_no_error(quiet_graft(graft_missing_taxa, s$toy$phy, mtl, total_richness_per_clades = s$totals, seed = 1))
  mtl_no_names <- mtl; mtl_no_names$Target_clade_name <- NULL
  expect_no_error(quiet_graft(graft_missing_taxa, s$toy$phy, mtl_no_names, total_richness_per_clades = tr, seed = 1))
  # The names do not change the grafting
  g_no_names <- quiet_graft(graft_missing_taxa, s$toy$phy, mtl_no_names, total_richness_per_clades = s$totals, seed = 1)
  expect_equal(g$edge.length, g_no_names$edge.length)

  # A different name for the same clade, or a name identifying two clades
  tr_wrong <- with_names(s$totals, c("clade_A", "clade_B"))
  expect_error(quiet_graft(graft_missing_taxa, s$toy$phy, mtl, total_richness_per_clades = tr_wrong, seed = 1), "different name")
  tr_dup <- with_names(s$totals, c("same", "same"))
  expect_error(quiet_graft(graft_missing_taxa, s$toy$phy, mtl, total_richness_per_clades = tr_dup, seed = 1), "must be unique")
  tr_len <- with_names(s$totals, "clade_B")
  expect_error(quiet_graft(graft_missing_taxa, s$toy$phy, mtl, total_richness_per_clades = tr_len, seed = 1), "one name per clade")
  # Taxa of the same clade with different names in the list
  mtl_conflict <- mtl; mtl_conflict$Target_clade_name <- c("clade_B", "other", "clade_A")
  expect_error(quiet_graft(graft_missing_taxa, s$toy$phy, mtl_conflict, total_richness_per_clades = tr, seed = 1), "share the same clade")
})

test_that("graft_missing_clades() checks the names of total_richness_per_clades against Sister_clade_name", 
{
  toy <- toy_setup()
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$C), Sister_clade_name = c("sis_B", "sis_C"))
  tr <- list(Sister_clade = list(toy$B, toy$C), Sister_clade_name = c("sis_B", "sis_C"), total_richness = c(toy$nB + 3, NA))
  expect_s3_class(quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = tr, seed = 1), "phylo")
  # Names in the list only, or in total_richness_per_clades only
  tr_no_names <- tr; tr_no_names$Sister_clade_name <- NULL
  expect_no_error(quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = tr_no_names, seed = 1))
  mcl_no_names <- mcl; mcl_no_names$Sister_clade_name <- NULL
  expect_no_error(quiet_graft(graft_missing_clades, toy$phy, mcl_no_names, total_richness_per_clades = tr, seed = 1))

  tr_wrong <- tr; tr_wrong$Sister_clade_name <- c("sis_B", "other")
  expect_error(quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = tr_wrong, seed = 1), "different name")
  tr_swapped <- tr; tr_swapped$Sister_clade_name <- c("sis_C", "sis_B")
  expect_error(quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = tr_swapped, seed = 1), "different name")
  # The names are optional with NA
  tr_na <- tr; tr_na$Sister_clade_name <- c(NA, "sis_C")
  expect_no_error(quiet_graft(graft_missing_clades, toy$phy, mcl, total_richness_per_clades = tr_na, seed = 1))
})

test_that("convert_nodes_to_tips() keeps the names of the clades of total_richness_per_clades", 
{
  s <- names_setup()
  node_B <- ape::getMRCA(s$toy$phy, s$toy$B)
  node_A <- ape::getMRCA(s$toy$phy, s$toy$A)
  mtl <- list(Taxon_name = s$taxa, Target_clade = list(X1 = node_B, X2 = node_B, X3 = node_A), Target_clade_name = c("clade_B", "clade_B", "clade_A"))
  tr <- list(Target_clade = list(node_B, node_A), Target_clade_name = c("clade_B", "clade_A"), total_richness = c(s$toy$nB + 5, s$toy$nA + 8))
  converted <- suppressMessages(convert_nodes_to_tips(mtl, s$toy$phy, total_richness_per_clades = tr, verbose = FALSE))
  expect_equal(converted$total_richness_per_clades$Target_clade_name, c("clade_B", "clade_A"))
  expect_equal(converted$total_richness_per_clades$total_richness, tr$total_richness)
  expect_equal(converted$missing_taxa_list$Target_clade_name, mtl$Target_clade_name)
  # The converted objects pass the checks, and the names are reported
  g <- quiet_graft(graft_missing_taxa, s$toy$phy, converted$missing_taxa_list, total_richness_per_clades = converted$total_richness_per_clades, seed = 1)
  info <- extract_sampling_info(g, converted$total_richness_per_clades)
  expect_equal(info$Clade_name, c("clade_B", "clade_A"))
})
