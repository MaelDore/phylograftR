# Target clades defined by many tips: the keys built from the tip names exceed 10,000 bytes
# (an environment cannot store such keys: "variable names are limited to 10000 bytes")

test_that("Large target clades (long keys) can be grafted and checked", {
  set.seed(1)
  phy <- ape::rcoal(400)
  phy$tip.label <- sprintf("a_long_species_name_for_the_tip_%04d", seq_along(phy$tip.label))
  phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
  all_tips <- phy$tip.label
  expect_gt(nchar(paste(sort(all_tips), collapse = "|")), 10000)

  missing_taxa_list <- list(Taxon_name = c("New_1", "New_2", "New_3"),
                            Target_clade = stats::setNames(rep(list(all_tips), 3), c("New_1", "New_2", "New_3")))

  for (method in c("uniform", "BD_clades"))
  {
    grafted <- suppressWarnings(graft_missing_taxa(phy, missing_taxa_list, method = method, seed = 1, verbose = FALSE))
    expect_s3_class(grafted, "phylo")
    expect_true(all(missing_taxa_list$Taxon_name %in% grafted$tip.label))
    expect_equal(ape::Ntip(grafted), 403)
  }

  # The grafting plan uses the same kind of keys
  expect_no_error(suppressWarnings(check_grafting_plan(phy, missing_taxa_list, verbose = FALSE)))
})
