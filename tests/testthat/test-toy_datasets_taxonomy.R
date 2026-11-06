
##### Test suite for the toy taxonomic tables #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


test_that("the toy taxonomic tables are consistent with the toy backbone", {
  data(toy_backbone, package = "phylograftR")
  data(toy_backbone_taxonomic_table, package = "phylograftR")
  data(toy_missing_taxa_taxonomic_table, package = "phylograftR")
  expect_s3_class(toy_backbone_taxonomic_table, "data.frame")
  expect_equal(names(toy_backbone_taxonomic_table), c("Taxa", "Family", "Subfamily", "Tribe", "Genus"))
  expect_equal(toy_backbone_taxonomic_table$Taxa, toy_backbone$tip.label)
  expect_equal(nrow(toy_missing_taxa_taxonomic_table), 10)
  expect_equal(names(toy_missing_taxa_taxonomic_table),
               c("Taxa", "Type", "Family", "Subfamily", "Tribe", "Genus", "Target_clade", "Sister_clade", "Excluded_subclades", "include_target_clade_stem", "exclude_subclades_stem"))
  expect_equal(table(toy_missing_taxa_taxonomic_table$Type)[c("clades", "taxa")], table(c(rep("taxa", 7), rep("clades", 3)))[c("clades", "taxa")])
  # The taxa to graft are absent from the backbone
  expect_false(any(toy_missing_taxa_taxonomic_table$Taxa %in% toy_backbone$tip.label))
  # Some tips are unannotated at the finest levels, but not at the broadest one
  expect_false(anyNA(toy_backbone_taxonomic_table$Family))
  expect_true(anyNA(toy_backbone_taxonomic_table$Genus))
})

test_that("the csv files are the datasets", {
  for (name in c("toy_backbone_taxonomic_table", "toy_missing_taxa_taxonomic_table"))
  {
    file <- system.file("extdata", paste0(name, ".csv"), package = "phylograftR")
    expect_true(nzchar(file))
    tab <- utils::read.csv(file, stringsAsFactors = FALSE, na.strings = "")
    expect_equal(tab, get(name, envir = asNamespace("phylograftR")) , ignore_attr = TRUE)
  }
})

test_that("the toy taxonomic tables build the toy missing taxa list and a list of missing clades without warning", {
  data(toy_backbone, package = "phylograftR")
  data(toy_missing_taxa_list, package = "phylograftR")
  data(toy_backbone_taxonomic_table, package = "phylograftR")
  data(toy_missing_taxa_taxonomic_table, package = "phylograftR")
  expect_no_warning(lists <- build_missing_taxa_list_from_taxonomy(toy_backbone_taxonomic_table, toy_missing_taxa_taxonomic_table, toy_backbone, verbose = FALSE))
  expect_equal(names(lists), c("taxa", "clades"))
  info <- attr(lists, "taxonomic_groups")
  expect_true(all(info$monophyletic))

  # Same clades, excluded subclades and stems as in toy_missing_taxa_list
  node_of <- function (tips) if (length(tips) == 1) match(tips, toy_backbone$tip.label) else ape::getMRCA(toy_backbone, tips)
  nodes <- function (x) unname(vapply(x, node_of, numeric(1)))
  taxa <- lists$taxa
  expect_equal(taxa$Taxon_name, toy_missing_taxa_list$Taxon_name)
  expect_equal(nodes(taxa$Target_clade), nodes(toy_missing_taxa_list$Target_clade))
  expect_equal(lapply(taxa$Excluded_subclades, nodes), lapply(toy_missing_taxa_list$Excluded_subclades, nodes), ignore_attr = TRUE)
  expect_equal(lapply(taxa$exclude_subclades_stem, unname), lapply(toy_missing_taxa_list$exclude_subclades_stem, unname), ignore_attr = TRUE)
  expect_equal(unname(taxa$include_target_clade_stem), unname(toy_missing_taxa_list$include_target_clade_stem))

  # The two taxa of the new genus are grafted together: same sister clade
  clades <- lists$clades
  expect_equal(clades$Taxon_name, c("Genus_d_sp1", "Genus_d_sp2", "Genus_e_sp1"))
  expect_equal(clades$Sister_clade$Genus_d_sp1, clades$Sister_clade$Genus_d_sp2)
  expect_equal(nodes(clades$Sister_clade), c(50, 50, 40))
  expect_equal(unname(clades$Sister_clade_name), c("Gen_1", "Gen_1", "Tribe_A1b"))
})

test_that("the lists built from the toy taxonomic tables can be grafted", {
  data(toy_backbone, package = "phylograftR")
  data(toy_backbone_taxonomic_table, package = "phylograftR")
  data(toy_missing_taxa_taxonomic_table, package = "phylograftR")
  lists <- build_missing_taxa_list_from_taxonomy(toy_backbone_taxonomic_table, toy_missing_taxa_taxonomic_table, toy_backbone, verbose = FALSE)
  g_taxa <- quiet_graft(graft_missing_taxa, toy_backbone, lists$taxa, method = "uniform", seed = 1)
  expect_valid_graft(g_taxa, toy_backbone, lists$taxa)
  g_clades <- quiet_graft(graft_missing_clades, toy_backbone, lists$clades, seed = 1)
  expect_true(all(lists$clades$Taxon_name %in% g_clades$tip.label))
  expect_true(ape::is.ultrametric(g_clades))
})
