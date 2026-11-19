##### Test suite for build_total_richness_per_clades() and import_total_richness_per_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Lists of taxa and missing clades built from the toy taxonomic tables. Taxonomic groups used below (clades of the toy backbone):
# Family_A (17 tips, 4 missing taxa: 3 within, 1 as sister of Tribe_A1b), Family_B (13 tips), Subfam_B1 (9 tips), Tribe_B2 (3 tips), Gen_2 (1 tip), Gen_1 (3 tips), Tribe_A1b (2 tips)

toy_lists <- function ()
{
  build_missing_taxa_list_from_taxonomy(backbone_taxonomic_table = phylograftR::toy_backbone_taxonomic_table,
                                        missing_taxa_taxonomic_table = phylograftR::toy_missing_taxa_taxonomic_table,
                                        phy = phylograftR::toy_backbone, verbose = FALSE)
}
toy_richness_table <- function () data.frame(Clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"), total_richness = c(25, NA, 14, 8), stringsAsFactors = FALSE)
build_tr <- function (data, lists = toy_lists(), ...) build_total_richness_per_clades(data = data, phy = phylograftR::toy_backbone, missing_taxa_list = lists, verbose = FALSE, ...)


### 1/ Output ####

test_that("one total_richness_per_clades is returned per type of taxa, with all the clades of the table and of the lists",
{
  lists <- toy_lists()
  r <- build_tr(toy_richness_table(), lists)
  expect_named(r, c("taxa", "clades"))
  expect_named(r$taxa, c("Target_clade", "Target_clade_name", "total_richness"))
  expect_named(r$clades, c("Sister_clade", "Sister_clade_name", "total_richness"))
  # Same clades, names and values in both lists
  expect_identical(r$taxa$Target_clade, r$clades$Sister_clade)
  expect_identical(r$taxa$Target_clade_name, r$clades$Sister_clade_name)
  expect_identical(r$taxa$total_richness, r$clades$total_richness)
  expect_type(r$taxa$total_richness, "integer")
  expect_identical(r$taxa$total_richness[1:4], c(25L, NA, 14L, 8L))
  # The 4 clades of the table, then the clades of the lists that are not in the table (NA richness)
  expect_length(r$taxa$Target_clade, 7)
  expect_true(all(is.na(r$taxa$total_richness[5:7])))
  expect_identical(r$taxa$Target_clade_name[1:4], c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"))
  # The tips are the tips of the lists
  expect_identical(r$taxa$Target_clade[[1]], lists$taxa$Target_clade$New_sp1)
  expect_true(all(unlist(lapply(lists$taxa$Target_clade, function (x) any(vapply(r$taxa$Target_clade, identical, logical(1), x))))))
  expect_true(all(unlist(lapply(lists$clades$Sister_clade, function (x) any(vapply(r$clades$Sister_clade, identical, logical(1), x))))))
})

test_that("the output follows the structure of missing_taxa_list: a single list for a single type",
{
  lists <- toy_lists()
  # The taxonomic groups of the lists are kept as an attribute of the object, so they are carried over to the single lists to resolve names
  lists_taxa <- lists$taxa; attr(lists_taxa, "taxonomic_groups") <- attr(lists, "taxonomic_groups")
  lists_clades <- lists$clades; attr(lists_clades, "taxonomic_groups") <- attr(lists, "taxonomic_groups")
  r_taxa <- build_tr(toy_richness_table(), lists_taxa)
  expect_named(r_taxa, c("Target_clade", "Target_clade_name", "total_richness"))
  r_clades <- build_tr(toy_richness_table(), lists_clades)
  expect_named(r_clades, c("Sister_clade", "Sister_clade_name", "total_richness"))
  # Same clades as with both lists if the table covers the clades of the list
  both <- build_tr(toy_richness_table(), lists)
  expect_identical(r_taxa$Target_clade[1:4], both$taxa$Target_clade[1:4])
  expect_length(r_taxa$Target_clade, 4 + 1) # 4 rows of the table + the target clade of the list that is not in the table (Gen_2)
  expect_length(r_clades$Sister_clade, 4 + 2) # + the sister clades of the list that are not in the table (Gen_1 and Tribe_A1b)
})

test_that("the richness table reports the richness and the sampling fraction of the clades before grafting",
{
  r <- build_tr(toy_richness_table())
  tab <- attr(r, "richness_table")
  expect_s3_class(tab, "data.frame")
  expect_named(tab, c("Clade_name", "Clade_tips", "total_richness_provided", "total_richness", "nb_backbone_species", "nb_missing_species", "nb_unsampled_species", "sampling_fraction"))
  expect_equal(nrow(tab), 7)
  fam_a <- tab[tab$Clade_name == "Family_A", ]
  expect_equal(fam_a$nb_backbone_species, 17)
  expect_equal(fam_a$nb_missing_species, 4)   # 3 taxa within Family_A, and the missing clade grafted on the stem of Tribe_A1b
  expect_equal(fam_a$nb_unsampled_species, 4)
  expect_equal(fam_a$sampling_fraction, 21 / 25)
  # Unknown richness derived from nested clades
  expect_true(is.na(tab$total_richness_provided[tab$Clade_name == "Family_B"]))
  expect_false(is.na(tab$total_richness[tab$Clade_name == "Family_B"]))
})

test_that("the richness table is returned by default and only on request",
{
  lists <- toy_lists()
  expect_false(is.null(attr(build_tr(toy_richness_table(), lists), "richness_table")))
  expect_false(is.null(attr(build_tr(toy_richness_table(), lists, return_richness_table = TRUE), "richness_table")))
  r <- build_tr(toy_richness_table(), lists, return_richness_table = FALSE)
  expect_null(attr(r, "richness_table"))
  # The rest of the output is the same
  with_table <- build_tr(toy_richness_table(), lists)
  attr(with_table, "richness_table") <- NULL
  expect_identical(r, with_table)
  # The message points to the table only if it is returned
  expect_message(build_total_richness_per_clades(toy_richness_table(), phylograftR::toy_backbone, lists, verbose = TRUE), "attr\\(, \"richness_table\"\\)")
  msg <- testthat::capture_messages(build_total_richness_per_clades(toy_richness_table(), phylograftR::toy_backbone, lists, return_richness_table = FALSE, verbose = TRUE))
  expect_false(any(grepl("richness_table", msg)))
  expect_error(build_tr(toy_richness_table(), lists, return_richness_table = NA), "'return_richness_table' must be TRUE or FALSE")
  expect_error(build_tr(toy_richness_table(), lists, return_richness_table = "yes"), "'return_richness_table' must be TRUE or FALSE")
})

test_that("the expected sampling fractions of the table can differ from the effective ones extracted after grafting",
{
  lists <- toy_lists()
  r <- build_tr(toy_richness_table(), lists)
  expected <- attr(r, "richness_table")
  g <- graft_missing_taxa(phylograftR::toy_backbone, lists$taxa, total_richness_per_clades = r$taxa, nb_simulations = 5, seed = 1, verbose = FALSE)
  # The effective richness and sampling fraction of the clades are extracted a posteriori, with the same total richness
  effective <- extract_sampling_info(g, total_richness_per_clades = r$taxa)
  expect_s3_class(effective, "data.frame")
  expect_true(all(c("Clade_name", "simulation", "sampling_fraction") %in% names(effective)))
  fam_a <- effective[effective$Clade_name == "Family_A", ]
  expect_equal(nrow(fam_a), 5)
  # The table counts all the taxa of both types in the clade, whereas only the taxa are grafted so far
  expect_false(isTRUE(all.equal(fam_a$sampling_fraction, rep(expected$sampling_fraction[expected$Clade_name == "Family_A"], 5))))
  # The effective sampling fraction of a clade can vary among simulations, as taxa are not grafted at the same place
  varies <- vapply(split(effective$sampling_fraction, effective$Clade_name), function (x) length(unique(round(x, 8))) > 1, logical(1))
  expect_true(any(varies))
})

test_that("the lists can be used in the grafting functions in sequence, with the same richness table",
{
  lists <- toy_lists()
  r <- build_tr(toy_richness_table(), lists)
  g1 <- graft_missing_taxa(phylograftR::toy_backbone, lists$taxa, total_richness_per_clades = r$taxa, seed = 1, verbose = FALSE)
  expect_true(all(lists$taxa$Taxon_name %in% g1$tip.label))
  g2 <- graft_missing_clades(g1, lists$clades, total_richness_per_clades = r$clades, seed = 1, verbose = FALSE)
  expect_true(all(c(lists$taxa$Taxon_name, lists$clades$Taxon_name) %in% g2$tip.label))
  # The other order works too
  g3 <- graft_missing_clades(phylograftR::toy_backbone, lists$clades, total_richness_per_clades = r$clades, seed = 1, verbose = FALSE)
  g4 <- graft_missing_taxa(g3, lists$taxa, total_richness_per_clades = r$taxa, seed = 1, verbose = FALSE)
  expect_true(all(c(lists$taxa$Taxon_name, lists$clades$Taxon_name) %in% g4$tip.label))
})


### 2/ Ways to write the clades ####

test_that("clades can be written with tips, node IDs or names, and give the tips used in the lists",
{
  lists <- toy_lists(); phy <- phylograftR::toy_backbone
  ref <- build_tr(toy_richness_table(), lists)
  node_A <- ape::getMRCA(phy, c("t26", "t19"))
  d <- data.frame(Clade = c(as.character(node_A), "t8,t22", "t8, t7", "t5,t12"), total_richness = c(25, NA, 14, 8))
  r <- build_tr(d, lists)
  expect_identical(r$taxa$Target_clade, ref$taxa$Target_clade)
  expect_identical(r$taxa$total_richness, ref$taxa$total_richness)
  # The node ID is converted to the tips of the lists, and the name comes from the lists
  expect_identical(r$taxa$Target_clade_name[1], "Family_A")
  # A pair of tips of another pair with the same MRCA node is also converted to the tips of the lists
  d2 <- d; d2$Clade[1] <- "t26, t18"
  expect_identical(build_tr(d2, lists)$taxa$Target_clade, ref$taxa$Target_clade)
  # Conversions are printed
  expect_message(build_total_richness_per_clades(d, phy, lists, verbose = TRUE), "node ID\\(s\\) converted into tips")
})

test_that("clade names are resolved with the names of the lists and with the taxonomic groups",
{
  lists <- toy_lists()
  # Remove the names of the clades in the lists: names are then resolved with the taxonomic groups of the taxonomy
  lists_no_names <- lists
  lists_no_names$taxa$Target_clade_name <- NULL
  lists_no_names$clades$Sister_clade_name <- NULL
  lists_no_names$taxa$Excluded_subclades <- lapply(lists_no_names$taxa$Excluded_subclades, function (x) if (is.null(x)) NULL else unname(x))
  expect_s3_class(attr(lists_no_names, "taxonomic_groups"), "data.frame")
  r <- build_tr(toy_richness_table(), lists_no_names)
  ref <- build_tr(toy_richness_table(), lists)
  expect_identical(r$taxa$Target_clade[1:4], ref$taxa$Target_clade[1:4])
  expect_identical(r$taxa$Target_clade_name[1:4], c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"))
  # A group that is not used by any taxon can be named
  d <- data.frame(Clade = c("Family_A", "Gen_3"), total_richness = c(25, 3))
  expect_no_error(build_tr(d, lists_no_names))
  # Without the lists names and without the taxonomic groups, the name is unknown
  attr(lists_no_names, "taxonomic_groups") <- NULL
  expect_error(build_tr(toy_richness_table(), lists_no_names), "'Clade' of row 1 \\('Family_A'\\) is not a tip of the backbone")
  # A name that does not exist
  expect_error(build_tr(data.frame(Clade = "Family_Z", total_richness = 5), lists), "'Clade' of row 1 \\('Family_Z'\\).*name of a clade")
})

test_that("clade names written in a second column are checked",
{
  lists <- toy_lists()
  d <- data.frame(Clade = c("t26,t19", "t8,t22"), Clade_name = c("Family_A", "other_B"), total_richness = c(25, NA))
  # The name of Family_B in the lists is different
  expect_error(build_tr(d, lists), "Inconsistent labels")
  d$Clade_name[2] <- "Family_B"
  r <- build_tr(d, lists)
  expect_identical(r$taxa$Target_clade_name[1:2], c("Family_A", "Family_B"))
  # A clade without name in the table takes the name of the lists
  d$Clade_name <- NA
  expect_identical(build_tr(d, lists)$taxa$Target_clade_name[1:2], c("Family_A", "Family_B"))
  # A name written in the cell and a different 'Clade_name'
  d2 <- data.frame(Clade = "Family_A", Clade_name = "Fam", total_richness = 25)
  expect_error(build_tr(d2, lists), "different 'Clade_name'")
  # A name used for two clades
  d3 <- data.frame(Clade = c("t26,t19", "t8,t22"), Clade_name = c("Family_B", "Family_B"), total_richness = c(25, NA))
  expect_error(build_tr(d3, lists), "Inconsistent labels")
})

test_that("legacy column names are accepted, with the same result",
{
  lists <- toy_lists()
  ref <- build_tr(toy_richness_table(), lists)
  d <- data.frame(Target_clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"), total_richness = c(25, NA, 14, 8))
  expect_identical(build_tr(d, lists), ref)
  d2 <- data.frame(Sister_clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"), Sister_clade_name = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"), total_richness = c(25, NA, 14, 8))
  expect_identical(build_tr(d2, lists), ref)
})


### 3/ Coverage ####

test_that("clades of the lists that are not in the table are added with an unknown richness, and reported",
{
  lists <- toy_lists()
  expect_message(build_total_richness_per_clades(toy_richness_table(), phylograftR::toy_backbone, lists, verbose = TRUE), "7 clade\\(s\\) in total richness: 3 with a value, 4 with an unknown richness \\(including 3 clade\\(s\\)")
  # With a table that covers everything, nothing is added
  all_groups <- attr(lists, "taxonomic_groups")
  d <- data.frame(Clade = c("Family_A", "Family_B", "Subfam_B1", "Gen_2", "Tribe_B2", "Gen_1", "Tribe_A1b"), total_richness = NA)
  expect_message(build_total_richness_per_clades(d, phylograftR::toy_backbone, lists, verbose = TRUE), "7 clade\\(s\\) in total richness: 0 with a value, 7 with an unknown richness\\. ")
  # No value at all: the grafting is the same as without richness (all NA is accepted by the grafting functions)
  r <- build_tr(d, lists)
  expect_true(all(is.na(r$taxa$total_richness)))
})


### 4/ Consistency of the richness with the full sequence of grafting ####

test_that("a richness lower than the species of the clade is rejected, counting the taxa of both types",
{
  lists <- toy_lists()
  # Family_A has 17 backbone species and 4 taxa to graft (3 as taxa, 1 as a missing clade): 21 species at least
  expect_error(build_tr(data.frame(Clade = "Family_A", total_richness = 20), lists), "not consistent with the species of the backbone and the taxa to graft.*lower than the number of backbone species plus missing species")
  expect_no_error(build_tr(data.frame(Clade = "Family_A", total_richness = 21), lists))
  # The grafting function of the taxa alone would accept 20 (17 + 3 taxa), as it does not know the missing clades
  r_taxa_only <- build_tr(data.frame(Clade = "t26,t19", total_richness = 20), lists$taxa)
  expect_equal(attr(r_taxa_only, "richness_table")$nb_missing_species[1], 3)
})

test_that("a clade must have at least as many unsampled species as the clades with a richness nested within it",
{
  lists <- toy_lists()
  # Family_B has 13 backbone species and 6 taxa to graft: with 20 species, it has 1 unsampled species. Subfam_B1 (nested) with 14 species has 2 unsampled species
  expect_error(build_tr(data.frame(Clade = c("Family_B", "Subfam_B1"), total_richness = c(20, 14)), lists), "too low to contain the clades")
  expect_no_error(build_tr(data.frame(Clade = c("Family_B", "Subfam_B1"), total_richness = c(21, 14)), lists))
})


### 5/ Checks of the inputs ####

test_that("build_total_richness_per_clades() checks the inputs",
{
  lists <- toy_lists(); d <- toy_richness_table(); phy <- phylograftR::toy_backbone
  expect_error(build_total_richness_per_clades(d, "phy", lists), "'phy' must be")
  expect_error(build_total_richness_per_clades(as.list(d), phy, lists), "'data' must be a data frame")
  expect_error(build_total_richness_per_clades(d[0, ], phy, lists), "no rows")
  expect_error(build_total_richness_per_clades(d, phy), "'missing_taxa_list' must be provided")
  expect_error(build_total_richness_per_clades(d, phy, "x"), "'missing_taxa_list' must be provided")
  expect_error(build_total_richness_per_clades(d, phy, lists, verbose = NA), "'verbose'")
  expect_error(build_total_richness_per_clades(d, phy, lists, return_richness_table = NULL), "'return_richness_table'")
  expect_error(build_tr(d[, "total_richness", drop = FALSE], lists), "must have a column 'Clade'")
  expect_error(build_tr(d[, "Clade", drop = FALSE], lists), "must have a column 'total_richness'")
  d2 <- d; d2$Target_clade <- d2$Clade
  expect_error(build_tr(d2, lists), "only one of the columns")
  expect_warning(build_tr(cbind(d, typo = 1), lists), "Unknown columns ignored: 'typo'")
  d3 <- d; d3$total_richness[2] <- "ten"
  expect_error(build_tr(d3, lists), "'total_richness' must be a non-negative whole number.*2 \\('ten'\\)")
  d3$total_richness[2] <- "10.5"
  expect_error(build_tr(d3, lists), "whole number")
  d3$total_richness[2] <- "-1"
  expect_error(build_tr(d3, lists), "whole number")
  d4 <- d; d4$Clade[3] <- NA
  expect_error(build_tr(d4, lists), "'Clade' is empty in the row\\(s\\): 3")
  d5 <- d; d5$Clade[2] <- "t8,,t22"
  expect_error(build_tr(d5, lists), "empty tip")
  d5$Clade[2] <- "t8,zzz"
  expect_error(build_tr(d5, lists), "missing from the backbone phylogeny \\(row 2\\): zzz")
  d6 <- rbind(d, data.frame(Clade = "t26, t18", total_richness = 25))
  expect_error(build_tr(d6, lists), "listed more than once.*rows 1, 5")
  d7 <- d; d7$Clade[2] <- "t8,t8"
  expect_error(build_tr(d7, lists), "repeated tips")
})

test_that("the missing_taxa_list must be valid and defined by tips",
{
  d <- toy_richness_table(); lists <- toy_lists()
  l2 <- lists$taxa; l2$Target_clade[[1]] <- 34
  expect_error(build_tr(d, l2), "defined by node IDs")
  both <- lists$taxa; both$Sister_clade <- both$Target_clade
  expect_error(build_tr(d, both), "has both 'Target_clade' and 'Sister_clade'")
  expect_error(build_tr(d, list(a = 1)), "must be a 'missing_taxa_list'")
  dup <- lists; dup$clades$Taxon_name[1] <- lists$taxa$Taxon_name[1]; names(dup$clades$Sister_clade)[1] <- lists$taxa$Taxon_name[1]
  expect_error(build_tr(d, dup), "must be different|unique|Duplicated")
})


### 6/ Import ####

test_that("import_total_richness_per_clades() reads a file, with the same result as the data.frame",
{
  lists <- toy_lists()
  f <- tempfile(fileext = ".csv")
  utils::write.csv(toy_richness_table(), f, row.names = FALSE)
  imported <- import_total_richness_per_clades(f, phy = phylograftR::toy_backbone, missing_taxa_list = lists, verbose = FALSE)
  expect_identical(imported, build_tr(toy_richness_table(), lists))
  # The argument to return the richness table is passed on
  expect_null(attr(import_total_richness_per_clades(f, phy = phylograftR::toy_backbone, missing_taxa_list = lists, return_richness_table = FALSE, verbose = FALSE), "richness_table"))
  # Text file with another separator
  f2 <- tempfile(fileext = ".txt")
  utils::write.table(toy_richness_table(), f2, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_identical(import_total_richness_per_clades(f2, phy = phylograftR::toy_backbone, missing_taxa_list = lists, verbose = FALSE), imported)
  unlink(c(f, f2))
  expect_error(import_total_richness_per_clades("no_file.csv", phylograftR::toy_backbone, lists), "File not found")
})

test_that("the .csv of the package is imported as the toy total richness dataset",
{
  f <- system.file("extdata", "toy_total_richness.csv", package = "phylograftR")
  expect_true(file.exists(f))
  imported <- import_total_richness_per_clades(f, phy = phylograftR::toy_backbone, missing_taxa_list = phylograftR::toy_missing_taxa_list, verbose = FALSE)
  expect_named(imported, c("Target_clade", "Target_clade_name", "total_richness"))
  expect_identical(imported$total_richness[1:5], phylograftR::toy_total_richness_per_clades$total_richness)
})
