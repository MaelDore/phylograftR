
##### Test suite for build_missing_taxa_list_from_taxonomy() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Small backbone: 10 tips, three genera (monophyletic) and two families
#   Fam1 = GenA (a1..a4) + GenB (b1..b3); Fam2 = GenC (c1..c3), so Fam2 and GenC have the same clade
small_taxonomy <- function ()
{
  phy <- ape::compute.brlen(ape::read.tree(text = "(((((a1,a2),a3),a4),((b1,b2),b3)),((c1,c2),c3));"))
  tips <- phy$tip.label
  bb <- data.frame(Taxa = tips,
                   Family = ifelse(grepl("^c", tips), "Fam2", "Fam1"),
                   Genus = c(a = "GenA", b = "GenB", c = "GenC")[substr(tips, 1, 1)],
                   stringsAsFactors = FALSE)
  list(phy = phy, bb = bb)
}

## Table of missing taxa from a list of columns (one element per taxon)
make_mt <- function (...) data.frame(..., stringsAsFactors = FALSE)

## Node of a clade written with tips
node_of <- function (phy, tips) if (length(tips) == 1) match(tips, phy$tip.label) else ape::getMRCA(phy, tips)

## Build quietly
build_tax <- function (bb, mt, phy, ...) build_missing_taxa_list_from_taxonomy(bb, mt, phy, verbose = FALSE, ...)


test_that("clades are retrieved from names of groups and from names of levels", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("x1", "x2", "x3"), Type = "Taxa", Genus = c("GenA", "GenB", "GenB"), Family = "Fam1",
                Target_clade = c("Genus", "GenA", "Family"))
  mtl <- build_tax(s$bb, mt, s$phy)
  expect_equal(names(mtl), c("Taxon_name", "Target_clade", "Target_clade_name"))
  expect_equal(mtl$Taxon_name, c("x1", "x2", "x3"))
  expect_equal(names(mtl$Target_clade), c("x1", "x2", "x3"))
  # The level 'Genus' of x1 is GenA: same clade as the group GenA (x2)
  expect_equal(mtl$Target_clade$x1, mtl$Target_clade$x2)
  expect_equal(node_of(s$phy, mtl$Target_clade$x1), node_of(s$phy, paste0("a", 1:4)))
  expect_equal(node_of(s$phy, mtl$Target_clade$x3), node_of(s$phy, c(paste0("a", 1:4), paste0("b", 1:3))))
  expect_equal(unname(mtl$Target_clade_name), c("GenA", "GenA", "Fam1"))
  expect_equal(names(mtl$Target_clade_name), c("x1", "x2", "x3"))
  # No excluded subclade: no item
  expect_null(mtl$Excluded_subclades)
  # The output is a valid list for the grafting functions
  expect_no_error(validate_missing_taxa_list(mtl, clade_item = "Target_clade"))
  # Tips used to write the clades are valid tips of the backbone
  expect_true(all(unlist(mtl$Target_clade) %in% s$phy$tip.label))
})

test_that("excluded subclades, stems, allowed states and the stem of the target clade are passed through", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("x1", "x2", "x3"), Type = "taxa", Target_clade = "Fam1",
                Excluded_subclades = c("GenA", "GenA; GenB", NA),
                exclude_subclades_stem = c("FALSE", "TRUE;FALSE", NA),
                include_target_clade_stem = c("TRUE", "FALSE", NA),
                allowed_states = c("S1,S2", NA, "S3"))
  mtl <- build_tax(s$bb, mt, s$phy)
  expect_equal(lengths(mtl$Excluded_subclades), c(x1 = 1, x2 = 2, x3 = 0))
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x2[[1]]), node_of(s$phy, paste0("a", 1:4)))
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x2[[2]]), node_of(s$phy, paste0("b", 1:3)))
  expect_equal(names(mtl$Excluded_subclades$x2), c("GenA", "GenB"))
  expect_equal(unname(mtl$exclude_subclades_stem$x1), FALSE)
  expect_equal(unname(mtl$exclude_subclades_stem$x2), c(TRUE, FALSE))
  expect_null(mtl$exclude_subclades_stem$x3)
  expect_equal(unname(mtl$include_target_clade_stem), c(TRUE, FALSE, NA))
  expect_equal(mtl$allowed_states$x1, c("S1", "S2"))
  expect_null(mtl$allowed_states$x2)
  # A single value of exclude_subclades_stem is used for all the groups
  mt$exclude_subclades_stem <- c("FALSE", "FALSE", NA)
  expect_equal(unname(build_tax(s$bb, mt, s$phy)$exclude_subclades_stem$x2), c(FALSE, FALSE))
  # Wrong number of values
  mt$exclude_subclades_stem <- c("FALSE", "TRUE;FALSE;TRUE", NA)
  expect_error(build_tax(s$bb, mt, s$phy), "one value per group.*x2")
  # Stem values without excluded subclades: ignored with a warning
  mt$exclude_subclades_stem <- c("FALSE", "TRUE", "TRUE")
  expect_warning(build_tax(s$bb, mt, s$phy), "without Excluded_subclades: ignored for x3")
})

test_that("sister clades are built for graft_missing_clades()", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("y1", "y2"), Type = "clades", Genus = c("GenC", "GenD"), Sister_clade = c("GenB", "Genus"))
  # y2 has Genus = GenD: absent from the backbone
  expect_error(build_tax(s$bb, mt, s$phy), "y2 \\(Sister_clade\\).*GenD.*absent from the backbone")
  mt$Genus <- c("GenC", "GenA")
  mtl <- build_tax(s$bb, mt, s$phy)
  expect_equal(names(mtl), c("Taxon_name", "Sister_clade", "Sister_clade_name"))
  expect_equal(unname(mtl$Sister_clade_name), c("GenB", "GenA"))
  expect_equal(node_of(s$phy, mtl$Sister_clade$y1), node_of(s$phy, paste0("b", 1:3)))
  expect_no_error(validate_missing_taxa_list(mtl, clade_item = "Sister_clade"))
  # allowed_states is also passed
  mt$allowed_states <- c("S1", NA)
  expect_equal(build_tax(s$bb, mt, s$phy)$allowed_states$y1, "S1")
})

test_that("both types return one list per type, a single type returns the list", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("x1", "y1", "x2"), Type = c("Taxa", "CLADES", " taxa "), Genus = c("GenA", "GenC", "GenB"),
                Target_clade = c("Genus", NA, "Genus"), Sister_clade = c(NA, "GenB", NA))
  res <- build_tax(s$bb, mt, s$phy)
  expect_equal(names(res), c("taxa", "clades"))
  expect_equal(res$taxa$Taxon_name, c("x1", "x2"))
  expect_equal(res$clades$Taxon_name, "y1")
  expect_equal(names(res$taxa), c("Taxon_name", "Target_clade", "Target_clade_name"))
  expect_equal(names(res$clades), c("Taxon_name", "Sister_clade", "Sister_clade_name"))
  # Rows of the other type are not used: no warning if the columns are blank for them, one if they are not
  mt$Sister_clade[1] <- "GenA"
  expect_warning(build_tax(s$bb, mt, s$phy), "'Sister_clade' is not used with Type = 'taxa': ignored for x1")
  mt$Sister_clade[1] <- NA
  mt$Target_clade[2] <- "GenA"
  expect_warning(build_tax(s$bb, mt, s$phy), "'Target_clade' is not used with Type = 'clades': ignored for y1")
  # Single type
  only_taxa <- build_tax(s$bb, mt[c(1, 3), ], s$phy)
  expect_equal(names(only_taxa)[1], "Taxon_name")
  expect_equal(only_taxa$Taxon_name, c("x1", "x2"))
  # Information on the groups
  info <- attr(res, "taxonomic_groups")
  expect_s3_class(info, "data.frame")
  expect_equal(names(info), c("Group", "Level", "n_members", "clade_node", "clade_tips", "n_tips_clade", "n_foreign", "n_unannotated", "monophyletic"))
  expect_equal(sort(info$Group), c("GenA", "GenB"))
  expect_equal(info$n_members[info$Group == "GenA"], 4)
  expect_equal(info$clade_node[info$Group == "GenB"], node_of(s$phy, paste0("b", 1:3)))
  expect_true(all(info$monophyletic))
  expect_equal(info$Level[info$Group == "GenB"], "Genus")
  expect_equal(attr(only_taxa, "taxonomic_groups")$Group, c("GenA", "GenB"))
})

test_that("groups of the same clade share their names", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("x1", "x2", "x3"), Type = "Taxa", Family = "Fam2", Genus = "GenC", Target_clade = c("Family", "GenC", "Fam1"),
                Excluded_subclades = c(NA, NA, "GenB"))
  expect_warning(mtl <- build_tax(s$bb, mt, s$phy), "Several groups have the same clade.*Fam2 / GenC \\(clade of Fam2, GenC\\)")
  # No warning when each clade has a single group
  expect_no_warning(build_tax(s$bb, mt[c(1, 3), ], s$phy))
  # Fam2 and GenC have the same clade (same MRCA node): same names, whatever the group used
  expect_equal(unname(mtl$Target_clade_name[1:2]), rep("Fam2 / GenC", 2))
  expect_equal(mtl$Target_clade$x1, mtl$Target_clade$x2)
  expect_equal(unname(mtl$Target_clade_name[3]), "Fam1")
})

test_that("non-monophyletic target groups: intruders are excluded, kept ('mrca') or rejected ('error')", {
  s <- small_taxonomy()
  # Fam1 loses b1 and b2, which belong to another family: Fam1 = GenA + b3, whose clade holds the pair (b1, b2)
  bb <- s$bb
  bb$Family[bb$Taxa %in% c("b1", "b2")] <- "FamX"
  mt <- make_mt(Taxa = c("x1", "x2"), Type = "Taxa", Target_clade = c("Fam1", "GenA"))
  expect_warning(mtl <- build_tax(bb, mt, s$phy), "Target groups are not monophyletic.*Fam1 \\(Family\\): 2 tip\\(s\\) of other groups among 7 tips")
  expect_equal(node_of(s$phy, mtl$Target_clade$x1), node_of(s$phy, c(paste0("a", 1:4), "b3")))
  expect_equal(lengths(mtl$Excluded_subclades), c(x1 = 1, x2 = 0))
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x1[[1]]), node_of(s$phy, c("b1", "b2")))
  expect_equal(names(mtl$Excluded_subclades$x1), "foreign(b1,b2)")
  # The stem of intruders is excluded
  expect_equal(unname(mtl$exclude_subclades_stem$x1), TRUE)
  expect_null(mtl$exclude_subclades_stem$x2)
  info <- attr(mtl, "taxonomic_groups")
  expect_equal(info$monophyletic[info$Group == "Fam1"], FALSE)
  expect_equal(info$n_foreign[info$Group == "Fam1"], 2)
  expect_no_error(validate_missing_taxa_list(mtl, clade_item = "Target_clade"))

  # 'mrca': the MRCA clade is used as such
  expect_warning(mrca <- build_tax(bb, mt, s$phy, non_monophyletic = "mrca"), "their clade \\(MRCA of their tips\\) is used and includes tips of other groups")
  expect_equal(mrca$Target_clade, mtl$Target_clade)
  expect_null(mrca$Excluded_subclades)
  # 'error'
  expect_error(build_tax(bb, mt, s$phy, non_monophyletic = "error"), "not monophyletic.*Target clades:.*Fam1 \\(Family\\)")
  expect_error(build_tax(bb, mt, s$phy, non_monophyletic = "none"), "should be one of")

  # Excluded subclades of the user are kept with their stem, and the intruders are added after them (stem excluded)
  mt$Excluded_subclades <- c("GenA", NA)
  mt$exclude_subclades_stem <- c("FALSE", NA)
  expect_warning(both <- build_tax(bb, mt, s$phy), "not monophyletic")
  expect_equal(lengths(both$Excluded_subclades)[["x1"]], 2)
  expect_equal(node_of(s$phy, both$Excluded_subclades$x1[[1]]), node_of(s$phy, paste0("a", 1:4)))
  expect_equal(unname(both$exclude_subclades_stem$x1), c(FALSE, TRUE))
  expect_equal(names(both$Excluded_subclades$x1), c("GenA", "foreign(b1,b2)"))

  # An intruder nested within a subclade excluded by the user is not added: GenB holds b1 and b2
  mt$Excluded_subclades <- c("GenB", NA)
  mt$exclude_subclades_stem <- NA
  expect_warning(nested <- build_tax(bb, mt, s$phy), "not monophyletic")
  expect_equal(lengths(nested$Excluded_subclades)[["x1"]], 1)
  expect_equal(node_of(s$phy, nested$Excluded_subclades$x1[[1]]), node_of(s$phy, paste0("b", 1:3)))
  expect_equal(names(nested$Excluded_subclades$x1), "GenB")
})

test_that("unannotated tips are compatible with every group", {
  s <- small_taxonomy()
  bb <- s$bb
  # a3 has no genus, and the table has no row for a2
  bb$Genus[bb$Taxa == "a3"] <- NA
  bb2 <- bb[bb$Taxa != "a2", ]
  mt <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "GenA")
  expect_no_warning(mtl <- build_tax(bb, mt, s$phy))
  info <- attr(mtl, "taxonomic_groups")
  expect_true(info$monophyletic)
  expect_equal(info$n_members, 3)
  expect_equal(info$n_unannotated, 1)
  expect_equal(info$n_tips_clade, 4)
  expect_null(mtl$Excluded_subclades)
  # Tips of the phylogeny absent from the table: warning, and they are never excluded
  mt2 <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "GenA")
  expect_warning(mtl2 <- build_tax(bb2, mt2, s$phy), "1 tips of 'phy' are absent from 'backbone_taxonomic_table'.*a2")
  expect_true(attr(mtl2, "taxonomic_groups")$monophyletic)
  expect_equal(attr(mtl2, "taxonomic_groups")$n_unannotated, 2) # a2 (absent) and a3 (no genus)
})

test_that("non-monophyletic excluded groups are split into their largest pure subclades", {
  s <- small_taxonomy()
  # a3 belongs to GenB: GenA = (a1, a2) + a4 in the clade (a1, a2, a3, a4)
  bb <- s$bb
  bb$Genus[bb$Taxa == "a3"] <- "GenB"
  mt <- make_mt(Taxa = c("x1", "x2"), Type = "Taxa", Target_clade = "Fam1", Excluded_subclades = c("GenA", "GenB"))
  expect_warning(mtl <- build_tax(bb, mt, s$phy), "Excluded groups are not monophyletic.*GenA \\(Genus\\)")
  # x1: GenA is replaced by the clade (a1, a2) and the tip a4. a3 is not excluded
  expect_equal(lengths(mtl$Excluded_subclades)[["x1"]], 2)
  expect_equal(names(mtl$Excluded_subclades$x1), c("GenA [1]", "GenA [2]"))
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x1[[1]]), node_of(s$phy, c("a1", "a2")))
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x1[[2]]), node_of(s$phy, "a4"))
  expect_false("a3" %in% unlist(mtl$Excluded_subclades$x1))
  # x2: GenB = a3 + (b1, b2, b3) is also split, into the tip a3 and the clade of b1 to b3
  expect_equal(lengths(mtl$Excluded_subclades)[["x2"]], 2)
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x2[[1]]), node_of(s$phy, "a3"))
  expect_equal(node_of(s$phy, mtl$Excluded_subclades$x2[[2]]), node_of(s$phy, paste0("b", 1:3)))
  # 'mrca': the clade of the MRCA is excluded
  expect_warning(mrca <- build_tax(bb, mt[1, ], s$phy, non_monophyletic = "mrca"), "Excluded groups are not monophyletic.*tips of other groups are excluded too")
  expect_equal(lengths(mrca$Excluded_subclades)[["x1"]], 1)
  expect_equal(node_of(s$phy, mrca$Excluded_subclades$x1[[1]]), node_of(s$phy, paste0("a", 1:4)))
  # 'error'
  expect_error(build_tax(bb, mt, s$phy, non_monophyletic = "error"), "Excluded subclades:")
})

test_that("non-monophyletic sister groups use the clade of the MRCA, with a warning", {
  s <- small_taxonomy()
  bb <- s$bb
  bb$Genus[bb$Taxa == "a3"] <- "GenB"
  mt <- make_mt(Taxa = "y1", Type = "clades", Sister_clade = "GenA")
  expect_warning(mtl <- build_tax(bb, mt, s$phy), "Sister groups are not monophyletic.*cannot be excluded from a sister clade.*GenA")
  expect_equal(node_of(s$phy, mtl$Sister_clade$y1), node_of(s$phy, paste0("a", 1:4)))
  expect_warning(build_tax(bb, mt, s$phy, non_monophyletic = "mrca"), "Sister groups are not monophyletic")
  expect_error(build_tax(bb, mt, s$phy, non_monophyletic = "error"), "Sister clades:")
})

test_that("a sister clade holding the whole backbone is an error", {
  s <- small_taxonomy()
  bb <- cbind(s$bb, Order = "Ord1")
  mt <- make_mt(Taxa = "y1", Type = "clades", Sister_clade = "Ord1")
  expect_error(build_tax(bb, mt, s$phy), "holds the whole backbone.*Ord1")
  # Allowed for a target clade
  mt <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "Ord1")
  expect_no_error(build_tax(bb, mt, s$phy))
})

test_that("groups that cannot be found are reported for all the taxa and clades", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("x1", "x2", "x3", "x4"), Type = "Taxa", Genus = c("GenA", NA, "GenZ", "GenA"),
                Target_clade = c("Foo", "Genus", "Genus", "GenA"), Excluded_subclades = c(NA, NA, NA, "GenB; Bar"))
  err <- tryCatch(build_tax(s$bb, mt, s$phy), error = function (e) conditionMessage(e))
  expect_match(err, "4 problem\\(s\\) found")
  expect_match(err, "x1 \\(Target_clade\\): 'Foo' is neither a taxonomic level nor a group")
  expect_match(err, "x2 \\(Target_clade\\): 'Genus' is a level, but the taxon has no value at this level")
  expect_match(err, "x3 \\(Target_clade\\): group 'GenZ' \\(Genus\\) is absent from the backbone")
  expect_match(err, "x4 \\(Excluded_subclades\\): 'Bar' is neither")
  # Level used but the column is absent from the table of the missing taxa
  mt2 <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "Genus")
  expect_error(build_tax(s$bb, mt2, s$phy), "'Genus' is a level, but 'missing_taxa_taxonomic_table' has no column 'Genus'")
  # Empty clade, and empty group in a list
  mt3 <- make_mt(Taxa = c("x1", "x2"), Type = "Taxa", Target_clade = c(NA, "GenA"), Excluded_subclades = c(NA, "GenB;;"))
  err3 <- tryCatch(build_tax(s$bb, mt3, s$phy), error = function (e) conditionMessage(e))
  expect_match(err3, "x1 \\(Target_clade\\): empty")
  expect_match(err3, "x2 \\(Excluded_subclades\\): empty group")
})

test_that("ambiguous names in the tables are rejected", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "GenA")
  # A level named like a tip
  bb <- s$bb; names(bb)[names(bb) == "Genus"] <- "a1"
  expect_error(build_tax(bb, mt, s$phy), "Ambiguous names.*levels that are also names of taxa: a1")
  # A group named like a tip
  bb <- s$bb; bb$Genus[bb$Taxa == "c1"] <- "b2"
  expect_error(build_tax(bb, mt, s$phy), "groups that are also names of taxa: b2")
  # A group named like a level
  bb <- s$bb; bb$Genus[bb$Taxa == "c1"] <- "Family"
  expect_error(build_tax(bb, mt, s$phy), "levels that are also names of groups: Family")
  # A group in two levels
  bb <- s$bb; bb$Genus[bb$Taxa == "c1"] <- "Fam1"
  expect_error(build_tax(bb, mt, s$phy), "several levels.*Fam1 \\(Family, Genus\\)")
  # A group of the table of the missing taxa that is also a group of another level
  mt2 <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "GenA", Genus = "Fam1")
  expect_error(build_tax(s$bb, mt2, s$phy), "several levels.*Fam1")
  # A group with a semicolon
  bb <- s$bb; bb$Genus[bb$Taxa == "c1"] <- "Gen;C"
  expect_error(build_tax(bb, mt, s$phy), "containing ';'")
  # Level names that are reserved
  bb <- s$bb; names(bb)[names(bb) == "Genus"] <- "Type"
  expect_error(build_tax(bb, mt, s$phy), "Names of taxonomic levels cannot be one of")
  # Duplicated columns
  bb <- s$bb; names(bb)[3] <- "Family"
  expect_error(build_tax(bb, mt, s$phy), "same name")
})

test_that("tables and phylogeny are checked", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "GenA")
  expect_error(build_missing_taxa_list_from_taxonomy(s$bb, mt, "not a tree", verbose = FALSE), "'phy' must be a phylogeny")
  expect_error(build_tax(as.matrix(s$bb), mt, s$phy), "'backbone_taxonomic_table' must be a data frame")
  expect_error(build_tax(s$bb, "mt", s$phy), "'missing_taxa_taxonomic_table' must be a data frame")
  expect_error(build_tax(s$bb[, -1], mt, s$phy), "must have a column 'Taxa'")
  expect_error(build_tax(s$bb[, 1, drop = FALSE], mt, s$phy), "at least one taxonomic level")
  # Taxa of the backbone absent from the phylogeny
  bb <- rbind(s$bb, data.frame(Taxa = "zz", Family = "Fam1", Genus = "GenA"))
  expect_error(build_tax(bb, mt, s$phy), "1 taxa of 'backbone_taxonomic_table' are not tips of 'phy': zz")
  bb <- rbind(s$bb, s$bb[1, ])
  expect_error(build_tax(bb, mt, s$phy), "must be unique in 'backbone_taxonomic_table'.*a1")
  # Taxa to graft
  expect_error(build_tax(s$bb, make_mt(Taxa = "a1", Type = "Taxa", Target_clade = "GenA"), s$phy), "absent from the backbone.*a1")
  expect_error(build_tax(s$bb, make_mt(Taxa = c("x", "x"), Type = "Taxa", Target_clade = "GenA"), s$phy), "must be unique in 'missing_taxa_taxonomic_table'.*x")
  expect_error(build_tax(s$bb, make_mt(Taxa = "x", Target_clade = "GenA"), s$phy), "must have a column 'Type'")
  expect_error(build_tax(s$bb, make_mt(Type = "Taxa", Target_clade = "GenA"), s$phy), "must have a column 'Taxa'")
  expect_error(build_tax(s$bb, make_mt(Taxa = c("x", "y"), Type = c("Taxa", "Genus"), Target_clade = "GenA"), s$phy), "'Type' must be 'Taxa' or 'clades'.*y \\('Genus'\\)")
  expect_error(build_tax(s$bb, make_mt(Taxa = "x", Type = NA, Target_clade = "GenA"), s$phy), "empty")
  expect_error(build_tax(s$bb, make_mt(Taxa = "x", Type = "Taxa", Sister_clade = "GenA"), s$phy), "needs a column 'Target_clade'")
  expect_error(build_tax(s$bb, make_mt(Taxa = "x", Type = "clades", Target_clade = "GenA"), s$phy), "needs a column 'Sister_clade'")
  # Unknown columns are ignored with a warning
  expect_warning(build_tax(s$bb, cbind(mt, Foo = "bar"), s$phy), "Unknown columns ignored.*'Foo'")
  # Blank and "NA" cells are missing values; factors are accepted
  bb <- s$bb
  bb$Genus[bb$Taxa == "a3"] <- ""
  bb$Family <- factor(bb$Family)
  info <- attr(build_tax(bb, mt, s$phy), "taxonomic_groups")
  expect_equal(info$n_unannotated, 1)
})

test_that("taxa of the missing table can have no taxonomic level at all", {
  s <- small_taxonomy()
  mt <- make_mt(Taxa = c("x1", "x2"), Type = "Taxa", Target_clade = c("GenA", "Fam2"))
  mtl <- build_tax(s$bb, mt, s$phy)
  expect_equal(unname(mtl$Target_clade_name), c("GenA", "Fam2"))
})

test_that("clades match the MRCA of random groups, and foreign tips are excluded as required", {
  phy <- make_toy_phy(n_tips = 40, seed = 77)
  n <- ape::Ntip(phy)
  set.seed(5)
  label <- sample(c("G1", "G2", "G3", "G4", NA), n, replace = TRUE, prob = c(1, 1, 1, 1, 0.5))
  bb <- data.frame(Taxa = sample(phy$tip.label), stringsAsFactors = FALSE)
  bb$Genus <- label[match(bb$Taxa, phy$tip.label)]
  mt <- make_mt(Taxa = paste0("x", 1:4), Type = "Taxa", Target_clade = paste0("G", 1:4))
  mtl <- suppressWarnings(build_tax(bb, mt, phy))
  parts <- ape::prop.part(phy)
  for (g in paste0("G", 1:4))
  {
    members <- bb$Taxa[which(bb$Genus == g)]
    node <- ape::getMRCA(phy, members)
    tips_clade <- phy$tip.label[parts[[node - n]]]
    # Clade of the target: same node as the MRCA of the members
    expect_equal(node_of(phy, mtl$Target_clade[[paste0("x", substr(g, 2, 2))]]), node)
    foreign <- tips_clade[!is.na(bb$Genus[match(tips_clade, bb$Taxa)]) & bb$Genus[match(tips_clade, bb$Taxa)] != g]
    excluded <- mtl$Excluded_subclades[[paste0("x", substr(g, 2, 2))]]
    excluded_tips <- unlist(lapply(excluded, function (tips) if (length(tips) == 1) tips else phy$tip.label[parts[[node_of(phy, tips) - n]]]))
    # All the foreign tips are excluded, and only them
    expect_setequal(excluded_tips, foreign)
    # Excluded subclades are disjoint, and each is a maximal pure subclade (its sister clade holds no more foreign tips only)
    expect_false(anyDuplicated(excluded_tips) > 0)
    info <- attr(mtl, "taxonomic_groups")
    expect_equal(info$n_foreign[info$Group == g], length(foreign))
    expect_equal(info$n_unannotated[info$Group == g], sum(is.na(bb$Genus[match(tips_clade, bb$Taxa)])))
  }
})

test_that("works with polytomies and with a phylogeny whose tips are not in the order of the table", {
  phy <- make_toy_polytomies(make_toy_phy(n_tips = 20, seed = 3), n_collapsed = 3)
  tips <- phy$tip.label
  bb <- data.frame(Taxa = rev(tips), Genus = rep(c("G1", "G2"), each = 10), stringsAsFactors = FALSE)
  mt <- make_mt(Taxa = "x1", Type = "Taxa", Target_clade = "G1")
  mtl <- suppressWarnings(build_tax(bb, mt, phy))
  members <- bb$Taxa[bb$Genus == "G1"]
  expect_equal(node_of(phy, mtl$Target_clade$x1), ape::getMRCA(phy, members))
})

test_that("the list is accepted by the grafting functions", {
  toy <- toy_setup()
  tipsB <- toy$tipsB; tipsC <- toy$tipsC
  bb <- data.frame(Taxa = toy$phy$tip.label, stringsAsFactors = FALSE)
  bb$Family <- ifelse(bb$Taxa %in% toy$tipsA, "FamA", "FamO")
  bb$Genus <- NA_character_
  bb$Genus[bb$Taxa %in% tipsB] <- "GenB"
  bb$Genus[bb$Taxa %in% tipsC] <- "GenC"
  bb$Genus[bb$Taxa %in% setdiff(toy$tipsA, c(tipsB, tipsC))] <- "GenR"
  # Group of A without one tip of B and one tip of C: not monophyletic
  bb$Grp <- NA_character_
  bb$Grp[bb$Taxa %in% toy$tipsA] <- "GrpA"
  bb$Grp[bb$Taxa %in% c(tipsB[1], tipsC[1])] <- "GrpX"
  mt <- make_mt(Taxa = c("X1", "X2", "X3", "X4", "Y1", "Y2"),
                Type = c("Taxa", "Taxa", "Taxa", "Taxa", "clades", "clades"),
                Family = "FamA", Genus = c("GenB", "GenB", "GenC", NA, "GenC", "GenB"),
                Target_clade = c("Genus", "FamA", "FamA", "GenC", NA, NA),
                Excluded_subclades = c(NA, "GenC", "GenB; GenC", NA, NA, NA),
                Sister_clade = c(NA, NA, NA, NA, "GenB", "GenC"))
  # GrpA is not monophyletic but it is not used here: no warning
  expect_no_warning(mtls <- build_tax(bb, mt, toy$phy))
  expect_equal(names(mtls), c("taxa", "clades"))
  expect_equal(names(attr(mtls, "taxonomic_groups")$monophyletic), NULL)

  g_taxa <- quiet_graft(graft_missing_taxa, toy$phy, mtls$taxa, method = "uniform", seed = 1)
  expect_valid_graft(g_taxa, toy$phy, mtls$taxa)
  g_clades <- quiet_graft(graft_missing_clades, toy$phy, mtls$clades, seed = 1)
  expect_true(all(mtls$clades$Taxon_name %in% g_clades$tip.label))
  expect_true(ape::is.ultrametric(g_clades))
  # The taxon X3 is grafted in A but neither in B nor in C
  parts <- ape::prop.part(g_taxa)
  n <- ape::Ntip(g_taxa)
  clade_tips <- function (tips) g_taxa$tip.label[parts[[ape::getMRCA(g_taxa, tips) - n]]]
  expect_false("X3" %in% clade_tips(tipsB))
  expect_false("X3" %in% clade_tips(tipsC))

  # With a non-monophyletic target (GrpA) the intruders (one tip of B and one tip of C) are excluded
  mt_r <- make_mt(Taxa = "R1", Type = "Taxa", Target_clade = "GrpA")
  expect_warning(mtl_r <- build_tax(bb, mt_r, toy$phy), "GrpA")
  expect_setequal(unlist(mtl_r$Excluded_subclades$R1), c(tipsB[1], tipsC[1]))
  g_r <- quiet_graft(graft_missing_taxa, toy$phy, mtl_r, method = "uniform", seed = 1)
  expect_valid_graft(g_r, toy$phy, mtl_r)
  # R1 is not grafted on the terminal branch of an intruder tip (it would then be sister to it)
  siblings <- function (phy, tip) { p <- phy$edge[match(match(tip, phy$tip.label), phy$edge[, 2]), 1]; phy$edge[phy$edge[, 1] == p, 2] }
  r1 <- match("R1", g_r$tip.label)
  expect_false(r1 %in% siblings(g_r, tipsB[1]))
  expect_false(r1 %in% siblings(g_r, tipsC[1]))

  # Names of the clades are consistent with total_richness_per_clades using the same names
  richness <- list(Target_clade = list(mtls$taxa$Target_clade$X1, mtls$taxa$Target_clade$X2, mtls$taxa$Target_clade$X4),
                   Target_clade_name = c("GenB", "FamA", "GenC"), total_richness = c(12, 60, 20))
  expect_no_error(quiet_graft(graft_missing_taxa, toy$phy, mtls$taxa, method = "uniform", total_richness_per_clades = richness, seed = 1))
})
