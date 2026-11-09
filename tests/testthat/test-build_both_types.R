##### Test suite for tables describing both types of taxa (build_missing_taxa_list() and import_missing_taxa_list()) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Toy backbone of the package (30 tips t1-t30). Clades: t26,t18 = clade A (17 tips) ; t8,t22 = clade B (13 tips) ; t5,t12 = clade C nested in B ;
# t10,t19 and t26,t4 = disjoint subclades nested in A

# Table with both types of taxa: 2 taxa with a target clade, 2 taxa grafted together as a missing clade
both_table <- function ()
{
  data.frame(Taxon_name = c("X1", "X2", "Y1", "Y2"),
             Type = c("taxa", "taxa", "clades", "clades"),
             Target_clade = c("t26,t18", "t8,t22", NA, NA),
             Sister_clade = c(NA, NA, "t5,t12", "t5,t12"),
             Excluded_subclades = c("t10,t19", NA, NA, NA),
             allowed_states = c("a", NA, "b", "b"),
             stringsAsFactors = FALSE)
}
phy_both <- function () phylograftR::toy_backbone


### 1/ Two lists are returned ####

test_that("build_missing_taxa_list() returns one list per type when the table has both types",
{
  lists <- build_missing_taxa_list(both_table(), phy_both(), verbose = FALSE)
  expect_type(lists, "list")
  expect_named(lists, c("taxa", "clades"))
  expect_named(lists$taxa, c("Taxon_name", "Target_clade", "Excluded_subclades", "allowed_states"))
  expect_named(lists$clades, c("Taxon_name", "Sister_clade", "allowed_states"))
  expect_identical(lists$taxa$Taxon_name, c("X1", "X2"))
  expect_identical(lists$clades$Taxon_name, c("Y1", "Y2"))
  expect_identical(lists$taxa$Target_clade, list(X1 = c("t26", "t18"), X2 = c("t8", "t22")))
  expect_identical(lists$clades$Sister_clade, list(Y1 = c("t5", "t12"), Y2 = c("t5", "t12")))
  expect_identical(lists$taxa$Excluded_subclades$X1, list(c("t10", "t19")))
  expect_null(lists$taxa$Excluded_subclades$X2)
  expect_identical(lists$clades$allowed_states, list(Y1 = "b", Y2 = "b"))
  # Each list is the same as the one built from the sub-table of its type, alone
  d <- both_table()
  expect_identical(lists$taxa, build_missing_taxa_list(d[1:2, c("Taxon_name", "Target_clade", "Excluded_subclades", "allowed_states")], phy_both(), verbose = FALSE))
  expect_identical(lists$clades, build_missing_taxa_list(d[3:4, c("Taxon_name", "Sister_clade", "allowed_states")], phy_both(), verbose = FALSE))
})

test_that("the type can be inferred from the clade column that is filled in each row",
{
  d <- both_table()
  expected <- build_missing_taxa_list(d, phy_both(), verbose = FALSE)
  d$Type <- NULL
  expect_identical(build_missing_taxa_list(d, phy_both(), verbose = FALSE), expected)
  # The case of 'Type' is ignored
  d <- both_table(); d$Type <- toupper(d$Type)
  expect_identical(build_missing_taxa_list(d, phy_both(), verbose = FALSE), expected)
})

test_that("the rows can be in any order",
{
  d <- both_table()[c(3, 1, 4, 2), ]
  lists <- build_missing_taxa_list(d, phy_both(), verbose = FALSE)
  expect_identical(lists$taxa$Taxon_name, c("X1", "X2"))
  expect_identical(lists$clades$Taxon_name, c("Y1", "Y2"))
})

test_that("optional columns empty for the taxa of a type are not kept in its list",
{
  d <- both_table(); d$Target_clade_name <- c("A", "B", NA, NA); d$Sister_clade_name <- c(NA, NA, "C", "C")
  lists <- build_missing_taxa_list(d, phy_both(), verbose = FALSE)
  expect_identical(lists$taxa$Target_clade_name, c(X1 = "A", X2 = "B"))
  expect_identical(lists$clades$Sister_clade_name, c(Y1 = "C", Y2 = "C"))
  expect_null(lists$clades$Excluded_subclades)
  d$allowed_states <- NA
  lists <- build_missing_taxa_list(d, phy_both(), verbose = FALSE)
  expect_null(lists$taxa$allowed_states)
  expect_null(lists$clades$allowed_states)
})


### 2/ A single list is still returned with a single type ####

test_that("a single list is returned when all the taxa have the same type, with or without 'Type'",
{
  d <- both_table()
  d1 <- d[1:2, ]
  single <- build_missing_taxa_list(d1, phy_both(), verbose = FALSE)
  expect_named(single, c("Taxon_name", "Target_clade", "Excluded_subclades", "allowed_states"))
  expect_identical(single, build_missing_taxa_list(d1[, c("Taxon_name", "Target_clade", "Excluded_subclades", "allowed_states")], phy_both(), verbose = FALSE))
  d2 <- d[3:4, ]
  single2 <- build_missing_taxa_list(d2, phy_both(), verbose = FALSE)
  expect_named(single2, c("Taxon_name", "Sister_clade", "allowed_states"))
  # Type provided, and a single clade column
  d3 <- d1[, c("Taxon_name", "Type", "Target_clade")]
  expect_identical(build_missing_taxa_list(d3, phy_both(), verbose = FALSE),
                   build_missing_taxa_list(d3[, c("Taxon_name", "Target_clade")], phy_both(), verbose = FALSE))
})

test_that("a table with a single clade column and without 'Type' behaves as before",
{
  d <- both_table()[1:2, c("Taxon_name", "Target_clade", "Excluded_subclades")]
  expect_named(build_missing_taxa_list(d, phy_both(), verbose = FALSE), c("Taxon_name", "Target_clade", "Excluded_subclades"))
  expect_error(build_missing_taxa_list(d[, c("Taxon_name", "Excluded_subclades")], phy_both(), verbose = FALSE), "'Target_clade'.*'Sister_clade'")
})


### 3/ Checks ####

test_that("the type of each row must be identifiable",
{
  d <- both_table(); d$Type <- NULL
  d2 <- d; d2$Sister_clade[1] <- "t5,t12"
  expect_error(build_missing_taxa_list(d2, phy_both(), verbose = FALSE), "both 'Target_clade' and 'Sister_clade' filled.*X1' \\(row 1\\)")
  d2 <- d; d2$Target_clade[2] <- NA
  expect_error(build_missing_taxa_list(d2, phy_both(), verbose = FALSE), "Neither 'Target_clade' nor 'Sister_clade'.*X2' \\(row 2\\)")
  d2 <- both_table(); d2$Type[3] <- "clade"
  expect_error(build_missing_taxa_list(d2, phy_both(), verbose = FALSE), "'Type' must be 'taxa' or 'clades'.*Y1' \\(row 3\\)")
  d2 <- both_table(); d2$Type[1] <- NA
  expect_error(build_missing_taxa_list(d2, phy_both(), verbose = FALSE), "'Type' must be 'taxa' or 'clades'.*X1' \\(row 1\\)")
  d2 <- both_table()[, names(both_table()) != "Sister_clade"]
  expect_error(build_missing_taxa_list(d2, phy_both(), verbose = FALSE), "needs a column 'Sister_clade'")
  d2 <- both_table()[, names(both_table()) != "Target_clade"]
  expect_error(build_missing_taxa_list(d2, phy_both(), verbose = FALSE), "needs a column 'Target_clade'")
})

test_that("a clade written in the column of the other type is ignored with a warning",
{
  d <- both_table(); d$Sister_clade[1] <- "t5,t12"
  expect_warning(lists <- build_missing_taxa_list(d, phy_both(), verbose = FALSE), "'Sister_clade' is not used with Type = 'taxa'.*X1' \\(row 1\\)")
  expect_identical(lists, build_missing_taxa_list(both_table(), phy_both(), verbose = FALSE))
})

test_that("columns of one type filled for the taxa of the other type are ignored with a warning",
{
  d <- both_table(); d$Excluded_subclades[3] <- "t5,t12"
  expect_warning(lists <- build_missing_taxa_list(d, phy_both(), verbose = FALSE), "not used with Type = 'clades' ignored: 'Excluded_subclades'")
  expect_null(lists$clades$Excluded_subclades)
  # Unknown columns are reported once
  d <- both_table(); d$typo <- 1
  expect_warning(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "Unknown columns ignored: 'typo'")
  w <- testthat::capture_warnings(build_missing_taxa_list(d, phy_both(), verbose = FALSE))
  expect_length(w, 1)
})

test_that("taxon names must be unique and absent from the backbone across both types",
{
  d <- both_table(); d$Taxon_name[3] <- "X1"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "must be unique. Duplicated: X1")
  d <- both_table(); d$Taxon_name[4] <- "t1"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "absent from the backbone.*t1")
  d <- both_table(); d$Taxon_name[2] <- NA
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "'Taxon_name' is empty in the row\\(s\\): 2")
})

test_that("messages refer to the rows of the original table",
{
  d <- both_table(); d$Sister_clade[4] <- "t5,zzz"; d$Sister_clade[3] <- "t5,zzz"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "Y1 \\(row 3\\)")
  d <- both_table(); d$Target_clade[2] <- "t8,zzz"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "X2 \\(row 2\\)")
})

test_that("invalid constraints of a type are caught as with a single type",
{
  # Excluded subclade identical to the target clade: an error, as with a single type
  d <- both_table(); d$Excluded_subclades[1] <- "t26,t18"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "identical to, or contain")
})


### 4/ Same clades, same tips, same labels in both lists ####

test_that("the same clade is written with the same tips in both lists when node IDs are used",
{
  phy <- phy_both()
  node_A <- ape::getMRCA(phy, c("t26", "t18"))
  d <- data.frame(Taxon_name = c("X1", "Y1"), Type = c("taxa", "clades"),
                  Target_clade = c("t26,t18", NA), Sister_clade = c(NA, as.character(node_A)),
                  stringsAsFactors = FALSE)
  lists <- build_missing_taxa_list(d, phy, verbose = FALSE)
  expect_identical(unname(lists$clades$Sister_clade), unname(lists$taxa$Target_clade))
  # And the other way round
  d <- data.frame(Taxon_name = c("X1", "Y1"), Type = c("taxa", "clades"),
                  Target_clade = c(as.character(node_A), NA), Sister_clade = c(NA, "t26,t18"),
                  stringsAsFactors = FALSE)
  lists <- build_missing_taxa_list(d, phy, verbose = FALSE)
  expect_identical(unname(lists$clades$Sister_clade), unname(lists$taxa$Target_clade))
})

test_that("clade labels are checked across the two lists",
{
  d <- data.frame(Taxon_name = c("X1", "Y1"), Type = c("taxa", "clades"),
                  Target_clade = c("t26,t18", NA), Target_clade_name = c("A", NA),
                  Sister_clade = c(NA, "t26,t18"), Sister_clade_name = c(NA, "A"),
                  stringsAsFactors = FALSE)
  expect_no_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE))
  d$Sister_clade_name[2] <- "other"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "Inconsistent labels in 'Target_clade_name and Sister_clade_name'")
  # The same label for two different clades
  d$Sister_clade_name[2] <- "A"; d$Sister_clade[2] <- "t8,t22"
  expect_error(build_missing_taxa_list(d, phy_both(), verbose = FALSE), "Inconsistent labels")
})


### 5/ Import and use ####

test_that("import_missing_taxa_list() returns two lists for a file with both types",
{
  f <- tempfile(fileext = ".csv")
  utils::write.csv(both_table(), f, row.names = FALSE)
  imported <- import_missing_taxa_list(f, phy = phy_both(), verbose = FALSE)
  expect_identical(imported, build_missing_taxa_list(both_table(), phy_both(), verbose = FALSE))
  expect_named(imported, c("taxa", "clades"))
  unlink(f)
})

test_that("both lists can be used in the grafting functions, one after the other",
{
  d <- both_table(); d$allowed_states <- NA
  lists <- build_missing_taxa_list(d, phy_both(), verbose = FALSE)
  g1 <- graft_missing_taxa(phy_both(), lists$taxa, method = "uniform", seed = 1, verbose = FALSE)
  expect_true(all(lists$taxa$Taxon_name %in% g1$tip.label))
  g2 <- graft_missing_clades(phy = g1, missing_taxa_list = lists$clades, seed = 1, verbose = FALSE)
  expect_true(all(c(lists$taxa$Taxon_name, lists$clades$Taxon_name) %in% g2$tip.label))
})
