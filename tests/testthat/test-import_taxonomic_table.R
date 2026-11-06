##### Test suite for import_taxonomic_table() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

backbone_path <- function () system.file("extdata", "toy_backbone_taxonomic_table.csv", package = "phylograftR")
missing_path <- function () system.file("extdata", "toy_missing_taxa_taxonomic_table.csv", package = "phylograftR")

# Use a new empty directory in the temporary directory of the session (removed at its end)
temp_dir <- function () { d <- tempfile("phylograftR_"); dir.create(d); d }

# Small tables
bb_table <- function ()
{
  data.frame(Taxa = c("t1", "t2", "t3"), Family = c("F1", "F1", NA), Genus = c("G1", NA, "G2"), stringsAsFactors = FALSE)
}
mt_table <- function ()
{
  data.frame(Taxa = c("X1", "X2"), Type = c("taxa", "Clades"), Family = c("F1", "F1"),
             Target_clade = c("Family", NA), Sister_clade = c(NA, "Genus"), stringsAsFactors = FALSE)
}

### 1/ Datasets ####

test_that("import_taxonomic_table() reads the .csv files of the package as the toy datasets",
{
  expect_true(file.exists(backbone_path()))
  expect_true(file.exists(missing_path()))
  bb <- import_taxonomic_table(backbone_path(), type = "backbone", phy = phylograftR::toy_backbone, verbose = FALSE)
  mt <- import_taxonomic_table(missing_path(), type = "missing_taxa", phy = phylograftR::toy_backbone, verbose = FALSE)
  expect_identical(bb, phylograftR::toy_backbone_taxonomic_table)
  # All columns are read as text, whereas the dataset stores include_target_clade_stem as logical
  toy_mt <- phylograftR::toy_missing_taxa_taxonomic_table
  expect_identical(names(mt), names(toy_mt))
  expect_true(all(vapply(mt, is.character, logical(1))))
  toy_mt[] <- lapply(toy_mt, as.character)
  expect_identical(mt, toy_mt)
  # The default type is the backbone
  expect_identical(import_taxonomic_table(backbone_path(), verbose = FALSE), bb)
})

test_that("tables imported from files give the same missing_taxa_lists as the toy datasets",
{
  phy <- phylograftR::toy_backbone
  bb <- import_taxonomic_table(backbone_path(), "backbone", verbose = FALSE)
  mt <- import_taxonomic_table(missing_path(), "missing_taxa", verbose = FALSE)
  from_files <- build_missing_taxa_list_from_taxonomy(bb, mt, phy, verbose = FALSE)
  from_datasets <- build_missing_taxa_list_from_taxonomy(phylograftR::toy_backbone_taxonomic_table,
                                                         phylograftR::toy_missing_taxa_taxonomic_table, phy, verbose = FALSE)
  expect_identical(from_files, from_datasets)
})

### 2/ Reading ####

test_that("import_taxonomic_table() reads .csv, .tsv, .txt and files with a given separator",
{
  dir <- temp_dir()
  d <- bb_table()
  utils::write.csv(d, file.path(dir, "t.csv"), row.names = FALSE, na = "")
  utils::write.table(d, file.path(dir, "t.tsv"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)
  utils::write.table(d, file.path(dir, "t.txt"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)
  utils::write.table(d, file.path(dir, "t.dat"), sep = "|", row.names = FALSE, na = "", quote = FALSE)
  for (f in c("t.csv", "t.tsv", "t.txt"))
  {
    x <- import_taxonomic_table(file.path(dir, f), verbose = FALSE)
    expect_identical(x, d, info = f)
  }
  expect_identical(import_taxonomic_table(file.path(dir, "t.dat"), sep = "|", verbose = FALSE), d)
  expect_error(import_taxonomic_table(file.path(dir, "t.dat"), verbose = FALSE), "Unknown file extension '.dat'")
  # Blank cells are NA, and columns are text
  x <- import_taxonomic_table(file.path(dir, "t.csv"), verbose = FALSE)
  expect_true(all(vapply(x, is.character, logical(1))))
  expect_true(is.na(x$Family[3]))
})

test_that("import_taxonomic_table() reads the missing taxa table with its constraints written inside quoted cells",
{
  d <- mt_table()
  d$Excluded_subclades <- c("Gen_1; Gen_2", NA)
  d$allowed_states <- c("blue, green", NA)
  path <- file.path(temp_dir(), "m.csv")
  utils::write.csv(d, path, row.names = FALSE, na = "")
  x <- import_taxonomic_table(path, type = "missing_taxa", verbose = FALSE)
  expect_identical(x, d)
  expect_identical(x$allowed_states[1], "blue, green")
})

test_that("import_taxonomic_table() strips spaces, drops empty rows, and does not change the case of the values",
{
  d <- data.frame(` Taxa ` = c(" t1", "t2 ", NA, NA), Family = c("F1 ", " f2", NA, ""), check.names = FALSE, stringsAsFactors = FALSE)
  path <- file.path(temp_dir(), "t.csv")
  utils::write.csv(d, path, row.names = FALSE, na = "")
  x <- import_taxonomic_table(path, verbose = FALSE)
  expect_identical(names(x), c("Taxa", "Family"))
  expect_identical(x$Taxa, c("t1", "t2"))
  expect_identical(x$Family, c("F1", "f2"))
  expect_identical(rownames(x), c("1", "2"))
  # Case of Type is kept: it is ignored later
  d <- mt_table()
  path <- file.path(temp_dir(), "m.csv")
  utils::write.csv(d, path, row.names = FALSE, na = "")
  expect_identical(import_taxonomic_table(path, "missing_taxa", verbose = FALSE)$Type, c("taxa", "Clades"))
})

test_that("import_taxonomic_table() reads Excel files",
{
  skip_if_not_installed("readxl")
  skip_if_not_installed("writexl")
  path <- file.path(temp_dir(), "t.xlsx")
  writexl::write_xlsx(list(backbone = bb_table(), missing = mt_table()), path)
  expect_identical(import_taxonomic_table(path, "backbone", sheet = "backbone", verbose = FALSE), bb_table())
  expect_identical(import_taxonomic_table(path, "missing_taxa", sheet = 2, verbose = FALSE), mt_table())
})

test_that("import_taxonomic_table() accepts a user-supplied reader",
{
  d <- bb_table()
  path <- file.path(temp_dir(), "t.rds")
  saveRDS(d, path)
  expect_identical(import_taxonomic_table(path, reader = readRDS, verbose = FALSE), d)
  reader <- function (file, keep) { x <- readRDS(file); x[keep, ] }
  expect_identical(import_taxonomic_table(path, reader = reader, keep = 1:2, verbose = FALSE)$Taxa, c("t1", "t2"))
  expect_error(import_taxonomic_table(path, reader = "readRDS"), "'reader' must be a function")
  expect_error(import_taxonomic_table(path, reader = function (file) "a"), "not a data frame")
})

### 3/ Checks ####

test_that("import_taxonomic_table() checks the file and the arguments",
{
  expect_error(import_taxonomic_table(file.path(tempdir(), "nothing.csv")), "File not found")
  expect_error(import_taxonomic_table(c("a.csv", "b.csv")), "single path")
  expect_error(import_taxonomic_table(backbone_path(), sep = ";;"), "single character")
  expect_error(import_taxonomic_table(backbone_path(), type = "other"), "'arg' should be one of")
  expect_error(import_taxonomic_table(backbone_path(), phy = "a"), "'phy' must be NULL or a phylogeny")
  expect_error(import_taxonomic_table(backbone_path(), verbose = NA), "'verbose' must be TRUE or FALSE")
  path <- file.path(temp_dir(), "empty.csv")
  writeLines("Taxa,Family", path)
  expect_error(import_taxonomic_table(path), "no rows")
  writeLines(c("Taxa,Family", ",", ","), path)
  expect_error(import_taxonomic_table(path), "no rows with values")
})

test_that("import_taxonomic_table() checks the columns and the names of the taxa",
{
  dir <- temp_dir()
  write_table <- function (d, name = "t.csv") { p <- file.path(dir, name); utils::write.csv(d, p, row.names = FALSE, na = ""); p }
  expect_error(import_taxonomic_table(write_table(data.frame(Tip = "t1", Family = "F1"))), "must have a column 'Taxa'")
  expect_error(import_taxonomic_table(write_table(data.frame(Taxa = "t1"))), "at least one taxonomic level")
  expect_error(import_taxonomic_table(write_table(data.frame(Taxa = c("t1", "t1"), Family = "F1"))), "'Taxa' must be unique.*t1")
  expect_error(import_taxonomic_table(write_table(data.frame(Taxa = c("t1", NA), Family = c("F1", "F2")))), "'Taxa' is empty in the row\\(s\\).*: 2")
  expect_error(import_taxonomic_table(write_table(data.frame(Taxa = "X1", Family = "F1")), type = "missing_taxa"), "must have a column 'Type'")
  d <- mt_table(); d$Type[2] <- "family"
  expect_error(import_taxonomic_table(write_table(d), type = "missing_taxa"), "'Type' must be 'taxa' or 'clades'.*X2")
  d <- mt_table(); d$Type[1] <- NA
  expect_error(import_taxonomic_table(write_table(d), type = "missing_taxa"), "Invalid for: X1")
  # Two columns with the same name
  p <- file.path(dir, "dup.csv"); writeLines(c("Taxa,Family,Family", "t1,F1,F2"), p)
  expect_error(import_taxonomic_table(p), "same name: Family")
})

test_that("import_taxonomic_table() checks the table against the tips of the phylogeny, when provided",
{
  phy <- phylograftR::toy_backbone
  dir <- temp_dir()
  d <- bb_table(); d$Taxa[3] <- "not_a_tip"
  p <- file.path(dir, "b.csv"); utils::write.csv(d, p, row.names = FALSE, na = "")
  expect_error(import_taxonomic_table(p, phy = phy), "1 taxa of the backbone table are not tips of 'phy': not_a_tip")
  expect_no_error(import_taxonomic_table(p, verbose = FALSE))
  # Taxa to graft cannot be tips
  d <- mt_table(); d$Taxa[1] <- "t1"
  p <- file.path(dir, "m.csv"); utils::write.csv(d, p, row.names = FALSE, na = "")
  expect_error(import_taxonomic_table(p, "missing_taxa", phy = phy), "must be absent from the backbone.*t1")
  expect_no_error(import_taxonomic_table(p, "missing_taxa", verbose = FALSE))
})

### 4/ Messages ####

test_that("import_taxonomic_table() summarizes the table if verbose = TRUE",
{
  expect_message(import_taxonomic_table(backbone_path()), "Backbone taxonomic table: 30 tips, with")
  expect_message(import_taxonomic_table(missing_path(), "missing_taxa"), "Missing taxa taxonomic table: 10 taxa to graft \\(7 as 'taxa', 3 as 'clades'\\)")
  expect_no_message(import_taxonomic_table(backbone_path(), verbose = FALSE))
})
