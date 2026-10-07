##### Test suite for import_missing_taxa_list() and build_missing_taxa_list() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Toy backbone of the package (30 tips t1-t30). Clades used below:
# t26,t18 = clade A (17 tips) ; t10,t19 = group_1 and t26,t4 = group_2, disjoint and nested in A
# t8,t22 = clade B (13 tips) ; t8,t7 = B1 (9 tips) nested in B ; t13,t27 = group_3 nested in B1 ; t5,t12 = clade C nested in B

# Function to load path to the toy_missing_taxa.csv file
csv_path <- function () system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")

## Minimal table, with the expected list for comparison
toy_table <- function ()
{
  data.frame(Taxon_name = c("X1", "X2", "X3"),
             Target_clade = c("t26,t18", "t26,t18", "t8,t22"),
             Excluded_subclades = c("t10,t19;t26,t4", "t10,t19", NA),
             Excluded_subclades_names = c("group_1;group_2", "group_1", NA),
             exclude_subclades_stem = c("TRUE;FALSE", NA, NA),
             stringsAsFactors = FALSE)
}

# Use a new empty directory in the temporary directory of the session (removed at its end)
temp_dir <- function () { d <- tempfile("phylograftR_"); dir.create(d); d }

# Perform the silent conversion of a table into a missing_taxa_list (i.e., build a missing_taxa_list)
build_quiet <- function (data, phy = phylograftR::toy_backbone) suppressWarnings(build_missing_taxa_list(data, phy))


### 1/ Datasets ####

test_that("import_missing_taxa_list() transforms the .csv of the package into the expected toy_missing_taxa_list dataset",
{
  expect_true(file.exists(csv_path()))
  imported <- import_missing_taxa_list(csv_path(), phy = phylograftR::toy_backbone)
  expect_identical(imported, phylograftR::toy_missing_taxa_list)
})

test_that("the dataset toy_missing_taxa_list has the format required by the grafting functions",
{
  l <- phylograftR::toy_missing_taxa_list
  taxa <- c("Genus_a_sp1", "Genus_a_sp2", "Genus_a_sp3", "Genus_b_sp1", "Genus_b_sp2", "Genus_c_sp1", "Genus_c_sp2")
  expect_identical(l$Taxon_name, taxa)
  # Every item is named after the taxa
  for (item in c("Target_clade", "Target_clade_name", "include_target_clade_stem", "Excluded_subclades", "exclude_subclades_stem"))
  {
    expect_identical(names(l[[item]]), taxa, info = item)
  }
  expect_identical(l$Target_clade$Genus_c_sp1, "t5")
  expect_identical(l$Target_clade$Genus_b_sp2, c("t8", "t7"))
  # Nested and named excluded subclades
  expect_identical(l$Excluded_subclades$Genus_a_sp1, list(group_1 = c("t10", "t19"), group_2 = c("t26", "t4")))
  expect_identical(l$Excluded_subclades$Genus_a_sp2, list(group_2 = c("t26", "t4"), tip_t30 = "t30"))
  expect_null(l$Excluded_subclades$Genus_a_sp3)
  expect_identical(l$exclude_subclades_stem$Genus_a_sp1, c(group_1 = TRUE, group_2 = FALSE))
  expect_identical(l$include_target_clade_stem[["Genus_a_sp3"]], FALSE)
  expect_true(is.na(l$include_target_clade_stem[["Genus_b_sp1"]]))
  expect_silent(phylograftR:::validate_missing_taxa_list(l, clade_item = "Target_clade"))
})

test_that("the datasets can be used to graft, with the constraints respected",
{
  phy <- phylograftR::toy_backbone
  for (method in c("uniform", "BD_clades"))
  {
    grafted <- quiet_graft(graft_missing_taxa, phy, phylograftR::toy_missing_taxa_list,
                           method = method, total_richness_per_clades = phylograftR::toy_total_richness_per_clades, seed = 1)
    expect_true(all(phylograftR::toy_missing_taxa_list$Taxon_name %in% grafted$tip.label))
    res <- check_grafting_constraints(grafted, phy, phylograftR::toy_missing_taxa_list, verbose = FALSE)
    expect_true(all(res$ok), info = method)
  }
})

test_that("the .csv of total richness matches the package dataset",
{
  path <- system.file("extdata", "toy_total_richness.csv", package = "phylograftR")
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  rebuilt <- list(Target_clade = strsplit(tab$Target_clade, ",", fixed = TRUE), total_richness = tab$total_richness)
  expect_equal(rebuilt, phylograftR::toy_total_richness_per_clades)
})

test_that("the backbone dataset is an ultrametric tree with the expected size",
{
  phy <- phylograftR::toy_backbone
  expect_s3_class(phy, "phylo")
  expect_identical(ape::Ntip(phy), 30L)
  expect_true(ape::is.ultrametric(phy))
  expect_equal(max(ape::node.depth.edgelength(phy)), 50)
})


### 2/ build_missing_taxa_list(): output ####

test_that("build_missing_taxa_list() splits the tips and the subclades, and names every item after the taxa",
{
  l <- build_quiet(toy_table())
  expect_identical(l$Taxon_name, c("X1", "X2", "X3"))
  expect_identical(l$Target_clade, list(X1 = c("t26", "t18"), X2 = c("t26", "t18"), X3 = c("t8", "t22")))
  expect_identical(l$Excluded_subclades$X1, list(group_1 = c("t10", "t19"), group_2 = c("t26", "t4")))
  expect_identical(l$Excluded_subclades$X2, list(group_1 = c("t10", "t19")))
  expect_null(l$Excluded_subclades$X3)
  expect_identical(names(l$Excluded_subclades), c("X1", "X2", "X3"))
  expect_identical(l$exclude_subclades_stem$X1, c(group_1 = TRUE, group_2 = FALSE))
  expect_null(l$exclude_subclades_stem$X2)
  expect_false("Excluded_subclades_names" %in% names(l))
})

test_that("build_missing_taxa_list() trims spaces and accepts spaces around separators",
{
  d <- toy_table()
  d$Target_clade <- c(" t26 , t18 ", "t26,t18", "t8, t22")
  d$Excluded_subclades <- c("t10, t19 ; t26 , t4", "t10,t19", NA)
  d$Excluded_subclades_names <- c(" group_1 ; group_2 ", "group_1", NA)
  l <- build_quiet(d)
  expect_identical(l, build_quiet(toy_table()))
})

test_that("the optional columns can be omitted when building a missing_taxa_list from a data.frame",
{
  d <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("t26,t18", "t5"), stringsAsFactors = FALSE)
  l <- build_quiet(d)
  expect_identical(names(l), c("Taxon_name", "Target_clade"))
  expect_identical(l$Target_clade, list(X1 = c("t26", "t18"), X2 = "t5"))
  # Excluded subclades without names: no inner names
  d$Excluded_subclades <- c("t10,t19;t26,t4", NA)
  l <- build_quiet(d)
  expect_null(names(l$Excluded_subclades$X1))
  expect_length(l$Excluded_subclades$X1, 2)
})

test_that("build_missing_taxa_list() keeps a single value of exclude_subclades_stem for the taxon, TRUE/FALSE are parsed whatever the case",
{
  d <- toy_table()
  d$exclude_subclades_stem <- c("true;False", "FALSE", NA)
  l <- build_quiet(d)
  expect_identical(l$exclude_subclades_stem$X1, c(group_1 = TRUE, group_2 = FALSE))
  expect_identical(unname(l$exclude_subclades_stem$X2), FALSE)
  d <- toy_table()
  d$include_target_clade_stem <- c("TRUE", "false", "")
  l <- build_quiet(d)
  expect_identical(unname(l$include_target_clade_stem), c(TRUE, FALSE, NA))
  expect_identical(names(l$include_target_clade_stem), c("X1", "X2", "X3"))
  d$include_target_clade_stem <- c("yes", "N", "1")
  expect_identical(unname(build_quiet(d)$include_target_clade_stem), c(TRUE, FALSE, TRUE))
})

test_that("build_missing_taxa_list() accepts the sister clades and yields a list usable by graft_missing_clades()",
{
  d <- data.frame(Taxon_name = c("X1", "X2"), Sister_clade = c("t8,t22", "t5"), Sister_clade_name = c("clade_B", "tip_t5"), stringsAsFactors = FALSE)
  l <- build_quiet(d)
  expect_identical(names(l), c("Taxon_name", "Sister_clade", "Sister_clade_name"))
  phy <- phylograftR::toy_backbone
  grafted <- quiet_graft(graft_missing_clades, phy, l, seed = 1)
  expect_true(all(c("X1", "X2") %in% grafted$tip.label))
})

test_that("the list with the optional label columns grafts like the list without them",
{
  l <- build_quiet(toy_table())
  grafted <- quiet_graft(graft_missing_taxa, phylograftR::toy_backbone, l, method = "uniform", seed = 5)
  expect_true(all(l$Taxon_name %in% grafted$tip.label))
  res <- check_grafting_constraints(grafted, phylograftR::toy_backbone, l, verbose = FALSE)
  expect_true(all(res$ok))
})


### 3/ build_missing_taxa_list(): checks ####

test_that("build_missing_taxa_list() checks the inputs",
{
  d <- toy_table()
  phy <- phylograftR::toy_backbone
  expect_error(build_missing_taxa_list(d, phy = "phy"), "'phy' must be")
  expect_error(build_missing_taxa_list(as.list(d), phy), "'data' must be a data frame")
  expect_error(build_missing_taxa_list(d[0, ], phy), "no rows")
  expect_error(build_missing_taxa_list(d[, names(d) != "Taxon_name"], phy), "'Taxon_name'")
  expect_error(build_missing_taxa_list(d[, names(d) != "Target_clade"], phy), "'Target_clade'.*'Sister_clade'")
  d2 <- d; d2$Sister_clade <- d2$Target_clade
  expect_error(build_missing_taxa_list(d2, phy), "both 'Target_clade' and 'Sister_clade'")
  d2 <- d; names(d2)[2] <- "Taxon_name"
  expect_error(build_missing_taxa_list(d2, phy), "same name")
  d2 <- d[, c("Taxon_name", "Target_clade", "Excluded_subclades_names")]
  expect_error(build_missing_taxa_list(d2, phy), "needs the column 'Excluded_subclades'")
})

test_that("build_missing_taxa_list() checks the taxon names",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Taxon_name[3] <- NA
  expect_error(build_missing_taxa_list(d, phy), "'Taxon_name' is empty")
  d <- toy_table(); d$Taxon_name[3] <- "X1"
  expect_error(build_missing_taxa_list(d, phy), "unique.*X1")
  d <- toy_table(); d$Taxon_name[3] <- "t1"
  expect_error(build_missing_taxa_list(d, phy), "absent from the backbone.*t1")
})

test_that("build_missing_taxa_list() warns about unknown columns and columns not used with the chosen clade column",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Comment <- "a"
  expect_warning(build_missing_taxa_list(d, phy), "Unknown columns ignored: 'Comment'")
  l <- suppressWarnings(build_missing_taxa_list(d, phy))
  expect_false("Comment" %in% names(l))
  d <- data.frame(Taxon_name = "X1", Sister_clade = "t8,t22", include_target_clade_stem = TRUE, stringsAsFactors = FALSE)
  expect_warning(build_missing_taxa_list(d, phy), "not used with 'Sister_clade'.*include_target_clade_stem")
})

test_that("build_missing_taxa_list() checks the tips",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Target_clade[1] <- "t26,t99"
  expect_error(build_missing_taxa_list(d, phy), "missing from the backbone.*t99")
  d <- toy_table(); d$Target_clade[1] <- "t26,,t18"
  expect_error(build_missing_taxa_list(d, phy), "empty tip")
  d <- toy_table(); d$Target_clade[1] <- "t26,t26"
  expect_error(build_missing_taxa_list(d, phy), "repeated tips")
  d <- toy_table(); d$Target_clade[2] <- NA
  expect_error(build_missing_taxa_list(d, phy), "'Target_clade' is empty for.*X2")
  d <- toy_table(); d$Excluded_subclades[1] <- "t10,t19;t98"
  expect_error(build_missing_taxa_list(d, phy), "missing from the backbone.*t98")
  d <- toy_table(); d$Excluded_subclades[1] <- "t10,t19;;t26,t4"
  expect_error(build_missing_taxa_list(d, phy), "empty subclade.*';'")
  # Separators other than commas: the tip is not found
  d <- toy_table(); d$Target_clade[1] <- "t26 t18"
  expect_error(build_missing_taxa_list(d, phy), "missing from the backbone")
})

test_that("build_missing_taxa_list() checks the logical columns",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$include_target_clade_stem <- c("maybe", NA, NA)
  expect_error(build_missing_taxa_list(d, phy), "must be TRUE, FALSE or NA, not 'maybe'")
  d <- toy_table(); d$exclude_subclades_stem[1] <- "TRUE;FALSE;TRUE"
  expect_error(build_missing_taxa_list(d, phy), "3 value\\(s\\) but the taxon has 2")
  d <- toy_table(); d$exclude_subclades_stem[1] <- "TRUE;maybe"
  expect_error(build_missing_taxa_list(d, phy), "TRUE, FALSE or NA")
  d <- toy_table(); d$exclude_subclades_stem[3] <- "TRUE"
  expect_warning(build_missing_taxa_list(d, phy), "without excluded subclades.*X3")
})

test_that("build_missing_taxa_list() checks the names of the excluded subclades of a taxon",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Excluded_subclades_names[1] <- "group_1"
  expect_error(build_missing_taxa_list(d, phy), "1 name\\(s\\) but the taxon has 2")
  d <- toy_table(); d$Excluded_subclades_names[1] <- "group_1;group_1"
  expect_error(build_missing_taxa_list(d, phy), "repeated names.*group_1")
  d <- toy_table(); d$Excluded_subclades_names[3] <- "group_3"
  expect_error(build_missing_taxa_list(d, phy), "'Excluded_subclades' is empty.*X3")
})

test_that("build_missing_taxa_list() manage labels properly: the same subclade must have the same label across taxa",
{
  phy <- phylograftR::toy_backbone
  # Same clade (same MRCA, tips written differently) with the same label is fine
  d <- toy_table(); d$Excluded_subclades[2] <- "t18,t19"  # MRCA of t18,t19 (node 40) is not the MRCA of t10,t19 (node 37)
  expect_error(build_missing_taxa_list(d, phy), "Inconsistent labels in 'Excluded_subclades_names'")
  # t19,t10,t30 has the same MRCA as t10,t19 (node 37): fine with the same label
  d <- toy_table(); d$Excluded_subclades[2] <- "t19,t10,t30"
  l <- build_quiet(d)
  expect_identical(l$Excluded_subclades$X2, list(group_1 = c("t19", "t10", "t30")))
  # The same clade with different labels
  d <- toy_table(); d$Excluded_subclades_names[2] <- "group_9"
  expect_error(build_missing_taxa_list(d, phy), "same clade must always have the same label.*different labels.*group_1.*group_9")
})

test_that("build_missing_taxa_list() manage labels properly: different subclades must have different labels across taxa",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Excluded_subclades[2] <- "t26,t4"
  # 'group_1' is now used for t26,t4 in X2 and for t10,t19 in X1
  expect_error(build_missing_taxa_list(d, phy), "different clades must have different labels.*'group_1' is used for 2 different clades")
})

test_that("build_missing_taxa_list() manage labels properly: the optional labels of the target clades follow the same rules",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Target_clade_name <- c("clade_A", "clade_A", "clade_B")
  l <- build_quiet(d)
  expect_identical(unname(l$Target_clade_name), c("clade_A", "clade_A", "clade_B"))
  d$Target_clade_name <- c("clade_A", "clade_B", "clade_B")
  expect_error(build_missing_taxa_list(d, phy), "Inconsistent labels in 'Target_clade_name'")
  d$Target_clade_name <- c("clade_A", "clade_Z", "clade_B")
  expect_error(build_missing_taxa_list(d, phy), "2 different labels.*clade_A.*clade_Z")
  # Missing labels are allowed
  d$Target_clade_name <- c("clade_A", NA, NA)
  expect_error(build_missing_taxa_list(d, phy), NA)
})

test_that("build_missing_taxa_list() checks the position of the excluded subclades relative to the target clade",
{
  phy <- phylograftR::toy_backbone
  # Subclade identical to the target clade
  d <- data.frame(Taxon_name = "X1", Target_clade = "t8,t7", Excluded_subclades = "t8,t7", stringsAsFactors = FALSE)
  expect_error(build_missing_taxa_list(d, phy), "identical to, or contain, the target clade.*X1")
  # Subclade containing the target clade
  d <- data.frame(Taxon_name = "X1", Target_clade = "t8,t7", Excluded_subclades = "t8,t22", stringsAsFactors = FALSE)
  expect_error(build_missing_taxa_list(d, phy), "identical to, or contain, the target clade")
  # Subclade outside the target clade: warning only
  d <- data.frame(Taxon_name = "X1", Target_clade = "t8,t7", Excluded_subclades = "t10,t19", stringsAsFactors = FALSE)
  expect_warning(build_missing_taxa_list(d, phy), "outside the target clade.*X1")
})


### 4/ import_missing_taxa_list(): formats ####

test_that("import_missing_taxa_list() reads .csv, .tsv, .txt and files with a given separator",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table()
  expected <- build_quiet(d)
  dir <- temp_dir()

  utils::write.csv(d, file.path(dir, "t.csv"), row.names = FALSE, na = "")
  expect_identical(suppressWarnings(import_missing_taxa_list(file.path(dir, "t.csv"), phy)), expected)

  utils::write.table(d, file.path(dir, "t.tsv"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)
  expect_identical(suppressWarnings(import_missing_taxa_list(file.path(dir, "t.tsv"), phy)), expected)

  utils::write.table(d, file.path(dir, "t.txt"), sep = "\t", row.names = FALSE, na = "", quote = FALSE)
  expect_identical(suppressWarnings(import_missing_taxa_list(file.path(dir, "t.txt"), phy)), expected)

  # Semicolons cannot be the field separator here, as they are used inside the cells: use a pipe
  utils::write.table(d, file.path(dir, "t.dat"), sep = "|", row.names = FALSE, na = "", quote = FALSE)
  expect_identical(suppressWarnings(import_missing_taxa_list(file.path(dir, "t.dat"), phy, sep = "|")), expected)
  expect_error(import_missing_taxa_list(file.path(dir, "t.dat"), phy), "Unknown file extension '.dat'")

  # Upper-case extension, and a .csv whose cells are quoted
  file.copy(file.path(dir, "t.csv"), file.path(dir, "T2.CSV"))
  expect_identical(suppressWarnings(import_missing_taxa_list(file.path(dir, "T2.CSV"), phy)), expected)
})

test_that("import_missing_taxa_list() reads Excel files",
{
  skip_if_not_installed("readxl")
  skip_if_not_installed("writexl")
  phy <- phylograftR::toy_backbone
  d <- toy_table()
  path <- file.path(temp_dir(), "t.xlsx")
  writexl::write_xlsx(d, path)
  expect_identical(suppressWarnings(import_missing_taxa_list(path, phy)), build_quiet(d))
  # Several sheets
  writexl::write_xlsx(list(first = d[1, ], second = d), path)
  expect_identical(suppressWarnings(import_missing_taxa_list(path, phy, sheet = "second")), build_quiet(d))
  expect_identical(names(import_missing_taxa_list(path, phy, sheet = 1)$Target_clade), "X1")
  # Logical columns written as booleans in Excel
  d$include_target_clade_stem <- c(TRUE, FALSE, NA)
  writexl::write_xlsx(d, path)
  l <- suppressWarnings(import_missing_taxa_list(path, phy))
  expect_identical(unname(l$include_target_clade_stem), c(TRUE, FALSE, NA))
})

test_that("import_missing_taxa_list() accepts a user-supplied reader",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table()
  path <- file.path(temp_dir(), "t.rds")
  saveRDS(d, path)
  expect_identical(suppressWarnings(import_missing_taxa_list(path, phy, reader = readRDS)), build_quiet(d))
  # Extra arguments are passed to the reader
  reader <- function (file, keep) { x <- readRDS(file); x[keep, ] }
  l <- suppressWarnings(import_missing_taxa_list(path, phy, reader = reader, keep = 1:2))
  expect_identical(l$Taxon_name, c("X1", "X2"))
  expect_error(import_missing_taxa_list(path, phy, reader = "readRDS"), "'reader' must be a function")
  expect_error(import_missing_taxa_list(path, phy, reader = function (file) "a"), "not a data frame")
})

test_that("import_missing_taxa_list() checks the file",
{
  phy <- phylograftR::toy_backbone
  expect_error(import_missing_taxa_list(file.path(tempdir(), "nothing.csv"), phy), "File not found")
  expect_error(import_missing_taxa_list(c("a.csv", "b.csv"), phy), "single path")
  expect_error(import_missing_taxa_list(csv_path(), phy, sep = ";;"), "single character")
  path <- file.path(temp_dir(), "empty.csv")
  writeLines("Taxon_name,Target_clade", path)
  expect_error(import_missing_taxa_list(path, phy), "no rows")
})

test_that("import_missing_taxa_list() reports the checks of build_missing_taxa_list()",
{
  phy <- phylograftR::toy_backbone
  d <- toy_table(); d$Excluded_subclades_names[2] <- "group_9"
  path <- file.path(temp_dir(), "t.csv")
  utils::write.csv(d, path, row.names = FALSE, na = "")
  expect_error(import_missing_taxa_list(path, phy), "same clade must always have the same label")
})
