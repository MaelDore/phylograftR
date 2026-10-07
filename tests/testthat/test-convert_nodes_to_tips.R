
##### Test suite for convert_nodes_to_tips() and the node IDs of build_missing_taxa_list() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Toy backbone of the package (30 tips t1-t30, tips are nodes 1-30, the root is node 31). Nodes used below:
# 32 = clade A (17 tips) ; 34 = clade B1 (9 tips) ; 37 = group_1 (8 tips, nested in A) ; 39 = group_2 (6 tips, nested in A)
# 50 = group_3 (t13, t17, t27; nested in B1) ; 53 = clade C (t5, t23, t12; nested in B) ; node 5 = tip t11

# Build a list with a mix of node IDs and tips for testing
mixed_list <- function ()
{
  list(Taxon_name = c("X1", "X2", "X3"),
       Target_clade = list(X1 = 34, X2 = c("t5", "t12"), X3 = "32"),
       include_target_clade_stem = c(X1 = TRUE, X2 = NA, X3 = FALSE),
       Excluded_subclades = list(X1 = list(group_3 = 50, tip_t7 = "t7"), X2 = NULL, X3 = list(37, c("t26", "t4"))),
       exclude_subclades_stem = list(X1 = c(group_3 = TRUE, tip_t7 = FALSE), X2 = NULL, X3 = TRUE))
}

quiet <- function (expr) suppressMessages(expr)


### 1/ Conversion of node IDs ####

test_that("convert_nodes_to_tips() converts node IDs, text of digits and tips into tips, keeping the other items unchanged",
{
  phy <- phylograftR::toy_backbone
  mtl <- mixed_list()
  res <- quiet(convert_nodes_to_tips(mtl, phy))
  expect_identical(names(res), names(mtl))
  expect_identical(res$Taxon_name, mtl$Taxon_name)
  expect_identical(names(res$Target_clade), c("X1", "X2", "X3"))
  expect_identical(res$Target_clade$X2, c("t5", "t12"))
  # The converted clades are the same clades as the node IDs
  nodes <- vapply(res$Target_clade, function (tips) get_clade_node(phy, tips), numeric(1))
  expect_identical(unname(nodes), c(34, ape::getMRCA(phy, c("t5", "t12")), 32))
  expect_true(all(vapply(res$Target_clade, is.character, logical(1))))
  # Excluded subclades, with their inner names
  expect_identical(names(res$Excluded_subclades$X1), c("group_3", "tip_t7"))
  expect_identical(res$Excluded_subclades$X1$tip_t7, "t7")
  expect_equal(get_clade_node(phy, res$Excluded_subclades$X1$group_3), 50)
  expect_null(res$Excluded_subclades$X2)
  expect_identical(names(res$Excluded_subclades), c("X1", "X2", "X3"))
  expect_identical(vapply(res$Excluded_subclades$X3, function (tips) get_clade_node(phy, tips), numeric(1)),
                   c(37, ape::getMRCA(phy, c("t26", "t4"))))
  # Other items unchanged
  expect_identical(res$include_target_clade_stem, mtl$include_target_clade_stem)
  expect_identical(res$exclude_subclades_stem, mtl$exclude_subclades_stem)
})

test_that("convert_nodes_to_tips() is consistent in converting a tip node into its tip label, and an internal node into one pair of tips corresponding to the MRCA, deterministically",
{
  phy <- phylograftR::toy_backbone
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = 5, X2 = 32))
  res <- quiet(convert_nodes_to_tips(mtl, phy))
  expect_identical(res$Target_clade$X1, phy$tip.label[5])
  expect_length(res$Target_clade$X2, 2)
  expect_equal(get_clade_node(phy, res$Target_clade$X2), 32)
  expect_identical(quiet(convert_nodes_to_tips(mtl, phy)), res)
  # All the internal nodes, including the root, are converted into a pair with the right MRCA
  for (node in seq(ape::Ntip(phy) + 1, ape::Ntip(phy) + phy$Nnode))
  {
    tips <- node_to_tips(phy, node)
    expect_length(tips, 2)
    expect_identical(as.numeric(ape::getMRCA(phy, tips)), as.numeric(node))
  }
})

test_that("node_to_tips() conversion is consistent with a polytomy",
{
  phy <- ape::read.tree(text = "((a1:1,a2:1,a3:1):1,(b1:1,b2:1):1);")
  tips <- node_to_tips(phy, 6)
  expect_identical(as.numeric(ape::getMRCA(phy, tips)), 6)
  expect_length(tips, 2)
})

test_that("convert_nodes_to_tips() produces a converted list that is usable in the grafting functions",
{
  phy <- phylograftR::toy_backbone
  res <- quiet(convert_nodes_to_tips(mixed_list(), phy))
  grafted <- quiet_graft(graft_missing_taxa, phy, res, method = "uniform", seed = 3)
  expect_true(all(c("X1", "X2", "X3") %in% grafted$tip.label))
  check <- check_grafting_constraints(grafted, phy, res, verbose = FALSE)
  expect_true(all(check$ok))
})

test_that("convert_nodes_to_tips() returns unchanged a list without node IDs, with a message",
{
  phy <- phylograftR::toy_backbone
  mtl <- phylograftR::toy_missing_taxa_list
  expect_message(res <- convert_nodes_to_tips(mtl, phy), "No node ID found")
  expect_identical(res, mtl)
  expect_silent(convert_nodes_to_tips(mtl, phy, verbose = FALSE))
})

test_that("convert_nodes_to_tips() converts sister clades",
{
  phy <- phylograftR::toy_backbone
  mtl <- list(Taxon_name = c("X1", "X2"), Sister_clade = list(X1 = 53, X2 = c("t8", "t22")))
  res <- quiet(convert_nodes_to_tips(mtl, phy))
  expect_equal(get_clade_node(phy, res$Sister_clade$X1), 53)
  expect_identical(res$Sister_clade$X2, c("t8", "t22"))
})

test_that("convert_nodes_to_tips() converts total_richness_per_clades, and the result has two items",
{
  phy <- phylograftR::toy_backbone
  tr <- list(Target_clade = list(32, c("t5", "t12"), "t30"), total_richness = c(25, 8, NA))
  res <- quiet(convert_nodes_to_tips(mixed_list(), phy, total_richness_per_clades = tr))
  expect_identical(names(res), c("missing_taxa_list", "total_richness_per_clades"))
  expect_equal(get_clade_node(phy, res$total_richness_per_clades$Target_clade[[1]]), 32)
  expect_identical(res$total_richness_per_clades$Target_clade[[2]], c("t5", "t12"))
  expect_identical(res$total_richness_per_clades$total_richness, c(25, 8, NA))
  expect_identical(res$missing_taxa_list, quiet(convert_nodes_to_tips(mixed_list(), phy)))
  # With sister clades
  mtl <- list(Taxon_name = "X1", Sister_clade = list(X1 = 53))
  tr <- list(Sister_clade = list(53, 34), total_richness = c(5, 12))
  res <- quiet(convert_nodes_to_tips(mtl, phy, total_richness_per_clades = tr))
  expect_equal(get_clade_node(phy, res$total_richness_per_clades$Sister_clade[[2]]), 34)
})


### 2/ Messages ####

test_that("convert_nodes_to_tips() prints the conversions, grouped by node ID, with where they are used",
{
  phy <- phylograftR::toy_backbone
  mtl <- list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = 32, X2 = 32, X3 = "t5"),
              Excluded_subclades = list(X1 = list(5), X2 = NULL, X3 = NULL))
  expect_message(convert_nodes_to_tips(mtl, phy), "2 node ID\\(s\\) converted into tips")
  msg <- paste(capture_messages(convert_nodes_to_tips(mtl, phy)), collapse = "")
  expect_match(msg, "node 32 \\(internal node, 17 tips\\) -> t[0-9]+, t[0-9]+ +\\[Target_clade: X1, X2\\]")
  expect_match(msg, "node 5 \\(tip node\\) -> t11 +\\[Excluded_subclades: X1\\]")
  expect_silent(convert_nodes_to_tips(mtl, phy, verbose = FALSE))
  # Long lists of taxa are truncated
  taxa <- paste0("X", 1:10)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(32), 10), taxa))
  msg <- paste(capture_messages(convert_nodes_to_tips(mtl, phy)), collapse = "")
  expect_match(msg, "\\.\\.\\. \\(6 more\\)")
})


### 3/ Checks ####

test_that("convert_nodes_to_tips() checks the inputs and produce the expected informative errors",
{
  phy <- phylograftR::toy_backbone
  mtl <- mixed_list()
  expect_error(convert_nodes_to_tips(mtl, phy = "phy"), "'phy' must be")
  expect_error(convert_nodes_to_tips("a", phy), "must be a list")
  expect_error(convert_nodes_to_tips(mtl, phy, verbose = NA), "'verbose'")
  m2 <- mtl; m2$Sister_clade <- m2$Target_clade
  expect_error(convert_nodes_to_tips(m2, phy), "both 'Target_clade' and 'Sister_clade'")
  m2 <- mtl; m2$Target_clade <- NULL
  expect_error(convert_nodes_to_tips(m2, phy), "must have an item 'Target_clade'")
  # Names of items
  m2 <- mtl; names(m2$Target_clade) <- c("X1", "X3", "X2")
  expect_error(convert_nodes_to_tips(m2, phy), "same name|same order|must match")
  m2 <- mtl; m2$Target_clade$X1 <- NULL
  expect_error(convert_nodes_to_tips(m2, phy), "same length")
  # Empty clade, subclades not in a list
  m2 <- mtl; m2$Target_clade["X1"] <- list(NULL)
  expect_error(convert_nodes_to_tips(m2, phy), "is empty for: X1")
  m2 <- mtl; m2$Excluded_subclades$X1 <- c("t1", "t2")
  expect_error(convert_nodes_to_tips(m2, phy), "'X1' must be a list of subclades")
  m2 <- mtl; m2$Excluded_subclades$X1 <- list(character(0))
  expect_error(convert_nodes_to_tips(m2, phy), "empty clade")
  m2 <- mtl; m2$Target_clade$X1 <- TRUE
  expect_error(convert_nodes_to_tips(m2, phy), "tips \\(character\\) or node IDs \\(numbers\\)")
  m2 <- mtl; m2$Target_clade$X1 <- list("t1", "t2")
  expect_error(convert_nodes_to_tips(m2, phy), "defined by a list")
})

test_that("convert_nodes_to_tips() detects a clade mixing tips and node IDs, or having several node IDs, as an error",
{
  phy <- phylograftR::toy_backbone
  m2 <- mixed_list(); m2$Target_clade$X2 <- c("t5", "34")
  expect_error(convert_nodes_to_tips(m2, phy), "'X2' mixes node IDs and tips")
  m2 <- mixed_list(); m2$Target_clade$X2 <- c(34, 35)
  expect_error(convert_nodes_to_tips(m2, phy), "'X2' has several node IDs")
  m2 <- mixed_list(); m2$Excluded_subclades$X1 <- list(c("t7", "50"))
  expect_error(convert_nodes_to_tips(m2, phy), "'Excluded_subclades' of 'X1' mixes")
})

test_that("convert_nodes_to_tips() detects invalid node IDs as errors",
{
  phy <- phylograftR::toy_backbone
  for (bad in list(0, 60, 999, 34.5, -3, NA_real_))
  {
    m2 <- mixed_list(); m2$Target_clade$X2 <- bad
    expect_error(convert_nodes_to_tips(m2, phy), "invalid node ID|missing value", info = as.character(bad))
  }
  # Valid limits: the first tip and the last internal node
  m2 <- mixed_list(); m2$Target_clade$X2 <- 1
  expect_error(quiet(convert_nodes_to_tips(m2, phy)), NA)
  m2$Target_clade$X2 <- 59
  expect_error(quiet(convert_nodes_to_tips(m2, phy)), NA)
})

test_that("convert_nodes_to_tips() reports tips absent from the backbone as error",
{
  phy <- phylograftR::toy_backbone
  m2 <- mixed_list(); m2$Target_clade$X2 <- c("t5", "t99")
  expect_error(convert_nodes_to_tips(m2, phy), "missing from the backbone.*t99")
  m2 <- mixed_list(); m2$Excluded_subclades$X1 <- list(50, "tZ")
  expect_error(convert_nodes_to_tips(m2, phy), "missing from the backbone.*tZ")
})

test_that("convert_nodes_to_tips() refuses node IDs when some tips of the backbone are only made of digits",
{
  phy <- phylograftR::toy_backbone
  phy$tip.label[1] <- "12"
  expect_error(convert_nodes_to_tips(mixed_list(), phy), "only made of digits.*12.*Rename these tips")
  # Without node IDs, numeric tip labels are not a problem here
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = c("t5", "t12")))
  expect_error(convert_nodes_to_tips(mtl, phy, verbose = FALSE), NA)
})


### 4/ Node IDs in tables: build_missing_taxa_list() and import_missing_taxa_list() ####

test_that("build_missing_taxa_list() converts node IDs, alone or mixed with tips, as the exported function",
{
  phy <- phylograftR::toy_backbone
  d <- data.frame(Taxon_name = c("X1", "X2", "X3"),
                  Target_clade = c("34", "t5, t12", "32"),
                  Excluded_subclades = c("50; t7", NA, "37; t26, t4"),
                  Excluded_subclades_names = c("group_3; tip_t7", NA, "group_1; group_2"),
                  exclude_subclades_stem = c("TRUE; FALSE", NA, "TRUE"),
                  stringsAsFactors = FALSE)
  res <- quiet(build_missing_taxa_list(d, phy))
  ref <- quiet(convert_nodes_to_tips(mixed_list(), phy))
  expect_identical(res$Target_clade, ref$Target_clade)
  # Same subclades, the names of the subclades being those of the table
  expect_identical(lapply(res$Excluded_subclades, unname), lapply(ref$Excluded_subclades, unname))
  expect_identical(names(res$Excluded_subclades$X3), c("group_1", "group_2"))
  # Equivalent to the table with the tips written by hand
  d_tips <- d
  d_tips$Target_clade <- c(paste(ref$Target_clade$X1, collapse = ","), "t5, t12", paste(ref$Target_clade$X3, collapse = ","))
  d_tips$Excluded_subclades <- c(paste(paste(ref$Excluded_subclades$X1$group_3, collapse = ","), "t7", sep = ";"), NA,
                                 paste(paste(ref$Excluded_subclades$X3[[1]], collapse = ","), "t26,t4", sep = ";"))
  expect_identical(quiet(build_missing_taxa_list(d_tips, phy)), res)
})

test_that("build_missing_taxa_list() prints the conversions unless verbose = FALSE",
{
  phy <- phylograftR::toy_backbone
  d <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("34", "t5,t12"), Excluded_subclades = c("50", NA), stringsAsFactors = FALSE)
  expect_message(build_missing_taxa_list(d, phy), "node 34 \\(internal node, 9 tips\\).*\\[Target_clade: X1\\]")
  expect_message(build_missing_taxa_list(d, phy), "node 50.*\\[Excluded_subclades: X1\\]")
  expect_silent(build_missing_taxa_list(d, phy, verbose = FALSE))
  # No message without node IDs
  expect_silent(build_missing_taxa_list(data.frame(Taxon_name = "X1", Target_clade = "t5,t12"), phy))
})

test_that("build_missing_taxa_list() checks the clades defined by node IDs",
{
  phy <- phylograftR::toy_backbone
  d <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("34", "t5,t12"), Excluded_subclades = c("50", NA), stringsAsFactors = FALSE)
  d2 <- d; d2$Target_clade[1] <- "t8, 34"
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "'X1' mixes node IDs and tips")
  d2 <- d; d2$Target_clade[1] <- "34, 35"
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "several node IDs")
  d2 <- d; d2$Target_clade[1] <- "999"
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "invalid node ID.*999")
  d2 <- d; d2$Excluded_subclades[1] <- "50; t7, 50"
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "mixes node IDs and tips")
  # Empty values are still reported by the checks of tips
  d2 <- d; d2$Excluded_subclades[1] <- "50;;t7"
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "empty subclade")
  d2 <- d; d2$Target_clade[1] <- "34,,t7"
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "empty tip")
  d2 <- d; d2$Target_clade[2] <- NA
  expect_error(build_missing_taxa_list(d2, phy, verbose = FALSE), "is empty for.*X2")
  # Numeric tip labels in the backbone
  phy2 <- phy; phy2$tip.label[1] <- "12"
  expect_error(build_missing_taxa_list(d, phy2, verbose = FALSE), "only made of digits")
})

test_that("the checks of the tips apply to the tips converted from node IDs",
{
  phy <- phylograftR::toy_backbone
  # Excluded subclade (node 39, in clade A) outside the target clade (node 34): warning refers to the converted tips
  d <- data.frame(Taxon_name = "X1", Target_clade = "34", Excluded_subclades = "39", stringsAsFactors = FALSE)
  expect_warning(quiet(build_missing_taxa_list(d, phy)), "outside the target clade.*X1.*target clade tips: t8, t7")
  # Excluded subclade identical to the target clade, defined by different means
  d <- data.frame(Taxon_name = "X1", Target_clade = "34", Excluded_subclades = "t8,t7", stringsAsFactors = FALSE)
  expect_error(build_missing_taxa_list(d, phy, verbose = FALSE), "identical to, or contain")
  # Labels: the same clade defined by a node ID in one row and by tips in another must have the same label
  d <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("34", "t8,t7"), Target_clade_name = c("B1", "B1"), stringsAsFactors = FALSE)
  expect_error(quiet(build_missing_taxa_list(d, phy)), NA)
  d$Target_clade_name <- c("B1", "other")
  expect_error(build_missing_taxa_list(d, phy, verbose = FALSE), "Inconsistent labels in 'Target_clade_name'")
})

test_that("build_missing_taxa_list() deals with sister clades defined by node IDs and gives a list usable by graft_missing_clades()",
{
  phy <- phylograftR::toy_backbone
  d <- data.frame(Taxon_name = c("X1", "X2"), Sister_clade = c("53", "t8,t22"), stringsAsFactors = FALSE)
  l <- quiet(build_missing_taxa_list(d, phy))
  expect_equal(get_clade_node(phy, l$Sister_clade$X1), 53)
  grafted <- quiet_graft(graft_missing_clades, phy, l, seed = 1)
  expect_true(all(c("X1", "X2") %in% grafted$tip.label))
})

test_that("import_missing_taxa_list() reads node IDs from .csv and Excel files (read as text)",
{
  phy <- phylograftR::toy_backbone
  d <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("34", "t5,t12"), Excluded_subclades = c("50", NA), stringsAsFactors = FALSE)
  expected <- quiet(build_missing_taxa_list(d, phy))
  dir <- tempfile("phylograftR_"); dir.create(dir)
  utils::write.csv(d, file.path(dir, "t.csv"), row.names = FALSE, na = "")
  expect_identical(quiet(import_missing_taxa_list(file.path(dir, "t.csv"), phy)), expected)
  expect_message(import_missing_taxa_list(file.path(dir, "t.csv"), phy), "node ID\\(s\\) converted")
  expect_silent(import_missing_taxa_list(file.path(dir, "t.csv"), phy, verbose = FALSE))
  skip_if_not_installed("readxl")
  skip_if_not_installed("writexl")
  # Node IDs written as text, or typed as numbers in Excel
  writexl::write_xlsx(data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("34", "t5,t12")), file.path(dir, "t.xlsx"))
  expect_identical(quiet(import_missing_taxa_list(file.path(dir, "t.xlsx"), phy))$Target_clade,
                   expected$Target_clade)
  d_pure <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c(34, 32))
  writexl::write_xlsx(d_pure, file.path(dir, "n.xlsx"))
  res <- quiet(import_missing_taxa_list(file.path(dir, "n.xlsx"), phy))
  expect_equal(get_clade_node(phy, res$Target_clade$X2), 32)
})


### 5/ Node IDs are refused by the grafting functions, with a pointer to convert_nodes_to_tips() ####

test_that("the grafting functions stop if clades are defined by node IDs, list some of them and point to convert_nodes_to_tips()",
{
  phy <- phylograftR::toy_backbone
  mtl <- mixed_list()
  # Target clades (number and text of digits) and excluded subclades
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE),
               "4 clade\\(s\\) of 'missing_taxa_list' are defined by node IDs, but tips are required")
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE), "'Target_clade' of 'X1': 34")
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE), "'Target_clade' of 'X3': \"32\"")
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE), "'Excluded_subclades' of 'X1' \\(subclade 'group_3'\\): 50")
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE), "'Excluded_subclades' of 'X3' \\(subclade #1\\): 37")
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE), "convert_nodes_to_tips\\(missing_taxa_list, phy\\)")
  expect_error(graft_missing_taxa(phy, mtl, method = "BD_clades", verbose = FALSE), "defined by node IDs")
  # Only the excluded subclades
  m2 <- list(Taxon_name = "X1", Target_clade = list(X1 = c("t8", "t22")), Excluded_subclades = list(X1 = list(50)))
  expect_error(graft_missing_taxa(phy, m2, method = "uniform", verbose = FALSE), "1 clade\\(s\\).*'Excluded_subclades' of 'X1'")
  # The other functions using the list
  expect_error(check_grafting_constraints(phy, phy, mtl, verbose = FALSE), "defined by node IDs")
  ms <- list(Taxon_name = c("X1", "X2"), Sister_clade = list(X1 = 53, X2 = c("t8", "t22")))
  expect_error(graft_missing_clades(phy, ms, verbose = FALSE), "'Sister_clade' of 'X1': 53")
  # The conversion fixes it
  fixed <- quiet(convert_nodes_to_tips(mtl, phy))
  expect_error(quiet_graft(graft_missing_taxa, phy, fixed, method = "uniform", seed = 1), NA)
  expect_error(quiet_graft(graft_missing_clades, phy, quiet(convert_nodes_to_tips(ms, phy)), seed = 1), NA)
})

test_that("inforamtive errors only lists some node IDs when there are many",
{
  phy <- phylograftR::toy_backbone
  taxa <- paste0("X", 1:8)
  mtl <- list(Taxon_name = taxa, Target_clade = stats::setNames(rep(list(32), 8), taxa))
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", verbose = FALSE), "8 clade\\(s\\).*\\.\\.\\. and 3 more")
})

test_that("is_node_id_clade() can distinguish between tips and a numeric-looking tip: only digits-only values are node IDs",
{
  expect_false(is_node_id_clade(c("t1", "t2")))
  expect_false(is_node_id_clade("sp_12"))
  expect_true(is_node_id_clade(12))
  expect_true(is_node_id_clade("12"))
  expect_true(is_node_id_clade(" 12 "))
  expect_false(is_node_id_clade(NULL))
})

test_that("total_richness_per_clades with node IDs is refused by the grafting functions",
{
  phy <- phylograftR::toy_backbone
  mtl <- phylograftR::toy_missing_taxa_list
  tr <- list(Target_clade = list(32, c("t8", "t22"), "34"), total_richness = c(25, NA, 14))
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", total_richness_per_clades = tr, verbose = FALSE),
               "2 clade\\(s\\) of 'total_richness_per_clades\\$Target_clade' are defined by node IDs.*clade 1: 32.*clade 3: 34.*convert_nodes_to_tips")
  expect_error(graft_missing_taxa(phy, mtl, method = "uniform", total_richness_per_clades = tr, verbose = FALSE),
               "AT THE SAME TIME as its associated 'missing_taxa_list'.*already converted.*template")
  # The conversion fixes it
  converted <- quiet(convert_nodes_to_tips(mtl, phy, total_richness_per_clades = tr))
  # The richness list must cover all the target clades of the missing taxa: use the dataset clades with a converted node
  tr2 <- phylograftR::toy_total_richness_per_clades
  tr2$Target_clade[[1]] <- 32
  fixed <- quiet(convert_nodes_to_tips(mtl, phy, total_richness_per_clades = tr2))
  expect_error(quiet_graft(graft_missing_taxa, phy, fixed$missing_taxa_list, method = "uniform", seed = 1,
                           total_richness_per_clades = fixed$total_richness_per_clades), NA)
})


### 6/ Same clade, same tips in the missing_taxa_list and total_richness_per_clades ####

test_that("convert_nodes_to_tips() converts a node ID into the tips already used for the same clade (by hand or converted), in the list and the richness",
{
  phy <- phylograftR::toy_backbone
  # Clade A (node 32) is written by hand with t26,t18 in X1: the node ID of X2 and of the richness use the same tips
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = c("t26", "t18"), X2 = 32),
              Excluded_subclades = list(X1 = list(37), X2 = NULL))
  tr <- list(Target_clade = list(32, 37, "t5"), total_richness = c(25, 10, NA))
  res <- quiet(convert_nodes_to_tips(mtl, phy, total_richness_per_clades = tr))
  expect_identical(res$missing_taxa_list$Target_clade$X2, c("t26", "t18"))
  expect_identical(res$total_richness_per_clades$Target_clade[[1]], c("t26", "t18"))
  # Node 37 is not defined by tips anywhere: default pair, the same in the list and the richness
  expect_identical(res$total_richness_per_clades$Target_clade[[2]], res$missing_taxa_list$Excluded_subclades$X1[[1]])
  expect_identical(res$total_richness_per_clades$Target_clade[[3]], "t5")
  expect_message(convert_nodes_to_tips(mtl, phy, total_richness_per_clades = tr), "tips already used in the list")
})

test_that("convert_nodes_to_tips() can use an already converted missing_taxa_list as a template for total_richness_per_clades",
{
  phy <- phylograftR::toy_backbone
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = c("t26", "t18"), X2 = 34))
  converted <- quiet(convert_nodes_to_tips(mtl, phy))
  tr <- list(Target_clade = list(32, 34), total_richness = c(25, 12))
  res <- quiet(convert_nodes_to_tips(converted, phy, total_richness_per_clades = tr))
  # The list is unchanged, the richness takes the tips of the list
  expect_identical(res$missing_taxa_list, converted)
  expect_identical(res$total_richness_per_clades$Target_clade[[1]], c("t26", "t18"))
  expect_identical(res$total_richness_per_clades$Target_clade[[2]], converted$Target_clade$X2)
  expect_identical(res$total_richness_per_clades$total_richness, c(25, 12))
  # Converting the two at the same time or in two steps gives the same tips
  together <- quiet(convert_nodes_to_tips(mtl, phy, total_richness_per_clades = tr))
  expect_identical(together$total_richness_per_clades, res$total_richness_per_clades)
  expect_identical(together$missing_taxa_list, converted)
})

test_that("build_missing_taxa_list() also uses the tips already written for the same clade",
{
  phy <- phylograftR::toy_backbone
  d <- data.frame(Taxon_name = c("X1", "X2"), Target_clade = c("t26, t18", "32"), stringsAsFactors = FALSE)
  res <- quiet(build_missing_taxa_list(d, phy))
  expect_identical(res$Target_clade$X2, c("t26", "t18"))
})
