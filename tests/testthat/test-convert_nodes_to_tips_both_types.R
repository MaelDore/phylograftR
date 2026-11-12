##### Test suite for convert_nodes_to_tips() with list(taxa, clades) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Toy backbone of the package (30 tips t1-t30). Nodes: 32 = clade A (17 tips) ; 34 = clade B1 (9 tips) ; 50 = group_3 (nested in B1)

both_nodes <- function ()
{
  list(taxa = list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = 34, X2 = c("t5", "t12")),
                   Excluded_subclades = list(X1 = list(group_3 = 50), X2 = NULL)),
       clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = "34", M2 = 32)))
}
quiet <- function (expr) suppressMessages(expr)


test_that("convert_nodes_to_tips() converts both types and returns a list with the same structure",
{
  phy <- phylograftR::toy_backbone
  inp <- both_nodes()
  res <- quiet(convert_nodes_to_tips(inp, phy))
  expect_identical(names(res), c("taxa", "clades"))
  expect_identical(names(res$taxa), names(inp$taxa)); expect_identical(names(res$clades), names(inp$clades))
  expect_true(all(vapply(c(res$taxa$Target_clade, res$clades$Sister_clade), is.character, logical(1))))
  expect_equal(unname(get_clade_node(phy, res$clades$Sister_clade$M2)), 32)
  # Each type is the same as converting it alone
  expect_identical(res$taxa$Taxon_name, inp$taxa$Taxon_name)
  expect_identical(res$taxa$Target_clade$X2, c("t5", "t12"))
})

test_that("the same clade is defined by the same tips in both types, whatever the order of definition",
{
  phy <- phylograftR::toy_backbone
  res <- quiet(convert_nodes_to_tips(both_nodes(), phy))
  expect_identical(res$taxa$Target_clade$X1, res$clades$Sister_clade$M1)
  # A clade already defined by hand by other tips in one type is used for the node of the other type
  inp <- both_nodes(); inp$clades$Sister_clade$M1 <- c("t1", "t2")
  node12 <- get_clade_node(phy, c("t1", "t2"))
  inp$taxa$Target_clade$X1 <- as.character(node12)
  res2 <- quiet(convert_nodes_to_tips(inp, phy))
  expect_identical(res2$taxa$Target_clade$X1, c("t1", "t2"))
})

test_that("the total richness of both types is converted with the same tips",
{
  phy <- phylograftR::toy_backbone
  richness <- list(taxa = list(Target_clade = list(34, c("t5", "t12")), total_richness = c(20, 8)),
                   clades = list(Sister_clade = list("34", 32), total_richness = c(15, 40)))
  res <- quiet(convert_nodes_to_tips(both_nodes(), phy, total_richness_per_clades = richness))
  expect_identical(names(res), c("missing_taxa_list", "total_richness_per_clades"))
  expect_identical(names(res$missing_taxa_list), c("taxa", "clades"))
  expect_identical(names(res$total_richness_per_clades), c("taxa", "clades"))
  expect_identical(res$total_richness_per_clades$taxa$Target_clade[[1]], res$missing_taxa_list$taxa$Target_clade$X1)
  expect_identical(res$total_richness_per_clades$clades$Sister_clade[[1]], res$missing_taxa_list$clades$Sister_clade$M1)
  expect_identical(res$total_richness_per_clades$clades$total_richness, c(15, 40))
  # One type of richness only
  one <- quiet(convert_nodes_to_tips(both_nodes(), phy, total_richness_per_clades = list(taxa = richness$taxa, clades = NULL)))
  expect_identical(names(one$total_richness_per_clades), "taxa")
})

test_that("a list with one non-NULL element is handled as that single type",
{
  phy <- phylograftR::toy_backbone
  inp <- both_nodes()
  a <- quiet(convert_nodes_to_tips(inp$taxa, phy))
  b <- quiet(convert_nodes_to_tips(list(taxa = inp$taxa, clades = NULL), phy))
  expect_identical(a, b)
})

test_that("invalid inputs with list(taxa, clades) stop with explicit messages",
{
  phy <- phylograftR::toy_backbone
  inp <- both_nodes()
  expect_error(convert_nodes_to_tips(list(taxa = NULL, clades = NULL), phy), "both are NULL")
  bad <- inp; bad$clades <- inp$taxa
  expect_error(convert_nodes_to_tips(bad, phy), "Sister_clade")
  expect_error(convert_nodes_to_tips(inp, phy, total_richness_per_clades = list(Target_clade = list(34), total_richness = 10)), "list\\(taxa")
  expect_error(convert_nodes_to_tips(list(taxa = inp$taxa, clades = NULL), phy,
                                     total_richness_per_clades = list(clades = list(Sister_clade = list(34), total_richness = 10))), "not in 'missing_taxa_list'")
  # An error in one type (mixed node IDs and tips) is reported
  bad2 <- inp; bad2$clades$Sister_clade$M1 <- c("t1", "34")
  expect_error(convert_nodes_to_tips(bad2, phy, verbose = FALSE), "mixes node IDs and tips")
})

test_that("the conversions of both types are printed together",
{
  phy <- phylograftR::toy_backbone
  msgs <- testthat::capture_messages(convert_nodes_to_tips(both_nodes(), phy))
  txt <- paste(msgs, collapse = "\n")
  expect_match(txt, "Target_clade"); expect_match(txt, "Sister_clade")
  expect_message(convert_nodes_to_tips(list(taxa = list(Taxon_name = "X1", Target_clade = list(X1 = c("t1", "t2"))),
                                           clades = list(Taxon_name = "M1", Sister_clade = list(M1 = c("t5", "t12")))), phy), "No node ID found")
})

test_that("the converted lists of both types can be used by the wrapper",
{
  toy <- toy_setup()
  node <- ape::getMRCA(toy$phy, toy$B)
  inp <- list(taxa = list(Taxon_name = "X1", Target_clade = list(X1 = ape::getMRCA(toy$phy, toy$A))),
              clades = list(Taxon_name = "M1", Sister_clade = list(M1 = node)))
  conv <- quiet(convert_nodes_to_tips(inp, toy$phy))
  g <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, conv, method = "uniform", seed = 3)
  expect_true(all(c("X1", "M1") %in% g$tip.label))
})
