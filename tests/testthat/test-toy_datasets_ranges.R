##### Test suite for the toy datasets with biogeographic ranges #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


# ---- Toy simmap of ranges ----

## Test that the toy simmap of ranges is a valid simmap of the toy backbone, with correctly coded ranges
test_that("toy_simmap_ranges is a valid simmap of toy_backbone with four ranges coded as in BioGeoBEARS",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_simmap_ranges, package = "phylograftR")
  expect_s3_class(toy_simmap_ranges, "simmap")
  expect_identical(toy_simmap_ranges$tip.label, toy_backbone$tip.label)
  expect_identical(toy_simmap_ranges$edge, toy_backbone$edge)
  expect_equal(toy_simmap_ranges$edge.length, toy_backbone$edge.length)
  expect_no_error(validate_simmap(toy_simmap_ranges, toy_backbone))
  expect_setequal(colnames(toy_simmap_ranges$mapped.edge), c("A", "B", "AB", "BC"))
  expect_no_error(validate_range_coding(colnames(toy_simmap_ranges$mapped.edge)))
  # All the ranges are found in the tips, and the root is in "AB"
  expect_setequal(unique(as.character(toy_simmap_ranges$states)), c("A", "B", "AB", "BC"))
  expect_identical(unname(toy_simmap_ranges$node.states[1, 1]), "AB")
  # Only changes of a single area within the branches: d (gain) and e (loss), so the states are connected as a chain A - AB - B - BC
  changes <- unlist(lapply(toy_simmap_ranges$maps, function (x) paste(head(names(x), -1), tail(names(x), -1), sep = ">")))
  expect_true(all(changes %in% c("A>AB", "AB>A", "B>AB", "AB>B", "B>BC", "BC>B")))
  # The toy simmap of ranges is not the toy simmap of colors
  expect_false(identical(colnames(toy_simmap_ranges$mapped.edge), colnames(phylograftR::toy_simmap$mapped.edge)))
  # The total durations per range and branch are consistent with the maps
  expect_equal(toy_simmap_ranges$mapped.edge, t(vapply(toy_simmap_ranges$maps, function (x) vapply(colnames(toy_simmap_ranges$mapped.edge), function (r) sum(x[names(x) == r]), numeric(1)), numeric(4))),
               ignore_attr = TRUE)
})

## Test the events represented in the toy simmap of ranges, with the syntax of BioGeoBEARS
test_that("toy_simmap_ranges represents anagenetic (d, e) and cladogenetic (y, s, v) events of BioGeoBEARS, as documented",
{
  data(toy_simmap_ranges, package = "phylograftR")
  n_tips <- ape::Ntip(toy_simmap_ranges)
  areas <- function (range) strsplit(range, "")[[1]]
  range_above <- function (node)
  {
    parent_edge <- match(node, toy_simmap_ranges$edge[, 2])
    if (is.na(parent_edge)) toy_simmap_ranges$node.states[1, 1] else tail(names(toy_simmap_ranges$maps[[parent_edge]]), 1)
  }

  # Anagenetic changes within the branches: d (gain of one area) or e (loss of one area)
  changes <- do.call(rbind, lapply(toy_simmap_ranges$maps, function (x) if (length(x) > 1) data.frame(from = head(names(x), -1), to = tail(names(x), -1))))
  n_gain <- sum(mapply(function (from, to) length(areas(to)) == length(areas(from)) + 1 && all(areas(from) %in% areas(to)), changes$from, changes$to))
  n_loss <- sum(mapply(function (from, to) length(areas(to)) == length(areas(from)) - 1 && all(areas(to) %in% areas(from)), changes$from, changes$to))
  expect_identical(c(d = n_gain, e = n_loss), c(d = 8L, e = 4L))
  expect_identical(n_gain + n_loss, nrow(changes))

  # Cladogenetic events at the nodes
  type <- vapply((n_tips + 1):(n_tips + toy_simmap_ranges$Nnode), function (node)
  {
    node_range <- range_above(node)
    descendants <- vapply(which(toy_simmap_ranges$edge[, 1] == node), function (edge) names(toy_simmap_ranges$maps[[edge]])[1], character(1))
    if (all(descendants == node_range)) { return(if (nchar(node_range) == 1) "y" else "range_copying") }
    if (sum(descendants == node_range) == 1) { return(if (nchar(node_range) > 1 && nchar(descendants[descendants != node_range]) == 1 && descendants[descendants != node_range] %in% areas(node_range)) "s" else "other") }
    if (!any(areas(descendants[1]) %in% areas(descendants[2])) && setequal(unlist(lapply(descendants, areas)), areas(node_range))) "v" else "other"
  }, character(1))
  expect_identical(as.list(table(type)[c("y", "range_copying", "s", "v")]), list(y = 13L, range_copying = 10L, s = 5L, v = 1L))
  expect_false("other" %in% type)
  # Cladogenetic events are branches starting in a range that differs from the range at the end of the parental branch
  expect_identical(sum(type %in% c("s", "v")), 6L)
  expect_gt(sum(type %in% c("s", "v")) / (sum(type %in% c("s", "v")) + nrow(changes)), 0.25)
  expect_lt(sum(type %in% c("s", "v")) / (sum(type %in% c("s", "v")) + nrow(changes)), 0.36)
})


# ---- Toy list with allowed ranges ----

## Test that the list is the import of the .csv file, with the same taxa as the toy list and correct ranges
test_that("toy_missing_taxa_list_with_allowed_ranges is the import of its .csv file and has valid ranges for toy_simmap_ranges",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_missing_taxa_list, package = "phylograftR")
  data(toy_simmap_ranges, package = "phylograftR")
  data(toy_missing_taxa_list_with_allowed_ranges, package = "phylograftR")
  csv <- system.file("extdata", "toy_missing_taxa_with_allowed_ranges.csv", package = "phylograftR")
  expect_true(nzchar(csv))
  expect_identical(import_missing_taxa_list(file = csv, phy = toy_backbone), toy_missing_taxa_list_with_allowed_ranges)

  # Same items as the toy list, with allowed ranges
  expect_identical(toy_missing_taxa_list_with_allowed_ranges[names(toy_missing_taxa_list)], toy_missing_taxa_list)
  allowed <- toy_missing_taxa_list_with_allowed_ranges$allowed_states
  expect_identical(names(allowed), toy_missing_taxa_list$Taxon_name)
  expect_identical(allowed[lengths(allowed) > 0],
                   list(New_sp1 = "A", New_sp4 = c("B", "BC"), New_sp5 = "B", New_sp7 = c("A", "BC")))
  # Ranges are states of the simmap, so nothing is ignored and no warning is raised
  expect_true(all(unlist(allowed) %in% colnames(toy_simmap_ranges$mapped.edge)))
  expect_no_warning(resolve_allowed_ranges(allowed_states_taxa = allowed, allowed_states = c("A", "B"), states = colnames(toy_simmap_ranges$mapped.edge),
                                           taxon_names = toy_missing_taxa_list$Taxon_name))
})

## Test that grafting with the toy datasets works without warning, respects the constraints, and includes subset sympatric speciation events
test_that("grafting the toy list with allowed ranges on the toy simmap of ranges is valid and gives subset sympatric speciation events",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_simmap_ranges, package = "phylograftR")
  data(toy_missing_taxa_list_with_allowed_ranges, package = "phylograftR")
  mtl <- toy_missing_taxa_list_with_allowed_ranges
  taxa <- mtl$Taxon_name

  for (seed in c(1234, 1:3))
  {
    expect_no_warning(grafted <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = mtl, simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                                                    treat_states_as_ranges = TRUE, verbose = FALSE, seed = seed))
    checks <- check_grafting_constraints(grafted, toy_backbone, mtl, allowed_states = c("A", "B"), treat_states_as_ranges = TRUE, verbose = FALSE)
    expect_true(attr(checks, "all_ok"))
    expect_true(all(checks$ranges_ok, na.rm = TRUE))
    # The ranges inherited are the allowed ranges of the taxa (global ranges for the taxa without taxon-level ones)
    inherited <- graft_states_in_simmap(grafted, taxa)
    expected <- Map(function (taxon) if (is.null(mtl$allowed_states[[taxon]])) c("A", "B") else mtl$allowed_states[[taxon]], taxa)
    expect_true(all(mapply(function (range, taxon) range %in% expected[[taxon]], inherited$state, inherited$taxon)))
  }

  # With the seed of the examples, a lineage inherits a range smaller than the state of its parent branch (subset sympatric speciation)
  grafted <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = mtl, simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                                treat_states_as_ranges = TRUE, verbose = FALSE, seed = 1234)
  inherited <- graft_states_in_simmap(grafted, taxa)
  expect_true(any(inherited$state != inherited$parent_state))
  expect_true(all(range_in_states_pairwise(inherited$state, inherited$parent_state)))

  # Matching the states exactly, the ranges "A" and "B" are not allowed on all the branches: a list of ranges requires treat_states_as_ranges = TRUE to use ranges as subsets
  expect_error(graft_missing_taxa(phy = toy_backbone, missing_taxa_list = mtl, allowed_states = c("A", "B"), verbose = FALSE, seed = 1234), "simmap")
})
