
##### Test suite for the toy datasets of two traits #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


taxa_states_in <- function (simmap, taxa)
{
  vapply(taxa, function (taxon) names(simmap$maps[[match(match(taxon, simmap$tip.label), simmap$edge[, 2])]])[1], character(1))
}


## Test that the second trait is a valid simmap of the toy backbone
test_that("toy_simmap_elevation is a valid simmap of toy_backbone with two states",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_simmap_elevation, package = "phylograftR")
  expect_s3_class(toy_simmap_elevation, "simmap")
  expect_identical(toy_simmap_elevation$tip.label, toy_backbone$tip.label)
  expect_identical(toy_simmap_elevation$edge, toy_backbone$edge)
  expect_equal(toy_simmap_elevation$edge.length, toy_backbone$edge.length)
  expect_no_error(validate_simmap(toy_simmap_elevation, toy_backbone))
  expect_setequal(colnames(toy_simmap_elevation$mapped.edge), c("high", "low"))
  # Balanced map: documented numbers
  share <- colSums(toy_simmap_elevation$mapped.edge) / sum(toy_simmap_elevation$mapped.edge)
  expect_equal(round(unname(share[c("low", "high")]), 2), c(0.54, 0.46))
  expect_equal(sum(lengths(toy_simmap_elevation$maps) - 1), 12)
  expect_equal(sum(lengths(toy_simmap_elevation$maps) > 1), 9)
  expect_lte(max(lengths(toy_simmap_elevation$maps)), 4)
  expect_equal(unname(table(factor(toy_simmap_elevation$states, levels = c("high", "low")))), unname(table(c(rep("high", 18), rep("low", 12)))), ignore_attr = TRUE)
})

## Test that the .csv and the dataset are the same, with the syntaxes of the several traits
test_that("toy_missing_taxa_list_with_two_traits is the import of its .csv, the toy list of missing taxa with allowed states of two traits",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_missing_taxa_list, package = "phylograftR")
  data(toy_missing_taxa_list_with_two_traits, package = "phylograftR")
  path <- system.file("extdata", "toy_missing_taxa_with_two_traits.csv", package = "phylograftR")
  expect_true(nzchar(path))
  imported <- import_missing_taxa_list(path, phy = toy_backbone, verbose = FALSE)
  expect_equal(imported, toy_missing_taxa_list_with_two_traits)
  # The same other items as toy_missing_taxa_list
  other_items <- setdiff(names(toy_missing_taxa_list), "allowed_states")
  expect_equal(toy_missing_taxa_list_with_two_traits[other_items], toy_missing_taxa_list[other_items])
  # The documented allowed states
  states <- toy_missing_taxa_list_with_two_traits$allowed_states
  expect_equal(states$Genus_a_sp1, list("red", NULL))
  expect_equal(states$Genus_b_sp1, list(c("blue", "green"), "high"))
  expect_equal(states$Genus_b_sp2, list(NULL, "low"))
  expect_equal(states$Genus_c_sp2, list(elevation = "high"))
  expect_null(states$Genus_a_sp2); expect_null(states$Genus_a_sp3); expect_null(states$Genus_c_sp1)
})

## Test the use of the datasets together, with both methods and without warnings
for (method in c("uniform", "BD_clades"))
{
  test_that(paste0("the toy datasets of two traits graft without warnings with ", method, ", and the taxa are in their allowed states"),
  {
    data(toy_backbone, package = "phylograftR")
    data(toy_simmap, package = "phylograftR")
    data(toy_simmap_elevation, package = "phylograftR")
    data(toy_missing_taxa_list_with_two_traits, package = "phylograftR")
    expect_no_warning(grafted <- graft_missing_taxa(toy_backbone, toy_missing_taxa_list_with_two_traits,
                                                    simmap = list(color = toy_simmap, elevation = toy_simmap_elevation),
                                                    allowed_states = list(c("blue", "green"), NULL), method = method, verbose = FALSE, seed = 1234))
    expect_equal(names(grafted$simmap), c("color", "elevation"))
    taxa <- toy_missing_taxa_list_with_two_traits$Taxon_name
    color <- taxa_states_in(grafted$simmap$color, taxa)
    elevation <- taxa_states_in(grafted$simmap$elevation, taxa)
    expect_equal(unname(color["Genus_a_sp1"]), "red")
    expect_true(color["Genus_b_sp1"] %in% c("blue", "green") && elevation["Genus_b_sp1"] == "high")
    expect_true(elevation["Genus_b_sp2"] == "low" && color["Genus_b_sp2"] %in% c("blue", "green"))
    expect_equal(unname(elevation["Genus_c_sp2"]), "high")
    # Taxa without taxon-level states: the global value for the color
    expect_true(all(color[c("Genus_a_sp2", "Genus_a_sp3", "Genus_c_sp1")] %in% c("blue", "green")))
  })
}


##### Dataset with the ranges and the elevation as traits #####

## Test that the .csv and the dataset are the same, with the syntaxes of the ranges and the elevation
test_that("toy_missing_taxa_list_with_ranges_and_elevation is the import of its .csv, the toy list of missing taxa with allowed states of ranges and elevation",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_missing_taxa_list, package = "phylograftR")
  data(toy_missing_taxa_list_with_ranges_and_elevation, package = "phylograftR")
  path <- system.file("extdata", "toy_missing_taxa_with_ranges_and_elevation.csv", package = "phylograftR")
  expect_true(nzchar(path))
  imported <- import_missing_taxa_list(path, phy = toy_backbone, verbose = FALSE)
  expect_equal(imported, toy_missing_taxa_list_with_ranges_and_elevation)
  # The same other items as toy_missing_taxa_list
  other_items <- setdiff(names(toy_missing_taxa_list), "allowed_states")
  expect_equal(toy_missing_taxa_list_with_ranges_and_elevation[other_items], toy_missing_taxa_list[other_items])
  # The documented allowed states
  states <- toy_missing_taxa_list_with_ranges_and_elevation$allowed_states
  expect_equal(states$Genus_a_sp1, list("A", "high"))
  expect_equal(states$Genus_b_sp1, list(c("B", "BC"), NULL))
  expect_equal(states$Genus_b_sp2, list(elevation = "low"))
  expect_equal(states$Genus_c_sp2, list(biogeo = c("A", "BC"), elevation = "high"))
  expect_null(states$Genus_a_sp2); expect_null(states$Genus_a_sp3); expect_null(states$Genus_c_sp1)
})

## Test the use of the datasets together, with both methods and without warnings
for (method in c("uniform", "BD_clades"))
{
  test_that(paste0("the toy datasets of ranges and elevation graft without warnings with ", method, ", respect the constraints and pass the check"),
  {
    data(toy_backbone, package = "phylograftR")
    data(toy_simmap_ranges, package = "phylograftR")
    data(toy_simmap_elevation, package = "phylograftR")
    data(toy_missing_taxa_list_with_ranges_and_elevation, package = "phylograftR")
    maps <- list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation)
    global_states <- list(c("A", "B"), NULL)
    ranges_flag <- c(biogeo = TRUE, elevation = FALSE)
    expect_no_warning(grafted <- graft_missing_taxa(toy_backbone, toy_missing_taxa_list_with_ranges_and_elevation,
                                                    simmap = maps, allowed_states = global_states, treat_states_as_ranges = ranges_flag,
                                                    method = method, verbose = FALSE, seed = 1234))
    expect_equal(names(grafted$simmap), c("biogeo", "elevation"))
    taxa <- toy_missing_taxa_list_with_ranges_and_elevation$Taxon_name
    biogeo <- taxa_states_in(grafted$simmap$biogeo, taxa)
    elevation <- taxa_states_in(grafted$simmap$elevation, taxa)
    # Taxon-level states of the documented list
    expect_equal(unname(biogeo["Genus_a_sp1"]), "A"); expect_equal(unname(elevation["Genus_a_sp1"]), "high")
    expect_true(biogeo["Genus_b_sp1"] %in% c("B", "BC"))
    expect_equal(unname(elevation["Genus_b_sp2"]), "low")
    expect_true(biogeo["Genus_c_sp2"] %in% c("A", "BC") && elevation["Genus_c_sp2"] == "high")
    # Taxa without taxon-level states: the global range ("A" or "B", possibly within a larger range of the branch)
    expect_true(all(grepl("A|B", biogeo[c("Genus_a_sp2", "Genus_a_sp3", "Genus_c_sp1")])))
    # The check of the constraints agrees, for both traits
    check <- check_grafting_constraints(grafted, toy_backbone, toy_missing_taxa_list_with_ranges_and_elevation,
                                        allowed_states = global_states, treat_states_as_ranges = ranges_flag, verbose = FALSE)
    expect_true(attr(check, "all_ok"))
    expect_true(all(check$ok))
    expect_true(all(c("states_ok_biogeo", "states_ok_elevation", "ranges_ok") %in% colnames(check)))
  })
}

## Test the use with the clades
test_that("the toy datasets of ranges and elevation can be used with graft_missing_clades() for a range trait and a normal trait",
{
  data(toy_backbone, package = "phylograftR")
  data(toy_simmap_ranges, package = "phylograftR")
  data(toy_simmap_elevation, package = "phylograftR")
  clades <- list(Taxon_name = c("M1", "M2", "M3", "M4", "M5"),
                 Sister_clade = list(M1 = c("t26", "t18"), M2 = c("t26", "t18"), M3 = c("t26", "t18"), M4 = "t14", M5 = "t14"))
  expect_no_warning(grafted <- graft_missing_clades(toy_backbone, clades,
                                                    simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
                                                    allowed_states = list(c("A", "B"), "low"),
                                                    treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
                                                    verbose = FALSE, seed = 1234))
  expect_equal(names(grafted$simmap), c("biogeo", "elevation"))
  expect_true(all(taxa_states_in(grafted$simmap$elevation, clades$Taxon_name) == "low"))
  expect_true(all(grepl("A|B", taxa_states_in(grafted$simmap$biogeo, clades$Taxon_name))))
})
