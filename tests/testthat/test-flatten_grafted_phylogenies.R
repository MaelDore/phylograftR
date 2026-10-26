
##### Test suite for flatten_grafted_phylogenies() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Create toy set for tests
flat_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}
flat_multi_simmap <- function (maps, names = NULL)
{
  names(maps) <- names
  structure(maps, class = c("multiSimmap", "multiPhylo"))
}
flat_graft <- function (toy, ..., names = c("hyp_a", "hyp_b", "hyp_c"), n_maps = 3)
{
  maps <- lapply(seq_len(n_maps), function (k) make_toy_simmap(toy$phy, seed = k))
  suppressWarnings(graft_missing_taxa(toy$phy, flat_missing_taxa(toy), simmap = flat_multi_simmap(maps, if (is.null(names)) NULL else names[seq_len(n_maps)]),
                                      allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4, ...))
}


test_that("flatten_grafted_phylogenies() flattens a list of maps with several simulations into one multiPhylo, with the names of the input",
{
  toy <- toy_setup()
  out <- flat_graft(toy, nb_simulations = 2)
  flat <- flatten_grafted_phylogenies(out)
  expect_s3_class(flat, "multiPhylo")
  expect_length(flat, 6)
  expect_equal(names(flat), c("hyp_a_sim_1", "hyp_a_sim_2", "hyp_b_sim_1", "hyp_b_sim_2", "hyp_c_sim_1", "hyp_c_sim_2"))
  # Order and content: the phylogenies are the ones of the input
  expect_identical(flat[[3]]$edge.length, out$hyp_b[[1]]$edge.length)
  expect_identical(flat[[6]]$tip.label, out$hyp_c[[2]]$tip.label)
  expect_s3_class(flat[[1]]$simmap, "simmap")
  # Attributes of each phylogeny
  expect_equal(attr(flat[[3]], "map"), "hyp_b")
  expect_equal(attr(flat[[3]], "seed"), attr(out$hyp_b[[1]], "seed"))

  table <- attr(flat, "phylogenies")
  expect_equal(names(table), c("name", "map", "simulation", "seed"))
  expect_equal(table$name, names(flat))
  expect_equal(table$map, rep(c("hyp_a", "hyp_b", "hyp_c"), each = 2))
  expect_equal(table$simulation, rep(1:2, 3))
  expect_equal(table$seed, unname(unlist(lapply(out, function (set) attr(set, "seeds")))))
  expect_equal(attr(flat, "seeds"), setNames(table$seed, names(flat)))
  expect_equal(nrow(attr(flat, "failed")), 0)
  expect_equal(names(attr(flat, "failed")), c("map", "seed", "reason"))
  expect_equal(nrow(attr(flat, "failed_maps")), 0)

  # The set works with functions that expect a multiPhylo
  heights <- sapply(flat, function (x) max(ape::node.depth.edgelength(x)))
  expect_length(heights, 6)
  expect_true(all(abs(heights - toy$root_age) < 1e-6))
})

test_that("flatten_grafted_phylogenies() flattens a list of phylogenies (one simulation per map), with the names of the maps",
{
  toy <- toy_setup()
  out <- flat_graft(toy)
  flat <- flatten_grafted_phylogenies(out)
  expect_s3_class(flat, "multiPhylo")
  expect_equal(names(flat), c("hyp_a", "hyp_b", "hyp_c"))
  expect_equal(attr(flat, "phylogenies")$simulation, c(1L, 1L, 1L))
  expect_equal(attr(flat, "phylogenies")$map, c("hyp_a", "hyp_b", "hyp_c"))
  # Without names in the multiSimmap: map_1 ...
  flat_unnamed <- flatten_grafted_phylogenies(flat_graft(toy, names = NULL, n_maps = 2))
  expect_equal(names(flat_unnamed), c("map_1", "map_2"))
})

test_that("flatten_grafted_phylogenies() selects maps and simulations, and removes the simmaps if requested",
{
  toy <- toy_setup()
  out <- flat_graft(toy, nb_simulations = 3)
  flat <- flatten_grafted_phylogenies(out, maps = c("hyp_c", "hyp_a"), simulations = c(1, 3), keep_simmap = FALSE)
  # Order of the input kept (a, then c), not the order of 'maps'
  expect_equal(names(flat), c("hyp_a_sim_1", "hyp_a_sim_3", "hyp_c_sim_1", "hyp_c_sim_3"))
  expect_true(all(vapply(flat, function (x) is.null(x$simmap), logical(1))))
  expect_equal(attr(flat, "phylogenies")$map, c("hyp_a", "hyp_a", "hyp_c", "hyp_c"))
  expect_equal(attr(flat, "phylogenies")$simulation, c(1L, 3L, 1L, 3L))
  expect_equal(names(attr(flat, "seeds")), names(flat))
  # The input is not modified
  expect_s3_class(out$hyp_a[[1]]$simmap, "simmap")
  # Only maps, or only simulations
  expect_length(flatten_grafted_phylogenies(out, maps = "hyp_b"), 3)
  expect_length(flatten_grafted_phylogenies(out, simulations = 2), 3)

  # Errors
  expect_error(flatten_grafted_phylogenies(out, maps = c("hyp_a", "unknown")), "not found: unknown.*Available maps: hyp_a, hyp_b, hyp_c")
  expect_error(flatten_grafted_phylogenies(out, simulations = 4), "No phylogeny left.*Simulations available: 1, 2, 3")
  expect_error(flatten_grafted_phylogenies(out, maps = 1), "'maps' must be NULL or a character vector")
  expect_error(flatten_grafted_phylogenies(out, simulations = 0), "'simulations' must be NULL or a vector of positive integers")
  expect_error(flatten_grafted_phylogenies(out, simulations = 1.5), "'simulations' must be NULL or a vector of positive integers")
  expect_error(flatten_grafted_phylogenies(out, keep_simmap = NA), "'keep_simmap' must be TRUE or FALSE")
  expect_error(flatten_grafted_phylogenies("a"), "must be a phylogeny")
})

test_that("flatten_grafted_phylogenies() merges the failed simulations and the failed maps, with their map",
{
  toy <- toy_setup()
  taxa <- paste0("Y", 1:2)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 2), taxa))
  maps <- list(make_young_state_simmap(toy$phy, 1e-4), make_young_state_simmap(toy$phy, 1), make_young_state_simmap(toy$phy, 40))
  out <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = flat_multi_simmap(maps, c("tiny", "narrow", "wide")), allowed_states = "B", method = "BD_clades",
                                             max_age_redraws = 1, nb_simulations = 6, verbose = FALSE, seed = 1))
  expect_false("tiny" %in% names(out))
  flat <- flatten_grafted_phylogenies(out)
  expect_equal(attr(flat, "failed_maps")$map, "tiny")
  failed <- attr(flat, "failed")
  expect_gt(nrow(failed), 0)
  expect_true(all(failed$map %in% c("narrow", "wide")))
  expect_equal(nrow(failed), sum(vapply(out, function (set) nrow(attr(set, "failed")), integer(1))))
  expect_equal(length(flat), sum(lengths(out)))
  # The names keep the original index of the simulation (gaps if a simulation failed)
  expect_equal(attr(flat, "phylogenies")$simulation, as.integer(sub(".*_sim_", "", names(flat))))
  # Selecting a map keeps its failures only
  only_narrow <- flatten_grafted_phylogenies(out, maps = "narrow")
  expect_true(all(attr(only_narrow, "failed")$map == "narrow"))
  expect_equal(nrow(attr(only_narrow, "failed_maps")), 0)
  expect_error(flatten_grafted_phylogenies(out, maps = "tiny"), "not found: tiny")
})

test_that("flatten_grafted_phylogenies() accepts a multiPhylo or a phylo (already flat), and the output without map",
{
  toy <- toy_setup()
  mtl <- flat_missing_taxa(toy)
  single <- graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, seed = 4)
  multi <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 3, verbose = FALSE, seed = 4)

  # A phylo stays a phylo: unchanged by default
  flat1 <- flatten_grafted_phylogenies(single)
  expect_s3_class(flat1, "phylo")
  expect_false(inherits(flat1, "multiPhylo"))
  expect_identical(flat1, single)
  expect_equal(attr(flat1, "seed"), 4L)
  # The only possible action is to remove the simmap
  with_simmap <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = make_toy_simmap(toy$phy), allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4))
  expect_s3_class(with_simmap$simmap, "simmap")
  expect_s3_class(flatten_grafted_phylogenies(with_simmap)$simmap, "simmap")
  no_simmap <- flatten_grafted_phylogenies(with_simmap, keep_simmap = FALSE)
  expect_s3_class(no_simmap, "phylo")
  expect_null(no_simmap$simmap)
  expect_equal(no_simmap$edge.length, with_simmap$edge.length)
  expect_equal(attr(no_simmap, "seed"), 4L)
  expect_s3_class(with_simmap$simmap, "simmap") # the input is not modified
  # Selections that cannot apply to a phylogeny without a map
  expect_error(flatten_grafted_phylogenies(single, maps = "map_1"), "'maps' cannot be used")
  expect_error(flatten_grafted_phylogenies(single, simulations = 2), "No phylogeny left")

  flat3 <- flatten_grafted_phylogenies(multi)
  expect_equal(names(flat3), c("sim_1", "sim_2", "sim_3"))
  expect_true(all(is.na(attr(flat3, "phylogenies")$map)))
  expect_equal(attr(flat3, "seeds"), attr(multi, "seeds"))
  expect_equal(flatten_grafted_phylogenies(multi, simulations = c(1, 3))[[2]]$edge.length, multi[[3]]$edge.length)
  expect_error(flatten_grafted_phylogenies(multi, maps = "map_1"), "'maps' cannot be used")

  # Simulations with gaps are kept as in the input
  gap <- multi[c(1, 3)]
  expect_equal(attr(flatten_grafted_phylogenies(gap), "phylogenies")$simulation, c(1L, 3L))

  # seed = NULL: no seeds
  no_seed <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 2, verbose = FALSE, seed = NULL)
  expect_null(attr(flatten_grafted_phylogenies(no_seed), "seeds"))
  expect_true(all(is.na(attr(flatten_grafted_phylogenies(no_seed), "phylogenies")$seed)))
})

test_that("flatten_grafted_phylogenies() makes duplicated names unique with a warning, and the result works with the other functions",
{
  toy <- toy_setup()
  mtl <- flat_missing_taxa(toy)
  multi <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 2, verbose = FALSE, seed = 4)
  twice <- list(a = multi, b = multi)
  names(twice[[1]]) <- c("same", "same")
  names(twice[[2]]) <- c("same", "x")
  expect_warning(flat <- flatten_grafted_phylogenies(twice), "same name")
  expect_equal(anyDuplicated(names(flat)), 0)

  # The flat set can be checked and updated by the other functions
  out <- flat_graft(toy, nb_simulations = 2)
  flat <- flatten_grafted_phylogenies(out)
  checks <- check_grafting_constraints(flat, toy$phy, mtl, verbose = FALSE)
  expect_s3_class(checks, "data.frame")
  expect_equal(nrow(checks), 6 * 5)
  expect_true(attr(checks, "all_ok"))
  expect_equal(unique(checks$map), c("hyp_a", "hyp_b", "hyp_c"))
  expect_equal(checks$simulation, rep(rep(1:2, each = 5), 3))
  updated <- identify_grafted(flat, mtl$Taxon_name)
  expect_equal(names(updated), names(flat))
  expect_equal(attr(updated, "phylogenies"), attr(flat, "phylogenies"))
})
