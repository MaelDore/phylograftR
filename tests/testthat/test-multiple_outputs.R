##### Test suite for check_grafting_constraints() and identify_grafted() with several grafted phylogenies #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


outputs_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}
as_multi_simmap_outputs <- function (maps, names = NULL)
{
  names(maps) <- names
  structure(maps, class = c("multiSimmap", "multiPhylo"))
}


## ---- list_grafted_phylogenies() ---- ##

test_that("list_grafted_phylogenies() lists the phylogenies, their maps, simulations and seeds",
{
  toy <- toy_setup()
  mtl <- outputs_missing_taxa(toy)
  single <- graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, seed = 4)
  multi <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 3, verbose = FALSE, seed = 4)

  l1 <- list_grafted_phylogenies(single)
  expect_true(is.na(l1$map)); expect_equal(l1$simulation, 1); expect_equal(l1$seed, 4L); expect_false(l1$is_list)

  l3 <- list_grafted_phylogenies(multi)
  expect_equal(l3$map, rep(NA_character_, 3)); expect_equal(l3$simulation, 1:3); expect_equal(l3$seed, unname(attr(multi, "seeds"))); expect_false(l3$is_list)

  # A single simmap without name: map_1, and no map (NA) without simmap
  one_map <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = make_toy_simmap(toy$phy), allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4))
  expect_equal(list_grafted_phylogenies(one_map)$map, "map_1")
  no_simmap_output <- one_map; no_simmap_output$simmap <- NULL
  expect_true(is.na(list_grafted_phylogenies(no_simmap_output)$map))

  # A multiPhylo with failed simulations: the index of the simulation is the one of the names
  gap <- multi[c(1, 3)]
  expect_equal(list_grafted_phylogenies(gap)$simulation, c(1L, 3L))

  # List of maps
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  out <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap_outputs(maps, c("hyp_a", "hyp_b")), allowed_states = c("A", "B", "C"),
                                             method = "uniform", nb_simulations = 2, verbose = FALSE, seed = 4))
  lm <- list_grafted_phylogenies(out)
  expect_true(lm$is_list)
  expect_equal(lm$map, c("hyp_a", "hyp_a", "hyp_b", "hyp_b"))
  expect_equal(lm$simulation, c(1L, 2L, 1L, 2L))
  out1 <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap_outputs(maps), allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4))
  lm1 <- list_grafted_phylogenies(out1)
  expect_equal(lm1$map, c("map_1", "map_2")); expect_equal(lm1$simulation, c(1L, 1L)); expect_true(lm1$is_list)

  expect_error(list_grafted_phylogenies("a"), "must be a phylogeny")
  expect_error(list_grafted_phylogenies(list(1, 2)), "must be a phylogeny")
  expect_error(list_grafted_phylogenies(list()), "must be a phylogeny")
})


## ---- identify_grafted() ---- ##

test_that("identify_grafted() updates every phylogeny of a multiPhylo or of a list of maps, and keeps the names and the attributes",
{
  toy <- toy_setup()
  mtl <- outputs_missing_taxa(toy)
  multi <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 3, verbose = FALSE, seed = 4)
  # Manipulate the phylogenies: the elements are not valid anymore
  pruned <- multi
  for (j in seq_along(pruned)) { pruned[[j]] <- ape::drop.tip(pruned[[j]], c("X2", toy$phy$tip.label[1])); attr(pruned[[j]], "seed") <- attr(multi[[j]], "seed") }
  expect_false(length(pruned[[1]]$grafted_tips) == ape::Ntip(pruned[[1]]))
  taxa <- setdiff(mtl$Taxon_name, "X2")
  updated <- identify_grafted(pruned, taxa)
  expect_s3_class(updated, "multiPhylo")
  expect_equal(names(updated), names(multi))
  expect_equal(attr(updated, "seeds"), attr(multi, "seeds"))
  expect_equal(attr(updated, "failed"), attr(multi, "failed"))
  for (j in seq_along(updated))
  {
    expect_equal(length(updated[[j]]$grafted_tips), ape::Ntip(updated[[j]]))
    expect_equal(length(updated[[j]]$grafted_edges), nrow(updated[[j]]$edge))
    expect_equal(sum(updated[[j]]$grafted_tips), length(taxa))
    expect_identical(updated[[j]]$grafted_tips, identify_grafted(pruned[[j]], taxa)$grafted_tips)
    expect_equal(attr(updated[[j]], "seed"), attr(multi[[j]], "seed"))
  }

  # List of maps
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  out <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap_outputs(maps, c("hyp_a", "hyp_b")), allowed_states = c("A", "B", "C"),
                                             method = "uniform", nb_simulations = 2, verbose = FALSE, seed = 4))
  out_pruned <- lapply(out, function (set) { for (j in seq_along(set)) { set[[j]] <- ape::drop.tip(set[[j]], "X2") }; set })
  attributes(out_pruned) <- attributes(out)
  updated_out <- identify_grafted(out_pruned, taxa)
  expect_equal(names(updated_out), c("hyp_a", "hyp_b"))
  expect_equal(attr(updated_out, "failed_maps"), attr(out, "failed_maps"))
  expect_s3_class(updated_out$hyp_b, "multiPhylo")
  expect_equal(sum(updated_out$hyp_b[[2]]$grafted_tips), length(taxa))
  # One simulation per map: list of phylogenies
  out1 <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap_outputs(maps), allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4))
  updated_out1 <- identify_grafted(out1, mtl$Taxon_name)
  expect_equal(names(updated_out1), c("map_1", "map_2"))
  expect_equal(attr(updated_out1$map_2, "map"), "map_2")

  # Warnings about absent taxa are sent once
  n_warnings <- 0
  withCallingHandlers(identify_grafted(multi, c(mtl$Taxon_name, "unknown")), warning = function (w) { n_warnings <<- n_warnings + 1; invokeRestart("muffleWarning") })
  expect_equal(n_warnings, 1)
  expect_warning(identify_grafted(multi, c(mtl$Taxon_name, "unknown")), "1 taxa of 'grafted_taxa' are not tips")

  # Invalid inputs
  expect_error(identify_grafted("a", "X1"), "'phy' must be an object of class 'phylo'")
  expect_error(identify_grafted(list(1, 2), "X1"), "'phy' must be an object of class 'phylo'")
})


## ---- check_grafting_constraints() ---- ##

test_that("check_grafting_constraints() with a single phylogeny keeps its results and adds the columns map, simulation and seed",
{
  toy <- toy_setup()
  mtl <- outputs_missing_taxa(toy)
  single <- graft_missing_taxa(toy$phy, mtl, method = "uniform", verbose = FALSE, seed = 4)
  res <- check_grafting_constraints(single, toy$phy, mtl, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(names(res)[1:4], c("map", "simulation", "seed", "Taxon_name"))
  expect_equal(nrow(res), 5)
  expect_true(all(is.na(res$map))); expect_equal(unique(res$simulation), 1); expect_equal(unique(res$seed), 4L)
  # With a single simmap without name: map_1
  with_simmap <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = make_toy_simmap(toy$phy), allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4))
  expect_equal(unique(check_grafting_constraints(with_simmap, toy$phy, mtl, verbose = FALSE)$map), "map_1")
  expect_true(attr(res, "all_ok"))
  expect_equal(attr(res, "phylogenies")$n_failed, 0)
  # Same verbose summary as before
  expect_output(check_grafting_constraints(single, toy$phy, mtl, verbose = TRUE), "5 taxa checked.*All constraints respected")
})

test_that("check_grafting_constraints() with a multiPhylo merges the simulations in one data.frame and gives an aggregated summary",
{
  toy <- toy_setup()
  mtl <- outputs_missing_taxa(toy)
  multi <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 3, verbose = FALSE, seed = 4)
  res <- check_grafting_constraints(multi, toy$phy, mtl, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 15)
  expect_equal(res$simulation, rep(1:3, each = 5))
  expect_equal(res$seed, rep(unname(attr(multi, "seeds")), each = 5))
  expect_true(all(is.na(res$map)))
  expect_equal(res$Taxon_name, rep(mtl$Taxon_name, 3))
  expect_true(attr(res, "all_ok"))
  expect_true(attr(res, "backbone_preserved"))
  expect_equal(nrow(attr(res, "phylogenies")), 3)
  # Same rows as the check of each phylogeny alone
  for (j in 1:3)
  {
    alone <- check_grafting_constraints(multi[[j]], toy$phy, mtl, verbose = FALSE)
    expect_equal(res[res$simulation == j, "ok"], alone$ok)
    expect_equal(res[res$simulation == j, "target_position"], alone$target_position)
  }
  out <- capture.output(check_grafting_constraints(multi, toy$phy, mtl, verbose = TRUE))
  expect_match(paste(out, collapse = "\n"), "3 phylogenies checked \\(no map\\), 5 taxa each.*respecting all the constraints: 3/3.*All constraints respected")
  # One summary only
  expect_equal(sum(grepl("phylogenies checked", out)), 1)
})

test_that("check_grafting_constraints() with a list of maps returns one data.frame per map, with its own simmap checked",
{
  toy <- toy_setup()
  mtl <- outputs_missing_taxa(toy)
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  out <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap_outputs(maps, c("hyp_a", "hyp_b")), allowed_states = c("A", "B", "C"),
                                             method = "uniform", nb_simulations = 2, verbose = FALSE, seed = 4))
  res <- check_grafting_constraints(out, toy$phy, mtl, allowed_states = c("A", "B"), verbose = FALSE)
  expect_type(res, "list")
  expect_equal(names(res), c("hyp_a", "hyp_b"))
  for (m in names(res))
  {
    expect_equal(unique(res[[m]]$map), m)
    expect_equal(res[[m]]$simulation, rep(1:2, each = 5))
    expect_false(anyNA(res[[m]]$states_ok))
    expect_equal(nrow(attr(res[[m]], "phylogenies")), 2)
  }
  expect_equal(attr(res, "all_ok"), all(vapply(res, function (r) isTRUE(attr(r, "all_ok")), logical(1))))

  # States checked in each phylogeny with its own simmap: a state outside the allowed ones is detected
  res_c <- check_grafting_constraints(out, toy$phy, mtl, allowed_states = "C", verbose = FALSE)
  expect_true(any(unlist(lapply(res_c, function (r) r$states_ok %in% FALSE))))

  # One simulation per map: still a list of data.frames (with the columns map, simulation = 1)
  out1 <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap_outputs(maps), allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, seed = 4))
  res1 <- check_grafting_constraints(out1, toy$phy, mtl, verbose = FALSE)
  expect_equal(names(res1), c("map_1", "map_2"))
  expect_equal(unique(res1$map_2$simulation), 1)
  expect_equal(unique(res1$map_2$map), "map_2")

  # Aggregated summary over the maps
  text <- paste(capture.output(check_grafting_constraints(out, toy$phy, mtl, verbose = TRUE)), collapse = "\n")
  expect_match(text, "4 phylogenies checked \\(2 map\\(s\\)\\), 5 taxa each")
  expect_match(text, "Allowed states violated")
})

test_that("check_grafting_constraints() with several phylogenies reports the violations, stops on failure, and names the phylogeny in the errors",
{
  toy <- toy_setup()
  mtl <- outputs_missing_taxa(toy)
  multi <- graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 3, verbose = FALSE, seed = 4)
  # Violate the constraints in the second phylogeny: the checked list requires another target clade for X3
  multi_bad <- multi
  mtl_other <- mtl
  mtl_other$Target_clade$X3 <- toy$C
  res <- check_grafting_constraints(multi, toy$phy, mtl_other, verbose = FALSE)
  expect_false(attr(res, "all_ok"))
  failing <- attr(res, "phylogenies")
  expect_true(any(!failing$all_ok))
  text <- paste(capture.output(check_grafting_constraints(multi, toy$phy, mtl_other, verbose = TRUE)), collapse = "\n")
  expect_match(text, "Phylogenies with violations:.*no map, simulation [0-9].*taxa failed")
  expect_match(text, "Some constraints NOT respected")
  expect_error(check_grafting_constraints(multi, toy$phy, mtl_other, stop_on_failure = TRUE, verbose = FALSE),
               "Some grafting constraints are not respected in [0-9]+ of 3 phylogenies:\n  - no map, simulation [0-9]")

  # Arguments on the states need the simmap of each phylogeny: the error names the phylogeny
  expect_error(check_grafting_constraints(multi, toy$phy, mtl, allowed_states = "A", verbose = FALSE), "Phylogeny, simulation 1: 'allowed_states' can only be used")
  expect_error(check_grafting_constraints("a", toy$phy, mtl), "must be a phylogeny")

  # The warnings about the inputs are sent once
  mtl_no_excluded <- mtl; mtl_no_excluded$Excluded_subclades <- NULL
  n_warnings <- 0
  withCallingHandlers(check_grafting_constraints(multi, toy$phy, mtl_no_excluded, verbose = FALSE),
                      warning = function (w) { n_warnings <<- n_warnings + 1; invokeRestart("muffleWarning") })
  expect_equal(n_warnings, 1)
})
