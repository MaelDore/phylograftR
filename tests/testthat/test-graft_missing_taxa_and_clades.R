##### Test suite for graft_missing_taxa_and_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## ---- Helpers ---- ##

# Lists of missing taxa and missing clades of the toy backbone, provided together
both_lists <- function (toy)
{
  list(taxa = list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(X1 = toy$A, X2 = toy$B, X3 = toy$s)),
       clades = list(Taxon_name = c("M1", "M2", "M3"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C)))
}
as_multi <- function (maps, names = NULL) { names(maps) <- names; structure(maps, class = c("multiSimmap", "multiPhylo")) }

# Same tree, same grafted flags, same maps (the attributes added by the runs are ignored)
same_grafted <- function (x, y)
{
  isTRUE(all.equal(ape::as.phylo(x)$edge, ape::as.phylo(y)$edge)) && isTRUE(all.equal(x$edge.length, y$edge.length, tolerance = 0)) &&
    identical(x$tip.label, y$tip.label) && identical(x$grafted_tips, y$grafted_tips) && identical(x$grafted_edges, y$grafted_edges) && identical(x$maps, y$maps)
}

# Manual chain: first step, then second step on the first grafted phylogeny with its updated simmap and the same seed
manual_chain <- function (toy, lists, order = c("taxa", "clades"), seed = 1234, simmap = NULL, method = "uniform", ...)
{
  run <- function (step, phy, simmap)
  {
    if (step == "taxa") { quiet_graft(graft_missing_taxa, phy, lists$taxa, method = method, simmap = simmap, seed = seed, ...) }
    else { quiet_graft(graft_missing_clades, phy, lists$clades, simmap = simmap, seed = seed, ...) }
  }
  g1 <- run(order[1], toy$phy, simmap)
  run(order[2], g1, if (is.null(simmap)) NULL else g1$simmap)
}


## ---- Same result as the manual chain ---- ##

test_that("a single simulation is the same as grafting the taxa and then the clades with the same seed",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", seed = 11)
  expect_s3_class(both, "phylo")
  expect_true(same_grafted(both, manual_chain(toy, lists, seed = 11)))
  expect_equal(sum(both$grafted_tips), 6)
  expect_equal(attr(both, "seed"), 11L)
  expect_equal(attr(both, "grafting_order"), c("taxa", "clades"))
  expect_true(ape::is.ultrametric(both))
  expect_no_error(graft_missing_taxa_and_clades(toy$phy, lists, verbose = FALSE)) # default method BD_clades
})

test_that("grafting_order = c('clades', 'taxa') grafts the clades first and is the same as the manual chain in this order",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, grafting_order = c("clades", "taxa"), method = "uniform", seed = 3)
  expect_true(same_grafted(both, manual_chain(toy, lists, order = c("clades", "taxa"), seed = 3)))
  expect_equal(attr(both, "grafting_order"), c("clades", "taxa"))
  expect_equal(sum(both$grafted_tips), 6)
  expect_false(same_grafted(both, quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", seed = 3)))
})

test_that("several simulations are the manual chain of each simulation, with the seed of the simulation for both steps",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", nb_simulations = 3, seed = 5)
  expect_s3_class(both, "multiPhylo")
  expect_equal(names(both), c("sim_1", "sim_2", "sim_3"))
  seeds <- attr(both, "seeds")
  expect_equal(seeds, derive_simulation_seeds(5, 3) |> stats::setNames(names(both)))
  expect_equal(nrow(attr(both, "failed")), 0)
  expect_equal(names(attr(both, "failed")), c("seed", "reason", "step"))
  for (i in 1:3)
  {
    expect_true(same_grafted(both[[i]], manual_chain(toy, lists, seed = seeds[[i]])))
    expect_equal(attr(both[[i]], "seed"), as.integer(seeds[[i]]))
    expect_equal(sum(both[[i]]$grafted_tips), 6)
  }
  # The first simulation is the single simulation with the same seed
  expect_true(same_grafted(both[[1]], quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", seed = 5)))
})

test_that("total richness of both types of taxa is passed to each step",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  richness <- list(taxa = list(Target_clade = list(toy$A, toy$B, toy$s), total_richness = c(toy$nA + 15, toy$nB + 5, NA)),
                   clades = list(Sister_clade = list(toy$B, toy$C), total_richness = c(toy$nB + 12, NA)))
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, total_richness_per_clades = richness, seed = 2, return_BD_model_summary = TRUE)
  manual <- manual_chain(toy, lists, seed = 2, method = "BD_clades")
  expect_false(same_grafted(both, manual))
  g1 <- quiet_graft(graft_missing_taxa, toy$phy, lists$taxa, total_richness_per_clades = richness$taxa, seed = 2)
  g2 <- quiet_graft(graft_missing_clades, g1, lists$clades, total_richness_per_clades = richness$clades, seed = 2)
  expect_true(same_grafted(both, g2))
  # Only one of the two elements can be provided
  expect_no_error(quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, total_richness_per_clades = list(clades = richness$clades), method = "uniform"))
})

test_that("the birth-death models of both steps are reported in one table with the column 'step'",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, return_BD_model_summary = TRUE, seed = 4)
  table <- attr(both, "BD_model_summary")
  expect_equal(names(table)[1], "step")
  expect_equal(attr(table, "type"), "taxa_and_clades")
  expect_setequal(unique(table$step), c("taxa", "clades"))
  # The rows follow the order of the steps
  expect_equal(table$step, c(sort(table$step[table$step == "taxa"]), table$step[table$step == "clades"]))
  expect_equal(nrow(table), length(unique(lists$taxa$Target_clade)) + length(unique(lists$clades$Sister_clade)))
  expect_s3_class(extract_BD_model_summary(both), "data.frame")
  # Grafting the clades first
  rev <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, return_BD_model_summary = TRUE, grafting_order = c("clades", "taxa"), seed = 4)
  expect_equal(attr(rev, "BD_model_summary")$step[1], "clades")
  # method = "uniform": no model for the taxa, and no error
  uni <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, return_BD_model_summary = TRUE, method = "uniform", seed = 4)
  expect_equal(unique(attr(uni, "BD_model_summary")$step), "clades")
  # Not requested: no attribute
  expect_null(attr(quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform"), "BD_model_summary"))
})


## ---- Stochastic maps ---- ##

test_that("P maps x N simulations give P x N phylograftes, each the manual chain with the map of the simulation for both steps",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  maps <- list(make_toy_simmap(toy$phy, seed = 1), make_toy_simmap(toy$phy, seed = 2))
  multi <- as_multi(maps, c("hyp_a", "hyp_b"))
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", simmap = multi, nb_simulations = 2, seed = 7)
  expect_equal(names(both), c("hyp_a", "hyp_b"))
  expect_equal(names(both$hyp_a), c("hyp_a_sim_1", "hyp_a_sim_2"))
  expect_s3_class(attr(both, "failed_maps"), "data.frame")
  seeds <- derive_simulation_seeds(7, 2)
  for (k in 1:2) for (i in 1:2)
  {
    expected <- manual_chain(toy, lists, seed = seeds[i], simmap = maps[[k]])
    expect_true(same_grafted(both[[k]][[i]], expected))
    expect_equal(attr(both[[k]][[i]], "map"), c("hyp_a", "hyp_b")[k])
    expect_s3_class(both[[k]][[i]]$simmap, "simmap")
    expect_equal(sum(both[[k]][[i]]$grafted_tips), 6)
  }
  # The same map and seed give the same phylogeny for both maps with identical maps (maps are compared with the same random draws)
  expect_equal(attr(both$hyp_a, "seeds"), attr(both$hyp_b, "seeds") |> stats::setNames(names(attr(both$hyp_a, "seeds"))))
})

test_that("return_simmap = FALSE removes the simmap from the output, but the second step still uses the updated simmap",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  map <- make_toy_simmap(toy$phy, seed = 1)
  with_map <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", simmap = map, seed = 8)
  without <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", simmap = map, return_simmap = FALSE, seed = 8)
  expect_s3_class(with_map$simmap, "simmap")
  expect_null(without$simmap)
  without$simmap <- with_map$simmap <- NULL
  expect_true(same_grafted(with_map, without))
})

test_that("several traits are paired and given to both steps",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  traits <- list(color = make_toy_simmap(toy$phy, seed = 1), other = make_toy_simmap(toy$phy, seed = 2))
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", simmap = traits, seed = 9)
  expect_s3_class(both, "phylo")
  expect_equal(names(both$simmap), c("color", "other"))
  expect_true(same_grafted(both, manual_chain(toy, lists, seed = 9, simmap = traits)))
})


## ---- Failures ---- ##

# Make a step fail (failure specific to a simulation) for some seeds
fail_step_for_seeds <- function (step, seeds)
{
  original <- if (step == "taxa") phylograftR::graft_missing_taxa else phylograftR::graft_missing_clades
  function (..., seed) { if (!is.null(seed) && seed %in% seeds) stop_simulation_failure("no usable age for seed ", seed) else original(..., seed = seed) }
}

test_that("a failure at one step skips the simulation, with the step in the table of failures",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  seeds <- derive_simulation_seeds(21, 4)
  local_mocked_bindings(graft_missing_clades = fail_step_for_seeds("clades", seeds[2]), graft_missing_taxa = fail_step_for_seeds("taxa", seeds[3]))
  expect_warning(both <- graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", nb_simulations = 4, seed = 21, verbose = FALSE), "Only 2 of the 4")
  expect_equal(names(both), c("sim_1", "sim_4"))
  failed <- attr(both, "failed")
  expect_equal(failed$seed, seeds[2:3])
  expect_equal(failed$step, c("clades", "taxa"))
  expect_match(failed$reason[1], "^Step 2/2 \\(missing clades\\): no usable age")
  expect_match(failed$reason[2], "^Step 1/2 \\(missing taxa\\): no usable age")
  expect_equal(unname(attr(both, "seeds")), seeds[c(1, 4)])
})

test_that("force_simulations redoes both steps with new seeds until the complete phylogenies requested are obtained",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  seeds <- derive_simulation_seeds(21, 10)
  local_mocked_bindings(graft_missing_clades = fail_step_for_seeds("clades", seeds[2]))
  both <- suppressWarnings(graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", nb_simulations = 3, force_simulations = TRUE, seed = 21, verbose = FALSE))
  expect_equal(names(both), c("sim_1", "sim_2", "sim_3"))
  expect_equal(unname(attr(both, "seeds")), seeds[c(1, 3, 4)])
  expect_equal(attr(both, "failed")$step, "clades")
  for (x in both) { expect_equal(sum(x$grafted_tips), 6) }
})

test_that("a single simulation that fails raises the error, with the step",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  local_mocked_bindings(graft_missing_clades = fail_step_for_seeds("clades", 5))
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", seed = 5, verbose = FALSE), "Step 2/2 \\(missing clades\\): no usable age for seed 5")
})

test_that("a map without any simulation is dropped with its failure, and the table of failures keeps the column 'step'",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  maps <- list(make_toy_simmap(toy$phy, seed = 1), make_toy_simmap(toy$phy, seed = 2))
  multi <- as_multi(maps, c("a", "b"))
  original <- phylograftR::graft_missing_taxa
  local_mocked_bindings(graft_missing_taxa = function (..., simmap) {
    if (!is.null(simmap) && identical(simmap$maps, maps[[2]]$maps)) stop_simulation_failure("map b cannot be used")
    original(..., simmap = simmap)
  })
  expect_warning(both <- graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", simmap = multi, nb_simulations = 2, seed = 1, verbose = FALSE), "Map 'b'")
  expect_equal(names(both), "a")
  expect_equal(attr(both, "failed_maps")$map, "b")
  expect_equal(names(attr(both$a, "failed")), c("seed", "reason", "step"))
  expect_equal(nrow(attr(both$a, "failed")), 0)
})


## ---- Temporary saves ---- ##

test_that("each step of each phylogeny has its own temporary file, deleted at the end of the run",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  dir <- tempfile("wrap_temp_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "temp.rds")
  kept <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", nb_simulations = 2, seed = 6, save_temp = TRUE, temp_file_path = path, clean_temp_files = FALSE)
  expect_setequal(list.files(dir), c("temp_step_taxa_sim_1.rds", "temp_step_taxa_sim_2.rds", "temp_step_clades_sim_1.rds", "temp_step_clades_sim_2.rds"))
  cleaned <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", nb_simulations = 2, seed = 6, save_temp = TRUE, temp_file_path = path, resume = TRUE)
  expect_length(list.files(dir), 0)
  expect_true(all(mapply(same_grafted, kept, cleaned)))
})

test_that("a run interrupted during the second step resumes from the temporary files, and gives the same phylogeny",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  dir <- tempfile("wrap_resume_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "temp.rds")
  reference <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", seed = 12)
  original <- phylograftR::graft_missing_clades
  local_mocked_bindings(graft_missing_clades = function (...) stop("crash"))
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", seed = 12, save_temp = TRUE, temp_file_path = path, verbose = FALSE), "crash")
  # The first step is complete and kept after the crash
  expect_equal(list.files(dir), "temp_step_taxa.rds")
  local_mocked_bindings(graft_missing_clades = original)
  resumed <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", seed = 12, save_temp = TRUE, temp_file_path = path, resume = TRUE)
  expect_true(same_grafted(resumed, reference))
  expect_length(list.files(dir), 0)
})


## ---- Inputs ---- ##

test_that("invalid inputs stop before any grafting, with informative messages",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  expect_error(graft_missing_taxa_and_clades("phy", lists), "'phy' must be an object of class 'phylo'")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, grafting_order = "taxa"), "'grafting_order' must be")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, grafting_order = c("taxa", "taxa")), "'grafting_order' must be")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists$taxa), "graft_missing_taxa\\(\\) instead")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists$clades), "graft_missing_clades\\(\\) instead")
  expect_error(graft_missing_taxa_and_clades(toy$phy, list(taxa = lists$taxa)), "two elements named 'taxa' and 'clades'")
  dup <- lists; dup$clades$Taxon_name[1] <- "X1"; names(dup$clades$Sister_clade)[1] <- "X1"
  expect_error(graft_missing_taxa_and_clades(toy$phy, dup), "unique across.*Duplicated names: X1")
  inbb <- lists; inbb$clades$Taxon_name[1] <- toy$phy$tip.label[1]; names(inbb$clades$Sister_clade)[1] <- toy$phy$tip.label[1]
  expect_error(graft_missing_taxa_and_clades(toy$phy, inbb), "absent from 'phy\\$tip.label'")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, total_richness_per_clades = list(Target_clade = list(toy$A), total_richness = 40)), "'total_richness_per_clades' must be NULL or a list")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, return_BD_model_summary = NA), "'return_BD_model_summary' must be TRUE or FALSE")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, return_simmap = "yes"), "'return_simmap' must be TRUE or FALSE")
  expect_error(graft_missing_taxa_and_clades(toy$phy, lists, nb_simulations = 0), "'nb_simulations' must be a single positive integer")
  nonultra <- toy$phy; nonultra$edge.length[1] <- nonultra$edge.length[1] * 2
  expect_error(graft_missing_taxa_and_clades(nonultra, lists, verbose = FALSE), "ultrametric")
})

test_that("verbose reports the two steps of each simulation",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  expect_output(graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", verbose = TRUE), "Step 1/2: grafting missing taxa.*Step 2/2: grafting missing clades")
  expect_output(graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", grafting_order = c("clades", "taxa"), verbose = TRUE), "Step 1/2: grafting missing clades.*Step 2/2: grafting missing taxa")
  expect_silent(graft_missing_taxa_and_clades(toy$phy, lists, method = "uniform", verbose = FALSE))
})

test_that("the downstream functions work on the output (flags of both steps)",
{
  toy <- toy_setup(); lists <- both_lists(toy)
  both <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", nb_simulations = 2, seed = 13)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_grafted(both, index = 1))
  info <- extract_sampling_info(both, list(Target_clade = list(toy$A), total_richness = toy$nA + 30))
  expect_equal(nrow(info), 2)
  expect_true(all(info$nb_grafted_species >= 1))
})


## ---- Parallel ---- ##

skip_if_no_workers_wrap <- function ()
{
  testthat::skip_on_cran(); testthat::skip_if_not_installed("future"); testthat::skip_if_not_installed("future.apply")
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  ok <- tryCatch(isTRUE(future::value(future::future(requireNamespace("phylograftR", quietly = TRUE)))), error = function (e) FALSE)
  if (!ok) { testthat::skip("the workers cannot load the package phylograftR (install it first)") }
}

test_that("running the simulations in parallel gives the same phylogenies as sequentially",
{
  skip_if_no_workers_wrap()
  toy <- toy_setup(); lists <- both_lists(toy)
  multi <- as_multi(list(make_toy_simmap(toy$phy, seed = 1), make_toy_simmap(toy$phy, seed = 2)), c("a", "b"))
  sequential <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", simmap = multi, nb_simulations = 2, seed = 3, parallel = FALSE)
  old <- future::plan(future::multisession, workers = 2); on.exit(future::plan(old), add = TRUE)
  parallel <- quiet_graft(graft_missing_taxa_and_clades, toy$phy, lists, method = "uniform", simmap = multi, nb_simulations = 2, seed = 3, parallel = TRUE)
  for (k in 1:2) for (i in 1:2) { expect_true(same_grafted(sequential[[k]][[i]], parallel[[k]][[i]])) }
})
