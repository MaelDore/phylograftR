
##### Test suite for the multiple simulations of graft_missing_taxa() and graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Toy lists of missing taxa / missing clades
multi_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}

multi_missing_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}

## Remove what depends on the seed only through the phylogeny, to compare two phylogenies
same_phylo <- function (x, y) { isTRUE(all.equal(x, y, tolerance = 0)) }

## Check that removing the grafted taxa gives back the backbone
expect_backbone_preserved_multi <- function (grafted, phy, taxa)
{
  pruned <- ape::drop.tip(grafted, taxa)
  expect_true(ape::all.equal.phylo(pruned, phy, use.edge.length = TRUE, tolerance = 1e-6))
}


## Run the two grafting functions in the same way
graft_functions <- list(
  taxa_BD = function (toy, ...) graft_missing_taxa(toy$phy, multi_missing_taxa(toy), method = "BD_clades", verbose = FALSE, ...),
  taxa_uniform = function (toy, ...) graft_missing_taxa(toy$phy, multi_missing_taxa(toy), method = "uniform", verbose = FALSE, ...),
  clades = function (toy, ...) graft_missing_clades(toy$phy, multi_missing_clades(toy), verbose = FALSE, ...))
taxa_of <- list(taxa_BD = paste0("X", 1:5), taxa_uniform = paste0("X", 1:5), clades = paste0("M", 1:4))


## Test the output
for (name in names(graft_functions))
{
  test_that(paste0("multiple simulations (", name, "): 1 simulation gives a phylo, several give a multiPhylo named sim_1 ... sim_N with the attributes seeds and failed"),
  {
    toy <- toy_setup()
    graft <- graft_functions[[name]]
    taxa <- taxa_of[[name]]

    single <- graft(toy, seed = 5)
    expect_s3_class(single, "phylo")
    expect_false(inherits(single, "multiPhylo"))
    expect_equal(same_phylo(single, graft(toy, nb_simulations = 1, seed = 5)), TRUE)

    multi <- graft(toy, nb_simulations = 4, seed = 5)
    expect_s3_class(multi, "multiPhylo")
    expect_equal(names(multi), paste0("sim_", 1:4))
    expect_equal(names(attr(multi, "seeds")), names(multi))
    expect_equal(unname(attr(multi, "seeds")[1]), 5)
    expect_equal(nrow(attr(multi, "failed")), 0)
    for (phy_i in multi)
    {
      expect_s3_class(phy_i, "phylo")
      expect_true(ape::is.ultrametric(phy_i))
      expect_true(all(taxa %in% phy_i$tip.label))
      expect_equal(sum(phy_i$grafted_tips), length(taxa))
      expect_backbone_preserved_multi(phy_i, toy$phy, taxa)
    }

    # The simulations are different phylogenies (stochastic), and the first one is the single run with the same seed
    expect_true(same_phylo(multi$sim_1, single))
    expect_false(same_phylo(multi$sim_1, multi$sim_2))
    expect_false(same_phylo(multi$sim_2, multi$sim_3))
  })

  test_that(paste0("multiple simulations (", name, "): each simulation is reproduced alone with its seed, and the first simulations of a longer run are the same"),
  {
    toy <- toy_setup()
    graft <- graft_functions[[name]]

    multi <- graft(toy, nb_simulations = 3, seed = 11)
    seeds <- attr(multi, "seeds")
    for (k in 1:3) { expect_true(same_phylo(multi[[k]], graft(toy, seed = seeds[[k]]))) }

    longer <- graft(toy, nb_simulations = 5, seed = 11)
    expect_equal(unname(attr(longer, "seeds")[1:3]), unname(seeds))
    for (k in 1:3) { expect_true(same_phylo(longer[[k]], multi[[k]])) }

    # Same call, same result
    again <- graft(toy, nb_simulations = 3, seed = 11)
    for (k in 1:3) { expect_true(same_phylo(again[[k]], multi[[k]])) }
  })

  test_that(paste0("multiple simulations (", name, "): the RNG of the session is restored, and seed = NULL follows the session RNG without 'seeds' attribute"),
  {
    toy <- toy_setup()
    graft <- graft_functions[[name]]

    set.seed(321)
    before <- .Random.seed
    invisible(graft(toy, nb_simulations = 3, seed = 5))
    expect_identical(.Random.seed, before)

    set.seed(321)
    a <- graft(toy, nb_simulations = 2, seed = NULL)
    set.seed(321)
    b <- graft(toy, nb_simulations = 2, seed = NULL)
    expect_true(same_phylo(a[[1]], b[[1]]))
    expect_true(same_phylo(a[[2]], b[[2]]))
    expect_false(same_phylo(a[[1]], a[[2]]))
    expect_null(attr(a, "seeds"))
  })

  test_that(paste0("multiple simulations (", name, "): invalid 'nb_simulations' and 'force_simulations' stop before any grafting"),
  {
    toy <- toy_setup()
    graft <- graft_functions[[name]]
    for (value in list(0, -2, 1.5, NA, "3", c(1, 2)))
    {
      expect_error(graft(toy, nb_simulations = value), "'nb_simulations' must be a single positive integer")
    }
    expect_error(graft(toy, force_simulations = NA), "'force_simulations' must be TRUE or FALSE")
    expect_error(graft(toy, nb_simulations = 2, force_simulations = "yes"), "'force_simulations' must be TRUE or FALSE")
  })

  test_that(paste0("multiple simulations (", name, "): the warnings about the inputs are sent once, and a generic error stops the run at once"),
  {
    toy <- toy_setup()
    graft <- graft_functions[[name]]
    # A generic error (invalid input) is raised as it is, whatever the number of simulations
    expect_error(graft(toy, nb_simulations = 3, ultrametric_tol = "a"))
    phy_bad <- toy$phy
    phy_bad$edge.length[1] <- phy_bad$edge.length[1] * 1.5
    toy_bad <- toy; toy_bad$phy <- phy_bad
    expect_error(graft(toy_bad, nb_simulations = 3), "ultrametric")
  })

  test_that(paste0("multiple simulations (", name, "): with 'save_temp', each simulation saves its own file with the suffix _sim_k"),
  {
    toy <- toy_setup()
    graft <- graft_functions[[name]]
    dir <- tempfile("sims"); dir.create(dir)
    path <- file.path(dir, "temp.rds")
    out <- graft(toy, nb_simulations = 2, save_temp = TRUE, save_temp_every = 1, temp_file_path = path, seed = 5)
    expect_setequal(list.files(dir), c("temp_sim_1.rds", "temp_sim_2.rds"))
    expect_s3_class(readRDS(file.path(dir, "temp_sim_2.rds")), "phylo")

    # With a single simulation, the path is used as it is
    dir1 <- tempfile("sim1"); dir.create(dir1)
    path1 <- file.path(dir1, "single.rds")
    out1 <- graft(toy, nb_simulations = 1, save_temp = TRUE, save_temp_every = 1, temp_file_path = path1, seed = 5)
    expect_equal(list.files(dir1), "single.rds")
    unlink(c(dir, dir1), recursive = TRUE)
  })
}


## Test the multiple simulations with a simmap
test_that("multiple simulations with a simmap: each phylogeny holds its own $simmap, consistent with its tree, and the constraints hold",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- multi_missing_taxa(toy)

  for (method in c("uniform", "BD_clades"))
  {
    multi <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = method, nb_simulations = 3, verbose = FALSE, seed = 2))
    expect_equal(names(multi), paste0("sim_", 1:3))
    for (k in 1:3)
    {
      expect_s3_class(multi[[k]]$simmap, "simmap")
      expect_equal(multi[[k]]$simmap$tip.label, multi[[k]]$tip.label)
      expect_equal(unname(attr(suppressWarnings(check_grafting_constraints(multi[[k]], toy$phy, mtl, verbose = FALSE)), "all_ok")), TRUE)
      single <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = method, verbose = FALSE, seed = attr(multi, "seeds")[[k]]))
      expect_true(same_phylo(multi[[k]], single))
    }
    expect_false(same_phylo(multi[[1]]$simmap, multi[[2]]$simmap))

    # Without the simmap in the output
    no_simmap <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B"), method = method, return_simmap = FALSE, nb_simulations = 2, verbose = FALSE, seed = 2))
    expect_null(no_simmap[[1]]$simmap)
    expect_null(no_simmap[[2]]$simmap)
  }

  mcl <- multi_missing_clades(toy)
  multi_clades <- suppressWarnings(graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B"), nb_simulations = 3, verbose = FALSE, seed = 2))
  for (k in 1:3) { expect_s3_class(multi_clades[[k]]$simmap, "simmap") }
  expect_true(same_phylo(multi_clades[[1]], suppressWarnings(graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = c("A", "B"), verbose = FALSE, seed = 2))))
})

test_that("multiple simulations: the warnings about the inputs (e.g., 'min_valid_fraction') are sent once for all the simulations",
{
  toy <- toy_setup()
  simmap <- make_young_state_simmap(toy$phy, 4)
  taxa <- paste0("Y", 1:3)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 3), taxa))
  n_warnings <- 0
  withCallingHandlers(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", nb_simulations = 3, verbose = FALSE, seed = 1),
                      warning = function (w) { n_warnings <<- n_warnings + 1; invokeRestart("muffleWarning") })
  expect_equal(n_warnings, 1)
})


## Test the failures specific to a simulation
test_that("graft_missing_taxa() with several simulations skips the simulations without usable divergence age, with a warning, and reports the seeds and reasons",
{
  toy <- toy_setup()
  # The state "B" only exists in a narrow window of the youngest ages: with one redraw, some simulations fail and some succeed
  simmap <- make_young_state_simmap(toy$phy, 1)
  taxa <- paste0("Y", 1:2)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 2), taxa))
  graft <- function (...) graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 1, verbose = FALSE, ...)

  warnings_sent <- character(0)
  out <- withCallingHandlers(graft(nb_simulations = 10, seed = 1),
                             warning = function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") })
  n_obtained <- length(out)
  expect_gt(n_obtained, 0)
  expect_lt(n_obtained, 10)
  expect_match(warnings_sent, "Only [0-9]+ of the 10 simulations requested were obtained", all = FALSE)
  expect_match(warnings_sent, "No usable divergence age", all = FALSE)

  # Names keep the index of the simulation: some are missing
  indices <- as.integer(sub("sim_", "", names(out)))
  expect_true(all(indices %in% 1:10))
  expect_true(any(diff(indices) > 1) || min(indices) > 1 || max(indices) < 10)
  failed <- attr(out, "failed")
  expect_equal(nrow(failed), 10 - n_obtained)
  expect_match(failed$reason, "No usable divergence age")
  seeds_all <- derive_simulation_seeds(1, 10)
  expect_equal(unname(c(attr(out, "seeds"), failed$seed)[order(match(c(attr(out, "seeds"), failed$seed), seeds_all))]), seeds_all)

  # A simulation obtained is the same as its single run, and a failed seed fails alone with the same error
  k <- 1
  single <- suppressWarnings(graft(seed = attr(out, "seeds")[[k]]))
  expect_true(same_phylo(out[[k]], single))
  expect_error(suppressWarnings(graft(seed = failed$seed[1])), "No usable divergence age")
})

test_that("graft_missing_taxa() with 'force_simulations' tries new seeds until the number of simulations requested is obtained",
{
  toy <- toy_setup()
  simmap <- make_young_state_simmap(toy$phy, 1)
  taxa <- paste0("Y", 1:2)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 2), taxa))
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 1, verbose = FALSE, ...))

  without_force <- graft(nb_simulations = 10, seed = 1)
  expect_lt(length(without_force), 10)

  forced <- graft(nb_simulations = 6, force_simulations = TRUE, seed = 1)
  expect_length(forced, 6)
  expect_equal(names(forced), paste0("sim_", 1:6))
  expect_gt(nrow(attr(forced, "failed")), 0)
  # The phylogenies of a run without 'force_simulations' are the first ones of a run with it
  expect_true(same_phylo(forced[[1]], without_force[[1]]))
  # Each phylogeny is reproduced alone with its seed
  expect_true(same_phylo(forced[[6]], graft(seed = attr(forced, "seeds")[[6]])))
  expect_false(any(attr(forced, "seeds") %in% attr(forced, "failed")$seed))
})

test_that("graft_missing_taxa() stops when no simulation is obtained, and stops on a generic failure without trying the other simulations",
{
  toy <- toy_setup()
  simmap <- make_young_state_simmap(toy$phy, 1e-4)
  taxa <- paste0("Y", 1:2)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 2), taxa))
  expect_error(suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 2,
                                                   nb_simulations = 3, verbose = FALSE, seed = 1)),
               "None of the 3 simulation\\(s\\) could be obtained.*Reasons:.*No usable divergence age")
  expect_error(suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 2,
                                                   nb_simulations = 2, force_simulations = TRUE, verbose = FALSE, seed = 1)),
               "None of the 20 simulation\\(s\\) could be obtained")
  # With a single simulation, the error is the usual one
  expect_error(suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = "B", method = "BD_clades", max_age_redraws = 2,
                                                   nb_simulations = 1, verbose = FALSE, seed = 1)),
               "No usable divergence age.*after 2 redraws")

  # A generic failure (taxon already in the backbone) stops the whole run before any simulation, with the usual message
  mtl_present <- mtl
  mtl_present$Taxon_name[1] <- toy$phy$tip.label[1]
  names(mtl_present$Target_clade)[1] <- toy$phy$tip.label[1]
  expect_error(graft_missing_taxa(toy$phy, mtl_present, nb_simulations = 3, verbose = FALSE, seed = 1), "backbone")
})

test_that("graft_missing_clades() with several simulations skips the simulations without usable divergence time, and 'force_simulations' reaches the number requested",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  # Only a small piece of the stem branch of B (at its tipward end) is in state "C"
  simmap <- set_branch_states(simmap, toy$B, c("A", "C"), c(1, 0.1))
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B))
  graft <- function (...) graft_missing_clades(toy$phy, mcl, simmap = simmap, allowed_states = "C", max_age_redraws = 1, verbose = FALSE, ...)

  warnings_sent <- character(0)
  out <- withCallingHandlers(graft(nb_simulations = 10, seed = 1),
                             warning = function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_gt(length(out), 0)
  expect_lt(length(out), 10)
  expect_match(warnings_sent, "Only [0-9]+ of the 10 simulations requested were obtained.*No usable divergence time", all = FALSE)
  expect_equal(nrow(attr(out, "failed")), 10 - length(out))
  expect_match(attr(out, "failed")$reason, "No usable divergence time")
  expect_true(same_phylo(out[[1]], suppressWarnings(graft(seed = attr(out, "seeds")[[1]]))))
  expect_error(suppressWarnings(graft(seed = attr(out, "failed")$seed[1])), "No usable divergence time")

  forced <- suppressWarnings(graft(nb_simulations = 5, force_simulations = TRUE, seed = 1))
  expect_equal(names(forced), paste0("sim_", 1:5))
  expect_gt(nrow(attr(forced, "failed")), 0)
  expect_true(same_phylo(forced[[5]], suppressWarnings(graft(seed = attr(forced, "seeds")[[5]]))))

  # None obtained
  simmap_none <- set_branch_states(make_toy_simmap(toy$phy), toy$B, c("A", "C"), c(1, 1e-7))
  expect_error(suppressWarnings(graft_missing_clades(toy$phy, mcl, simmap = simmap_none, allowed_states = "C", max_age_redraws = 1, nb_simulations = 3, verbose = FALSE, seed = 1)),
               "None of the 3 simulation\\(s\\) could be obtained.*No usable divergence time")
})
