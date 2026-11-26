
##### Test suite for the use of several stochastic maps ('multiSimmap') to constrain grafting in graft_missing_taxa() and graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Toy lists of missing taxa / missing clades, and multiSimmap
maps_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}
maps_missing_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}
as_multi_simmap <- function (maps, names = NULL)
{
  names(maps) <- names
  structure(maps, class = c("multiSimmap", "multiPhylo"))
}
## Remove the attribute 'map' of a phylogeny, to compare it with the one of a run with a single map
drop_attr_map <- function (phy) { attr(phy, "map") <- NULL; phy }
same_phylo <- function (x, y) { isTRUE(all.equal(x, y, tolerance = 0)) }


## ---- Helpers ---- ##

test_that("check_simmap_input() lists the maps, names them, and rejects what is not a map or a multiSimmap",
{
  expect_equal(check_simmap_input(NULL), list(maps = list(NULL), names = NULL, multiple = FALSE, traits = NULL))
  toy <- toy_setup()
  one <- make_toy_simmap(toy$phy)
  single <- check_simmap_input(one)
  expect_false(single$multiple)
  expect_null(single$names)
  expect_length(single$maps, 1)

  # No names: map_1, ..., map_P
  multi <- check_simmap_input(as_multi_simmap(list(one, one, one)))
  expect_true(multi$multiple)
  expect_equal(multi$names, c("map_1", "map_2", "map_3"))
  expect_length(multi$maps, 3)
  expect_equal(check_simmap_input(as_multi_simmap(list(one, one), c("", "")))$names, c("map_1", "map_2"))
  # One map in a multiSimmap: the output is still a list
  expect_true(check_simmap_input(as_multi_simmap(list(one)))$multiple)
  # Names are kept
  expect_equal(check_simmap_input(as_multi_simmap(list(one, one), c("hypothesis_a", "hypothesis_b")))$names, c("hypothesis_a", "hypothesis_b"))

  # Errors
  expect_error(check_simmap_input(as_multi_simmap(list(one, one), c("a", "a"))), "different names.*Duplicated names: a")
  expect_error(check_simmap_input(as_multi_simmap(list(one, one), c("a", ""))), "named and others are not.*2")
  expect_error(check_simmap_input(as_multi_simmap(list())), "without any map")
  expect_error(check_simmap_input(list(one)), "list of 1 element.*several traits")
  expect_error(check_simmap_input(toy$phy), "must be a stochastic character map of class 'simmap'")
  expect_error(check_simmap_input("map"), "must be a stochastic character map of class 'simmap'")
})

test_that("with_map_context() adds the name of the map to the errors and the warnings, only if there is a name",
{
  expect_error(with_map_context(stop("problem"), NULL), "^problem$")
  expect_error(with_map_context(stop("problem"), "map_2"), "^Map 'map_2': problem$")
  expect_warning(with_map_context(warning("careful"), "map_2"), "^Map 'map_2': careful$")
  expect_warning(with_map_context(warning("careful"), NULL), "^careful$")
  expect_equal(with_map_context(3, "map_1"), 3)
})

test_that("temp_file_path_simulation() adds the map and the simulation suffixes",
{
  expect_equal(temp_file_path_simulation("a/b.rds", 2, multiple = FALSE, map_name = "m1"), "a/b_map_m1.rds")
  expect_equal(temp_file_path_simulation("a/b.rds", 2, multiple = TRUE, map_name = "m1"), "a/b_map_m1_sim_2.rds")
  expect_equal(temp_file_path_simulation("a/b.rds", 2, multiple = TRUE, map_name = "hyp A/1"), "a/b_map_hyp_A_1_sim_2.rds")
  expect_equal(temp_file_path_simulation(NULL, 1, multiple = FALSE, map_name = "m1"), "./temp_phy_from_grafting_map_m1.rds")
  expect_equal(temp_file_path_simulation("a/b.rds", 2, multiple = TRUE), "a/b_sim_2.rds")
  # The prefix "map_" is not repeated in the name of the maps that start with it
  expect_equal(temp_file_path_simulation("a/b.rds", 2, multiple = TRUE, map_name = "map_1"), "a/b_map_1_sim_2.rds")
  expect_equal(temp_file_path_simulation("a/b.rds", 2, multiple = FALSE, map_name = "map_A"), "a/b_map_A.rds")
})

## Simulation function for run_maps(): fails for the maps named in 'failing_maps'
fake_map_simulation <- function (failing_maps = character(0))
{
  function (tracker, missing_taxa_list, map_name, sim_seed, sim_index, verbose_sim)
  {
    if (!is.null(map_name) && map_name %in% failing_maps) { stop_simulation_failure("no usable age for the map ", map_name) }
    phy <- ape::rtree(4)
    phy$info <- c(tracker = tracker, index = sim_index)
    phy
  }
}
fake_prepared <- function (n) lapply(seq_len(n), function (k) list(tracker = k, missing_taxa_list = list()))

test_that("run_maps() without multiSimmap is run_simulations() on the single map",
{
  out <- run_maps(fake_map_simulation(), fake_prepared(1), map_names = NULL, multiple = FALSE, nb_simulations = 2, force_simulations = FALSE, seed = 4)
  expect_s3_class(out, "multiPhylo")
  expect_equal(names(out), c("sim_1", "sim_2"))
  out1 <- run_maps(fake_map_simulation(), fake_prepared(1), map_names = NULL, multiple = FALSE, nb_simulations = 1, force_simulations = FALSE, seed = 4)
  expect_s3_class(out1, "phylo")
  expect_null(attr(out1, "map"))
  expect_equal(attr(out1, "seed"), 4L)
})

test_that("run_maps() with a multiSimmap returns a list named after the maps: phylogenies for 1 simulation, multiPhylo named <map>_sim_k for several",
{
  out1 <- run_maps(fake_map_simulation(), fake_prepared(3), map_names = c("a", "b", "c"), multiple = TRUE, nb_simulations = 1, force_simulations = FALSE, seed = 4)
  expect_type(out1, "list")
  expect_equal(names(out1), c("a", "b", "c"))
  for (k in 1:3)
  {
    expect_s3_class(out1[[k]], "phylo")
    expect_equal(attr(out1[[k]], "map"), c("a", "b", "c")[k])
    expect_equal(attr(out1[[k]], "seed"), 4L) # same seeds in each map
    expect_equal(unname(out1[[k]]$info["tracker"]), k) # each map has its own tracker
  }
  expect_equal(nrow(attr(out1, "failed_maps")), 0)

  out3 <- run_maps(fake_map_simulation(), fake_prepared(2), map_names = c("a", "b"), multiple = TRUE, nb_simulations = 3, force_simulations = FALSE, seed = 4)
  expect_equal(names(out3), c("a", "b"))
  for (name in c("a", "b"))
  {
    expect_s3_class(out3[[name]], "multiPhylo")
    expect_equal(names(out3[[name]]), paste0(name, "_sim_", 1:3))
    expect_equal(names(attr(out3[[name]], "seeds")), names(out3[[name]]))
    expect_equal(unname(attr(out3[[name]], "seeds")), derive_simulation_seeds(4, 3))
    expect_equal(attr(out3[[name]][[2]], "map"), name)
    expect_equal(nrow(attr(out3[[name]], "failed")), 0)
  }
  expect_equal(attr(out3$a, "seeds"), setNames(attr(out3$b, "seeds"), names(out3$a)))
})

test_that("run_maps() drops the maps without any simulation with a warning, and stops if no map gives a simulation",
{
  expect_warning(out <- run_maps(fake_map_simulation("b"), fake_prepared(3), map_names = c("a", "b", "c"), multiple = TRUE, nb_simulations = 2, force_simulations = FALSE, seed = 4),
                 "Map 'b': none of the 2 simulation\\(s\\) could be obtained, the map is dropped from the output.*no usable age for the map b")
  expect_equal(names(out), c("a", "c"))
  expect_equal(attr(out, "failed_maps")$map, "b")
  expect_match(attr(out, "failed_maps")$reason, "no usable age for the map b")

  # With one simulation per map
  expect_warning(out1 <- run_maps(fake_map_simulation("a"), fake_prepared(2), map_names = c("a", "b"), multiple = TRUE, nb_simulations = 1, force_simulations = FALSE, seed = 4),
                 "Map 'a': none of the 1 simulation\\(s\\) could be obtained")
  expect_equal(names(out1), "b")

  expect_error(suppressWarnings(run_maps(fake_map_simulation(c("a", "b")), fake_prepared(2), map_names = c("a", "b"), multiple = TRUE, nb_simulations = 2, force_simulations = FALSE, seed = 4)),
               "None of the 2 maps gave a simulation.*a: no usable age for the map a.*b: no usable age for the map b")

  # A partial failure in a map is reported for this map, with its name
  partial <- function (tracker, missing_taxa_list, map_name, sim_seed, sim_index, verbose_sim)
  {
    if (map_name == "b" && sim_index == 2) { stop_simulation_failure("failed once") }
    ape::rtree(4)
  }
  expect_warning(out <- run_maps(partial, fake_prepared(2), map_names = c("a", "b"), multiple = TRUE, nb_simulations = 3, force_simulations = FALSE, seed = 4),
                 "Map 'b': Only 2 of the 3 simulations requested were obtained")
  expect_equal(names(out$b), c("b_sim_1", "b_sim_3"))
  expect_equal(names(out$a), c("a_sim_1", "a_sim_2", "a_sim_3"))
})

test_that("run_maps() warns when maps x simulations exceeds 1000",
{
  expect_no_warning(run_maps(fake_map_simulation(), fake_prepared(2), map_names = c("a", "b"), multiple = TRUE, nb_simulations = 5, force_simulations = FALSE, seed = 4))
  expect_warning(run_maps(fake_map_simulation(), fake_prepared(11), map_names = paste0("m", 1:11), multiple = TRUE, nb_simulations = 100, force_simulations = FALSE, seed = 4),
                 "11 maps x 100 simulations = 1100 grafted phylogenies.*more than 1000")
})


## ---- Grafting functions ---- ##

graft_with_maps <- list(
  taxa_BD = function (toy, simmap, ...) graft_missing_taxa(toy$phy, maps_missing_taxa(toy), simmap = simmap, allowed_states = c("A", "B", "C"), method = "BD_clades", verbose = FALSE, ...),
  taxa_uniform = function (toy, simmap, ...) graft_missing_taxa(toy$phy, maps_missing_taxa(toy), simmap = simmap, allowed_states = c("A", "B", "C"), method = "uniform", verbose = FALSE, ...),
  clades = function (toy, simmap, ...) graft_missing_clades(toy$phy, maps_missing_clades(toy), simmap = simmap, allowed_states = c("A", "B", "C"), verbose = FALSE, ...))

for (name in names(graft_with_maps))
{
  test_that(paste0("multiSimmap (", name, "): 1 simulation per map gives a list of phylogenies named after the maps, each the same as the single run with its map"),
  {
    toy <- toy_setup()
    graft <- function (...) suppressWarnings(graft_with_maps[[name]](toy, ...))
    maps <- lapply(1:3, function (k) make_toy_simmap(toy$phy, seed = k))

    # Without names
    out <- graft(simmap = as_multi_simmap(maps), seed = 5)
    expect_type(out, "list")
    expect_false(inherits(out, "phylo"))
    expect_equal(names(out), c("map_1", "map_2", "map_3"))
    expect_equal(nrow(attr(out, "failed_maps")), 0)
    for (k in 1:3)
    {
      expect_s3_class(out[[k]], "phylo")
      expect_s3_class(out[[k]]$simmap, "simmap")
      expect_equal(attr(out[[k]], "map"), paste0("map_", k))
      expect_equal(attr(out[[k]], "seed"), 5L)
      # Same as a run with this map alone, and the same seed
      single <- graft(simmap = maps[[k]], seed = 5)
      expect_true(same_phylo(drop_attr_map(out[[k]]), single))
      if (grepl("^taxa", name)) { expect_true(isTRUE(attr(suppressWarnings(check_grafting_constraints(out[[k]], toy$phy, maps_missing_taxa(toy), verbose = FALSE)), "all_ok"))) }
    }
    # The maps are different hypotheses: the grafted phylogenies are not the same
    expect_false(identical(out$map_1$simmap$maps, out$map_2$simmap$maps))

    # With names
    named <- graft(simmap = as_multi_simmap(maps, c("hyp_a", "hyp_b", "hyp_c")), seed = 5)
    expect_equal(names(named), c("hyp_a", "hyp_b", "hyp_c"))
    expect_equal(attr(named$hyp_b, "map"), "hyp_b")
    expect_true(same_phylo(named$hyp_c, out$map_3 |> (\(x) { attr(x, "map") <- "hyp_c"; x })()))

    # A multiSimmap with one map still gives a list
    one <- graft(simmap = as_multi_simmap(maps[1]), seed = 5)
    expect_equal(names(one), "map_1")

    # Without the simmaps in the output
    no_simmap <- graft(simmap = as_multi_simmap(maps), return_simmap = FALSE, seed = 5)
    expect_null(no_simmap$map_1$simmap)
  })

  test_that(paste0("multiSimmap (", name, "): several simulations per map give a list of multiPhylo named <map>_sim_k, each simulation reproducible with its map and seed"),
  {
    toy <- toy_setup()
    graft <- function (...) suppressWarnings(graft_with_maps[[name]](toy, ...))
    maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
    out <- graft(simmap = as_multi_simmap(maps, c("hyp_a", "hyp_b")), nb_simulations = 3, seed = 5)
    expect_equal(names(out), c("hyp_a", "hyp_b"))
    for (k in 1:2)
    {
      map_name <- c("hyp_a", "hyp_b")[k]
      expect_s3_class(out[[k]], "multiPhylo")
      expect_equal(names(out[[k]]), paste0(map_name, "_sim_", 1:3))
      seeds <- attr(out[[k]], "seeds")
      expect_equal(names(seeds), names(out[[k]]))
      for (j in 1:3)
      {
        expect_s3_class(out[[k]][[j]]$simmap, "simmap")
        expect_equal(attr(out[[k]][[j]], "map"), map_name)
        expect_equal(attr(out[[k]][[j]], "seed"), seeds[[j]])
        single <- graft(simmap = maps[[k]], seed = seeds[[j]])
        expect_true(same_phylo(drop_attr_map(out[[k]][[j]]), single))
      }
    }
    # The same seeds in each map
    expect_equal(unname(attr(out$hyp_a, "seeds")), unname(attr(out$hyp_b, "seeds")))
  })

  test_that(paste0("multiSimmap (", name, "): a state absent from one map, or a plain list of maps, stops before any simulation, naming the map"),
  {
    toy <- toy_setup()
    maps <- list(make_toy_simmap(toy$phy, seed = 1), make_young_state_simmap(toy$phy, 4), make_toy_simmap(toy$phy, seed = 2)) # The second map only has the states "A" and "B"
    run <- function (map_list, ..., allowed_states = c("A", "B", "C")) {
      if (grepl("^taxa", name)) graft_missing_taxa(toy$phy, maps_missing_taxa(toy), simmap = map_list, allowed_states = allowed_states, method = if (name == "taxa_BD") "BD_clades" else "uniform", verbose = FALSE, ...)
      else graft_missing_clades(toy$phy, maps_missing_clades(toy), simmap = map_list, allowed_states = allowed_states, verbose = FALSE, ...)
    }
    calls <- 0
    original <- run_simulations
    local_mocked_bindings(run_simulations = function (...) { calls <<- calls + 1; original(...) })
    expect_error(run(as_multi_simmap(maps, c("first", "young", "third")), nb_simulations = 5), "Map 'young':.*'C'")
    expect_error(run(as_multi_simmap(maps)), "Map 'map_2':.*'C'")
    expect_equal(calls, 0)
    expect_error(run(list(maps[[1]], maps[[3]])), "must be given for each trait") # A plain list is the maps of several traits
    expect_error(run(as_multi_simmap(maps[c(1, 3)], c("a", "a"))), "different names")
  })

  test_that(paste0("multiSimmap (", name, "): with 'save_temp', each map and simulation saves its own file"),
  {
    toy <- toy_setup()
    maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
    dir <- tempfile("maps"); dir.create(dir)
    path <- file.path(dir, "temp.rds")
    out <- suppressWarnings(graft_with_maps[[name]](toy, simmap = as_multi_simmap(maps, c("a", "b")), nb_simulations = 2, save_temp = TRUE, save_temp_every = 1, clean_temp_files = FALSE, temp_file_path = path, seed = 5))
    expect_setequal(list.files(dir), c("temp_map_a_sim_1.rds", "temp_map_a_sim_2.rds", "temp_map_b_sim_1.rds", "temp_map_b_sim_2.rds"))
    unlink(dir, recursive = TRUE)
  })
}

## ---- Failures per map ---- ##

test_that("graft_missing_taxa() with a multiSimmap drops the map without any simulation, with a warning, and keeps the others",
{
  toy <- toy_setup()
  taxa <- paste0("Y", 1:2)
  mtl <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(toy$A), 2), taxa))
  # In the first map, the state "B" only exists in a tiny window of the youngest ages: no divergence age can be used
  maps <- list(make_young_state_simmap(toy$phy, 1e-4), make_young_state_simmap(toy$phy, 40))
  graft <- function (map_list = maps, ...) graft_missing_taxa(toy$phy, mtl, simmap = as_multi_simmap(map_list, c("tiny", "wide")), allowed_states = "B", method = "BD_clades", max_age_redraws = 2, verbose = FALSE, ...)

  warnings_sent <- character(0)
  out <- withCallingHandlers(graft(nb_simulations = 2, seed = 1),
                             warning = function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_equal(names(out), "wide")
  expect_s3_class(out$wide, "multiPhylo")
  expect_match(warnings_sent, "Map 'tiny': none of the 2 simulation\\(s\\) could be obtained, the map is dropped.*No usable divergence age", all = FALSE)
  expect_equal(attr(out, "failed_maps")$map, "tiny")
  expect_match(attr(out, "failed_maps")$reason, "No usable divergence age")

  # With one simulation per map
  out1 <- suppressWarnings(graft(nb_simulations = 1, seed = 1))
  expect_equal(names(out1), "wide")
  expect_s3_class(out1$wide, "phylo")

  # No map gives a simulation
  both_tiny <- list(make_young_state_simmap(toy$phy, 1e-4), make_young_state_simmap(toy$phy, 1e-5))
  expect_error(suppressWarnings(graft(both_tiny, nb_simulations = 2, seed = 1)), "None of the 2 maps gave a simulation.*tiny.*wide")
})

test_that("graft_missing_clades() with a multiSimmap drops the map without any simulation, with a warning, and keeps the others",
{
  toy <- toy_setup()
  mcl <- list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = toy$B, M2 = toy$B))
  bad <- set_branch_states(make_toy_simmap(toy$phy), toy$B, c("A", "C"), c(1, 1e-7))
  good <- set_branch_states(make_toy_simmap(toy$phy), toy$B, c("A", "C"), c(1, 1))
  graft <- function (...) graft_missing_clades(toy$phy, mcl, simmap = as_multi_simmap(list(bad, good), c("bad", "good")), allowed_states = "C", max_age_redraws = 2, verbose = FALSE, ...)
  warnings_sent <- character(0)
  out <- withCallingHandlers(graft(nb_simulations = 2, seed = 1),
                             warning = function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") })
  # The warnings about the inputs name the map too (once)
  expect_equal(sum(grepl("^Map 'bad': 1 missing clade\\(s\\) have less than 10%", warnings_sent)), 1)
  expect_match(warnings_sent, "Map 'bad': none of the 2 simulation\\(s\\) could be obtained.*No usable divergence time", all = FALSE)
  expect_equal(names(out), "good")
  expect_equal(attr(out, "failed_maps")$map, "bad")
})

test_that("the progress messages state the map and the simulation",
{
  toy <- toy_setup()
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  lines <- capture.output(invisible(suppressWarnings(graft_missing_clades(toy$phy, maps_missing_clades(toy), simmap = as_multi_simmap(maps, c("a", "b")), allowed_states = c("A", "B"),
                                                                           nb_simulations = 2, verbose = TRUE, seed = 3))))
  expect_true(any(grepl("Map a .*Start grafting process", lines)))
  expect_true(any(grepl("Map b - Simulation n[^ ]*2 - Start grafting process", lines)))
})
