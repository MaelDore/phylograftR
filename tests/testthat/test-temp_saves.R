
##### Test suite for the temporary saves and the resume of the grafting functions #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Toy lists of missing taxa / missing clades, and multiSimmap (same as in test-multiple_simmaps.R)
ts_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}
ts_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}
ts_multi_simmap <- function (maps, names = NULL)
{
  class(maps) <- c("multiSimmap", "multiPhylo")
  if (!is.null(names)) { names(maps) <- names }
  maps
}

## Run an expression in which the grafting 'crashes' after 'n' grafts (the function that enforces the ultrametry of each new tip or clade is called once per graft)
crash_after <- function (n, expr)
{
  count <- 0
  original <- enforce_ultrametry_new_tips
  testthat::with_mocked_bindings(expr, enforce_ultrametry_new_tips = function (...)
  {
    count <<- count + 1
    if (count > n) { stop("simulated crash") }
    original(...)
  }, .package = "phylograftR")
}

## Temporary directory and path
new_temp_path <- function (name = "temp.rds")
{
  dir <- tempfile("tempsaves")
  dir.create(dir)
  file.path(dir, name)
}

## Skip if the packages are missing, or if the workers cannot load the package phylograftR
skip_if_no_workers_ts <- function ()
{
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  ok <- tryCatch(isTRUE(future::value(future::future(requireNamespace("phylograftR", quietly = TRUE)))), error = function (e) FALSE)
  if (!ok) { testthat::skip("the workers cannot load the package phylograftR (install it first)") }
}


# ---- Checks of the arguments ----

test_that("check_save_temp() validates save_temp, save_temp_every and resume", {
  expect_true(check_save_temp(TRUE, 5, FALSE))
  expect_true(check_save_temp(TRUE, NULL, FALSE))
  expect_true(check_save_temp(FALSE, 10, FALSE))
  expect_true(check_save_temp(TRUE, 1, TRUE))

  expect_error(check_save_temp(NA, 10, FALSE), "'save_temp'")
  expect_error(check_save_temp("yes", 10, FALSE), "'save_temp'")
  expect_error(check_save_temp(TRUE, 10, NA), "'resume'")
  expect_error(check_save_temp(TRUE, 0, FALSE), "positive whole number")
  expect_error(check_save_temp(TRUE, -3, FALSE), "positive whole number")
  expect_error(check_save_temp(TRUE, 2.5, FALSE), "positive whole number")
  expect_error(check_save_temp(TRUE, c(2, 3), FALSE), "positive whole number")
  expect_error(check_save_temp(TRUE, NA, FALSE), "positive whole number")
  expect_error(check_save_temp(TRUE, Inf, FALSE), "positive whole number")
  expect_error(check_save_temp(TRUE, "5", FALSE), "positive whole number")
  expect_error(check_save_temp(FALSE, 10, TRUE), "needs 'save_temp = TRUE'")
})

test_that("The grafting functions check save_temp_every, before any grafting", {
  toy <- toy_setup()
  expect_error(graft_missing_taxa(toy$phy, ts_taxa(toy), save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 0, verbose = FALSE), "positive whole number")
  expect_error(graft_missing_clades(toy$phy, ts_clades(toy), save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1.5, verbose = FALSE), "positive whole number")
  expect_error(graft_missing_taxa(toy$phy, ts_taxa(toy), resume = TRUE, verbose = FALSE), "needs 'save_temp = TRUE'")
  expect_error(graft_missing_clades(toy$phy, ts_clades(toy), resume = TRUE, verbose = FALSE), "needs 'save_temp = TRUE'")
})

test_that("A save_temp_every larger than the number of taxa / groups gives one warning, NULL does not", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  path <- new_temp_path()

  # Uniform: 5 taxa
  expect_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 6, temp_file_path = path, verbose = FALSE),
                 "larger than the number of taxa to graft \\(5\\).*once fully grafted")
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 5, temp_file_path = path, verbose = FALSE))
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = NULL, temp_file_path = path, verbose = FALSE))
  # Not sent without save_temp
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "uniform", save_temp = FALSE, save_temp_every = 100, verbose = FALSE))

  # BD_clades: the unit is the group of taxa (same target clade, stem availability and allowed states)
  n_groups <- length(group_and_order_clades(toy$phy, mtl$Target_clade, include_target_clade_stem = rep(TRUE, 5))$group_clade_tips)
  expect_warning(graft_missing_taxa(toy$phy, mtl, method = "BD_clades", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = n_groups + 1, temp_file_path = path, verbose = FALSE),
                 paste0("larger than the number of groups of taxa to graft \\(", n_groups, "\\)"))
  expect_no_warning(graft_missing_taxa(toy$phy, mtl, method = "BD_clades", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = n_groups, temp_file_path = path, verbose = FALSE))

  # Missing clades
  mcl <- ts_clades(toy)
  n_clade_groups <- length(group_and_order_clades(toy$phy, mcl$Sister_clade)$group_clade_tips)
  expect_warning(graft_missing_clades(toy$phy, mcl, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = n_clade_groups + 1, temp_file_path = path, verbose = FALSE),
                 paste0("larger than the number of missing clades to graft \\(", n_clade_groups, "\\)"))
  expect_no_warning(graft_missing_clades(toy$phy, mcl, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = NULL, temp_file_path = path, verbose = FALSE))
  unlink(dirname(path), recursive = TRUE)
})

test_that("The warning on save_temp_every is sent once, even with several simulations and maps", {
  toy <- toy_setup()
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  path <- new_temp_path()
  warnings_sent <- character(0)
  withCallingHandlers(
    graft_missing_taxa(toy$phy, ts_taxa(toy), simmap = ts_multi_simmap(maps, c("a", "b")), method = "uniform", nb_simulations = 2, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 50,
                       temp_file_path = path, verbose = FALSE, parallel = FALSE),
    warning = function (w) { warnings_sent <<- c(warnings_sent, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_equal(sum(grepl("save_temp_every", warnings_sent)), 1)
  unlink(dirname(path), recursive = TRUE)
})


# ---- Names of the files ----

test_that("The path is made absolute and the default is in the working directory", {
  expect_equal(resolve_temp_file_path(NULL), file.path(normalizePath(".", winslash = "/"), "temp_phy_from_grafting.rds"))
  expect_equal(resolve_temp_file_path("temp.rds"), file.path(normalizePath(".", winslash = "/"), "temp.rds"))
  absolute <- file.path(normalizePath(tempdir(), winslash = "/"), "x.rds")
  expect_equal(resolve_temp_file_path(absolute), absolute)
})

test_that("temp_map_tags() keeps distinct the maps that give the same file name once the prefix 'map_' is not repeated", {
  tags <- temp_map_tags(c("map_1", "1", "x"))
  expect_equal(temp_file_path_simulation("d/t.rds", 1, multiple = FALSE, map_name = "x", map_tag = tags[["x"]]), "d/t_map_x.rds")
  paths <- vapply(c("map_1", "1"), function (m) temp_file_path_simulation("d/t.rds", 1, multiple = FALSE, map_name = m, map_tag = tags[[m]]), character(1))
  expect_false(anyDuplicated(paths) > 0)
  expect_equal(unname(temp_file_path_simulation("d/t.rds", 1, multiple = FALSE, map_name = "map_2", map_tag = temp_map_tags("map_2")[["map_2"]])), "d/t_map_2.rds")
})

test_that("temp_map_tags() sanitises the names of the maps and keeps them distinct", {
  expect_null(temp_map_tags(NULL))
  expect_equal(unname(temp_map_tags(c("a", "hyp A/1"))), c("a", "hyp_A_1"))
  # Two names with the same tag
  tags <- temp_map_tags(c("a b", "a/b", "c"))
  expect_equal(unname(tags), c("a_b-1", "a_b-2", "c"))
  expect_equal(names(tags), c("a b", "a/b", "c"))
  # Used in the name of the file
  expect_equal(temp_file_path_simulation("d/t.rds", 2, multiple = TRUE, map_name = "a b", map_tag = tags[["a b"]]), "d/t_map_a_b-1_sim_2.rds")
})

test_that("Maps with the same sanitised name get distinct files, map first then simulation", {
  toy <- toy_setup()
  maps <- lapply(1:3, function (k) make_toy_simmap(toy$phy, seed = k))
  path <- new_temp_path()
  out <- graft_missing_taxa(toy$phy, ts_taxa(toy), simmap = ts_multi_simmap(maps, c("a b", "a/b", "c")), method = "uniform", nb_simulations = 2,
                            save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = NULL, temp_file_path = path, verbose = FALSE, parallel = FALSE)
  expect_setequal(list.files(dirname(path)), c("temp_map_a_b-1_sim_1.rds", "temp_map_a_b-1_sim_2.rds", "temp_map_a_b-2_sim_1.rds", "temp_map_a_b-2_sim_2.rds",
                                               "temp_map_c_sim_1.rds", "temp_map_c_sim_2.rds"))
  expect_named(out, c("a b", "a/b", "c"))
  unlink(dirname(path), recursive = TRUE)
})


# ---- Content of the files ----

test_that("Without save_temp_every, only the completed phylogeny is saved, and with it the partial progress is saved first", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)

  # Only completed (NULL): nothing is saved if the run is interrupted
  path <- new_temp_path()
  res <- try(crash_after(3, graft_missing_taxa(toy$phy, mtl, method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = NULL, temp_file_path = path, verbose = FALSE)), silent = TRUE)
  expect_s3_class(res, "try-error")
  expect_length(list.files(dirname(path)), 0)

  # Progress saved every 2 taxa: after 3 grafts the file holds the progress after 2 taxa, as a plain phylo
  res <- try(crash_after(3, graft_missing_taxa(toy$phy, mtl, method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 2, temp_file_path = path, verbose = FALSE)), silent = TRUE)
  expect_s3_class(res, "try-error")
  partial <- readRDS(path)
  expect_s3_class(partial, "phylo")
  expect_equal(ape::Ntip(partial), ape::Ntip(toy$phy) + 2)
  state <- attr(partial, "phylograftR_temp")
  expect_equal(state$status, "partial")
  expect_equal(state$progress, list(done = 2, total = 5, unit = "taxa"))
  expect_false(any(grepl("^bb", partial$node.label)))
  unlink(dirname(path), recursive = TRUE)
})

test_that("The completed phylogeny replaces the progress in the file, for any save_temp_every", {
  toy <- toy_setup()
  for (every in list(1, 2, NULL, 50))
  {
    path <- new_temp_path()
    out <- suppressWarnings(graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = every, temp_file_path = path, verbose = FALSE))
    saved <- readRDS(path)
    expect_equal(attr(saved, "phylograftR_temp")$status, "complete")
    attr(saved, "phylograftR_temp") <- NULL
    # The attribute "seed" is added to the output after the save
    attr(out, "seed") <- NULL
    expect_identical(saved, out)
    expect_false(any(grepl("\\.part$", list.files(dirname(path)))))
    unlink(dirname(path), recursive = TRUE)
  }
})

test_that("graft_missing_clades() saves its progress and the completed phylogeny the same way", {
  toy <- toy_setup()
  mcl <- ts_clades(toy)
  path <- new_temp_path()
  res <- try(crash_after(2, graft_missing_clades(toy$phy, mcl, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, verbose = FALSE)), silent = TRUE)
  expect_s3_class(res, "try-error")
  state <- attr(readRDS(path), "phylograftR_temp")
  expect_equal(state$status, "partial")
  expect_equal(state$progress$unit, "missing clades")
  expect_equal(state$progress$done, 2)
  out <- graft_missing_clades(toy$phy, mcl, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, verbose = FALSE)
  expect_equal(attr(readRDS(path), "phylograftR_temp")$status, "complete")
  unlink(dirname(path), recursive = TRUE)
})

test_that("A failing write stops with an informative message and never leaves a side file", {
  path <- file.path(tempfile("nodir"), "x.rds")
  expect_error(suppressWarnings(write_temp_file(1, path)), "could not be written")
})


# ---- Resume: the resumed phylogenies are identical to the uninterrupted ones ----

# Crash after 'n_crash' grafts, then resume with the same arguments: the result must be identical to the uninterrupted run
expect_resume_identical <- function (graft, n_crash, every, ...)
{
  path <- new_temp_path()
  reference <- graft(save_temp = FALSE, ...)
  crashed <- try(crash_after(n_crash, graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = every, temp_file_path = path, ...)), silent = TRUE)
  expect_s3_class(crashed, "try-error")
  resumed <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = every, temp_file_path = path, resume = TRUE, ...)
  expect_identical(resumed, reference)
  unlink(dirname(path), recursive = TRUE)
  invisible(resumed)
}

test_that("graft_missing_taxa() resumed after an interruption gives the same phylogeny, uniform and BD_clades", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  for (method in c("uniform", "BD_clades"))
  {
    graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = method, seed = 11, verbose = FALSE, ...))
    expect_resume_identical(graft, n_crash = 3, every = 1)
    expect_resume_identical(graft, n_crash = 4, every = 2)
  }
})

test_that("Resuming works with the total richness of the clades (sampling fractions of the BD models)", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  total <- list(Target_clade = list(toy$A, toy$B, toy$C, toy$s), total_richness = c(toy$nA + 12, toy$nB + 4, NA, NA))
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "BD_clades", total_richness_per_clades = total, seed = 3, verbose = FALSE, ...))
  resumed <- expect_resume_identical(graft, n_crash = 3, every = 1)
  # The realized counts of the clades are the ones of the original run
  expect_identical(extract_sampling_info(resumed, total), extract_sampling_info(graft(save_temp = FALSE), total))
})

test_that("Resuming works with a simmap: updated simmap and states are identical", {
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy, seed = 2)
  mtl <- ts_taxa(toy)
  for (method in c("uniform", "BD_clades"))
  {
    graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = simmap, allowed_states = c("A", "B", "C"), method = method, seed = 4, verbose = FALSE, ...))
    resumed <- expect_resume_identical(graft, n_crash = 3, every = 1)
    expect_s3_class(resumed$simmap, "simmap")
    expect_setequal(resumed$grafted_tips[resumed$grafted_tips] |> names(), mtl$Taxon_name)
  }
})

test_that("graft_missing_clades() resumed after an interruption gives the same phylogeny, with and without simmap", {
  toy <- toy_setup()
  mcl <- ts_clades(toy)
  graft <- function (...) suppressWarnings(graft_missing_clades(toy$phy, mcl, seed = 6, verbose = FALSE, ...))
  expect_resume_identical(graft, n_crash = 2, every = 1)
  simmap <- make_toy_simmap(toy$phy, seed = 3)
  graft_states <- function (...) suppressWarnings(graft_missing_clades(toy$phy, mcl, simmap = simmap, seed = 6, verbose = FALSE, ...))
  resumed <- expect_resume_identical(graft_states, n_crash = 2, every = 1)
  expect_s3_class(resumed$simmap, "simmap")
})

test_that("Resuming from a save at the very end of the loop only runs the final steps", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 2, verbose = FALSE, ...))
  path <- new_temp_path()
  reference <- graft(save_temp = FALSE)
  # All the 5 grafts are done (the crash happens at the following call, the final check of the ultrametry that is not mocked here): save the progress by hand after 5 taxa
  res <- try(crash_after(5, graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 5, temp_file_path = path)), silent = TRUE)
  # The run completed (no further call): replace the complete file by a partial one saved after the last taxon
  partial_done <- attr(readRDS(path), "phylograftR_temp")$status
  expect_true(partial_done %in% c("complete", "partial"))
  resumed <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 5, temp_file_path = path, resume = TRUE)
  expect_identical(resumed, reference)
  unlink(dirname(path), recursive = TRUE)
})

test_that("Resuming a run with several simulations and several maps: complete files are reused, unfinished ones continued, missing ones started", {
  toy <- toy_setup()
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  mtl <- ts_taxa(toy)
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = ts_multi_simmap(maps, c("a", "b")), allowed_states = c("A", "B", "C"),
                                                              method = "uniform", nb_simulations = 2, seed = 9, verbose = FALSE, parallel = FALSE, ...))
  path <- new_temp_path()
  reference <- graft(save_temp = FALSE)
  # 4 tasks of 5 grafts: crash during the second task (map a, simulation 2)
  crashed <- try(crash_after(7, graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path)), silent = TRUE)
  expect_s3_class(crashed, "try-error")
  status <- vapply(list.files(dirname(path), full.names = TRUE), function (f) attr(readRDS(f), "phylograftR_temp")$status, character(1))
  expect_equal(unname(status[basename(names(status)) == "temp_map_a_sim_1.rds"]), "complete")
  expect_equal(unname(status[basename(names(status)) == "temp_map_a_sim_2.rds"]), "partial")
  expect_false("temp_map_b_sim_1.rds" %in% basename(names(status)))

  resumed <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE)
  expect_identical(resumed, reference)
  expect_setequal(list.files(dirname(path)), c("temp_map_a_sim_1.rds", "temp_map_a_sim_2.rds", "temp_map_b_sim_1.rds", "temp_map_b_sim_2.rds"))
  expect_true(all(vapply(list.files(dirname(path), full.names = TRUE), function (f) attr(readRDS(f), "phylograftR_temp")$status == "complete", logical(1))))

  # A second resume reuses everything: no graft is done
  count <- 0
  again <- testthat::with_mocked_bindings(graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE),
                                          enforce_ultrametry_new_tips = function (...) { count <<- count + 1; stop("must not be called") }, .package = "phylograftR")
  expect_identical(again, reference)
  expect_equal(count, 0)
  unlink(dirname(path), recursive = TRUE)
})

test_that("Resuming with force_simulations uses the attempt numbers in the names of the files", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "uniform", nb_simulations = 2, force_simulations = TRUE, seed = 5, verbose = FALSE, parallel = FALSE, ...))
  path <- new_temp_path()
  reference <- graft(save_temp = FALSE)
  graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path)
  expect_true(all(c("temp_sim_1.rds", "temp_sim_2.rds") %in% list.files(dirname(path))))
  resumed <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE)
  expect_identical(resumed, reference)
  unlink(dirname(path), recursive = TRUE)
})

test_that("Resuming without any file starts from scratch, and without a seed the run still completes", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  path <- new_temp_path()
  graft <- function (seed, ...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = seed, verbose = FALSE, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 2, temp_file_path = path, ...))
  fresh <- graft(seed = 8, resume = TRUE)
  expect_identical(drop_simmap(fresh), drop_simmap(graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 8, verbose = FALSE)))
  unlink(path)

  # No seed: the state of the random number generator is not saved, a resumed run is valid but not reproducible
  res <- try(crash_after(3, graft(seed = NULL)), silent = TRUE)
  expect_s3_class(res, "try-error")
  expect_null(attr(readRDS(path), "phylograftR_temp")$rng)
  resumed <- graft(seed = NULL, resume = TRUE)
  expect_equal(ape::Ntip(resumed), ape::Ntip(toy$phy) + 5)
  expect_setequal(resumed$tip.label[resumed$grafted_tips], mtl$Taxon_name)
  unlink(dirname(path), recursive = TRUE)
})

test_that("Resume restores the random number generator of the session after the run", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  path <- new_temp_path()
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 8, verbose = FALSE, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, ...))
  try(crash_after(3, graft()), silent = TRUE)
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  graft(resume = TRUE)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  unlink(dirname(path), recursive = TRUE)
})


# ---- Resume: inputs that differ, old or corrupted files ----

test_that("Resuming with different inputs stops with an informative message, and resume = FALSE overwrites", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  path <- new_temp_path()
  graft <- function (mtl = mtl, seed = 8, method = "uniform", ...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = method, seed = seed, verbose = FALSE, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, ...))
  mtl_main <- ts_taxa(toy)
  try(crash_after(3, graft(mtl_main)), silent = TRUE)

  expect_error(graft(mtl_main, seed = 9, resume = TRUE), "different inputs \\(other arguments")
  expect_error(graft(mtl_main, seed = 9, resume = TRUE), "Seed of the simulation.*previous = [0-9]+, new = [0-9]+")
  expect_error(graft(mtl_main, method = "BD_clades", resume = TRUE), "different inputs \\(other arguments")
  expect_error(graft(mtl_main, method = "BD_clades", resume = TRUE), "Argument 'method': previous = uniform, new = BD_clades")
  mtl_other <- mtl_main
  mtl_other$Taxon_name[1] <- "Z1"
  names(mtl_other$Target_clade)[1] <- "Z1"; names(mtl_other$Excluded_subclades)[1] <- "Z1"
  expect_error(graft(mtl_other, resume = TRUE), "different inputs \\(the list of missing taxa")
  expect_error(graft(mtl_other, resume = TRUE), "only in the interrupted run \\(1\\): 'X1'")
  expect_error(graft(mtl_other, resume = TRUE), "only in the current call \\(1\\): 'Z1'")
  mtl_constraint <- mtl_main
  mtl_constraint$include_target_clade_stem <- stats::setNames(rep(FALSE, length(mtl_main$Taxon_name)), mtl_main$Taxon_name)
  expect_error(graft(mtl_constraint, resume = TRUE), "items of the list that differ: include_target_clade_stem")
  other_phy <- toy$phy
  other_phy$edge.length[1] <- other_phy$edge.length[1] * 0.5
  expect_error(suppressWarnings(graft_missing_taxa(other_phy, mtl_main, method = "uniform", seed = 8, verbose = FALSE, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE)),
               "(different inputs|not ultrametric)")
  simmap <- make_toy_simmap(toy$phy, seed = 1)
  expect_error(graft(mtl_main, simmap = simmap, allowed_states = c("A", "B", "C"), resume = TRUE), "different inputs \\(.*the simmap")

  # Without 'resume', the files are overwritten as before
  expect_no_error(graft(mtl_main, seed = 9, resume = FALSE))
  expect_equal(attr(readRDS(path), "phylograftR_temp")$status, "complete")
  unlink(dirname(path), recursive = TRUE)
})

test_that("check_temp_signature() describes what differs, with the previous and the new values", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  sig <- function (phy = toy$phy, mtl = mtl, tracker = NULL, options = list(method = "uniform", min_tips = 5, sim_seed = 1, richness = NULL))
    temp_signature(phy = phy, missing_taxa_list = mtl, tracker = tracker, options = options)
  expect_true(check_temp_signature(sig(mtl = mtl), sig(mtl = mtl), "f.rds"))
  # Backbone: tips and branch lengths
  smaller <- ape::drop.tip(toy$phy, toy$phy$tip.label[1])
  expect_error(check_temp_signature(sig(mtl = mtl), sig(phy = smaller, mtl = mtl), "f.rds"), "Backbone: number of tips, previous = [0-9]+, new = [0-9]+")
  longer <- toy$phy; longer$edge.length <- longer$edge.length * 2
  expect_error(check_temp_signature(sig(mtl = mtl), sig(phy = longer, mtl = mtl), "f.rds"), "Backbone: different branch lengths")
  expect_error(check_temp_signature(sig(mtl = mtl), sig(phy = longer, mtl = mtl), "f.rds"), "age of the root, previous")
  # Taxa: different constraints for one taxon
  mtl2 <- mtl; mtl2$Target_clade$X3 <- toy$A
  expect_error(check_temp_signature(sig(mtl = mtl), sig(mtl = mtl2), "f.rds"), "different constraints for 1 taxa: 'X3'")
  # Arguments, with their previous and new values
  expect_error(check_temp_signature(sig(mtl = mtl), sig(mtl = mtl, options = list(method = "BD_clades", min_tips = 8, sim_seed = 1, richness = NULL)), "f.rds"),
               "Argument 'method': previous = uniform, new = BD_clades.*Argument 'min_tips': previous = 5, new = 8")
  # A file saved without summary (older format) still stops, with the generic message
  old <- sig(mtl = mtl); old$details <- NULL
  expect_error(check_temp_signature(old, sig(mtl = mtl, options = list(method = "BD_clades", min_tips = 5, sim_seed = 1, richness = NULL)), "f.rds"), "different inputs \\(other arguments")
})

test_that("A file that was not saved by this version, or is corrupted, is ignored with a warning", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  path <- new_temp_path()
  graft <- function (...) graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 8, verbose = FALSE, save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE, ...)
  reference <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 8, verbose = FALSE))

  # Old format: a plain phylo
  saveRDS(toy$phy, path)
  expect_warning(out <- graft(), "does not hold the state needed to resume")
  expect_identical(out, reference)

  # Corrupted
  writeLines("not an rds file", path)
  expect_warning(out <- graft(), "could not be read")
  expect_identical(out, reference)
  unlink(dirname(path), recursive = TRUE)
})

test_that("read_resume_state() returns NULL without a file, the complete phylogeny, or the state of a partial file", {
  toy <- toy_setup()
  path <- new_temp_path()
  signature <- list(backbone = "a", taxa = "b", states = NULL, options = "c")
  expect_null(read_resume_state(path, signature))

  phy <- toy$phy
  phy$node.label <- paste0("bb", seq_len(phy$Nnode))
  labelled <- label_backbone_nodes(toy$phy)
  save_temp_progress(path = path, phy = labelled$phy, labelled_backbone = labelled, signature = signature, sim_seed = 1, done = 0, total = 5, unit = "taxa")
  resumed <- read_resume_state(path, signature)
  expect_equal(resumed$status, "partial")
  expect_identical(resumed$phy$node.label, labelled$phy$node.label)
  expect_null(attr(resumed$phy, "phylograftR_temp"))
  expect_false(is.null(resumed$state$rng))

  save_temp_complete(path, toy$phy, signature)
  expect_equal(read_resume_state(path, signature)$status, "complete")
  expect_error(read_resume_state(path, list(backbone = "x", taxa = "b", states = NULL, options = "c")), "the backbone phylogeny")
  unlink(dirname(path), recursive = TRUE)
})


# ---- Parallel computation ----

test_that("The files of a parallel run are written by the workers, and a resume in parallel reuses the completed ones", {
  skip_if_no_workers_ts()
  toy <- toy_setup()
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  mtl <- ts_taxa(toy)
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = ts_multi_simmap(maps, c("a", "b")), allowed_states = c("A", "B", "C"),
                                                              method = "uniform", nb_simulations = 2, seed = 9, verbose = FALSE, ...))
  reference <- graft(save_temp = FALSE, parallel = FALSE)

  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  # Relative path: made absolute in the main session
  withr_dir <- tempfile("parallel"); dir.create(withr_dir)
  path <- file.path(withr_dir, "temp.rds")
  out <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 2, temp_file_path = path, parallel = TRUE)
  expect_identical(out, reference)
  files <- list.files(withr_dir)
  expect_setequal(files, c("temp_map_a_sim_1.rds", "temp_map_a_sim_2.rds", "temp_map_b_sim_1.rds", "temp_map_b_sim_2.rds"))

  # Delete one completed file: a resume in parallel recomputes only this one, and reuses the others
  unlink(file.path(withr_dir, "temp_map_b_sim_1.rds"))
  resumed <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 2, temp_file_path = path, resume = TRUE, parallel = TRUE)
  expect_identical(resumed, reference)
  expect_setequal(list.files(withr_dir), files)
  unlink(withr_dir, recursive = TRUE)
})


# ---- Cleaning of the temporary files, and order of the arguments ----

test_that("check_save_temp() validates clean_temp_files", {
  expect_true(check_save_temp(TRUE, NULL, FALSE, clean_temp_files = FALSE))
  expect_error(check_save_temp(TRUE, NULL, FALSE, clean_temp_files = NA), "'clean_temp_files'")
  expect_error(check_save_temp(TRUE, NULL, FALSE, clean_temp_files = "yes"), "'clean_temp_files'")
  toy <- toy_setup()
  expect_error(graft_missing_taxa(toy$phy, ts_taxa(toy), save_temp = TRUE, clean_temp_files = 2, verbose = FALSE), "'clean_temp_files'")
  expect_error(graft_missing_clades(toy$phy, ts_clades(toy), save_temp = TRUE, clean_temp_files = NULL, verbose = FALSE), "'clean_temp_files'")
})

test_that("The temporary files are deleted by default once the run is completed, and only these files", {
  toy <- toy_setup()
  maps <- lapply(1:2, function (k) make_toy_simmap(toy$phy, seed = k))
  path <- new_temp_path()
  other <- file.path(dirname(path), "other.rds")
  saveRDS(1, other)
  # Single phylogeny
  expect_formals <- formals(graft_missing_taxa)
  expect_true(expect_formals$clean_temp_files)
  out <- graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", save_temp = TRUE, save_temp_every = 1, temp_file_path = path, verbose = FALSE)
  expect_s3_class(out, "phylo")
  expect_equal(list.files(dirname(path)), "other.rds")
  # Simulations x maps, taxa and clades
  out <- suppressWarnings(graft_missing_taxa(toy$phy, ts_taxa(toy), simmap = ts_multi_simmap(maps, c("a", "b")), allowed_states = c("A", "B", "C"), method = "uniform",
                                             nb_simulations = 2, save_temp = TRUE, temp_file_path = path, verbose = FALSE, parallel = FALSE))
  expect_equal(list.files(dirname(path)), "other.rds")
  out <- suppressWarnings(graft_missing_clades(toy$phy, ts_clades(toy), simmap = ts_multi_simmap(maps, c("a", "b")), nb_simulations = 2, save_temp = TRUE,
                                               save_temp_every = 1, temp_file_path = path, verbose = FALSE, parallel = FALSE))
  expect_equal(list.files(dirname(path)), "other.rds")
  # With force_simulations (up to 10 attempts per simulation)
  out <- graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", nb_simulations = 2, force_simulations = TRUE, save_temp = TRUE, temp_file_path = path, verbose = FALSE, parallel = FALSE)
  expect_equal(list.files(dirname(path)), "other.rds")
  # Nothing is created without save_temp
  graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", temp_file_path = path, verbose = FALSE)
  expect_equal(list.files(dirname(path)), "other.rds")
  unlink(dirname(path), recursive = TRUE)
})

test_that("clean_temp_files = FALSE keeps the files, and the message of the deletion is printed with verbose", {
  toy <- toy_setup()
  path <- new_temp_path()
  graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", save_temp = TRUE, clean_temp_files = FALSE, temp_file_path = path, verbose = FALSE)
  expect_true(file.exists(path))
  graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", save_temp = TRUE, nb_simulations = 2, temp_file_path = path, verbose = FALSE)
  expect_true(file.exists(path))  # The file of the first run is not one of this run (suffixes _sim_k): not touched
  expect_output(graft_missing_taxa(toy$phy, ts_taxa(toy), method = "uniform", save_temp = TRUE, temp_file_path = path, verbose = TRUE), "temporary file\\(s\\) deleted")
  expect_false(file.exists(path))
  unlink(dirname(path), recursive = TRUE)
})

test_that("A crashed run keeps its files, and the resumed run that completes deletes them", {
  toy <- toy_setup()
  mtl <- ts_taxa(toy)
  path <- new_temp_path()
  graft <- function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 3, verbose = FALSE, save_temp = TRUE, save_temp_every = 1, temp_file_path = path, ...))
  reference <- graft_missing_taxa(toy$phy, mtl, method = "uniform", seed = 3, verbose = FALSE)
  crashed <- try(crash_after(2, graft()), silent = TRUE)
  expect_s3_class(crashed, "try-error")
  expect_true(file.exists(path))
  resumed <- graft(resume = TRUE)
  expect_identical(resumed, reference)
  expect_length(list.files(dirname(path)), 0)
  unlink(dirname(path), recursive = TRUE)
})

test_that("temp_files_of_run() lists the files of all the maps and attempts", {
  tags <- temp_map_tags(c("a", "b"))
  expect_equal(basename(temp_files_of_run("/d/t.rds", NULL, NULL, 1, FALSE))[1], "t.rds")
  files <- temp_files_of_run("/d/t.rds", c("a", "b"), tags, 2, FALSE)
  expect_equal(basename(files[1:4]), c("t_map_a_sim_1.rds", "t_map_a_sim_2.rds", "t_map_b_sim_1.rds", "t_map_b_sim_2.rds"))
  expect_equal(basename(files[5:8]), paste0(basename(files[1:4]), ".part"))
  expect_length(temp_files_of_run("/d/t.rds", NULL, NULL, 2, TRUE), 2 * 20)
})

test_that("remove_temp_files() warns about the files it cannot delete", {
  path <- new_temp_path()
  saveRDS(1, path)
  with_mocked_bindings(
    expect_warning(remove_temp_files(path, NULL, NULL, 1, FALSE), "could not be deleted"),
    file.remove = function (...) FALSE, .package = "base")
  unlink(dirname(path), recursive = TRUE)
})

test_that("The arguments follow the blocks of the Details: inputs, method, constraints, richness, attributes, multi-run, saving and resuming, technical, verbose", {
  common_tail <- c("seed", "nb_simulations", "force_simulations", "parallel", "save_temp", "save_temp_every", "temp_file_path", "clean_temp_files", "resume", "ultrametric_tol", "verbose")
  expect_equal(names(formals(graft_missing_taxa)),
               c("phy", "missing_taxa_list", "method", "min_tips", "return_BD_model_summary", "include_target_clade_stem", "use_excluded_subclades", "exclude_subclades_stem", "total_richness_per_clades",
                 "simmap", "allowed_states", "treat_states_as_ranges", "return_simmap", "min_valid_fraction", "max_age_redraws", common_tail))
  expect_equal(names(formals(graft_missing_clades)),
               c("phy", "missing_taxa_list", "min_tips", "return_BD_model_summary", "total_richness_per_clades",
                 "simmap", "allowed_states", "treat_states_as_ranges", "return_simmap", "min_valid_fraction", "max_age_redraws", common_tail))
  for (f in list(graft_missing_taxa, graft_missing_clades))
  {
    expect_null(formals(f)$save_temp_every)
    expect_false(formals(f)$save_temp)
    expect_false(formals(f)$resume)
    expect_true(formals(f)$clean_temp_files)
  }
})
