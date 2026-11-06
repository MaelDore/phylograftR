
##### Test suite for the summary of the birth-death models (return_BD_model_summary and extract_BD_model_summary()) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Expected columns (without the optional Clade_name)
bd_columns <- c("Taxa", "n_taxa", "initial_clade_node", "initial_clade_tips", "initial_clade_n_tips", "initial_clade_crown_age", "enlarged",
                "fitted_clade_node", "fitted_clade_tips", "fitted_clade_n_tips", "fitted_clade_n_backbone_tips", "fitted_clade_crown_age", "min_tips",
                "b", "d", "d_fitted", "d_adjustment", "net_diversification", "relative_extinction", "rho", "n_pending", "n_unsampled")

## Crash after 'n' grafts (the function that enforces the ultrametry of each new tip or clade is called once per graft)
bd_crash_after <- function (n, expr)
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


test_that("graft_missing_taxa() reports one BD model per target clade, only if requested", {
  toy <- toy_setup()
  # B and C are disjoint and have at least 3 tips (B: 3 to 6 tips): B is smaller than min_tips = 8 so the model is fitted on an enclosing clade
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$s),
              Target_clade_name = c("clade_A", "clade_A", "clade_B", "tip_s"))

  default <- quiet_graft(graft_missing_taxa, toy$phy, mtl, seed = 1)
  expect_null(attr(default, "BD_model_summary"))
  expect_error(extract_BD_model_summary(default), "return_BD_model_summary = TRUE")

  g <- quiet_graft(graft_missing_taxa, toy$phy, mtl, min_tips = 8, return_BD_model_summary = TRUE, seed = 1)
  # The summary does not change the grafting
  expect_equal(g[c("edge", "edge.length", "tip.label")], default[c("edge", "edge.length", "tip.label")])

  info <- extract_BD_model_summary(g)
  expect_s3_class(info, "data.frame")
  expect_equal(names(info), c("map", "simulation", "seed", "Clade_name", bd_columns))
  expect_equal(attr(info, "type"), "taxa")
  expect_equal(nrow(info), 3) # One per target clade (A, B and the tip s), not per group nor per taxon
  expect_true(all(is.na(info$map)))
  expect_equal(info$simulation, rep(1, 3))
  # Rows in the order of the taxa of the list
  expect_equal(info$Clade_name, c("clade_A", "clade_B", "tip_s"))
  expect_equal(info$Taxa, c("X1, X2", "X3", "X4"))
  expect_equal(info$n_taxa, c(2, 1, 1))

  # The initial clades are located in the backbone
  expect_equal(info$initial_clade_node, c(toy$node_A, toy$node_B, match(toy$s, toy$phy$tip.label)))
  expect_equal(info$initial_clade_n_tips, c(toy$nA, toy$nB, 1))
  expect_equal(info$initial_clade_crown_age, c(toy$root_age - phytools::nodeheight(toy$phy, toy$node_A), toy$root_age - phytools::nodeheight(toy$phy, toy$node_B), NA))
  # A and B have fewer tips than min_tips = 8, but A has more: the clades are enlarged when they are too small
  expect_equal(info$enlarged, info$fitted_clade_node != info$initial_clade_node)
  expect_true(info$enlarged[3]) # A tip is always enlarged
  expect_true(all(info$fitted_clade_n_backbone_tips >= 8))
  expect_true(all(info$fitted_clade_n_tips >= info$fitted_clade_n_backbone_tips))
  # The fitted clades contain the initial ones
  for (k in 1:3)
  {
    fitted <- ape::extract.clade(toy$phy, info$fitted_clade_node[k])$tip.label
    initial <- if (info$initial_clade_n_tips[k] == 1) toy$s else ape::extract.clade(toy$phy, info$initial_clade_node[k])$tip.label
    expect_true(all(initial %in% fitted))
    expect_equal(length(fitted), info$fitted_clade_n_backbone_tips[k])
  }
  expect_equal(info$initial_clade_tips[2], paste(phylograftR:::node_to_tips(toy$phy, toy$node_B), collapse = ", "))

  # Rates
  expect_true(all(info$b > 0))
  expect_equal(info$net_diversification, info$b - info$d)
  expect_equal(info$relative_extinction, info$d / info$b)
  expect_true(all(info$d >= 0 & info$d <= 0.99 * info$b))
  expect_true(all(info$d_adjustment %in% c("none", "set_to_0", "capped")))
  expect_equal(info$d[info$d_adjustment == "none"], info$d_fitted[info$d_adjustment == "none"])
  # Sampling: nothing unsampled, the taxa still to graft are counted
  expect_equal(info$n_unsampled, rep(0, 3))
  expect_equal(info$rho, info$fitted_clade_n_tips / (info$fitted_clade_n_tips + info$n_pending))
  expect_true(all(info$n_pending >= 0))
})

test_that("The model of a clade that is not enlarged is the one fitted by fit_bd_on_clade() with its sampling fraction", {
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$A))
  total <- list(Target_clade = list(toy$A), total_richness = toy$nA + 2 + 5)
  g <- quiet_graft(graft_missing_taxa, toy$phy, mtl, total_richness_per_clades = total, min_tips = 5, return_BD_model_summary = TRUE, seed = 1)
  info <- extract_BD_model_summary(g)
  expect_equal(nrow(info), 1)
  expect_false(info$enlarged)
  expect_equal(info$n_unsampled, 5)
  expect_equal(info$n_pending, 2)
  expect_equal(info$rho, toy$nA / (toy$nA + 2 + 5))
  # Same rates as the direct fit on the backbone (the first fit is made before any graft of this clade)
  richness <- build_richness_info(toy$phy, mtl$Taxon_name, mtl$Target_clade, "target", total)
  fit <- fit_bd_on_clade(toy$phy, toy$node_A, min_tips = 5, richness_info = richness)
  expect_equal(c(info$b, info$d, info$rho), c(fit$b, fit$d, fit$rho))
})

test_that("graft_missing_clades() reports one BD model per missing clade", {
  toy <- toy_setup()
  mcl <- list(Taxon_name = c("M1", "M2", "M3", "M4"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s),
              Sister_clade_name = c("sis_B", "sis_B", "sis_C", "sis_s"))
  expect_null(attr(quiet_graft(graft_missing_clades, toy$phy, mcl, seed = 1), "BD_model_summary"))
  g <- quiet_graft(graft_missing_clades, toy$phy, mcl, return_BD_model_summary = TRUE, seed = 1)
  info <- extract_BD_model_summary(g)
  expect_equal(names(info), c("map", "simulation", "seed", "Clade_name", bd_columns))
  expect_equal(attr(info, "type"), "clades")
  expect_equal(nrow(info), 3)
  expect_equal(info$Clade_name, c("sis_B", "sis_C", "sis_s"))
  expect_equal(info$Taxa, c("M1, M2", "M3", "M4"))
  expect_equal(info$n_taxa, c(2, 1, 1))
  expect_equal(info$initial_clade_n_tips, c(toy$nB, toy$nC, 1))
  expect_true(info$enlarged[3])
  expect_true(all(info$fitted_clade_n_backbone_tips >= 5))
  # Without names, no Clade_name column
  mcl_no_names <- mcl; mcl_no_names$Sister_clade_name <- NULL
  info_no_names <- extract_BD_model_summary(quiet_graft(graft_missing_clades, toy$phy, mcl_no_names, return_BD_model_summary = TRUE, seed = 1))
  expect_equal(names(info_no_names), c("map", "simulation", "seed", bd_columns))
  expect_equal(info_no_names, info[, names(info) != "Clade_name"], ignore_attr = TRUE)
})

test_that("The BD summary has one table per simulation and map, and the same models in all of them (different rates)", {
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$C))
  g <- quiet_graft(graft_missing_taxa, toy$phy, mtl, return_BD_model_summary = TRUE, nb_simulations = 3, seed = 4)
  info <- extract_BD_model_summary(g)
  expect_equal(nrow(info), 6)
  expect_equal(info$simulation, rep(1:3, each = 2))
  expect_false(anyNA(info$seed))
  expect_equal(info$seed, rep(unname(attr(g, "seeds")), each = 2))
  # The models are fitted on the same clades in all the simulations
  expect_equal(info$fitted_clade_node[info$simulation == 2], info$fitted_clade_node[info$simulation == 1])
  # Subsetting a phylogeny keeps its table
  expect_equal(nrow(attr(g[[2]], "BD_model_summary")), 2)
  # Flattening keeps the tables
  expect_equal(nrow(extract_BD_model_summary(flatten_grafted_phylogenies(g))), 6)

  # Several maps
  maps <- structure(list(m1 = make_toy_simmap(toy$phy, seed = 1), m2 = make_toy_simmap(toy$phy, seed = 2)), class = c("multiSimmap", "multiPhylo"))
  g_maps <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = maps, return_BD_model_summary = TRUE, nb_simulations = 2, seed = 4)
  info_maps <- extract_BD_model_summary(g_maps)
  expect_equal(nrow(info_maps), 2 * 2 * 2)
  expect_equal(unique(info_maps$map), c("m1", "m2"))
  expect_equal(info_maps$simulation[info_maps$map == "m2"], rep(1:2, each = 2))
})

test_that("return_BD_model_summary is checked, and requires BD_clades", {
  toy <- toy_setup()
  mtl <- list(Taxon_name = "X1", Target_clade = list(X1 = toy$A))
  mcl <- list(Taxon_name = "M1", Sister_clade = list(M1 = toy$B))
  expect_error(graft_missing_taxa(toy$phy, mtl, method = "uniform", return_BD_model_summary = TRUE, verbose = FALSE), "requires method = 'BD_clades'")
  expect_error(graft_missing_taxa(toy$phy, mtl, return_BD_model_summary = NA, verbose = FALSE), "must be TRUE or FALSE")
  expect_error(graft_missing_clades(toy$phy, mcl, return_BD_model_summary = "yes", verbose = FALSE), "must be TRUE or FALSE")
  expect_no_error(graft_missing_taxa(toy$phy, mtl, method = "uniform", return_BD_model_summary = FALSE, verbose = FALSE))
  expect_error(extract_BD_model_summary("not a phylogeny"))
})

test_that("The BD summary is the same after resuming an interrupted run", {
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"), Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s))
  mcl <- list(Taxon_name = c("M1", "M2", "M3", "M4"), Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
  expect_resume <- function (graft)
  {
    dir <- tempfile("bdsummary"); dir.create(dir); path <- file.path(dir, "temp.rds")
    on.exit(unlink(dir, recursive = TRUE))
    reference <- graft(save_temp = FALSE)
    crashed <- try(bd_crash_after(2, graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path)), silent = TRUE)
    expect_s3_class(crashed, "try-error")
    resumed <- graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE)
    expect_identical(attr(resumed, "BD_model_summary"), attr(reference, "BD_model_summary"))
    expect_identical(resumed, reference)
    # A resume asking for another summary setting is refused
    expect_error(graft(save_temp = TRUE, clean_temp_files = FALSE, save_temp_every = 1, temp_file_path = path, resume = TRUE, return_BD_model_summary = FALSE))
  }
  expect_resume(function (...) suppressWarnings(graft_missing_taxa(toy$phy, mtl, return_BD_model_summary = TRUE, seed = 11, verbose = FALSE, ...)))
  expect_resume(function (...) suppressWarnings(graft_missing_clades(toy$phy, mcl, return_BD_model_summary = TRUE, seed = 11, verbose = FALSE, ...)))
})
