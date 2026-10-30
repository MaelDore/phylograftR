
##### Test suite for plot_grafted() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Toy grafted trees, with and without simmap
plot_setup <- function ()
{
  toy <- toy_setup()
  mtl <- list(Taxon_name = c("X1", "X2"), Target_clade = list(X1 = toy$A, X2 = toy$B))
  simmap <- make_toy_simmap(toy$phy, seed = 3)
  list(toy = toy, mtl = mtl, simmap = simmap,
       g = quiet_graft(graft_missing_taxa, toy$phy, mtl),
       g_sm = quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap))
}

## Run a plot on a null device and return what is plotted (the plotting functions are called with the arguments they receive)
plot_quiet <- function (expr)
{
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(expr)
}


test_that("dash_map() cuts the pieces in dashes and gaps, preserving the total length",
{
  # Long piece: alternating dashes and gaps, starting and ending with a dash
  out <- phylograftR:::dash_map(c(A = 10), dash_length = 1, gap_state = ".gap")
  expect_equal(sum(out), 10)
  expect_equal(names(out)[1], "A")
  expect_equal(names(out)[length(out)], "A")
  expect_equal(length(out) %% 2, 1)
  expect_equal(names(out), rep(c("A", ".gap"), length.out = length(out)))
  expect_equal(length(out), 2 * floor((10 / 1 - 1) / 2) + 1)
  expect_true(all(abs(out - out[1]) < 1e-12)) # Regular dashes

  # Short piece (< 3 dashes): kept as it is
  expect_equal(phylograftR:::dash_map(c(A = 2.9), dash_length = 1, gap_state = ".gap"), c(A = 2.9))
  # Exactly 3 dashes: one gap
  expect_equal(names(phylograftR:::dash_map(c(A = 3), dash_length = 1, gap_state = ".gap")), c("A", ".gap", "A"))

  # Several pieces: each is cut separately and the order of the states is kept
  out <- phylograftR:::dash_map(c(B = 0.5, A = 6, C = 1), dash_length = 1, gap_state = ".gap")
  expect_equal(sum(out), 7.5)
  expect_equal(unique(names(out)[names(out) != ".gap"]), c("B", "A", "C"))
  expect_equal(sum(out[names(out) == "B"]), 0.5)
  expect_equal(sum(out[names(out) == "C"]), 1)
  expect_equal(sum(out[names(out) %in% c("A", ".gap")]), 6)
})


test_that("plot_grafted() plots plain phylogenies with the requested elements colored, and returns NULL invisibly",
{
  s <- plot_setup()
  for (w in c("both", "edges", "tips"))
  {
    expect_no_error(res <- plot_quiet(withVisible(plot_grafted(s$g, which = w, cex = 0.7))))
    expect_null(res$value)
    expect_false(res$visible)
  }
  expect_no_error(plot_quiet(plot_grafted(s$g, color = "blue", main = "Title")))
  # Arguments are passed to plot.phylo
  expect_no_error(plot_quiet(plot_grafted(s$g, type = "fan", show.tip.label = FALSE)))
  # The simmap is ignored if not requested
  expect_no_error(plot_quiet(plot_grafted(s$g_sm, use_simmap = FALSE, cex = 0.7)))
})


test_that("plot_grafted() sends the right colors to plot.phylo",
{
  s <- plot_setup()
  captured <- NULL
  # Replace temporarily the plot method for phylo to capture its arguments
  graphics_ns <- asNamespace("graphics")
  registerS3method("plot", "phylo", function (x, ...) { captured <<- list(x = x, ...); invisible(NULL) }, envir = graphics_ns)
  plot_grafted(s$g, color = "blue")
  expect_equal(captured$edge.color == "blue", unname(s$g$grafted_edges))
  expect_equal(captured$tip.color == "blue", unname(s$g$grafted_tips))
  expect_true(all(captured$edge.color[!s$g$grafted_edges] == "black"))

  plot_grafted(s$g, which = "edges")
  expect_true(all(captured$tip.color == "black"))
  expect_true(any(captured$edge.color == "red"))
  plot_grafted(s$g, which = "tips")
  expect_true(all(captured$edge.color == "black"))
  expect_true(any(captured$tip.color == "red"))

  # User colors override the defaults
  plot_grafted(s$g, edge.color = "green")
  expect_equal(captured$edge.color, "green")
  registerS3method("plot", "phylo", ape::plot.phylo, envir = graphics_ns) # Restore the original method
})


test_that("plot_grafted() plots simmaps, with dashed grafted edges and colored grafted tips",
{
  s <- plot_setup()
  for (w in c("both", "edges", "tips"))
  {
    expect_no_error(plot_quiet(plot_grafted(s$g_sm, which = w, fsize = 0.7, lwd = 2)))
  }
  expect_no_error(plot_quiet(plot_grafted(s$g_sm, simmap_colors = c(A = "lightblue", B = "orange", C = "grey"), color = "purple")))
  # No labels
  expect_no_error(plot_quiet(plot_grafted(s$g_sm, ftype = "off")))
  # Font types
  for (ft in c("reg", "b", "i", "bi")) { expect_no_error(plot_quiet(plot_grafted(s$g_sm, ftype = ft))) }
  # Explicit xlim
  expect_no_error(plot_quiet(plot_grafted(s$g_sm, xlim = c(0, 200))))

  # The simmap sent to plotSimmap: the grafted edges are dashed, and the other ones are untouched
  captured <- NULL
  original_plotSimmap <- phytools::plotSimmap
  testthat::local_mocked_bindings(plotSimmap = function (tree, ...) { captured <<- tree; original_plotSimmap(tree = tree, ...) }, .package = "phytools")
  # The grafted edges of the toy tree are short (they stay solid): lengthen one of them
  g_long <- s$g_sm
  grafted <- which(g_long$grafted_edges)
  other <- which(!g_long$grafted_edges)
  g_long$simmap$maps[[grafted[1]]] <- c(A = 4, B = 6)
  g_long$simmap$edge.length[grafted[1]] <- 10
  plot_quiet(plot_grafted(g_long, which = "both", ftype = "off"))
  expect_equal(captured$maps[other], g_long$simmap$maps[other])
  expect_equal(vapply(captured$maps, sum, numeric(1)), vapply(g_long$simmap$maps, sum, numeric(1)))
  expect_true(".phylograftR_gap" %in% names(captured$maps[[grafted[1]]]))
  expect_equal(captured$maps[[grafted[2]]], g_long$simmap$maps[[grafted[2]]]) # Short edge: solid
  expect_null(captured$mapped.edge)

  # With edges not requested, nothing is dashed
  plot_quiet(plot_grafted(g_long, which = "tips", ftype = "off"))
  expect_equal(captured$maps, g_long$simmap$maps)
})


test_that("plot_grafted() selects the phylogeny and the trait",
{
  s <- plot_setup()
  g_multi <- quiet_graft(graft_missing_taxa, s$toy$phy, s$mtl, nb_simulations = 3, seed = 1)
  g_multi_sm <- quiet_graft(graft_missing_taxa, s$toy$phy, s$mtl, simmap = s$simmap, nb_simulations = 3, seed = 1)
  expect_no_error(plot_quiet(plot_grafted(g_multi, index = 3)))
  expect_no_error(plot_quiet(plot_grafted(g_multi_sm, index = 2)))
  expect_error(plot_quiet(plot_grafted(g_multi, index = 4)), "only 3 phylogenies")

  # Several traits
  traits <- list(biogeo = s$simmap, elev = make_toy_simmap(s$toy$phy, seed = 2))
  g_t <- quiet_graft(graft_missing_taxa, s$toy$phy, s$mtl, simmap = traits)
  expect_no_error(plot_quiet(plot_grafted(g_t, trait = 1)))
  expect_no_error(plot_quiet(plot_grafted(g_t, trait = "elev")))
  expect_error(plot_quiet(plot_grafted(g_t, trait = 3)), "'trait' must be")
  expect_error(plot_quiet(plot_grafted(g_t, trait = "none")), "'trait' must be")
  expect_error(plot_quiet(plot_grafted(g_t, trait = c(1, 2))), "'trait' must be")

  captured <- NULL
  original_plotSimmap <- phytools::plotSimmap
  testthat::local_mocked_bindings(plotSimmap = function (tree, ...) { captured <<- tree; original_plotSimmap(tree = tree, ...) }, .package = "phytools")
  plot_quiet(plot_grafted(g_t, trait = "elev", which = "tips", ftype = "off"))
  expect_equal(captured$maps, g_t$simmap$elev$maps)
  plot_quiet(plot_grafted(g_t, trait = 1, which = "tips"))
  expect_equal(captured$maps, g_t$simmap$biogeo$maps)
})


test_that("plot_grafted() checks its inputs",
{
  s <- plot_setup()
  expect_error(plot_quiet(plot_grafted(s$g, color = "notacolor")), "'color' must be")
  expect_error(plot_quiet(plot_grafted(s$g, color = c("red", "blue"))), "'color' must be")
  expect_error(plot_quiet(plot_grafted(s$g, color = NA_character_)), "'color' must be")
  expect_error(plot_quiet(plot_grafted(s$g, which = "branches")), "should be one of")
  expect_error(plot_quiet(plot_grafted(s$g, use_simmap = NA)), "'use_simmap' must be")
  expect_error(plot_quiet(plot_grafted(s$g, index = 0)), "'index' must be")
  expect_error(plot_quiet(plot_grafted(s$g, index = 1.5)), "'index' must be")
  expect_error(plot_quiet(plot_grafted(s$g, index = c(1, 2))), "'index' must be")

  # Phylogeny without the grafted elements
  no_info <- s$g; no_info$grafted_tips <- NULL; no_info$grafted_edges <- NULL
  expect_error(plot_quiet(plot_grafted(no_info)), "identify_grafted")

  # Elements that do not match the tree (tips removed)
  dropped <- ape::drop.tip(s$g, s$g$tip.label[1])
  dropped$grafted_tips <- s$g$grafted_tips
  dropped$grafted_edges <- s$g$grafted_edges
  expect_error(plot_quiet(plot_grafted(dropped)), "do not match")

  # Simmap that does not match the tree
  bad_sm <- s$g_sm
  bad_sm$simmap$edge <- bad_sm$simmap$edge[nrow(bad_sm$simmap$edge):1, ]
  expect_error(plot_quiet(plot_grafted(bad_sm)), "does not match")
  expect_no_error(plot_quiet(plot_grafted(bad_sm, use_simmap = FALSE)))

  # Simmap colors
  expect_error(plot_quiet(plot_grafted(s$g_sm, simmap_colors = c("red", "blue", "green"))), "named character vector")
  expect_error(plot_quiet(plot_grafted(s$g_sm, simmap_colors = c(A = "red"))), "no color for the states")

  # Direction and type are not supported for simmaps
  expect_error(plot_quiet(plot_grafted(s$g_sm, direction = "upwards")), "rightwards phylograms")
  expect_error(plot_quiet(plot_grafted(s$g_sm, type = "fan")), "rightwards phylograms")
  expect_no_error(plot_quiet(plot_grafted(s$g_sm, direction = "rightwards", type = "phylogram")))
})
