##### Function to plot grafted phylogenies #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Plot a grafted phylogeny with the grafted taxa highlighted ####

#' @title Plot a grafted phylogeny with the grafted taxa highlighted
#'
#' @description Plot a phylogeny returned by the grafting functions, with the grafted branches and/or the grafted tips highlighted.
#'   If the phylogeny has a `$simmap` (grafting with a `simmap`), the updated stochastic character map is plotted instead of the plain phylogeny:
#'   the colors of the states are kept, the grafted branches are drawn dashed lines, and the grafted tips are labelled with the highlight color.
#'
#'   The function uses the plotting functions of `ape` and `phytools` (`plot.phylo()`, `phytools::plotSimmap()` and `ape::tiplabels()`):
#'   it only prepares the colors and the dashes, and other plotting arguments are passed to these functions.
#'
#' @param phy_grafted A grafted phylogeny (class `"phylo"`, with the elements `$grafted_tips` and `$grafted_edges`, see [`phylograftR::identify_grafted()`][identify_grafted]),
#'   a set of grafted phylogenies (class `"multiPhylo"`, e.g., several simulations), or a list of them with one element per map (output of the grafting functions with a `"multiSimmap"`).
#'   **Only one phylogeny is plotted**: see `index`.
#' @param which Character string. What to highlight: `"both"` (default) for the grafted branches and tips, `"edges"` for the grafted branches only, or `"tips"` for the grafted tips only.
#' @param color Character string. Color used to highlight the grafted elements. Default is `"red"`.
#'   For a plain phylogeny, the grafted branches and the grafted tip labels are drawn in this color.
#'   For a simmap, the branches keep the colors of their states and the grafted tip labels are drawn in this color.
#' @param use_simmap Logical. If `TRUE` (default) and the phylogeny has a `$simmap`, plot the updated simmap instead of the plain phylogeny. If there is no `$simmap`, the plain phylogeny is plotted.
#' @param simmap_colors Named character vector. Colors of the states of the simmap: the names are the states. If `NULL` (default), the current default palette is used.
#'   Only used if the simmap is plotted.
#' @param index Integer. Position of the phylogeny to plot among the phylogenies of `phy_grafted`. Default is `1`.
#'   For a `"multiPhylo"`, it is the position of the simulation. For a list of maps, the phylogenies of the first map come first, then those of the second map, etc.
#' @param trait Integer or character string. The trait to plot when the phylogeny has the simmaps of several traits (position or name of the trait). Default is `1`. Ignored with a single trait.
#' @param ... Other arguments passed to `plot()` (see `ape::plot.phylo()`) for a plain phylogeny, or to `phytools::plotSimmap()` for a simmap
#'   (e.g., `cex`, `show.tip.label`, `edge.width` for the plain phylogeny; `fsize`, `ftype`, `lwd` for the simmap).
#'   A title can be added with `main` (and `sub`).
#'
#' @details
#' ### Plain phylogeny
#'
#' The grafted branches (`which = "edges"` or `"both"`) and the labels of the grafted tips (`which = "tips"` or `"both"`) are drawn in `color`, the rest in black.
#' This is the equivalent of passing `edge.color` and `tip.color` to `plot()`, which can still be overridden through `...`.
#'
#' ### Simmap
#'
#' A branch of a simmap is drawn in the color of its state, so the grafted branches cannot be highlighted by a color without hiding their state.
#' They are **drawn as dashed lines** instead, with `which = "edges"` or `"both"`.
#' The labels of the grafted tips, with `which = "tips"` or `"both"`, are drawn in `color`.
#' Only rightwards phylograms (`direction = "rightwards"`) are supported for simmaps.
#'
#' @importFrom ape tiplabels
#' @importFrom phytools plotSimmap nodeHeights
#' @importFrom graphics par strwidth
#' @importFrom grDevices col2rgb hcl.colors
#' @importFrom utils modifyList
#'
#' @export
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa based on target clades and excluded subclades.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades of backbone clades.
#'  * [`phylograftR::identify_grafted()`][identify_grafted] to identify the grafted tips and branches of a phylogeny.
#'
#' @return `NULL`, invisibly. The function is called for its plot.
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' # List the missing taxa to graft
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")))
#'
#' # --- Plot a grafted phylogeny --- #
#'
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   verbose = FALSE, seed = 1234)
#'
#' # Grafted branches and tips are highlighted in red
#' plot_grafted(grafted_phy, cex = 0.7)
#'
#' # Highlight only the tips, in green
#' plot_grafted(grafted_phy, which = "tips", color = "limegreen", cex = 0.7)
#'
#' # --- Plot a grafted phylogeny with its simmap --- #
#'
#' # Simulate a simmap with three states on the backbone
#' Q <- matrix(0.04, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' diag(Q) <- -0.08
#' simmap <- phytools::sim.history(phy, Q, anc = "A", message = FALSE)
#'
#' # Graft the taxa, updating the simmap
#' grafted_phy_states <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                          simmap = simmap, verbose = FALSE, seed = 1234)
#'
#' # The simmap is plotted: grafted branches are dashed, grafted tips are labelled in red
#' plot_grafted(grafted_phy_states, simmap_colors = c(A = "lightblue", B = "orange", C = "grey"),
#'              fsize = 0.7, lwd = 3)
#'
#' # Plot the phylogeny without its simmap
#' plot_grafted(grafted_phy_states, use_simmap = FALSE, cex = 0.7)
#'
#' # --- Plot one of several simulations --- #
#'
#' grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                    nb_simulations = 3, verbose = FALSE, seed = 1234)
#' plot_grafted(grafted_phys, index = 2, cex = 0.7, main = "Simulation 2")
#'
#'

plot_grafted <- function (phy_grafted, # Output of the grafting functions: phylo, multiPhylo, or list of maps
                          which = c("both", "edges", "tips"), # What to highlight
                          color = "red", # Highlight color
                          use_simmap = TRUE, # Plot the simmap (if available)?
                          simmap_colors = NULL, # Named colors of the states of the simmap
                          index = 1, # Which phylogeny to plot
                          trait = 1, # Which trait to plot (simmaps of several traits)
                          ...) # Passed to plot() or phytools::plotSimmap()
{
  ### ---- Check input validity ---- ###

  which <- match.arg(which)
  if (!is.character(color) || length(color) != 1 || is.na(color) || inherits(try(grDevices::col2rgb(color), silent = TRUE), "try-error"))
  {
    stop("'color' must be a single valid color")
  }
  if (!is.logical(use_simmap) || length(use_simmap) != 1 || is.na(use_simmap)) { stop("'use_simmap' must be TRUE or FALSE") }
  if (!is.numeric(index) || length(index) != 1 || is.na(index) || index < 1 || index != round(index)) { stop("'index' must be a single positive integer") }

  ## Select the phylogeny
  listed <- list_grafted_phylogenies(phy_grafted, arg = "phy_grafted")
  if (index > length(listed$phylos)) { stop("'index' is ", index, " but 'phy_grafted' has only ", length(listed$phylos), " phylogenies") }
  phy <- listed$phylos[[index]]

  ## Grafted tips and branches
  if (is.null(phy$grafted_tips) || is.null(phy$grafted_edges))
  {
    stop("The phylogeny has no '$grafted_tips' and '$grafted_edges' elements. Use identify_grafted() to add them.")
  }
  if (length(phy$grafted_tips) != ape::Ntip(phy) || length(phy$grafted_edges) != nrow(phy$edge))
  {
    stop("'$grafted_tips' and '$grafted_edges' do not match the phylogeny (e.g., tips were removed). Use identify_grafted() to update them.")
  }
  highlight_edges <- which %in% c("both", "edges")
  highlight_tips <- which %in% c("both", "tips")
  dots <- list(...)

  ## The simmap to plot, if any
  simmap <- NULL
  if (use_simmap && !is.null(phy$simmap))
  {
    simmap <- phy$simmap
    if (!inherits(simmap, "simmap")) # Several traits: one simmap per trait
    {
      if (length(trait) != 1 || (is.numeric(trait) && (is.na(trait) || trait < 1 || trait > length(simmap) || trait != round(trait))) ||
          (is.character(trait) && !(trait %in% names(simmap))))
      {
        stop("'trait' must be the position or the name of one of the traits: ", paste(if (is.null(names(simmap))) seq_along(simmap) else names(simmap), collapse = ", "))
      }
      simmap <- simmap[[trait]]
    }
    if (!identical(simmap$edge, phy$edge))
    {
      stop("The '$simmap' of the phylogeny does not match the phylogeny (e.g., tips were removed). Use 'use_simmap = FALSE' to plot the phylogeny.")
    }
  }

  ### ---- Plain phylogeny ---- ###

  if (is.null(simmap))
  {
    edge_color <- rep("black", nrow(phy$edge))
    tip_color <- rep("black", ape::Ntip(phy))
    if (highlight_edges) { edge_color[phy$grafted_edges] <- color }
    if (highlight_tips) { tip_color[phy$grafted_tips] <- color }
    args <- utils::modifyList(list(x = phy, edge.color = edge_color, tip.color = tip_color), dots)
    do.call(graphics::plot, args)
    return(invisible(NULL))
  }

  ### ---- Simmap ---- ###

  if (!is.null(dots$direction) && dots$direction != "rightwards" || !is.null(dots$type) && dots$type != "phylogram")
  {
    stop("Only rightwards phylograms are supported for simmaps: do not use 'direction' or 'type' (or use 'use_simmap = FALSE').")
  }

  ## Colors of the states
  gap_state <- ".phylograftR_gap"
  states <- sort(unique(unlist(lapply(simmap$maps, names), use.names = FALSE)))
  if (is.null(simmap_colors))
  {
    simmap_colors <- stats::setNames(grDevices::hcl.colors(length(states), "Dark 2"), states)
  } else {
    if (!is.character(simmap_colors) || is.null(names(simmap_colors)))
    {
      stop("'simmap_colors' must be a named character vector (names = states)")
    }
    if (!all(states %in% names(simmap_colors)))
    {
      stop("'simmap_colors' has no color for the states: ", paste(setdiff(states, names(simmap_colors)), collapse = ", "))
    }
  }
  plot_colors <- c(simmap_colors[states], stats::setNames("transparent", gap_state))

  ## Dashed grafted branches: alternating pieces in the state of the branch and in a transparent state
  root_age <- max(phytools::nodeHeights(simmap)[, 2])
  if (highlight_edges)
  {
    simmap$maps[phy$grafted_edges] <- lapply(simmap$maps[phy$grafted_edges], dash_map, dash_length = 0.005 * root_age, gap_state = gap_state)
  }
  simmap$mapped.edge <- NULL

  ## Tip labels: drawn separately to color the grafted ones
  label_size <- if (is.null(dots$fsize)) 1 else dots$fsize
  label_type <- if (is.null(dots$ftype)) "reg" else dots$ftype
  draw_labels <- label_type != "off"
  dots[c("fsize", "ftype", "colors")] <- NULL
  # plotSimmap() has no title: draw it separately, and leave room for it in the margins
  plot_title <- dots[intersect(c("main", "sub"), names(dots))]
  dots[c("main", "sub")] <- NULL
  if (is.null(dots$mar) && length(plot_title) > 0)
  {
    dots$mar <- c(if (is.null(plot_title$sub)) 0.1 else 1.5, 0.1, if (is.null(plot_title$main)) 0.1 else 1.5, 0.1)
  }
  # Room for the labels at the right of the tree (the labels are not known by plotSimmap)
  if (draw_labels && is.null(dots$xlim))
  {
    label_width <- max(graphics::strwidth(simmap$tip.label, units = "inches", cex = label_size, font = c(reg = 1, b = 2, i = 3, bi = 4)[[label_type]]))
    label_fraction <- min(label_width / graphics::par("pin")[1], 0.6)
    dots$xlim <- c(0, root_age / (1 - label_fraction - 0.03))
  }
  old_par <- graphics::par(lend = "butt") # Dashes as long as the gaps
  on.exit(graphics::par(old_par), add = TRUE)
  do.call(phytools::plotSimmap, c(list(tree = simmap, colors = plot_colors, ftype = "off", fsize = label_size), dots))
  if (length(plot_title) > 0) { do.call(graphics::title, plot_title) }
  if (draw_labels)
  {
    ape::tiplabels(text = simmap$tip.label, tip = seq_len(ape::Ntip(simmap)), frame = "none", adj = c(0, 0.5), offset = 0.01 * root_age, cex = label_size,
                   font = c(reg = 1, b = 2, i = 3, bi = 4)[[label_type]], col = ifelse(highlight_tips & phy$grafted_tips, color, "black"))
  }
  invisible(NULL)
}


#' @title Cut the map of a branch in dashes
#'
#' @description Cut each piece of the map of a branch (one state) in alternating pieces of the same state and of a transparent state (the gap), starting and ending with a dash.
#'   The gaps are twice as long as the dashes. A piece that is too short to have a gap (shorter than three dashes) is kept as it is.
#'
#' @param map Named numeric vector. The pieces of the branch from the rootward node to the tipward node: the durations, named after the states.
#' @param dash_length Numeric. Approximate length of the dashes (the gaps are twice as long).
#' @param gap_state Character string. Name of the transparent state.
#'
#' @noRd
#'
#' @return A named numeric vector, with the same total duration as `map`.
#'

dash_map <- function (map, dash_length, gap_state)
{
  out <- lapply(seq_along(map), function (k)
  {
    length_k <- map[[k]]
    if (length_k < 3 * dash_length) { return(map[k]) }
    # n_gaps gaps (twice as long as the dashes) between n_gaps + 1 dashes: length = (n_gaps + 1) * dash + n_gaps * 2 * dash
    n_gaps <- max(1, floor((length_k / dash_length - 1) / 3))
    unit <- length_k / (3 * n_gaps + 1)
    n_pieces <- 2 * n_gaps + 1
    stats::setNames(rep(c(unit, 2 * unit), length.out = n_pieces), rep(c(names(map)[k], gap_state), length.out = n_pieces))
  })
  do.call(c, out)
}
