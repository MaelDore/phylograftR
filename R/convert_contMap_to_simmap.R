
##### Convert continuous trait maps (contMap) into discrete stochastic maps (simmap) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - convert_contMap_to_simmap(): converts one contMap, or a list of contMaps, into a simmap, or a multiSimmap, with discrete states
#  - Internal helpers: check_contMaps(), extract_contMap_values(), resolve_contMap_breaks(), discretise_edge()


##### 1/ Exported function #####

#' @title Convert contMap(s) of a continuous trait into simmap(s) with discrete states
#'
#' @description Converts the continuous trait history(ies) mapped on a phylogeny into a stochastic character map of class `"simmap"`
#'   with discrete states obtained by cutting the trait values into bins (e.g., `"low"`, `"medium"`, `"high"`).
#'
#'   The output can be used as a `simmap` in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and
#'   [`phylograftR::graft_missing_clades()`][graft_missing_clades] to constrain the grafting of missing taxa
#'   by the continuous trait (e.g., a taxon known to have a high body size can only be grafted on branches in the `"high"` state).
#'
#'   The function accepts a single contMap, typically built with [`phytools::contMap()`][phytools::contMap],
#'   or a list of contMaps typically obtained from continuous stochastic mapping approach such as the one featured in the `contsimmap` R package.
#'   A single `contMap` gives a single `simmap`. A list of `contMap`s, typically covering the uncertainty of the reconstruction
#'   (e.g., one `contMap` per stochastic history of a continuous trait), gives a `multiSimmap`, with one `simmap` per `contMap`.
#'
#' @param contMap An object of class `"contMap"`, or a list of such objects.
#'   All the `contMap`s of a list must be built on the same phylogeny (same tips, topology, branch lengths and order of the edges).
#'   A list is named after the names of the `contMap`s, with `map_1`, `map_2`, ... if no names are provided.
#' @param breaks Numeric vector. The interior cut points between the bins, in the units of the trait. Values lower than the first
#'   cut point belong to the first state, and so on. They must be provided in increasing order and lie within the range of the trait values (see Details).
#'   Exactly one of `breaks`, `n_bins` and `probs` arguments must be provided.
#' @param n_bins Integer, at least 2. The number of bins of equal width over the range of the trait values observed at the tips. Exactly one of `breaks`, `n_bins` and `probs` arguments must be provided.
#' @param probs Numeric vector with values strictly between 0 and 1. The quantile levels, computed on the trait values observed at the tips, at which to cut the bins
#'   (e.g., `c(1/3, 2/3)` for three bins of equal numbers of tips). Exactly one of `breaks`, `n_bins` and `probs` arguments must be provided.
#' @param states Character vector. The names of the states, from the bin with the lowest values to the bin with the highest values.
#'   Its length must match the number of bins: the number of interior cut points + 1.
#' @param min_segment Numeric between 0 (included) and 1 (excluded). The minimum duration of a segment of a branch, as a proportion of the depth of the tree (root age). Default is `0.01`.
#'   Shorter segments are merged into their longest neighbour on the branch (see Details for why this may matter). Use `0` to keep every crossing of a cut point.
#' @param verbose Logical. Whether to print the progress of the conversion (every 10% of the `contMap`s) when a list of several `contMap`s is provided. Default is `TRUE`.
#'
#' @details
#' ### Principle
#'
#' A `contMap` stores the trait history as a series of short segments along each branch,
#' with ancestral node values typically estimated from a maximum likelihood reconstruction under a Brownian motion (BM) model
#' and linearly interpolated along the branches.
#'
#' During discretization, each segment is assigned to the bin its value falls in, and consecutive segments of the same bin are merged.
#' A value equal to a cut point belongs to the bin above it. The result is a `simmap` (as in [`phytools::make.simmap()`][phytools::make.simmap]) where the states are the bins.
#'
#' ### Cutting the trait values
#'
#' Provide exactly one of the three arguments:
#'   * `breaks`: the cut points are chosen by the user.
#'   * `n_bins`: uses bins of equal width over the range of the tip values.
#'   * `probs`: uses bins defined by quantiles of the tip values.
#'
#' `n_bins` and `probs` are computed on the tip values (averaged across the `contMap`s of a list), which are the observed data.
#' The values of the internal branches are reconstructed and are expected to shrink towards the intermediate values under a BM process
#' (or towards one specific direction if a directional model was used to fit the `contMap`s), because ancestral values are averaged over their descendants.
#' As a consequence, bins defined from the tips will contain a smaller share of the internal branches for the extreme states than for the intermediate states.
#'
#' `breaks` must be interior cut points: they must lie within the range of the trait values, i.e.,
#' between the minimum and the maximum of the values of all the `contMap`s. A cut point beyond this range is an error.
#' A cut point equal to the minimum or the maximum is ignored with a warning, because it would define an empty bin.
#'
#' The number of `states` must be equal to the number of bins (the number of cut points + 1, after the cut points equal to the boundaries were ignored).
#'
#' ### Short segments
#'
#' When using `contMap`s derived from continuous stochastic mapping processes such as the one featured in the `contsimmap` R package,
#' trait values can wander around a cut point: the trait evolutionary history can crosses the threshold several times over short distances,
#' resulting in tiny segments alternating states.
#' Segments shorter than `min_segment` times the depth of the tree are merged, from the shortest to the longest, into their longest neighbour on the same branch.
#' A branch shorter than this minimum cannot change state: it has only one state, the most represented along the branch
#' (for a terminal branch, the state of the tip, so that the states of the tips always match the tip data).
#' On a terminal branch, the last segment (at the tip) is never merged away to ensure tip states matched with observed tip data.
#' Setting `min_segment = 0` keeps all the crossings.
#'
#' This is not an issue for `contMap` featuring linear interpolation of ancestral trait values as in [`phytools::contMap()`][phytools::contMap],
#' since trait values can only vary monotonously along branches (i.e., decrease, increase, or remain stable, but not a combination of two),
#' thus cannot cross each cut point multiple times.
#'
#' ### Cautions for interpretation
#'
#' * The resulting `simmap`s is/are a discretization of the reconstructed trait to constrain the grafting of missing taxa.
#'   It is not a stochastic map of a discrete trait: the number of transitions depends on the resolution of the `contMap` and on `min_segment`
#'   and should not be interpreted as a rate of evolution of discrete states.
#'
#' * Nodes with reconstructed values close to a cut point have uncertain states. Providing a list of `contMap`s
#'   covering the uncertainty (see below) is the way to account for it, as each grafted phylogeny can then use a different `simmap`.
#'
#' * As for the use of `simmap` derived from discrete states, constraining the grafting by a trait that will then be analyzed
#'  (e.g., the association of the trait with diversification) biases the results towards the constraint. Interpret such analyses accordingly.
#'
#' ### Accounting for uncertainty in trait estimates: continuous stochastic mapping
#'
#' A single `contMap` is a point estimate. To cover the uncertainty of the reconstruction, provide a list of `contMap`s, each one built from one possible history of the trait.
#' For example, the `contsimmap` R package samples histories of a continuous trait on a phylogeny, and the function [`deepSTRAPP::convert_contsimmap_to_contMaps()`][deepSTRAPP::convert_contsimmap_to_contMaps]
#' converts them into a list of `contMap`s (see Examples). The output of this function is then a `multiSimmap`, that can be used in the grafting functions like any `multiSimmap`,
#' to account for uncertainty in the ancestral trait reconstructions.
#'
#'
#' @return Depending on the nature of the input:
#' * If `contMap` is a single `contMap`: an object of class `"simmap"`.
#' * If `contMap` is a list of `contMap`s: an object of class `"multiSimmap"` in which each element is a `"simmap"` ordered to match the initial `contMap`s.
#'
#' The interior cut points used to discretize the trait values are stored in the attribute `"breaks"` of the output.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] that use the `simmap` to constrain the grafting.
#'
#' @importFrom ape Ntip
#' @importFrom stats quantile
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Load backbone phylogeny
#' data(toy_backbone)
#'
#' # Simulate a continuous trait on the backbone phylogeny under a Brownian Motion
#' # and interpolate its ML estimates with phytools::contMap()
#' tip_data <- phytools::fastBM(toy_backbone,
#'                           a = 0,  # Ancestral value at the root
#'                           sig2 = 1) # Evolutionary rate as instantaneous variance of the BM process
#' my_contMap <- phytools::contMap(tree = toy_backbone, x = tip_data, plot = FALSE)
#' plot(my_contMap, legend = FALSE)
#'
#' # --- Discretize using quantiles of tip values --- #
#'
#' # Convert into simmap of three states using quantiles of the tip values
#' my_simmap_1 <- convert_contMap_to_simmap(my_contMap, probs = c(1/3, 2/3),
#'                                        states = c("low", "medium", "high"))
#' attr(my_simmap_1, "breaks")
#'
#' # Check the distribution of states across tips
#' tip_states <- phytools::getStates(my_simmap_1, type = "tips")
#' table(tip_states)
#'
#' # Plot the resulting simmap
#' plot(my_simmap_1, colors = c(low = "blue", medium = "grey", high = "red"),
#'      fsize = 0.7, ftype = "i")
#'
#' # --- Discretize using cut points --- #
#'
#' # Use two interior cut points defined by the user to obtain three states
#' my_simmap_2 <- convert_contMap_to_simmap(my_contMap, breaks = c(-1, 1),
#'                                          states = c("low", "medium", "high"))
#'
#' # Plot the resulting simmap
#' plot(my_simmap_2, colors = c(low = "blue", medium = "grey", high = "red"),
#'      fsize = 0.7, ftype = "i")
#'
#' # --- Discretize using bins of equal width --- #
#'
#' # Use bins of equal width to assign states
#' my_simmap_3 <- convert_contMap_to_simmap(my_contMap, n_bins = 3,
#'                                          states = c("low", "medium", "high"))
#'
#' # Plot the resulting simmap
#' plot(my_simmap_3, colors = c(low = "blue", medium = "grey", high = "red"),
#'      fsize = 0.7, ftype = "i")
#'
#'
#' # --- Discretize multiple continuous stochastic maps --- #
#'
#' \dontrun{ # Requires suggested packaged
#' ## The suggested R packages 'deepSTRAPP' and 'contsimmap' are needed for this example to work.
#' if (requireNamespace("deepSTRAPP", quietly = TRUE) && requireNamespace("contsimmap", quietly = TRUE))
#' {
#'  # Use phytools::anc.ML() to extract the evolutionary rate = sigma²
#'  BM_fit <- phytools::anc.ML(tree = rescaled_phylo,
#'                             x = tip_data,
#'                             model = "BM")
#'  # Extract evolutionary rate = sigma²
#'  sigma2 <- BM_fit$sig2
#'
#'  # Simulate 100 continuous stochastic maps for the same trait data and BM model
#'  my_contsimmap <- contsimmap::make.contsimmap(
#'      tree = toy_backbone,
#'      trait.data = tip_data,
#'      Xsig2 = list(sigma2),
#'      nsim = 100)
#'  dim(my_contsimmap) # 58 edges, 1 trait, 100 simulations
#'
#'  # Convert contsimmap output into a list of contMaps
#'  my_contMaps_list <- deepSTRAPP::convert_contsimmap_to_contMaps(
#'      contsimmap = my_contsimmap, verbose = TRUE)
#'
#'  # Convert multiple contMaps at once into corresponding simmaps
#'  # using quantiles of the tip values
#'  my_multiSimmap <- convert_contMap_to_simmap(contMaps_list, probs = c(1/3, 2/3),
#'                                              states = c("low", "medium", "high"))
#'
#'  # Plot the resulting simmaps describing different evolutionary histories simulated
#'  # from the same evolutionary process conditioned to the tip data and phylogeny,
#'  # therefore accounting for uncertainty in ancestral trait estimates.
#'  plot(my_multiSimmap[[1]], colors = c(low = "blue", medium = "grey", high = "red"),
#'      main = "simmap n\u00B01", fsize = 0.7, ftype = "i")
#'  plot(my_multiSimmap[[2]], colors = c(low = "blue", medium = "grey", high = "red"),
#'       main = "simmap n\u00B02", fsize = 0.7, ftype = "i")
#'  plot(my_multiSimmap[[3]], colors = c(low = "blue", medium = "grey", high = "red"),
#'      main = "simmap n\u00B03", fsize = 0.7, ftype = "i")
#' }
#' }
#'

convert_contMap_to_simmap <- function (contMap, # One contMap or a list of contMaps
                                       breaks = NULL, # Interior cut points
                                       n_bins = NULL, # Number of bins of equal width over the range of the tip values
                                       probs = NULL, # Quantile levels of the tip values
                                       states, # Names of the states, from the lowest to the highest values
                                       min_segment = 0.01, # Minimum duration of a segment, as a proportion of the root age
                                       verbose = TRUE) # Print the progress across the contMaps
{
  ## 1/ Check the inputs ####

  # contMap(s)
  if (missing(contMap) || is.null(contMap))
  {
    stop("'contMap' must be an object of class 'contMap', or a list of 'contMap' objects.")
  }
  if (inherits(contMap, "multiContMap"))
  {
    stop("'contMap' cannot be a 'multiContMap' (phytools uses it for the maps of several traits).\n",
         "Provide a plain list of 'contMap' objects instead.")
  }
  is_single <- inherits(contMap, "contMap")
  if (is_single)
  {
    contMaps <- list(contMap)
  } else {
    if (!is.list(contMap) || length(contMap) == 0 || !all(vapply(contMap, inherits, logical(1), what = "contMap")))
    {
      stop("'contMap' must be an object of class 'contMap', or a list of 'contMap' objects.")
    }
    contMaps <- unclass(contMap)
  }
  check_contMaps(contMaps)

  # Which binning argument
  provided <- c(breaks = !is.null(breaks), n_bins = !is.null(n_bins), probs = !is.null(probs))
  if (sum(provided) != 1)
  {
    stop("Exactly one of 'breaks', 'n_bins' and 'probs' must be provided (", sum(provided), " provided).")
  }

  # states
  if (missing(states) || !is.character(states) || length(states) < 2 || anyNA(states) || any(!nzchar(states)) || anyDuplicated(states) > 0)
  {
    stop("'states' must be a character vector of at least two unique non-empty names, from the lowest to the highest trait values.")
  }

  # min_segment
  if (!is.numeric(min_segment) || length(min_segment) != 1 || is.na(min_segment) || min_segment < 0 || min_segment >= 1)
  {
    stop("'min_segment' must be a single number between 0 (included) and 1 (excluded).")
  }

  # verbose
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose))
  {
    stop("'verbose' must be TRUE or FALSE.")
  }

  ## 2/ Extract the values and define the bins ####

  n_maps <- length(contMaps)
  all_values <- lapply(contMaps, extract_contMap_values)

  # Tip values (averaged across maps)
  tip_values <- rowMeans(do.call(cbind, lapply(all_values, function (x) x$tips)))

  # Range of the trait values, pooled across maps
  range_values <- range(unlist(lapply(contMaps, function (x) x$lims)))

  breaks <- resolve_contMap_breaks(breaks = breaks, n_bins = n_bins, probs = probs, tip_values = tip_values, range_values = range_values)

  if (length(states) != length(breaks) + 1)
  {
    stop("'states' must have one name per bin, from the lowest to the highest values: ", length(breaks) + 1, " bins (", length(breaks),
         " interior cut point(s) + 1) but ", length(states), " state(s) provided.")
  }

  ## 3/ Convert each contMap ####

  tree <- contMaps[[1]]$tree
  depth <- max(ape::node.depth.edgelength(tree))
  min_len <- min_segment * depth
  tip_edges <- tree$edge[, 2] <= ape::Ntip(tree)
  n_edges <- nrow(tree$edge)
  n_states <- length(states)

  simmaps <- vector("list", n_maps)
  last_pct <- 0
  for (i in seq_len(n_maps))
  {
    maps_i <- contMaps[[i]]$tree$maps
    values_i <- all_values[[i]]$values
    new_maps <- vector("list", n_edges)
    mapped_edge <- matrix(0, nrow = n_edges, ncol = n_states,
                          dimnames = list(paste(tree$edge[, 1], tree$edge[, 2], sep = ","), states))
    for (e in seq_len(n_edges))
    {
      seg <- discretise_edge(lens = as.numeric(maps_i[[e]]), idx = findInterval(values_i[[e]], breaks) + 1L, min_len = min_len, terminal = tip_edges[e])
      new_map <- seg$lens
      names(new_map) <- states[seg$idx]
      new_maps[[e]] <- new_map
      for (j in seq_along(seg$idx)) { mapped_edge[e, seg$idx[j]] <- mapped_edge[e, seg$idx[j]] + seg$lens[j] }
    }
    simmap_i <- contMaps[[i]]$tree
    simmap_i$maps <- new_maps
    simmap_i$mapped.edge <- mapped_edge
    class(simmap_i) <- c("simmap", "phylo")
    simmaps[[i]] <- simmap_i

    if (verbose && n_maps > 1)
    {
      pct <- floor(10 * i / n_maps) * 10
      if (pct > last_pct)
      {
        message("Converting contMaps: ", pct, "% (", i, "/", n_maps, ")")
        last_pct <- pct
      }
    }
  }

  ## 4/ Format the output ####

  if (is_single)
  {
    output <- simmaps[[1]]
  } else {
    map_names <- names(contMap)
    if (is.null(map_names)) { map_names <- rep("", n_maps) }
    map_names[is.na(map_names) | map_names == ""] <- paste0("map_", seq_len(n_maps))[is.na(map_names) | map_names == ""]
    names(simmaps) <- map_names
    output <- simmaps
    class(output) <- c("multiSimmap", "multiPhylo")
  }
  attr(output, "breaks") <- breaks

  return(output)
}


##### 2/ Internal helpers #####

#' @title Check the contMaps and their consistency
#'
#' @description Check that each `contMap` is well formed, and that all share the same tree (tips, topology, order of the edges, branch lengths).
#'
#' @param contMaps List of objects of class `"contMap"`.
#' @param tol Numeric. Relative tolerance (to the longest branch) for the comparison of the branch lengths. Default is `1e-6`.
#'
#' @noRd
#'
#' @return Invisibly `TRUE`. Stops with an informative message otherwise.
#'

check_contMaps <- function (contMaps, tol = 1e-6)
{
  for (i in seq_along(contMaps))
  {
    who <- if (length(contMaps) > 1) paste0("contMap ", i, " of the list") else "'contMap'"
    cm <- contMaps[[i]]
    if (!inherits(cm$tree, "phylo") || is.null(cm$tree$maps) || is.null(cm$cols) || is.null(cm$lims))
    {
      stop("The ", who, " must have a '$tree' (class 'simmap'), '$cols' and '$lims', as returned by phytools::contMap().")
    }
    if (length(cm$cols) < 2 || length(cm$lims) != 2 || !is.numeric(cm$lims) || anyNA(cm$lims) || !all(is.finite(cm$lims)) || cm$lims[2] <= cm$lims[1])
    {
      stop("The ", who, " has invalid '$cols' or '$lims': the range of its trait values must be a finite interval.")
    }
    idx <- suppressWarnings(as.numeric(unlist(lapply(cm$tree$maps, names), use.names = FALSE)))
    if (anyNA(idx) || any(idx < 0) || any(idx > length(cm$cols) - 1))
    {
      stop("The segments of the ", who, " must be named after the colour indices ('0' to 'length(cols) - 1'), as returned by phytools::contMap().")
    }
    if (i > 1)
    {
      ref <- contMaps[[1]]$tree
      tr <- cm$tree
      if (!identical(tr$tip.label, ref$tip.label))
      {
        stop("The ", who, " has different tips (or in a different order) than the first contMap. All the contMaps must be built on the same phylogeny.")
      }
      if (!identical(unname(tr$edge), unname(ref$edge)))
      {
        stop("The ", who, " has a different topology (or a different order of the edges) than the first contMap. All the contMaps must be built on the same phylogeny.")
      }
      if (any(abs(tr$edge.length - ref$edge.length) > tol * max(ref$edge.length)))
      {
        stop("The ", who, " has different branch lengths than the first contMap. All the contMaps must be built on the same phylogeny.")
      }
    }
  }
  invisible(TRUE)
}


#' @title Extract the trait values of a contMap
#'
#' @description Convert the colour indices of the segments of a `contMap` into trait values (the index `k` of `length(cols)` colours is
#'   the value `lims[1] + k / (length(cols) - 1) * diff(lims)`), and estimate the values at the tips.
#'   The value of a segment is its value at its middle, so the value at the tip is extrapolated from the last two segments of the terminal branch.
#'
#' @param contMap Object of class `"contMap"`.
#'
#' @noRd
#'
#' @return A list with 2 elements:
#'  * `values` List with, for each branch (row of `$tree$edge`), the numeric vector of the values of its segments.
#'  * `tips` Named numeric vector of the estimated values at the tips (in the order of the tips).
#'

extract_contMap_values <- function (contMap)
{
  tree <- contMap$tree
  lims <- contMap$lims
  n_cols <- length(contMap$cols)
  values <- lapply(tree$maps, function (seg) lims[1] + as.numeric(names(seg)) / (n_cols - 1) * (lims[2] - lims[1]))
  tip_edges <- match(seq_along(tree$tip.label), tree$edge[, 2])
  tips <- vapply(tip_edges, function (e) {
    v <- values[[e]]
    s <- as.numeric(tree$maps[[e]])
    k <- length(s)
    if (k < 2) { return(v[k]) }
    v[k] + (v[k] - v[k - 1]) * s[k] / (s[k] + s[k - 1])
  }, numeric(1))
  names(tips) <- tree$tip.label
  list(values = values, tips = pmin(pmax(tips, lims[1]), lims[2]))
}


#' @title Resolve the interior cut points
#'
#' @description Compute the cut points from `breaks`, `n_bins` or `probs` (only one is not `NULL`), and check them.
#'
#' @param breaks,n_bins,probs See [`convert_contMap_to_simmap()`][convert_contMap_to_simmap].
#' @param tip_values Numeric vector. Trait values at the tips.
#' @param range_values Numeric vector of length 2. Range of the trait values (pooled across the contMaps).
#'
#' @noRd
#'
#' @return Sorted numeric vector of the interior cut points.
#'

resolve_contMap_breaks <- function (breaks = NULL, n_bins = NULL, probs = NULL, tip_values, range_values)
{
  if (!is.null(breaks))
  {
    if (!is.numeric(breaks) || length(breaks) == 0 || anyNA(breaks) || !all(is.finite(breaks)))
    {
      stop("'breaks' must be a numeric vector of finite cut points.")
    }
    if (anyDuplicated(breaks) > 0)
    {
      stop("'breaks' must not contain duplicated values.")
    }
    breaks <- sort(breaks)
    tol <- sqrt(.Machine$double.eps) * diff(range_values)
    outside <- breaks < range_values[1] - tol | breaks > range_values[2] + tol
    if (any(outside))
    {
      stop("'breaks' must lie within the range of the trait values [", signif(range_values[1], 6), ", ", signif(range_values[2], 6), "].\n",
           "Cut point(s) outside: ", paste(signif(breaks[outside], 6), collapse = ", "), ".")
    }
    on_boundary <- !outside & (abs(breaks - range_values[1]) <= tol | abs(breaks - range_values[2]) <= tol)
    if (any(on_boundary))
    {
      warning("'breaks' equal to the minimum or the maximum of the trait values (", paste(signif(breaks[on_boundary], 6), collapse = ", "),
              ") were ignored: only interior cut points are expected, as a cut point on a boundary defines an empty bin.\n",
              "This is particularly important for lists of contMaps derived from stochastic histories, as ancestral values can vary beyond the values of the tips.",
              call. = FALSE)
      breaks <- breaks[!on_boundary]
    }
    if (length(breaks) == 0)
    {
      stop("No interior cut point left in 'breaks'. At least one cut point strictly within the range of the trait values is required.")
    }
  } else if (!is.null(n_bins)) {
    if (!is.numeric(n_bins) || length(n_bins) != 1 || is.na(n_bins) || n_bins < 2 || n_bins != round(n_bins))
    {
      stop("'n_bins' must be a single whole number, at least 2.")
    }
    tip_range <- range(tip_values)
    if (tip_range[1] == tip_range[2])
    {
      stop("The trait values at the tips are all equal: bins cannot be defined.")
    }
    breaks <- seq(tip_range[1], tip_range[2], length.out = n_bins + 1)[2:n_bins]
  } else {
    if (!is.numeric(probs) || length(probs) == 0 || anyNA(probs) || any(probs <= 0) || any(probs >= 1) || anyDuplicated(probs) > 0)
    {
      stop("'probs' must be a numeric vector of unique values strictly between 0 and 1.")
    }
    breaks <- unname(stats::quantile(tip_values, probs = sort(probs), names = FALSE))
    if (anyDuplicated(breaks) > 0 || any(breaks <= range_values[1]) || any(breaks >= range_values[2]))
    {
      stop("The quantiles of the tip values in 'probs' do not define distinct interior cut points (many tips share the same value?). Use 'breaks' or fewer 'probs'.")
    }
  }
  breaks
}


#' @title Discretise the segments of one branch
#'
#' @description Merge consecutive segments of the same state, then segments shorter than `min_len` into their longest neighbour (shortest first).
#'   A branch shorter than `min_len` gets a single state (the one with the longest total duration, or the state of the tip for a terminal branch).
#'   On a terminal branch, the segment at the tip is never merged away.
#'
#' @param lens Numeric vector. Durations of the segments, from the root to the tip.
#' @param idx Integer vector. State (bin) index of each segment.
#' @param min_len Numeric. Minimum duration of a segment.
#' @param terminal Logical. Whether the branch is terminal.
#'
#' @noRd
#'
#' @return A list with `lens` (durations) and `idx` (state indices) of the new segments.
#'

discretise_edge <- function (lens, idx, min_len, terminal)
{
  merge_equal <- function (lens, idx)
  {
    if (length(idx) < 2) { return(list(lens = lens, idx = idx)) }
    r <- rle(idx)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    list(lens = vapply(seq_along(ends), function (k) sum(lens[starts[k]:ends[k]]), numeric(1)), idx = r$values)
  }

  seg <- merge_equal(lens, idx)
  if (min_len > 0 && length(seg$idx) > 1)
  {
    total <- sum(lens)
    if (total < min_len)
    {
      if (terminal) {
        state <- idx[length(idx)]
      } else {
        by_state <- tapply(lens, idx, sum)
        state <- as.integer(names(by_state)[which.max(by_state)])
      }
      return(list(lens = total, idx = state))
    }
    repeat
    {
      k <- length(seg$idx)
      if (k < 2) { break }
      mergeable <- seg$lens < min_len
      if (terminal) { mergeable[k] <- FALSE }
      if (!any(mergeable)) { break }
      i <- which(mergeable)[which.min(seg$lens[mergeable])]
      neighbours <- c(if (i > 1) i - 1, if (i < k) i + 1)
      j <- neighbours[which.max(seg$lens[neighbours])]
      seg$lens[j] <- seg$lens[j] + seg$lens[i]
      seg <- merge_equal(seg$lens[-i], seg$idx[-i])
    }
  }
  seg
}
