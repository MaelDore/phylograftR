
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
#'   **A single `contMap` gives a single `simmap`.** A list of `contMap`s, typically covering the uncertainty of the reconstruction
#'   (e.g., one `contMap` per stochastic history of a continuous trait), gives a **`multiSimmap`**, with **one `simmap` per `contMap`**.
#'
#' @param contMap An object of class `"contMap"`, or a list of such objects.
#'   **All the `contMap`s of a list must be built on the same phylogeny** (same tips, topology, branch lengths and order of the edges).
#'   A list is named after the names of the `contMap`s, with `map_1`, `map_2`, ... if no names are provided.
#' @param breaks Numeric vector. The interior cut points between the bins, in the units of the trait. Values lower than the first
#'   cut point belong to the first state, and so on. They must be provided in increasing order, be different values, and lie within the range of the trait values (see Details).
#'   If not, a warning is sent and the cut points are sorted and/or the duplicated values removed.
#'   Exactly one of `breaks`, `n_bins` and `probs` arguments must be provided.
#' @param n_bins Integer, at least 2. The number of bins of equal width over the range of the trait values at the tips.
#'   Exactly one of `breaks`, `n_bins` and `probs` arguments must be provided.
#' @param probs Numeric vector with values strictly between 0 and 1. The quantile levels, computed on the trait values observed at the tips, at which to cut the bins.
#'   They must be provided in increasing order and be different values, otherwise a warning is sent and they are sorted and/or the duplicated values removed
#'   (e.g., `c(1/3, 2/3)` for three bins of equal numbers of tips). Exactly one of `breaks`, `n_bins` and `probs` arguments must be provided.
#' @param states Character vector. The names of the states, from the bin with the lowest values to the bin with the highest values.
#'   Its length must match the number of bins: **the number of interior cut points + 1**.
#' @param tip_values Numeric vector (optional). The **observed trait values at the tips**, named after the tips of the phylogeny (or provided in the order of the tips if unnamed).
#'   If provided, they are used instead of the values read from the `contMap`(s) to **compute the cut points** (`n_bins` and `probs`) and to **set the state of each tip**
#'   (the last segment of each terminal branch), so that the tips are in the state of their observed value. The same values are used for all the `contMap`s of a list.
#'   They must lie within the range of the trait values of the `contMap`(s). If `NULL`, the tip values are read from the `contMap`(s), where they are only approximated (see Details).
#' @param min_segment_proportion Numeric between 0 (included) and 1 (excluded). The minimum duration of a segment of a branch, as a proportion of the depth of the tree (root age). Default is `0.01`.
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
#' **A value equal to a cut point belongs to the bin above it.** The result is a `simmap` (as in [`phytools::make.simmap()`][phytools::make.simmap]) where the states are the bins.
#'
#' ### Cutting the trait values
#'
#' **Provide exactly one of the three arguments:**
#'   * `breaks`: the cut points are chosen by the user.
#'   * `n_bins`: uses bins of equal width over the range of the tip values.
#'   * `probs`: uses bins defined by quantiles of the tip values.
#'
#' `n_bins` and `probs` are computed on the tip values (averaged across the `contMap`s of a list).
#' By default, they are read from the `contMap`s, which store the trait values on a finite scale of colors (from 0 to 1000 in `phytools::contMap()`), so they approximate the observed data.
#' **Provide the observed values with `tip_values` to use them instead** (see below).
#'
#' ### Tip values
#'
#' A `contMap` stores the trait history as a series of short segments along each branch, with the value of each segment interpolated from the estimated values at the nodes.
#' The value associated to a segment is a single value for the whole segment, on a finite scale of colors. As a result, the value of the last segment of a terminal branch
#' is **nearly, but not identical to the observed value at the tip**. With the default settings, the function extrapolates the tip value from the last segments, which is only an approximate.
#' The difference is usually small compared to the range of the trait, but it can move the tip of a taxon with a value close to a cut point into the neighboring state,
#' or shift the cut points obtained from `n_bins` or `probs`.
#'
#' To avoid this, provide the **observed tip values** with `tip_values`: they are then used to compute the cut points, and **each tip is set in the state of its observed value**,
#' so that the tips of the `simmap`s are in the same states as the observed data converted with [`phylograftR::convert_values_to_states()`][convert_values_to_states].
#' This is recommended whenever the tip data used to build the `contMap`(s) are available.
#' When the observed state of a tip differs from the state of the interpolated segments just before it, the terminal branch ends with a short segment in the state of the observed value
#' (the segment at the tip is never merged away, see `min_segment_proportion`).
#'
#' The values of the internal branches are reconstructed and are expected to shrink towards the intermediate values under a BM process
#' (or towards one specific direction if a directional model was used to fit the `contMap`s), because ancestral values are averaged over their descendants.
#' As a consequence, bins defined from the tips will contain a smaller share of the internal branches for the extreme states than for the intermediate states.
#'
#' `breaks` must be interior cut points: they must lie within the range of the trait values, i.e.,
#' between the minimum and the maximum of the values of all the `contMap`s. **A cut point beyond this range is an error.**
#' A cut point equal to the minimum or the maximum is ignored with a warning, because it would define an empty bin.
#'
#' The cut points must be in increasing order and different values. This is checked for the three arguments, including the cut points derived from `n_bins` and `probs`
#' (e.g., many tips sharing the same value can give identical quantiles). Unsorted cut points are sorted, and duplicated values are removed, with an informative warning.
#'
#' **The number of `states` must be equal to the number of bins** (the number of cut points + 1, after the cut points equal to the boundaries were ignored).
#' Removing cut points reduces the number of bins, so `states` must match the number of bins left.
#'
#' ### Short segments
#'
#' When using `contMap`s derived from continuous stochastic mapping processes such as the one featured in the `contsimmap` R package,
#' trait values can wander around a cut point: the trait evolutionary history can crosses the threshold several times over short distances,
#' resulting in tiny segments alternating states.
#' **Segments shorter than `min_segment_proportion` times the depth of the tree are merged**, from the shortest to the longest, into their longest neighbour on the same branch.
#' A branch shorter than the minimum length gets a single state: the one with the longest total duration, or the state of the tip for a terminal branch.
#' On a terminal branch, the segment at the tip is never merged away to ensure tip states matched with observed tip data.
#' Setting `min_segment_proportion = 0` keeps all the crossings.
#'
#' This is not an issue for `contMap` featuring linear interpolation of ancestral trait values as in [`phytools::contMap()`][phytools::contMap],
#' since trait values can only vary monotonously along branches (i.e., decrease, increase, or remain stable, but not a combination of two),
#' thus cannot cross each cut point multiple times.
#'
#' ### Cautions for interpretation
#'
#' * The resulting `simmap`s is/are a discretization of the reconstructed trait to constrain the grafting of missing taxa.
#'   **It is not a stochastic map of a discrete trait**: the number of transitions depends on the resolution of the `contMap` and on `min_segment_proportion`
#'   and should not be interpreted as a rate of evolution of discrete states.
#'
#' * Nodes with reconstructed values close to a cut point have uncertain states. Providing a list of `contMap`s
#'   covering the uncertainty (see below) is the way to account for it, as each grafted phylogeny can then use a different `simmap`.
#'
#' * The states depend on the cut points. A taxon whose allowed states are defined by a threshold can have different eligible branches
#'  if the threshold is moved. Ideally, you would want to compare several sets of cut points before relying on the results (see `vignette("continuous_data")`).
#'
#' * As for the use of `simmap` derived from discrete states, constraining the grafting by a trait that will then be analyzed
#'  (e.g., the association of the trait with diversification) biases the results towards the constraint. Interpret such analyses accordingly.
#'
#' ### Accounting for uncertainty in trait estimates: continuous stochastic mapping
#'
#' A single `contMap` is a point estimate. To cover the uncertainty of the reconstruction, provide a list of `contMap`s, each one built from one possible history of the trait.
#' For example, the `contsimmap` R package samples histories of a continuous trait on a phylogeny, and the function [`deepSTRAPP::convert_contsimmap_to_contMaps()`][deepSTRAPP::convert_contsimmap_to_contMaps]
#' converts them into a list of `contMap`s (see Examples). The output of this function is then a **`multiSimmap`**, that can be used in the grafting functions like any `multiSimmap`,
#' to account for uncertainty in the ancestral trait reconstructions.
#'
#' The function `deepSTRAPP::prepare_trait_data()` is the **recommended way** to obtain such a list of `contMap`s: with `run_stochastic_maps = TRUE`, it fits the evolutionary model(s) (selected by AICc),
#' runs the continuous stochastic mapping with `contsimmap::make.contsimmap()` and converts the maps into `contMap`s in a single call, with additional options to fine-tune the process (see Examples).
#' The example below uses `contsimmap` directly to credit the authors of the continuous stochastic mapping approach and to show the steps usde internally by `deepSTRAPP::prepare_trait_data()`.
#'
#'
#' @return Depending on the nature of the input:
#' * If `contMap` is a single `contMap`: an object of class `"simmap"`.
#' * If `contMap` is a list of `contMap`s: an object of class `"multiSimmap"` in which each element is a `"simmap"` ordered to match the initial `contMap`s.
#'
#' The interior cut points used to discretize the trait values are stored in the attribute **`"breaks"`** of the output,
#' and the names of the states (from the lowest to the highest values) in its attribute **`"states"`**.
#' They allow to convert other values, e.g., the trait values of the missing taxa, into the same states
#' with [`phylograftR::convert_values_to_states()`][convert_values_to_states].
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::convert_values_to_states()`][convert_values_to_states] to convert continuous values of the missing taxa into the states of the `simmap`.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa based on target clades and excluded subclades, using the states to constrain the grafting.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades, using the states to constrain the grafting.
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
#' # --- Use the observed tip values --- #
#'
#' # The contMap only approximates the values at the tips (interpolation within segments).
#' # Provide the observed tip values to compute the cut points and to set the state of each tip from the data
#' my_simmap_1b <- convert_contMap_to_simmap(my_contMap, probs = c(1/3, 2/3),
#'                                           states = c("low", "medium", "high"),
#'                                           tip_values = tip_data)
#' # Cut points computed from the observed tip values
#' attr(my_simmap_1b, "breaks")
#' quantile(tip_data, probs = c(1/3, 2/3))
#'
#' # Tips are in the state of their observed values
#' tip_states_b <- phytools::getStates(my_simmap_1b, type = "tips")
#' all(tip_states_b == unlist(convert_values_to_states(tip_data, simmap = my_simmap_1b))[names(tip_states_b)])
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
#' \dontrun{ # Requires one package not available on CRAN
#' ## The R packages 'deepSTRAPP' (on CRAN) and 'contsimmap' (not on CRAN)
#' ## are needed for this example to work.
#' if (requireNamespace("deepSTRAPP", quietly = TRUE)
#'     && requireNamespace("contsimmap", quietly = TRUE))
#' {
#'  # Use phytools::anc.ML() to extract the evolutionary rate = sigma²
#'  BM_fit <- phytools::anc.ML(tree = toy_backbone,
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
#'  my_multiSimmap <- convert_contMap_to_simmap(my_contMaps_list, probs = c(1/3, 2/3),
#'                                              states = c("low", "medium", "high"),
#'                                              tip_values = tip_data) # Use the observed tip values
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
#'
#'  # --- Alternative: use deepSTRAPP::prepare_trait_data() (recommended) --- #
#'
#'  # The steps above (fit of the BM model, simulation of the continuous stochastic maps
#'  # with contsimmap::make.contsimmap(), and conversion into a list of contMaps)
#'  # can be run at once with deepSTRAPP::prepare_trait_data().
#'
#'  # The function offers additional options to fine-tune the process, such as the comparison of several
#'  # evolutionary models (selected by AICc), the resolution of the contMaps, and the number of simulations.
#'  # It can also output directly both the interpolated point estimate as in phytools::contMap(),
#'  # and the set of stochastic continuous maps simulated with contsimmap::make.contsimmap().
#'
#'  # Not to run: it replaces the steps above and returns the same type of output (a list of contMaps).
#'  # my_trait_data <- deepSTRAPP::prepare_trait_data(
#'  #     tip_data = tip_data,
#'  #     trait_data_type = "continuous",
#'  #     phylo = toy_backbone,
#'  #     seed = 1234, # Set seed for reproducibility
#'  #     evolutionary_models = "BM", # Evolutionary model(s) to fit
#'  #     run_stochastic_maps = TRUE, # Run continuous stochastic mapping to account for uncertainty
#'  #     nb_simulations = 100, # Number of stochastic maps
#'  #     plot_map = FALSE)
#'  #
#'  # my_contMap <- my_trait_data$contMap # Point estimate: a single contMap
#'  # my_contMaps_list <- my_trait_data$contMaps # List of contMaps from the stochastic maps
#'  #
#'  # my_multiSimmap <- convert_contMap_to_simmap(my_contMaps_list, probs = c(1/3, 2/3),
#'  #                                             states = c("low", "medium", "high"),
#'  #                                             tip_values = tip_data)
#' }
#' }
#'

convert_contMap_to_simmap <- function (contMap, # One contMap or a list of contMaps
                                       breaks = NULL, # Interior cut points
                                       n_bins = NULL, # Number of bins of equal width over the range of the tip values
                                       probs = NULL, # Quantile levels of the tip values
                                       states, # Names of the states, from the lowest to the highest values
                                       tip_values = NULL, # Observed trait values at the tips (optional)
                                       min_segment_proportion = 0.01, # Minimum duration of a segment, as a proportion of the root age
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

  # tip_values: observed values at the tips, reordered to match the tips of the phylogeny
  tip_labels <- contMaps[[1]]$tree$tip.label
  if (!is.null(tip_values))
  {
    if (!is.numeric(tip_values) || anyNA(tip_values) || !all(is.finite(tip_values)))
    {
      stop("'tip_values' must be a numeric vector of finite values (no NA).")
    }
    if (length(tip_values) != length(tip_labels))
    {
      stop("'tip_values' must have one value per tip of the phylogeny (", length(tip_labels), "), but ", length(tip_values), " value(s) were provided.")
    }
    if (is.null(names(tip_values)))
    {
      # Unnamed: assumed to be in the order of the tips of the phylogeny
      names(tip_values) <- tip_labels
    } else {
      if (anyDuplicated(names(tip_values)) > 0 || !setequal(names(tip_values), tip_labels))
      {
        stop("The names of 'tip_values' must be the tip labels of the contMap(s), each used once.")
      }
      tip_values <- tip_values[tip_labels]
    }
    # The observed values must lie within the range of the values of the contMap(s) (small tolerance for rounding)
    lims_all <- range(unlist(lapply(contMaps, function (x) x$lims)))
    tol_lims <- 1e-6 * diff(lims_all)
    if (any(tip_values < lims_all[1] - tol_lims | tip_values > lims_all[2] + tol_lims))
    {
      stop("Some 'tip_values' are outside the range of the trait values of the contMap(s) (", format(lims_all[1]), " to ", format(lims_all[2]), ").\n",
           "'tip_values' must be the observed data used to build the contMap(s).")
    }
  }

  # min_segment_proportion
  if (!is.numeric(min_segment_proportion) || length(min_segment_proportion) != 1 || is.na(min_segment_proportion) || min_segment_proportion < 0 || min_segment_proportion >= 1)
  {
    stop("'min_segment_proportion' must be a single number between 0 (included) and 1 (excluded).")
  }

  # verbose
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose))
  {
    stop("'verbose' must be TRUE or FALSE.")
  }

  ## 2/ Extract the values and define the bins ####

  # Get number of maps
  n_maps <- length(contMaps)
  # Extract all values along all branches of all maps (requires back-conversion into the trait scale)
  all_values <- lapply(contMaps, extract_contMap_values)
  # Tip values used to compute the cut points: the observed values if provided, or else the values read from the contMaps (averaged across maps), which are only approximate
  tip_values_cutpoints <- if (!is.null(tip_values)) tip_values else rowMeans(do.call(cbind, lapply(all_values, function (x) x$tips)))
  # Extract the max range of the trait values across all maps based on '$lims'
  range_values <- range(unlist(lapply(contMaps, function (x) x$lims)))

  # Interior cut points: from 'breaks', or computed from the tip values ('n_bins', 'probs').
  # Always sorted, distinct, and strictly within 'range_values' (warnings or errors otherwise)
  breaks <- resolve_contMap_breaks(breaks = breaks, n_bins = n_bins, probs = probs, tip_values = tip_values_cutpoints, range_values = range_values)

  # Check that one state per bin is provided: k interior cut points define k + 1 bins.
  # This is checked after the cut points are resolved, as some may have been removed (duplicates, boundaries)
  if (length(states) != length(breaks) + 1)
  {
    stop("'states' must have one name per bin, from the lowest to the highest values: ", length(breaks) + 1, " bins (", length(breaks),
         " interior cut point(s) + 1) but ", length(states), " state(s) provided.")
  }

  ## 3/ Convert each contMap into a simmap ####

  ## Compute tree-level stuff
  # All the contMaps share the same tree (checked in check_contMaps()), so the quantities that depend on the tree are computed once.

  tree <- contMaps[[1]]$tree
  # Depth of the tree = root age for an ultrametric tree = largest distance from the root to a tip
  depth <- max(ape::node.depth.edgelength(tree))
  # Minimum accepted duration of a segment, in the units of the branch lengths (e.g., Myr)
  min_len <- min_segment_proportion * depth
  # Identify which edges are terminal (end on a tip): their last segment is kept to preserve the tip state
  tip_edges <- tree$edge[, 2] <= ape::Ntip(tree)
  n_edges <- nrow(tree$edge)
  n_states <- length(states)
  # State (bin index) of the observed value of each tip, if provided. In the order of the tips, i.e., indexed by the tipward node of the terminal edges
  tip_bins <- if (!is.null(tip_values)) findInterval(tip_values, breaks) + 1L else NULL

  # Initiate the list of simmaps
  simmaps <- vector("list", n_maps)
  # Initiate percentage of progress
  last_pct <- 0

  ## Loop per map
  for (i in seq_len(n_maps))
  {
    # Original segments of the contMap: for each edge, a named vector with the duration of each segment from the root to the tip.
    # The names are the colour indices, not the states
    maps_i <- contMaps[[i]]$tree$maps
    # Trait values of the same segments: list (one element per edge) of numeric vectors of the same lengths as 'maps_i'
    values_i <- all_values[[i]]$values

    # Initiate elements of the new simmap: '$maps' (segments per edge) and '$mapped.edge' (total duration per edge and state, as in phytools::make.simmap())
    new_maps <- vector("list", n_edges)
    mapped_edge <- matrix(0, nrow = n_edges, ncol = n_states,
                          dimnames = list(paste(tree$edge[, 1], tree$edge[, 2], sep = ","), states))

    # Loop per edge
    for (e in seq_len(n_edges))
    {
      # 1. Assign each segment to a bin: findInterval() counts the cut points <= value, so values below the first cut point give 0,
      #    between the first and second give 1, ..., above the last give k. The + 1L turns it into a bin index from 1 to k + 1
      #    (a value equal to a cut point goes to the bin above it).
      idx_e <- findInterval(values_i[[e]], breaks) + 1L
      # 2. With observed tip values, the last segment of a terminal edge (at the tip) is set to the state of the observed value of the tip,
      #    as the value of this segment only approximates it (this segment is kept as it is by discretise_edge())
      if (!is.null(tip_bins) && tip_edges[e]) { idx_e[length(idx_e)] <- tip_bins[tree$edge[e, 2]] }
      # 3. discretise_edge() merges the neighbour segments of the same bin and removes the segments shorter than 'min_len'
      seg <- discretise_edge(lens = as.numeric(maps_i[[e]]), idx = idx_e, min_len = min_len, terminal = tip_edges[e])

      # New map of the edge: named vector of durations, named after the states (as in a simmap)
      new_map <- seg$lens
      names(new_map) <- states[seg$idx]
      new_maps[[e]] <- new_map

      # Sum the durations per state in '$mapped.edge' (a state can occur in several segments of the same edge)
      for (j in seq_along(seg$idx)) { mapped_edge[e, seg$idx[j]] <- mapped_edge[e, seg$idx[j]] + seg$lens[j] }
    }

    # Build the simmap: the tree of the contMap with the new maps elements.
    # The tree structure (edges, branch lengths, tips) is unchanged
    simmap_i <- contMaps[[i]]$tree
    simmap_i$maps <- new_maps
    simmap_i$mapped.edge <- mapped_edge
    class(simmap_i) <- c("simmap", "phylo")
    simmaps[[i]] <- simmap_i

    # Progress: message each time a new multiple of 10% of the maps is reached (only for lists of several contMaps)
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
    # One contMap in: one simmap out
    output <- simmaps[[1]]
  } else {
    # A list in: a multiSimmap out, named after the list. Missing or empty names are replaced by map_1, map_2, ... (position in the list)
    map_names <- names(contMap)
    if (is.null(map_names)) { map_names <- rep("", n_maps) }
    map_names[is.na(map_names) | map_names == ""] <- paste0("map_", seq_len(n_maps))[is.na(map_names) | map_names == ""]
    names(simmaps) <- map_names
    output <- simmaps
    class(output) <- c("multiSimmap", "multiPhylo")
  }
  # In both cases, keep the cut points and the names of the states, so that other values (e.g., of the missing taxa) can be binned the same way (see convert_values_to_states())
  attr(output, "breaks") <- breaks
  attr(output, "states") <- states

  return(output)
}


##### 2/ Internal helpers #####

#' @title Check the contMaps and their consistency
#'
#' @description Check that each `contMap` has a valid structure, and that all share the same tree (tips, topology, order of the edges, branch lengths).
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
#' @description Convert the colour indices of the segments of a `contMap` into trait values (the index `k` of `length(cols)` colors is
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
  lims <- contMap$lims # Range of the trait values mapped onto the colors
  n_cols <- length(contMap$cols) # Number of colors = number of trait value levels (res + 1 in phytools::contMap())

  # Color index k (0 to n_cols - 1, stored as the name of the segment) -> trait value: linear scale from lims[1] (k = 0) to lims[2] (k = n_cols - 1).
  # The resolution of the values is therefore diff(lims) / (n_cols - 1)
  values <- lapply(tree$maps, function (seg) lims[1] + as.numeric(names(seg)) / (n_cols - 1) * (lims[2] - lims[1]))

  # Tip values. The value of a segment is the value at its middle, so the last segment of a terminal edge gives the value
  # half a segment before the tip. The tip value is extrapolated linearly from the last two segments:
  # slope between the middles of the last two segments = (v[k] - v[k - 1]) / ((s[k] + s[k - 1]) / 2),
  # applied over the half-length of the last segment (s[k] / 2) = (v[k] - v[k - 1]) * s[k] / (s[k] + s[k - 1]).
  # With a single segment, its value is used as it is.
  tip_edges <- match(seq_along(tree$tip.label), tree$edge[, 2]) # Terminal edge of each tip, in the order of the tips
  tips <- vapply(tip_edges, function (e) {
    v <- values[[e]] # Values of the segments of the edge
    s <- as.numeric(tree$maps[[e]]) # Their durations
    k <- length(s) # Index of the last segment (at the tip)
    if (k < 2) { return(v[k]) }
    v[k] + (v[k] - v[k - 1]) * s[k] / (s[k] + s[k - 1])
  }, numeric(1))
  names(tips) <- tree$tip.label

  # Values are clipped to the range of the contMap (the extrapolation can overshoot it)
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
  # Source of the cut points, used in the messages
  if (!is.null(breaks))
  {
    if (!is.numeric(breaks) || length(breaks) == 0 || anyNA(breaks) || !all(is.finite(breaks)))
    {
      stop("'breaks' must be a numeric vector of finite cut points.")
    }
    source <- "'breaks'"
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
    source <- "the cut points derived from 'n_bins'"
  } else {
    if (!is.numeric(probs) || length(probs) == 0 || anyNA(probs) || any(probs <= 0) || any(probs >= 1))
    {
      stop("'probs' must be a numeric vector of values strictly between 0 and 1.")
    }
    if (is.unsorted(probs, strictly = FALSE))
    {
      warning("'probs' must be provided in increasing order: they were sorted. The 'states' are assigned from the lowest to the highest values, i.e., following the sorted 'probs'.", call. = FALSE)
      probs <- sort(probs)
    }
    if (anyDuplicated(probs) > 0)
    {
      warning("'probs' must be different values: the duplicated values were removed.", call. = FALSE)
      probs <- unique(probs)
    }
    breaks <- unname(stats::quantile(tip_values, probs = probs, names = FALSE))
    source <- "the cut points derived from 'probs'"
  }

  # Increasing order and different values (for the three arguments)
  if (is.unsorted(breaks, strictly = FALSE))
  {
    warning("The cut points must be provided in increasing order: ", source, " were sorted. The 'states' are assigned from the lowest to the highest values, i.e., following the sorted cut points.", call. = FALSE)
    breaks <- sort(breaks)
  }
  if (anyDuplicated(breaks) > 0)
  {
    warning("The cut points must be different values: duplicated values among ", source, " (", paste(signif(breaks[duplicated(breaks)], 6), collapse = ", "),
            ") were removed", if (source == "the cut points derived from 'probs'") " (many tips share the same value?)" else "",
            ". Fewer bins than expected are defined: check the number of 'states'.", call. = FALSE)
    breaks <- unique(breaks)
  }

  # Interior cut points
  tol <- sqrt(.Machine$double.eps) * diff(range_values)
  outside <- breaks < range_values[1] - tol | breaks > range_values[2] + tol
  if (any(outside))
  {
    stop("The cut points must lie within the range of the trait values [", signif(range_values[1], 6), ", ", signif(range_values[2], 6), "].\n",
         "Cut point(s) outside: ", paste(signif(breaks[outside], 6), collapse = ", "), ".")
  }
  on_boundary <- abs(breaks - range_values[1]) <= tol | abs(breaks - range_values[2]) <= tol
  if (any(on_boundary))
  {
    warning("Cut points equal to the minimum or the maximum of the trait values (", paste(signif(breaks[on_boundary], 6), collapse = ", "),
            ") were ignored: only interior cut points are expected, as a cut point on a boundary defines an empty bin.\n",
            "This is particularly important for lists of contMaps derived from stochastic histories, as ancestral values can vary beyond the values of the tips.",
            call. = FALSE)
    breaks <- breaks[!on_boundary]
  }
  if (length(breaks) == 0)
  {
    stop("No interior cut point left. At least one cut point strictly within the range of the trait values is required.")
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
  # Merge consecutive segments with the same state index into a single segment (durations summed).
  # rle() gives the runs of identical consecutive states: e.g., idx = c(1, 1, 2, 2, 2, 1) -> values c(1, 2, 1), lengths c(2, 3, 1)
  merge_equal <- function (lens, idx)
  {
    if (length(idx) < 2) { return(list(lens = lens, idx = idx)) }
    r <- rle(idx)
    ends <- cumsum(r$lengths) # Position of the last segment of each run
    starts <- ends - r$lengths + 1 # Position of the first segment of each run
    list(lens = vapply(seq_along(ends), function (k) sum(lens[starts[k]:ends[k]]), numeric(1)), idx = r$values)
  }

  # Step 1: Merge the neighbours with the same state. This is all that is done if min_len = 0 (every crossing of a cut point is kept)
  seg <- merge_equal(lens, idx)

  # Step 2: Remove the segments that are too short (only if there is a minimum, and several segments left)
  if (min_len > 0 && length(seg$idx) > 1)
  {
    total <- sum(lens) # Length of the edge

    # Case A: the whole edge is shorter than the minimum: it cannot change state, so it gets a single segment
    if (total < min_len)
    {
      if (terminal) {
        # Terminal edge: the state of the tip (last original segment), so that the tip keeps the state of its data
        state <- idx[length(idx)]
      } else {
        # Internal edge: the state with the longest total duration along the edge (the lowest state index in case of a tie)
        by_state <- tapply(lens, idx, sum)
        state <- as.integer(names(by_state)[which.max(by_state)])
      }
      return(list(lens = total, idx = state))
    }

    # Case B: the edge is long enough to have several segments, but some of them may be too short.
    # Merge them one at a time, the shortest first, until no short segment is left to merge.
    # Merging only one segment per iteration matters: the neighbours it is merged into become longer, and may not need to be merged themselves.
    repeat
    {
      k <- length(seg$idx) # Current number of segments
      if (k < 2) { break } # A single segment cannot be merged with anything

      # Segments shorter than the minimum, except the last segment (at the tip) of a terminal edge, which is protected:
      # other segments can be merged into it, but it is never absorbed, so the tip keeps its state
      mergeable <- seg$lens < min_len
      if (terminal) { mergeable[k] <- FALSE }
      if (!any(mergeable)) { break } # All remaining segments are long enough (or protected): done

      # The shortest mergeable segment ...
      i <- which(mergeable)[which.min(seg$lens[mergeable])]
      # ... is absorbed by its longest neighbour on the edge (the one rootward in case of a tie, as it comes first in 'neighbours').
      # A segment at the end of the edge has a single neighbour
      neighbours <- c(if (i > 1) i - 1, if (i < k) i + 1)
      j <- neighbours[which.max(seg$lens[neighbours])]
      seg$lens[j] <- seg$lens[j] + seg$lens[i] # The neighbour is extended over the absorbed segment (the state of the neighbour is kept)
      # Remove the absorbed segment, and merge the neighbours that now have the same state (e.g., A-B-A: if B is absorbed, A and A are merged)
      seg <- merge_equal(seg$lens[-i], seg$idx[-i])
    }
  }
  seg
}
