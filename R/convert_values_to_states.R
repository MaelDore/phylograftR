
##### Function to convert continuous trait values into the states of a discretized trait #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### 1/ Main function ####

#' @title Convert continuous trait values into the states of a discretized trait
#'
#' @description Assigns each trait value to a state, using the same cut points and states as a `simmap` obtained with
#'   [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap].
#'   This is typically used to convert the continuous trait values of the missing taxa as states of the `simmap` used to constrain the grafting,
#'   e.g., to fill the `allowed_states` of the taxa.
#'
#' @param values Numeric vector. The trait values to convert, in the same units as the trait used to build the `simmap`. If named, the names are kept in the output (e.g., the names of the missing taxa).
#'   `NA` values give `NA` states.
#' @param simmap A `"simmap"` or `"multiSimmap"` produced by [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap],
#'   which carries the cut points and the names of the states in its attributes `"breaks"` and `"states"`.
#'   Alternatively, leave it as `NULL` and provide `breaks` and `states` (e.g., for a `simmap` built from tip values binned by yourself).
#' @param breaks Numeric vector. Only used if `simmap = NULL`. The interior cut points between the states, in increasing order and different from each other.
#' @param states Character vector. Only used if `simmap = NULL`. The names of the states, from the lowest to the highest values.
#'   Must have one more element than `breaks`, because N interior cuts (`breaks`) produce N+1 groups.
#'
#' @details Each state covers the values from its lower cut point (included) up to its upper cut point (excluded): a value equal to a cut point belongs to the state above it.
#'   The first and last states are open-ended, so values below the first cut point or above the last one are assigned to the lowest and highest states.
#'
#'   Each value is assigned to a unique state. If you wish to include uncertainty by assigning multiple states to values close to cut points (as possible for `allowed_states`),
#'   you must edit the output after the conversion.
#'
#' @return A character vector with the state of each value, named as in `values`.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap] to obtain the `simmap` and its cut points from a continuous trait map (`contMap`).
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] that use the states to constrain the grafting.
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Load backbone phylogeny
#' data(toy_backbone)
#'
#' # Simulate a continuous trait on the backbone and interpolate it with phytools::contMap()
#' tip_data <- phytools::fastBM(toy_backbone, a = 0, sig2 = 1)
#' my_contMap <- phytools::contMap(tree = toy_backbone, x = tip_data, plot = FALSE)
#'
#' # Discretize the trait into three states using quantiles
#' my_simmap <- convert_contMap_to_simmap(my_contMap, probs = c(1/3, 2/3),
#'                                        states = c("low", "medium", "high"))
#' # Inspect the trait values used as breaks
#' attr(my_simmap, "breaks")
#'
#' # Define continuous trait values of three missing taxa
#' missing_values <- c(Taxon_1 = -2, Taxon_2 = 0.1, Taxon_3 = 3)
#' # Convert continuous values into the same states used in the converted simmap
#' convert_values_to_states(missing_values, simmap = my_simmap)
#'
#' # Same conversion with the cut points and states given directly
#' convert_values_to_states(missing_values, breaks = attr(my_simmap, "breaks"),
#'                          states = c("low", "medium", "high"))
#'

convert_values_to_states <- function (values, # Numeric trait values to convert, named after the taxa if wanted
                                      simmap = NULL, # simmap or multiSimmap from convert_contMap_to_simmap(), with attributes "breaks" and "states"
                                      breaks = NULL, # Interior cut points, if simmap = NULL
                                      states = NULL) # Names of the states from the lowest to the highest values, if simmap = NULL
{
  ### ---- Checks ---- ###

  if (!is.numeric(values)) { stop("'values' must be a numeric vector of trait values.") }
  if (any(is.infinite(values))) { stop("'values' must not contain infinite values.") }

  if (!is.null(simmap))
  {
    # The cut points and the states come from the simmap
    if (!is.null(breaks) || !is.null(states)) { stop("Provide either 'simmap', or 'breaks' and 'states', not both.") }
    if (!inherits(simmap, c("simmap", "multiSimmap"))) { stop("'simmap' must be a 'simmap' or a 'multiSimmap'.") }
    breaks <- attr(simmap, "breaks")
    states <- attr(simmap, "states")
    if (is.null(breaks) || is.null(states))
    {
      stop("'simmap' has no attributes \"breaks\" and \"states\": it was not produced by convert_contMap_to_simmap() (or by a version that stores them). ",
           "Provide 'breaks' and 'states' instead.")
    }
  } else {
    if (is.null(breaks) || is.null(states)) { stop("Provide either 'simmap', or both 'breaks' and 'states'.") }
  }

  if (!is.numeric(breaks) || length(breaks) < 1 || anyNA(breaks) || any(!is.finite(breaks))) { stop("'breaks' must be a numeric vector of finite cut points.") }
  if (is.unsorted(breaks, strictly = TRUE)) { stop("'breaks' must be in increasing order and have different values.") }
  if (!is.character(states) || anyNA(states) || any(states == "") || anyDuplicated(states) > 0) { stop("'states' must be a character vector of unique, non-empty names.") }
  if (length(states) != length(breaks) + 1)
  {
    stop("'states' must have one more element than 'breaks': ", length(breaks), " cut point(s) define ", length(breaks) + 1, " states, but ", length(states), " state name(s) were provided.")
  }

  ### ---- Conversion ---- ###

  # findInterval() counts the cut points lower than or equal to each value: 0 below the first cut point, up to length(breaks) above the last one
  # (so a value equal to a cut point goes to the state above it). The state is the (count + 1)th, from the lowest values; NA values stay NA
  output <- states[findInterval(values, breaks) + 1L]
  names(output) <- names(values)

  return(output)
}
