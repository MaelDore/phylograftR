
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
#'   **A `tolerance` can be applied to the values to account for their uncertainty, or for the uncertainty of the cut points**:
#'   each value is then assigned to all the states covered by the interval `[value - tolerance, value + tolerance]`.
#'
#' @param values Numeric vector. The trait values to convert, in the same units as the trait used to build the `simmap`. If named, the names are kept in the output (e.g., the names of the missing taxa).
#'   `NA` values give `NA` states.
#' @param simmap A `"simmap"` or `"multiSimmap"` produced by [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap],
#'   which carries the cut points and the names of the states in its attributes `"breaks"` and `"states"`.
#'   Alternatively, leave it as `NULL` and provide `breaks` and `states` (e.g., for a `simmap` built from tip values binned by yourself).
#' @param breaks Numeric vector. Only used if `simmap = NULL`. The interior cut points between the states, in increasing order and different from each other.
#' @param states Character vector. Only used if `simmap = NULL`. The names of the states, from the lowest to the highest values.
#'   **Must have one more element than `breaks`**, because N interior cuts (`breaks`) produce N+1 groups.
#' @param tolerance Numeric, zero or positive. Default is `0`: each value is assigned to a single state.
#'   With a tolerance, each value is assigned to **all the states covered by the interval `[value - tolerance, value + tolerance]`**,
#'   in the units of the trait. Values that are farther than `tolerance` from every cut point are still assigned to a single state.
#'   It can be a single number used for all the values, or a vector with one tolerance per value (e.g., the measurement error of each missing taxon).
#'   See 'Details' for more information.
#'
#' @details Each state covers the values from its lower cut point (included) up to its upper cut point (excluded): **a value equal to a cut point belongs to the state above it.**
#'   The first and last states are open-ended, so values below the first cut point or above the last one are assigned to the lowest and highest states.
#'
#'   ### Tolerance
#'
#'   **By default (`tolerance = 0`), each value is assigned to a unique state.**
#'   To include uncertainty, `tolerance` assigns **several states** to the values that are close to a cut point, as accepted by `allowed_states` in the grafting functions:
#'   a value gets the state of `value - tolerance`, the state of `value + tolerance`, and all the states in between.
#'   The same rule applies at the edges of the interval as for the values: if `value + tolerance` is equal to a cut point, the state above this cut point is included
#'   (and if `value - tolerance` is equal to a cut point, the state above it is the first state of the value).
#'   A tolerance makes the constraint of the grafting weaker: the taxa can be grafted on a wider set of branches (see [`phylograftR::check_grafting_plan()`][check_grafting_plan]
#'   to evaluate how much, with the `valid_fraction`).
#'
#' @return **A list of character vectors, with one element per value**, named as in `values` (unnamed if `values` is not named).
#'   Each element contains the states of the value, ordered from the lowest to the highest values: a single state by default,
#'   and several states for the values close to a cut point if a `tolerance` is used. The state of a `NA` value is `NA`.
#'
#'   The output can be used directly as the `allowed_states` item of a list of missing taxa, if `values` is named after the taxa.
#'   To get a simple character vector with a single state per value, e.g., to compare the states of tips or nodes, use `unlist()`
#'   (**only when `tolerance = 0`**, as the elements have then one state each. With a tolerance, the elements can have several states).
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap] to obtain the `simmap` and its cut points from a continuous trait map (`contMap`).
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa based on target clades and excluded subclades, using the states to constrain the grafting.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades, using the states to constrain the grafting.
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
#' # The output is a list with the states of each value
#' my_states <- convert_values_to_states(missing_values, simmap = my_simmap)
#' my_states
#'
#' # To get a character vector with a single state per value, use unlist() (if tolerance = 0)
#' unlist(my_states)
#'
#' # Same conversion with the cut points and states given directly
#' convert_values_to_states(missing_values, breaks = attr(my_simmap, "breaks"),
#'                          states = c("low", "medium", "high"))
#'
#' # --- Account for the uncertainty of the values with a tolerance --- #
#'
#' # Taxon_2 has a value just above the second cut point, between the states "medium" and "high"
#' close_values <- c(Taxon_1 = -2, Taxon_2 = attr(my_simmap, "breaks")[2] + 0.1, Taxon_3 = 3)
#' convert_values_to_states(close_values, simmap = my_simmap)
#'
#' # With a tolerance of 0.3 trait units, values close to a cut point get several states
#' allowed_states <- convert_values_to_states(close_values, simmap = my_simmap, tolerance = 0.3)
#' allowed_states
#'
#' # One tolerance per taxon, e.g., the measurement error of each value (no tolerance for Taxon_1 and Taxon_3)
#' convert_values_to_states(close_values, simmap = my_simmap, tolerance = c(0, 0.3, 0))
#'
#' # The output can be used directly as the 'allowed_states' of a list of missing taxa
#' missing_taxa_list <- list(Taxon_name = names(close_values),
#'                           Target_clade = list(Taxon_1 = "t1", Taxon_2 = c("t1", "t2"), Taxon_3 = "t1"))
#' missing_taxa_list$allowed_states <- allowed_states
#'

convert_values_to_states <- function (values, # Numeric trait values to convert, named after the taxa if wanted
                                      simmap = NULL, # simmap or multiSimmap from convert_contMap_to_simmap(), with attributes "breaks" and "states"
                                      breaks = NULL, # Interior cut points, if simmap = NULL
                                      states = NULL, # Names of the states from the lowest to the highest values, if simmap = NULL
                                      tolerance = 0) # Tolerance around the values: all the states covered by [value - tolerance, value + tolerance] are returned
{
  ### ---- Checks ---- ###

  if (!is.numeric(values)) { stop("'values' must be a numeric vector of trait values.") }
  if (any(is.infinite(values))) { stop("'values' must not contain infinite values.") }
  if (!is.numeric(tolerance) || length(tolerance) < 1 || anyNA(tolerance) || any(!is.finite(tolerance)) || any(tolerance < 0))
  {
    stop("'tolerance' must be a numeric vector of finite values that are zero or positive.")
  }
  if (length(tolerance) != 1 && length(tolerance) != length(values))
  {
    stop("'tolerance' must have a single value, or one value per element of 'values' (", length(values), "), but it has ", length(tolerance), " values.")
  }

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
  # (so a value equal to a cut point goes to the state above it). The state is the (count + 1)th, from the lowest values.
  # With a tolerance: all the states from the state of (value - tolerance) to the state of (value + tolerance), ordered from the lowest to the highest values.
  # The states are in the order of the cut points, so the states between the two are the consecutive ones. Without tolerance, both are the same state.
  # NA values give NA. The output is always a list with one element per value
  tolerance <- rep_len(tolerance, length(values))
  first_state <- findInterval(values - tolerance, breaks) + 1L
  last_state <- findInterval(values + tolerance, breaks) + 1L
  output <- lapply(seq_along(values), function (i)
  {
    if (is.na(values[i])) { return(NA_character_) }
    states[first_state[i]:last_state[i]]
  })
  names(output) <- names(values)

  return(output)
}
