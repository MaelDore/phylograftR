
##### Internal helpers dealing with biogeographic ranges used as states of a simmap #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Coding, areas and matching of biogeographic ranges #####

# Rationale
#   The states of a simmap can be biogeographic ranges, coded as in BioGeoBEARS: one CAPITAL letter per area, and a multi-area range is
#   the string of the letters of its areas, each area once, in alphabetical order (e.g., "A", "AB", "ABC").
#   With 'treat_states_as_areas = TRUE', a range is a SET of areas: an allowed range R matches a state S of the simmap if R is contained in S
#   (e.g., "A" matches "A", "AB" and "ABC"; "AB" matches "AB" and "ABC", but not "A").
#   The helpers below only deal with the coding and the matching of ranges. They do not use phylogenies or simmaps.


#' @title Split ranges in areas
#'
#' @description Split biogeographic ranges in their areas, one letter per area (e.g., `"AB"` gives `c("A", "B")`).
#'   No check of the coding is done: see `check_range_coding()`.
#'
#' @param ranges Character vector. Ranges (or states).
#'
#' @noRd
#'
#' @return A list of character vectors (the areas of each range), in the order of `ranges`.
#'

split_range <- function (ranges)
{
  strsplit(ranges, split = "", fixed = TRUE)
}


#' @title Check the coding of biogeographic ranges
#'
#' @description Check that ranges are coded as in BioGeoBEARS: only capital letters (one per area), each area once, and the areas in alphabetical order.
#'   The reasons why a range is not correctly coded are all reported (e.g., lowercase letters and a repeated area).
#'
#' @param ranges Character vector. Ranges (or states) to check. Duplicated values are checked once.
#'
#' @noRd
#'
#' @return A data.frame with one row per unique range and 3 columns:
#'  * `range` Character. The range.
#'  * `valid` Logical. Is the range correctly coded?
#'  * `reason` Character. Why the range is not correctly coded. Empty string if it is.
#'

check_range_coding <- function (ranges)
{
  ranges <- unique(ranges)
  get_reason <- function (range)
  {
    if (is.na(range) || !nzchar(range)) { return("empty or NA") }
    areas <- strsplit(range, split = "", fixed = TRUE)[[1]]
    reasons <- character(0)
    lowercase <- areas[grepl("[a-z]", areas)]
    not_letters <- areas[!grepl("[A-Za-z]", areas)]
    if (length(lowercase) > 0) { reasons <- c(reasons, paste0("lowercase letters (", paste(unique(lowercase), collapse = ", "), ")")) }
    if (length(not_letters) > 0) { reasons <- c(reasons, paste0("characters that are not letters (", paste0("'", unique(not_letters), "'", collapse = ", "), ")")) }
    repeated <- unique(areas[duplicated(areas)])
    if (length(repeated) > 0) { reasons <- c(reasons, paste0("repeated areas (", paste(repeated, collapse = ", "), ")")) }
    if (length(lowercase) == 0 && length(not_letters) == 0)
    {
      expected <- sort(unique(areas), method = "radix")
      if (!identical(unique(areas), expected))
      {
        reasons <- c(reasons, paste0("areas not in alphabetical order (expected '", paste(expected, collapse = ""), "')"))
      }
    }
    paste(reasons, collapse = "; ")
  }
  reason <- vapply(ranges, get_reason, character(1), USE.NAMES = FALSE)
  data.frame(range = ranges, valid = !nzchar(reason), reason = reason, stringsAsFactors = FALSE)
}


#' @title Stop if ranges are not coded as in BioGeoBEARS
#'
#' @description Strict validation used with `treat_states_as_areas = TRUE`: stops with an informative message listing each incorrect range
#'   and the reason why it is not correctly coded (see `check_range_coding()`).
#'
#' @param ranges Character vector. Ranges (or states) to validate.
#' @param who Character string. Describes where the ranges come from, used in the message (e.g., "the simmap", "the global 'allowed_states' argument").
#'
#' @noRd
#'
#' @return `ranges` (invisibly).
#'

validate_range_coding <- function (ranges, who = "the simmap")
{
  if (!is.character(ranges) || length(ranges) == 0)
  {
    stop("The states in ", who, " must be a non-empty character vector.")
  }
  checked <- check_range_coding(ranges)
  incorrect <- checked[!checked$valid, , drop = FALSE]
  if (nrow(incorrect) > 0)
  {
    stop("With 'treat_states_as_areas = TRUE', the states in ", who, " must be biogeographic ranges coded as in BioGeoBEARS: ",
         "one capital letter per area, each area once, and the areas of a multi-area range in alphabetical order (e.g., 'A', 'AB', 'ABC').\n",
         "Incorrect states: ", paste0("'", incorrect$range, "' (", incorrect$reason, ")", collapse = "; "), ".")
  }
  invisible(ranges)
}


#' @title Get the areas of a set of ranges
#'
#' @description Get the unique areas found in ranges (e.g., `c("A", "AB", "BC")` gives `c("A", "B", "C")`).
#'
#' @param ranges Character vector. Ranges (or states).
#'
#' @noRd
#'
#' @return A character vector of the unique areas, in alphabetical order.
#'

get_range_areas <- function (ranges)
{
  sort(unique(unlist(split_range(ranges))), method = "radix")
}


#' @title Check that the areas of the allowed ranges exist in the simmap
#'
#' @description Stops with an informative message if some areas of the allowed ranges are not in any state of the simmap (e.g., typo),
#'   listing the areas available. Same logic as `check_allowed_states()`, at the level of the areas.
#'
#' @param allowed_ranges Character vector. Allowed ranges to check.
#' @param states Character vector. States of the simmap.
#' @param who Character string. Describes where the ranges come from, used in the message (e.g., "the global 'allowed_states' argument").
#'
#' @noRd
#'
#' @return `allowed_ranges` (invisibly).
#'

check_allowed_areas <- function (allowed_ranges, states, who = "'allowed_states'")
{
  unknown <- setdiff(get_range_areas(allowed_ranges), get_range_areas(states))
  if (length(unknown) > 0)
  {
    stop("The areas of the states in ", who, " are not in the simmap: ", paste0("'", unknown, "'", collapse = ", "), ".\n",
         "Areas available in the simmap: ", paste0("'", get_range_areas(states), "'", collapse = ", "), ".")
  }
  invisible(allowed_ranges)
}


#' @title Is a range contained in states?
#'
#' @description Test if all the areas of a range are in each of a set of states (R is contained in S).
#'
#' @param range Character string. One range (e.g., `"AB"`).
#' @param states Character vector. States to test (e.g., the states of a simmap).
#'
#' @noRd
#'
#' @return A logical vector, in the order of `states`: `TRUE` if `range` is contained in the state.
#'   A range contained in itself is `TRUE`. For example, `"A"` is contained in `"A"`, `"AB"`, `"ABC"`, but not in `"B"` or `"BC"`.
#'

range_in_states <- function (range, states)
{
  areas <- split_range(range)[[1]]
  vapply(split_range(states), function (state_areas) all(areas %in% state_areas), logical(1))
}


#' @title Expand the allowed ranges into acceptable states
#'
#' @description Match each allowed range to the states of the simmap that contain it (see `range_in_states()`).
#'   A state is acceptable if it contains at least one allowed range. For each acceptable state, the allowed ranges it contains are the
#'   candidate ranges that a lineage grafted in this state can inherit.
#'
#' @param allowed_ranges Character vector. Allowed ranges.
#' @param states Character vector. States of the simmap.
#'
#' @noRd
#'
#' @return A list with 4 items:
#'  * `acceptable_states` Character vector. States of the simmap that contain at least one allowed range, in the order of `states`.
#'  * `candidates` Named list (named after `acceptable_states`) of character vectors: the allowed ranges contained in each acceptable state,
#'    in the order of `allowed_ranges`.
#'  * `inherited_ranges` Character vector. All the ranges that can be inherited (unique candidates, in the order of `allowed_ranges`).
#'  * `unmatched_ranges` Character vector. Allowed ranges contained in no state of the simmap (e.g., `"AB"` if the simmap has only `"A"` and `"B"`).
#'

expand_allowed_ranges <- function (allowed_ranges, states)
{
  allowed_ranges <- unique(allowed_ranges)
  states <- unique(states)
  contained <- matrix(vapply(allowed_ranges, range_in_states, logical(length(states)), states = states),
                      nrow = length(states), ncol = length(allowed_ranges), dimnames = list(states, allowed_ranges))
  candidates <- lapply(seq_along(states), function (i) allowed_ranges[contained[i, ]])
  names(candidates) <- states
  candidates <- candidates[lengths(candidates) > 0]
  list(acceptable_states = as.character(names(candidates)),
       candidates = candidates,
       inherited_ranges = allowed_ranges[allowed_ranges %in% unlist(candidates, use.names = FALSE)],
       unmatched_ranges = allowed_ranges[colSums(contained) == 0])
}
