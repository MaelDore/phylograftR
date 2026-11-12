##### Internal helpers to handle the lists of both types of taxa (taxa and missing clades) #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


#' @title Detect a list of both types of taxa: `list(taxa, clades)`
#'
#' @description The builders return a `missing_taxa_list` (or a `total_richness_per_clades`) for each type of taxa when both types are provided:
#'   a list with the elements `taxa` (with `$Target_clade`, grafted with `graft_missing_taxa()`) and `clades` (with `$Sister_clade`, grafted with `graft_missing_clades()`).
#'   This function detects such a list, so the checking and extracting functions can handle both types in the same call.
#'   Elements that are `NULL` are dropped: a list with a single element is then handled as this single type.
#'
#' @param x The input given by the user (`missing_taxa_list` or `total_richness_per_clades`).
#' @param arg Name of the argument, for the messages.
#'
#' @noRd
#'
#' @return A list with the elements `typed` (whether `x` is a list named after the types), `both` (whether the two types are present)
#'   and `lists` (the elements of `x` that are not `NULL`, named `"taxa"` and/or `"clades"`, in this order; `NULL` if `typed = FALSE`).
#'
split_typed_list <- function (x, arg)
{
  type_names <- c("taxa", "clades")
  typed <- is.list(x) && length(x) > 0 && !is.null(names(x)) && all(names(x) %in% type_names) && !anyDuplicated(names(x))
  if (!typed) { return(list(typed = FALSE, both = FALSE, lists = NULL)) }

  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0) { stop("'", arg, "' has the elements 'taxa' and 'clades', but both are NULL.", call. = FALSE) }
  x <- x[intersect(type_names, names(x))]

  # Each element must hold the clade item of its type
  clade_items <- c(taxa = "Target_clade", clades = "Sister_clade")
  for (type in names(x))
  {
    item <- clade_items[[type]]
    other <- setdiff(clade_items, item)
    if (!is.list(x[[type]]) || is.null(x[[type]][[item]]) || !is.null(x[[type]][[other]]))
    {
      stop("'", arg, "$", type, "' must be a list with the item '", item, "' (and not '", other, "'): ",
           if (type == "taxa") "the taxa grafted with graft_missing_taxa()" else "the missing clades grafted with graft_missing_clades()", ".", call. = FALSE)
    }
  }
  list(typed = TRUE, both = length(x) == 2, lists = x)
}


#' @title Print the title of the part of a check made for one type of taxa
#'
#' @param type `"taxa"` or `"clades"`.
#'
#' @noRd
#'
#' @return `NULL` (invisibly).
#'
print_type_title <- function (type)
{
  cat(paste0("=== ", if (type == "taxa") "Missing taxa (graft_missing_taxa())" else "Missing clades (graft_missing_clades())", " ===\n"))
  invisible(NULL)
}
