
##### Internal helpers used to handle attributeèbased grafting with multiple traits #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Rationale: the grafting machinery follows ONE simmap.
# To use several traits at once as constraints, the maps of the traits are merged into a composite simmap:
# each branch is cut at the union of the change points of all the traits, and each piece gets a composite state "<state trait 1>|<state trait 2>|...".
# Only the combinations present on the backbone are states of the composite simmap. The allowed states of a taxon (one set of states per trait) are
# resolved into the composite states that belong to all the sets (a state of every trait is allowed).
# Once the grating is achieved on the composite map and all composite states tracked accordingly, the updated composite simmap is split back into one simmap per trait.


# Separator of the states of the traits in the composite states such as "<state trait 1>|<state trait 2>|..."
trait_state_separator <- "|"


#' @title Get the maps of a `"multiSimmap"` and their names
#'
#' @param simmap Object of class `"multiSimmap"`.
#' @param who Character string. Describes the object in the messages.
#'
#' @noRd
#'
#' @return A list with the elements `maps` (list of maps without names) and `names` (character vector, `NULL` if the maps have no names).
#'   The function stops if the maps are partly named or if the names are duplicated.
#'

get_multisimmap_maps <- function (simmap, who = "'simmap'")
{
  maps <- unclass(simmap)
  attributes(maps) <- list(names = names(simmap))
  if (length(maps) == 0) { stop(who, " is a 'multiSimmap' without any map.") }
  map_names <- names(maps)
  if (is.null(map_names) || all(is.na(map_names) | map_names == "")) { return(list(maps = unname(maps), names = NULL)) }
  unnamed <- which(is.na(map_names) | map_names == "")
  if (length(unnamed) > 0)
  {
    stop("Some maps of ", who, " are named and others are not (maps n°", paste(unnamed, collapse = ", "), "): name all the maps or none of them.")
  }
  duplicated_names <- unique(map_names[duplicated(map_names)])
  if (length(duplicated_names) > 0)
  {
    stop("The maps of ", who, " must have different names. Duplicated names: ", paste(duplicated_names, collapse = ", "), ".")
  }
  list(maps = unname(maps), names = map_names)
}


#' @title Check a list of maps of several traits
#'
#' @description Checks the structure of a plain list of maps given as `simmap`: one element per trait, each a `"simmap"` or a `"multiSimmap"`.
#'   All traits must have the same number of maps (a single `"simmap"` counts as one map): the maps are paired (map 1 of every trait together, then map 2, ...).
#'   The maps are only listed here: they are validated against the backbone when they are merged (see `merge_trait_simmaps()`).
#'
#' @param simmap Plain list (named after the traits or not) of `"simmap"` and `"multiSimmap"` objects.
#'
#' @noRd
#'
#' @return As `check_simmap_input()`: a list with the elements `maps` (list with one element per map of class `"phylograftR_trait_maps"`:
#'   a list of the simmaps of the traits named after the traits), `names`, `multiple` and `traits` (names of the traits).
#'

check_trait_maps <- function (simmap)
{
  n_traits <- length(simmap)
  if (n_traits < 2)
  {
    stop("'simmap' is a list of ", n_traits, " element(s): a list of maps is used to provide the maps of several traits (at least 2).\n",
         "For a single trait, provide the 'simmap' (or the 'multiSimmap') itself, not a list.")
  }
  trait_names <- names(simmap)
  if (is.null(trait_names) || all(is.na(trait_names) | trait_names == ""))
  {
    trait_names <- paste0("trait_", seq_len(n_traits))
  } else {
    unnamed <- which(is.na(trait_names) | trait_names == "")
    if (length(unnamed) > 0)
    {
      stop("Some traits of 'simmap' are named and others are not (traits n°", paste(unnamed, collapse = ", "), "): name all the traits or none of them.")
    }
    duplicated_names <- unique(trait_names[duplicated(trait_names)])
    if (length(duplicated_names) > 0)
    {
      stop("The traits of 'simmap' must have different names. Duplicated names: ", paste(duplicated_names, collapse = ", "), ".")
    }
    forbidden <- grepl("[;,=]", trait_names) | grepl(trait_state_separator, trait_names, fixed = TRUE)
    if (any(forbidden))
    {
      stop("The names of the traits of 'simmap' cannot contain ';', ',', '=' or '", trait_state_separator, "': ", paste0("'", trait_names[forbidden], "'", collapse = ", "), ".")
    }
  }

  # Each trait: a simmap (one map) or a multiSimmap
  is_multi <- vapply(simmap, inherits, logical(1), what = "multiSimmap")
  is_single <- vapply(simmap, function (x) inherits(x, "simmap") && !inherits(x, "multiSimmap"), logical(1))
  invalid <- which(!is_multi & !is_single)
  if (length(invalid) > 0)
  {
    stop("Each element of 'simmap' (the maps of one trait) must be a stochastic character map ('simmap') or a set of maps ('multiSimmap'). ",
         "Not valid: ", paste0("'", trait_names[invalid], "'", collapse = ", "), ".")
  }

  trait_maps <- vector("list", n_traits)
  trait_map_names <- vector("list", n_traits)
  for (t in seq_len(n_traits))
  {
    if (is_multi[t])
    {
      got <- get_multisimmap_maps(simmap[[t]], who = paste0("trait '", trait_names[t], "'"))
      trait_maps[[t]] <- got$maps
      trait_map_names[t] <- list(got$names)
    } else {
      trait_maps[[t]] <- list(simmap[[t]])
    }
  }

  # The maps of the traits are paired: same number of maps for every trait
  n_maps_trait <- lengths(trait_maps)
  if (length(unique(n_maps_trait)) > 1)
  {
    stop("All the traits must have the same number of maps, since the maps are paired (map 1 of every trait, then map 2, ...). Number of maps per trait: ",
         paste0("'", trait_names, "' = ", n_maps_trait, collapse = ", "), ".")
  }
  n_maps <- n_maps_trait[1]
  multiple <- any(is_multi)

  # Names of the maps: the names of the 'multiSimmap' if they have some (the same for all the traits), otherwise map_1, map_2, ...
  map_names <- NULL
  if (multiple)
  {
    named <- which(is_multi & !vapply(trait_map_names, is.null, logical(1)))
    if (length(named) > 0)
    {
      if (length(named) < sum(is_multi) || any(vapply(trait_map_names[named], function (x) !identical(x, trait_map_names[[named[1]]]), logical(1))))
      {
        stop("The maps of the traits are paired, so their names must be the same for every trait ('multiSimmap'), or absent for all: the maps of 'simmap' are named differently across the traits.")
      }
      map_names <- trait_map_names[[named[1]]]
    } else {
      map_names <- paste0("map_", seq_len(n_maps))
    }
  }

  maps <- lapply(seq_len(n_maps), function (k)
  {
    structure(stats::setNames(lapply(trait_maps, function (x) x[[k]]), trait_names), class = "phylograftR_trait_maps")
  })
  list(maps = maps, names = map_names, multiple = multiple, traits = trait_names)
}


#' @title Create a simmap from maps
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param maps List with one named vector per branch (in the order of `phy$edge`): durations of the pieces from the rootward to the tipward node (names = states).
#' @param states Character vector. States of the simmap (columns of `mapped.edge`). The states of the maps that are not in `states` are added.
#'
#' @noRd
#'
#' @return An object of class `c("simmap", "phylo")`.
#'

create_simmap <- function (phy, maps, states)
{
  n_edges <- nrow(phy$edge)
  names(maps) <- NULL
  states <- unique(c(states, unlist(lapply(maps, names), use.names = FALSE)))
  mapped_edge <- matrix(0, nrow = n_edges, ncol = length(states), dimnames = list(paste(phy$edge[, 1], phy$edge[, 2], sep = ","), states))
  cell <- (match(unlist(lapply(maps, names), use.names = FALSE), states) - 1) * n_edges + rep(seq_len(n_edges), lengths(maps))
  total <- tapply(unlist(maps, use.names = FALSE), cell, sum)
  mapped_edge[as.integer(names(total))] <- total
  simmap <- phy
  simmap$maps <- maps
  simmap$mapped.edge <- mapped_edge
  class(simmap) <- c("simmap", "phylo")
  simmap
}


#' @title Merge the simmaps of several traits into a composite simmap
#'
#' @description Validates the simmap of each trait against the backbone (`validate_simmap()`), cuts each branch at the union of the change points of the traits
#'   and gives each piece the composite state `"<state trait 1>|<state trait 2>|..."`. Adjacent pieces with the same composite state are merged.
#'   Change points closer than `tol` (relative to the longest branch) are considered identical.
#'
#' @param phy Backbone phylogeny of class `"phylo"` (the simmap is built on it, with its branches in the same order).
#' @param trait_maps List of `"simmap"` objects named after the traits (one map per trait).
#' @param tol Numeric. Relative tolerance. Default is `1e-6`.
#'
#' @noRd
#'
#' @return A list with the elements `simmap` (the composite simmap, class `"simmap"`) and `traits` (a list with the names of the traits `names`, their states `states`
#'   (list of character vectors) and the separator `sep`).
#'

merge_trait_simmaps <- function (phy, trait_maps, tol = 1e-6)
{
  trait_names <- names(trait_maps)
  n_traits <- length(trait_maps)
  sep <- trait_state_separator

  checked <- lapply(seq_len(n_traits), function (t)
  {
    tryCatch(validate_simmap(simmap = trait_maps[[t]], phy = phy, tol = tol),
             error = function (e) stop("Trait '", trait_names[t], "': ", conditionMessage(e), call. = FALSE))
  })
  for (t in seq_len(n_traits))
  {
    bad <- grep(sep, checked[[t]]$states, fixed = TRUE, value = TRUE)
    if (length(bad) > 0)
    {
      stop("Trait '", trait_names[t], "': the states cannot contain '", sep, "' (separator of the states of the traits in the merged map): ", paste0("'", bad, "'", collapse = ", "), ".")
    }
  }

  tol_abs <- tol * max(phy$edge.length)
  edge_pieces <- lapply(seq_len(n_traits), function (t) trait_maps[[t]]$maps[checked[[t]]$edge_map])
  n_edges <- nrow(phy$edge)

  merge_edge <- function (e)
  {
    pieces <- lapply(edge_pieces, function (x) x[[e]])
    if (all(lengths(pieces) == 1))
    {
      # Fast path: a single state for every trait on the branch
      return(stats::setNames(phy$edge.length[e], paste(vapply(pieces, function (x) names(x), character(1)), collapse = sep)))
    }
    ends <- lapply(pieces, cumsum)
    cuts <- c(0, sort(unlist(ends, use.names = FALSE)))
    cuts <- cuts[c(TRUE, diff(cuts) > tol_abs)][-1] # Change points of any trait (identical ones merged); the last is the end of the branch
    starts <- c(0, cuts[-length(cuts)])
    duration <- cuts - starts
    middle <- starts + duration / 2
    state_trait <- lapply(seq_len(n_traits), function (t)
    {
      pieces_start <- c(0, ends[[t]][-length(ends[[t]])])
      names(pieces[[t]])[findInterval(middle, pieces_start)]
    })
    composite <- do.call(paste, c(state_trait, sep = sep))
    run <- cumsum(c(TRUE, composite[-1] != composite[-length(composite)]))
    duration <- as.numeric(tapply(duration, factor(run, levels = unique(run)), sum))
    stats::setNames(duration * phy$edge.length[e] / sum(duration), composite[!duplicated(run)])
  }
  maps <- lapply(seq_len(n_edges), merge_edge)

  list(simmap = create_simmap(phy = phy, maps = maps, states = character(0)),
       traits = list(names = trait_names, states = stats::setNames(lapply(checked, function (x) x$states), trait_names), sep = sep))
}


#' @title Split a composite simmap into one simmap per trait
#'
#' @description Inverse of `merge_trait_simmaps()`, applied to the updated composite simmap of the grafted tree:
#'   for each trait, the state of each piece is the component of its composite state, and adjacent pieces of a branch in the same state are merged.
#'   A grafted branch is entirely in the composite state it inherited, so each trait inherits the state of the lineage it was grafted on.
#'   The states of each trait are those of its original simmap.
#'
#' @param simmap Composite simmap (class `"simmap"`), e.g., the output of `build_updated_simmap()`.
#' @param traits The element `traits` returned by `merge_trait_simmaps()`.
#'
#' @noRd
#'
#' @return A list named after the traits, with the simmap of each trait (class `"simmap"`).
#'

split_trait_simmap <- function (simmap, traits)
{
  phy <- simmap
  phy$maps <- NULL
  phy$mapped.edge <- NULL
  class(phy) <- "phylo"
  n_edges <- length(simmap$maps)
  edge_of_piece <- rep(seq_len(n_edges), lengths(simmap$maps))
  composite <- unlist(lapply(simmap$maps, names), use.names = FALSE)
  duration <- unlist(simmap$maps, use.names = FALSE)
  parts <- strsplit(composite, traits$sep, fixed = TRUE)

  result <- lapply(seq_along(traits$names), function (t)
  {
    state <- vapply(parts, function (p) p[t], character(1))
    maps <- lapply(split(seq_along(state), factor(edge_of_piece, levels = seq_len(n_edges))), function (i)
    {
      run <- cumsum(c(TRUE, state[i][-1] != state[i][-length(i)]))
      stats::setNames(as.numeric(tapply(duration[i], factor(run, levels = unique(run)), sum)), state[i][!duplicated(run)])
    })
    create_simmap(phy = phy, maps = maps, states = traits$states[[t]])
  })
  names(result) <- traits$names
  result
}


#' @title Parse the allowed states of several traits given as text
#'
#' @description Converts a text such as `"A1,A2;;C3"` (the states of different traits separated by `;`, the states of a trait by `,`) into a list with one element per trait.
#'   A trait without states (empty field) gives `NULL` (the global value is used, see `resolve_trait_allowed_states()`).
#'   The traits can be named: `"color=A1,A2; elevation=C3"` (all the fields or none).
#'
#' @param cell Character string.
#' @param who Character string. Describes where the text comes from, for the messages.
#'
#' @noRd
#'
#' @return A list with one element per field (character vector of states, or `NULL`), named after the traits if they are named.
#'

parse_trait_states_cell <- function (cell, who = "'allowed_states'")
{
  fields <- trimws(strsplit(paste0(cell, ";"), ";", fixed = TRUE)[[1]]) # The last field is kept even if empty
  has_name <- grepl("=", fields, fixed = TRUE)
  if (any(has_name) && !all(has_name | fields == ""))
  {
    stop(who, " mixes named and unnamed traits (name all the traits, 'trait=state1,state2', or none): ", cell)
  }
  trait_name <- rep("", length(fields))
  if (any(has_name))
  {
    trait_name[has_name] <- trimws(sub("=.*$", "", fields[has_name]))
    fields[has_name] <- sub("^[^=]*=", "", fields[has_name])
    if (any(trait_name[has_name] == "")) { stop(who, " has a trait without name before '=': ", cell) }
    if (anyDuplicated(trait_name[has_name])) { stop(who, " has repeated traits (", paste(unique(trait_name[has_name][duplicated(trait_name[has_name])]), collapse = ", "), "): ", cell) }
    fields <- fields[has_name]
    trait_name <- trait_name[has_name]
  }
  states <- lapply(fields, function (field)
  {
    if (field == "") { return(NULL) }
    values <- trimws(strsplit(field, ",", fixed = TRUE)[[1]])
    if (any(values == "")) { stop(who, " has an empty state (check the ',' separators): ", cell) }
    if (anyDuplicated(values)) { stop(who, " has repeated states (", paste(unique(values[duplicated(values)]), collapse = ", "), "): ", cell) }
    values
  })
  if (any(nzchar(trait_name))) { names(states) <- trait_name }
  states
}


#' @title Allowed states of each trait: normalise the input
#'
#' @description Converts the value of `allowed_states` (global, or of a taxon) into a list with one element per trait (in the order of the traits):
#'   a character vector of states, or `NULL` (no states: the global value, or all the states of the trait, is used).
#'   The input is a list with one element per trait, in the order of the traits (e.g., `list(c("A1", "A2"), NULL, "C3")`)
#'   or named after the traits (traits absent from the list have no states), or a text parsed by `parse_trait_states_cell()`.
#'
#' @param value `NULL`, list, or character string with `;` or `=`.
#' @param trait_names Character vector. Names of the traits.
#' @param who Character string. Describes where the value comes from, for the messages.
#'
#' @noRd
#'
#' @return A list named after the traits.
#'

normalize_trait_states <- function (value, trait_names, who)
{
  n_traits <- length(trait_names)
  empty <- stats::setNames(vector("list", n_traits), trait_names)
  if (is.null(value) || length(value) == 0 || (is.atomic(value) && length(value) == 1 && is.na(value))) { return(empty) }
  if (is.character(value))
  {
    if (length(value) == 1 && grepl("[;=]", value))
    {
      value <- parse_trait_states_cell(value, who = who)
    } else {
      stop("The maps of ", n_traits, " traits are used, so the allowed states of ", who, " must be given for each trait: a list with one element per trait in the order of the traits (",
           paste(trait_names, collapse = ", "), "), e.g., list(c('state_1', 'state_2'), NULL, 'state_3') ('NULL': no restriction on the trait), ",
           "or a text such as 'state_1,state_2;;state_3'.")
    }
  }
  if (!is.list(value)) { stop("The allowed states of ", who, " must be a list with one element per trait.") }
  value_names <- names(value)
  if (!is.null(value_names) && any(value_names != ""))
  {
    if (any(is.na(value_names) | value_names == "")) { stop("The allowed states of ", who, " name some traits but not all of them.") }
    if (anyDuplicated(value_names)) { stop("The allowed states of ", who, " repeat traits: ", paste(unique(value_names[duplicated(value_names)]), collapse = ", "), ".") }
    unknown <- setdiff(value_names, trait_names)
    if (length(unknown) > 0)
    {
      stop("The allowed states of ", who, " refer to traits that are not in 'simmap': ", paste0("'", unknown, "'", collapse = ", "), ". Traits: ", paste0("'", trait_names, "'", collapse = ", "), ".")
    }
    empty[value_names] <- value
    value <- empty
  } else {
    if (length(value) != n_traits)
    {
      stop("The allowed states of ", who, " have ", length(value), " element(s) but there are ", n_traits, " traits (", paste(trait_names, collapse = ", "), "): ",
           "give one element per trait in this order ('NULL' for no restriction), or name the elements after the traits.")
    }
    names(value) <- trait_names
  }
  lapply(value, function (x)
  {
    if (length(x) == 0 || (length(x) == 1 && is.na(x))) { return(NULL) }
    if (!is.character(x) || anyNA(x)) { stop("The states of each trait in ", who, " must be character vectors without NA (or NULL).") }
    unique(x)
  })
}


#' @title Resolve the allowed states of each trait for each taxon
#'
#' @description For each taxon, the allowed states of each trait are the states of `missing_taxa_list$allowed_states` for this trait if provided,
#'   otherwise the global `allowed_states` for this trait, otherwise all the states of the trait (no restriction). A taxon-level value replaces the global one (no intersection), trait by trait.
#'   The function stops if a state is not in the states of its trait. Used by `resolve_trait_allowed_states()` (grafting, on the composite simmap) and by `check_grafting_constraints()` (check, trait by trait).
#'
#' @param traits List with `names` (trait names) and `states` (list of the states of each trait, named after the traits).
#' @inheritParams resolve_trait_allowed_states
#'
#' @noRd
#'
#' @return A list named after `taxon_names`, each element a list named after the traits with the allowed states (character vector) of the trait.
#'

resolve_trait_state_sets <- function (allowed_states_taxa = NULL, allowed_states = NULL, traits, taxon_names)
{
  trait_names <- traits$names
  n_traits <- length(trait_names)
  n_taxa <- length(taxon_names)

  check_trait_states <- function (spec, who)
  {
    for (t in seq_len(n_traits))
    {
      if (!is.null(spec[[t]])) { check_allowed_states(spec[[t]], traits$states[[t]], who = paste0(who, " (trait '", trait_names[t], "')")) }
    }
  }

  global <- normalize_trait_states(allowed_states, trait_names, who = "the global 'allowed_states' argument")
  check_trait_states(global, "the global 'allowed_states' argument")

  if (!is.null(allowed_states_taxa) && !is.list(allowed_states_taxa))
  {
    stop("'missing_taxa_list$allowed_states' must be a list (named after 'Taxon_name').")
  }
  if (!is.null(allowed_states_taxa) && length(allowed_states_taxa) != n_taxa)
  {
    stop("'missing_taxa_list$allowed_states' must have the same length as 'Taxon_name' (", length(allowed_states_taxa), " vs ", n_taxa, ").")
  }

  sets_all <- vector("list", n_taxa)
  names(sets_all) <- taxon_names
  for (i in seq_len(n_taxa))
  {
    who <- paste0("'missing_taxa_list$allowed_states' of taxon '", taxon_names[i], "'")
    spec <- normalize_trait_states(if (is.null(allowed_states_taxa)) NULL else allowed_states_taxa[[i]], trait_names, who = who)
    check_trait_states(spec, who)
    sets_all[[i]] <- stats::setNames(lapply(seq_len(n_traits), function (t)
    {
      if (!is.null(spec[[t]])) spec[[t]] else if (!is.null(global[[t]])) global[[t]] else traits$states[[t]]
    }), trait_names)
  }
  sets_all
}


#' @title Resolve the allowed states of each taxon with the maps of several traits
#'
#' @description Equivalent of `resolve_allowed_states()` for a composite simmap. For each taxon, the allowed states of each trait are the states of `missing_taxa_list$allowed_states`
#'   for this trait if provided, otherwise the global `allowed_states` for this trait, otherwise all the states of the trait (no restriction).
#'   A taxon-level value replaces the global one (no intersection), trait by trait.
#'   The allowed composite states are the composite states observed on the backbone whose state of each trait is allowed.
#'   The function stops if a state is not in the map of its trait, or if no composite state is allowed for a taxon.
#'
#' @param allowed_states_taxa `NULL`, or list (named after `Taxon_name`) of the allowed states of each taxon (see `normalize_trait_states()`): `missing_taxa_list$allowed_states`.
#' @param allowed_states `NULL`, or list with one element per trait (see `normalize_trait_states()`). The global `allowed_states` argument.
#' @param traits The element `traits` returned by `merge_trait_simmaps()`.
#' @param composite_states Character vector. Composite states of the merged simmap.
#' @param taxon_names Character vector. `missing_taxa_list$Taxon_name`.
#'
#' @noRd
#'
#' @return A list named after `taxon_names`, with the allowed composite states (character vector) of each taxon.
#'

resolve_trait_allowed_states <- function (allowed_states_taxa = NULL, allowed_states = NULL, traits, composite_states, taxon_names)
{
  trait_names <- traits$names
  n_traits <- length(trait_names)
  n_taxa <- length(taxon_names)
  sets_all <- resolve_trait_state_sets(allowed_states_taxa = allowed_states_taxa, allowed_states = allowed_states, traits = traits, taxon_names = taxon_names)

  # State of each trait in each composite state
  parts <- matrix(unlist(strsplit(composite_states, traits$sep, fixed = TRUE)), ncol = n_traits, byrow = TRUE)

  resolved <- vector("list", n_taxa)
  names(resolved) <- taxon_names
  for (i in seq_len(n_taxa))
  {
    sets <- sets_all[[i]]
    allowed <- Reduce(`&`, lapply(seq_len(n_traits), function (t) parts[, t] %in% sets[[t]]))
    if (!any(allowed))
    {
      stop("No combination of the allowed states of taxon '", taxon_names[i], "' is found on the maps: ",
           paste0(trait_names, " = ", vapply(sets, paste, character(1), collapse = ", "), collapse = "; "), ".\n",
           "Combinations of states found on the maps (", paste(trait_names, collapse = " | "), "): ",
           paste0("'", utils::head(composite_states, 10), "'", collapse = ", "), if (length(composite_states) > 10) ", ..." else "", ".")
    }
    resolved[[i]] <- composite_states[allowed]
  }
  resolved
}
