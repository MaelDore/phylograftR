
##### Temporary saves of the grafting functions: checks, files, saved state and resume #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - check_save_temp(), resolve_temp_file_path(), temp_map_tags(), warn_save_temp_every(): checks of the arguments and names of the files
#  - write_temp_file(), read_temp_file(): atomic writing and safe reading of the temporary files
#  - temp_signature(), check_temp_signature(): detect that a file was not produced with the inputs of the current call
#  - save_temp_progress(), save_temp_complete(), read_resume_state(): content of the files (see below), and how to resume from them
#
# A temporary file is the phylogeny grafted so far, as a normal 'phylo' (original node labels restored, so it can be read and plotted),
# with the attribute "phylograftR_temp": a list with
#  - status: "partial" (grafting in progress) or "complete" (the phylogeny is the final output of the simulation, with $simmap, grafted_tips, ...)
#  - signature: fingerprints of the inputs, to check that a resume uses the same inputs
#  - for a partial file only: progress (done, total, unit), node_labels (labels followed during grafting), tracker (states followed during grafting),
#    bd_cache (birth-death models already fitted), counters, rng (state of the random number generator, to continue the same random draws)
# One file per task (simulation x map): see temp_file_path_simulation().


##### 1/ Checks of the arguments #####

#' @title Check the arguments of the temporary saves
#'
#' @param save_temp,resume Logicals.
#' @param save_temp_every `NULL` (default of the grafting functions: only the completed phylogenies are saved), or a single positive whole number.
#' @param clean_temp_files Single logical. Whether the temporary files are deleted once the run is completed.
#'
#' @noRd
#'
#' @return `TRUE` (invisibly). Stops with an informative message otherwise.
#'

check_save_temp <- function (save_temp, save_temp_every, resume, clean_temp_files = TRUE)
{
  if (!is.logical(save_temp) || length(save_temp) != 1 || is.na(save_temp))
  {
    stop("'save_temp' must be TRUE or FALSE.")
  }
  if (!is.logical(resume) || length(resume) != 1 || is.na(resume))
  {
    stop("'resume' must be TRUE or FALSE.")
  }
  if (!is.null(save_temp_every) && (!is.numeric(save_temp_every) || length(save_temp_every) != 1 || !is.finite(save_temp_every) || save_temp_every < 1 || save_temp_every != round(save_temp_every)))
  {
    stop("'save_temp_every' must be a single positive whole number (number of taxa, or of groups of taxa, grafted between two temporary saves), ",
         "or NULL to only save the phylogenies once fully grafted.")
  }
  if (!is.logical(clean_temp_files) || length(clean_temp_files) != 1 || is.na(clean_temp_files))
  {
    stop("'clean_temp_files' must be TRUE or FALSE.")
  }
  if (resume && !save_temp)
  {
    stop("'resume = TRUE' needs 'save_temp = TRUE': the grafting is resumed from the temporary files of a previous run with the same arguments, ",
         "and continues to save its progress in these files.")
  }
  invisible(TRUE)
}


#' @title Absolute path of the temporary file
#'
#' @description The path is made absolute in the main session, so the workers of a parallel run (whose working directory can differ) write to the same place.
#'   `NULL` gives the default file name in the working directory.
#'
#' @param temp_file_path `NULL` or path ending with `".rds"` (checked by `check_temp_file_path()`).
#'
#' @noRd
#'
#' @return A path.
#'

resolve_temp_file_path <- function (temp_file_path)
{
  if (is.null(temp_file_path)) { temp_file_path <- "temp_phy_from_grafting.rds" }
  file.path(normalizePath(dirname(temp_file_path), winslash = "/", mustWork = FALSE), basename(temp_file_path))
}


#' @title Tags of the maps used in the names of the temporary files
#'
#' @description The name of a map with the characters other than letters, digits, `_`, `.` and `-` replaced by `_`.
#'   Maps whose names give the same tag (e.g., `"a b"` and `"a/b"`) get their position added (`"a_b-1"`, `"a_b-2"`), so they never share a file.
#'
#' @param map_names `NULL` (no `multiSimmap`) or character vector.
#'
#' @noRd
#'
#' @return `NULL`, or character vector of tags named after the maps.
#'

temp_map_tags <- function (map_names)
{
  if (is.null(map_names)) { return(NULL) }
  tags <- gsub("[^A-Za-z0-9_.-]", "_", map_names)
  duplicated_tags <- tags %in% tags[duplicated(tags)]
  tags[duplicated_tags] <- paste0(tags[duplicated_tags], "-", which(duplicated_tags))
  stats::setNames(tags, map_names)
}


#' @title Warn if `save_temp_every` is larger than the number of units to graft
#'
#' @description If `save_temp_every` is larger than the number of taxa (or groups of taxa), no partial save happens: the phylogeny is only saved once fully grafted.
#'   Sent once, before any simulation. `save_temp_every = NULL` is the way to ask for it without a warning.
#'
#' @param save_temp_every `NULL` or single positive whole number.
#' @param n_units Integer vector. Number of units to graft, for each map.
#' @param unit Character. Name of the unit, e.g., `"taxa"`.
#' @param map_names `NULL` or names of the maps.
#'
#' @noRd
#'
#' @return `NULL` (invisibly). Sends a warning if `save_temp_every` is larger than the number of units.
#'

warn_save_temp_every <- function (save_temp_every, n_units, unit, map_names = NULL)
{
  if (is.null(save_temp_every)) { return(invisible(NULL)) }
  too_large <- save_temp_every > n_units
  if (!any(too_large)) { return(invisible(NULL)) }
  which_maps <- if (!is.null(map_names) && length(unique(n_units)) > 1) paste0(" for the map(s) ", paste0("'", map_names[too_large], "'", collapse = ", ")) else ""
  largest <- max(n_units[too_large])
  warning("'save_temp_every' (", save_temp_every, ") is larger than the number of ", unit, " to graft (", paste(unique(n_units[too_large]), collapse = " or "), ")", which_maps, ": ",
          "no partial phylogeny is saved, each phylogeny is only saved once fully grafted.\n",
          "Use 'save_temp_every = NULL' to ask for it without this warning, or a value up to ", largest, " (the number of ", unit, ") to save progress along the way.",
          call. = FALSE)
  invisible(NULL)
}


##### 2/ Reading and writing the files #####

#' @title Write a temporary file
#'
#' @description Written to a side file, then renamed, so a crash during the writing never leaves a corrupted temporary file in place of a good one.
#'
#' @param object Object to save.
#' @param path Path of the `.rds` file.
#'
#' @noRd
#'
#' @return `TRUE` (invisibly).
#'

write_temp_file <- function (object, path)
{
  side_file <- paste0(path, ".part")
  tryCatch(saveRDS(object = object, file = side_file),
           error = function (e) stop("The temporary file '", path, "' could not be written: ", conditionMessage(e), call. = FALSE))
  renamed <- suppressWarnings(file.rename(side_file, path))
  if (!renamed)
  {
    # Some systems do not replace an existing file when renaming
    renamed <- file.copy(side_file, path, overwrite = TRUE)
    unlink(side_file)
  }
  if (!renamed) { stop("The temporary file '", path, "' could not be written.", call. = FALSE) }
  invisible(TRUE)
}


#' @title Read a temporary file
#'
#' @param path Path of the `.rds` file.
#'
#' @noRd
#'
#' @return The object saved, or `NULL` if the file does not exist or cannot be read (with a warning in the second case).
#'

read_temp_file <- function (path)
{
  if (!file.exists(path)) { return(NULL) }
  tryCatch(readRDS(path), error = function (e)
  {
    warning("The temporary file '", path, "' could not be read and is ignored (", conditionMessage(e), "): this simulation starts from scratch.", call. = FALSE)
    NULL
  })
}


##### 3/ Content of the files #####

#' @title Fingerprint of an object
#'
#' @description MD5 sum of the serialized object, to compare the inputs of a run with the ones of the run that produced a temporary file without storing them.
#'
#' @param x Any object.
#'
#' @importFrom tools md5sum
#'
#' @noRd
#'
#' @return Character string.
#'

fingerprint <- function (x)
{
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(x, file = file, compress = FALSE, version = 2)
  unname(tools::md5sum(file))
}


#' @title Signature of the inputs of a simulation
#'
#' @description Fingerprints of the inputs that determine the simulation, stored in the temporary file and compared when resuming.
#'
#' @param phy Backbone phylogeny, as provided by the user.
#' @param missing_taxa_list The list of the taxa of the map, after resolution of the arguments and before the random order.
#' @param tracker Initial tracker of the states of the map (`NULL` without simmap).
#' @param options List of the other arguments and settings that change the result (including the seed and the index of the simulation).
#'
#' @noRd
#'
#' @return Named list of fingerprints.
#'

temp_signature <- function (phy, missing_taxa_list, tracker, options)
{
  list(backbone = fingerprint(list(phy$edge, phy$edge.length, phy$tip.label, phy$node.label)),
       taxa = fingerprint(missing_taxa_list),
       states = if (is.null(tracker)) NULL else fingerprint(tracker),
       options = fingerprint(options))
}


#' @title Check that a temporary file was produced with the inputs of the current call
#'
#' @param stored,current Signatures (see `temp_signature()`).
#' @param path Path of the file, for the message.
#'
#' @noRd
#'
#' @return `TRUE` (invisibly). Stops with an informative message if the signatures differ.
#'

check_temp_signature <- function (stored, current, path)
{
  differs <- names(current)[!vapply(names(current), function (name) identical(stored[[name]], current[[name]]), logical(1))]
  if (length(differs) == 0) { return(invisible(TRUE)) }
  labels <- c(backbone = "the backbone phylogeny", taxa = "the list of missing taxa", states = "the simmap", options = "other arguments (method, seed, number of simulations or maps, ...)")
  stop("Cannot resume from the temporary file '", path, "': it was produced with different inputs (", paste(labels[differs], collapse = ", "), ").\n",
       "Use the same inputs as in the interrupted run, or run with 'resume = FALSE' to overwrite the temporary files.", call. = FALSE)
}


#' @title State of the random number generator
#'
#' @noRd
#'
#' @return The `.Random.seed` of the session (`NULL` if none).
#'

get_rng_state <- function ()
{
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) { get(".Random.seed", envir = globalenv(), inherits = FALSE) } else { NULL }
}


#' @title Set the state of the random number generator
#'
#' @param state Output of `get_rng_state()`. Nothing is done if `NULL`.
#'
#' @noRd
#'
#' @return `NULL` (invisibly).
#'

set_rng_state <- function (state)
{
  if (!is.null(state)) { assign(".Random.seed", state, envir = globalenv()) }
  invisible(NULL)
}


#' @title Save the progress of a simulation in its temporary file
#'
#' @description Saves the phylogeny grafted so far (original node labels restored, as in the final output) with the state needed to resume.
#'
#' @param path Path of the temporary file of the simulation.
#' @param phy Phylogeny grafted so far, with the labels used to follow the nodes of the backbone.
#' @param labelled_backbone Output of `label_backbone_nodes()`.
#' @param signature Output of `temp_signature()`.
#' @param sim_seed Seed of the simulation (`NULL` if none: the state of the random number generator is then not saved).
#' @param done,total,unit Number of units grafted so far, total number of units, and name of the unit.
#' @param tracker,bd_cache,counters State followed during grafting (`NULL` if not used).
#'
#' @noRd
#'
#' @return `TRUE` (invisibly).
#'

save_temp_progress <- function (path, phy, labelled_backbone, signature, sim_seed, done, total, unit, tracker = NULL, bd_cache = NULL, counters = NULL)
{
  state <- list(version = 1L, status = "partial", signature = signature,
                progress = list(done = done, total = total, unit = unit),
                node_labels = phy$node.label, tracker = tracker, bd_cache = bd_cache, counters = counters,
                rng = if (is.null(sim_seed)) NULL else get_rng_state())
  saved <- restore_node_labels(phy = phy, labelled = labelled_backbone)
  attr(saved, "phylograftR_temp") <- state
  write_temp_file(saved, path)
}


#' @title Save the completed phylogeny of a simulation in its temporary file
#'
#' @description Replaces the partial file by the final output of the simulation (status `"complete"`).
#'
#' @param path Path of the temporary file of the simulation.
#' @param phy Final output of the simulation.
#' @param signature Output of `temp_signature()`.
#'
#' @noRd
#'
#' @return `TRUE` (invisibly).
#'

save_temp_complete <- function (path, phy, signature)
{
  attr(phy, "phylograftR_temp") <- list(version = 1L, status = "complete", signature = signature)
  write_temp_file(phy, path)
}


#' @title Read the temporary file of a simulation to resume it
#'
#' @param path Path of the temporary file of the simulation.
#' @param signature Signature of the inputs of the current call (see `temp_signature()`).
#'
#' @noRd
#'
#' @return `NULL` if there is no usable file (the simulation starts from scratch). Otherwise a list with `status` (`"complete"` or `"partial"`),
#'   `phy` (the final phylogeny, or the phylogeny grafted so far with the labels followed during grafting) and, for a partial file, `state`.
#'

read_resume_state <- function (path, signature)
{
  saved <- read_temp_file(path)
  if (is.null(saved)) { return(NULL) }
  state <- attr(saved, "phylograftR_temp")
  if (!is.list(state) || is.null(state$status) || is.null(state$signature))
  {
    warning("The temporary file '", path, "' does not hold the state needed to resume (it was not saved by this version of the function) and is ignored: this simulation starts from scratch.", call. = FALSE)
    return(NULL)
  }
  check_temp_signature(stored = state$signature, current = signature, path = path)
  attr(saved, "phylograftR_temp") <- NULL
  if (identical(state$status, "complete")) { return(list(status = "complete", phy = saved)) }
  saved$node.label <- state$node_labels
  list(status = "partial", phy = saved, state = state)
}


#' @title Paths of all the temporary files a run could write
#'
#' @description Rebuilds, from the arguments of the grafting function, the paths of the files of the simulations: one per map (if a `"multiSimmap"` is used) and per attempt
#'   (up to `nb_simulations`, or `10 * nb_simulations` with `force_simulations`, as in `run_attempts()`), following `temp_file_path_simulation()`.
#'
#' @param temp_file_path Absolute path template, see `resolve_temp_file_path()`.
#' @param map_names `NULL`, or names of the maps of the `"multiSimmap"`.
#' @param map_tags Named character vector of the tags of the maps, see `temp_map_tags()`.
#' @param nb_simulations,force_simulations Arguments of the grafting function.
#'
#' @return Character vector of paths (the `.part` files of the atomic writes are added). The files do not need to exist.
#'
#' @noRd
temp_files_of_run <- function (temp_file_path, map_names, map_tags, nb_simulations, force_simulations)
{
  multiple <- nb_simulations > 1
  max_attempts <- if (force_simulations) 10 * nb_simulations else nb_simulations
  maps <- if (is.null(map_names)) list(NULL) else as.list(map_names)
  paths <- unlist(lapply(maps, function (map_name)
  {
    vapply(seq_len(max_attempts), function (k)
    {
      temp_file_path_simulation(temp_file_path = temp_file_path, index = k, multiple = multiple, map_name = map_name,
                                map_tag = if (is.null(map_name)) NULL else map_tags[[map_name]])
    }, character(1))
  }))
  paths <- unique(paths)
  c(paths, paste0(paths, ".part"))
}


#' @title Delete the temporary files of a completed run
#'
#' @description Deletes the files returned by `temp_files_of_run()` that exist. Nothing else is touched.
#'   A warning lists the files that could not be deleted.
#'
#' @inheritParams temp_files_of_run
#'
#' @return The paths of the files deleted (invisibly).
#'
#' @noRd
remove_temp_files <- function (temp_file_path, map_names, map_tags, nb_simulations, force_simulations)
{
  paths <- temp_files_of_run(temp_file_path = temp_file_path, map_names = map_names, map_tags = map_tags,
                             nb_simulations = nb_simulations, force_simulations = force_simulations)
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) { return(invisible(character(0))) }
  removed <- suppressWarnings(file.remove(paths))
  if (any(!removed))
  {
    warning("Some temporary files could not be deleted: ", paste0("'", paths[!removed], "'", collapse = ", "), ".", call. = FALSE)
  }
  invisible(paths[removed])
}
