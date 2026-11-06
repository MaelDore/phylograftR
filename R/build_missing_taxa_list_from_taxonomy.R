##### Function to build a missing_taxa_list from taxonomic tables #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Build a missing_taxa_list from taxonomic tables ####

#' @title Build a `missing_taxa_list` from taxonomic tables
#'
#' @description Build a ready-to-be-used `missing_taxa_list` object based on taxonomic information from the backbone phylogeny and the taxa to graft.
#'   The `missing_taxa_list` can then be used as input for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'
#'   Target clades, sister clades and excluded subclades are identified with taxonomic names (e.g., `"Homo"`) or with taxonomic levels (e.g., `"Genus"`) instead of tips or node IDs.
#'   The function retrieves the MRCA node of each taxonomic group in the backbone phylogeny, checks that the group is monophyletic,
#'   and writes each clade as a pair of tips, as expected by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list].
#'
#' @param backbone_taxonomic_table Data frame. Taxonomy of the tips of the backbone phylogeny, with one row per tip.
#'   The first column must be named `Taxa` and holds the tip labels of `phy`.
#'   All the other columns are taxonomic levels (e.g., `Family`, `Subfamily`, `Tribe`, `Genus`, but they can be anything), usually from the broadest to the narrowest.
#'   The values of a level are the names of the taxonomic groups (e.g., `"Homo"`). Levels do not have to be informed for all the tips: blank cells and `NA` are missing values.
#'   See Details for the constraints on the names of the levels and groups.
#' @param missing_taxa_taxonomic_table Data frame. Description of the taxa to graft, with one row per taxon.
#'   Columns `Taxa` (names of the taxa to graft) and `Type` (`"Taxa"` or `"clades"`) are required, as well as the column of the clade to use for grafting (`Target_clade` or `Sister_clade`).
#'   Entries can be either the name of the taxonomic group (e.g., `"Homo"`) or the taxonomic level of the taxa to use (e.g., `Genus`).
#'   Next columns provide the taxonomic levels of the taxa to graft using the same names as in `backbone_taxonomic_table` (e.g., `Family`, `Subfamily`, `Tribe`, `Genus`, but they can be anything).
#'   Other optional columns provide additional information on constraints to use for grafting: `Excluded_subclades`, `include_target_clade_stem`, `exclude_subclades_stem`, and `allowed_states`,
#'   as defined in `missing_taxa_list` argument of the grafting functions ([`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades]).
#'   See 'Missing taxa taxonomic table' in Details section for an extended description.
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which the taxa will be grafted.
#' @param non_monophyletic Character string. What to do when a taxonomic group is not monophyletic in `phy`, i.e., when the clade defined by the MRCA of its tips
#'   also holds tips of other groups of the same level. Must be one of `"exclude"` (default), `"mrca"` or `"error"`. See 'Non-monophyletic groups' in Details.
#' @param verbose Logical. Whether to print a summary of the conversion. Default is `TRUE`.
#'
#' @details
#' **Backbone taxonomic table**
#'
#' Names of taxonomic levels (the column names; e.g., `"Genus"`) and names of taxonomic groups (their values; e.g., `"Homo"`) must be unambiguous:
#'
#' * A level name cannot be the name of a tip, of a taxon to graft, or of a group.
#' * A group name belongs to a single level: the same name cannot be found in two columns (e.g., a Family cannot be named the same as a Genus).
#' * Level names and group names cannot contain `";"`. Level names cannot be `"Taxa"`, `"Type"` or the names of the other columns of `missing_taxa_taxonomic_table`.
#' * All the tips of `backbone_taxonomic_table` must be in `phy`. Tips of `phy` absent from the table are allowed (with a warning).
#'
#' **Missing taxa taxonomic table**
#'
#' * `Taxa` (Required): Names of the taxa to graft. Names must be unique and absent from `phy$tip.label`.
#' * `Type` (Required): `"taxa"` for taxa grafted with [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] within a target clade,
#'   or `"clades"` for taxa grafted with [`phylograftR::graft_missing_clades()`][graft_missing_clades] as a sister of a clade.
#' * Taxonomic levels (Optional): Same columns as in `backbone_taxonomic_table`. They are used to retrieve the taxonomic group of the taxon when `Target_clade`, `Sister_clade`,
#'   and/or `Excluded_subclades` (see below) use a taxonomic level (e.g., `"Genus"`) rather than directly a taxonomic group (e.g., `"Homo"`) as identification.
#'   Columns that are not levels of the backbone, nor one of the columns listed here, are ignored with a warning.
#' * `Target_clade` (Required for `Type = "taxa"`) or `Sister_clade` (Required for `Type = "clades"`): A single taxonomic group, written either as the name of the group
#'   (e.g., `"Homo"`), or as the name of a taxonomic level (e.g., `"Genus"`). If a taxonomic level is provided, the group will be identified as recorded in the column of the same name
#'   in `missing_taxa_taxonomic_table` (e.g., `"Genus"` with `Genus = "Camponotus"` gives the clade of the genus `Camponotus`).
#'   It is an error if the group is absent from the backbone, or if the taxon has no value at the requested level.
#' * `Excluded_subclades` (Optional, `Type = "taxa"`): Subclades of the target clade where the taxon is known not to belong, separated by semicolons (e.g., `"Pheidole; Camponotus"`).
#'   Each one is a group or a level name, as in `Target_clade`.
#' * `$include_target_clade_stem` (Optional, `Type = "taxa"`): Logical. Whether the stem branch of the target clade is a valid placement to graft this taxa (`NA` uses the global argument set in the grafting function).
#' * `$exclude_subclades_stem` (Optional, `Type = "taxa"`): Logical. Whether the stem branch of the excluded clade(s) is/are a exlcuded from grafting too. Separated by semicolons if multiple `Excluded_subclades` are provided
#'   (`NA` uses the global argument set in the grafting function).
#' * `$allowed_states` (Optional): For each taxon, the states of a `simmap` in which it can be grafted
#'   (see the `simmap` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades]), separated by commas.
#'   It is accepted with both `Target_clade` and `Sister_clade`. (`NA` uses the global argument set in the grafting function).
#'
#' **Non-monophyletic groups**
#'
#' What to do when a taxonomic group is not monophyletic in `phy`, i.e., when the clade defined by the MRCA of its tips also holds tips of other groups of the same level?
#'
#' * `non_monophyletic = "exclude"` (default): For target clades, intruders (i.e., tips from other groups descending from the MRCA) are added to the
#'   `Excluded_subclades` of the taxon (with their stem, `exclude_subclades_stem = TRUE`), so the taxon is never grafted among intruders.
#'   Intruders nested within an excluded subclade provided by the user are not added.
#'   Intruders cannot be removed from a sister clade: the clade of the MRCA is used with a warning (as for `"mrca"`).
#' * `non_monophyletic = "mrca"`: The clade of the MRCA is used as such, with a warning, including its intruders.
#' * `non_monophyletic = "error"`: Stop if one group used is not monophyletic.
#'
#' @return A `missing_taxa_list` object as returned by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]:
#'   a list ready-to-be-used as input of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] if all the taxa have `Type = "Taxa"`,
#'   or of [`phylograftR::graft_missing_clades()`][graft_missing_clades] if they all have `Type = "clades"`.
#'   If both types are present, a list with two elements, `taxa` and `clades`, each holding the `missing_taxa_list` of one type.
#'
#'   The attribute `"taxonomic_groups"` is a data frame with one row per taxonomic group used, with the columns `Group`, `Level`, `n_members` (number of tips of the group),
#'   `clade_node` (MRCA node ID in `phy`), `clade_tips` (pair of tips defining the clade), `n_tips_clade` (number of tips of the clade), `n_foreign` (tips of other groups in the clade),
#'   `n_unannotated` (tips without value at this level in the clade), and `monophyletic`.
#'
#' @seealso Other functions from phylograftR:
#' * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]: to build a `missing_taxa_list` object from a data.frame loaded within R,
#' * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list]: to build a `missing_taxa_list` object from an external table,
#' * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa], [`phylograftR::graft_missing_clades()`][graft_missing_clades]: to use the list for grafting.
#'
#' @importFrom ape reorder.phylo
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' ## Backbone phylogeny with 10 tips
#' phy <- ape::read.tree(text = "(((((a1,a2),a3),a4),((b1,b2),b3)),((c1,c2),c3));")
#' phy <- ape::compute.brlen(phy)
#'
#' ## Taxonomy of the backbone
#' # Fam1 = GenA + GenB; Fam2 = GenC
#' backbone_taxonomic_table <- data.frame(
#'   Taxa = phy$tip.label,
#'   Family = ifelse(grepl("^c", phy$tip.label), "Fam2", "Fam1"),
#'   Genus = unname(c(a = "GenA", b = "GenB", c = "GenC")[substr(phy$tip.label, 1, 1)]))
#' backbone_taxonomic_table
#'
#' ## Taxa to graft
#' # 'a_new' is grafted in its genus, 'b_new' in its family except in the genus GenA,
#' # and 'd_new' is grafted as the sister of the genus GenC.
#' missing_taxa_taxonomic_table <- data.frame(
#'   Taxa = c("a_new", "b_new", "d_new"),
#'   Type = c("taxa", "taxa", "clades"),
#'   Family = c("Fam1", "Fam1", "Fam2"),
#'   Genus = c("GenA", "GenB", "GenD"),
#'   Target_clade = c("Genus", "Family", NA),
#'   Excluded_subclades = c(NA, "GenA", NA),
#'   Sister_clade = c(NA, NA, "GenC"))
#'
#' ## Build the missing taxa lists associated to these taxonomic tables
#' missing_taxa_lists <- build_missing_taxa_list_from_taxonomy(
#'   backbone_taxonomic_table = backbone_taxonomic_table,
#'   missing_taxa_taxonomic_table = missing_taxa_taxonomic_table,
#'   phy = phy)
#'
#' # Since both types of taxa ti graft where provided, the output contains one list per grafting function
#' names(missing_taxa_lists)
#' print(missing_taxa_lists$taxa) # Missing_taxa_list for graft_missing_taxa()
#' print(missing_taxa_lists$clades) # Missing_taxa_list for graft_missing_clades()
#'
#' # Information on the groups used
#' attr(missing_taxa_lists, "taxonomic_groups")
#'

build_missing_taxa_list_from_taxonomy <- function (backbone_taxonomic_table, # Data frame, one row per tip of the backbone
                                                   missing_taxa_taxonomic_table, # Data frame, one row per taxon to graft
                                                   phy, # Backbone phylogeny
                                                   non_monophyletic = c("exclude", "mrca", "error"),
                                                   verbose = TRUE)
{
  ## 1.1/ Check the inputs ####

  non_monophyletic <- match.arg(non_monophyletic)
  if (!inherits(phy, "phylo")) { stop("'phy' must be a phylogeny of class 'phylo'") }

  bb <- clean_taxonomic_table(backbone_taxonomic_table, "backbone_taxonomic_table")
  mt <- clean_taxonomic_table(missing_taxa_taxonomic_table, "missing_taxa_taxonomic_table")

  # Columns recognized in the table of the missing taxa
  fields <- c("Target_clade", "Sister_clade", "Excluded_subclades", "include_target_clade_stem", "exclude_subclades_stem", "allowed_states")
  reserved <- c("Taxa", "Type", fields)

  ## 1.2/ Backbone table ####

  if (!("Taxa" %in% names(bb))) { stop("'backbone_taxonomic_table' must have a column 'Taxa' with the tip labels of 'phy'") }
  levels_bb <- setdiff(names(bb), "Taxa")
  if (length(levels_bb) == 0) { stop("'backbone_taxonomic_table' must have at least one taxonomic level (a column other than 'Taxa')") }
  if (nrow(bb) == 0) { stop("'backbone_taxonomic_table' has no rows") }
  if (any(levels_bb %in% reserved)) { stop("Names of taxonomic levels cannot be one of: ", paste(reserved, collapse = ", "), ". Found: ", paste(intersect(levels_bb, reserved), collapse = ", ")) }
  if (any(grepl(";", levels_bb, fixed = TRUE))) { stop("Names of taxonomic levels cannot contain ';'. Found: ", paste(grep(";", levels_bb, fixed = TRUE, value = TRUE), collapse = ", ")) }
  if (anyNA(bb$Taxa)) { stop("'Taxa' is empty in the row(s) of 'backbone_taxonomic_table': ", paste(which(is.na(bb$Taxa)), collapse = ", ")) }
  if (anyDuplicated(bb$Taxa)) { stop("'Taxa' must be unique in 'backbone_taxonomic_table'. Duplicated: ", paste(unique(bb$Taxa[duplicated(bb$Taxa)]), collapse = ", ")) }
  not_in_phy <- setdiff(bb$Taxa, phy$tip.label)
  if (length(not_in_phy) > 0)
  {
    stop(length(not_in_phy), " taxa of 'backbone_taxonomic_table' are not tips of 'phy': ", format_names(not_in_phy))
  }
  not_in_table <- setdiff(phy$tip.label, bb$Taxa)
  if (length(not_in_table) > 0)
  {
    warning(length(not_in_table), " tips of 'phy' are absent from 'backbone_taxonomic_table': they are considered unannotated and are never excluded from a clade. ",
            format_names(not_in_table), call. = FALSE)
  }

  ## 1.3/ Table of the missing taxa ####

  for (column in c("Taxa", "Type")) { if (!(column %in% names(mt))) { stop("'missing_taxa_taxonomic_table' must have a column '", column, "'") } }
  if (nrow(mt) == 0) { stop("'missing_taxa_taxonomic_table' has no rows") }
  unknown <- setdiff(names(mt), c("Taxa", "Type", levels_bb, fields))
  if (length(unknown) > 0)
  {
    warning("Unknown columns ignored in 'missing_taxa_taxonomic_table': ", paste0("'", unknown, "'", collapse = ", "), ".\n",
            "Recognized columns: Taxa, Type, the taxonomic levels of the backbone (", paste(levels_bb, collapse = ", "), "), ", paste(fields, collapse = ", "), ".", call. = FALSE)
  }
  levels_mt <- intersect(levels_bb, names(mt))
  n_missing <- nrow(mt)
  row_id <- seq_len(n_missing)
  taxa_mt <- mt$Taxa
  if (anyNA(taxa_mt)) { stop("'Taxa' is empty in the row(s) of 'missing_taxa_taxonomic_table': ", paste(row_id[is.na(taxa_mt)], collapse = ", ")) }
  if (anyDuplicated(taxa_mt)) { stop("'Taxa' must be unique in 'missing_taxa_taxonomic_table'. Duplicated: ", paste(unique(taxa_mt[duplicated(taxa_mt)]), collapse = ", ")) }
  in_backbone <- taxa_mt %in% c(phy$tip.label, bb$Taxa)
  if (any(in_backbone))
  {
    stop("'Taxa' of 'missing_taxa_taxonomic_table' must be absent from the backbone (these are taxa to graft). Present in the backbone: ", format_names(taxa_mt[in_backbone]))
  }

  # Type of each taxon
  type <- tolower(mt$Type)
  bad_type <- is.na(type) | !(type %in% c("taxa", "clades"))
  if (any(bad_type))
  {
    stop("'Type' must be 'Taxa' or 'clades' (case is ignored). Invalid for: ",
         format_names(paste0(taxa_mt[bad_type], " (", ifelse(is.na(mt$Type[bad_type]), "empty", paste0("'", mt$Type[bad_type], "'")), ")")))
  }
  for (what in c("taxa", "clades"))
  {
    needed <- if (what == "taxa") "Target_clade" else "Sister_clade"
    if (any(type == what) && !(needed %in% names(mt))) { stop("'missing_taxa_taxonomic_table' needs a column '", needed, "' for the rows with Type = '", what, "'") }
  }

  # Columns not used by the type of the row: warn and ignore
  irrelevant <- list(taxa = c("Sister_clade"), clades = c("Target_clade", "Excluded_subclades", "include_target_clade_stem", "exclude_subclades_stem"))
  for (what in names(irrelevant))
  {
    for (column in intersect(irrelevant[[what]], names(mt)))
    {
      filled <- which(type == what & !is.na(mt[[column]]))
      if (length(filled) > 0)
      {
        warning("'", column, "' is not used with Type = '", what, "': ignored for ", format_names(taxa_mt[filled]), ".", call. = FALSE)
      }
    }
  }

  ## 1.4/ Names of levels and groups ####

  all_taxa_names <- c(bb$Taxa, phy$tip.label, taxa_mt)
  group_values <- rbind(
    do.call(rbind, lapply(levels_bb, function (l) data.frame(value = bb[[l]], level = l, stringsAsFactors = FALSE))),
    do.call(rbind, lapply(levels_mt, function (l) data.frame(value = mt[[l]], level = l, stringsAsFactors = FALSE))))
  group_values <- unique(group_values[!is.na(group_values$value), , drop = FALSE])
  problems <- character(0)
  if (any(levels_bb %in% all_taxa_names)) { problems <- c(problems, paste0("Names of levels that are also names of taxa: ", format_names(intersect(levels_bb, all_taxa_names)))) }
  if (any(levels_bb %in% group_values$value)) { problems <- c(problems, paste0("Names of levels that are also names of groups: ", format_names(intersect(levels_bb, group_values$value)))) }
  if (any(group_values$value %in% all_taxa_names)) { problems <- c(problems, paste0("Names of groups that are also names of taxa: ", format_names(intersect(group_values$value, all_taxa_names)))) }
  if (any(grepl(";", group_values$value, fixed = TRUE))) { problems <- c(problems, paste0("Names of groups containing ';': ", format_names(grep(";", group_values$value, fixed = TRUE, value = TRUE)))) }
  multi_level <- unique(group_values$value[duplicated(group_values$value)])
  if (length(multi_level) > 0)
  {
    detail <- vapply(multi_level, function (g) paste0(g, " (", paste(group_values$level[group_values$value == g], collapse = ", "), ")"), character(1))
    problems <- c(problems, paste0("Names of groups found in several levels (a group must belong to a single level): ", format_names(detail)))
  }
  if (length(problems) > 0)
  {
    stop("Ambiguous names in the taxonomic tables. A clade is written with the name of a group or of a level, so these names must be distinct from each other and from the names of taxa.\n",
         paste0("  - ", problems, collapse = "\n"))
  }

  ## 1.5/ Index the backbone ####

  index <- index_tree_ranges(phy)
  tip_id <- match(bb$Taxa, phy$tip.label)
  # For each level: the group of each tip, ordered by position in the tree (NA for unannotated tips)
  group_by_pos <- lapply(stats::setNames(levels_bb, levels_bb), function (l)
  {
    v <- rep(NA_character_, index$n_tips)
    v[index$pos[tip_id]] <- bb[[l]]
    v
  })

  ## 2/ Resolve the names written in the clades ####

  # A cell is the name of a group or of a level. A level means 'the group of the taxon at this level'.
  # Returns the pair (level, group), or the description of the problem
  problem_rows <- list()
  resolve_token <- function (token, i, column)
  {
    level <- NA_character_
    group <- NA_character_
    problem <- NA_character_
    if (token %in% levels_bb)
    {
      level <- token
      if (!(token %in% names(mt))) { problem <- paste0("'", token, "' is a level, but 'missing_taxa_taxonomic_table' has no column '", token, "'") }
      else if (is.na(mt[[token]][i])) { problem <- paste0("'", token, "' is a level, but the taxon has no value at this level") }
      else { group <- mt[[token]][i] }
    }
    else if (token %in% group_values$value)
    {
      group <- token
      level <- group_values$level[match(token, group_values$value)]
    }
    else
    {
      problem <- paste0("'", token, "' is neither a taxonomic level nor a group of the backbone")
    }
    if (is.na(problem) && !any(group_by_pos[[level]] == group, na.rm = TRUE))
    {
      problem <- paste0("group '", group, "' (", level, ") is absent from the backbone")
    }
    if (!is.na(problem))
    {
      problem_rows[[length(problem_rows) + 1]] <<- data.frame(Taxa = taxa_mt[i], Column = column, Problem = problem, stringsAsFactors = FALSE)
      return(NULL)
    }
    list(level = level, group = group, key = paste(level, group, sep = "\r"))
  }

  rows_taxa <- which(type == "taxa")
  rows_clades <- which(type == "clades")
  resolve_single <- function (rows, column)
  {
    res <- vector("list", n_missing)
    for (i in rows)
    {
      if (is.na(mt[[column]][i]))
      {
        problem_rows[[length(problem_rows) + 1]] <<- data.frame(Taxa = taxa_mt[i], Column = column, Problem = "empty", stringsAsFactors = FALSE)
        next
      }
      res[[i]] <- resolve_token(mt[[column]][i], i, column)
    }
    res
  }
  target_res <- if (length(rows_taxa) > 0) resolve_single(rows_taxa, "Target_clade") else vector("list", n_missing)
  sister_res <- if (length(rows_clades) > 0) resolve_single(rows_clades, "Sister_clade") else vector("list", n_missing)

  excluded_res <- vector("list", n_missing)
  stem_values <- vector("list", n_missing)
  stem_ignored <- character(0)
  stem_wrong <- character(0)
  if ("Excluded_subclades" %in% names(mt))
  {
    for (i in rows_taxa)
    {
      cell <- mt$Excluded_subclades[i]
      n_user <- 0
      if (!is.na(cell))
      {
        tokens <- trimws(strsplit(cell, ";", fixed = TRUE)[[1]])
        if (length(tokens) == 0 || any(tokens == ""))
        {
          problem_rows[[length(problem_rows) + 1]] <- data.frame(Taxa = taxa_mt[i], Column = "Excluded_subclades",
                                                                 Problem = paste0("empty group (check the ';' separators): ", cell), stringsAsFactors = FALSE)
          next
        }
        resolved <- lapply(tokens, resolve_token, i = i, column = "Excluded_subclades")
        if (any(vapply(resolved, is.null, logical(1)))) { next }
        excluded_res[[i]] <- resolved
        n_user <- length(tokens)
      }
      # Stem of the excluded subclades: one value per subclade or a single one
      if ("exclude_subclades_stem" %in% names(mt) && !is.na(mt$exclude_subclades_stem[i]))
      {
        values <- trimws(strsplit(mt$exclude_subclades_stem[i], ";", fixed = TRUE)[[1]])
        if (n_user == 0) { stem_ignored <- c(stem_ignored, taxa_mt[i]) }
        else if (!(length(values) %in% c(1, n_user))) { stem_wrong <- c(stem_wrong, paste0(taxa_mt[i], " (", length(values), " value(s) for ", n_user, " groups)")) }
        else { stem_values[[i]] <- if (length(values) == 1) rep(values, n_user) else values }
      }
    }
  }
  else if ("exclude_subclades_stem" %in% names(mt))
  {
    stem_ignored <- taxa_mt[rows_taxa[!is.na(mt$exclude_subclades_stem[rows_taxa])]]
  }
  if (length(stem_ignored) > 0) { warning("'exclude_subclades_stem' is provided for taxa without Excluded_subclades: ignored for ", format_names(stem_ignored), ".", call. = FALSE) }
  if (length(stem_wrong) > 0)
  {
    stop("'exclude_subclades_stem' needs one value per group of 'Excluded_subclades' (separated by ';'), or a single value for all of them. Problem for: ", format_names(stem_wrong))
  }

  if (length(problem_rows) > 0)
  {
    problems <- do.call(rbind, problem_rows)
    n_show <- min(nrow(problems), 15)
    stop(nrow(problems), " problem(s) found when retrieving the taxonomic groups of the clades:\n",
         paste0("  - ", problems$Taxa[seq_len(n_show)], " (", problems$Column[seq_len(n_show)], "): ", problems$Problem[seq_len(n_show)], collapse = "\n"),
         if (nrow(problems) > n_show) paste0("\n  ... and ", nrow(problems) - n_show, " more") else "",
         "\nA clade is written with the name of a group of 'backbone_taxonomic_table', or with the name of a level (the taxon must then have a value at this level).")
  }

  ## 3/ Clades of the taxonomic groups ####

  used <- c(target_res, sister_res, unlist(excluded_res, recursive = FALSE))
  used <- used[!vapply(used, is.null, logical(1))]
  used_keys <- unique(vapply(used, function (x) x$key, character(1)))
  groups_info <- lapply(stats::setNames(used_keys, used_keys), function (key)
  {
    parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
    describe_group(index, group_by_pos[[parts[1]]], level = parts[1], group = parts[2])
  })

  # Role of each group: target clade (target), sister clade (sister), excluded subclade (excluded)
  keys_of <- function (res) unique(vapply(res[!vapply(res, is.null, logical(1))], function (x) x$key, character(1)))
  role_keys <- list(target = keys_of(target_res), sister = keys_of(sister_res), excluded = keys_of(unlist(excluded_res, recursive = FALSE)))
  nonmono_role <- lapply(role_keys, function (keys) keys[!vapply(keys, function (k) groups_info[[k]]$monophyletic, logical(1))])

  # Sister clade with the whole backbone: no stem
  whole_tree <- role_keys$sister[vapply(role_keys$sister, function (k) groups_info[[k]]$node == index$root, logical(1))]
  if (length(whole_tree) > 0)
  {
    stop("The clade of a sister group holds the whole backbone, so there is no stem branch to graft on: ",
         format_names(vapply(whole_tree, function (k) paste0(groups_info[[k]]$group, " (", groups_info[[k]]$level, ")"), character(1))))
  }

  describe_nonmono <- function (keys)
  {
    paste0("  - ", vapply(keys, function (k) { g <- groups_info[[k]]; paste0(g$group, " (", g$level, "): ", g$n_foreign, " tip(s) of other groups among ", g$n_tips_clade, " tips") }, character(1)), collapse = "\n")
  }
  n_nonmono <- vapply(nonmono_role, length, integer(1))
  if (any(n_nonmono > 0) && non_monophyletic == "error")
  {
    stop("Some taxonomic groups used are not monophyletic in 'phy' (the clade of their MRCA holds tips of other groups):\n",
         paste0(vapply(names(nonmono_role)[n_nonmono > 0], function (r) paste0(c(target = "Target clades", sister = "Sister clades", excluded = "Excluded subclades")[r], ":\n", describe_nonmono(nonmono_role[[r]])), character(1)), collapse = "\n"),
         "\nUse 'non_monophyletic' = \"exclude\" or \"mrca\" to build the list anyway.")
  }
  if (n_nonmono["target"] > 0)
  {
    warning(if (non_monophyletic == "exclude")
              "Target groups are not monophyletic: the largest subclades of other groups are added to Excluded_subclades (stem included) to keep them out of the target clade:\n"
            else "Target groups are not monophyletic: their clade (MRCA of their tips) is used and includes tips of other groups:\n",
            describe_nonmono(nonmono_role$target), call. = FALSE)
  }
  if (n_nonmono["sister"] > 0)
  {
    warning("Sister groups are not monophyletic: their clade (MRCA of their tips) is used and includes tips of other groups (they cannot be excluded from a sister clade):\n",
            describe_nonmono(nonmono_role$sister), call. = FALSE)
  }
  if (n_nonmono["excluded"] > 0)
  {
    warning(if (non_monophyletic == "exclude")
              "Excluded groups are not monophyletic: they are replaced by their largest subclades with only tips of the group, so tips of other groups are not excluded:\n"
            else "Excluded groups are not monophyletic: their clade (MRCA of their tips) is used, so tips of other groups are excluded too:\n",
            describe_nonmono(nonmono_role$excluded), call. = FALSE)
  }

  # Nodes excluded for each group of Excluded_subclades (several nodes for non-monophyletic groups in 'exclude' mode), with their labels
  excluded_nodes <- lapply(stats::setNames(role_keys$excluded, role_keys$excluded), function (k)
  {
    g <- groups_info[[k]]
    if (non_monophyletic == "exclude" && !g$monophyletic)
    {
      is_member <- !is.na(group_by_pos[[g$level]]) & group_by_pos[[g$level]] == g$group
      nodes <- maximal_subclades(index, g$node, is_member)
      nodes[order(index$lo[nodes])]
    } else g$node
  })

  # Intruders of the target groups in 'exclude' mode
  intruders <- list()
  if (non_monophyletic == "exclude")
  {
    intruders <- lapply(stats::setNames(nonmono_role$target, nonmono_role$target), function (k)
    {
      g <- groups_info[[k]]
      is_foreign <- !is.na(group_by_pos[[g$level]]) & group_by_pos[[g$level]] != g$group
      nodes <- maximal_subclades(index, g$node, is_foreign)
      nodes[order(index$lo[nodes])]
    })
  }

  ## 4/ Names of the clades ####

  # A node can be the clade of several groups (e.g., a genus with a single species, and its family): all the names are kept, joined by " / ",
  # so the same node always has the same name, and different nodes have different names.
  {
    contrib_node <- numeric(0)
    contrib_name <- character(0)
    for (k in c(role_keys$target, role_keys$sister))
    {
      contrib_node <- c(contrib_node, groups_info[[k]]$node)
      contrib_name <- c(contrib_name, groups_info[[k]]$group)
    }
    for (k in role_keys$excluded)
    {
      nodes <- excluded_nodes[[k]]
      contrib_node <- c(contrib_node, nodes)
      contrib_name <- c(contrib_name, if (length(nodes) > 1) paste0(groups_info[[k]]$group, " [", seq_along(nodes), "]") else groups_info[[k]]$group)
    }
    node_names <- vapply(split(contrib_name, contrib_node), function (x) paste(unique(x), collapse = " / "), character(1))
    # Several groups for the same clade: report the names that were joined
    shared <- split(contrib_name, contrib_node)
    shared <- shared[vapply(shared, function (x) length(unique(x)) > 1, logical(1))]
    if (length(shared) > 0)
    {
      warning("Several groups have the same clade in 'phy' (e.g., a genus with a single species and its family). The clade has the same name in all the taxa, ",
              "made of the names of all these groups joined by ' / ', so some names are not the group written for the taxon:\n",
              paste0("  - ", vapply(shared, function (x) paste0(paste(unique(x), collapse = " / "), " (clade of ", paste(unique(x), collapse = ", "), ")"), character(1)), collapse = "\n"),
              "\nUse the same names in 'total_richness_per_clades'.", call. = FALSE)
    }
  }
  name_of_nodes <- function (nodes)
  {
    out <- unname(node_names[as.character(nodes)])
    unnamed <- is.na(out)
    if (any(unnamed)) { out[unnamed] <- paste0("foreign(", node_tips_text(index, nodes[unnamed], phy$tip.label), ")") }
    out
  }

  ## 5/ Tables for build_missing_taxa_list() ####

  # Tips of the clades
  tips_text <- function (nodes) node_tips_text(index, nodes, phy$tip.label)
  bad_tips <- phy$tip.label[grepl("[,;]", phy$tip.label)]
  if (length(bad_tips) > 0) { stop("Tip labels cannot contain ',' or ';' to be written in a clade. Found: ", format_names(bad_tips)) }

  build_table <- function (rows, what)
  {
    clade_item <- if (what == "taxa") "Target_clade" else "Sister_clade"
    res <- if (what == "taxa") target_res else sister_res
    group_info_of <- function (i) groups_info[[res[[i]]$key]]
    tbl <- data.frame(Taxon_name = taxa_mt[rows], stringsAsFactors = FALSE)
    clade_nodes <- vapply(rows, function (i) group_info_of(i)$node, numeric(1))
    tbl[[clade_item]] <- tips_text(clade_nodes)
    tbl[[paste0(clade_item, "_name")]] <- name_of_nodes(clade_nodes)
    if (what == "taxa")
    {
      if ("include_target_clade_stem" %in% names(mt)) { tbl$include_target_clade_stem <- mt$include_target_clade_stem[rows] }

      # Excluded subclades: the ones of the user, then the intruders of the target group
      excl_tips <- excl_names <- excl_stem <- rep(NA_character_, length(rows))
      for (r in seq_along(rows))
      {
        i <- rows[r]
        user_nodes <- numeric(0)
        user_stem <- character(0)
        for (j in seq_along(excluded_res[[i]]))
        {
          nodes <- excluded_nodes[[excluded_res[[i]][[j]]$key]]
          user_nodes <- c(user_nodes, nodes)
          user_stem <- c(user_stem, rep(if (is.null(stem_values[[i]])) "NA" else stem_values[[i]][j], length(nodes)))
        }
        intr_nodes <- if (res[[i]]$key %in% names(intruders)) intruders[[res[[i]]$key]] else numeric(0)

        # Keep the largest clades: remove user subclades nested within another one, then the intruders nested within a user subclade, and vice versa
        if (length(user_nodes) > 1)
        {
          lo <- index$lo[user_nodes]; hi <- index$hi[user_nodes]
          inside <- outer(lo, lo, ">=") & outer(hi, hi, "<=") # [a, b]: a within b
          diag(inside) <- FALSE
          same <- inside & t(inside)
          drop <- apply(inside & !same, 1, any) | apply(same & lower.tri(same), 1, any)
          user_nodes <- user_nodes[!drop]; user_stem <- user_stem[!drop]
        }
        if (length(intr_nodes) > 0 && length(user_nodes) > 0)
        {
          within_user <- vapply(seq_along(intr_nodes), function (x) any(index$lo[intr_nodes[x]] >= index$lo[user_nodes] & index$hi[intr_nodes[x]] <= index$hi[user_nodes]), logical(1))
          intr_nodes <- intr_nodes[!within_user]
          if (length(intr_nodes) > 0)
          {
            contains_intruder <- vapply(seq_along(user_nodes), function (x) any(index$lo[user_nodes[x]] >= index$lo[intr_nodes] & index$hi[user_nodes[x]] <= index$hi[intr_nodes]), logical(1))
            user_nodes <- user_nodes[!contains_intruder]; user_stem <- user_stem[!contains_intruder]
          }
        }
        nodes <- c(user_nodes, intr_nodes)
        if (length(nodes) == 0) { next }
        stems <- c(user_stem, rep("TRUE", length(intr_nodes)))
        excl_tips[r] <- paste(tips_text(nodes), collapse = ";")
        excl_names[r] <- paste(name_of_nodes(nodes), collapse = ";")
        if (!all(stems == "NA")) { excl_stem[r] <- paste(stems, collapse = ";") }
      }
      if (any(!is.na(excl_tips)))
      {
        tbl$Excluded_subclades <- excl_tips
        tbl$Excluded_subclades_names <- excl_names
        if (any(!is.na(excl_stem))) { tbl$exclude_subclades_stem <- excl_stem }
      }
    }
    if ("allowed_states" %in% names(mt)) { tbl$allowed_states <- mt$allowed_states[rows] }
    tbl
  }

  ## 6/ Build the lists ####

  build_list <- function (what, rows)
  {
    tbl <- build_table(rows, what)
    tryCatch(build_missing_taxa_list(tbl, phy = phy, verbose = FALSE),
             error = function (e) stop("The clades retrieved from the taxonomy are not valid:\n", conditionMessage(e), call. = FALSE))
  }
  lists <- list()
  if (length(rows_taxa) > 0) { lists$taxa <- build_list("taxa", rows_taxa) }
  if (length(rows_clades) > 0) { lists$clades <- build_list("clades", rows_clades) }
  result <- if (length(lists) == 1) lists[[1]] else lists

  # Information on the groups
  info_df <- do.call(rbind, lapply(groups_info, function (g)
  {
    data.frame(Group = g$group, Level = g$level, n_members = g$n_members, clade_node = g$node, clade_tips = tips_text(g$node), n_tips_clade = g$n_tips_clade,
               n_foreign = g$n_foreign, n_unannotated = g$n_unannotated, monophyletic = g$monophyletic, stringsAsFactors = FALSE)
  }))
  if (is.null(info_df))
  {
    info_df <- data.frame(Group = character(0), Level = character(0), n_members = integer(0), clade_node = numeric(0), clade_tips = character(0),
                          n_tips_clade = integer(0), n_foreign = integer(0), n_unannotated = integer(0), monophyletic = logical(0))
  }
  rownames(info_df) <- NULL
  attr(result, "taxonomic_groups") <- info_df

  if (verbose)
  {
    message(length(rows_taxa), " taxa for graft_missing_taxa() and ", length(rows_clades), " for graft_missing_clades() built from ", nrow(info_df), " taxonomic group(s), ",
            sum(!info_df$monophyletic), " not monophyletic. See attr(, \"taxonomic_groups\").")
  }
  return(result)
}


### 2/ Internal helpers ####

#' @title Clean a taxonomic table
#'
#' @description Check that the input is a data frame, trim the names of the columns, and convert every cell to trimmed text, with blank cells and `"NA"` as `NA`.
#'
#' @param x Data frame.
#' @param arg Character string. Name of the argument, for the messages.
#'
#' @noRd
#'
#' @return Data frame of character columns.
#'

clean_taxonomic_table <- function (x, arg)
{
  if (!is.data.frame(x)) { stop("'", arg, "' must be a data frame") }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  names(x) <- trimws(names(x))
  if (anyDuplicated(names(x))) { stop("Some columns of '", arg, "' have the same name: ", paste(unique(names(x)[duplicated(names(x))]), collapse = ", ")) }
  x[] <- lapply(x, function (v) { v <- trimws(as.character(v)); v[!is.na(v) & v %in% c("", "NA")] <- NA; v })
  x
}


#' @title Format a vector of names for messages
#'
#' @param x Character vector.
#' @param n Integer. Maximum number of names shown.
#'
#' @noRd
#'
#' @return Character string.
#'

format_names <- function (x, n = 10)
{
  x <- unique(x)
  paste0(paste(utils::head(x, n), collapse = ", "), if (length(x) > n) paste0(", ... and ", length(x) - n, " more") else "")
}


#' @title Index the tips of each clade of a phylogeny
#'
#' @description Order the tips so that each clade is a range of positions. Each node gets the first (`lo`) and last (`hi`) positions of its tips.
#'   Then the MRCA of a set of tips is the smallest node whose range covers their positions, and the tips of a clade are the tips in this range.
#'
#' @param phy Phylogeny of class `"phylo"`.
#'
#' @noRd
#'
#' @return A list with `n_tips`, `root`, `pos` (position of each tip ID), `tip_order` (tip ID of each position), `parent` (parent of each node, 0 for the root),
#'   and `lo`, `hi` (first and last position of the tips of each node).
#'

index_tree_ranges <- function (phy)
{
  n_tips <- length(phy$tip.label)
  n_nodes <- n_tips + phy$Nnode
  phy <- ape::reorder.phylo(phy, "cladewise")
  edge <- phy$edge
  tip_order <- edge[edge[, 2] <= n_tips, 2]
  pos <- integer(n_tips)
  pos[tip_order] <- seq_len(n_tips)
  parent <- integer(n_nodes)
  parent[edge[, 2]] <- edge[, 1]
  root <- setdiff(edge[, 1], edge[, 2])
  lo <- rep(Inf, n_nodes)
  hi <- rep(-Inf, n_nodes)
  lo[seq_len(n_tips)] <- pos
  hi[seq_len(n_tips)] <- pos
  # Descendants come before their ancestors in the reversed cladewise order
  for (e in rev(seq_len(nrow(edge))))
  {
    p <- edge[e, 1]
    child <- edge[e, 2]
    if (lo[child] < lo[p]) { lo[p] <- lo[child] }
    if (hi[child] > hi[p]) { hi[p] <- hi[child] }
  }
  list(n_tips = n_tips, root = root, pos = pos, tip_order = tip_order, parent = parent, lo = as.integer(lo), hi = as.integer(hi))
}


#' @title Pair of tips defining each clade
#'
#' @param index Output of `index_tree_ranges()`.
#' @param nodes Integer vector of node IDs.
#' @param tip_labels Character vector. Tip labels of the phylogeny.
#'
#' @noRd
#'
#' @return Character vector with one text per node: the tip itself for a tip, and the first and last tips of the clade separated by a comma otherwise.
#'

node_tips_text <- function (index, nodes, tip_labels)
{
  if (length(nodes) == 0) { return(character(0)) }
  first <- tip_labels[index$tip_order[index$lo[nodes]]]
  last <- tip_labels[index$tip_order[index$hi[nodes]]]
  ifelse(index$lo[nodes] == index$hi[nodes], first, paste(first, last, sep = ","))
}


#' @title Clade of a taxonomic group
#'
#' @param index Output of `index_tree_ranges()`.
#' @param group_pos Character vector. Group of each tip at the level, by position in the tree (`NA` if unannotated).
#' @param level,group Character strings.
#'
#' @noRd
#'
#' @return A list with the MRCA node of the group and counts of the tips of the clade.
#'

describe_group <- function (index, group_pos, level, group)
{
  member_pos <- which(group_pos == group)
  first <- min(member_pos)
  last <- max(member_pos)
  candidates <- which(index$lo <= first & index$hi >= last)
  node <- candidates[which.min(index$hi[candidates] - index$lo[candidates])]
  in_clade <- group_pos[index$lo[node]:index$hi[node]]
  n_foreign <- sum(!is.na(in_clade) & in_clade != group)
  list(level = level, group = group, node = node, n_members = length(member_pos), n_tips_clade = length(in_clade),
       n_foreign = n_foreign, n_unannotated = sum(is.na(in_clade)), monophyletic = n_foreign == 0)
}


#' @title Largest subclades with only tips that satisfy a condition
#'
#' @param index Output of `index_tree_ranges()`.
#' @param node Integer. Node of the clade where to search.
#' @param good Logical vector. Whether each tip (by position in the tree) satisfies the condition.
#'
#' @noRd
#'
#' @return Integer vector of nodes of the subclades whose tips all satisfy the condition, and whose parent does not.
#'

maximal_subclades <- function (index, node, good)
{
  cs <- c(0, cumsum(good))
  within <- which(index$lo >= index$lo[node] & index$hi <= index$hi[node])
  all_good <- (cs[index$hi + 1] - cs[index$lo]) == (index$hi - index$lo + 1)
  keep <- all_good[within] & (within == node | !all_good[pmax(index$parent[within], 1)])
  within[keep]
}
