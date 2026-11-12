##### Functions to build total_richness_per_clades from an internal data.frame #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Build total_richness_per_clades from an internal data.frame ####

#' @title Build the `total_richness_per_clades` object from a table of clades
#'
#' @description Build the `total_richness_per_clades` object used by [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#'   and [`phylograftR::graft_missing_clades()`][graft_missing_clades] to account for the species that are neither in the backbone nor in the `missing_taxa_list` (i.e., unsampled species),
#'   from **a single table** with one row per clade (e.g., read from a .csv file) and the total number of species of the clade.
#'
#'   The clades can be written with tips, node IDs, or **names** (see Details), and are matched with the clades of the `missing_taxa_list`
#'   object built to be used in the same call to a grafting function: this ensures that **the tips and the names of a clade are the same in all the objects**.
#'   The function checks that the richness is **consistent with the full sequence of grafting**, before any grafting.
#'
#'   The output follows the structure of `missing_taxa_list`:
#'   if it contains both types of taxa (as returned by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]
#'   and [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy]), a list for each grafting function is returned.
#'
#'   Use [`phylograftR::import_total_richness_per_clades()`][import_total_richness_per_clades] to read the table directly from a file.
#'
#' @param data A data.frame, with one row per clade.
#'   Columns are:
#'   * `Clade`: a monophyletic group from the backbone phylogeny defined by tips, a node ID, or a name (see Details), from which the total richness is known.
#'   * `total_richness`: Total number of species of the clade. It includes the backbone species, **all the taxa to graft that belong to the clade**, and the unsampled species that are known to exist but will not be grafted.
#'     Use `NA` if unknown.
#'   * `Clade_name` (Optional): Name of the clade, used to check matching with other inputs, and to report warnings in downstream functions. If the clade is written with a name in `Clade`, it is its name.
#'
#'   The columns `Target_clade` and `Target_clade_name`, or `Sister_clade` and `Sister_clade_name`, are accepted instead of `Clade` and `Clade_name` (e.g., for tables prepared for a single grafting function).
#'   All the columns are read as text. Blank cells and `NA` are missing values.
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which the taxa will be grafted. It is used to check the tips and clades, and to convert node IDs into tips.
#' @param missing_taxa_list A `missing_taxa_list`, or the list of two `missing_taxa_list` (`$taxa` and `$clades`) returned by
#'   [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list], [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list]
#'   or [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy].
#'   It provides the tips and the names used for the clades when they are only named, and is used to check the validity of the richness values provided.
#' @param return_richness_table Logical. Whether to return the expected richness and sampling fraction of each clade in the attribute `"richness_table"` of the output.
#'   It is computed before any grafting, from the backbone, the `missing_taxa_list` and the total richness. It is a tool to **examine a priori the sampling fraction
#'   expected in each clade**. The actual **effective sampling fractions of the grafted phylogenies**, accessible with [`phylograftR::extract_sampling_info()`][extract_sampling_info],
#'   may differ (see 'Richness table' in Details). Default is `TRUE`.
#' @param verbose Logical. Whether to print the conversions of node IDs into tips (if any), and a summary of the table. Default is `TRUE`.
#'
#' @details
#' ### Identifying clades
#'
#' Clades can be defined by tips, a node ID, or a name:
#'  * A pair of tips defines the clade as their MRCA node (e.g., `"t26, t18"`), separated by commas. A single tip defines the "clade" as the terminal branch of this tip.
#'    **Tip labels must include at least one other character than digits.**
#'  *  A node ID of `phy` (digits only, e.g., `"34"`) also defines the clade, as in [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list].
#'  *  A single value that is not a tip of `phy` is a **clade name** (e.g., `"Family_A"`). Names are matched with the names of the clades of `missing_taxa_list`
#'    (`$Target_clade_name`, `$Sister_clade_name`, and the names of the excluded subclades) and, for lists built from taxonomic tables, with all the taxonomic groups of the backbone taxonomy.
#'    The taxonomic groups are stored in the attribute `"taxonomic_groups"` of the object returned by [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy]:
#'    they are only available if `missing_taxa_list` is this object itself, not one of its elements (`$taxa` or `$clades`) alone.
#'
#' **Clades are matched by their MRCA node in `phy`**: a clade can only be listed once, and
#' the tips written in the output are those already used in `missing_taxa_list` for the same clade, so the tips match in all the objects.
#' The grafting functions require **all the clades used in `missing_taxa_list` to be listed**:
#' the clades that are not in the table are added with an unknown richness (`NA`), and their number is reported.
#'
#' ### Unsampled species
#'
#' For a focal clade with a `total_richness` value provided: `total_richness = nb_backbone + nb_missing + nb_unsampled`,
#' with `nb_missing` the number of lineages to graft within the focal clade, and `nb_unsampled` the number of species in this clade that are not represented in the backbone, nor will be grafted,
#' but must be accounted in the sampling fraction for when fitting diversification models.
#'
#' For a clade with no `total_richness` value (`NA` or `NULL`), `nb_unsampled` are estimated as the sum of `nb_unsampled` over the largest clades nested within it with a `total_richness` value provided.
#' If there is none, `nb_unsampled = 0`. The rationale is that even if `nb_unsampled` cannot be derived from its `total_richness`
#' (as it is not provided for this clade), it has to at least account for the `nb_unsampled` species that have being identified in its own nested clades.
#'
#' ### Full sequence of grafting
#'
#' The same table is valid for the grafting of the taxa within target clades and for the grafting of the missing clades as sister to backbone clades, in either order.
#' With the recommended sequence (taxa first, then missing clades on the grafted phylogeny):
#'  * When the taxa are grafted, the taxa to graft as missing clades are not in the backbone nor in the list, so they are counted among the unsampled species of the clades.
#'  * When the missing clades are grafted, the taxa grafted before are tips of the backbone, and the unsampled species are the species still unaccounted for.
#'
#' The function checks the richness against all the taxa to graft at once: **the function stops if the total richness of a clade is lower than the number of species it contains**
#' (backbone species and all the taxa to graft within it, or within clades nested within it),
#' or if a clade has fewer unsampled species than the clades with a total richness nested within it.
#'
#' ### Richness table
#'
#' With `return_richness_table = TRUE`, the output has the attribute `"richness_table"`: a data.frame that reports, **before any grafting**,
#' the richness and the sampling fraction expected in each clade (see Value for details on its columns). It can be used to **examine a priori the sampling fraction of the clades**
#' and to spot errors in the total richness: e.g., an unlikely low sampling fraction, or clades with an unknown richness (`NA`) for which unsampled species were derived from nested clades, or assumed to be none.
#' It can also be used to compare the sampling of different clades, or to decide which clades need a better estimate of their richness.
#'
#' **The sampling fractions of the table are expectations, and the effective sampling fractions of the grafted phylogenies may differ.**
#' The table assumes that all the taxa to graft end up in the clade where they belong: a taxon counts for all the clades that are, or enclose, its target clade,
#' and a missing clade counts for the clades that enclose its sister clade (it is grafted on the stem of the sister clade, so not within the sister clade itself).
#' However, where taxa are actually grafted varies among the grafted phylogenies (and among the stochastic maps, if any):
#' a taxon can be grafted on the stem of its target clade, which is outside of the clade, or within a clade nested in its target clade, which changes the richness of these clades.
#' The actual effective richness and sampling fraction of each clade, in each grafted phylogeny, are extracted **a posteriori** with [`phylograftR::extract_sampling_info()`][extract_sampling_info]
#' on the output of the grafting functions, using the same `total_richness_per_clades`.
#'
#' @return If `data` contains one type of taxa, a `total_richness_per_clades` list;
#'   if it contains both types (`$taxa` and `$clades`), a list with two `total_richness_per_clades` named `$taxa` (for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]) and `$clades`
#'   (for [`phylograftR::graft_missing_clades()`][graft_missing_clades]).
#'
#'   Each `total_richness_per_clades` list contains all the clades of the table, and has the following items:
#'  * `$Target_clade` (or `$Sister_clade`, for the list of missing clades) List of character vectors of tips defining the clades.
#'  * `$Target_clade_name` (or `$Sister_clade_name`) Only if some clades have a name. Character vector (`NA` if a clade has no name).
#'  * `$total_richness` Integer vector of the total richness of the clades (`NA` if unknown).
#'
#'  If `return_richness_table = TRUE`, the output has the attribute `"richness_table"`, a data.frame with one row per clade, and the columns:
#'  * `Clade_name`: Name of the clade (`NA` if none).
#'  * `Clade_tips`: Tips defining the clade, as in the output.
#'  * `total_richness_provided`: Total richness of the table (`NA` if unknown).
#'  * `total_richness`: Total richness as used for the grafting, i.e., the value provided, or, if unknown, the species known in the clade (backbone and taxa to graft) plus the unsampled species derived from the clades nested within it (none if there are no such clades).
#'  * `nb_backbone_species`: Number of tips of the backbone in the clade.
#'  * `nb_missing_species`: Number of taxa to graft expected in the clade, of both types (see 'Richness table' in Details).
#'  * `nb_unsampled_species`: Number of unsampled species, i.e., species that are neither in the backbone nor to graft (`total_richness` - `nb_backbone_species` - `nb_missing_species`).
#'  * `sampling_fraction`: Expected sampling fraction of the clade in the grafted phylogenies, `(nb_backbone_species + nb_missing_species) / total_richness`. It is computed **a priori**,
#'    and the effective sampling fraction may differ after grafting: use [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract it **a posteriori**.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::import_total_richness_per_clades()`][import_total_richness_per_clades] to read the table directly from an external file.
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] and [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to prepare the `missing_taxa_list`.
#'  * [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy] to prepare the `missing_taxa_list` from taxonomic tables.
#'  * [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] to convert the node IDs of a `total_richness_per_clades` written 'by hand'.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft the taxa using the richness.
#'  * [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract the richness and sampling fraction of the clades after grafting.
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Load the backbone, the taxonomic tables, and the toy lists of taxa and missing clades
#' data(toy_backbone)
#' data(toy_backbone_taxonomic_table)
#' data(toy_missing_taxa_taxonomic_table)
#'
#' # Build the missing taxa list for both types of missing lineages (taxa and clades)
#' lists <- build_missing_taxa_list_from_taxonomy(
#'   backbone_taxonomic_table = toy_backbone_taxonomic_table,
#'   missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
#'   phy = toy_backbone, verbose = FALSE)
#' names(lists)
#'
#' # Inform the total richness of some clades, written with their taxonomic names.
#' # Use NA if the total richness is unknown
#' richness_table <- data.frame(Clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"),
#'                              total_richness = c(25, NA, 14, 8))
#'
#' # Use the the function to validate this richness table
#' # against the backbone tree and missing taxon lists
#' richness_per_clades <- build_total_richness_per_clades(data = richness_table, phy = toy_backbone,
#'                                                        missing_taxa_list = lists)
#'
#' # Output = one list for each grafting function,
#' # with the clades of both lists and the same tips as in 'lists'
#' str(richness_per_clades, 2)
#'
#' # Expected richness and sampling fraction of each clade, examined before any grafting
#' attr(richness_per_clades, "richness_table")
#'
#' # Use the lists to build the fully grafted phylogenies in the recommended sequence:
#' # the missing taxa first, and then the missing clades
#' grafted_taxa <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = lists$taxa,
#'                                    total_richness_per_clades = richness_per_clades$taxa,
#'                                    seed = 1234, verbose = FALSE)
#' grafted_phy <- graft_missing_clades(phy = grafted_taxa, missing_taxa_list = lists$clades,
#'                                     total_richness_per_clades = richness_per_clades$clades,
#'                                     seed = 1234, verbose = FALSE)
#'


build_total_richness_per_clades <- function (data, # Data frame, one row per clade
                                             phy, # Backbone phylogeny
                                             missing_taxa_list, # Missing taxa list, or list of two missing taxa lists ($taxa and $clades)
                                             return_richness_table = TRUE, # Return the expected richness and sampling fraction of each clade as an attribute
                                             verbose = TRUE) # Print the conversions of node IDs into tips
{
  ## 1/ Check the inputs ####

  if (!inherits(phy, "phylo")) { stop("'phy' must be a phylogeny of class 'phylo'") }
  if (!is.data.frame(data)) { stop("'data' must be a data frame") }
  if (nrow(data) == 0) { stop("'data' has no rows") }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) { stop("'verbose' must be TRUE or FALSE") }
  if (!is.logical(return_richness_table) || length(return_richness_table) != 1 || is.na(return_richness_table)) { stop("'return_richness_table' must be TRUE or FALSE") }
  if (missing(missing_taxa_list) || !is.list(missing_taxa_list))
  {
    stop("'missing_taxa_list' must be provided: a 'missing_taxa_list', or the list of the two 'missing_taxa_list' ($taxa and $clades) built from the same table or taxonomic tables")
  }

  ## 1.1/ Structure of the missing_taxa_list ####

  # A single missing_taxa_list (with Target_clade or Sister_clade), or a list with the two ($taxa and $clades)
  is_pair <- is.null(missing_taxa_list$Taxon_name) && all(c("taxa", "clades") %in% names(missing_taxa_list))
  lists <- list()
  if (is_pair)
  {
    lists$taxa <- missing_taxa_list$taxa
    lists$clades <- missing_taxa_list$clades
  } else {
    has_target <- !is.null(missing_taxa_list$Target_clade)
    has_sister <- !is.null(missing_taxa_list$Sister_clade)
    if (has_target && has_sister) { stop("'missing_taxa_list' has both 'Target_clade' and 'Sister_clade': use one list for each type of taxa ($taxa and $clades)") }
    if (!has_target && !has_sister) { stop("'missing_taxa_list' must be a 'missing_taxa_list' (with 'Target_clade' or 'Sister_clade'), or a list of two ('$taxa' and '$clades')") }
    lists[[if (has_target) "taxa" else "clades"]] <- missing_taxa_list
  }
  items <- c(taxa = "Target_clade", clades = "Sister_clade")
  for (type in names(lists))
  {
    if (!is.list(lists[[type]])) { stop("'missing_taxa_list$", type, "' must be a 'missing_taxa_list'") }
    validate_missing_taxa_list(lists[[type]], clade_item = items[[type]])
  }
  all_taxa <- unlist(lapply(lists, function (x) x$Taxon_name), use.names = FALSE)
  if (anyDuplicated(all_taxa)) { stop("The taxa of the two lists must be different. Duplicated: ", paste(unique(all_taxa[duplicated(all_taxa)]), collapse = ", ")) }
  absent <- unique(unlist(lapply(names(lists), function (type) unlist(lists[[type]][[items[[type]]]]))))
  absent <- setdiff(absent, phy$tip.label)
  if (length(absent) > 0) { stop("Some tips of 'missing_taxa_list' are missing from the backbone phylogeny: ", paste(utils::head(absent, 10), collapse = ", ")) }

  ## 1.2/ Clades and names of the missing_taxa_list ####

  # Clades of the lists (target and sister clades, and excluded subclades), to harmonize the tips used for each clade
  known <- list()
  pool_tips <- list(); pool_names <- character(0); pool_nodes <- numeric(0)
  add_to_pool <- function (tips, name)
  {
    if (is.null(name) || is.na(name)) { return(invisible(NULL)) }
    pool_tips[[length(pool_tips) + 1]] <<- tips
    pool_names <<- c(pool_names, name)
    pool_nodes <<- c(pool_nodes, get_clade_node(phy, tips))
  }
  for (type in names(lists))
  {
    list_i <- lists[[type]]
    item <- items[[type]]
    clades_i <- unname(list_i[[item]])
    known <- c(known, clades_i)
    clade_names_i <- list_i[[paste0(item, "_name")]]
    if (!is.null(clade_names_i))
    {
      for (k in seq_along(clades_i)) { add_to_pool(clades_i[[k]], unname(clade_names_i)[k]) }
    }
    for (subclades in list_i$Excluded_subclades)
    {
      known <- c(known, unname(subclades))
      for (k in seq_along(subclades)) { add_to_pool(subclades[[k]], names(subclades)[k]) }
    }
  }
  lookup <- build_tips_lookup(phy, known)
  if (length(pool_names) > 0)
  {
    # Consistent names across the lists: the same clade always has the same name
    check_clade_labels(labels = pool_names, nodes = pool_nodes, tips = pool_tips, taxon_name = rep("missing_taxa_list", length(pool_names)),
                       what = "the names of the clades of 'missing_taxa_list'")
  }
  # Taxonomic groups of the lists built from taxonomic tables: names of all the groups of the backbone taxonomy, with lower priority
  groups <- attr(missing_taxa_list, "taxonomic_groups")
  group_tips <- list()
  if (is.data.frame(groups) && all(c("Group", "clade_tips") %in% names(groups)))
  {
    for (g in seq_len(nrow(groups)))
    {
      tips_g <- trimws(strsplit(as.character(groups$clade_tips[g]), ",", fixed = TRUE)[[1]])
      if (length(tips_g) > 0 && all(tips_g %in% phy$tip.label)) { group_tips[[as.character(groups$Group[g])]] <- tips_g }
    }
  }

  ## 2/ Check the table ####

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  names(data) <- trimws(names(data))
  if (anyDuplicated(names(data))) { stop("Some columns have the same name: ", paste(unique(names(data)[duplicated(names(data))]), collapse = ", ")) }
  data[] <- lapply(data, function (x) { x <- trimws(as.character(x)); x[!is.na(x) & x %in% c("", "NA")] <- NA; x })

  clade_columns <- intersect(c("Clade", "Target_clade", "Sister_clade"), names(data))
  if (length(clade_columns) == 0) { stop("The table must have a column 'Clade' (or 'Target_clade', or 'Sister_clade')") }
  if (length(clade_columns) > 1) { stop("The table must have only one of the columns 'Clade', 'Target_clade' and 'Sister_clade'. Found: ", paste(clade_columns, collapse = ", ")) }
  clade_column <- clade_columns
  name_column <- paste0(clade_column, "_name")
  if (!("total_richness" %in% names(data))) { stop("The table must have a column 'total_richness'") }
  unknown <- setdiff(names(data), c(clade_column, name_column, "total_richness"))
  if (length(unknown) > 0)
  {
    warning("Unknown columns ignored: ", paste0("'", unknown, "'", collapse = ", "), ".\n",
            "Recognized columns: ", paste(c(clade_column, name_column, "total_richness"), collapse = ", "), ".", call. = FALSE)
  }
  n_rows <- nrow(data)

  # Total richness: whole numbers (NA allowed)
  richness_text <- data$total_richness
  richness <- suppressWarnings(as.numeric(richness_text))
  bad_richness <- which(!is.na(richness_text) & (is.na(richness) | richness < 0 | richness != round(richness)))
  if (length(bad_richness) > 0)
  {
    stop("'total_richness' must be a non-negative whole number (or NA) in the row(s): ",
         paste0(bad_richness, " ('", richness_text[bad_richness], "')", collapse = ", "))
  }
  empty_clades <- which(is.na(data[[clade_column]]))
  if (length(empty_clades) > 0) { stop("'", clade_column, "' is empty in the row(s): ", paste(empty_clades, collapse = ", ")) }

  ## 3/ Identify the clades ####

  # A clade is written with tips, with a node ID, or with a name. Node IDs are converted into tips (reusing the tips of the lists)
  recorder <- new.env(parent = emptyenv())
  recorder$rows <- list()
  clade_tips <- vector("list", n_rows)
  name_from_cell <- rep(NA_character_, n_rows)
  for (i in seq_len(n_rows))
  {
    tokens <- trimws(strsplit(data[[clade_column]][i], ",", fixed = TRUE)[[1]])
    if (length(tokens) == 0 || any(tokens == "")) { stop("'", clade_column, "' of row ", i, " has an empty tip (check the ',' separators): ", data[[clade_column]][i]) }
    if (anyDuplicated(tokens)) { stop("'", clade_column, "' of row ", i, " has repeated tips: ", data[[clade_column]][i]) }
    if (length(tokens) == 1 && !grepl("^[0-9]+$", tokens) && !(tokens %in% phy$tip.label))
    {
      # A name: from the clades of the lists, then from the taxonomic groups
      node_name <- if (tokens %in% pool_names) pool_nodes[match(tokens, pool_names)] else NA
      if (!is.na(node_name))
      {
        tips <- lookup[[as.character(node_name)]]
        if (is.null(tips)) { tips <- pool_tips[[match(tokens, pool_names)]] }
      } else if (!is.null(group_tips[[tokens]])) {
        tips <- group_tips[[tokens]]
        node_group <- as.character(get_clade_node(phy, tips))
        if (!is.null(lookup[[node_group]])) { tips <- lookup[[node_group]] }
      } else {
        stop("'", clade_column, "' of row ", i, " ('", tokens, "') is not a tip of the backbone, a node ID, nor the name of a clade of 'missing_taxa_list'",
             if (length(group_tips) > 0) " or of its taxonomic groups" else "", ".\n",
             "Check the spelling. A clade is written with tips (e.g., \"t26, t18\"), a node ID (e.g., \"34\"), or a name.")
      }
      clade_tips[[i]] <- tips
      name_from_cell[i] <- tokens
    } else {
      tips <- resolve_clade_tokens(tokens, phy = phy, item = clade_column, who = paste0("row ", i), recorder = recorder, lookup = lookup)
      absent_tips <- setdiff(tips, phy$tip.label)
      if (length(absent_tips) > 0)
      {
        stop("Some tips of '", clade_column, "' are missing from the backbone phylogeny (row ", i, "): ", paste(absent_tips, collapse = ", "),
             "\nCheck the spelling of the tips, or the separators.")
      }
      clade_tips[[i]] <- tips
    }
  }
  report_node_conversions(recorder$rows, phy = phy, verbose = verbose)
  clade_nodes <- vapply(clade_tips, function (tips) get_clade_node(phy, tips), numeric(1))

  # A clade is listed once (clades are matched by their MRCA node)
  if (anyDuplicated(clade_nodes))
  {
    dup <- unique(clade_nodes[duplicated(clade_nodes)])
    stop("Some clades are listed more than once in the table (clades are matched by their MRCA node):\n",
           paste0("  - rows ", vapply(dup, function (nd) paste(which(clade_nodes == nd), collapse = ", "), character(1)),
                  " (tips: ", vapply(dup, function (nd) paste(clade_tips[[which(clade_nodes == nd)[1]]], collapse = ", "), character(1)), ")", collapse = "\n"))
  }

  ## 4/ Names of the clades ####

  # Name written in the table, or used as the clade
  clade_names <- if (name_column %in% names(data)) data[[name_column]] else rep(NA_character_, n_rows)
  conflict <- which(!is.na(clade_names) & !is.na(name_from_cell) & clade_names != name_from_cell)
  if (length(conflict) > 0)
  {
    stop("The clade is written with a name in '", clade_column, "' but has a different '", name_column, "' in the row(s): ",
         paste0(conflict, " ('", name_from_cell[conflict], "' vs '", clade_names[conflict], "')", collapse = ", "))
  }
  clade_names[is.na(clade_names)] <- name_from_cell[is.na(clade_names)]
  # Names must be consistent in the table and with the names of 'missing_taxa_list'
  check_clade_labels(labels = c(clade_names, pool_names), nodes = c(clade_nodes, pool_nodes), tips = c(clade_tips, pool_tips),
                     taxon_name = c(paste0("row ", seq_len(n_rows)), rep("missing_taxa_list", length(pool_names))),
                     what = paste0(name_column, " and the names of missing_taxa_list"))
  # Names of the clades of the lists are used when the table has none
  for (i in which(is.na(clade_names)))
  {
    if (clade_nodes[i] %in% pool_nodes) { clade_names[i] <- pool_names[match(clade_nodes[i], pool_nodes)] }
  }
  # Same tips as in the lists
  for (i in seq_len(n_rows))
  {
    reused <- lookup[[as.character(clade_nodes[i])]]
    if (!is.null(reused)) { clade_tips[[i]] <- reused }
  }

  ## 5/ Add the clades of the lists that are not in the table ####

  # Every clade of the lists must be listed: they are added with an unknown total richness
  needed_tips <- list(); needed_names <- character(0)
  for (type in names(lists))
  {
    list_i <- lists[[type]]
    item <- items[[type]]
    clades_i <- unname(list_i[[item]])
    names_i <- list_i[[paste0(item, "_name")]]
    for (k in seq_along(clades_i))
    {
      needed_tips[[length(needed_tips) + 1]] <- clades_i[[k]]
      needed_names <- c(needed_names, if (is.null(names_i)) NA_character_ else unname(names_i)[k])
    }
  }
  needed_nodes <- vapply(needed_tips, function (tips) get_clade_node(phy, tips), numeric(1))
  to_add <- which(!(needed_nodes %in% clade_nodes) & !duplicated(needed_nodes))
  n_added <- length(to_add)
  provided_rows <- seq_len(n_rows)
  if (n_added > 0)
  {
    clade_tips <- c(clade_tips, needed_tips[to_add])
    clade_nodes <- c(clade_nodes, needed_nodes[to_add])
    clade_names <- c(clade_names, needed_names[to_add])
    richness <- c(richness, rep(NA_real_, n_added))
  }
  n_clades <- length(clade_tips)

  ## 6/ Check the richness against the full sequence of grafting ####

  # The taxa to graft as a missing clade are placed on the stem of their sister clade: they belong to the parent of this clade.
  # For the check, all the taxa are therefore described as taxa with a target clade (the sister clades are replaced by their parent).
  # The clades used for the check without being in the table are added with an unknown richness.
  all_clades <- list()
  if (!is.null(lists$taxa)) { all_clades <- c(all_clades, unname(lists$taxa$Target_clade)) }
  check_tips <- clade_tips; check_values <- richness
  if (!is.null(lists$clades))
  {
    for (tips in unname(lists$clades$Sister_clade))
    {
      node <- get_clade_node(phy, tips)
      parent <- phy$edge[match(node, phy$edge[, 2]), 1]
      if (is.na(parent)) { stop("A sister clade of 'missing_taxa_list' is the root of the backbone, so it has no stem branch to graft on") }
      parent_tips <- node_to_tips(phy, parent)
      all_clades <- c(all_clades, list(parent_tips))
      if (!(parent %in% vapply(check_tips, function (x) get_clade_node(phy, x), numeric(1))))
      {
        check_tips <- c(check_tips, list(parent_tips)); check_values <- c(check_values, NA_real_)
      }
    }
  }
  info <- tryCatch(
    build_richness_info(phy = phy, missing_taxa_names = all_taxa, missing_clades = all_clades, mode = "target",
                        total_richness_per_clades = list(Target_clade = check_tips, total_richness = check_values)),
    error = function (e)
    {
      stop("The total richness is not consistent with the species of the backbone and the taxa to graft ",
           "(all the taxa to graft are counted in the clades they belong to, whatever their type):\n", conditionMessage(e), call. = FALSE)
    })
  richness_table <- info$table[seq_len(n_clades), , drop = FALSE]
  richness_table <- data.frame(Clade_name = clade_names,
                               Clade_tips = vapply(clade_tips, paste, character(1), collapse = ", "),
                               total_richness_provided = richness,
                               total_richness = richness_table$total_richness,
                               nb_backbone_species = richness_table$nb_backbone_species,
                               nb_missing_species = richness_table$nb_missing_species,
                               nb_unsampled_species = richness_table$nb_unsampled_species,
                               sampling_fraction = richness_table$sampling_fraction,
                               stringsAsFactors = FALSE)

  ## 7/ Build the outputs ####

  values <- if (all(is.na(richness) | richness <= .Machine$integer.max)) as.integer(richness) else richness
  build_one <- function (item)
  {
    out <- list()
    out[[item]] <- clade_tips
    if (any(!is.na(clade_names))) { out[[paste0(item, "_name")]] <- clade_names }
    out$total_richness <- values
    out
  }
  result <- lapply(lists, function (x) NULL)
  for (type in names(lists)) { result[[type]] <- build_one(items[[type]]) }
  if (!is_pair) { result <- result[[1]] }
  if (return_richness_table) { attr(result, "richness_table") <- richness_table }

  if (verbose)
  {
    message(n_clades, " clade(s) in total richness: ", sum(!is.na(richness)), " with a value, ", sum(is.na(richness)), " with an unknown richness",
            if (n_added > 0) paste0(" (including ", n_added, " clade(s) of 'missing_taxa_list' added to the ", n_rows, " of the table)") else "",
            if (return_richness_table) ". See attr(, \"richness_table\") for the expected sampling fraction of the clades." else ".")
  }
  return(result)
}
