##### Functions to build a missing_taxa_list from an internal data.frame #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Build a missing_taxa_list from an internal data.frame ####

#' @title Build a `missing_taxa_list` from an internal data.frame of missing taxa
#'
#' @description Build a valid `missing_taxa_list` from a data.frame describing the taxa to graft (one row per taxon, e.g., read from a .csv file),
#'   for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'   Lists of tips (or node IDs, see Details) are written in the cells of the table, and the function parses them and names every item after the taxa (as required by the grafting functions),
#'   and runs all the checks of validity against the backbone phylogeny.
#'
#'   Use [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read the table directly from a file.
#'
#'   Clades are defined by tips or node IDs, which allows for fine-tuning the constraints of each taxon (see 'Fine-tuning the clades' in Details).
#'   For large phylogenies with neatly defined taxonomic levels, taxonomic tables are likely easier to prepare:
#'   see [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy] and [`phylograftR::import_taxonomic_table()`][import_taxonomic_table].
#'
#' @param data A data frame, with one row per taxon to graft.
#'   Its columns are the items of `missing_taxa_list` as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'   See Details for a full description. All the columns are read as text.  Blank cells and `NA` are missing values.
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which the taxa will be grafted. It is used to check the tips and clades, and to convert node IDs into tips.
#' @param verbose Logical. Whether to print the conversions of node IDs into tips (if any). Default is `TRUE`.
#'
#' @details
#' The columns of `data` are the items of `missing_taxa_list` as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#' Column names must be written exactly as below. Tips are separated by commas (`,`), and the subclades by semicolons (`;`). Spaces around separators are ignored.
#' Values are `TRUE` or `FALSE` (also `T`, `F`, `yes`, `no`, `1`, `0`), or `NA`.
#'
#'  * `Taxon_name` (Required): Names of the taxa to graft. Names must be unique and absent from `phy$tip.label`.
#'  * `Target_clade` or `Sister_clade` (Required): Tip labels defining a target clade per taxon, i.e., the clade where the taxon is known to belong (e.g., its genus);
#'    or a sister clade, i.e., the clade on whose stem branch the taxon is grafted, as its sister.
#'    A pair of tips defines the clade as their MRCA node (e.g., `"t26, t18"`).
#'    A single tip defines the "clade" as the terminal branch of this tip (e.g., `"t15"`).
#'    A node ID of `phy` (digits only, e.g., `"34"`) also defines the clade, see 'Node IDs' below.
#'    Use `Sister_clade` to prepare a list for [`phylograftR::graft_missing_clades()`][graft_missing_clades]:
#'    the columns specific to [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] below are then ignored, with a warning.
#'  * `include_target_clade_stem`: (Optional) Logical values. For each taxon, whether the stem branch of its target clade is available for grafting.
#'    `NA` values and blank cells use the global `include_target_clade_stem` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#'  * `Excluded_subclades`: (Optional) For each taxon, tip labels defining the subclades of its target clade where the taxon is known not to belong.
#'   Each subclade is defined by a pair of tips, a single tip, or a node ID.
#'   Ex: `"t10, t19; 50; t7"` for three subclades, the two first being clades with multiple species defined by a pair of tips (`"t10, t19"`), and an internal node ID (`"50"`).
#'   The last one is a single tip define by its tip label (`"t7"`). Use blank cells for taxa without excluded subclades.
#'  * `exclude_subclades_stem`: (Optional) Logical values. For each taxon, whether the stem of each of its excluded subclades is also excluded,
#'    with one value per subclade of `Excluded_subclades` separated by semicolons, or a single value used for all of them.
#'    Ex: `"FALSE; TRUE; TRUE"`.
#'    `NA` values and blank cells use the global `exclude_subclades_stem` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#'  * `allowed_states`: (Optional) For each taxon, the states of a `simmap` in which it can be grafted (see the `simmap` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#'    or [`phylograftR::graft_missing_clades()`][graft_missing_clades]), separated by commas. It is accepted with both `Target_clade` and `Sister_clade`.
#'    Ex: `"blue, green"`. Blank cells and `NA` have no taxon-level value, so the global `allowed_states` argument of the grafting function is used for the taxon.
#'    States are matched exactly with the states of the `simmap` (e.g., case matters) when grafting, and so cannot contain commas.
#'    If the `simmap` argument is a list of the maps of several traits, give the states of each trait separated by semicolons, in the order of the traits (e.g., `"A1, A2;; C3"` for
#'    three traits: no states for the second trait, so the global `allowed_states` of the trait, or all its states, are used), or name the traits (e.g., `"color=A1,A2; elevation=C3"`).
#'    The names of the states are not checked here as the `simmap` is not provided: the grafting function checks them.
#'
#' Optional columns that do not match items in `missing_taxa_list` can be used to provide names of clades:
#'  * `Target_clade_name` or `Sister_clade_name` (Optional): Label of the target (or sister) clade.
#'    It is only used to check that the clades are defined consistently (see 'Labels of the clades' section below).
#'  * `Excluded_subclades_names` (Optional): Labels of the excluded subclades, separated by semicolons, in the same order as in `Excluded_subclades`.
#'    Ex: `"group_1; group_2; tip_t7"`. There must be one label per subclade.
#'    If provided, they are used as names of the subclades in the `Excluded_subclades` item of the list, and in the messages of the grafting functions.
#'
#'  Other columns are ignored, with a warning.
#'
#' ### Node IDs
#'
#' Clades (`Target_clade`, `Sister_clade`, and each subclade of `Excluded_subclades`) are typically defined by (pairs of) tips.
#' However, the function also accepts clades defined by the ID of a node of `phy` instead of tips, in the same table, and even in the same cell
#' (e.g., `"50; t5, t12; t7"` for three excluded subclades defined respectively by an internal node, a pair of tips, a terminal tip label).
#' Node IDs are integers: `1` to N, with N being the number of tips (`ape::Ntip(phy)`) correspond to terminal nodes, and thus define single-tip clades.
#' Following IDs are the internal nodes, as shown by `ape::nodelabels()`, and correspond to clades with multiple descending species.
#'
#' A value made only of digits is a node ID, so **tip labels must include at least one other character than digits** to be distinguinshed.
#' A clade is defined by tips or by one node ID, not both, and one node ID defines the clade alone.
#' Node IDs are converted into tips before any other check, as done by [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips],
#' and the conversions are printed if `verbose = TRUE`.
#'
#' The warnings and errors that follow, and the grafting functions, therefore use the tips that are traced back to your node IDs in the messages.
#' If a clade is already defined by tips in the table, a node ID of the same clade is converted into the same tips.
#' Convert `total_richness_per_clades` (if you use it) at the same time as the list, with [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips],
#' so the tips of the clades match between both objects.
#'
#' ### Fine-tuning the clades
#'
#' Defining clades with tips or node IDs gives full control on the constraints of each taxon. **Any node of `phy` can be used** as a `Target_clade`,
#' as a `Sister_clade`, or as an excluded subclade, simply by writing its node ID (e.g., `"34"`), or two tips whose MRCA is the node (e.g., `"t8, t7"`).
#' The node does not need to match a taxonomic group: it can be a clade without a name, a clade that groups several taxonomic units (e.g., two neighboring genera),
#' a part of a taxonomic group only, or a clade you identified from the phylogeny itself or from another source (e.g., a molecular study).
#' Clades can be nested at any depth, and each taxon has its own set of excluded subclades, with or without their stem.
#'
#' This is the best route to fine-tune the grafting of a limited number of taxa, but the clades have to be written for each taxon.
#' When the backbone is large, and taxonomic levels (e.g., Family, Subfamily, Tribe, Genus) are defined for its tips and for the taxa to graft, taxonomic tables are the go-to option:
#' clades are retrieved from their taxonomic names, which ease the attribution of clades across many taxa.
#' See [`phylograftR::import_taxonomic_table()`][import_taxonomic_table] and [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy].
#'
#' Both routes return the same `missing_taxa_list`, so they can be combined: build a list from taxonomic tables, and then adjust the clades of a few taxa with node IDs or tips.
#'
#' ### Labels of the clades
#'
#' Labels (`Target_clade_name`, `Sister_clade_name`, and `Excluded_subclades_names`) are optional but, when used, they identify the clades
#' across all the rows of the table: the same clade must always have the same label, and different clades must have different labels.
#' A clade is identified by the MRCA node of its tips in `phy`, so the same clade can be defined by different pairs of tips, as long as they share this MRCA node.
#' The function stops if a label is used for different clades, or if a clade is given different labels, and lists the labels and tips in conflict.
#' This catches typos in the tips or in the labels. Subclades without a label are not checked.
#'
#' ### Validity checks
#'
#' The function stops with an informative message if:
#'  * Required columns are missing, or one `Taxon_name` is missing, duplicated, or already in `phy`.
#'  * A clade has no tips, or tips are repeated, or are absent from `phy`.
#'  * The number of values in `exclude_subclades_stem`, or of labels in `Excluded_subclades_names`, does not match the number of excluded subclades.
#'  * Labels are inconsistent (see above).
#'  * A cell of `allowed_states` has an empty state (check the commas) or repeated states.
#'  * An excluded subclade is identical to, or contains, the target clade of its taxon (no branch would remain).
#'    A warning is issued if it is outside of the target clade (the exclusion has no effect).
#'
#' The messages give the names of the problematic taxa and the numbers of the rows in the table.
#'
#' @return A `missing_taxa_list`: a list, ready to be used in the grafting functions, with the following items:
#'  * `$Taxon_name` Character vector.
#'  * `$Target_clade` (or `$Sister_clade`) List of character vectors of tips.
#'  * `$Target_clade_name` (or `$Sister_clade_name`) Only if the column was provided. Character vector of labels (`NA` if none). Not used for the grafting, but checked against the names of `total_richness_per_clades`, if provided (see the grafting functions).
#'  * `$include_target_clade_stem` Only if the column was provided. Logical vector (`NA` uses the global argument).
#'  * `$Excluded_subclades` Only if the column was provided. List (one element per taxon) of lists of character vectors of tips. `NULL` for taxa with no excluded subclades.
#'    The subclades are named after their labels if `Excluded_subclades_names` was provided.
#'  * `$exclude_subclades_stem` Only if the column was provided. List (one element per taxon) of logical vectors (one value per excluded subclade), or `NULL`.
#'  * `$allowed_states` Only if the column was provided. List (one element per taxon) of character vectors of states, or `NULL` for taxa without a taxon-level value (the global `allowed_states` argument is used).
#'
#'  Every item is named after the taxa (`$Taxon_name`), in the same order.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read the table from a file,
#'  * [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy] and [`phylograftR::import_taxonomic_table()`][import_taxonomic_table] to define the clades with taxonomic names instead of tips or node IDs,
#'  * [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] to convert the node IDs of a list built 'by hand' into pairs of tips defining MRCA,
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] that use the list for grafting.
#'
#' @importFrom ape getMRCA
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Load a table of missing taxa as a data.frame in R
#' file_path <- system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")
#' missing_taxa_table <- read.csv(file_path)
#'
#' # Convert it into a 'missing_taxa_list', checking it against the backbone phylogeny
#' missing_taxa_list <- build_missing_taxa_list(data = missing_taxa_table, phy = toy_backbone)
#' # Items are named after the taxa, tips are split, and excluded subclades are named
#' str(missing_taxa_list, 2)
#'
#' # A table can also be built in R and then converted to a valid 'missing_taxa_list' object
#' my_table <- data.frame(Taxon_name = c("M1", "M2"),
#'                        Target_clade = c("t26, t18", "t5"),
#'                        Excluded_subclades = c("t10, t19; t26, t4", NA),
#'                        Excluded_subclades_names = c("group_1; group_2", NA))
#' my_missing_taxa_list <- build_missing_taxa_list(data = my_table, phy = toy_backbone)
#' str(my_missing_taxa_list, 2)
#'
#' # Clades can also be defined by node IDs (digits only), alone or mixed with tips.
#' # Node IDs are converted into tips, and the conversions are printed.
#' ape::plot.phylo(toy_backbone, cex = 0.7, no.margin = TRUE)
#' ape::nodelabels(cex = 0.6, frame = "none", col = "blue")
#' my_table_nodes <- data.frame(Taxon_name = c("M1", "M2"),
#'                              Target_clade = c("34", "t5, t12"),
#'                              Excluded_subclades = c("50", NA))
#' str(build_missing_taxa_list(data = my_table_nodes, phy = toy_backbone), 2)
#'

build_missing_taxa_list <- function (data, # Data frame, one row per taxon
                                     phy, # Backbone phylogeny
                                     verbose = TRUE) # Print the conversions of node IDs into tips
{
  ## 2.1/ Check the inputs ####

  if (!inherits(phy, "phylo")) { stop("'phy' must be a phylogeny of class 'phylo'") }
  if (!is.data.frame(data)) { stop("'data' must be a data frame") }
  if (nrow(data) == 0) { stop("'data' has no rows") }

  ## 2.1.1/ Check columns ####

  # Columns recognized with each type of clade column. 'Target_clade' goes with the items of graft_missing_taxa(), 'Sister_clade' with those of graft_missing_clades()
  columns_taxa <- c("Taxon_name", "Target_clade", "Target_clade_name", "include_target_clade_stem",
                    "Excluded_subclades", "Excluded_subclades_names", "exclude_subclades_stem", "allowed_states")
  columns_clades <- c("Taxon_name", "Sister_clade", "Sister_clade_name", "allowed_states")
  all_columns <- unique(c(columns_taxa, columns_clades))

  # Clean the table: names of the columns without extra spaces, and every cell as trimmed text. Blank cells and "NA" are missing values (NA)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  names(data) <- trimws(names(data))
  if (anyDuplicated(names(data))) { stop("Some columns have the same name: ", paste(unique(names(data)[duplicated(names(data))]), collapse = ", ")) }
  data[] <- lapply(data, function (x) { x <- trimws(as.character(x)); x[!is.na(x) & x %in% c("", "NA")] <- NA; x })

  if (!("Taxon_name" %in% names(data))) { stop("The table must have a column 'Taxon_name'") }

  # Exactly one of 'Target_clade' and 'Sister_clade' is required: it decides which list is built (for graft_missing_taxa() or graft_missing_clades())
  has_target <- "Target_clade" %in% names(data)
  has_sister <- "Sister_clade" %in% names(data)
  if (has_target && has_sister) { stop("The table has both 'Target_clade' and 'Sister_clade': use only one, depending on the grafting function (graft_missing_taxa() or graft_missing_clades())") }
  if (!has_target && !has_sister) { stop("The table must have a column 'Target_clade' (for graft_missing_taxa()) or 'Sister_clade' (for graft_missing_clades())") }
  clade_item <- if (has_target) "Target_clade" else "Sister_clade"
  label_item <- paste0(clade_item, "_name")
  allowed_columns <- if (has_target) columns_taxa else columns_clades

  # Columns not recognized (typos) or not used with the clade column provided: warn, and keep only the usable columns
  unknown <- setdiff(names(data), all_columns)
  not_used <- setdiff(intersect(names(data), all_columns), allowed_columns)
  if (length(unknown) > 0)
  {
    warning("Unknown columns ignored: ", paste0("'", unknown, "'", collapse = ", "), ".\n",
            "Recognized columns: ", paste(allowed_columns, collapse = ", "), ".", call. = FALSE)
  }
  if (length(not_used) > 0)
  {
    warning("Columns not used with '", clade_item, "' ignored: ", paste0("'", not_used, "'", collapse = ", "), ".", call. = FALSE)
  }
  data <- data[, intersect(names(data), allowed_columns), drop = FALSE]

  if ("Excluded_subclades_names" %in% names(data) && !("Excluded_subclades" %in% names(data)))
  {
    stop("The column 'Excluded_subclades_names' needs the column 'Excluded_subclades'")
  }

  ## 2.1.2/ Check Taxa ####

  n_taxa <- nrow(data)
  row_id <- seq_len(n_taxa)
  taxon_name <- data$Taxon_name
  if (anyNA(taxon_name)) { stop("'Taxon_name' is empty in the row(s): ", paste(row_id[is.na(taxon_name)], collapse = ", ")) }
  if (anyDuplicated(taxon_name)) { stop("'Taxon_name' must be unique. Duplicated: ", paste(unique(taxon_name[duplicated(taxon_name)]), collapse = ", ")) }
  in_phy <- taxon_name %in% phy$tip.label
  if (any(in_phy))
  {
    stop("'Taxon_name' must be absent from the backbone phylogeny (these are taxa to graft). Present in 'phy': ", paste(taxon_name[in_phy], collapse = ", "))
  }

  # Label of a taxon in messages
  label_taxon <- function (i) paste0("'", taxon_name[i], "' (row ", i, ")")

  ## 2.1.3/ Convert the node IDs into tips ####

  # GOAL: A clade can be written in a cell with tips ("t8, t7") or with a node ID of 'phy' ("34"). Here, every node ID is replaced by tips in the TEXT of the cell,
  # so the rest of the function (parsing and checks, from 2.2) only ever sees tips. Cells that already contain tips are left unchanged.
  # Cells with problems (e.g., empty tips like "t1,,t2") are left as they are: they are reported by the checks of the tips in the next sections.

  # Example for a cell of 'Excluded_subclades': "50; t5, t12; t7"
  #   split on ";"  ->  "50"  |  "t5, t12"  |  "t7"
  #   node ID "50" is replaced by a pair of tips  ->  "t13,t27"  |  "t5,t12"  |  "t7"
  #   cell becomes: "t13,t27;t5,t12;t7"

  # 'recorder' is a small environment (i.e., a list that can be modified from inside the functions below) where each conversion is noted
  # (node ID, tips used, column, taxon). It is used to print the conversions at the end of this section.
  recorder <- new.env(parent = emptyenv())
  recorder$rows <- list()

  # Step A. Find the clades that are ALREADY written with tips in the table, to build a 'lookup' (node of the clade -> tips used to write it).
  # It ensures the same clade is always written with the same tips: if a row writes "t26, t18" (MRCA node 32) and another row writes "32",
  # the node ID 32 is converted into "t26, t18" (the tips already used) and not into another pair of tips with the same MRCA node.

  # This internal function splits the text of a cell into vectors of tips (here, node IDs are just pieces of text that are skipped when building the lookup).
  split_cell <- function (cell, split_subclades)
  {
    if (is.na(cell)) { return(list()) }
    # Excluded_subclades: several clades per cell, separated by ";". Target_clade/Sister_clade: one clade per cell
    pieces <- if (split_subclades) trimws(strsplit(cell, ";", fixed = TRUE)[[1]]) else cell
    # Each clade is a vector of tips separated by ","
    lapply(pieces, function (piece) trimws(strsplit(piece, ",", fixed = TRUE)[[1]]))
  }

  # Record all the clades of the table, as a flat list of vectors of tips: target (or sister) clades, then excluded subclades
  known <- unlist(lapply(seq_len(n_taxa), function (i) split_cell(data[[clade_item]][i], FALSE)), recursive = FALSE)
  if ("Excluded_subclades" %in% names(data))
  {
    known <- c(known, unlist(lapply(seq_len(n_taxa), function (i) split_cell(data$Excluded_subclades[i], TRUE)), recursive = FALSE))
  }
  # Keep the clades written with tips (valid in 'phy'), and associate each to its MRCA node
  lookup <- build_tips_lookup(phy, known)

  # Step B. Build a function to rewrite one cell, replacing its node IDs by tips.

  # Internal function to parse text of a cell and replace its node IDs by tips
  # 'cell': text of the cell. 'item': name of the column, and 'who': name of the taxon (both for the messages and the record of the conversions).
  # 'split_subclades': TRUE if the cell has several clades separated by ";" ('Excluded_subclades'), FALSE if it has one clade ('Target_clade' or 'Sister_clade').
  convert_cell <- function (cell, item, who, split_subclades)
  {
    if (is.na(cell)) { return(cell) }
    pieces <- if (split_subclades) trimws(strsplit(cell, ";", fixed = TRUE)[[1]]) else cell
    converted <- vapply(pieces, function (piece)
    {
      # 'piece' is one clade, e.g., "t5, t12" or "50". Its tokens are the tips, or the node ID
      tokens <- trimws(strsplit(piece, ",", fixed = TRUE)[[1]])
      if (length(tokens) == 0 || any(tokens == "")) { return(piece) } # Left unchanged: reported later by parse_tips()
      # resolve_clade_tokens() returns the tips unchanged if the clade is written with tips; and converts it if the clade is a node ID.
      # It stops if tips and node IDs are mixed in the same clade, or if the node ID is invalid.
      # The tips already used for this node ('lookup') are preferred over the default pair of tips of the node.
      paste(resolve_clade_tokens(tokens, phy = phy, item = item, who = who, recorder = recorder, cell = piece, lookup = lookup), collapse = ",")
    }, character(1), USE.NAMES = FALSE)
    # Reassemble the clades of the cell with ";"
    paste(converted, collapse = ";")
  }

  # Step C. Apply it to every row of the columns defining clades: the text of these columns is replaced
  data[[clade_item]] <- vapply(seq_len(n_taxa), function (i) convert_cell(data[[clade_item]][i], clade_item, taxon_name[i], FALSE), character(1))
  if ("Excluded_subclades" %in% names(data))
  {
    data$Excluded_subclades <- vapply(seq_len(n_taxa), function (i) convert_cell(data$Excluded_subclades[i], "Excluded_subclades", taxon_name[i], TRUE), character(1))
  }

  # Print the conversions (nothing is printed if there was no node ID, or if verbose = FALSE)
  report_node_conversions(recorder$rows, phy = phy, verbose = verbose)

  ## 2.2/ Clades of the taxa: Target_clade or Sister_clade ####

  # From here, all the clades are written with tips: parse the text of each cell into a vector of tips.
  # We now check that the tips are in 'phy', and identify the node of each clade (tip node for one tip, MRCA node for two or more tips)
  clade_cells <- data[[clade_item]]
  if (anyNA(clade_cells)) { stop("'", clade_item, "' is empty for: ", paste(label_taxon(which(is.na(clade_cells))), collapse = ", ")) }
  clade_tips <- lapply(seq_len(n_taxa), function (i) parse_tips(clade_cells[i], label_taxon(i), clade_item))
  check_tips_in_phy(clade_tips, phy, taxon_name, clade_item)
  clade_nodes <- vapply(clade_tips, function (tips) get_clade_node(phy, tips), numeric(1))

  # Optional labels of the clades: the same clade (same MRCA node) must have the same label across taxa, and different clades different labels
  clade_labels <- NULL
  if (label_item %in% names(data))
  {
    clade_labels <- data[[label_item]]
    check_clade_labels(labels = clade_labels, nodes = clade_nodes, tips = clade_tips, taxon_name = taxon_name, what = label_item)
  }

  ## 2.3/ Stem of the target clade ####

  # One logical value per taxon (NA = use the global argument of graft_missing_taxa())
  include_stem <- NULL
  if ("include_target_clade_stem" %in% names(data))
  {
    include_stem <- vapply(seq_len(n_taxa), function (i) parse_logical(data$include_target_clade_stem[i], label_taxon(i), "include_target_clade_stem"), logical(1))
  }

  ## 2.4/ Excluded subclades and their stem ####

  # 'excluded' will be a list with one element per taxon: a list of vectors of tips (one per excluded subclade), or NULL for a taxon with no excluded subclade.
  # 'exclude_stem' will be a list with one element per taxon: a logical vector (one value per excluded subclade), or NULL.
  excluded <- NULL
  exclude_stem <- NULL
  if ("Excluded_subclades" %in% names(data))
  {
    excluded <- vector("list", n_taxa)
    # Loop on the taxa that have excluded subclades
    for (i in which(!is.na(data$Excluded_subclades)))
    {
      # Split the cell into subclades (separated by ";"), then each subclade into tips (separated by ",")
      subclades <- trimws(strsplit(data$Excluded_subclades[i], ";", fixed = TRUE)[[1]])
      if (any(subclades == "")) { stop("'Excluded_subclades' of ", label_taxon(i), " has an empty subclade (check the ';' separators): ", data$Excluded_subclades[i]) }
      subclades_tips <- lapply(subclades, function (cell) parse_tips(cell, label_taxon(i), "Excluded_subclades"))

      # Optional names of the subclades: one name per subclade
      if ("Excluded_subclades_names" %in% names(data) && !is.na(data$Excluded_subclades_names[i]))
      {
        subclade_names <- trimws(strsplit(data$Excluded_subclades_names[i], ";", fixed = TRUE)[[1]])
        if (length(subclade_names) != length(subclades_tips) || any(subclade_names == ""))
        {
          stop("'Excluded_subclades_names' of ", label_taxon(i), " has ", length(subclade_names), " name(s) but the taxon has ", length(subclades_tips),
               " excluded subclade(s): provide one non-empty name per subclade, separated by ';'.\n",
               "  Excluded_subclades: ", data$Excluded_subclades[i], "\n  Excluded_subclades_names: ", data$Excluded_subclades_names[i])
        }
        if (anyDuplicated(subclade_names))
        {
          stop("'Excluded_subclades_names' of ", label_taxon(i), " has repeated names (", paste(unique(subclade_names[duplicated(subclade_names)]), collapse = ", "),
               "): different subclades of a taxon must have different names.")
        }
        names(subclades_tips) <- subclade_names
      }
      excluded[[i]] <- subclades_tips
    }
    # Taxa with names but without subclades
    if ("Excluded_subclades_names" %in% names(data))
    {
      orphans <- which(is.na(data$Excluded_subclades) & !is.na(data$Excluded_subclades_names))
      if (length(orphans) > 0) { stop("'Excluded_subclades_names' is provided but 'Excluded_subclades' is empty for: ", paste(label_taxon(orphans), collapse = ", ")) }
    }

    # Tips of the excluded subclades must be in the backbone
    check_tips_in_phy(lapply(excluded, function (x) unique(unlist(x))), phy, taxon_name, "Excluded_subclades")

    # Labels of the excluded subclades: the same clade (same MRCA node) must have the same label, and different clades different labels
    if ("Excluded_subclades_names" %in% names(data))
    {
      flat_tips <- unlist(lapply(seq_len(n_taxa), function (i) lapply(excluded[[i]], identity)), recursive = FALSE)
      flat_labels <- unlist(lapply(excluded, names))
      flat_taxon <- rep(taxon_name, vapply(excluded, length, integer(1)))
      if (length(flat_labels) > 0)
      {
        flat_nodes <- vapply(flat_tips, function (tips) get_clade_node(phy, tips), numeric(1))
        check_clade_labels(labels = flat_labels, nodes = flat_nodes, tips = unname(flat_tips), taxon_name = flat_taxon, what = "Excluded_subclades_names")
      }
    }

    # Stem of the excluded subclades: one value per subclade, separated by ";", or a single value used for all the subclades of the taxon
    if ("exclude_subclades_stem" %in% names(data))
    {
      exclude_stem <- vector("list", n_taxa)
      no_subclades <- character(0)
      for (i in which(!is.na(data$exclude_subclades_stem)))
      {
        if (length(excluded[[i]]) == 0) { no_subclades <- c(no_subclades, taxon_name[i]); next }
        values <- trimws(strsplit(data$exclude_subclades_stem[i], ";", fixed = TRUE)[[1]])
        if (!(length(values) %in% c(1, length(excluded[[i]]))))
        {
          stop("'exclude_subclades_stem' of ", label_taxon(i), " has ", length(values), " value(s) but the taxon has ", length(excluded[[i]]),
               " excluded subclade(s): provide one value per subclade, separated by ';', or a single value for all of them.\n",
               "  Excluded_subclades: ", data$Excluded_subclades[i], "\n  exclude_subclades_stem: ", data$exclude_subclades_stem[i])
        }
        parsed <- vapply(values, function (v) parse_logical(v, label_taxon(i), "exclude_subclades_stem"), logical(1))
        exclude_stem[[i]] <- stats::setNames(unname(parsed), if (length(parsed) == length(excluded[[i]])) names(excluded[[i]]) else NULL)
      }
      if (length(no_subclades) > 0)
      {
        warning("'exclude_subclades_stem' is provided for taxa without excluded subclades: ignored for ", paste(no_subclades, collapse = ", "), ".", call. = FALSE)
      }
    }
  }

  ## 2.5/ Position of the excluded subclades relative to the target clade ####

  # Only for target clades: the excluded subclades must be nested within the target clade of their taxon.
  # The check uses the same function as check_grafting_constraints(), on a temporary list ('mtl_check') with the same items as a 'missing_taxa_list'
  # to ensure that the position of the excluded subclades is valid relative to their target clade.
  if (!is.null(excluded) && has_target)
  {
    # The stem of a subclade is excluded unless set to FALSE (the global default of graft_missing_taxa())
    mtl_check <- list(Taxon_name = taxon_name, Target_clade = clade_tips, Excluded_subclades = excluded,
                      exclude_subclades_stem = lapply(seq_len(n_taxa), function (i)
                      {
                        value <- exclude_stem[[i]]
                        if (is.null(value)) TRUE else { value[is.na(value)] <- TRUE; value }
                      }))
    problems <- check_excluded_subclades_position(phy = phy, missing_taxa_list = mtl_check)
    format_problems <- function (df)
    {
      n_show <- min(nrow(df), 10)
      paste0(paste0("  - ", df$Taxon_name[seq_len(n_show)], " (target clade tips: ", df$Target_clade[seq_len(n_show)],
                    "; excluded subclade: ", df$Excluded_subclade[seq_len(n_show)], ")", collapse = "\n"),
             if (nrow(df) > n_show) paste0("\n  ... and ", nrow(df) - n_show, " more") else "")
    }
    wrong <- problems[problems$relation %in% c("identical", "contains_target"), ]
    if (nrow(wrong) > 0)
    {
      stop(nrow(wrong), " excluded subclade(s) are identical to, or contain, the target clade of their taxon: no branch would remain for grafting.\n",
           format_problems(wrong), "\nExcluded subclades must be nested within their target clade. Check the tips in 'Target_clade' and 'Excluded_subclades'.")
    }
    outside <- problems[problems$relation == "disjoint", ]
    if (nrow(outside) > 0)
    {
      warning(nrow(outside), " excluded subclade(s) are outside the target clade of their taxon: these exclusions have no effect.\n",
              format_problems(outside), "\nCheck the tips in 'Target_clade' and 'Excluded_subclades'.", call. = FALSE)
    }
  }

  ## 2.6/ Allowed states ####

  # 'allowed_states' will be a list with one element per taxon: a vector of states, or NULL for a taxon without taxon-level value (the global argument is then used).
  # The names of the states cannot be checked here (the simmap is not provided): the grafting functions check them.
  allowed_states <- NULL
  if ("allowed_states" %in% names(data))
  {
    allowed_states <- vector("list", n_taxa)
    for (i in which(!is.na(data$allowed_states)))
    {
      # Several traits ("A1,A2; ;C3" or "color = A1,A2 ; elevation = C3"): a list with the states of each trait
      allowed_states[[i]] <- if (grepl("[;=]", data$allowed_states[i])) parse_trait_states_cell(data$allowed_states[i], paste0("'allowed_states' of ", label_taxon(i))) else parse_states(data$allowed_states[i], label_taxon(i), "allowed_states")
    }
  }

  ## 2.7/ Build the list, with all the items named after the taxa ####

  # Every item is named after the taxa (same names, same order), as required by the grafting functions.
  # Items absent from the table are not added. Items that may contain NULL (excluded subclades and their stem) are assigned with 'list()'
  # so the NULL elements (taxa without excluded subclades, or without allowed states) are kept in the list.
  missing_taxa_list <- list(Taxon_name = taxon_name)
  missing_taxa_list[[clade_item]] <- stats::setNames(clade_tips, taxon_name)
  if (!is.null(clade_labels)) { missing_taxa_list[[label_item]] <- stats::setNames(clade_labels, taxon_name) }
  if (!is.null(include_stem)) { missing_taxa_list$include_target_clade_stem <- stats::setNames(include_stem, taxon_name) }
  if (!is.null(excluded)) { missing_taxa_list["Excluded_subclades"] <- list(stats::setNames(excluded, taxon_name)) }
  if (!is.null(exclude_stem)) { missing_taxa_list["exclude_subclades_stem"] <- list(stats::setNames(exclude_stem, taxon_name)) }
  if (!is.null(allowed_states)) { missing_taxa_list["allowed_states"] <- list(stats::setNames(allowed_states, taxon_name)) }

  # Final check of the structure, as done by the grafting functions
  validate_missing_taxa_list(missing_taxa_list, clade_item = clade_item)

  return(missing_taxa_list)
}


### 2/ Internal helpers ####

#' @title Split the tips of a cell
#'
#' @description Split a text with tips separated by commas, remove the spaces around them, and check that there is at least one tip and no repeated tip.
#'
#' @param cell Character string. Tips separated by commas.
#' @param taxon_label Character string. Taxon, for the messages.
#' @param column Character string. Name of the column, for the messages.
#'
#' @noRd
#'
#' @return Character vector of tips.
#'

parse_tips <- function (cell, taxon_label, column)
{
  tips <- trimws(strsplit(cell, ",", fixed = TRUE)[[1]])
  if (length(tips) == 0 || any(tips == ""))
  {
    stop("'", column, "' of ", taxon_label, " has an empty tip (check the ',' separators): ", cell)
  }
  if (anyDuplicated(tips))
  {
    stop("'", column, "' of ", taxon_label, " has repeated tips (", paste(unique(tips[duplicated(tips)]), collapse = ", "), "): ", cell)
  }
  return(tips)
}


#' @title Split the states of a cell
#'
#' @description Split a text with states separated by commas, remove the spaces around them, and check that there is at least one state and no repeated state.
#'
#' @param cell Character string. States separated by commas.
#' @param taxon_label Character string. Taxon, for the messages.
#' @param column Character string. Name of the column, for the messages.
#'
#' @noRd
#'
#' @return Character vector of states.
#'

parse_states <- function (cell, taxon_label, column)
{
  states <- trimws(strsplit(cell, ",", fixed = TRUE)[[1]])
  if (length(states) == 0 || any(states == ""))
  {
    stop("'", column, "' of ", taxon_label, " has an empty state (check the ',' separators): ", cell)
  }
  if (anyDuplicated(states))
  {
    stop("'", column, "' of ", taxon_label, " has repeated states (", paste(unique(states[duplicated(states)]), collapse = ", "), "): ", cell)
  }
  return(states)
}


#' @title Check that the tips of clades are in the backbone
#'
#' @param tips_list List of character vectors of tips (one per taxon, can be `NULL`).
#' @param phy Phylogeny of class `"phylo"`.
#' @param taxon_name Character vector. Name of the taxa.
#' @param column Character string. Name of the column, for the messages.
#'
#' @noRd
#'
#' @return Nothing. Stops if tips are missing from `phy`.
#'

check_tips_in_phy <- function (tips_list, phy, taxon_name, column)
{
  absent <- lapply(tips_list, function (tips) setdiff(tips, phy$tip.label))
  has_absent <- which(vapply(absent, length, integer(1)) > 0)
  if (length(has_absent) > 0)
  {
    n_show <- min(length(has_absent), 10)
    stop("Some tips of '", column, "' are missing from the backbone phylogeny (check the spelling of the tips, or the separators):\n",
         paste0("  - ", taxon_name[has_absent[seq_len(n_show)]], " (row ", has_absent[seq_len(n_show)], "): ",
                vapply(absent[has_absent[seq_len(n_show)]], paste, character(1), collapse = ", "), collapse = "\n"),
         if (length(has_absent) > n_show) paste0("\n  ... and ", length(has_absent) - n_show, " more") else "")
  }
  invisible(NULL)
}


#' @title Parse a logical value written in a cell
#'
#' @param value Character string. `TRUE`, `FALSE`, `T`, `F`, `yes`, `no`, `1`, `0`, or `NA`.
#' @param taxon_label Character string. Taxon, for the messages.
#' @param column Character string. Name of the column, for the messages.
#'
#' @noRd
#'
#' @return A logical value.
#'

parse_logical <- function (value, taxon_label, column)
{
  if (is.na(value)) { return(NA) }
  v <- toupper(trimws(value))
  if (v %in% c("TRUE", "T", "YES", "Y", "1")) { return(TRUE) }
  if (v %in% c("FALSE", "F", "NO", "N", "0")) { return(FALSE) }
  if (v == "NA") { return(NA) }
  stop("'", column, "' of ", taxon_label, " must be TRUE, FALSE or NA, not '", value, "'")
}


#' @title Check that clade labels and clades match one-to-one
#'
#' @description The same clade (same MRCA node) must always have the same label, and different clades must have different labels.
#'   Clades without a label (`NA`) are not checked.
#'
#' @param labels Character vector. Label of each clade (`NA` if none).
#' @param nodes Numeric vector. MRCA node of each clade in the backbone.
#' @param tips List of character vectors. Tips defining each clade.
#' @param taxon_name Character vector. Taxon of each clade, for the messages.
#' @param what Character string. Name of the column of the labels, for the messages.
#'
#' @noRd
#'
#' @return Nothing. Stops if labels and clades are inconsistent.
#'

check_clade_labels <- function (labels, nodes, tips, taxon_name, what)
{
  keep <- !is.na(labels)
  if (!any(keep)) { return(invisible(NULL)) }
  labels <- labels[keep]; nodes <- nodes[keep]; tips <- tips[keep]; taxon_name <- taxon_name[keep]
  tips_text <- vapply(tips, paste, character(1), collapse = ", ")
  problems <- character(0)

  # A label used for several clades
  for (lab in unique(labels))
  {
    idx <- which(labels == lab)
    if (length(unique(nodes[idx])) > 1)
    {
      per_clade <- vapply(unique(nodes[idx]), function (nd)
      {
        k <- idx[nodes[idx] == nd]
        paste0("tips (", paste(unique(tips_text[k]), collapse = " | "), ") in ", paste(unique(taxon_name[k]), collapse = ", "))
      }, character(1))
      problems <- c(problems, paste0("  - the label '", lab, "' is used for ", length(per_clade), " different clades: ", paste(per_clade, collapse = "; ")))
    }
  }
  # A clade with several labels
  for (nd in unique(nodes))
  {
    idx <- which(nodes == nd)
    if (length(unique(labels[idx])) > 1)
    {
      problems <- c(problems, paste0("  - the clade with tips (", paste(unique(tips_text[idx]), collapse = " | "), ") has ",
                                     length(unique(labels[idx])), " different labels: ", paste0("'", unique(labels[idx]), "'", collapse = ", "),
                                     " (in ", paste(unique(taxon_name[idx]), collapse = ", "), ")"))
    }
  }

  if (length(problems) > 0)
  {
    n_show <- min(length(problems), 10)
    stop("Inconsistent labels in '", what, "': ", length(problems), " problem(s).\n",
         "The same clade must always have the same label, and different clades must have different labels. ",
         "A clade is identified by the MRCA node of its tips in the backbone phylogeny.\n",
         paste(problems[seq_len(n_show)], collapse = "\n"),
         if (length(problems) > n_show) paste0("\n  ... and ", length(problems) - n_show, " more") else "",
         "\nCheck the spelling of the labels and the tips of these clades.", call. = FALSE)
  }
  invisible(NULL)
}
