##### Functions to import and build a missing_taxa_list from a table #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Import a table (.csv, .tsv, .xlsx, ...) as a missing_taxa_list ####

#' @title Import a table of missing taxa as a `missing_taxa_list` object
#'
#' @description Read a table (e.g., a .csv file, or an Excel sheet) describing the taxa to graft, and convert it into a valid `missing_taxa_list`
#'   object ready-to-be-used as input for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'
#'   The table is converted by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list], which also runs all the checks of validity
#'   against the backbone phylogeny to ensure the `missing_taxa_list` provides information that are valid to process the grafting.
#'
#' @param file Character string. Path to the file with the table.
#'   The format is chosen from the extension of the file: `".csv"` (comma-separated), `".tsv"`, `".tab"` or `".txt"` (tab-separated),
#'   `".xlsx"` or `".xls"` (Excel, requires the package `readxl`). Use `sep` or `reader` for other formats.
#'
#'   The table must record one row per taxon to graft.
#'   Its columns are the items of `missing_taxa_list` as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'   See Details for a full description. All the columns are read as text.
#'   Blank cells and `NA` are missing values.
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which the taxa will be grafted. It is used to check the validity of the tips and clades provided.
#' @param sep Character string. Field separator of a text file. Default is `NULL`, to choose it from the extension of the file:
#'   `","` for `.csv`, and a tab for `.tsv`, `.tab` and `.txt`. Commas and semicolons are also used inside the cells, to separate the tips and the subclades:
#'   to distinguish from the field separator, the cells containing tips and the subclades data must be quoted.
#' @param sheet Character string or integer. Only used for Excel files. Name or position of the sheet to read. Default is `1`.
#' @param reader Function. Optional function used to read `file`, to import any other format (e.g., `.ods`, or an export of an online spreadsheet).
#'   It must take the path of the file as first argument and return a data frame. If provided, it is used whatever the extension of the file.
#'   Default is `NULL`.
#' @param verbose Logical. Whether to print the conversions of node IDs into tips (if any). Default is `TRUE`.
#' @param ... Further arguments passed to the function reading the file (`utils::read.table()` for text files, `readxl::read_excel()` for Excel files, or `reader`).
#'
#' @details
#' The columns of `file` are the items of `missing_taxa_list` as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#' Column names must be written exactly as below. Tips are separated by commas (`,`), and the subclades by semicolons (`;`). Spaces around separators are ignored.
#' Values are `TRUE` or `FALSE` (also accepts `T`, `F`, `yes`, `no`, `1`, `0`), or `NA`.
#'
#' * `Taxon_name` (Required): Names of the taxa to graft. Names must be unique and absent from `phy$tip.label`.
#' * `Target_clade` or `Sister_clade` (Required): Tip labels defining a target clade per taxon, i.e., the clade where the taxon is known to belong (e.g., its genus);
#'   or a sister clade, i.e., the clade on whose stem branch the taxon is grafted, as its sister.
#'   A pair of tips defines the clade as their MRCA node (e.g., `"t26, t18"`).
#'   A single tip defines the "clade" as the terminal branch of this tip (e.g., `"t15"`).
#'   A node ID of `phy` (digits only, e.g., `"34"`) also defines the clade, see 'Node IDs' below.
#'   Use `Sister_clade` to prepare a list for [`phylograftR::graft_missing_clades()`][graft_missing_clades]:
#'   the columns specific to [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] below (`Target_clade_name`, `include_target_clade_stem`, `Excluded_subclades`, `Excluded_subclades_names` and `exclude_subclades_stem`) are ignored, with a warning.
#' * `include_target_clade_stem`: (Optional) Logical values. For each taxon, whether the stem branch of its target clade is available for grafting.
#'   `NA` values and blank cells use the global `include_target_clade_stem` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#' * `Excluded_subclades`: (Optional) For each taxon, tip labels defining the subclades of its target clade where the taxon is known not to belong.
#'   Each subclade is defined by a pair of tips, a single tip, or a node ID.
#'   Ex: `"t10, t19; 50; t7"` for three subclades, the two first being clades with multiple species defined by a pair of tips (`"t10, t19"`), and an internal node ID (`"50"`).
#'   The last one is a single tip define by its tip label (`"t7"`). Use blank cells for taxa without excluded subclades.
#' * `exclude_subclades_stem`: (Optional) Logical values. For each taxon, whether the stem of each of its excluded subclades is also excluded,
#'   with one value per subclade of `Excluded_subclades` separated by semicolons, or a single value used for all of them.
#'   Ex: `"FALSE; TRUE; TRUE"`.
#'   `NA` values and blank cells use the global `exclude_subclades_stem` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#'
#' Optional columns that do not match items in `missing_taxa_list` can be used to provide names of clades:
#'
#' * `Target_clade_name` or `Sister_clade_name` (Optional): Label of the target (or sister) clade.
#'   It is only used to check that the clades are defined consistently (see [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]).
#' * `Excluded_subclades_names` (Optional): Labels of the excluded subclades, separated by semicolons, in the same order as in `Excluded_subclades`.
#'   Ex: `"group_1; group_2; tip_t7"`. There must be one label per subclade.
#'   If provided, they are used as names of the subclades in the `Excluded_subclades` item of the list, and in the messages of the grafting functions.
#'
#' Other columns are ignored, with a warning.
#'
#' ### Node IDs
#'
#' Clades (`Target_clade`, `Sister_clade`, and each subclade of `Excluded_subclades`) are typically defined by (pairs of) tips.
#' However, the function also accepts clades defined by the ID of a node of `phy` instead of tips, in the same table, and even in the same cell
#' (e.g., `"50; t5, t12; t7"` for three excluded subclades defined respectively by an internal node, a pair of tips, a terminal tip label).
#' Node IDs are integers: `1` to N, with N being the number of tips (`ape::Ntip(phy)`), correspond to terminal nodes, and thus define single-tip clades.
#' Following IDs are the internal nodes, as shown by `ape::nodelabels()`, and correspond to clades with multiple descending species.
#'
#' A value made only of digits is a node ID, so **tip labels must include at least one other character than digits** to be distinguished.
#' A clade is defined by tips or by one node ID, not both, and one node ID defines the clade alone.
#' Node IDs are converted into tips before any other check, as done by [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips],
#' and the conversions are printed if `verbose = TRUE`.
#'
#' The warnings and errors that follow, and the grafting functions, therefore use the tips that are traced back to your node IDs in the messages.
#' If a clade is already defined by tips in the table, a node ID of the same clade is converted into the same tips.
#' Convert `total_richness_per_clades` (if you use it) at the same time as the list, with [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips],
#' so the tips of the clades match between both objects.
#'
#' @return A `missing_taxa_list`: a list, ready to be used in the grafting functions, with the following items:
#'  * `$Taxon_name` Character vector.
#'  * `$Target_clade` (or `$Sister_clade`) List of character vectors of tips.
#'  * `$Target_clade_name` (or `$Sister_clade_name`) Only if the column was provided. Character vector of labels (`NA` if none). Ignored by the grafting functions.
#'  * `$include_target_clade_stem` Only if the column was provided. Logical vector (`NA` uses the global argument).
#'  * `$Excluded_subclades` Only if the column was provided. List (one element per taxon) of lists of character vectors of tips. `NULL` for taxa with no excluded subclades.
#'    The subclades are named after their labels if `Excluded_subclades_names` was provided.
#'  * `$exclude_subclades_stem` Only if the column was provided. List (one element per taxon) of logical vectors (one value per excluded subclade), or `NULL`.
#'
#'  Every item is named after the taxa (`$Taxon_name`), in the same order.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] to convert a data frame already loaded in R,
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] that use the list for grafting.
#'
#' @importFrom utils read.table
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Path to an example table of missing taxa, provided with the package
#' file_path <- system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")
#' # Look at the table: one row per taxon, and the fields of the items of 'missing_taxa_list'
#' read.csv(file_path)
#'
#' # Load the associated backbone phylogeny
#' data(toy_backbone)
#'
#' # Import it as a 'missing_taxa_list', checking it against the backbone phylogeny
#' missing_taxa_list <- import_missing_taxa_list(file = file_path, phy = toy_backbone)
#' str(missing_taxa_list, max.level = 1)
#'
#' # The same list is provided as a dataset
#' data(toy_missing_taxa_list)
#' identical(missing_taxa_list, toy_missing_taxa_list)
#'
#' # Use it to graft the taxa
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = missing_taxa_list,
#'                                   method = "uniform", verbose = FALSE, seed = 1234)
#'
#' # Plot the pruned phylogeny with the updated colors
#' plot(grafted_phy, cex = 0.7,
#'      edge.color = c("black", "red")[grafted_phy$grafted_edges + 1],
#'      tip.color = c("black", "red")[grafted_phy$grafted_tips + 1])
#'

import_missing_taxa_list <- function (file, # Path to the table
                                      phy, # Backbone phylogeny
                                      sep = NULL, # Field separator of text files. By default, from the extension
                                      sheet = 1, # Sheet of Excel files
                                      reader = NULL, # Optional function to read other formats
                                      verbose = TRUE, # Print the conversions of node IDs into tips
                                      ...)
{
  ## 1.1/ Check the inputs ####

  if (!is.character(file) || length(file) != 1 || is.na(file)) { stop("'file' must be a single path") }
  if (!file.exists(file)) { stop("File not found: ", file) }
  if (!is.null(reader) && !is.function(reader)) { stop("'reader' must be a function that takes the path of the file and returns a data frame") }
  if (!is.null(sep) && (!is.character(sep) || length(sep) != 1 || is.na(sep) || nchar(sep) != 1)) { stop("'sep' must be a single character") }

  ## 1.2/ Read the table ####

  extension <- tolower(sub(".*\\.", "", basename(file)))

  if (!is.null(reader))
  {
    data <- reader(file, ...)
  } else if (extension %in% c("xlsx", "xls"))
  {
    if (!requireNamespace("readxl", quietly = TRUE))
    {
      stop("The package 'readxl' is required to read Excel files. Install it with install.packages(\"readxl\"), ",
           "or export the table as a .csv file.")
    }
    data <- readxl::read_excel(path = file, sheet = sheet, col_types = "text", ...)
  } else if (extension %in% c("csv", "tsv", "tab", "txt") || !is.null(sep))
  {
    if (is.null(sep)) { sep <- if (extension == "csv") "," else "\t" }
    # Everything is read as text (the cells contain lists of tips), blank cells are missing values, and the names of the columns are kept as they are
    data <- utils::read.table(file = file, header = TRUE, sep = sep, quote = "\"", colClasses = "character",
                              na.strings = c("", "NA"), strip.white = TRUE, check.names = FALSE, comment.char = "",
                              fill = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8", ...)
  } else {
    stop("Unknown file extension '.", extension, "'. Use a '.csv', '.tsv', '.tab', '.txt', '.xlsx' or '.xls' file, ",
         "or provide 'sep' (for a text file) or 'reader' (a function that reads the file).")
  }

  if (!is.data.frame(data)) { stop("The table read from '", file, "' is not a data frame") }
  if (nrow(data) == 0) { stop("The table read from '", file, "' has no rows") }

  ## 1.3/ Convert it ####

  build_missing_taxa_list(data = data, phy = phy, verbose = verbose)
}


### 2/ Build a missing_taxa_list from a data frame ####

#' @title Build a `missing_taxa_list` from a table of missing taxa
#'
#' @description Build a valid `missing_taxa_list` from a data frame describing the taxa to graft (one row per taxon, e.g., read from a .csv file),
#'   for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'   Lists of tips (or node IDs, see Details) are written in the cells of the table, and the function parses them and names every item after the taxa (as required by the grafting functions),
#'   and runs all the checks of validity against the backbone phylogeny.
#'
#'   Use [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read the table directly from a file.
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
#'  * An excluded subclade is identical to, or contains, the target clade of its taxon (no branch would remain).
#'    A warning is issued if it is outside of the target clade (the exclusion has no effect).
#'
#' The messages give the names of the problematic taxa and the numbers of the rows in the table.
#'
#' @return A `missing_taxa_list`: a list, ready to be used in the grafting functions, with the following items:
#'  * `$Taxon_name` Character vector.
#'  * `$Target_clade` (or `$Sister_clade`) List of character vectors of tips.
#'  * `$Target_clade_name` (or `$Sister_clade_name`) Only if the column was provided. Character vector of labels (`NA` if none). Ignored by the grafting functions.
#'  * `$include_target_clade_stem` Only if the column was provided. Logical vector (`NA` uses the global argument).
#'  * `$Excluded_subclades` Only if the column was provided. List (one element per taxon) of lists of character vectors of tips. `NULL` for taxa with no excluded subclades.
#'    The subclades are named after their labels if `Excluded_subclades_names` was provided.
#'  * `$exclude_subclades_stem` Only if the column was provided. List (one element per taxon) of logical vectors (one value per excluded subclade), or `NULL`.
#'
#'  Every item is named after the taxa (`$Taxon_name`), in the same order.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read the table from a file,
#'  * [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] to convert the node IDs of a list built 'by hand',
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

  ## 2.1.1/ Columns ####

  columns_taxa <- c("Taxon_name", "Target_clade", "Target_clade_name", "include_target_clade_stem",
                    "Excluded_subclades", "Excluded_subclades_names", "exclude_subclades_stem")
  columns_clades <- c("Taxon_name", "Sister_clade", "Sister_clade_name")
  all_columns <- unique(c(columns_taxa, columns_clades))

  # Names of the columns without extra spaces. All the cells are read as trimmed text. Blank cells are missing values
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  names(data) <- trimws(names(data))
  if (anyDuplicated(names(data))) { stop("Some columns have the same name: ", paste(unique(names(data)[duplicated(names(data))]), collapse = ", ")) }
  data[] <- lapply(data, function (x) { x <- trimws(as.character(x)); x[!is.na(x) & x %in% c("", "NA")] <- NA; x })

  if (!("Taxon_name" %in% names(data))) { stop("The table must have a column 'Taxon_name'") }

  has_target <- "Target_clade" %in% names(data)
  has_sister <- "Sister_clade" %in% names(data)
  if (has_target && has_sister) { stop("The table has both 'Target_clade' and 'Sister_clade': use only one, depending on the grafting function (graft_missing_taxa() or graft_missing_clades())") }
  if (!has_target && !has_sister) { stop("The table must have a column 'Target_clade' (for graft_missing_taxa()) or 'Sister_clade' (for graft_missing_clades())") }
  clade_item <- if (has_target) "Target_clade" else "Sister_clade"
  label_item <- paste0(clade_item, "_name")
  allowed_columns <- if (has_target) columns_taxa else columns_clades

  # Columns not recognized (typos) or not used with the clade column provided
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

  ## 2.1.2/ Taxa ####

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

  # Clades defined by a node ID (digits only) are converted into tips, before any check based on tips.
  # Clades defined by tips are left unchanged. Empty tips or subclades are left to the checks of the tips below.

  # Use an internal environment to record conversion such as a node ID is converted into the tips already used for the same clade
  recorder <- new.env(parent = emptyenv())
  recorder$rows <- list()

  # Function to parse cell information
  split_cell <- function (cell, split_subclades)
  {
    if (is.na(cell)) { return(list()) }
    pieces <- if (split_subclades) trimws(strsplit(cell, ";", fixed = TRUE)[[1]]) else cell
    lapply(pieces, function (piece) trimws(strsplit(piece, ",", fixed = TRUE)[[1]]))
  }
  known <- unlist(lapply(seq_len(n_taxa), function (i) split_cell(data[[clade_item]][i], FALSE)), recursive = FALSE)
  if ("Excluded_subclades" %in% names(data))
  {
    known <- c(known, unlist(lapply(seq_len(n_taxa), function (i) split_cell(data$Excluded_subclades[i], TRUE)), recursive = FALSE))
  }
  lookup <- build_tips_lookup(phy, known)
  convert_cell <- function (cell, item, who, split_subclades)
  {
    if (is.na(cell)) { return(cell) }
    pieces <- if (split_subclades) trimws(strsplit(cell, ";", fixed = TRUE)[[1]]) else cell
    converted <- vapply(pieces, function (piece)
    {
      tokens <- trimws(strsplit(piece, ",", fixed = TRUE)[[1]])
      if (length(tokens) == 0 || any(tokens == "")) { return(piece) }
      paste(resolve_clade_tokens(tokens, phy = phy, item = item, who = who, recorder = recorder, cell = piece, lookup = lookup), collapse = ",")
    }, character(1), USE.NAMES = FALSE)
    paste(converted, collapse = ";")
  }
  data[[clade_item]] <- vapply(seq_len(n_taxa), function (i) convert_cell(data[[clade_item]][i], clade_item, taxon_name[i], FALSE), character(1))
  if ("Excluded_subclades" %in% names(data))
  {
    data$Excluded_subclades <- vapply(seq_len(n_taxa), function (i) convert_cell(data$Excluded_subclades[i], "Excluded_subclades", taxon_name[i], TRUE), character(1))
  }
  report_node_conversions(recorder$rows, phy = phy, verbose = verbose)

  ## 2.2/ Clades of the taxa: Target_clade or Sister_clade ####

  clade_cells <- data[[clade_item]]
  if (anyNA(clade_cells)) { stop("'", clade_item, "' is empty for: ", paste(label_taxon(which(is.na(clade_cells))), collapse = ", ")) }
  clade_tips <- lapply(seq_len(n_taxa), function (i) parse_tips(clade_cells[i], label_taxon(i), clade_item))
  check_tips_in_phy(clade_tips, phy, taxon_name, clade_item)
  clade_nodes <- vapply(clade_tips, function (tips) get_clade_node(phy, tips), numeric(1))

  # Optional labels of the clades
  clade_labels <- NULL
  if (label_item %in% names(data))
  {
    clade_labels <- data[[label_item]]
    check_clade_labels(labels = clade_labels, nodes = clade_nodes, tips = clade_tips, taxon_name = taxon_name, what = label_item)
  }

  ## 2.3/ Stem of the target clade ####

  include_stem <- NULL
  if ("include_target_clade_stem" %in% names(data))
  {
    include_stem <- vapply(seq_len(n_taxa), function (i) parse_logical(data$include_target_clade_stem[i], label_taxon(i), "include_target_clade_stem"), logical(1))
  }

  ## 2.4/ Excluded subclades and their stem ####

  excluded <- NULL
  exclude_stem <- NULL
  if ("Excluded_subclades" %in% names(data))
  {
    excluded <- vector("list", n_taxa)
    for (i in which(!is.na(data$Excluded_subclades)))
    {
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

    # Tips in the backbone
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

    # Stem of the excluded subclades
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

  ## 2.6/ Build the list, with all the items named after the taxa ####

  missing_taxa_list <- list(Taxon_name = taxon_name)
  missing_taxa_list[[clade_item]] <- stats::setNames(clade_tips, taxon_name)
  if (!is.null(clade_labels)) { missing_taxa_list[[label_item]] <- stats::setNames(clade_labels, taxon_name) }
  if (!is.null(include_stem)) { missing_taxa_list$include_target_clade_stem <- stats::setNames(include_stem, taxon_name) }
  if (!is.null(excluded)) { missing_taxa_list["Excluded_subclades"] <- list(stats::setNames(excluded, taxon_name)) }
  if (!is.null(exclude_stem)) { missing_taxa_list["exclude_subclades_stem"] <- list(stats::setNames(exclude_stem, taxon_name)) }

  # Final check of the structure, as done by the grafting functions
  validate_missing_taxa_list(missing_taxa_list, clade_item = clade_item)

  return(missing_taxa_list)
}


### 3/ Internal helpers ####

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
