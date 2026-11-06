##### Function to import a missing_taxa_list from an external table #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Import a table (.csv, .tsv, .xlsx, ...) as a missing_taxa_list ####

#' @title Import an external table of missing taxa as a `missing_taxa_list` object
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
#' * `allowed_states`: (Optional) For each taxon, the states of a `simmap` in which it can be grafted (see the `simmap` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#'   or [`phylograftR::graft_missing_clades()`][graft_missing_clades]), separated by commas. It is accepted with both `Target_clade` and `Sister_clade`.
#'   Ex: `"blue, green"`. Blank cells and `NA` have no taxon-level value, so the global `allowed_states` argument of the grafting function is used for the taxon.
#'   States are matched exactly with the states of the `simmap` (e.g., case matters) when grafting, and so cannot contain commas.
#'   If the `simmap` argument is a list of the maps of several traits, give the states of each trait separated by semicolons, in the order of the traits (e.g., `"A1, A2; ; C3"` for
#'   three traits: no states for the second trait, so the global `allowed_states` of the trait, or all its states, are used), or name the traits (e.g., `"color = A1, A2 ; elevation = C3"`).
#'   The names of the states are not checked here as the `simmap` is not provided: the grafting function checks them.
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
#' plot_grafted(grafted_phy, cex = 0.7)
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
