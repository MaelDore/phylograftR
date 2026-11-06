##### Function to import taxonomic tables from external files #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Import a taxonomic table (.csv, .tsv, .xlsx, ...) as a data frame ####

#' @title Import an external taxonomic table as a data frame
#'
#' @description Read a taxonomic table from a file (e.g., a .csv file, or an Excel sheet), check its format, and return it as a data frame
#'   ready-to-be-used as input of [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy].
#'
#'   It is the mirror of [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] for taxonomic tables.
#'   Two tables are used to build a `missing_taxa_list` from taxonomy: the taxonomy of the tips of the backbone phylogeny (`type = "backbone"`),
#'   and the description of the taxa to graft (`type = "missing_taxa"`). Each one is read with a call of the function.
#'
#'   Taxonomic tables are the go-to option for large phylogenies with neatly defined taxonomic levels (see Details).
#'
#' @param file Character string. Path to the file with the table.
#'   The format is chosen from the extension of the file: `".csv"` (comma-separated), `".tsv"`, `".tab"` or `".txt"` (tab-separated),
#'   `".xlsx"` or `".xls"` (Excel, requires the package `readxl`). Use `sep` or `reader` for other formats.
#'   All the columns are read as text. Blank cells and `NA` are missing values.
#' @param type Character string. The type of table to read:
#'   * `"backbone"` for the taxonomy of the tips of the backbone phylogeny (one row per tip),
#'   * `"missing_taxa"` for the description of the taxa to graft (one row per taxon).
#'   See Details for the columns of each one.
#' @param phy An object of class `"phylo"`, or `NULL`. Optional backbone phylogeny, used to check the table against its tips:
#'   all the tips of a `"backbone"` table must be tips of `phy`, and no taxon of a `"missing_taxa"` table can already be a tip of `phy`.
#'   Default is `NULL`, with no check against the phylogeny (it is anyway done by [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy]).
#' @param sep Character string. Field separator of a text file. Default is `NULL`, to choose it from the extension of the file:
#'   `","` for `.csv`, and a tab for `.tsv`, `.tab` and `.txt`. Semicolons are also used inside the cells of the table of the missing taxa (to separate the excluded subclades, and their stems),
#'   so do not use them as field separator, and quote the cells that contain commas (e.g., `allowed_states`).
#' @param sheet Character string or integer. Only used for Excel files. Name or position of the sheet to read. Default is `1`.
#'   Both taxonomic tables can be stored in two sheets of the same Excel file, and read with two calls.
#' @param reader Function. Optional function used to read `file`, to import any other format (e.g., `.ods`, or an export of an online spreadsheet).
#'   It must take the path of the file as first argument and return a data frame. If provided, it is used whatever the extension of the file.
#'   Default is `NULL`.
#' @param verbose Logical. Whether to print a summary of the table. Default is `TRUE`.
#' @param ... Further arguments passed to the function reading the file (`utils::read.table()` for text files, `readxl::read_excel()` for Excel files, or `reader`).
#'
#' @details
#' **When to use taxonomic tables**
#'
#' Taxonomic tables are the go-to option for large phylogenies with neatly defined taxonomic levels (e.g., Family, Subfamily, Tribe, Genus). The taxonomy of the backbone is provided once,
#' and each taxon to graft only needs its own taxonomic levels: the constraints (target clade, sister clade, excluded subclades) are written with the names of the groups or of the taxonomic levels to use,
#' and the clades are retrieved from the phylogeny. The same group is therefore always the same clade, whatever the number of taxa, which avoids the typos and
#' inconsistencies of clades written taxon by taxon, and the tables can be taken as such from a taxonomic database.
#'
#' The limit is that only the taxonomic groups of the backbone table can be used as clades. To fine-tune the grafting of some taxa with any other node of the phylogeny
#' (e.g., a clade without a name, a part of a group, or a group of several genera), define their clades by tips or node IDs with
#' [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] or [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] (see 'Fine-tuning the clades' in their Details sections).
#' Both routes return the same `missing_taxa_list`, so a list built from taxonomic tables can be adjusted afterwards.#'
#'
#' **Backbone table** (`type = "backbone"`)
#'
#' One row per tip of the backbone phylogeny. The column `Taxa` (Required) holds the tip labels, and must have unique values.
#' All the other columns are taxonomic levels (e.g., `Family`, `Subfamily`, `Tribe`, `Genus`, but they can be anything), usually from the broadest to the narrowest.
#' At least one level is required. The values of a level are the names of the taxonomic groups (e.g., `"Homo"` for a `"Genus"`).
#' Levels do not have to be informed for all the tips: blank cells and `NA` are missing values.
#'
#' **Missing taxa table** (`type = "missing_taxa"`)
#'
#' One row per taxon to graft. The columns `Taxa` (Required; unique names of the taxa to graft) and `Type` (Required; `"taxa"` or `"clades"`) are required.
#' The next columns are the taxonomic levels of the taxa to graft, with the same names as in the backbone table, and the addtional constraints:
#' * `Target_clade` (Required for `Type = "taxa"`) or `Sister_clade` (Required for `Type = "clades"`): A single taxonomic group, written either as the name of the group
#'   (e.g., `"Homo"`), or as the name of a taxonomic level (e.g., `"Genus"`). If a taxonomic level is provided, the group will be identified as recorded in the column of the same name
#'   in `missing_taxa_taxonomic_table` (e.g., `"Genus"` with `Genus = "Camponotus"` gives the clade of the genus `Camponotus`).
#'   It is an error if the group is absent from the backbone, or if the taxon has no value at the requested level.
#' * `Excluded_subclades` (Optional, `Type = "taxa"`): Subclades of the target clade where the taxon is known not to belong, separated by semicolons (e.g., `"Pheidole; Camponotus"`).
#'   Each one is a group or a level name, as in `Target_clade`.
#' * `$include_target_clade_stem` (Optional, `Type = "taxa"`): Logical. Whether the stem branch of the target clade is a valid placement to graft this taxa (`NA` or empty cell uses the global argument set in the grafting function).
#' * `$exclude_subclades_stem` (Optional, `Type = "taxa"`): Logical. Whether the stem branch of the excluded clade(s) is/are a excluded from grafting too. Separated by semicolons if multiple `Excluded_subclades` are provided
#'   (`NA` or empty cell uses the global argument set in the grafting function).
#' * `$allowed_states` (Optional): For each taxon, the states of a `simmap` in which it can be grafted
#'   (see the `simmap` argument of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades]), separated by commas.
#'   If multiple traits are used, each set of states must be separated by semicolons (e.g., `"blue, red; low"`).
#'   It is accepted with both `Target_clade` and `Sister_clade`. (`NA` or empty cell uses the global argument set in the grafting function).
#'
#' **Checks**
#'
#' The function checks that the required columns are present, that names of columns are unique, that `Taxa` is never empty and always unique,
#' and, for the table of the missing taxa, that `Type` is valid. Names of columns and values are stripped from the spaces around them.
#' Rows with only blank cells (e.g., the trailing rows of a spreadsheet) are dropped.
#'
#' The other checks (e.g., names of levels and groups, consistency with the backbone, and the groups used as clades) are done when the clades are retrieved,
#' by [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy].
#'
#' @return A data frame with only character columns, ready to be used as `backbone_taxonomic_table` (if `type = "backbone"`)
#'   or `missing_taxa_taxonomic_table` (if `type = "missing_taxa"`) in [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy].
#'   Blank cells are `NA`.
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy] that uses the tables to build a `missing_taxa_list`,
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read a table of the missing taxa with clades defined by tips or node IDs,
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] that use the list for grafting.
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Paths to example tables, provided with the package
#' backbone_path <- system.file("extdata", "toy_backbone_taxonomic_table.csv", package = "phylograftR")
#' missing_taxa_path <- system.file("extdata", "toy_missing_taxa_taxonomic_table.csv", package = "phylograftR")
#'
#' # Load the associated backbone phylogeny
#' data(toy_backbone)
#'
#' # Import the taxonomy of the backbone: one row per tip, one column per taxonomic level
#' backbone_taxonomic_table <- import_taxonomic_table(file = backbone_path, type = "backbone", phy = toy_backbone)
#' head(backbone_taxonomic_table)
#'
#' # Import the description of the taxa to graft: one row per taxon
#' missing_taxa_taxonomic_table <- import_taxonomic_table(file = missing_taxa_path, type = "missing_taxa", phy = toy_backbone)
#' missing_taxa_taxonomic_table
#'
#' # The tables are the datasets provided with the package
#' data(toy_backbone_taxonomic_table)
#' data(toy_missing_taxa_taxonomic_table)
#' identical(backbone_taxonomic_table, toy_backbone_taxonomic_table)
#' # (the logical column include_target_clade_stem of the dataset is read as text from the file)
#' all(mapply(function (x, y) identical(as.character(x), y),
#'            toy_missing_taxa_taxonomic_table, missing_taxa_taxonomic_table))
#'
#' # Use them to build the lists of missing taxa
#' missing_taxa_lists <- build_missing_taxa_list_from_taxonomy(
#'   backbone_taxonomic_table = backbone_taxonomic_table,
#'   missing_taxa_taxonomic_table = missing_taxa_taxonomic_table,
#'   phy = toy_backbone, verbose = FALSE)
#' names(missing_taxa_lists)
#'
#' # Use them to graft the taxa
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = missing_taxa_lists$taxa,
#'                                   verbose = FALSE, seed = 1234)
#'

import_taxonomic_table <- function (file, # Path to the table
                                    type, # Table to read
                                    phy = NULL, # Optional backbone phylogeny
                                    sep = NULL, # Field separator of text files. By default, from the extension
                                    sheet = 1, # Sheet of Excel files
                                    reader = NULL, # Optional function to read other formats
                                    verbose = TRUE, # Print a summary of the table
                                    ...)
{
  ## 1.1/ Check the inputs ####

  type <- match.arg(type)
  if (!is.null(phy) && !inherits(phy, "phylo")) { stop("'phy' must be NULL or a phylogeny of class 'phylo'") }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) { stop("'verbose' must be TRUE or FALSE") }

  ## 1.2/ Read the table ####

  data <- read_table_file(file = file, sep = sep, sheet = sheet, reader = reader, ...)

  # Strip the spaces, and turn blank cells into NA. Same cleaning as when the table is provided directly to build_missing_taxa_list_from_taxonomy()
  arg <- if (type == "backbone") "backbone_taxonomic_table" else "missing_taxa_taxonomic_table"
  data <- clean_taxonomic_table(data, arg)
  # Drop the rows without any value, e.g., the trailing empty rows of a spreadsheet
  empty_rows <- rowSums(!is.na(data)) == 0
  data <- data[!empty_rows, , drop = FALSE]
  rownames(data) <- NULL

  ## 1.3/ Check the table ####

  required <- if (type == "backbone") "Taxa" else c("Taxa", "Type")
  for (column in required)
  {
    if (!(column %in% names(data))) { stop("The table read from '", file, "' must have a column '", column, "'. Columns found: ", format_names(names(data))) }
  }
  if (nrow(data) == 0) { stop("The table read from '", file, "' has no rows with values") }
  if (anyNA(data$Taxa)) { stop("'Taxa' is empty in the row(s) of '", file, "': ", paste(which(is.na(data$Taxa)), collapse = ", ")) }
  if (anyDuplicated(data$Taxa)) { stop("'Taxa' must be unique in '", file, "'. Duplicated: ", format_names(unique(data$Taxa[duplicated(data$Taxa)]))) }

  if (type == "backbone")
  {
    if (ncol(data) < 2) { stop("The backbone table read from '", file, "' must have at least one taxonomic level (a column other than 'Taxa')") }
    if (!is.null(phy))
    {
      not_in_phy <- setdiff(data$Taxa, phy$tip.label)
      if (length(not_in_phy) > 0) { stop(length(not_in_phy), " taxa of the backbone table are not tips of 'phy': ", format_names(not_in_phy)) }
    }
  } else {
    invalid_type <- is.na(data$Type) | !(tolower(data$Type) %in% c("taxa", "clades"))
    if (any(invalid_type))
    {
      stop("'Type' must be 'taxa' or 'clades' (case is ignored). Invalid for: ", format_names(data$Taxa[invalid_type]))
    }
    if (!is.null(phy))
    {
      in_phy <- intersect(data$Taxa, phy$tip.label)
      if (length(in_phy) > 0) { stop("Taxa to graft must be absent from the backbone, but are tips of 'phy': ", format_names(in_phy)) }
    }
  }

  ## 1.4/ Summarize it ####

  if (verbose)
  {
    levels_found <- setdiff(names(data), required)
    if (type == "backbone")
    {
      message("Backbone taxonomic table: ", nrow(data), " tips, with ", length(levels_found), " taxonomic level(s): ", paste(levels_found, collapse = ", "))
    } else {
      message("Missing taxa taxonomic table: ", nrow(data), " taxa to graft (", sum(tolower(data$Type) == "taxa"), " as 'taxa', ",
              sum(tolower(data$Type) == "clades"), " as 'clades'), with ", length(levels_found), " other column(s): ", paste(levels_found, collapse = ", "))
    }
  }

  return(data)
}
