##### Function to import total_richness_per_clades from an external table #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Import a table (.csv, .tsv, .xlsx, ...) as total_richness_per_clades ####

#' @title Import an external table of clades as `total_richness_per_clades` objects
#'
#' @description Read a table (e.g., a .csv file, or an Excel sheet) with the total number of species per clades, and convert it into the `total_richness_per_clades` objects
#'   ready-to-be-used as input for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and/or [`phylograftR::graft_missing_clades()`][graft_missing_clades].
#'
#'   The table is converted by [`phylograftR::build_total_richness_per_clades()`][build_total_richness_per_clades], which matches the clades with those of `missing_taxa_list`,
#'   and checks that the richness is consistent with the full sequence of grafting.
#'
#' @param file Character string. Path to the file with the table.
#'   The format is chosen from the extension of the file: `".csv"` (comma-separated), `".tsv"`, `".tab"` or `".txt"` (tab-separated),
#'   `".xlsx"` or `".xls"` (Excel, requires the package `readxl`). Use `sep` or `reader` for other formats.
#'
#'   The table must record one row per clade, with the columns `Clade` (tips, node ID, or name), `total_richness`, and, optionally, `Clade_name`.
#'   * `Clade`: a monophyletic group from the backbone phylogeny defined by tips, a node ID, or a name (see Details), from which the total richness is known.
#'   * `total_richness`: Total number of species of the clade. It includes the backbone species, **all the taxa to graft that belong to the clade**, and the unsampled species that are known to exist but will not be grafted.
#'     Use `NA` if unknown.
#'   * `Clade_name` (Optional): Name of the clade, used to check matching with other inputs, and to report warnings in downstream functions. If the clade is written with a name in `Clade`, it is its name.
#'
#'   The columns `Target_clade` and `Target_clade_name`, or `Sister_clade` and `Sister_clade_name`, are accepted instead of `Clade` and `Clade_name` (e.g., for tables prepared for a single grafting function).
#'   All the columns are read as text. Blank cells and `NA` are missing values.
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which the taxa will be grafted. It is used to check the validity of the tips and clades provided.
#' @param missing_taxa_list A `missing_taxa_list`, or the list of two `missing_taxa_list` (`$taxa` and `$clades`) built for the same backbone phylogeny,
#'   e.g., with [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] or [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy].
#'   It defines the type(s) of the output, and is used to harmonize the clades and check the validity of the richness values provided.
#' @param sep Character string. Field separator of a text file. Default is `NULL`, to choose it from the extension of the file:
#'   `","` for `.csv`, and a tab for `.tsv`, `.tab` and `.txt`. Commas are also used inside the cells to separate the tips:
#'   to distinguish from the field separator, **the cells containing tips must be quoted.**
#' @param sheet Character string or integer. Only used for Excel files. Name or position of the sheet to read. Default is `1`.
#' @param reader Function. Optional function used to read `file`, to import any other format (e.g., `.ods`, or an export of an online spreadsheet).
#'   It must take the path of the file as first argument and return a data frame. If provided, it is used whatever the extension of the file.
#'   Default is `NULL`.
#' @param verbose Logical. Whether to print the conversions of node IDs into tips (if any), and a summary of the table. Default is `TRUE`.
#' @param ... Further arguments passed to the function reading the file (`utils::read.table()` for text files, `readxl::read_excel()` for Excel files, or `reader`).
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
#' The resulting richness and sampling fraction of each clade, before any grafting, are returned in the attribute `"richness_table"`.
#'
#' @return If `file` contains one type of taxa, a `total_richness_per_clades` list;
#'   if it contains both types (`$taxa` and `$clades`), a list with two `total_richness_per_clades` named `$taxa` (for [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]) and `$clades`
#'   (for [`phylograftR::graft_missing_clades()`][graft_missing_clades]).
#'
#'   Each `total_richness_per_clades` list contains all the clades of the table, and has the following items:
#'  * `$Target_clade` (or `$Sister_clade`, for the list of missing clades) List of character vectors of tips defining the clades.
#'  * `$Target_clade_name` (or `$Sister_clade_name`) Only if some clades have a name. Character vector (`NA` if a clade has no name).
#'  * `$total_richness` Integer vector of the total richness of the clades (`NA` if unknown).
#'
#'  The output has the attribute `"richness_table"`, a data.frame with one row per clade: `Clade_name`, `Clade_tips`, `total_richness_provided`
#'  (the value of the table), `total_richness` (as used for the grafting, derived from the nested clades if unknown), `nb_backbone_species`,
#'  `nb_missing_species` (all the taxa to graft in the clade, of both types), `nb_unsampled_species`, and `sampling_fraction` (`(nb_backbone_species + nb_missing_species) / total_richness`).
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::build_total_richness_per_clades()`][build_total_richness_per_clades] to convert a data.frame of clades already loaded in R.
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to read the table of missing taxa from an external file.
#'  * [`phylograftR::import_taxonomic_table()`][import_taxonomic_table] to read taxonomic tables from external files.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft the taxa using the richness.
#'
#' @export
#'
#' @author Maël Doré
#'
#' @examples
#' # Path to an example table of the total richness of clades, provided with the package
#' file_path <- system.file("extdata", "toy_total_richness.csv", package = "phylograftR")
#' read.csv(file_path)
#'
#' # Load the backbone, and the list of missing taxa associated with the table
#' data(toy_backbone)
#' data(toy_missing_taxa_list)
#'
#' # Import it as a 'total_richness_per_clades' object, checked against the backbone and the list
#' total_richness_per_clades <- import_total_richness_per_clades(file = file_path, phy = toy_backbone,
#'                                                               missing_taxa_list = toy_missing_taxa_list)
#' str(total_richness_per_clades, max.level = 1)
#'
#' # Use it to graft the taxa
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
#'                                   total_richness_per_clades = total_richness_per_clades,
#'                                   seed = 1234, verbose = FALSE)
#'

import_total_richness_per_clades <- function (file, # Path to the table
                                              phy, # Backbone phylogeny
                                              missing_taxa_list, # Missing taxa list, or list of two missing taxa lists ($taxa and $clades)
                                              sep = NULL, # Field separator of text files. By default, from the extension
                                              sheet = 1, # Sheet of Excel files
                                              reader = NULL, # Optional function to read other formats
                                              verbose = TRUE, # Print the conversions of node IDs into tips
                                              ...)
{
  ## 1.1/ Read the table ####

  # Checks the inputs, and reads the file with the reader matching its extension (or 'reader')
  data <- read_table_file(file = file, sep = sep, sheet = sheet, reader = reader, ...)

  ## 1.2/ Convert it ####

  build_total_richness_per_clades(data = data, phy = phy, missing_taxa_list = missing_taxa_list, verbose = verbose)
}
