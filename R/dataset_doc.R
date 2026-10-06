
##### Documentation of the datasets #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - A backbone phylogeny for grafting
#  - A list of missing taxa to be grafted
#  - A table of total richness per clades to inform sampling fractions with unsampled species


##### 1/ Toy backbone phylogeny #####

#' @title Toy backbone phylogeny
#'
#' @description Simulated ultrametric phylogeny with 30 tips (`t1` to `t30`) and a root age of 50, used in the examples and vignettes.
#' Generated with `ape::rcoal()` (seed 1234) and rescaled. See `data-raw/toy_data.R`.
#'
#' @usage data(toy_backbone)
#' @format An object of class `"phylo"` with 30 tips and 29 internal nodes.
#'
#' @name toy_backbone
#'
#' @seealso [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list]
#'
#' @examples
#' # Import the phylogeny
#' data(toy_backbone)
#'
#' # Plot the phylogeny
#' ape::plot.phylo(toy_backbone, cex = 0.7)
#'
"toy_backbone"


##### 2/ Toy missing taxa list #####

#' @title Toy list of missing taxa
#'
#' @description Valid `missing_taxa_list` for [`phylograftR::toy_backbone`][toy_backbone],
#'   with 7 taxa to graft by [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#'
#'   It is the output of [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] applied to the table
#'   `system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")`: open that file to see how to fill the table.
#'   It shows target clades defined by one or two tips, excluded subclades (named or not, with one tip or several), and the stem settings.
#'
#' @usage data(toy_missing_taxa_list)
#' @format A named list with the items `$Taxon_name`, `$Target_clade`, `$Target_clade_name`, `$include_target_clade_stem`,
#'   `$Excluded_subclades` and `$exclude_subclades_stem`.
#'   See [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] for the description of the items.
#'
#' @seealso [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list], [`phylograftR::convert_missing_taxa_list()`][convert_missing_taxa_list]
#'
#' @examples
#' # Import the list from loaded datasets
#' data(toy_missing_taxa_list)
#' str(toy_missing_taxa_list, max.level = 1)
#'
#' # Compare to the initial '.csv' file
#' system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")
#'
"toy_missing_taxa_list"


##### 3/ Toy total richness data #####

#' @title Toy total richness of the target clades
#'
#' @description List usable as `total_richness_per_clades` with [`phylograftR::toy_backbone`][toy_backbone] and [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list].
#' It is built from the table `system.file("extdata", "toy_total_richness.csv", package = "phylograftR")`, with the tips of each clade separated by commas:
#'
#' ```

#' ```
#'
#' @usage data(toy_total_richness_per_clades)
#' @format A list with the items `$Target_clade` (list of character vectors of tips) and `$total_richness` (numeric vector, `NA` if unknown).
#'
#' @seealso [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#' @examples
#' # Import the data
#' data(toy_total_richness_per_clades)
#' str(toy_total_richness_per_clades)
#'
#' # Compare to the initial '.csv' file
#' system.file("extdata", "toy_total_richness.csv", package = "phylograftR")
#'
#' # Rebuild the data from the '.csv' table
#' tab <- read.csv(system.file("extdata", "toy_total_richness.csv", package = "phylograftR"))
#' total_richness_per_clades <- list(
#'   Target_clade = strsplit(tab$Target_clade, ",", fixed = TRUE),
#'   total_richness = tab$total_richness)
#'
"toy_total_richness_per_clades"
