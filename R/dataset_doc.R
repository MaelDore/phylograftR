
##### Documentation of the datasets #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - A backbone phylogeny for grafting
#  - A list of missing taxa to be grafted, with and without allowed_states
#  - A table of total richness per clades to inform sampling fractions with unsampled species
#  - A stochastic character map (simmap) of the backbone phylogeny, mapping a trait or range with three states


##### 1/ Toy backbone phylogeny #####

#' @title Toy backbone phylogeny
#'
#' @description Simulated ultrametric phylogeny with 30 tips (`t1` to `t30`) and a root age of 50, used in the examples and vignettes.
#'
#' Generated with `ape::rcoal()` (seed 1234) and rescaled. Internal branches shorter than 1 (2% of the root age) were lengthened by moving
#' their rootward node towards the root, so that no node looks like a polytomy on a plot (the topology, the root age and the ultrametry are unchanged).
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
#' @format A named list with 7 taxa and the following items. Every item is named after the taxa (`$Taxon_name`), in the same order.
#'  * `$Taxon_name` Character vector. Names of the taxa to graft.
#'  * `$Target_clade` List of character vectors of tips defining the target clade of each taxon (one tip: the terminal branch of this tip).
#'  * `$Target_clade_name` Character vector. Labels of the target clades. Ignored by the grafting functions.
#'  * `$include_target_clade_stem` Logical vector. Whether the stem of the target clade of each taxon is available for grafting (`NA` uses the global argument).
#'  * `$Excluded_subclades` List (one element per taxon) of named lists of character vectors of tips defining the excluded subclades, or `NULL`.
#'  * `$exclude_subclades_stem` List (one element per taxon) of logical vectors (one value per excluded subclade), or `NULL`.
#'
#'   See [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] for the use of the items in the grafting.
#'
#' @seealso [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list], [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]
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


##### 2b/ Toy missing taxa list with allowed states #####

#' @title Toy list of missing taxa, with allowed states
#'
#' @description Same as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], with the extra item `$allowed_states`: the states of
#'   [`phylograftR::toy_simmap`][toy_simmap] in which each taxon can be grafted. It is used to showcase the attribute-based exclusion
#'   of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] (see its examples).
#'
#'   It is the output of [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] applied to the table
#'   `system.file("extdata", "toy_missing_taxa_with_allowed_states.csv", package = "phylograftR")`.
#'   Four taxa have a taxon-level value (one or several states) for `$allowed_states` and three have none (`NULL`),
#'   so they use the global `allowed_states` argument of the grafting function.
#'   Such a list requires a `simmap` to graft: use [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list] to graft without a `simmap` instead.
#'
#' @usage data(toy_missing_taxa_list_with_allowed_states)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of character vectors of states of [`phylograftR::toy_simmap`][toy_simmap]
#'    in which the taxon can be grafted, or `NULL` to use the global `allowed_states` argument. `"Genus_a_sp1"`: `"red"`;
#'    `"Genus_b_sp1"`: `"blue"` and `"green"`; `"Genus_b_sp2"`: `"green"`; `"Genus_c_sp2"`: `"blue"` and `"red"`; the others: `NULL`.
#'
#' @seealso [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], [`phylograftR::toy_simmap`][toy_simmap],
#'   [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list]
#'
#' @examples
#' # Import the list
#' data(toy_missing_taxa_list_with_allowed_states)
#' toy_missing_taxa_list_with_allowed_states$allowed_states
#'
#' # Compare to the initial '.csv' file
#' system.file("extdata", "toy_missing_taxa_with_allowed_states.csv", package = "phylograftR")
#'
#' # Graft with the toy simmap
#' # Taxa with taxon-level states follow their specific constraints.
#' # Taxa without taxon-level states follow the global argument and can only be grafted on "blue" or "green" branches
#' data(toy_backbone)
#' data(toy_simmap)
#'
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone,
#'                                   missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
#'                                   simmap = toy_simmap, allowed_states = c("blue", "green"),
#'                                   verbose = FALSE, seed = 1234)
#'
"toy_missing_taxa_list_with_allowed_states"


##### 3/ Toy total richness data #####

#' @title Toy total richness of the target clades
#'
#' @description List usable as `total_richness_per_clades` with [`phylograftR::toy_backbone`][toy_backbone] and [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list].
#' It is built from the table `system.file("extdata", "toy_total_richness.csv", package = "phylograftR")`, with the tips of each clade separated by commas (see the examples).
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


##### 4/ Toy simmap #####

#' @title Toy stochastic character map of the toy backbone phylogeny
#'
#' @description Stochastic character map (simmap) of [`phylograftR::toy_backbone`][toy_backbone] with three states, `"blue"`, `"green"` and `"red"`,
#'   named after colors to ease their visualization. It can be used as `simmap` to restrict grafting to branches in some states
#'   (see the examples of the grafting functions).
#'
#'   It was simulated with `phytools::sim.history()` from a root in state `"green"`, with equal transition rates between the three states
#'   (0.04 per unit of time, i.e., one transition every 12.5 time units on average on a lineage, for a root age of 50).
#'   Among the first 1000 seeds, the one used (50) is the first giving a balanced map:
#'   * The three states cover similar shares of the total branch length (about 34% `"green"`, 32% `"red"` and 34% `"blue"`).
#'   * A moderate number of transitions: 14 changes of state, along 12 of the 58 branches (no branch has more than 3 pieces).
#'   * Each state is found in 8 to 11 of the 30 tips.
#'
#' @usage data(toy_simmap)
#' @format An object of classes `"simmap"` and `"phylo"` (see `phytools::make.simmap()`), with the same tips, topology and branch lengths as
#'   [`phylograftR::toy_backbone`][toy_backbone], and the additional elements `$maps` (for each branch, the duration spent in each state,
#'   from the rootward node to the tipward node), `$mapped.edge` (total duration per state and branch), `$node.states` and `$states` (states of the nodes and tips).
#'
#' @seealso [`phylograftR::toy_backbone`][toy_backbone], [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
#'
#' @examples
#' # Import the data
#' data(toy_simmap)
#'
#' # Plot the map, using less flashy colors than the ones of the states
#' state_colors <- c(blue = "dodgerblue", green = "limegreen", red = "firebrick")
#' phytools::plotSimmap(toy_simmap, colors = state_colors, fsize = 0.7, ftype = "i")
#'
#' # Duration spent in each state (sum on all the branches)
#' colSums(toy_simmap$mapped.edge)
#'
"toy_simmap"
