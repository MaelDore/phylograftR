
##### Documentation of the datasets #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - A backbone phylogeny for grafting
#  - A list of missing taxa to be grafted, with and without allowed states or allowed ranges
#  - Taxonomic tables of the backbone and of the missing taxa, to build a list of missing taxa from taxonomic names
#  - A table of total richness per clades to inform sampling fractions with unsampled species
#  - A stochastic character map (simmap) of the backbone phylogeny, mapping a trait with three states
#  - A stochastic character map (simmap) of the backbone phylogeny, mapping biogeographic ranges with four states
#  - Tip data associated to all the simmaps: colors, elevation, ranges.
#  - Tip values fo a continuous trait associated with the contMap(s).


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
#' @seealso Other datasets from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list] for the list of missing taxa to graft on the backbone.
#'  * [`phylograftR::toy_taxonomic_tables`][toy_taxonomic_tables] for the taxonomic tables of the backbone and of the missing taxa.
#'  * [`phylograftR::toy_simmap`][toy_simmap] for the stochastic character map of colors.
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
#'   Clades and excluded subclades are named after the taxonomic groups of [`phylograftR::toy_backbone_taxonomic_table`][toy_taxonomic_tables] (e.g., `"Family_A"`).
#'
#' @usage data(toy_missing_taxa_list)
#' @format A named list with 7 taxa and the following items. **Every item is named after the taxa (`$Taxon_name`), in the same order.**
#'  * `$Taxon_name` Character vector. Names of the taxa to graft.
#'  * `$Target_clade` List of character vectors of tips defining the target clade of each taxon (one tip: the terminal branch of this tip).
#'  * `$Target_clade_name` Character vector. Labels of the target clades. Not used for the grafting, but checked against the names of `total_richness_per_clades`, if provided.
#'  * `$include_target_clade_stem` Logical vector. Whether the stem of the target clade of each taxon is available for grafting (`NA` uses the global argument).
#'  * `$Excluded_subclades` List (one element per taxon) of named lists of character vectors of tips defining the excluded subclades, or `NULL`.
#'  * `$exclude_subclades_stem` List (one element per taxon) of logical vectors (one value per excluded subclade), or `NULL`.
#'
#'   See [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] for the use of the items in the grafting.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_backbone`][toy_backbone] for the backbone phylogeny.
#'  * [`phylograftR::toy_total_richness_per_clades`][toy_total_richness_per_clades] for the total richness of the target clades.
#'  * [`phylograftR::toy_taxonomic_tables`][toy_taxonomic_tables] for the taxonomic tables of the backbone and of the missing taxa.
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to import such a list from an external table.
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] to build such a list from a data.frame.
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
#'    in which the taxon can be grafted, or `NULL` to use the global `allowed_states` argument. `"New_sp1"`: `"red"`;
#'    `"New_sp4"`: `"blue"` and `"green"`; `"New_sp5"`: `"green"`; `"New_sp7"`: `"blue"` and `"red"`; the others: `NULL`.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list] for the list of missing taxa to graft on the backbone, without allowed states.
#'  * [`phylograftR::toy_simmap`][toy_simmap] for the stochastic character map of colors.
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to import such a list from an external table.
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] to build such a list from a data.frame.
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
#' # Taxa without taxon-level states follow the global argument
#' # and can only be grafted on "blue" or "green" branches.
#' data(toy_backbone)
#' data(toy_simmap)
#'
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
#'     simmap = toy_simmap, allowed_states = c("blue", "green"),
#'     verbose = FALSE, seed = 1234)
#'
#' # Plot the updated simmap
#' plot_grafted(grafted_phy, fsize = 0.5, ftype = "i",
#'              simmap_colors = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))
#'
#' # Inspect final states of each grafted taxon
#' taxon_names <- toy_missing_taxa_list_with_allowed_states$Taxon_name
#' tip_edges <- match(match(taxon_names, grafted_phy$tip.label), grafted_phy$edge[, 2])
#' taxon_states <- setNames(object = vapply(grafted_phy$simmap$maps[tip_edges],
#'                                          function (x) names(x)[1], character(1)),
#'                          nm = taxon_names)
#' # Should match with the allowed ranges.
#' # NULL default to "blue" or "green" based on the 'allowed_states' argument.
#' taxon_states
#' toy_missing_taxa_list_with_allowed_states$allowed_states
#'
#'
"toy_missing_taxa_list_with_allowed_states"


##### 2c/ Toy missing taxa list with allowed ranges #####

#' @title Toy list of missing taxa, with allowed ranges
#'
#' @description Same as [`phylograftR::toy_missing_taxa_list_with_allowed_states`][toy_missing_taxa_list_with_allowed_states], with biogeographic ranges as allowed states:
#'   the `$allowed_states` item holds ranges of [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges], coded as in BioGeoBEARS, in which each taxon can be grafted.
#'   It is used to showcase the grafting with biogeographic ranges as states, with `treat_states_as_ranges = TRUE`
#'   (see the Details section "Biogeographic ranges as states" of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#'
#'   It is the output of [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] applied to the table
#'   `system.file("extdata", "toy_missing_taxa_with_allowed_ranges.csv", package = "phylograftR")`.
#'   Four taxa have a taxon-level value (one or several ranges) for `$allowed_states` and three have none (`NULL`),
#'   so they use the global `allowed_states` argument of the grafting function.
#'   All the allowed ranges are states of the toy simmap, so grafting does not raise warnings. Some are contained in other states of the simmap
#'   (`"A"` in `"AB"`, `"B"` in `"AB"` and `"BC"`): a lineage grafted in a branch in `"AB"` with the allowed range `"A"` inherits `"A"`,
#'   which is a subset sympatric speciation event, `(AB) -> (AB), (A)`.
#'   **Such a list requires a `simmap` to graft.**
#'
#' @usage data(toy_missing_taxa_list_with_allowed_ranges)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of character vectors of ranges of [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges]
#'    in which the taxon can be grafted, or `NULL` to use the global `allowed_states` argument. `"New_sp1"`: `"A"`;
#'    `"New_sp4"`: `"B"` and `"BC"`; `"New_sp5"`: `"B"`; `"New_sp7"`: `"A"` and `"BC"`; the others: `NULL`.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list_with_allowed_states`][toy_missing_taxa_list_with_allowed_states] for the same list, with allowed states of colors.
#'  * [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges] for the stochastic character map of biogeographic ranges.
#'  * [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] to import such a list from an external table.
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] to build such a list from a data.frame.
#'
#' @examples
#' # Import the list
#' data(toy_missing_taxa_list_with_allowed_ranges)
#' toy_missing_taxa_list_with_allowed_ranges$allowed_states
#'
#' # Compare to the initial '.csv' file
#' system.file("extdata", "toy_missing_taxa_with_allowed_ranges.csv", package = "phylograftR")
#'
#' # Graft with the toy simmap of ranges, treating the states as ranges
#' # Taxa with taxon-level ranges follow their specific constraints.
#' # Taxa without taxon-level ranges follow the global argument
#' # and can only be grafted on branches in a range containing "A" or "B".
#' data(toy_backbone)
#' data(toy_simmap_ranges)
#'
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone,
#'                                   missing_taxa_list = toy_missing_taxa_list_with_allowed_ranges,
#'                                   simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
#'                                   treat_states_as_ranges = TRUE,
#'                                   verbose = FALSE, seed = 1234)
#'
#' # Plot the updated simmap
#' plot_grafted(grafted_phy,
#'              simmap_colors = c("A" = "gold", "B" = "dodgerblue",
#'                                "AB" = "limegreen", "BC" = "firebrick2"))
#'
#' # Inspect range inherited by each grafted taxon (state of its terminal branch in the updated simmap)
#' taxon_names <- toy_missing_taxa_list_with_allowed_ranges$Taxon_name
#' tip_edges <- match(match(taxon_names, grafted_phy$tip.label), grafted_phy$edge[, 2])
#' taxon_ranges <- setNames(object = vapply(grafted_phy$simmap$maps[tip_edges],
#'                                          function (x) names(x)[1], character(1)),
#'                          nm = taxon_names)
#' # Should match with the allowed ranges. NULL default to "A" or "B"
#' # based on the 'allowed_ranges' argument.
#' taxon_ranges
#' toy_missing_taxa_list_with_allowed_ranges$allowed_states
#'
"toy_missing_taxa_list_with_allowed_ranges"


##### 3/ Toy total richness data #####

#' @title Toy total richness of the target clades
#'
#' @description *List usable as `total_richness_per_clades` with [`phylograftR::toy_backbone`][toy_backbone] and [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list].
#' It is built from the table `system.file("extdata", "toy_total_richness.csv", package = "phylograftR")`, with the tips of each clade separated by commas (see the examples).
#'
#' @usage data(toy_total_richness_per_clades)
#' @format A list with the items `$Target_clade` (list of character vectors of tips), `$Target_clade_name` (character vector of the names of the clades,
#'   the same as in [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list]) and `$total_richness` (numeric vector, `NA` if unknown).
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list] for the list of missing taxa to graft on the backbone.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade, using this dataset.
#'  * [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract the richness and the sampling fraction of the clades of grafted phylogenies.
#'
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
#'   Target_clade_name = tab$Target_clade_name,
#'   total_richness = tab$total_richness)
#'
"toy_total_richness_per_clades"


##### 4/ Toy taxonomic tables #####

#' @title Toy taxonomic tables
#'
#' @description Tables to build a `missing_taxa_list` from taxonomic names with
#'   [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy], consistent with [`phylograftR::toy_backbone`][toy_backbone].
#'
#'   `toy_backbone_taxonomic_table` gives the taxonomy of the 30 tips of the backbone, with four levels (`Family`, `Subfamily`, `Tribe`, `Genus`).
#'   **All the groups are monophyletic.** Some tips have no group at the finest levels (`NA`), as often in real data.
#'
#'   `toy_missing_taxa_taxonomic_table` describes 10 taxa to graft. The first seven (**`Type = "taxa"`**) are the taxa of
#'   [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], written with taxonomic names: the clades, excluded subclades and stem settings are the same.
#'   The last three (**`Type = "clades"`**) are taxa to graft as sister clades with [`phylograftR::graft_missing_clades()`][graft_missing_clades]:
#'   `Genus_a_sp1` and `Genus_a_sp2` belong to a new genus sister to the genus `Gen_1`, so they are grafted together as a missing clade, and `Genus_b_sp1` is sister to its own tribe.
#'   The target and sister clades are written with the name of a group (e.g., `"Family_B"`) or of a level (e.g., `"Family"`, meaning the family of the taxon itself).
#'   Both tables are also provided as `.csv` files in `system.file("extdata", package = "phylograftR")`
#'   (`toy_backbone_taxonomic_table.csv` and `toy_missing_taxa_taxonomic_table.csv`).
#'
#'   **The clades are named after their taxonomic groups** (e.g., `"Family_A"`), as in [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list]
#'   and [`phylograftR::toy_total_richness_per_clades`][toy_total_richness_per_clades], so the lists built from the taxonomy can be used with this richness list.
#'
#' @usage
#' data(toy_backbone_taxonomic_table)
#' data(toy_missing_taxa_taxonomic_table)
#'
#' @format
#' `toy_backbone_taxonomic_table` is a data frame with 30 rows (one per tip of the backbone) and the columns `Taxa` (tip labels), `Family`, `Subfamily`, `Tribe` and `Genus`.
#'
#' `toy_missing_taxa_taxonomic_table` is a data frame with 10 rows (one per taxon to graft) and the columns `Taxa`, `Type`, the taxonomic levels `Family`, `Subfamily`, `Tribe` and `Genus`
#' of the taxa (`NA` if unknown), `Target_clade`, `Sister_clade`, `Excluded_subclades`, `include_target_clade_stem` and `exclude_subclades_stem`.
#'
#' @name toy_taxonomic_tables
#' @aliases toy_backbone_taxonomic_table toy_missing_taxa_taxonomic_table
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_backbone`][toy_backbone] for the backbone phylogeny.
#'  * [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list] for the list of missing taxa to graft on the backbone.
#'  * [`phylograftR::import_taxonomic_table()`][import_taxonomic_table] to read such tables from files.
#'  * [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy] to build `missing_taxa_list` objects from these tables.
#'
#' @examples
#' data(toy_backbone)
#' data(toy_backbone_taxonomic_table)
#' data(toy_missing_taxa_taxonomic_table)
#'
#' head(toy_backbone_taxonomic_table)
#' toy_missing_taxa_taxonomic_table
#'
#' # Build the lists: one for graft_missing_taxa() and one for graft_missing_clades()
#' toy_lists <- build_missing_taxa_list_from_taxonomy(
#'   backbone_taxonomic_table = toy_backbone_taxonomic_table,
#'   missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
#'   phy = toy_backbone)
#' str(toy_lists$taxa, max.level = 1)
#' str(toy_lists$clades, max.level = 1)
#'
#' # Compare to the initial '.csv' files
#' system.file("extdata", "toy_backbone_taxonomic_table.csv", package = "phylograftR")
#' system.file("extdata", "toy_missing_taxa_taxonomic_table.csv", package = "phylograftR")
#'
NULL

#' @rdname toy_taxonomic_tables
"toy_backbone_taxonomic_table"

#' @rdname toy_taxonomic_tables
"toy_missing_taxa_taxonomic_table"


##### 5/ Toy simmap #####

#' @title Toy stochastic character map of the toy backbone phylogeny
#'
#' @description Stochastic character map (simmap) of [`phylograftR::toy_backbone`][toy_backbone] with **three states**, `"blue"`, `"green"` and `"red"`,
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
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_backbone`][toy_backbone] for the backbone phylogeny.
#'  * [`phylograftR::toy_missing_taxa_list_with_allowed_states`][toy_missing_taxa_list_with_allowed_states] for the list of missing taxa to graft, with associated allowed states of colors.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade based on color states, using this dataset.
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


##### 5b/ Toy simmap of biogeographic ranges #####

#' @title Toy stochastic character map of biogeographic ranges of the toy backbone phylogeny
#'
#' @description Stochastic character map (simmap) of [`phylograftR::toy_backbone`][toy_backbone] with **four states that are biogeographic ranges**
#'   coded as in BioGeoBEARS: `"A"`, `"B"`, `"AB"` and `"BC"` (areas A, B and C; the area C is only found in `"BC"`).
#'   It can be used as `simmap` with `treat_states_as_ranges = TRUE` to graft taxa on branches in ranges that contain their allowed ranges
#'   (see the Details section "Biogeographic ranges as states" of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and the examples of
#'   [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`][toy_missing_taxa_list_with_allowed_ranges]).
#'
#'   It was simulated with `phytools::sim.history()` from a root in range `"AB"`, with equal rates (0.04 per unit of time) for the anagenetic changes of range
#'   that gain or lose a single area (`"A"` <-> `"AB"`, `"B"` <-> `"AB"` and `"B"` <-> `"BC"`). A change of two areas at once (`"AB"` <-> `"BC"`) is not allowed.
#'   As a simulated history only has anagenetic changes, some of them were then moved to the parental node of their branch, to also represent cladogenetic events
#'   (see the Details section). Among the first 20000 seeds, the one used (3598) is the first giving a balanced map:
#'   * The four ranges cover similar shares of the total branch length (about 30% `"B"`, 26% `"AB"`, 25% `"A"` and 20% `"BC"`).
#'   * A moderate number of changes of range: 10 of the 58 branches include a change of range (no branch has more than 3 pieces).
#'   * Each range is found in 5 to 10 of the 30 tips.
#'
#' @details
#' Events represented in the map, with the syntax of BioGeoBEARS.
#'
#' **Anagenetic changes** (along the branches, 12 in total):
#' * 8 `"d"` (range expansion, gain of one area): 5 `"A"` -> `"AB"` and 3 `"B"` -> `"BC"`.
#' * 4 `"e"` (range contraction, loss of one area): 3 `"AB"` -> `"A"` and 1 `"AB"` -> `"B"`.
#'
#' **Cladogenetic events** (at the 29 nodes, from the range of the lineage above the node and the ranges at the start of the two descendant branches):
#' * 5 `"s"` (subset sympatry), such as `(AB) -> (AB), (A)`: one descendant keeps the widespread range, the other starts in a single area within it.
#'   They are the events that are produced by the grafting functions when a taxon inherits an allowed range smaller than the range of the lineage it is grafted on
#'   (subset sympatric speciation, see [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#' * 1 `"v"` (vicariance), `(AB) -> (A), (B)`: the widespread range is split between the two descendants.
#' * 13 `"y"` (narrow sympatry): both descendants inherit the single-area range of the lineage, such as `(A) -> (A), (A)`.
#' * 10 nodes where both descendants inherit the same widespread range (widespread sympatry): `(AB) -> (AB), (AB)` or `(BC) -> (BC), (BC)`
#'   (allowed in BAYAREALIKE-type models, but not in DEC).
#'
#' @usage data(toy_simmap_ranges)
#' @format An object of classes `"simmap"` and `"phylo"` (see `phytools::make.simmap()`), with the same tips, topology and branch lengths as
#'   [`phylograftR::toy_backbone`][toy_backbone], and the additional elements `$maps` (for each branch, the duration spent in each range,
#'   from the rootward node to the tipward node), `$mapped.edge` (total duration per range and branch), `$node.states` and `$states` (ranges of the nodes and tips).
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`][toy_missing_taxa_list_with_allowed_ranges] for the list of missing taxa to graft, with associated allowed biogeographic ranges.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade based on geographic ranges, using this dataset.
#'
#' @examples
#' # Import the data
#' data(toy_simmap_ranges)
#'
#' # Plot the map
#' range_colors <- c(A = "gold", B = "dodgerblue", AB = "limegreen", BC = "firebrick2")
#' phytools::plotSimmap(toy_simmap_ranges, colors = range_colors, fsize = 0.7, ftype = "i")
#'
#' # Duration spent in each range (sum on all the branches)
#' colSums(toy_simmap_ranges$mapped.edge)
#'
"toy_simmap_ranges"


##### 6/ Toy datasets for a second trait: simmap and list of missing taxa with allowed states of elevation #####

#' @title Toy stochastic character map of elevation states
#'
#' @description Stochastic character map (simmap) of [`phylograftR::toy_backbone`][toy_backbone] with **two elevation states**, `"high"` and `"low"`,
#'   to combine with [`phylograftR::toy_simmap`][toy_simmap] (color) or [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges] (biogeographic ranges)
#'   as examples of attribute-based constraints applied on two traits at once.
#'
#'   It was simulated with `phytools::sim.history()` from a root in state `"low"`, with equal transition rates between the two states (0.04 per unit of time).
#'   Among the first 1000 seeds, the one used (94) is the first giving a map that is balanced and compatible with the constraints of
#'   [`phylograftR::toy_missing_taxa_list_with_two_traits`][toy_missing_taxa_list_with_two_traits] used with [`phylograftR::toy_simmap`][toy_simmap]
#'   (at least one branch in the allowed states of each taxon):
#'   * The two states cover similar shares of the total branch length (54% `"low"` and 46% `"high"`).
#'   * A moderate number of transitions: 12 changes of state, along 9 of the 58 branches (no branch has more than 4 pieces).
#'   * `"high"` is found in 18 of the 30 tips and `"low"` in 12.
#'
#' @usage data(toy_simmap_elevation)
#' @format An object of classes `"simmap"` and `"phylo"` (see `phytools::make.simmap()`), with the same tips, topology and branch lengths as
#'   [`phylograftR::toy_backbone`][toy_backbone], and the additional elements `$maps` (for each branch, the duration spent in each state,
#'   from the rootward node to the tipward node), `$mapped.edge` (total duration per state and branch), `$node.states` and `$states` (states of the nodes and tips).
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_simmap`][toy_simmap] for the stochastic character map of colors.
#'  * [`phylograftR::toy_missing_taxa_list_with_two_traits`][toy_missing_taxa_list_with_two_traits] for the list of missing taxa to graft, with allowed states of colors and elevation.
#'  * [`phylograftR::toy_missing_taxa_list_with_ranges_and_elevation`][toy_missing_taxa_list_with_ranges_and_elevation] for the list of missing taxa to graft, with allowed ranges and elevation.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade based on elevation, using this dataset.
#'
#' @examples
#' # Import the data
#' data(toy_simmap_elevation)
#'
#' # Plot the map
#' plot(toy_simmap_elevation, colors = c(high = "darkorange", low = "steelblue"),
#'      fsize = 0.7, ftype = "i")
#'
#' # Duration spent in each state (sum on all the branches)
#' colSums(toy_simmap_elevation$mapped.edge)
#'
"toy_simmap_elevation"


#' @title Toy list of missing taxa, with allowed states of two traits: colors and elevation
#'
#' @description Same as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], with the extra item `$allowed_states` covering the states of two traits,
#'   [`phylograftR::toy_simmap`][toy_simmap] (`"blue"`, `"green"`, `"red"`) and [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation] (`"high"`, `"low"`),
#'   in which each taxon can be grafted. It is used to showcase the grafting constrained by several traits at once
#'   (see Several traits in Details section of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#'
#'   It is the output of [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] applied to the table
#'   `system.file("extdata", "toy_missing_taxa_with_two_traits.csv", package = "phylograftR")`.
#'   In the column `allowed_states`, **the states of the traits are separated by semicolons**
#'   (in the order of the maps provided in `simmap`: the color, then the elevation), **and the states of a trait by commas**. Open that file to see the different syntaxes:
#'   * `"red;"`: the color is `"red"`, and the elevation has no taxon-level value (the global value is used).
#'   * `"blue, green;high"`: color `"blue"` or `"green"`, and elevation `"high"`.
#'   * `";low"`: no taxon-level value for the color, elevation `"low"`.
#'   * `"elevation=high"`: traits named, so the color can be omitted (it has no taxon-level value).
#'   * Blank cell: no taxon-level value for any trait, so the global `allowed_states` argument is used for both.
#'
#' @usage data(toy_missing_taxa_list_with_two_traits)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of lists with the states of each trait (one element per trait, `NULL` for a trait without taxon-level value),
#'    or `NULL` to use the global `allowed_states` argument for all the traits. `"New_sp1"`: `list("red", NULL)`; `"New_sp4"`: `list(c("blue", "green"), "high")`;
#'    `"New_sp5"`: `list(NULL, "low")`; `"New_sp7"`: `list(elevation = "high")`; the others: `NULL`.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list_with_allowed_states`][toy_missing_taxa_list_with_allowed_states] for the same list, with allowed states of colors.
#'  * [`phylograftR::toy_simmap`][toy_simmap] for the stochastic character map of colors.
#'  * [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation] for the stochastic character map of elevation.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade based on colors and elevation, using this dataset.
#'
#' @examples
#' # Import the list
#' data(toy_missing_taxa_list_with_two_traits)
#' toy_missing_taxa_list_with_two_traits$allowed_states
#'
#' # Compare to the initial '.csv' file
#' system.file("extdata", "toy_missing_taxa_with_two_traits.csv", package = "phylograftR")
#'
#' # Graft with the maps of the two traits. The global 'allowed_states' has one element per trait:
#' # grafting on "blue" or "green" branches for the color, and no restriction for the elevation (NULL).
#' data(toy_backbone)
#' data(toy_simmap)
#' data(toy_simmap_elevation)
#'
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_two_traits,
#'     simmap = list(color = toy_simmap, elevation = toy_simmap_elevation),
#'     allowed_states = list(c("blue", "green"), NULL),
#'     verbose = FALSE, seed = 1234)
#'
#' # One updated simmap per trait
#' names(grafted_phy$simmap)
#'
#' # State of each grafted taxon for each trait: should match the allowed states
#' taxon_names <- toy_missing_taxa_list_with_two_traits$Taxon_name
#' taxon_states <- sapply(grafted_phy$simmap, function (simmap)
#' {
#'   tip_edges <- match(match(taxon_names, simmap$tip.label), simmap$edge[, 2])
#'   setNames(vapply(simmap$maps[tip_edges], function (x) names(x)[1], character(1)), taxon_names)
#' })
#' taxon_states
#' toy_missing_taxa_list_with_two_traits$allowed_states
#'
"toy_missing_taxa_list_with_two_traits"


#' @title Toy list of missing taxa, with allowed states of two traits: biogeographic ranges and elevation
#'
#' @description Same as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], with the extra item `$allowed_states` covering the states of two traits,
#'   the biogeographic ranges of [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges] (`"A"`, `"B"`, `"AB"`, `"BC"`) and the elevation of
#'   [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation] (`"high"`, `"low"`), in which each taxon can be grafted.
#'   It is used to showcase the grafting constrained by a range trait and a 'normal' trait at once,
#'   with `treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE)`
#'   (see the Details section "Biogeographic ranges as states" of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#'
#'   It is the output of [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list] applied to the table
#'   `system.file("extdata", "toy_missing_taxa_with_ranges_and_elevation.csv", package = "phylograftR")`.
#'   In the column `allowed_states`, **the states of the traits are separated by semicolons**
#'   (in the order of the maps provided in `simmap`: the ranges, then the elevation), **and the states of a trait by commas**. Open that file to see the different syntaxes:
#'   * `"A;high"`: the range is `"A"`, and the elevation is `"high"`.
#'   * `"B, BC;"`: the range is `"B"` or `"BC"`, and the elevation has no taxon-level value (the global value is used).
#'   * `"elevation=low"`: traits named, so the ranges can be omitted (they have no taxon-level value).
#'   * `"biogeo=A, BC; elevation=high"`: all the traits named, in any order.
#'   * Blank cell: no taxon-level value for any trait, so the global `allowed_states` argument is used for both.
#'
#'   **The names of the traits (`biogeo` and `elevation`) must match the names of the maps provided in `simmap`** to use the syntax with names.
#'
#' @usage data(toy_missing_taxa_list_with_ranges_and_elevation)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of lists with the states of each trait (one element per trait, `NULL` for a trait without taxon-level value),
#'    or `NULL` to use the global `allowed_states` argument for all the traits. `"New_sp1"`: `list("A", "high")`; `"New_sp4"`: `list(c("B", "BC"), NULL)`;
#'    `"New_sp5"`: `list(elevation = "low")`; `"New_sp7"`: `list(biogeo = c("A", "BC"), elevation = "high")`; the others: `NULL`.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`][toy_missing_taxa_list_with_allowed_ranges] for the same missing taxa list, with allowed biogeographic ranges.
#'  * [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges] for the stochastic character map of biogeographic ranges.
#'  * [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation] for the stochastic character map of elevation.
#'  * [`phylograftR::toy_missing_taxa_list_with_two_traits`][toy_missing_taxa_list_with_two_traits] for the list with allowed states of colors and elevation.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade based on ranges and elevation, using this dataset.
#'
#' @examples
#' # Import the list
#' data(toy_missing_taxa_list_with_ranges_and_elevation)
#' toy_missing_taxa_list_with_ranges_and_elevation$allowed_states
#'
#' # Compare to the initial '.csv' file
#' system.file("extdata", "toy_missing_taxa_with_ranges_and_elevation.csv", package = "phylograftR")
#'
#' # Graft with the maps of the two traits. The first trait is treated as ranges,
#' # the second as 'normal', non-biogeographic states.
#' # Taxa without taxon-level states can be grafted on branches
#' # in a range containing "A" or "B" (global value),
#' # whatever the elevation (NULL).
#'
#' # Load datasets
#' data(toy_backbone)
#' data(toy_simmap_ranges)
#' data(toy_simmap_elevation)
#'
#' # Graft with constraints from elevation states and biogeographic ranges
#' grafted_phy <- graft_missing_taxa(phy = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
#'     simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
#'     allowed_states = list(c("A", "B"), NULL),
#'     treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
#'     verbose = FALSE, seed = 1234)
#'
#' # Result stores one updated simmap per trait
#' names(grafted_phy$simmap)
#'
#' # State of each grafted taxon for each trait: should match the allowed states
#' taxon_names <- toy_missing_taxa_list_with_ranges_and_elevation$Taxon_name
#' taxon_states <- sapply(grafted_phy$simmap, function (simmap)
#' {
#'   tip_edges <- match(match(taxon_names, simmap$tip.label), simmap$edge[, 2])
#'   setNames(vapply(simmap$maps[tip_edges], function (x) names(x)[1], character(1)), taxon_names)
#' })
#' taxon_states
#' toy_missing_taxa_list_with_ranges_and_elevation$allowed_states
#'
#' # Check that all the grafted taxa respect their constraints, trait by trait
#' check_grafting_constraints(
#'     phy_grafted = grafted_phy,
#'     phy_backbone = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
#'     allowed_states = list(c("A", "B"), NULL),
#'     treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE))
#'
"toy_missing_taxa_list_with_ranges_and_elevation"


##### 7/ Toy datasets for a continuous trait: tip values, contMap and list of contMaps #####

#' @title Toy tip values of a continuous trait
#'
#' @description Values of a continuous trait for the 30 tips of [`phylograftR::toy_backbone`][toy_backbone].
#'   They are the data used to build [`phylograftR::toy_contMap`][toy_contMap] and [`phylograftR::toy_contMaps`][toy_contMaps],
#'   examples of **continuous traits discretized into states** to constrain the grafting (see [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap]).
#'
#'   They were simulated with `phytools::fastBM()` under a Brownian Motion (ancestral value of 0 at the root, evolutionary rate `sig2 = 1`), with seed 37.
#'   Among the first 150 seeds, this is one of the seeds for which the discretization of the contMap into three states using the terciles of the tip values
#'   gives states that cover similar shares of the branches (35% `"low"`, 43% `"medium"` and 22% `"high"`), and
#'   constraints that differ among the 7 toy missing taxa of [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list].
#'
#' @usage data(toy_tip_values)
#' @format A named numeric vector of length 30, named after the tips of [`phylograftR::toy_backbone`][toy_backbone] and in the same order.
#'   The values range from about -7.5 to 4.2.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_contMap`][toy_contMap] for the contMap of this trait (point estimate).
#'  * [`phylograftR::toy_contMaps`][toy_contMaps] for 10 stochastic histories of this trait (list of contMaps).
#'  * [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap] to discretize a contMap into a simmap, and
#'    [`phylograftR::convert_values_to_states()`][convert_values_to_states] to convert continuous values into the same states.
#'
#' @examples
#' # Import the data
#' data(toy_tip_values)
#'
#' # Distribution of the values
#' summary(toy_tip_values)
#'
#' # Terciles of the tip values: cut points used to define three states
#' quantile(toy_tip_values, probs = c(1/3, 2/3))
#'
"toy_tip_values"


#' @title Toy contMap of a continuous trait
#'
#' @description A `contMap` of [`phylograftR::toy_backbone`][toy_backbone] for the continuous trait of [`phylograftR::toy_tip_values`][toy_tip_values],
#'   built with `phytools::contMap()`. It is a **point estimate**: the ancestral values are estimated by maximum likelihood under a Brownian Motion
#'   and **linearly interpolated along the branches**, so it does not cover the uncertainty of the reconstruction.
#'
#'   It is an example of input of [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap], which discretizes the trait values into states
#'   (e.g., `"low"`, `"medium"`, `"high"`) to obtain a `simmap` that can **constrain the grafting**. The resulting `simmap` is a **constraint, not a reconstruction of a discrete trait**.
#'
#' @usage data(toy_contMap)
#' @format An object of class `"contMap"` (see `phytools::contMap()`) with the elements `$tree` (the backbone as a `"simmap"`, in which the states are
#'   the color levels of the trait values), `$cols` (the color palette) and `$lims` (range of the trait values, about -7.5 to 4.2).
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_tip_values`][toy_tip_values] for the data used to build the contMap.
#'  * [`phylograftR::toy_contMaps`][toy_contMaps] for 10 stochastic histories of the trait (list of contMaps).
#'  * [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap] to discretize a contMap into a simmap.
#'
#' @examples
#' # Import the data
#' data(toy_contMap)
#' data(toy_tip_values)
#'
#' # Plot the contMap
#' plot(toy_contMap, legend = FALSE, fsize = 0.7)
#'
#' # Discretize into a simmap of three states using the terciles of the tip values.
#' # Provide the observed tip values: the contMap only approximates them (interpolation within segments)
#' my_simmap <- convert_contMap_to_simmap(toy_contMap, probs = c(1/3, 2/3),
#'                                        states = c("low", "medium", "high"),
#'                                        tip_values = toy_tip_values)
#' attr(my_simmap, "breaks")
#' plot(my_simmap, colors = c(low = "blue", medium = "grey", high = "red"),
#'      fsize = 0.7, ftype = "i")
#'
"toy_contMap"


#' @title Toy list of contMaps: stochastic histories of a continuous trait
#'
#' @description List of 10 `contMap`s of [`phylograftR::toy_backbone`][toy_backbone], each describing **one possible history of the continuous
#'   trait** of [`phylograftR::toy_tip_values`][toy_tip_values] (continuous stochastic maps).
#'   Providing several histories accounts for the **uncertainty** of the reconstruction of the ancestral trait values, unlike [`phylograftR::toy_contMap`][toy_contMap].
#'
#'   It is the type of input that you need to **constrain the grafting with a continuous trait while accounting for the uncertainty**:
#'   the first step is to **discretize the list into a `multiSimmap`** with [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap]
#'   (see the examples), and to provide it as `simmap` to the grafting functions.
#'
#'   It was built with the script of the example of [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap]:
#'   1. 100 continuous stochastic maps were simulated with `contsimmap::make.contsimmap()`, conditioned on the tip values, the backbone and a Brownian Motion
#'   with the evolutionary rate estimated by `phytools::anc.ML()` (seed 1234).
#'   2. They were converted into a list of `contMap`s with `deepSTRAPP::convert_contsimmap_to_contMaps()`.
#'   3. Only the **first 10 maps** were kept, to keep the package light. 10 maps are enough for illustration, but **actual analyses should use more maps (100 or more)**.
#'
#'   The script uses `contsimmap` directly, to credit the authors of the continuous stochastic mapping approach.
#'   For your own data, `deepSTRAPP::prepare_trait_data()` is **recommended**: with `run_stochastic_maps = TRUE`, it fits the evolutionary model(s), runs the continuous stochastic mapping
#'   and converts the maps into a list of `contMap`s in a single call, with additional options to fine-tune the process.
#'   It can also output directly both the interpolated point estimate as in `phytools::contMap()`,
#'   and the set of stochastic continuous maps simulated with `contsimmap::make.contsimmap()`.
#'
#'   **Once discretized, these maps are a constraint for the grafting. They are not stochastic maps of a discrete trait:**
#'   the number of transitions between states depends on the cut points, on the resolution of the continuous maps and on `min_segment_proportion`,
#'   and must not be used to compute transition rates or counts.
#'
#' @usage data(toy_contMaps)
#' @format A list of 10 objects of class `"contMap"` (see `phytools::contMap()`), each with the elements:
#'   * `$tree`: the backbone as a `"simmap"`, in which the states are the color levels of the trait values from 0 to 1000.
#'   * `$cols`: the color palette across values discretized from 0 to 1000.
#'   * `$lims`: the range of the trait values of the map, which differs among the maps.
#' The topology and branch lengths of the trees are the same as [`phylograftR::toy_backbone`][toy_backbone].
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_tip_values`][toy_tip_values] and [`phylograftR::toy_contMap`][toy_contMap] for the data used and the point estimate of the trait history.
#'  * [`phylograftR::convert_contMap_to_simmap()`][convert_contMap_to_simmap] to discretize contMaps into simmaps.
#'  * [`phylograftR::convert_values_to_states()`][convert_values_to_states] to convert the continuous values of the missing taxa into the states of these maps.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa with the maps as constraints: each map is used for a different simulation.
#'
#' @examples
#' # Import the data
#' data(toy_contMaps)
#' data(toy_tip_values)
#' length(toy_contMaps)
#'
#' # Discretize the list of contMaps into a multiSimmap of three states, using the terciles of the observed tip values as cut points
#' my_multiSimmap <- convert_contMap_to_simmap(toy_contMaps, probs = c(1/3, 2/3),
#'                                             states = c("low", "medium", "high"),
#'                                             tip_values = toy_tip_values, verbose = FALSE)
#'
#' # Cut points and states used for the discretization
#' attr(my_multiSimmap, "breaks")
#' attr(my_multiSimmap, "states")
#'
#' # Plot the first two maps: two possible histories of the same trait
#' par(mfrow = c(1, 2))
#' state_colors <- c(low = "blue", medium = "grey", high = "red")
#' plot(my_multiSimmap[[1]], colors = state_colors, fsize = 0.6, ftype = "i")
#' plot(my_multiSimmap[[2]], colors = state_colors, fsize = 0.6, ftype = "i")
#' par(mfrow = c(1, 1))
#'
"toy_contMaps"

##### 8/ Toy tip data and stochastic maps inferred from them #####

#' @title Toy tip data of three species attributes
#'
#' @description States of three **species attributes** for the 30 tips of [`phylograftR::toy_backbone`][toy_backbone]: a color, an elevation and a biogeographic range.
#'   They are the **tip data** from which stochastic maps are inferred in `vignette("prepare_stochastic_maps")`
#'   (see [`phylograftR::toy_multiSimmap`][toy_multiSimmap], [`phylograftR::toy_multiSimmap_elevation`][toy_multiSimmap_elevation] and
#'   [`phylograftR::toy_multiSimmap_ranges`][toy_multiSimmap_ranges]).
#'
#'   They are the states found at the tips of [`phylograftR::toy_simmap`][toy_simmap], [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation]
#'   and [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges], so the datasets are consistent with all the other toy datasets.
#'
#' @usage data(toy_tip_data)
#' @format A data.frame with 30 rows (the tips of [`phylograftR::toy_backbone`][toy_backbone], as row names and in the same order) and three character columns:
#'   * `$color`: `"blue"`, `"green"` or `"red"`.
#'   * `$elevation`: `"high"` or `"low"`.
#'   * `$range`: biogeographic ranges coded as in `BioGeoBEARS`, i.e., `"A"`, `"B"`, `"AB"` or `"BC"`.
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_multiSimmap`][toy_multiSimmap], [`phylograftR::toy_multiSimmap_elevation`][toy_multiSimmap_elevation] and
#'    [`phylograftR::toy_multiSimmap_ranges`][toy_multiSimmap_ranges] for the stochastic maps inferred from these data.
#'  * [`phylograftR::toy_backbone`][toy_backbone] for the backbone phylogeny.
#'
#' @examples
#' # Import the data
#' data(toy_tip_data)
#' head(toy_tip_data)
#'
#' # Named vector of tip states of one attribute, as required by deepSTRAPP::prepare_trait_data()
#' tip_colors <- setNames(toy_tip_data$color, row.names(toy_tip_data))
#' table(tip_colors)
#'
"toy_tip_data"

#' @title Toy multiSimmap of a categorical trait: stochastic maps of the colors
#'
#' @description Set of 10 stochastic character maps (`multiSimmap`) of [`phylograftR::toy_backbone`][toy_backbone] for the **colors** of
#'   [`phylograftR::toy_tip_data`][toy_tip_data] (`"blue"`, `"green"` and `"red"`). Each map is **one possible history** of the color along the branches,
#'   so that together they account for the **uncertainty** in the reconstruction of the ancestral states, unlike a single `simmap`.
#'   All the maps have the **same states at the tips** (the observed ones) and differ in the internal branches.
#'
#'   They were inferred from the tip data with `deepSTRAPP::prepare_trait_data()` (see `vignette("prepare_stochastic_maps")`):
#'   the models of equal rates (`"ER"`), symmetrical rates (`"SYM"`) and all rates different (`"ARD"`) were compared with their AICc.
#'   The **best model** was the equal-rates model (`"ER"`), under which the maps were simulated with `phytools::make.simmap()` (seed 1234).
#'   100 maps were simulated, and only the first 10 are kept to keep the package light. 10 maps are enough for illustration,
#'   but **actual analyses should use more maps (100 or more)**.
#'
#' @usage data(toy_multiSimmap)
#' @format An object of classes `"multiSimmap"` and `"multiPhylo"`: a list of 10 objects of class `"simmap"` (see [`phylograftR::toy_simmap`][toy_simmap] for their structure),
#'   with the same tips, topology and branch lengths as [`phylograftR::toy_backbone`][toy_backbone].
#'   The class `"multiSimmap"` is **what tells the grafting functions that the maps are alternative histories of the same trait** (see `vignette("prepare_stochastic_maps")`).
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_tip_data`][toy_tip_data] for the tip data used to build the maps.
#'  * [`phylograftR::toy_simmap`][toy_simmap] for a single map of the same states.
#'  * [`phylograftR::toy_multiSimmap_elevation`][toy_multiSimmap_elevation] to combine with these maps as the maps of two traits.
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa with the maps as constraints: each map is used for a different simulation.
#'
#' @examples
#' # Import the data
#' data(toy_multiSimmap)
#' class(toy_multiSimmap)
#' length(toy_multiSimmap)
#'
#' # Plot the first two maps: two possible histories of the same trait
#' state_colors <- c(blue = "dodgerblue", green = "limegreen", red = "firebrick")
#' par(mfrow = c(1, 2))
#' phytools::plotSimmap(toy_multiSimmap[[1]], colors = state_colors, fsize = 0.6, ftype = "i")
#' phytools::plotSimmap(toy_multiSimmap[[2]], colors = state_colors, fsize = 0.6, ftype = "i")
#' par(mfrow = c(1, 1))
#'
"toy_multiSimmap"

#' @title Toy multiSimmap of a binary trait: stochastic maps of the elevation
#'
#' @description Set of 10 stochastic character maps (`multiSimmap`) of [`phylograftR::toy_backbone`][toy_backbone] for the **elevation** of
#'   [`phylograftR::toy_tip_data`][toy_tip_data] (`"high"` and `"low"`), inferred like [`phylograftR::toy_multiSimmap`][toy_multiSimmap]
#'   (models `"ER"` and `"ARD"` compared with their AICc; the **best model** was the equal-rates model `"ER"`; seed 2345; 100 maps simulated, the first 10 kept).
#'
#'   It has the same number of maps as [`phylograftR::toy_multiSimmap`][toy_multiSimmap], so that the two can be **combined as the maps of two traits**:
#'   the maps of the traits are paired (map 1 of each trait, then map 2 of each trait, and so on).
#'
#' @usage data(toy_multiSimmap_elevation)
#' @format An object of classes `"multiSimmap"` and `"multiPhylo"`: a list of 10 objects of class `"simmap"`,
#'   with the same tips, topology and branch lengths as [`phylograftR::toy_backbone`][toy_backbone].
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_tip_data`][toy_tip_data] for the tip data used to build the maps.
#'  * [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation] for a single map of the same states.
#'  * [`phylograftR::toy_multiSimmap`][toy_multiSimmap] for the maps of the colors.
#'
#' @examples
#' # Import the data
#' data(toy_multiSimmap)
#' data(toy_multiSimmap_elevation)
#'
#' # Maps of two traits: a list with one element per trait, each a multiSimmap with the same number of maps
#' maps_two_traits <- list(color = toy_multiSimmap, elevation = toy_multiSimmap_elevation)
#' sapply(maps_two_traits, length)
#'
"toy_multiSimmap_elevation"

#' @title Toy multiSimmap of biogeographic ranges
#'
#' @description Set of 10 stochastic character maps (`multiSimmap`) of [`phylograftR::toy_backbone`][toy_backbone] for the **biogeographic ranges** of
#'   [`phylograftR::toy_tip_data`][toy_tip_data], coded as in `BioGeoBEARS` (areas `"A"`, `"B"` and `"C"`; ranges of at most two areas).
#'   Each map is one possible biogeographic history, with anagenetic changes of range along the branches and cladogenetic events at the nodes.
#'
#'   They were inferred from the tip data with the **DEC** model of `BioGeoBEARS` (the only model fit), through `deepSTRAPP::prepare_trait_data()` (see `vignette("prepare_stochastic_maps")`):
#'   biogeographic stochastic mapping (BSM) was run, and its output converted into stochastic maps (seed 3456).
#'   100 maps were simulated, and only the first 10 are kept to keep the package light. 10 maps are enough for illustration,
#'   but **actual analyses should use more maps (100 or more)**.
#'
#'   The ranges of the maps include the ranges of the tips and may include **other ranges only found at the internal nodes** (e.g., `"C"` or `"AC"`).
#'   Use `treat_states_as_ranges = TRUE` to use the areas as allowed states in the grafting functions (see `vignette("biogeographic_data")`).
#'
#' @usage data(toy_multiSimmap_ranges)
#' @format An object of classes `"multiSimmap"` and `"multiPhylo"`: a list of 10 objects of class `"simmap"`,
#'   with the same tips, topology and branch lengths as [`phylograftR::toy_backbone`][toy_backbone].
#'
#' @seealso Other datasets and functions from phylograftR:
#'  * [`phylograftR::toy_tip_data`][toy_tip_data] for the tip data used to build the maps.
#'  * [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges] for a single map of ranges.
#'  * [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`][toy_missing_taxa_list_with_allowed_ranges] for missing taxa with allowed ranges.
#'
#' @examples
#' # Import the data
#' data(toy_multiSimmap_ranges)
#'
#' # Ranges found in the first map
#' colnames(toy_multiSimmap_ranges[[1]]$mapped.edge)
#'
"toy_multiSimmap_ranges"
