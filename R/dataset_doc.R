
##### Documentation of the datasets #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Includes:
#  - A backbone phylogeny for grafting
#  - A list of missing taxa to be grafted, with and without allowed states or allowed ranges
#  - A table of total richness per clades to inform sampling fractions with unsampled species
#  - A stochastic character map (simmap) of the backbone phylogeny, mapping a trait with three states
#  - A stochastic character map (simmap) of the backbone phylogeny, mapping biogeographic ranges with four states


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
#' plot(grafted_phy$simmap, fsize = 0.5, ftype = "i",
#'      colors = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))
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
#'   Such a list requires a `simmap` to graft.
#'
#' @usage data(toy_missing_taxa_list_with_allowed_ranges)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of character vectors of ranges of [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges]
#'    in which the taxon can be grafted, or `NULL` to use the global `allowed_states` argument. `"Genus_a_sp1"`: `"A"`;
#'    `"Genus_b_sp1"`: `"B"` and `"BC"`; `"Genus_b_sp2"`: `"B"`; `"Genus_c_sp2"`: `"A"` and `"BC"`; the others: `NULL`.
#'
#' @seealso [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges], [`phylograftR::toy_missing_taxa_list_with_allowed_states`][toy_missing_taxa_list_with_allowed_states],
#'   [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list]
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
#' plot(grafted_phy$simmap,
#'      col = c("A" = "gold", "B" = "dodgerblue", "AB" = "limegreen", "BC" = "firebrick2"))
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


##### 4b/ Toy simmap of biogeographic ranges #####

#' @title Toy stochastic character map of biogeographic ranges of the toy backbone phylogeny
#'
#' @description Stochastic character map (simmap) of [`phylograftR::toy_backbone`][toy_backbone] with four states that are biogeographic ranges
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
#' Anagenetic changes (along the branches, 12 in total):
#' * 8 `"d"` (range expansion, gain of one area): 5 `"A"` -> `"AB"` and 3 `"B"` -> `"BC"`.
#' * 4 `"e"` (range contraction, loss of one area): 3 `"AB"` -> `"A"` and 1 `"AB"` -> `"B"`.
#'
#' Cladogenetic events (at the 29 nodes, from the range of the lineage above the node and the ranges at the start of the two descendant branches):
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
#' @seealso [`phylograftR::toy_simmap`][toy_simmap] for a simmap with three states that are not ranges,
#'   [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`][toy_missing_taxa_list_with_allowed_ranges], [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]
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


##### 5/ Toy datasets for a second trait: simmap and list of missing taxa with allowed states of elevation #####

#' @title Toy stochastic character map of elevation states
#'
#' @description Stochastic character map (simmap) of [`phylograftR::toy_backbone`][toy_backbone] with two elevation states, `"high"` and `"low"`,
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
#' @seealso Associated missing_taxa_list object to combine with colors: [`phylograftR::toy_missing_taxa_list_with_two_traits`][toy_missing_taxa_list_with_two_traits].
#' Simmap of colors: [`phylograftR::toy_simmap`][toy_simmap].
#' Grafting function: [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#'
#' @examples
#' # Import the data
#' data(toy_simmap_elevation)
#'
#' # Plot the map
#' plot(toy_simmap_elevation, colors = c(high = "darkorange", low = "steelblue"), fsize = 0.7, ftype = "i")
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
#'   In the column `allowed_states`, the states of the traits are separated by semicolons
#'   (in the order of the maps provided in `simmap`: the color, then the elevation), and the states of a trait by commas. Open that file to see the different syntaxes:
#'   * `"red;"`: the color is `"red"`, and the elevation has no taxon-level value (the global value is used).
#'   * `"blue, green;high"`: color `"blue"` or `"green"`, and elevation `"high"`.
#'   * `";low"`: no taxon-level value for the color, elevation `"low"`.
#'   * `"elevation=high"`: traits named, so the color can be omitted (it has no taxon-level value).
#'   * Blank cell: no taxon-level value for any trait, so the global `allowed_states` argument is used for both.
#'
#' @usage data(toy_missing_taxa_list_with_two_traits)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of lists with the states of each trait (one element per trait, `NULL` for a trait without taxon-level value),
#'    or `NULL` to use the global `allowed_states` argument for all the traits. `"Genus_a_sp1"`: `list("red", NULL)`; `"Genus_b_sp1"`: `list(c("blue", "green"), "high")`;
#'    `"Genus_b_sp2"`: `list(NULL, "low")`; `"Genus_c_sp2"`: `list(elevation = "high")`; the others: `NULL`.
#'
#' @seealso Simmap for colors: [`phylograftR::toy_simmap`][toy_simmap]. Simmap for elevation: [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation].
#'   Missing taxa list for colors only: [`phylograftR::toy_missing_taxa_list_with_allowed_states`][toy_missing_taxa_list_with_allowed_states].
#'   Grafting function: [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list].
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
#'   In the column `allowed_states`, the states of the traits are separated by semicolons
#'   (in the order of the maps provided in `simmap`: the ranges, then the elevation), and the states of a trait by commas. Open that file to see the different syntaxes:
#'   * `"A;high"`: the range is `"A"`, and the elevation is `"high"`.
#'   * `"B, BC;"`: the range is `"B"` or `"BC"`, and the elevation has no taxon-level value (the global value is used).
#'   * `"elevation=low"`: traits named, so the ranges can be omitted (they have no taxon-level value).
#'   * `"biogeo=A, BC; elevation=high"`: all the traits named, in any order.
#'   * Blank cell: no taxon-level value for any trait, so the global `allowed_states` argument is used for both.
#'
#'   The names of the traits (`biogeo` and `elevation`) must match the names of the maps provided in `simmap` to use the syntax with names.
#'
#' @usage data(toy_missing_taxa_list_with_ranges_and_elevation)
#' @format A named list with the same 7 taxa and items as [`phylograftR::toy_missing_taxa_list`][toy_missing_taxa_list], and the extra item:
#'  * `$allowed_states` List (one element per taxon, named after the taxa) of lists with the states of each trait (one element per trait, `NULL` for a trait without taxon-level value),
#'    or `NULL` to use the global `allowed_states` argument for all the traits. `"Genus_a_sp1"`: `list("A", "high")`; `"Genus_b_sp1"`: `list(c("B", "BC"), NULL)`;
#'    `"Genus_b_sp2"`: `list(elevation = "low")`; `"Genus_c_sp2"`: `list(biogeo = c("A", "BC"), elevation = "high")`; the others: `NULL`.
#'
#' @seealso Simmap for ranges: [`phylograftR::toy_simmap_ranges`][toy_simmap_ranges]. Simmap for elevation: [`phylograftR::toy_simmap_elevation`][toy_simmap_elevation].
#'   Missing taxa list with ranges alone: [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`][toy_missing_taxa_list_with_allowed_ranges].
#'   Missing taxa list with two non-biogeographic traits: [`phylograftR::toy_missing_taxa_list_with_two_traits`][toy_missing_taxa_list_with_two_traits].
#'   Import function: [`phylograftR::import_missing_taxa_list()`][import_missing_taxa_list].
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
#' # Taxa without taxon-level states can be grafted on branches in a range containing "A" or "B" (global value),
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
#' check_grafting_constraints(grafted_phy,
#'     missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
#'     allowed_states = list(c("A", "B"), NULL),
#'     treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE))
#'
"toy_missing_taxa_list_with_ranges_and_elevation"
