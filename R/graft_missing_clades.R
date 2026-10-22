##### Function to graft missing clades #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


#' @title Graft missing clades on a time-calibrated phylogeny
#'
#' @description This function grafts missing taxa onto an ultrametric, time-calibrated backbone phylogeny as clades sister to clades of the backbone.
#'   Missing taxa sharing the same sister clade are grafted together as ONE missing clade, on the stem branch of their sister clade.
#'
#'   The divergence time of the missing clade is drawn from a birth-death model fitted on the sister clade,
#'   and the missing clade is grown under the same model, conditioned on its number of taxa and its stem age.
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which taxa are grafted.
#'   It must be rooted, binary (no polytomies) and ultrametric (i.e., be time-calibrated).
#' @param missing_taxa_list A list with the information on the taxa to graft. Items:
#'   * `$Taxon_name`: Character vector of the names of the taxa to graft. Names must be unique and absent from `phy$tip.label`.
#'   * `$Sister_clade`: List (named after `$Taxon_name`) of character vectors. For each taxon, the tip labels defining the clade of the backbone
#'     that is sister to the missing clade. Taxa with the same sister clade are grafted together as one missing clade.
#'     A pair of tips defines the clade as their MRCA node. A single tip defines the "clade" as the terminal branch of this tip.
#'     Clades, including the excluded subclades (see below), must be defined by tips: node IDs are not accepted (the function stops, and lists some of them).
#'     Use [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] to convert node IDs into pairs of tips.
#'   * `$allowed_states`: (Optional) List (named after `$Taxon_name`) of character vectors. For each taxon, the states of the `simmap` in which its missing clade
#'     can be grafted on the stem branch of its sister clade. It overrides the global argument `allowed_states` for this taxon. `NULL` or `NA` use the global argument.
#'     Only used with a `simmap`. Taxa with the same sister clade are grafted together as one missing clade, so they must have the same allowed states
#'     (the function stops otherwise, and lists the clades concerned).
#' @param simmap (Optional) A single stochastic character map of the backbone, of class `"simmap"` (e.g., from `phytools::make.simmap()` or `phytools::sim.history()`),
#'   mapping a species attribute (e.g., a trait or a range) with discrete states along the branches of `phy`.
#'   It must have the same tips, topology and branch lengths as `phy`. With a `simmap`, the missing clades can only be grafted on the parts of the stem branches
#'   of their sister clades that are in their allowed states (see `allowed_states`), and each missing clade inherits the state of the point of the stem branch
#'   where it is grafted, along all its branches, so the states can be followed alongside the phylogeny.
#'
#'   See 'Species attributes' in Details section. Default is `NULL` (no restriction for grafting based on states).
#' @param allowed_states Character vector. The states of `simmap` in which grafting is allowed. This is the default for all taxa:
#'   taxa with a value in `missing_taxa_list$allowed_states` use this specific value instead. States are matched exactly with the states of the `simmap`.
#'   Taxa with no allowed states (neither a global nor a taxon-level value) can be grafted in all the states of the `simmap`.
#'   The function stops if some states are not in the `simmap`, or if `allowed_states` is provided without a `simmap`.
#'   Default is `NULL` (all states allowed).
#' @param treat_states_as_ranges Logical. Only used with a `simmap`. Whether the states of the `simmap` are biogeographic ranges to be treated as sets of areas, rather than as unrelated states.
#'   Ranges must be coded as in `BioGeoBEARS`: one capital letter per area, and the letters of a multi-area range in alphabetical order (e.g., `"A"`, `"AB"`, `"ABC"`).
#'   If `TRUE`, an allowed range enables grafting on branch pieces in any range that contains it: `"A"` can be grafted on `"A"`, `"AB"` and `"ABC"`,
#'   whereas `"AB"` can be grafted on `"AB"` and `"ABC"` but not `"A"`.
#'
#'   See 'Biogeographic ranges as states' in Details section. Default is `FALSE` (states are matched exactly).
#' @param return_simmap Logical. Only used with a `simmap`, and silently ignored without it. Whether to add the updated `simmap` to the output, in `$simmap` (see 'Value').
#'   Default is `TRUE`.
#' @param min_valid_fraction Numeric between 0 and 1. Only used with a `simmap`. Threshold of the fraction of the time window of the divergence time of a missing clade
#'   (the stem branch of its sister clade) that is in its allowed states. A warning lists the missing clades below the threshold for inspection,
#'   as their divergence times are concentrated on a small part of the stem branch. The fraction is computed on the backbone, before grafting. Default is `0.1`.
#' @param max_age_redraws Integer. Only used with a `simmap`. Maximum number of times the divergence time of a missing clade can be drawn again when it falls
#'   in a part of the stem branch that is not in an allowed state. The function stops with an informative message if this limit is reached. Default is `1000`.
#' @param min_tips Integer. Minimum number of tips of the clade used to fit the birth-death model.
#'   If the sister clade has fewer tips, the smallest enclosing clade with at least `min_tips` tips is used instead. Must be at least 3. Default is `5`.
#' @param total_richness_per_clades (Optional) List with the total number of species of the clades, to account for species that are neither in the backbone
#'   nor in `missing_taxa_list` (i.e., unsampled species) when fitting the birth-death models. Items:
#'   * `$Sister_clade`: List of character vectors of tips, defining clades as in `missing_taxa_list$Sister_clade`.
#'     Node IDs are not accepted: convert them with [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] in the same call as `missing_taxa_list`,
#'     so the tips of the clades match between both objects.
#'     It must cover all the unique sister clades of `missing_taxa_list`. It can include other clades, e.g., enclosing clades.
#'   * `$total_richness`: Vector of positive integers, same length as `$Sister_clade`. Total number of species of each clade,
#'     i.e., backbone species, missing species nested within the clade, and unsampled species. Use `NA` if unknown: the unsampled species
#'     are then assumed to be those of the clades with an informed total richness nested within the clade (none if there are no such clades).
#'
#'   With `NULL` (default), no `sampling_fraction_df` is returned and missing species to graft and backbone species
#'   are assumed to cover all species in the focal group (no unsampled species).
#'   If the list is provided with only `NA` in `$total_richness`, the fits are identical (no unsampled species in any clade, `total_richness` = backbone + missing species),
#'   but the `sampling_fraction_df` table is returned (with `total_richness_provided = FALSE`).
#'
#'   See 'Sampling fraction and unsampled species' in Details section.
#' @param ultrametric_tol Numeric. Tolerance used by `ape::is.ultrametric()` for the ultrametricity checks of the backbone and of the final phylogeny.
#'   Default is `.Machine$double.eps^0.5`.
#' @param verbose Logical. Whether to display progress. Messages print every 10 missing clades. Default is `TRUE`.
#' @param save_temp Logical. Whether to save a temporary version of the grafted phylogeny during the process, to recover from a crash on long runs. Default is `FALSE`.
#' @param save_temp_every Integer. Frequency of the temporary saves, in number of missing clades grafted. Only used if `save_temp = TRUE`. Default is `10`.
#' @param temp_file_path Character string. Path of the `.rds` file used for temporary saves. It must end with `.rds` and its directory must exist and be writable. Default is `NULL`, which saves in `"./temp_phy_from_grafting.rds"`.
#' @param seed Integer. Seed used for reproducibility. The seed only applies within the function: the state of the random number generator of the R session
#'   is restored on exit. Set to `NULL` to use the random number generator of the session. Default is `1234`.
#'
#' @return The function returns an object of class `"phylo"`: the backbone phylogeny with the grafted taxa added as new tips.
#'   Backbone branch lengths are untouched, and the grafted phylogeny is ultrametric.
#'
#'   In addition to the standard elements of a `"phylo"` object, the phylogeny has the following elements:
#'   * `grafted_tips` Named logical vector with one value per tip, in the order of `$tip.label`. `TRUE` for grafted tips. Names are the tip labels.
#'   * `grafted_edges` Named logical vector with one value per branch, in the order of the rows of `$edge`. `TRUE` for grafted branches,
#'     i.e., branches with only grafted taxa descending from them. Names are the IDs of the tipward nodes of the branches (`$edge[, 2]`).
#'     Use them to color the phylogeny on plots (see examples).
#'
#'   These elements are not updated by functions that modify the phylogeny (e.g., `ape::drop.tip()`, which removes tips but keeps the other elements as they are,
#'   and `ape::write.tree()`, which removes them). Use [`phylograftR::identify_grafted()`][identify_grafted] to recompute `grafted_tips` and `grafted_edges` on a manipulated phylogeny.
#'
#'   Optional elements can be added to the output based on what optional argument was provided:
#'
#'   * `sampling_fraction_df` Only if `total_richness_per_clades` is provided (a list with only `NA` values also returns it, assuming no unsampled species).
#'     Data frame with one row per clade of `total_richness_per_clades`:
#'     * `Clade_tips` Character string. Tips used to identify the clade.
#'     * `total_richness` Integer. Derived when `NA` was provided.
#'     * `total_richness_provided` Logical. Whether `total_richness` of the clade was provided as input in `total_richness_per_clades`
#'     * `nb_backbone_species` Integer. Number of species initially present on the backbone phylogeny.
#'     * `nb_missing_species` Integer. Number of species that were missing from the backbone phylogeny and are now included in the grafted phylogeny
#'     * `nb_unsampled_species` Integer. Number of species that are neither on the backbone, nor grafted as missing. Derived from `total_richness`.
#'     * `sampling_fraction` Numeric. Fraction of the clade richness present in the grafted phylogeny as `(nb_backbone_species + nb_missing_species) / total_richness`.
#'       This can be used to account for clade-specific sampling fraction in downstream analyses.
#'   * `simmap` Only if a `simmap` is provided and `return_simmap = TRUE`. The updated stochastic character map of the grafted phylogeny (class `"simmap"`):
#'     the pieces of the branches of the backbone keep their states (cut where taxa were grafted), and every grafted clades is entirely in the state of the point where it was grafted.
#'     This element can be plotted as any `simmap`, e.g., `plot(grafted_phy$simmap)`, or used to extract the states of the tips with `phytools::getStates(grafted_phy$simmap, type = "tips")`.
#'     It has the same exact topology and branch lengths as the grafted phylogeny.
#'
#' @details
#' ### Stop conditions
#'
#' The function stops, before any grafting, if the backbone is not ultrametric, or if a sister clade is the whole tree (it has no stem branch to graft on).
#'
#' ### Grafting process
#'
#' Groups of missing taxa sharing the same sister clade are ordered by the size of the sister clade, so nested sister clades are processed first.
#' For each group:
#' 1. A birth-death model is fitted with `phytools::fit.bd()` on the sister clade, or on the smallest enclosing clade with at least `min_tips` tips.
#'    The model is fitted on the crown group of this clade, with the sampling fraction of this crown group (see Sampling fraction and unsampled species).
#' 2. The divergence time of the missing clade is drawn on the stem branch of the sister clade with `TreeSim::corsim()`:
#'    between the crown age of the sister clade and its stem age (or between the present and the stem age for a single-tip sister clade).
#' 3. A single missing taxon is grafted as a new tip at this divergence time. For several taxa, a missing clade is grown with `TreeSim::sim.bd.taxa.age()`
#'    under the fitted model, conditioned on the number of taxa and on its stem age (the drawn divergence time), and grafted on the stem branch of the sister clade.
#'
#' ### Species attributes: allowed states of a simmap
#'
#' A `simmap` describes the states of a discrete attribute (e.g., a trait or a range) along the branches of the backbone: each branch is cut in pieces, each piece being in one state.
#' Missing clades are grafted on the stem branch of a clade of the backbone, and the function assumes that the user considers this stem branch a valid place for them.
#' Restricting the grafting to some states of the stem branch is therefore only useful when the user has a strong hypothesis about a change of state
#' along the stem branch that should affect where the missing clade can be grafted. For most uses, the `simmap` is not needed.
#'
#' For example, in biogeography: a stem fossil is found in range A, and the sister clade is in range A, while the missing clade is in range B.
#' If the shift from range B to range A on the stem branch of the sister clade is estimated to have occurred before the age of the fossil,
#' the missing clade can be restricted to the part of the stem branch before the shift, with `allowed_states = "B"`.
#'
#' With a `simmap`, a missing clade is only grafted on the pieces of the stem branch of its sister clade that are in its allowed states
#' (`missing_taxa_list$allowed_states` of its taxa if provided, otherwise the global argument `allowed_states`).
#' The divergence time of the missing clade is drawn with `TreeSim::corsim()` as without a `simmap`: if it falls in a part of the stem branch that is not in an allowed state,
#' it is drawn again, up to `max_age_redraws` times.
#'
#' The `simmap` can also be provided so that the states are updated alongside the phylogeny: the whole missing clade (all its taxa and branches)
#' inherits the state of the point of the stem branch where it is grafted. The updated `simmap` is returned in `$simmap`, if `return_simmap = TRUE`).
#' The allowed states are not required: without any allowed states, all the states are allowed, grafting is not restricted by the states,
#' but the simmap will be updated, assuming grafted clades inherit the state of the point where they are grafted.
#'
#' The function stops, before any grafting, if the taxa of a same missing clade have different allowed states, or if the stem branch of a sister clade has no piece in the allowed states
#' (the clades are listed), and warns about the missing clades with a small share of their stem branch in allowed states (see `min_valid_fraction`).
#'
#' ### Biogeographic ranges as states
#'
#' With `treat_states_as_ranges = TRUE`, the states of the `simmap` are biogeographic ranges, coded as in `BioGeoBEARS`
#' (one capital letter per area, and the letters of a multi-area range in alphabetical order, e.g., `"AB"`),
#' and the allowed states are the ranges a missing clade can have.
#'
#' An allowed range matches a state of the `simmap` if it is contained in it: a missing clade allowed in `"A"` can be grafted on the parts of the stem branch in `"A"`, `"AB"` or `"ABC"`,
#' whereas a missing clade allowed in `"AB"` only can be grafted on the parts of the stem branch in `"AB"` or `"ABC"`.
#'
#' When multiple allowed states are provided, the inherited range is drawn randomly from the allowed ranges enclosed within the range of the stem branch piece where it is grafted.
#' When grafting on an edge with a larger range than the inherited range (e.g., `"A"` on `"AB"`), this results in a subset sympatry event such as (AB) -> (AB), (A).
#'
#' ### Enclosing clades
#'
#' When a sister clade contains fewer tips than `min_tips`, an alternative higher clade is selected to fit the birth-death model.
#' This 'enclosing clade' is identified by walking rootward the tree node by node, starting from the sister clade, until a clade meets the `min_tips` condition.
#' This ensures that the diversification dynamics modeled are based on enough tips to be relevant, while being closely associated with the sister clade.
#'
#' ### Sampling fraction
#'
#' The birth-death model is fitted on the crown group of the sister clade or its enclosing clade, so the sampling fraction is:
#' `n_present / (n_present + n_pending + n_unsampled)`, with `n_present` the tips of the fitted clade in the current tree
#' (backbone tips and taxa already grafted, so missing clades grafted in nested sister clades processed earlier are not missing anymore),
#' `n_pending` the taxa of `missing_taxa_list` still to graft whose sister clade is strictly nested within the fitted clade,
#' and `n_unsampled` the unsampled species derived from `total_richness_per_clades` (0 without it).
#'
#' Taxa that are sister to the fitted clade itself are outside of its crown group, so they are not counted:
#' when the sister clade itself is fitted (i.e. have at least `min_tips` tips) and nothing is pending nor unsampled within it, the sampling fraction is 1.
#' When an enclosing clade is fitted, the taxa to graft on the sister clade (nested in the enclosing clade) are counted.
#' The divergence time of each missing clade is still drawn with `TreeSim::corsim(missing = 1)`: unsampled species only influence the fitted rates,
#' they do not take a branch.
#'
#' ### Unsampled species
#'
#' For a clade with a `total_richness` value provided: `nb_unsampled = total_richness - nb_backbone - nb_missing`, with `nb_missing` the number of taxa to graft
#' whose sister clade is strictly nested within the clade. For a clade with no `total_richness` value (`NA`, or enclosing a sister clade without being listed),
#' `nb_unsampled` is the sum of `nb_unsampled` over the largest clades nested within it with a `total_richness` value provided. If there is none, `nb_unsampled = 0`.
#' The rationale is that even if `nb_unsampled` cannot be derived from its `total_richness` (as it is not provided for this clade), it has to at least account for the `nb_unsampled`
#' species that have being identified in its own nested clades.
#'
#' If `total_richness_per_clades = NULL`, `nb_unsampled = 0` and no `sampling_fraction_df` is returned.
#' To obtain the `sampling_fraction_df` without informing any total richness, provide `total_richness_per_clades` with the clades and only `NA` in `$total_richness`:
#' `nb_unsampled = 0` and `total_richness` equals the backbone plus missing species in all clades (the sampling fraction is 1).
#'
#' The function stops, before any grafting, if `total_richness_per_clades` does not cover all the sister clades,
#' if a `total_richness` provided is lower than the number of backbone species plus missing species of the clade,
#' or if a clade has fewer unsampled species than the sum of the clades with an informed `total_richness` nested within it.
#'
#'
#' @importFrom ape bind.tree branching.times extract.clade is.ultrametric
#'
#' @importFrom phytools bind.tip nodeHeights nodeheight
#'
#' @importFrom TreeSim corsim sim.bd.taxa.age
#'
#' @export
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::identify_grafted()`][identify_grafted] to identify the grafted tips and branches of a phylogeny,
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade with constraints (target clades and excluded subclades),
#'  * [`phylograftR::check_grafting_constraints()`][check_grafting_constraints] to check a posteriori that grafted taxa respect their constraints.
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' # List the missing taxa to graft
#' missing_clades_list <- list(
#'   Taxon_name = c("M1", "M2", "M3", "M4"),
#'   # Sister_clade is named after Taxon_name (same names, same order)
#'   # M1 and M2 share the same sister clade: they are grafted together as one missing clade
#'   Sister_clade = list(M1 = c("t8", "t22"), M2 = c("t8", "t22"), M3 = c("t26", "t18"), M4 = "t5"))
#'
#' # --- Perform grafting --- #
#'
#' # Graft the missing clades
#' grafted_phy <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
#'                                     verbose = FALSE, seed = 1234)
#' plot(grafted_phy, cex = 0.7)
#'
#' # Color the grafted tips and branches in red based on elements returned with the phylogeny
#' plot(grafted_phy, cex = 0.7,
#'      edge.color = c("black", "red")[grafted_phy$grafted_edges + 1],
#'      tip.color = c("black", "red")[grafted_phy$grafted_tips + 1])
#'
#' # --- Account for unsampled species --- #
#'
#' # Account for unsampled species with the total richness of the sister clades.
#' # All the unique sister clades must be listed (NA = unknown richness)
#' total_richness_per_clades <- list(
#'   Sister_clade = list(c("t8", "t22"), c("t26", "t18"), "t5"),
#'   total_richness = c(20, NA, 3))
#' grafted_phy_richness <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
#'                                              total_richness_per_clades = total_richness_per_clades,
#'                                              verbose = FALSE, seed = 1234)
#' grafted_phy_richness$sampling_fraction_df
#'
#' # --- Apply attribute-based exclusions --- #
#'
#' # Restrict the grafting to some states of a stochastic character map (simmap) of the backbone
#' # The toy simmap has three states, "blue", "green" and "red"
#' data(toy_backbone)
#' data(toy_simmap)
#'
#' # List the missing taxa to graft
#' missing_clades_list_states <- list(
#'    Taxon_name = c("M1", "M2", "M3", "M4", "M5"),
#'    Sister_clade = list(M1 = c("t26", "t18"), M2 = c("t26", "t18"), M3 = c("t26", "t18"),
#'                        M4 = "t14", M5 = "t14"))
#'
#' # The simmap shows that the targeted stem branches recorded state shifts.
#' plot(toy_simmap, col = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))
#'
#' # We can use 'allowed_states' to limit grafting on the "green" portions of the stem branches
#'    # The stem branch from the "bottom" clade shifts from "green" to "red"
#'    # The terminal (stem) branch of "t14" shift from "green" to "red" too
#' grafted_phy_states <- graft_missing_clades(phy = toy_backbone,
#'                                            missing_taxa_list = missing_clades_list_states,
#'                                            simmap = toy_simmap, allowed_states = "green",
#'                                            verbose = FALSE, seed = 1234)
#'
#' # Plot output phylogeny with colored grafts
#' plot(grafted_phy_states, cex = 0.7,
#'      edge.color = c("black", "red")[grafted_phy_states$grafted_edges + 1],
#'      tip.color = c("black", "red")[grafted_phy_states$grafted_tips + 1])
#'
#' # The updated simmap is in the output: the new grafted clades branch from the "green" segments
#' # of the stem branches of their sister clades, and each newly grafted clades are entirely "green"
#' plot(grafted_phy_states$simmap, col = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))
#'
#'


graft_missing_clades <- function (phy, # Backbone phylogeny
                                  missing_taxa_list, # List with Taxon_name = c(), and items named after Taxon_name: Sister_clade = list(X1 = c(), ...), allowed_states = list(X1 = c(), ...) (optional)
                                  simmap = NULL, # NULL, or stochastic character map (class "simmap") of the backbone: restricts grafting to the parts of the stem branches in allowed states
                                  allowed_states = NULL, # States of the simmap where grafting is allowed. Default for all taxa, overridden by missing_taxa_list$allowed_states. NULL: all states
                                  treat_states_as_ranges = FALSE, # Are the states of the simmap biogeographic ranges to be treated as sets of areas? An allowed range matches the states that contain it. (only with a simmap)
                                  return_simmap = TRUE, # Add the updated simmap to the output ($simmap)? Only with a simmap
                                  min_valid_fraction = 0.1, # Warn for the clades with a smaller fraction of their stem branch in allowed states (only with a simmap)
                                  max_age_redraws = 1000, # Maximum number of redraws of the divergence time of a missing clade when it falls outside the allowed states (only with a simmap)
                                  min_tips = 5, # Minimum number of tips of the clade used to fit the birth-death model (the enclosing clade is used if the sister clade is smaller)
                                  total_richness_per_clades = NULL, # NULL, or list with Sister_clade = list(c(), ...) and total_richness = c(...) (NA allowed). Total number of species of the clades accounting for the unsampled species that are neither present in the backbone, nor grafted.
                                  ultrametric_tol = .Machine$double.eps^0.5, # Tolerance used by ape::is.ultrametric() for the ultrametricity checks
                                  verbose = TRUE,
                                  save_temp = FALSE, # Save a temporary phylogeny every 'save_temp_every' clades
                                  save_temp_every = 10,
                                  temp_file_path = NULL, # Default = "./temp_phy_from_grafting.rds"
                                  seed = 1234)
{
  ### ---- Check input validity ---- ###
  {
    ## The backbone must be ultrametric (error otherwise). It is never modified: making it ultrametric is left to the user
    check_backbone_ultrametric(phy = phy, tol = ultrametric_tol)

    ## The backbone must be binary (error otherwise): polytomies are not handled, and resolving them is left to the user
    check_backbone_binary(phy = phy)

    ## Check 'min_tips'
    if (!is.numeric(min_tips) || length(min_tips) != 1 || min_tips < 3)
    {
      stop("'min_tips' must be a single number >= 3")
    }
    if (length(phy$tip.label) < min_tips)
    {
      stop("Backbone phylogeny too small to fit a birth-death model (< ", min_tips, " tips).\n",
           "Lower 'min_tips' or provide a larger backbone phylogeny.")
    }

    ## Check the structure of missing_taxa_list: Taxon_name, and Sister_clade (same length and names identical to Taxon_name, in the same order)
    validate_missing_taxa_list(missing_taxa_list, clade_item = "Sister_clade", optional_items = "allowed_states")
    n_missing <- length(missing_taxa_list$Taxon_name)

    ## Species attributes: 'simmap', 'allowed_states', 'min_valid_fraction' and 'max_age_redraws'
    # The states of each taxon are resolved below, once the simmap is validated and the backbone labelled
    use_states <- !is.null(simmap)
    if (!is.numeric(min_valid_fraction) || length(min_valid_fraction) != 1 || is.na(min_valid_fraction) || min_valid_fraction < 0 || min_valid_fraction > 1)
    {
      stop("'min_valid_fraction' must be a single number between 0 and 1")
    }
    if (!is.logical(return_simmap) || length(return_simmap) != 1 || is.na(return_simmap))
    {
      stop("'return_simmap' must be TRUE or FALSE")
    }
    if (!is.logical(treat_states_as_ranges) || length(treat_states_as_ranges) != 1 || is.na(treat_states_as_ranges))
    {
      stop("'treat_states_as_ranges' must be TRUE or FALSE")
    }
    if (!is.numeric(max_age_redraws) || length(max_age_redraws) != 1 || is.na(max_age_redraws) || max_age_redraws < 1 || max_age_redraws != round(max_age_redraws))
    {
      stop("'max_age_redraws' must be a single positive integer")
    }
    if (use_states)
    {
      if (inherits(simmap, "multiSimmap"))
      {
        stop("'simmap' must be a single stochastic character map (class 'simmap'), not a list of maps (class 'multiSimmap').")
      }
      warn_if_simmap_unused(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states, return_simmap = return_simmap)
    } else {
      if (treat_states_as_ranges)
      {
        stop("'treat_states_as_ranges = TRUE' requires a 'simmap', whose states are the biogeographic ranges of the branches of the backbone.")
      }
      states_in_list <- any(vapply(missing_taxa_list$allowed_states, function (x) length(x) > 0 && !all(is.na(x)), logical(1)))
      if (!is.null(allowed_states) || states_in_list)
      {
        stop("'allowed_states' (or 'missing_taxa_list$allowed_states') requires a 'simmap', which gives the states of the branches of the backbone.")
      }
    }

    ## Check that taxa names are truly missing from the backbone
    if (any(missing_taxa_list$Taxon_name %in% phy$tip.label))
      stop("Taxon_name must be absent from the backbone phylogeny. Present: ",
           paste(intersect(missing_taxa_list$Taxon_name, phy$tip.label), collapse = ", "))

    ## Check that all tips used to define sister clades are present in the backbone
    all_sister_clade_tips <- unique(unlist(missing_taxa_list$Sister_clade))
    if (!all(all_sister_clade_tips %in% phy$tip.label))
    {
      stop("Some tip labels provided to locate MRCA nodes of sister clades are missing from the backbone phylogeny\n",
           "Missing tips: ", paste(setdiff(all_sister_clade_tips, phy$tip.label), collapse = ", "))
    }

    ## temp_file_path
    # If provided, must be a single path ending with ".rds", in an existing and writable directory
    check_temp_file_path(temp_file_path = temp_file_path)

    ## Richness of the clades (unsampled species)
    # Validates 'total_richness_per_clades' if provided, and derives the number of unsampled species of each node of the backbone.
    # Always computed: it also provides the taxa still to graft to include in the sampling fraction of the BD models
    richness_info <- build_richness_info(phy = phy, missing_taxa_names = missing_taxa_list$Taxon_name,
                                         missing_clades = missing_taxa_list$Sister_clade,
                                         mode = "sister", total_richness_per_clades = total_richness_per_clades)
  }

  ## Set seed if provided, for the duration of this function only (the global RNG state is restored on exit)
  set_local_seed(seed)

  ## Label the internal nodes of the backbone, to identify them in the grafted tree (node IDs change after each graft) and restore the original node labels at the end
  labelled_backbone <- label_backbone_nodes(phy)
  phy <- labelled_backbone$phy

  ## With a simmap: validate it, start following the states of the branches through grafting, and resolve the allowed states of each taxon
  # The resolved list (one character vector per taxon, no NULL) replaces the item in the list, so it follows the taxa when they are reordered
  tracker <- NULL
  if (use_states)
  {
    tracker <- start_state_tracking(phy = phy, simmap = simmap)
    if (treat_states_as_ranges)
    {
      # Allowed ranges are resolved into the acceptable states (the states containing an allowed range), used to restrict grafting, and the candidate ranges of each state
      resolved <- resolve_allowed_ranges(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                         states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
      missing_taxa_list$allowed_states <- resolved$allowed_states # Allowed states adjusted to the ones observed in the simmap
      missing_taxa_list$range_candidates <- resolved$range_candidates # List of ranges that can be inherited by the taxa, based on which allowed state it is grafted on
    } else {
      missing_taxa_list$allowed_states <- resolve_allowed_states(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                                 states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
    }
  }

  ## Randomize the order of taxa so it has no influence on the pattern of grafting
  random_order <- sample.int(n_missing)
  missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[random_order] })

  # Root age: constant during grafting (all clades' tips reach the present) => computed once
  root_age <- max(phytools::nodeHeights(tree = phy)[, 2])

  if (verbose) cat(paste0(Sys.time(), " - Start identification and ordering of sister-clades\n"))

  # Identify sister-clades and order by clade size
  # Since nested clades are always smaller than their higher parents, this is equivalent to order by nesting scheme
  groups <- group_and_order_clades(phy, missing_taxa_list$Sister_clade)
  n_groups <- length(groups$group_clade_tips)

  # Stop before grafting if a sister clade is the whole tree: it has no stem to graft on
  root_node <- setdiff(phy$edge[, 1], phy$edge[, 2])
  is_whole_tree <- vapply(groups$group_clade_tips, function (tips) get_clade_node(phy, tips) == root_node, logical(1))
  if (any(is_whole_tree))
  {
    stop("Sister clade defined by ", paste(groups$group_clade_tips[[which(is_whole_tree)[1]]], collapse = ", "),
         " is the whole tree: it has no stem to graft on")
  }

  ## With a simmap: check before grafting that the stem branch of each sister clade can receive its missing clade
  if (use_states)
  {
    # Taxa grafted together as one missing clade must share their allowed states
    states_of_group <- lapply(seq_len(n_groups), function (g) { vapply(missing_taxa_list$allowed_states[groups$group_ID == g], function (x) paste(sort(unique(x)), collapse = ", "), character(1)) })
    inconsistent <- which(vapply(states_of_group, function (x) length(unique(x)) > 1, logical(1)))
    if (length(inconsistent) > 0)
    {
      n_show <- min(length(inconsistent), 10)
      stop(length(inconsistent), " missing clade(s) have taxa with different allowed states: taxa with the same sister clade are grafted together as one missing clade, so they must share their allowed states.\n",
           paste0("  - sister clade tips: ", vapply(groups$group_clade_tips[inconsistent[seq_len(n_show)]], paste, character(1), collapse = ", "), "; ",
                  vapply(inconsistent[seq_len(n_show)], function (g) { x <- states_of_group[[g]][order(names(states_of_group[[g]]))]; paste0(names(x), " (", x, ")", collapse = ", ") }, character(1)), collapse = "\n"),
           if (length(inconsistent) > n_show) paste0("\n  ... and ", length(inconsistent) - n_show, " more") else "",
           "\nAdjust 'allowed_states' (global or in missing_taxa_list), or give these taxa different sister clades.")
    }

    # The stem branch of the sister clade must have a piece in the allowed states, and the share of its length in these states is its share of the time window of the divergence time
    group_allowed_states <- lapply(seq_len(n_groups), function (g) missing_taxa_list$allowed_states[[which(groups$group_ID == g)[1]]])
    backbone_edge_ages <- root_age - phytools::nodeHeights(phy)
    stem_fraction <- vapply(seq_len(n_groups), function (g)
    {
      stem_node <- get_clade_node(phy, groups$group_clade_tips[[g]])
      stem_length <- phy$edge.length[match(stem_node, phy$edge[, 2])]
      stem_segments <- get_allowed_segments_tracked(tracker = tracker, phy = phy, nodes = stem_node, allowed_states = group_allowed_states[[g]], edge_ages = backbone_edge_ages)
      sum(stem_segments$to - stem_segments$from) / stem_length
    }, numeric(1))
    no_state_groups <- which(stem_fraction == 0)
    if (length(no_state_groups) > 0)
    {
      n_show <- min(length(no_state_groups), 10)
      stop(length(no_state_groups), " missing clade(s) have no piece of the stem branch of their sister clade in their allowed states.\n",
           paste0("  - sister clade tips: ", vapply(groups$group_clade_tips[no_state_groups[seq_len(n_show)]], paste, character(1), collapse = ", "),
                  "; allowed states: ", vapply(group_allowed_states[no_state_groups[seq_len(n_show)]], paste, character(1), collapse = ", "), collapse = "\n"),
           if (length(no_state_groups) > n_show) paste0("\n  ... and ", length(no_state_groups) - n_show, " more") else "",
           "\nAdjust 'allowed_states' (global or in missing_taxa_list) for these taxa.")
    }
    low_fraction <- which(stem_fraction < min_valid_fraction)
    if (length(low_fraction) > 0)
    {
      n_show <- min(length(low_fraction), 10)
      low_groups <- low_fraction[order(stem_fraction[low_fraction])][seq_len(n_show)]
      warning(length(low_fraction), " missing clade(s) have less than ", format(100 * min_valid_fraction), "% of the stem branch of their sister clade in their allowed states ",
              "('min_valid_fraction' = ", min_valid_fraction, "): their divergence times are concentrated on a small part of the stem branch.\n",
              paste0("  - sister clade tips: ", vapply(groups$group_clade_tips[low_groups], paste, character(1), collapse = ", "),
                     "; ", sprintf("%.2f%%", 100 * stem_fraction[low_groups]), " in allowed states", collapse = "\n"),
              if (length(low_fraction) > n_show) paste0("\n  ... and ", length(low_fraction) - n_show, " more") else "",
              call. = FALSE)
    }
  }

  if (verbose) cat(paste0(Sys.time(), " - Start grafting process per missing clades\n"))

  # Initiate number of taxa and missing clades grafted
  nb_taxa_grafted <- 0
  nb_missing_clades_grafted <- 0

  # Loop per ordered clades
  for (i in groups$group_order)
  {
    # Extract ID of taxa in group n°i
    taxa_ID_i <- which(groups$group_ID == i)
    # Extract associated data
    missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[taxa_ID_i] })
    # Extract the tips defining the sister clade
    sister_clade_tips_i <- groups$group_clade_tips[[i]]
    # Get number of taxa to include in this missing clade
    nb_missing_taxa <- length(taxa_ID_i)

    # Identify MRCA node and stem node in grafted tree
    # Node is recomputed from tips because MRCA node IDs change after each graft
    MRCA_node_i <- get_clade_node(phy, sister_clade_tips_i)
    stem_node_i <- get_parent_node(phy, MRCA_node_i)
    is_single_tip <- length(sister_clade_tips_i) == 1

    # Stop for special case where the sister clade is the whole tree
    if (is.na(stem_node_i))
    {
      stop("Sister clade defined by ", paste(sister_clade_tips_i, collapse = ", "), " is the whole tree: it has no stem to graft on")
    }

    # Extract stem age
    stem_age_i <- root_age - phytools::nodeheight(tree = phy, node = stem_node_i)
    upper_age_i <- stem_age_i

    ## Deal with different cases
    if (is_single_tip)
    {
      ## Case for single tips

      # Only its stem branch is available,
      # so only the stem age is provided to TreeSim::corsim()
      divergence_times_i <- stem_age_i
      # Lower bound is set to zero as the target branch is a tip
      lower_age_i <- 0

    } else {

      ## Case for true sister-clades

      sister_clade_tree_i <- ape::extract.clade(phy = phy, node = MRCA_node_i)

      # Get branching times within the sister clade
      divergence_times_i <- ape::branching.times(sister_clade_tree_i)
      # Use the crown age of the sister clade as the lower bound to simulate the divergence age
      # to force the divergence age to be on the stem branch of the sister-clade
      lower_age_i <- max(divergence_times_i)
    }

    ## Fit Birth-Death on the sister-clade, or on the smallest enclosing clade with >= min_tips tips
    # Stop at the root if reached
    # Sampling fraction of the crown group of the fitted clade: accounts for the taxa still to graft and the unsampled species nested in it
    BD_i <- fit_bd_on_clade(phy = phy, node = MRCA_node_i, min_tips = min_tips, richness_info = richness_info)

    # Compute crown age of the sister-clade Age of the sister clade's MRCA node and
    MRCA_age_i <- root_age - phytools::nodeheight(tree = phy, node = MRCA_node_i)
    # Extract stem branch length
    stem_branch_length_i <- phy$edge.length[match(MRCA_node_i, phy$edge[, 2])]

    # With a simmap: pieces of the stem branch in the allowed states of the missing clade
    allowed_states_i <- NULL
    if (use_states)
    {
      allowed_states_i <- missing_taxa_list_i$allowed_states[[1]] # Same states for all the taxa of the missing clade (checked before grafting)
      edge_ages_i <- root_age - phytools::nodeHeights(phy)
      stem_segments_i <- get_allowed_segments_tracked(tracker = tracker, phy = phy, nodes = MRCA_node_i, allowed_states = allowed_states_i, edge_ages = edge_ages_i)
    }

    ## Simulate the missing divergence time with TreeSim::corsim
    # With a simmap, the time is drawn again if it falls in a part of the stem branch that is not in an allowed state:
    # a truncated draw by rejection, which is the exact distribution of the divergence time conditioned on the allowed states
    n_redraws <- 0
    repeat
    {
      new_divergence_times_i <- TreeSim::corsim(
          x = divergence_times_i, lambda = BD_i$b, mu = BD_i$d,
          missing = 1, # The branch leading to the missing clade
          told = upper_age_i, # The stem age
          tyoung = lower_age_i) # The crown age (or 0 for tips)
      # Extract the divergence times for the missing clade
      divergence_time_missing_clade_i <- extract_new_divergence_times(all_times = new_divergence_times_i, old_times = divergence_times_i)

      if (!use_states) { break }
      usable <- nrow(get_segments_crossing_age(phy = phy, segments = stem_segments_i, age = divergence_time_missing_clade_i, root_age = root_age, edge_ages = edge_ages_i)) > 0
      if (usable) { break }

      n_redraws <- n_redraws + 1
      if (n_redraws > max_age_redraws)
      {
        stop("No usable divergence time was found for the missing clade (taxa: ", paste(utils::head(missing_taxa_list_i$Taxon_name, 10), collapse = ", "),
             if (nb_missing_taxa > 10) ", ..." else "", ") on the stem branch of its sister clade (tips: ", paste(sister_clade_tips_i, collapse = ", "),
             ") after ", max_age_redraws, " redraws ('max_age_redraws').\n",
             "Its allowed states: ", paste(allowed_states_i, collapse = ", "), ".\n",
             "Time window of the divergence time (stem branch of the sister clade): from age ", signif(lower_age_i, 4), " to age ", signif(upper_age_i, 4), ". ",
             "Share of this window in the allowed states: ", sprintf("%.2f%%", 100 * sum(stem_segments_i$to - stem_segments_i$from) / stem_branch_length_i), ".\n",
             "The parts of the stem branch in these states are too small compared to the time window of the divergence time.\n",
             "Increase 'max_age_redraws', or extend the allowed states.")
      }
    }

    # Compute position of the graft on the stem branch
    position_i <- min(max(0, divergence_time_missing_clade_i - MRCA_age_i), stem_branch_length_i)

    ## Graft missing clades
    if (nb_missing_taxa == 1)
    {
      ## Case where the missing 'clade' is a single tip

      # Graft the unique tip
      phy_before_graft <- phy
      phy <- phytools::bind.tip(tree = phy, tip.label = missing_taxa_list_i$Taxon_name,
                                edge.length = divergence_time_missing_clade_i,
                                where = MRCA_node_i, position = position_i)
      # With a simmap: the new taxon takes the state of the point where it was grafted, and the new node is labelled to keep following the branches
      if (use_states)
      {
        tracked <- track_graft(tracker = tracker, phy_before = phy_before_graft, phy_after = phy, node = MRCA_node_i, position = position_i,
                               tip_label = missing_taxa_list_i$Taxon_name, allowed_states = allowed_states_i, edge_ages = edge_ages_i)
        phy <- tracked$phy
        tracker <- tracked$tracker
      }
    } else {

      ## Case where the missing 'clade' is a true clade (nb_missing_taxa > 1)

      ## Grow a tree with the requested number of missing taxa under the same process as the sister-clade,
      # conditioned on BOTH the number of taxa and the stem age (= the divergence time drawn with TreeSim::corsim)
      missing_clade_i <- TreeSim::sim.bd.taxa.age(n = nb_missing_taxa, numbsim = 1,
                                                  lambda = BD_i$b, mu = BD_i$d,
                                                  frac = 1, age = divergence_time_missing_clade_i,
                                                  mrca = FALSE)[[1]] # mrca = FALSE: 'age' is the stem age (time since origin)
      missing_clade_i$tip.label <- missing_taxa_list_i$Taxon_name

      # Graft the generated subclade
      phy_before_graft <- phy
      phy <- ape::bind.tree(x = phy, y = missing_clade_i, where = MRCA_node_i, position = position_i)
      # With a simmap: the whole missing clade takes the state of the point where it was grafted, and the new nodes are labelled to keep following the branches
      if (use_states)
      {
        tracked <- track_graft_subtree(tracker = tracker, phy_before = phy_before_graft, phy_after = phy, node = MRCA_node_i, position = position_i,
                                       tip_labels = missing_taxa_list_i$Taxon_name, allowed_states = allowed_states_i, edge_ages = edge_ages_i)
        phy <- tracked$phy
        tracker <- tracked$tracker
      }
    }

    # Enforce ultrametricity for the new terminal branches only (backbone branches are left untouched)
    phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = missing_taxa_list_i$Taxon_name, root_age = root_age)

    # Update counts of grafted taxa and grafted missing clades
    nb_taxa_grafted <- nb_taxa_grafted + nb_missing_taxa
    nb_missing_clades_grafted <- nb_missing_clades_grafted + 1

    # Save temporary phylo if requested
    if (save_temp && (nb_missing_clades_grafted %% save_temp_every == 0))
    {
      if (is.null(temp_file_path))
      {
        temp_file_path <- "./temp_phy_from_grafting.rds"
      }
      # Save the phylogeny with original node labels restored
      saveRDS(object = restore_node_labels(phy = phy, labelled = labelled_backbone), file = file.path(temp_file_path))
    }

    if (verbose && (nb_missing_clades_grafted %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Missing taxa grafted within clades n\u00b0", nb_missing_clades_grafted, "/", n_groups,
                 " - Taxa n\u00b0", nb_taxa_grafted, "/", n_missing, "\n"))
    }
  }

  # Final check of ultrametricity (can take a long time for large phylogenies)
  is_ultrametric <- ape::is.ultrametric(phy, tol = ultrametric_tol)
  if (!is_ultrametric)
  {
    if (verbose) { cat(paste0(Sys.time(), " - Ensure the ultrametricity of the final tree\n") )}
    # Force ultrametry to avoid issue with bad rounding
    phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = missing_taxa_list$Taxon_name, root_age = root_age)
  }

  # Restore the original node labels (none if the backbone had none, and an empty string for the nodes created by grafting otherwise)
  phy_tracked <- phy # Save the tracked phylogeny to use to build the updated simmap below
  phy <- restore_node_labels(phy = phy_tracked, labelled = labelled_backbone)

  # With a simmap: the updated simmap of the grafted tree, built from the states followed during grafting (before the grafted elements are added)
  updated_simmap <- NULL
  if (use_states && return_simmap) { updated_simmap <- build_updated_simmap(tracker = tracker, phy_tracked = phy_tracked, phy = phy) }

  # Add the elements identifying the grafted tips and branches, and the sampling fraction table of the clades if 'total_richness_per_clades' was provided
  phy <- identify_grafted(phy = phy, grafted_taxa = missing_taxa_list$Taxon_name)
  if (!is.null(richness_info$table)) phy$sampling_fraction_df <- richness_info$table
  if (!is.null(updated_simmap)) phy$simmap <- updated_simmap

  return(phy)
}
