##### Functions to graft missing taxa #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Main function to graft missing taxa ####

#' @title Graft missing taxa on a time-calibrated phylogeny
#'
#' @description This function grafts missing taxa onto an ultrametric, time-calibrated backbone phylogeny.
#'   The position of each taxon is constrained by its taxonomic assignment (the **target clade**, defined in `Target_clade`),
#'   and optionally by taxonomic exclusions (the **excluded subclades**, defined in `Excluded_subclades`),
#'   and/or exclusions based on species attributes (the **allowed states**, defined in `allowed_states` based on the `simmap`).
#'
#'   The divergence times of the grafted branches are either drawn uniformly along the branches (`method = "uniform"`),
#'   or drawn from a birth-death model fitted on the target clade (`method = "BD_clades"`, default).
#'   **The backbone is never modified: only new tips are added.**
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which taxa are grafted.
#'   It must be rooted, binary (no polytomies) and ultrametric (i.e., be time-calibrated).
#'   It can be the output of a previous grafting: the grafted tips then accumulate (see 'Successive graftings').
#' @param missing_taxa_list A list with the information on the taxa to graft. **Every item provided must be named after `$Taxon_name`.** Items:
#'   * `$Taxon_name`: Character vector of the names of the taxa to graft. **Names must be unique and absent from `phy$tip.label`.**
#'   * `$Target_clade`: List (named after `$Taxon_name`) of character vectors. For each taxon, the tip labels defining its **target clade**,
#'     i.e., the clade where the taxon is known to belong (e.g., its genus). A pair of tips defines the clade as their MRCA node.
#'     A single tip defines the "clade" as the terminal branch of this tip.
#'     **Clades, including the excluded subclades (see below), must be defined by tips: node IDs are not accepted** (the function stops, and lists some of them).
#'     Use [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] to convert node IDs into pairs of tips.
#'   * `$Target_clade_name`: (Optional) Character vector (named after `$Taxon_name`) of the names of the target clades (`NA` if a clade has no name).
#'     **Taxa sharing the same target clade must have the same name.** The names are not used for the grafting, but, if `total_richness_per_clades` is provided,
#'     they are checked against its `$Target_clade_name` (a clade must have the same name in both objects).
#'   * `$include_target_clade_stem`: (Optional) Logical vector (named after `$Taxon_name`). For each taxon, whether the stem branch of its target clade
#'     is available for grafting. It overrides the global argument `include_target_clade_stem` for this taxon. `NA` values use the global argument.
#'     With `NULL` (default, or no item), the global argument applies to all taxa. Ignored for single-tip target clades (their terminal branch is their only branch,
#'     so the stem is always used; a warning lists the taxa with `FALSE`) and for target clades that are the whole tree.
#'   * `$Excluded_subclades`: (Optional) List (named after `$Taxon_name`) of lists of character vectors. For each taxon, the **excluded subclades** of its target clade
#'     where the taxon is known not to belong. Each subclade is defined by a pair of tips (or a single tip, which excludes its terminal branch, see `exclude_subclades_stem`).
#'     Use `NULL` for taxa with no excluded subclades.
#'   * `$exclude_subclades_stem`: (Optional) List (named after `$Taxon_name`) of logical vectors. For each taxon, whether the stem of each of its excluded subclades
#'     is also excluded, with one value per subclade of `$Excluded_subclades` (a single value is used for all the subclades of the taxon; `NULL` uses `exclude_subclades_stem`).
#'     It overrides the global argument `exclude_subclades_stem` for these subclades. `NA` values use `exclude_subclades_stem`.
#'     A logical vector with one value per taxon is also accepted (each value is used for all the excluded subclades of the taxon).
#'     For example, `Excluded_subclades = list(X1 = list(c("t1", "t2"), "t7"), X2 = NULL)` with `exclude_subclades_stem = list(X1 = c(FALSE, TRUE), X2 = NULL)`.
#'     The excluded subclades of a taxon can also be named (e.g., `list(X1 = list(grp_A = c("t1", "t2"), grp_B = "t7"))`): these optional names are used in messages,
#'     and must match the names of the values of `exclude_subclades_stem` of the taxon if both are named.
#'   * `$allowed_states`: (Optional) List (named after `$Taxon_name`) of character vectors. For each taxon, the states of the `simmap` in which it can be grafted.
#'     It **overrides the global argument** `allowed_states` for this taxon. `NULL` or `NA` use the global argument.
#'     Only used with a `simmap`. For example, `allowed_states = list(X1 = c("blue", "green"), X2 = NULL)`.
#'     With the maps of several traits, each element is a list with the states of each trait (`NULL`: no taxon-level value for the trait), e.g.,
#'     `list(X1 = list(c("blue", "green"), "high"), X2 = list(NULL, "low"))` (see 'Several traits' in Details section).
#' @param method Character string. Method used to draw the divergence times and the branches used for grafting.
#'   Must be one of **`"uniform"`** or **`"BD_clades"`** (default). See Details.
#' @param min_tips Integer. Only used for `method = "BD_clades"`. Minimum number of tips of the clade used to fit the birth-death model.
#'   If the target clade has fewer tips, the smallest enclosing clade with at least `min_tips` tips is used instead. **Must be at least 3.** Default is `5`.
#' @param return_BD_model_summary Logical. Only used for `method = "BD_clades"`. Whether to report the birth-death models fitted during the grafting.
#'   If `TRUE`, each grafted phylogeny has the attribute `"BD_model_summary"`: a data.frame with one row per model (i.e., per target clade),
#'   summarizing diversification rates and sampling fraction.
#'   It can be extracted a posteriori with [`phylograftR::extract_BD_model_summary()`][extract_BD_model_summary]. Default is `FALSE`.
#' @param include_target_clade_stem Logical. Whether the stem branch of target clades is available for grafting. This is the default for all taxa:
#'   taxa with a value in `missing_taxa_list$include_target_clade_stem` use this specific instead.
#'   If `TRUE`, **a taxon can be grafted on the stem branch, and becomes sister to the target clade**.
#'   Taxa already grafted on the stem can receive further grafts (see 'stem region' in Details section).
#'   Ignored for single-tip target clades as their terminal branch is their only branch, and is therefore always available for grafting.
#'   Ignored if the target clade is the whole tree: it has no stem. A warning is issued if the stem is requested for such a taxon.
#'   Default is `TRUE`.
#' @param use_excluded_subclades Logical. Whether to avoid grafting within the excluded subclades provided in `missing_taxa_list$Excluded_subclades`. Default is `TRUE`.
#' @param exclude_subclades_stem Logical. Whether the stem of the excluded subclades is also excluded. This is the default for all the excluded subclades:
#'   those with a value in `missing_taxa_list$exclude_subclades_stem` use this specific value instead.
#'   If `TRUE`, **a taxon cannot be grafted on the stem of an excluded subclade** (it cannot be sister to it): the subclade is excluded entirely, crown and stem.
#'   If `FALSE`, **only the branches within the subclade are excluded**, and the taxon can be sister to it.
#'   For a single-tip excluded subclade, the terminal branch is the only branch of the tip, and it is its stem: the tip excludes nothing if `FALSE`
#'   (a warning lists these cases). Only used if `use_excluded_subclades = TRUE`. Default is `TRUE`.
#' @param total_richness_per_clades (Optional) List with the total number of species of the clades, to account for species that are neither in the backbone
#'   nor in `missing_taxa_list` (unsampled species). Items:
#'   * `$Target_clade`: List of character vectors of tips, defining clades as in `missing_taxa_list$Target_clade`.
#'     **Node IDs are not accepted**: convert them with [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips] in the same call as `missing_taxa_list`,
#'     so the tips of the clades match between both objects.
#'     The easiest way to obtain this list is to build it from a table of clades with [`phylograftR::build_total_richness_per_clades()`][build_total_richness_per_clades]
#'     or [`phylograftR::import_total_richness_per_clades()`][import_total_richness_per_clades]: the clades can be written with names, tips or node IDs,
#'     they are matched with the clades of `missing_taxa_list`, and the richness is checked against all the taxa to graft.
#'     **It must cover all the unique target clades of `missing_taxa_list`.** It can include other clades, e.g., enclosing clades.
#'   * `$Target_clade_name`: (Optional) Character vector of the names of the clades, same length as `$Target_clade` (`NA` if a clade has no name).
#'     The names must be unique. If `missing_taxa_list$Target_clade_name` is provided, **a clade must have the same name in both objects** (the function stops otherwise).
#'     They are not used for the grafting, but they are reported by [`phylograftR::extract_sampling_info()`][extract_sampling_info].
#'   * `$total_richness`: Vector of positive integers, same length as `$Target_clade`. Total number of species of each clade,
#'     i.e., backbone species, missing species nested within the clade, and unsampled species. Use `NA` if unknown: the unsampled species
#'     are then assumed to be those of the clades with an informed total richness nested within the clade (none if there are no such clades).
#'
#'   With `NULL` (default), no unsampled species are accounted for.
#'
#'   If the list is provided with only `NA` in `$total_richness`, no unsampled species are assumed in any clade
#'   (`total_richness` = backbone + missing species): the fits are identical to those with `NULL`.
#'
#'   To obtain the resulting richness and the sampling fraction of the clades in each grafted phylogeny,
#'   use [`phylograftR::extract_sampling_info()`][extract_sampling_info] on the output, with the same `total_richness_per_clades`.
#'
#'   See 'Sampling fraction and unsampled species' in Details section.
#' @param simmap (Optional) A stochastic character map of the backbone, of class `"simmap"` (e.g., from `phytools::make.simmap()` or `phytools::sim.history()`),
#'   or a set of maps of the same trait, of class `"multiSimmap"` (e.g., from `phytools::make.simmap(nsim = 10)`), used as alternative hypotheses of the evolution of a species attribute
#'   (e.g., a trait or a range) with discrete states along the branches of `phy`.
#'
#'   Alternatively, you can provide a list of `"simmap"` and/or `"multiSimmap"` objects, one per trait, that contains the maps of several traits used together as constraints.
#'   **All maps must have the same tips, topology and branch lengths as `phy`.**
#'   With a `simmap`, **missing taxa can only be grafted on the parts of the branches that are in their allowed states** (see `allowed_states`).
#'   Each grafted taxa inherits the state(s) of the point of the branch piece where it is grafted.
#'
#'   See 'Species attributes', 'Multiple stochastic maps' and 'Several traits' in Details section. Default is `NULL` (no restriction for grafting based on states).
#' @param allowed_states Character vector. The states of `simmap` in which grafting is allowed. This is the default for all taxa:
#'   taxa with a value in `missing_taxa_list$allowed_states` use this specific value instead. States are matched exactly with the states of the `simmap`, unless `treat_states_as_ranges = TRUE`.
#'   **Taxa with no allowed states (neither a global nor a taxon-level value) can be grafted in all the states of the `simmap`.**
#'   The function stops if some states are not in the `simmap`, or if `allowed_states` is provided without a `simmap`.
#'   With `treat_states_as_ranges = TRUE`, a range that is not a state of the `simmap` is accepted if it is contained in some of them (e.g., `"A"` can be allowed if `"AB"` is found on the `simmap`).
#'
#'   With the maps of several traits (a list in `simmap`), you must provide a list with one element per trait, in the order of the traits or named after them.
#'   In each element, the states allowed for each trait are named, with `NULL` for no restriction on a trait, e.g., `list(c("blue", "green"), NULL, "high")`
#'   (see 'Several traits' in Details section).
#'   Default is `NULL` (all states allowed).
#' @param treat_states_as_ranges Logical. Only used with a `simmap`. Whether the states of the `simmap` are biogeographic ranges to be treated as sets of areas, rather than as unrelated states.
#'   **Ranges must be coded as in `BioGeoBEARS`**: one capital letter per area, and the letters of a multi-area range in alphabetical order (e.g., `"A"`, `"AB"`, `"ABC"`).
#'   If `TRUE`, an allowed range enables grafting on branch pieces in any range that contains it: `"A"` can be grafted on `"A"`, `"AB"` and `"ABC"`,
#'   whereas `"AB"` can be grafted on `"AB"` and `"ABC"` but not `"A"`.
#'
#'   If stochastic maps of several traits are provided (a list in `simmap`), you must provide a logical vector with one value per trait,
#'   with at most one `TRUE` to define the trait treated as biogeographic ranges, the other traits being matched exactly.
#'   The vector is positional (e.g., `c(FALSE, TRUE, FALSE)`: the second trait is the range trait) or named after the traits, in any order (e.g., `c(biogeo = TRUE, color = FALSE, elevation = FALSE)`).
#'   A single `FALSE` is accepted with several traits (no range trait), but a single `TRUE` stops.
#'
#'   See 'Biogeographic ranges as states' and 'Several traits' in Details section. Default is `FALSE` (states are matched exactly).
#' @param return_simmap Logical. Only used with a `simmap`, and silently ignored without it. Whether to add the updated `simmap` to the output, in `$simmap` (see 'Value').
#'   Default is `TRUE`.
#' @param min_valid_fraction Numeric between 0 and 1. Only used with a `simmap`. Threshold below which a warning lists the taxa for inspection,
#'   as their grafting positions are concentrated on a small part of their target clade. The fraction is computed on the backbone, before grafting:
#'   * With `method = "uniform"`, it is the fraction of the total length of the available branches of the taxon (its target clade, minus the excluded subclades)
#'     that is in its allowed states.
#'   * With `method = "BD_clades"`, it is the fraction of the time window of the divergence times of the taxon (from the present to the stem age of its target clade,
#'     or to its crown age if the stem is not available) where at least one available branch is in an allowed state.
#'   Default is `0.1`.
#' @param max_age_redraws Integer. Only used with `method = "BD_clades"`. Maximum number of times the divergence ages of a group of taxa
#'   (same target clade, same stem availability and, with a `simmap`, same allowed states) can be redrawn, when some of the ages drawn cannot be used because no available branch
#'   (in an allowed state with a `simmap`) crosses them. The function stops with an informative message if this limit is reached. See 'Methods' and 'Species attributes' in Details section. Default is `1000`.
#' @param seed Integer. Seed used for reproducibility. The seed only applies within the function: the state of the random number generator of the R session
#'   is restored on exit. Set to `NULL` to use the random number generator of the session. Default is `1234`.
#'   With several simulations, it is the master seed: the first simulation uses `seed`, and the seeds of the other simulations are derived from it (see 'Multiple simulations').
#' @param nb_simulations Integer. Number of simulations, i.e., of grafted phylogenies to produce in the same run, from the same backbone and the same `missing_taxa_list`.
#'   The grafting is stochastic, so **ideally downstream analyses should be run on a set of grafted phylogenies rather than on a single one**. Default is `1`.
#'   See the section 'Multiple simulations' for the output, the seeds and the failures.
#' @param force_simulations Logical. Whether to try new seeds until `nb_simulations` grafted phylogenies are obtained, up to `10 * nb_simulations` attempts
#'   (see 'Multiple simulations' in Details section). If `FALSE`, the simulations that fail are skipped. Default is `FALSE`.
#' @param parallel Logical or `NULL`. Whether to run the simulations in parallel, with the packages `future` and `future.apply` (**they must be installed separately from phylograftR**).
#'   Default is `NULL`: the simulations are run in parallel if the user set a `future` plan with several workers beforehand (e.g., `future::plan(future::multisession, workers = 4)`), and sequentially otherwise.
#'   `FALSE` forces a sequential run, even if a plan is set. `TRUE` execute the `future` plan, but stops with an informative error if the packages are not installed or no plan is set.
#'   See 'Parallel computation' in Details section for recommendations on whether parallelization is useful in your case.
#' @param save_temp Logical. Whether to save the grafted phylogenies during the process, to recover from a crash or an interruption of long runs (see 'Temporary saves and resuming' in the Details section).
#'   Each phylogeny (one per simulation and per stochastic map) has its own temporary file. Its progress is saved every `save_temp_every` taxa or groups of taxa (missing clades in `graft_missing_clades()`), if not `NULL`, and the phylogeny is always saved once fully grafted. Default is `FALSE`.
#' @param save_temp_every Integer or `NULL`. Frequency of the saves of the progress, in number of taxa for `method = "uniform"`, and groups of taxa for `method = "BD_clades"`
#'   (the taxa grafted together because they share the same target clade, availability of its stem and allowed states). Only used if `save_temp = TRUE`.
#'   Default is `NULL`: save each phylogeny only once fully grafted.
#' @param temp_file_path Character string. Path of the `.rds` file used for the temporary saves. It must end with `.rds` and its directory must exist and be writable.
#'   Default is `NULL`, which saves in `"temp_phy_from_grafting.rds"` in the working directory.
#'   With several simulations and/or several stochastic maps, the path is a template:
#'   the suffixes `"_map_<name>"` (several maps; `"map_"` is not repeated if the name of the map already starts with it, e.g., `"temp_map_1.rds"` for a map named `"map_1"`) and `"_sim_<k>"` (several simulations) are added before `".rds"`, in this order (e.g., `"temp_map_A_sim_2.rds"`).
#'   See 'Temporary saves and resuming' in the Details section.
#' @param clean_temp_files Logical. Whether to delete all the temporary files of the run once it is completed without stopping/crashing. Only used if `save_temp = TRUE`.
#'   If the run crashes or is interrupted, nothing is deleted so it can be resumed. Use `FALSE` to keep the files. Default is `TRUE`.
#' @param resume Logical. Whether to resume an interrupted run from its temporary files, which **must come from a run with the exact same arguments** (see 'Temporary saves and resuming').
#'   Needs `save_temp = TRUE`. Default is `FALSE`.
#' @param ultrametric_tol Numeric. Tolerance used by `ape::is.ultrametric()` for the ultrametricity checks of the backbone and of the final phylogeny.
#'   Default is `.Machine$double.eps^0.5`.
#' @param verbose Logical. Whether to display progress. Messages print every 100 taxa for `method = "uniform"`, and every 10 target clades for `method = "BD_clades"`.
#'   Default is `TRUE`.
#'
#' @return The function returns an object of class `"phylo"`: the backbone phylogeny with the grafted taxa added as new tips.
#'   Backbone branch lengths are untouched, and the grafted phylogeny is ultrametric.
#'
#'   In addition to the standard elements of a `"phylo"` object, the phylogeny has the following elements:
#'   * `grafted_tips` Named logical vector with one value per tip, in the order of `$tip.label`. `TRUE` for grafted tips (including the tips grafted by a previous grafting if `phy` was already grafted). Names are the tip labels.
#'   * `grafted_edges` Named logical vector with one value per branch, in the order of the rows of `$edge`. `TRUE` for grafted branches,
#'     i.e., branches with only grafted taxa below them. Names are the IDs of the tipward nodes of the branches (`$edge[, 2]`).
#'     Use them to color the phylogeny on plots (see examples).
#'
#'   **These elements are not updated by functions that modify the phylogeny** (e.g., `ape::drop.tip()`, which removes tips but keeps the other elements as they are,
#'   and `ape::write.tree()`, which removes them). Use [`phylograftR::identify_grafted()`][identify_grafted] to recompute `grafted_tips` and `grafted_edges` on a manipulated phylogeny.
#'
#'   Optional elements can be added to the output based on what optional argument was provided:
#'
#'   * `"BD_model_summary"` (attribute). Only if `return_BD_model_summary = TRUE`. Data.frame with one row per birth-death model fitted (one per target clade),
#'     summarizing information about the clade used for model fitting, diversification rates, and sampling fraction.
#'     It can be extracted a posteriori with [`phylograftR::extract_BD_model_summary()`][extract_BD_model_summary].
#'   * `simmap` Only if a `simmap` is provided and `return_simmap = TRUE`. The updated stochastic character map of the grafted phylogeny (class `"simmap"`):
#'     the pieces of the branches of the backbone keep their states (cut where taxa were grafted), and the grafted branch is entirely in the state of the point where it was grafted.
#'     With `treat_states_as_ranges = TRUE`, the grafted branch is entirely in the range it inherited (see 'Biogeographic ranges as states' in Details section).
#'     This element can be plotted as any `simmap`, e.g., `plot(grafted_phy$simmap)`, or used to extract the states of the tips with `phytools::getStates(grafted_phy$simmap, type = "tips")`.
#'     It has the same exact topology and branch lengths as the grafted phylogeny.
#'     With the maps of several traits, `$simmap` is a list of such simmaps, one per trait, named after the traits (see 'Several traits' in Details section).
#'
#'   With `nb_simulations > 1`, the function returns an object of class `"multiPhylo"` with one grafted phylogeny per simulation, named `"sim_1"`, `"sim_2"`, ...
#'   Each phylogeny has the elements described above (including `$simmap` if a `simmap` was provided and `return_simmap = TRUE`).
#'   The `"multiPhylo"` object has the attributes `"seeds"` (named vector of the seeds of the phylogenies, to reproduce each of them alone with `nb_simulations = 1`;
#'   `NULL` if `seed = NULL`) and `"failed"` (data.frame with the seed and the reason of the simulations that failed, 0 row if none).
#'
#'   With multiple stochastic maps (i.e., a `"multiSimmap"`) provided as `simmap`, the function returns a list with one element per map (named after the maps):
#'   a phylogeny if `nb_simulations = 1`, a `"multiPhylo"` otherwise (see 'Multiple stochastic maps').
#'   The list has the attribute `"failed_maps"` (data.frame with the `map` and the `reason` of the maps without any simulation, dropped from the output).
#'   The phylogenies have the attributes `"seed"` (their seed, unless `seed = NULL`) and, with several maps, `"map"`, to ensure reproducibility of each simulation.
#'
#' @details
#' ### Stop conditions
#'
#' The function **stops, before any grafting,** if the backbone is not ultrametric, or is not binary (i.e., all nodes have only two descendant lineages).
#' It also stops if some taxa have no branch available for grafting. The taxa concerned are listed in the error message,
#' so the constraints (`include_target_clade_stem`, `Excluded_subclades`, `exclude_subclades_stem` or, with a `simmap`, the allowed states) can be adjusted.
#' With `method = "BD_clades"`, the function can also stop during grafting if no usable divergence age is found after `max_age_redraws` redraws (see 'Methods').
#'
#' ### Target clades, stem regions and excluded subclades
#'
#' Each taxon is grafted within its target clade. The branches available for grafting are all the branches of the target clade
#' (descending from its MRCA node), plus its stem region if the stem is available for this taxon (`missing_taxa_list$include_target_clade_stem` of the taxon if provided, otherwise the global argument `include_target_clade_stem`).
#' **Grafting on the stem branch places the taxon sister to the target clade.**
#' Taxa of the same target clade can therefore differ in their stem availability: the clade is the same (same MRCA node, same richness information),
#' but the branches and the time window available for the divergence times of these taxa differ.
#'
#' **The stem region is not restricted to the original stem branch.** A taxon grafted on the stem splits it in two pieces, and becomes itself part of the stem region:
#' later taxa can be grafted on both pieces of the stem, and on the terminal branch of the earlier taxon (the two taxa are then sister, and together sister to the target clade).
#' More generally, the stem region includes all the branches of the 'missing-only subtrees' hanging from the stem, i.e., subtrees containing only taxa of `missing_taxa_list$Taxon_name`,
#' whatever their own target clade. This mirrors what happens within the crown group, where missing taxa can already be grafted on the branches of other missing taxa.
#'
#' **Branches inside excluded subclades are removed from the pool of available branches.** By default (`exclude_subclades_stem = TRUE`),
#' the stem region of the excluded subclades is removed as well: the taxon can be neither nested within the excluded subclade,
#' nor sister to it (nor sister to the missing taxa already grafted on its stem).
#' This is the only way to exclude a single tip, as its terminal branch is its only branch.
#' With `exclude_subclades_stem = FALSE` (globally, or for a given subclade in `missing_taxa_list$exclude_subclades_stem`), the stem of the excluded subclade stays available:
#' the taxon can be sister to the excluded subclade, but never nested within it. A single tip then excludes nothing.
#' The setting is specific to each excluded subclade, so a taxon can have excluded subclades with different settings.
#'
#' **Excluded subclades are expected to be nested within their target clade**. The function checks their position before grafting:
#' a warning is issued if an excluded subclade is outside the target clade (the exclusion then has no effect).
#' The function stops if an excluded subclade contains the whole target clade, or is identical to it: no branch (or only the stem branch) would be left for grafting.
#' To graft taxa specifically on the stem branch of a clade or of a single tip, use [`phylograftR::graft_missing_clades()`][graft_missing_clades] instead.
#' A warning lists the single-tip excluded subclades with a stem that is not excluded, as they exclude nothing.
#'
#' ### Methods
#'
#' * **`"uniform"`**: Taxa are grafted one by one, per target clade from the smallest to the largest, so nested clades are processed first (see below).
#'   Within each target clade, taxa are grafted in random order. The branch used for grafting is drawn among the available branches
#'   with probabilities proportional to branch lengths, and the position on the branch is drawn from a uniform distribution.
#'   This is fast, but does not account for diversification dynamics.
#'   Note that under a birth-death process, the probability of a side lineage splitting off a branch at age t follows the exponential function e^(-rt),
#'   with r the net diversification rate (speciation minus extinction).
#'   A uniform draw implicitly assumes a flat density (r = 0), and gives too much weight to old divergence times compared to the expectation from a birth-death process.
#'
#' * **`"BD_clades"`**: Taxa are grafted per target clade, from the smallest to the largest, so nested clades are processed first (see below).
#'   For each target clade, a birth-death model is fitted with `phytools::fit.bd()` on the clade, or on the smallest enclosing clade with at least `min_tips` tips.
#'   The sampling fraction accounts for the number of missing taxa provided for grafting in this clade.
#'   The divergence times of all taxa missing in the target clade are then jointly drawn with `TreeSim::corsim()`, between the present and the stem age of the clade
#'   (or its crown age if the stem is not available for these taxa).
#'   Taxa are grafted from the oldest to the youngest divergence time following this procedure:
#'     1. A taxon is grafted at its age if an available branch in an allowed state crosses this age in the current tree.
#'     2. The ages for which no such valid branch piece exists are discarded.
#'     3. The ages of the taxa left are drawn again jointly, and grafted or discarded in the same way.
#'     This is repeated until all the taxa of the group are grafted, up to `max_age_redraws` times.
#'   Taxa of the same target clade with different stem availability (or, with a `simmap`, different allowed states) are drawn in separate groups,
#'   and the birth-death model is fitted only once per target clade.
#'   **`"BD_clades"` is the recommended grafting method**, but it can be time consuming for large phylogenies.
#'
#' With both methods, **nested clades are processed first**, so the taxa of a larger clade are grafted on a tree that already
#' includes the taxa of the clades nested within it, and can be grafted on their branches.
#' For taxa of the same clade, those without access to the stem come first.
#' For `"BD_clades"`, it also ensures that the divergence times of nested clades are accounted for when drawing the divergence times of the larger clade.
#' Taxa of the same target clade are grafted in random order, so they are exchangeable.
#'
#' ### Enclosing clades
#'
#' Under the `"BD_clades"` method, when a target clade contains fewer tips than `min_tips`, an alternative higher clade is selected to fit the birth-death model.
#' This 'enclosing clade' is identified by walking rootward the tree node by node, starting from the target clade, until a clade meets the `min_tips` condition.
#' This ensures that the diversification dynamics modeled are based on enough tips to be relevant, while being closely associated with the target clade.
#'
#' ### Sampling fraction and unsampled species
#'
#' Under `method = "BD_clades"`, the birth-death model of a clade is fitted with the sampling fraction
#' `n_present / (n_present + n_pending + n_unsampled)`, with `n_present` the tips of the fitted clade in the current tree (backbone tips and taxa already grafted),
#' `n_pending` the taxa of `missing_taxa_list` still to graft in this clade, and `n_unsampled` the unsampled species derived from `total_richness_per_clades`.
#' Taxa grafted in nested clades processed earlier are already in the tree, so they are not counted as missing anymore.
#' **Unsampled species are only used to fit the rates: they are not grafted**, and do not compete with the grafted taxa for the branches.
#'
#' For a clade with a `total_richness` value provided: `nb_unsampled = total_richness - nb_backbone - nb_missing`, with `nb_missing` the number of taxa to graft
#' whose target clade is this clade or is nested within it. For a clade with no `total_richness` value (`NA`, or enclosing a target clade without being listed),
#' `nb_unsampled` is the sum over the largest clades nested within it with a `total_richness` value provided.
#' As such, the unsampled species of nested clades are summed across disjoint clades but never counted twice within nested clades.
#'
#' If `total_richness_per_clades = NULL`, `nb_unsampled = 0`: the clades have no unsampled species (the same as providing only `NA` in `$total_richness`).
#'
#' The species actually present in each listed clade, and the resulting sampling fractions for downstream analyses, may differ among grafted phylogenies,
#' as a taxon can be grafted on the stem of its target clade, or within a clade nested in it, which changes the richness of these clades.
#' Use [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract them a posteriori.
#'
#' The function stops, before any grafting, if `total_richness_per_clades` does not cover all the target clades, if a `total_richness` provided is lower than
#' the number of backbone species plus species to graft of the clade, or if a clade has fewer unsampled species than the sum of the clades with
#' a `total_richness` provided nested within it.
#' With `method = "uniform"`, no diversification model is used: `total_richness_per_clades` is not used during the grafting (a warning is issued).
#' It is only checked, and can be provided later to [`phylograftR::extract_sampling_info()`][extract_sampling_info].
#'
#' ### Species attributes: allowed states of a simmap
#'
#' A `simmap` describes the states of a discrete attribute (e.g., a trait or a range) along the branches of the backbone: each branch is cut in pieces, each piece being in one state.
#' With a `simmap`, **each taxon can only be grafted on the pieces of the available branches that are in its allowed states**
#' (`missing_taxa_list$allowed_states` of the taxon if provided, otherwise the global argument `allowed_states`).
#' This combines with the other constraints: the available branches are those of the target clade, minus the excluded subclades, restricted to the allowed states.
#' If all states are allowed (all states listed in `allowed_states`, or no allowed states at all), grafting is not restricted by the states,
#' and constraints are only based on target clades and excluded subclades.
#'
#' **A grafted taxon is assumed to be in the state of the point of the backbone where it was grafted, along its whole branch**,
#' and so are the taxa grafted on its branch later. Taxa grafted later can therefore be grafted on the branches of taxa grafted earlier,
#' if their state is allowed.
#'
#' The states are updated alongside the phylogeny, and the updated `simmap` is returned in `$simmap` if `return_simmap = TRUE`.
#' A `simmap` can thus be provided with `allowed_states = NULL`, only to get the states of the grafted taxa and branches:
#' grafting is then not restricted by the states and gives the same result as without `simmap` (for the same `seed`).
#'
#' With `method = "uniform"`, the draw is uniform over the available pieces (i.e., portion of branches in the allowed states), weighted by their length:
#' every time point in an allowed state has the same chance to be selected.
#'
#' With `method = "BD_clades"`, the divergence ages of the taxa of a group (same target clade, same stem availability and same allowed states) are drawn jointly
#' with `TreeSim::corsim()`, as without a `simmap`: the diversification model is not modified by the states.
#' The same grafting procedure, described above without a `simmap`, is then used except that
#' an age can only be used if an available branch in an allowed state crosses it (without a `simmap`, if an available branch crosses it).
#' The ages are used from the oldest to the youngest, each for one taxon of the group.
#'
#' **Cautions for interpretation: risks of circularity.**
#'
#' Constraining the grafting by a trait that will then be analyzed (e.g., the association of the trait with diversification) may bias the results towards the constraint.
#' The grafted taxa are placed only where the trait allows it, so they are in the allowed states by construction. This can create or inflate an association between the trait and the diversification dynamics,
#' affect the rates of transitions between states, or inflate the phylogenetic signal of the trait. It also reduces the uncertainty about the state of the grafted lineages since they always inherit a unique state.
#' It is recommended to acknowledge these potential biases when interpreting such analyses.
#'
#' To **limit the risk of circularity**, multiple approaches are possible:
#' * constrain the grafting by traits that are not the object of downstream analysis;
#' * compare the results with and without the constraint (or on the backbone only);
#' * keep the constraint as loose as the information allows (e.g., allowed states that include neighboring states);
#' * use several maps and simulations to propagate the uncertainty in ancestral state estimates.
#'
#' ### Biogeographic ranges as states
#'
#' With `treat_states_as_ranges = TRUE`, the states of the `simmap` are biogeographic ranges, coded as in `BioGeoBEARS`
#' (one capital letter per area, and the letters of a multi-area range in alphabetical order, e.g., `"AB"`),
#' and the allowed states are the ranges a missing clade can have.
#'
#' **An allowed range matches a state of the `simmap` if it is contained in it**: a missing taxon allowed in `"A"` can be grafted on branch pieces in `"A"`, `"AB"` or `"ABC"`,
#' whereas a missing taxon allowed in `"AB"` only can be grafted on branch pieces in `"AB"` or `"ABC"`.
#'
#' When multiple allowed states are provided, the inherited range is drawn randomly from the allowed ranges enclosed within the range of the branch piece where it is grafted.
#' When grafting on an edge with a larger range than the inherited range (e.g., `"A"` on `"AB"`), this results in a subset sympatric speciation event such as (AB) -> (AB), (A).
#' In the returned `simmap`, the branches of the grafted lineage are in their inherited range, whereas the branch it is grafted on continues in its own range.
#'
#' ### Multiple simulations
#'
#' With `nb_simulations > 1`, the grafting is simulated `nb_simulations` times from the same backbone and the same `missing_taxa_list`.
#' Each simulation starts again from the backbone: the birth-death models are fitted on the tree of the simulation (they depend on the taxa already grafted),
#' and the simulations do not share anything but the inputs. The phylogenies are returned in a `"multiPhylo"` object (a single `"phylo"` if `nb_simulations = 1`),
#' and the memory used grows with the number of simulations.
#'
#' **The first simulation uses `seed`**, so it is the same phylogeny as a run with `nb_simulations = 1` and the same `seed`.
#' The seeds of the other simulations are drawn from a stream started by `seed`, and returned in the attribute `"seeds"`:
#' a simulation can be reproduced alone with its seed, and the first simulations of a longer run are the same as the ones of a shorter run.
#' With `seed = NULL`, the simulations follow the random number generator of the session.
#'
#' A failure that is specific to a simulation, because of the random draws (no usable divergence age or time after `max_age_redraws` redraws), skips this simulation with a warning:
#' the name of a phylogeny keeps the index of its simulation (e.g., `"sim_3"` is missing if the third simulation failed). The function reports the number of phylogenies obtained
#' if it differs from `nb_simulations`, and stops if no simulation succeeded. With `force_simulations = TRUE`, new seeds are tried until `nb_simulations` phylogenies are obtained
#' (up to `10 * nb_simulations` attempts): the phylogenies are then named `"sim_1"` to `"sim_N"` without gaps, and the seeds of the failed attempts are in the attribute `"failed"`.
#' Any other error (invalid inputs, no branch available for a taxon, ...) concerns all the simulations and stops the whole run.
#'
#' With `save_temp = TRUE`, each simulation saves its own temporary file (see 'Temporary saves and resuming').
#'
#' ### Multiple stochastic maps
#'
#' When multiple stochastic maps of the same trait are provided as a `"multiSimmap"`, the grafting is run with each of its maps.
#' The maps are alternative hypotheses of the evolution of the same trait (or range), used as alternative constraints.
#' **All the maps and the simulations are crossed**: each map is used for `nb_simulations` simulations, so `P` maps give `P * nb_simulations` grafted phylogenies
#' This can lead to a quick inflation of computation time, thus a warning is sent if the total number of grafted phylogenies (`P * nb_simulations`) is more than **1000**.
#'
#' The output is a list with one element per stochastic map, named after the maps (the names of the `"multiSimmap"` if it has some, otherwise `"map_1"`, ..., `"map_P"`).
#' With `nb_simulations = 1`, each element is a grafted phylogeny. With several simulations, each element is a `"multiPhylo"` (see 'Multiple simulations')
#' with phylogenies named `"<map>_sim_1"`, `"<map>_sim_2"`, ... and the attributes `"seeds"` and `"failed"` of the map.
#' Each phylogeny can be reproduced alone with its map (as a single `simmap`) and its seed.
#'
#' A map for which no simulation could be obtained is dropped from the output with a warning, and listed in the attribute `"failed_maps"` of the list
#' (data.frame with the columns `map` and `reason`). The function stops if no map gives any phylogeny.
#'
#' With `save_temp = TRUE`, each map and each simulation saves its own temporary file (see 'Temporary saves and resuming').
#'
#' ### Several traits
#'
#' The maps of several traits can be used at once to constrain grafting. `simmap` is then a list with one element per trait, each a `"simmap"` or a `"multiSimmap"`
#' (e.g., `simmap = list(color = color_simmap, elevation = elevation_simmap)`). The elements can be named after the traits (all or none, with different names that cannot contain `;`, `,`, `=` or `|`),
#' otherwise the traits are named `"trait_1"`, `"trait_2"`, ... All the maps must have the same tips, topology and branch lengths as the backbone `phy`.
#'
#' **A taxon can only be grafted on the parts of the branches that are in one of its allowed states, for every trait.**
#' The allowed states are provided for each trait, in the order of the list of maps (or named after the traits):
#' the global argument `allowed_states` is a list with one element per trait, e.g., `list(c("blue", "green"), NULL, "high")` for three traits,
#' and `missing_taxa_list$allowed_states` has such a list for each taxon.
#' A trait with `NULL` (or an empty vector, or `NA`) has no taxon-level value, so the global `allowed_states` of the trait is used for the taxon
#' (a taxon-level value replaces the global one, trait by trait, as with a single trait).
#' A trait with no states, neither at the taxon level nor globally, has all its states allowed (no restriction).
#' Other traits can still be restricted: `list(NULL, "high")` only restricts the second trait.
#'
#' With `return_simmap = TRUE`, `$simmap` is a list of updated simmaps, one per trait, named after the traits.
#' In each simmap, a grafted branch is entirely in the state of each trait at the point where it was grafted (each trait is inherited from the lineage it was grafted on).
#'
#' With a list of `"multiSimmap"` provided as `simmap`, the maps of the traits are paired: the grafting is run with the map 1 of every trait, then the map 2 of every trait, and so on.
#' **All the traits must therefore have the same number of maps** to be associated in a given grafting process.
#'
#' One of the traits can be treated as biogeographic ranges combined with any number of other traits.
#' Use a logical vector in `treat_states_as_ranges` to define which of the trait is to be treated as biogeographic ranges for inheritance.
#' The range logic is the one of 'Biogeographic ranges as states' for the range trait:
#' its allowed values for grafting are any ranges that encompass the allowed range.
#' For instance, a missing taxon allowed in `"A;high"` can be grafted on branch pieces in `"A;high"`, `"AB;high"` or `"ABC;high"`,
#' whereas a missing taxon allowed in `"AB;high"` only can be grafted on branch pieces in `"AB;high"` or `"ABC;high"`.
#'
#' The other traits are matched exactly, and are inherited as they are from the point where the taxon is grafted: with `"AB"` and `"high"` at the graft point and the allowed range `"A"`,
#' the taxon inherits `"A"` and `"high"`. **Only one trait can be treated as ranges** to prevent inflation of the number of combined allowed states to consider.
#'
#' ### Parallel computation
#'
#' The simulations based on the different stochastic maps and seeds are independent, so they can be run in parallel with the package `future` and `future.apply`.
#' **Both packages must be installed independently from phylograftR**.
#'
#' The grafting function does not manage the workers: **you must set a plan beforehand**, e.g., `future::plan(future::multisession, workers = 4)`,
#' and reset it afterwards with `future::plan(future::sequential)`. With the default `parallel = NULL`, the function uses the plan if there is one, `parallel = FALSE` ignores it.
#' All the `P * nb_simulations` simulations (with `P` stochastic maps per trait) are run as one set of tasks, so the workers are used evenly.
#'
#' The output is exactly the same as without parallelization.
#' The messages, warnings and printed output of the simulations are relayed to the session when the tasks are collected (so they appear after the simulations, not live),
#' and with `verbose = TRUE` one line is printed when each simulation is done.
#'
#' Setting up the workers and sending the data to them has a time cost, so **parallelization only pays off with large backbone phylogenies, or many simulations or maps**.
#' As a rule of thumb, use parallel computation when the sequential run would last more than a minute or so, i.e., for large backbone phylogenies (500 tips or more) and several grafted phylogenies
#' (`nb_maps * nb_simulations`), or for 50 or more grafted phylogenies on small backbones.
#'
#' ### Temporary saves and resuming
#'
#' With `save_temp = TRUE`, each grafted phylogeny (one per simulation and per stochastic map) is saved in its own `.rds` file, built from `temp_file_path` by adding
#' the suffixes `"_map_<name>"` when a `"multiSimmap"` is used, then `"_sim_<k>"` when `nb_simulations > 1`.
#' For instance, with `temp_file_path = "dir/temp.rds"`:
#' * a single phylogeny uses `"dir/temp.rds"`;
#' * several simulations use `"dir/temp_sim_1.rds"`, `"dir/temp_sim_2.rds"`, ...;
#' * several maps use `"dir/temp_map_A.rds"`, `"dir/temp_map_B.rds"`, ...;
#' * several simulations across several maps use both suffixes as `"dir/temp_map_A_sim_1.rds"`;
#' * the prefix `"map_"` is not repeated if the name of the map already starts with it: maps named `"map_1"`, `"map_2"`, ... (default names of a `"multiSimmap"` without names) use `"dir/temp_map_1_sim_1.rds"`, `"dir/temp_map_2_sim_1.rds"`, ...
#'
#' With `force_simulations = TRUE`, `<k>` is the number of the attempt, not necessarily the ID to which the grafted phylogeny will be attributed to in the final output.
#' This is required to identify the appropriate seed used to produce each attempt; attempts that failed leave a file too for tractability.
#'
#' Temporary files can be read with `readRDS()` and used to plot the grafted phylogeny in its current state.
#' The object is a regular `"phylo"` with the attribute `"phylograftR_temp"` that stores what is needed to resume a run.
#' Once the phylogeny is fully grafted, the file is replaced by the complete output of the simulation (with `$simmap`, `$grafted_tips`, etc.).
#'
#' With `clean_temp_files = TRUE` (default), all the temporary files of the run are deleted once the function returns without crashing.
#' **A run that crashes or is interrupted never deletes its files**, so it can be resumed with `resume = TRUE`.
#' Use `clean_temp_files = FALSE` to keep the files after a successful run.
#'
#' With `resume = TRUE`, the function must be called again with **exactly the same arguments** as the interrupted run.
#' If the inputs differ (backbone, list of missing taxa, simmap or other arguments), the function stops and the message lists what differs, with the previous and the new values.
#' For each phylogeny, its temporary file is searched: a complete phylogeny is reused as it is, an unfinished one is continued from its last save, and a phylogeny without file is grafted from scratch.
#'
#' **Resuming is not equivalent to grafting the taxa left on the temporary phylogeny used as a new backbone**, with a `missing_taxa_list` restricted to these taxa:
#' the taxa already grafted would be treated as backbone species, so the sampling fractions of the birth-death models, the stem of the target clades,
#' the identification of the grafted tips and the updated `simmap` would be wrong, and a `simmap` of the initial backbone would not match the temporary phylogeny.
#' Use `resume = TRUE` instead.
#'
#' ### Successive graftings
#'
#' To graft the missing taxa and the missing clades on the same backbone, use [`graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades],
#' which chains both steps with a single set of inputs and a single seed, including with several simulations and/or stochastic maps.
#'
#' The output of this function can also be used as the `phy` of another grafting (e.g., to chain the steps by hand for a single grafted phylogeny).
#'
#' **The grafted tips and branches accumulate**: if `phy` has the `grafted_tips` element of a previous grafting, its flagged tips are kept as grafted
#' (they are still treated as backbone species by the grafting itself, e.g., to count the species used for the sampling fractions),
#' and the new taxa are added to them.
#' `grafted_tips` and `grafted_edges` of the output thus identify **all** the grafted taxa, not only those of the last step.
#' A warning is sent if `grafted_tips` does not match the tips of `phy` (e.g., after `ape::drop.tip()`):
#' use [`identify_grafted()`][identify_grafted] to recompute them first.
#'
#' If a `simmap` is used, provide the updated `simmap` of the previous step (`phy$simmap`, obtained with `return_simmap = TRUE`) to keep the same trait history.
#'
#'
#' @importFrom ape Ntip branching.times extract.clade getMRCA is.ultrametric
#'
#' @importFrom phytools bind.tip nodeHeights nodeheight
#'
#' @importFrom TreeSim corsim
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list] to build a `missing_taxa_list` from a data.frame loaded in R, with clades defined by tips or node IDs.
#'  * [`phylograftR::build_missing_taxa_list_from_taxonomy()`][build_missing_taxa_list_from_taxonomy] to build a `missing_taxa_list` from taxonomic tables, with clades defined by taxonomic names.
#'  * [`phylograftR::build_total_richness_per_clades()`][build_total_richness_per_clades] and [`phylograftR::import_total_richness_per_clades()`][import_total_richness_per_clades] to prepare `total_richness_per_clades` for the unsampled species.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades of backbone clades.
#'  * [`phylograftR::graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades] to graft missing taxa and missing clades in one run, with a single set of inputs.
#'  * [`phylograftR::check_grafting_plan()`][check_grafting_plan] to check a priori the feasibility of a grafting plan, before grafting.
#'  * [`phylograftR::check_grafting_constraints()`][check_grafting_constraints] to check a posteriori that grafted taxa respect their constraints.
#'  * [`phylograftR::extract_grafting_info()`][extract_grafting_info] to extract a posteriori information on the location of grafted lineages.
#'  * [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract a posteriori the richness and the sampling fraction of the clades in each grafted phylogeny.
#'  * [`phylograftR::extract_BD_model_summary()`][extract_BD_model_summary] to extract a posteriori information on the birth-death models fitted during grafting (with `return_BD_model_summary = TRUE`).
#'  * [`phylograftR::plot_grafted()`][plot_grafted] to plot grafted phylogenies.
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
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   # Target clades defined by pairs of tips (MRCA), or by a single tip
#'   # All the other items are named after Taxon_name (same names, same order)
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")),
#'   # Stem of the target clade available for grafting?
#'   # NA uses the global 'include_target_clade_stem' argument (TRUE by default)
#'   # X4 cannot be grafted on the stem
#'   include_target_clade_stem = c(X1 = NA, X2 = NA, X3 = NA, X4 = FALSE),
#'   # Excluded subclades defined by pairs of tips. Only X1 has an excluded subclade
#'   # (the optional name of the excluded subclade, here 'clade_A', is used in messages)
#'   Excluded_subclades = list(X1 = list(clade_A = c("t10", "t19")), X2 = NULL, X3 = NULL, X4 = NULL),
#'   # Stem of the excluded subclades also excluded?
#'   # NA / NULL use the global 'exclude_subclades_stem' argument (TRUE by default).
#'   # Here the stem of the excluded subclade of X1 stays available,
#'   # so X1 can be sister to it, but never nested within it.
#'   exclude_subclades_stem = list(X1 = FALSE, X2 = NULL, X3 = NULL, X4 = NULL))
#'
#' # --- Perform grafting --- #
#'
#' # Graft with the default method ("BD_clades")
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   method = "BD_clades", verbose = FALSE, seed = 1234)
#' plot(grafted_phy, cex = 0.7)
#'
#' # Highlight the grafted tips and branches in red
#' plot_grafted(grafted_phy, cex = 0.7)
#'
#' # Graft with the "uniform" method, without using the excluded subclades
#' grafted_phy_uniform <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                           use_excluded_subclades = FALSE,
#'                                           method = "uniform", verbose = FALSE, seed = 1234)
#'
#' # Check that the constraints are respected
#' check_grafting_constraints(phy_grafted = grafted_phy, phy_backbone = phy,
#'                            missing_taxa_list = missing_taxa_list)
#'
#'
#' # --- Produce several grafted phylogenies --- #
#'
#' # The grafting is stochastic: several simulations can be run at once from the same backbone and list
#' # (a "multiPhylo" named sim_1, sim_2, ...; sim_1 is the same as a single run with the same seed)
#' grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                    method = "uniform", nb_simulations = 3,
#'                                    verbose = FALSE, seed = 1234)
#' names(grafted_phys)
#' # Seeds of the simulations, to reproduce each phylogeny alone
#' attr(grafted_phys, "seeds")
#'
#' # Plot phylos
#' old_par <- par()["mfrow"]
#' par(mfrow = c(1, 3))
#' plot_grafted(grafted_phys, index = 1, cex = 0.7, main = "Sim #1") # sim_1
#' plot_grafted(grafted_phys, index = 2, cex = 0.7, main = "Sim #2") # sim_2
#' plot_grafted(grafted_phys, index = 3, cex = 0.7, main = "Sim #3") # sim_3
#' par(old_par)
#'
#'
#' # --- Account for unsampled species --- #
#'
#' # Account for unsampled species with the total richness of the clades.
#' # All the unique target clades must be listed
#' # (NA = unknown richness, i.e., no unsampled species beyond nested clades)
#' total_richness_per_clades <- list(
#'   Target_clade = list(c("t26", "t18"), c("t8", "t22"), "t5"),
#'   total_richness = c(25, NA, 3))
#' grafted_phy_richness <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                            total_richness_per_clades = total_richness_per_clades,
#'                                            method = "BD_clades", verbose = FALSE, seed = 1234)
#'
#' # Inspect the resulting clade richness and sampling fractions
#' extract_sampling_info(phy_grafted = grafted_phy_richness,
#'                       total_richness_per_clades = total_richness_per_clades)
#'
#'
#' # --- Apply attribute-based exclusions --- #
#'
#' # Restrict the grafting to some states of a stochastic character map (simmap) of the backbone
#' # The toy simmap has three states, "blue", "green" and "red"
#' data(toy_backbone)
#' data(toy_simmap)
#' data(toy_missing_taxa_list_with_allowed_states)
#'
#' # This missing_taxa_list has the item 'allowed_states',
#' # with the states where each taxon can be grafted
#' # (NULL for the taxa without a specific value).
#' toy_missing_taxa_list_with_allowed_states$allowed_states
#'
#' # Perform grafting informed with the simmap and associated 'allowed_states'.
#' # The taxa without taxon-level states use the global 'allowed_states' argument:
#' # here they can only be grafted on "blue" or "green" sections of the branches
#' grafted_phy_states <- graft_missing_taxa(phy = toy_backbone,
#'    missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
#'    simmap = toy_simmap, allowed_states = c("blue", "green"),
#'    method = "BD_clades", verbose = FALSE, seed = 1234)
#'
#' # Plot output phylogeny with colored grafts
#' plot_grafted(grafted_phy_states, use_simmap = FALSE, cex = 0.7)
#'
#' # The updated simmap is in the output:
#' # the grafted taxa are in the state of the point where they were grafted.
#' # Grafted branches are dashed lines, and grafted tips are labelled in red.
#' plot_grafted(grafted_phy_states, fsize = 0.5, ftype = "i",
#'              simmap_colors = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))
#'
#' # States of the grafted taxa: "New_sp1" and "New_sp7" can be "red"
#' # thanks to their taxon-level states, but the taxa using the global states cannot.
#' all_states <- phytools::getStates(grafted_phy_states$simmap, type = "tips")
#' all_states[toy_missing_taxa_list_with_allowed_states$Taxon_name]
#'
#'
#' # --- Use several stochastic maps (alternative hypotheses of trait evolution) as constraints --- #
#'
#' # A 'multiSimmap' gives several stochastic maps of the same trait:
#' # the grafting is run with each of them.
#' # Here, 3 maps are simulated from the states of the tips of the toy simmap.
#' tip_states <- phytools::getStates(toy_simmap, type = "tips")
#'
#' set.seed(118)
#' toy_multiSimmap <- phytools::make.simmap(toy_backbone, tip_states,
#'                                          model = "ER", nsim = 3, message = FALSE)
#'
#' # Graft using constraints from all 3 simmaps X 2 simulations per simmap
#' grafted_phys_maps <- graft_missing_taxa(phy = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
#'     simmap = toy_multiSimmap, allowed_states = c("blue", "green"),
#'     method = "BD_clades", nb_simulations = 2, verbose = FALSE, seed = 1234)
#'
#' # Inspect output:
#' # A list with one element per simmap, each a multiPhylo with the simulations of the map
#' names(grafted_phys_maps)
#' names(grafted_phys_maps$map_2)
#' # The map and the seed of a phylogeny
#' attr(grafted_phys_maps$map_2[[1]], "map")
#' attr(grafted_phys_maps$map_2[[1]], "seed")
#'
#'
#' # --- Constraints from several traits at once --- #
#'
#' # Load the simmap for the second trait: 'toy_simmap_elevation' has two states: "high" and "low".
#' data(toy_simmap_elevation)
#' # Load missing taxa list with '$allowed_states' covering both color and elevation traits
#' data(toy_missing_taxa_list_with_two_traits)
#' # NULL means no restriction for this taxon, for the given trait.
#' toy_missing_taxa_list_with_two_traits$allowed_states
#'
#' # Global 'allowed_states': "blue" or "green" for colors, no restriction on the elevation.
#' # Each taxon can only be grafted on branch pieces in allowed states for both traits.
#' grafted_phy_traits <- graft_missing_taxa(phy = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_two_traits,
#'     simmap = list(color = toy_simmap, elevation = toy_simmap_elevation),
#'     allowed_states = list(c("blue", "green"), NULL),
#'     verbose = FALSE, seed = 1234)
#'
#' # '$simmap' is a list of updated simmaps, one per trait
#' names(grafted_phy_traits$simmap)
#'
#'
#' # --- Constraints from biogeographic ranges and another trait at once --- #
#'
#' # One of the traits can be treated as biogeographic ranges:
#' # the allowed ranges do not have to match exactly the range of target stem branch pieces,
#' # any allowed range comprised in a target range is valid (e.g., "A" can be grafted on "AB").
#' # Here, the first trait holds ranges.
#'
#' # Load datasets
#' data(toy_simmap_ranges) # Ranges simmap
#' data(toy_missing_taxa_list_with_ranges_and_elevation) # Taxon-level constraints
#' toy_missing_taxa_list_with_ranges_and_elevation$allowed_states
#'
#' # Global 'allowed_states': ranges containing "A" or "B", no restriction on the elevation.
#' grafted_phy_ranges_traits <- graft_missing_taxa(phy = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
#'     simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
#'     allowed_states = list(c("A", "B"), NULL),
#'     treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
#'     verbose = FALSE, seed = 1234)
#'
#' # '$simmap' is a list of updated simmaps, one per trait: ranges, and elevation
#' names(grafted_phy_ranges_traits$simmap)
#'
#' # Check the constraints of both traits (columns 'states_ok_<trait>' and 'ranges_ok')
#' check_grafting_constraints(phy_grafted = grafted_phy_ranges_traits, phy_backbone = toy_backbone,
#'     missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
#'     allowed_states = list(c("A", "B"), NULL),
#'     treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE))
#'
#'


graft_missing_taxa <- function (phy, # Backbone phylogeny
                                missing_taxa_list, # List with Taxon_name = c(), and items named after Taxon_name: Target_clade = list(X1 = c(), ...), include_target_clade_stem = c(X1 = TRUE, ...) (optional), Excluded_subclades = list(X1 = list(c(), c()), X2 = NULL, ...), exclude_subclades_stem = list(X1 = c(TRUE, NA), X2 = NULL, ...) (optional)
                                method = "BD_clades", # "uniform" or "BD_clades"
                                min_tips = 5, # Minimum number of tips of the clade used to fit the birth-death model for "BD_clades" (the enclosing clade is used if the focal clade is smaller)
                                return_BD_model_summary = FALSE, # Report the birth-death models fitted during grafting in the attribute "BD_model_summary" of the phylogenies? Only with "BD_clades"
                                include_target_clade_stem = TRUE, # Include the stem of the target clades when grafting? Default for all taxa, overridden by missing_taxa_list$include_target_clade_stem
                                use_excluded_subclades = TRUE, # Avoid grafting within Excluded_subclades?
                                exclude_subclades_stem = TRUE, # Also exclude the stem of the excluded subclades? Default for all subclades, overridden by missing_taxa_list$exclude_subclades_stem
                                total_richness_per_clades = NULL, # NULL, or list with Target_clade = list(c(), ...) and total_richness = c(...) (NA allowed). Total number of species of the clades
                                simmap = NULL, # NULL, or stochastic character map(s) (class "simmap" or "multiSimmap", or a list of them for several traits)  of the backbone: restricts grafting to the branch segments in allowed states
                                allowed_states = NULL, # States of the simmap where grafting is allowed. Default for all taxa, overridden by missing_taxa_list$allowed_states. NULL: all states
                                treat_states_as_ranges = FALSE, # Are the states of the simmap biogeographic ranges to be treated as sets of areas? An allowed range matches the states that contain it. (only with a simmap)
                                return_simmap = TRUE, # Add the updated simmap to the output ($simmap)? Only with a simmap
                                min_valid_fraction = 0.1, # Warn for the taxa with a smaller share of valid branches in allowed states: fraction of their available branch length ('uniform') or of their time window ('BD_clades') (only with a simmap)
                                max_age_redraws = 1000, # Maximum number of redraws of the divergence ages of a group of taxa when some ages cannot be used because of the states (only with a simmap and 'BD_clades')
                                seed = 1234,
                                nb_simulations = 1, # Number of grafted phylogenies to produce in the same run (a multiPhylo if > 1)
                                force_simulations = FALSE, # If TRUE, try new seeds until nb_simulations phylogenies are obtained (up to 10 x nb_simulations attempts)
                                parallel = NULL, # NULL: run in parallel if a 'future' plan with several workers is set; FALSE: sequential; TRUE: parallel (stops if not possible). Packages 'future' and 'future.apply' (suggested)
                                save_temp = FALSE, # Save the grafted phylogenies during the process (partial every 'save_temp_every' taxa or groups of taxa, and complete once fully grafted)
                                save_temp_every = NULL, # NULL: only save the complete phylogenies
                                temp_file_path = NULL, # Default = "temp_phy_from_grafting.rds" in the working directory
                                clean_temp_files = TRUE, # Delete all the temporary files once the run is completed without crashing?
                                resume = FALSE, # Resume an interrupted run from its temporary files (same arguments)
                                ultrametric_tol = .Machine$double.eps^0.5, # Tolerance used by ape::is.ultrametric() for the ultrametricity checks
                                verbose = TRUE) # To print progress every 100 taxa or every 10 target clades
{
  ### ---- Check input validity ---- ###
  {
    ## Check 'nb_simulations' and 'force_simulations'
    check_nb_simulations(nb_simulations = nb_simulations, force_simulations = force_simulations)
    check_parallel(parallel)
    check_save_temp(save_temp = save_temp, save_temp_every = save_temp_every, resume = resume, clean_temp_files = clean_temp_files)

    ## Check 'method'
    if (!(method %in% c("uniform", "BD_clades")))
    {
      stop("'method' must be either 'uniform' or 'BD_clades'")
    }

    ## The backbone must be ultrametric (error otherwise). It is never modified: making it ultrametric is left to the user
    check_backbone_ultrametric(phy = phy, tol = ultrametric_tol)

    ## The backbone must be binary (error otherwise): polytomies are not handled, and resolving them is left to the user
    check_backbone_binary(phy = phy)

    ## Check 'min_tips'
    if (!is.numeric(min_tips) || length(min_tips) != 1 || min_tips < 3)
    {
      stop("'min_tips' must be a single number >= 3")
    }

    ## Check 'return_BD_model_summary'
    if (!is.logical(return_BD_model_summary) || length(return_BD_model_summary) != 1 || is.na(return_BD_model_summary))
    {
      stop("'return_BD_model_summary' must be TRUE or FALSE")
    }
    if (return_BD_model_summary && method == "uniform")
    {
      stop("'return_BD_model_summary = TRUE' requires method = 'BD_clades': no birth-death model is fitted with method = 'uniform'")
    }

    ## The backbone tree must have at least 'min_tips' tips, so a BD model can be fit
    if (method == "BD_clades" && length(phy$tip.label) < min_tips)
    {
      stop("Backbone phylogeny too small to fit a birth-death model (< ", min_tips, " tips).\n",
           "Choose method = 'uniform', lower 'min_tips', or provide a larger backbone phylogeny.")
    }

    ## Check the structure of missing_taxa_list: Taxon_name, and the items provided (same length and names identical to Taxon_name, in the same order)
    validate_missing_taxa_list(missing_taxa_list, clade_item = "Target_clade")
    n_missing <- length(missing_taxa_list$Taxon_name)

    ## Stem availability of each taxon: 'missing_taxa_list$include_target_clade_stem' of the taxon if provided (not NA), otherwise the global 'include_target_clade_stem' argument
    # The resolved vector (no NA) replaces the item in the list, so it follows the taxa when they are reordered
    stem_requested <- missing_taxa_list$include_target_clade_stem # Values provided by the user, before resolution
    missing_taxa_list$include_target_clade_stem <- resolve_include_target_clade_stem(
      include_target_clade_stem_taxa = stem_requested, include_target_clade_stem = include_target_clade_stem, n_taxa = n_missing,
      taxon_names = missing_taxa_list$Taxon_name)

    ## Stem exclusion of each excluded subclade: 'missing_taxa_list$exclude_subclades_stem' of the subclade if provided (not NA), otherwise the global 'exclude_subclades_stem'
    # The resolved list (one logical vector per taxon, no NA) replaces the item in the list, so it follows the taxa when they are reordered
    missing_taxa_list$exclude_subclades_stem <- resolve_exclude_subclades_stem(
      exclude_subclades_stem_taxa = missing_taxa_list$exclude_subclades_stem, exclude_subclades_stem = exclude_subclades_stem,
      excluded_subclades = missing_taxa_list$Excluded_subclades, n_taxa = n_missing,
      taxon_names = missing_taxa_list$Taxon_name)

    ## Species attributes: 'simmap', 'allowed_states' and 'min_valid_fraction'
    # The states of each taxon are resolved below, once the simmap is validated and the backbone labelled
    use_states <- !is.null(simmap)
    simmap_input <- check_simmap_input(simmap) # The maps (a single one, or the maps of a 'multiSimmap'), their names and whether the output is a list of maps
    if (!is.numeric(min_valid_fraction) || length(min_valid_fraction) != 1 || is.na(min_valid_fraction) || min_valid_fraction < 0 || min_valid_fraction > 1)
    {
      stop("'min_valid_fraction' must be a single number between 0 and 1")
    }
    if (!is.logical(return_simmap) || length(return_simmap) != 1 || is.na(return_simmap))
    {
      stop("'return_simmap' must be TRUE or FALSE")
    }
    # 'treat_states_as_ranges': TRUE or FALSE for a single map (or no map), with the maps of several traits a logical vector with one value per trait (one TRUE at most: the range trait)
    range_trait <- NULL # Position of the trait treated as ranges (several traits)
    if (is.null(simmap_input$traits))
    {
      if (!is.logical(treat_states_as_ranges) || length(treat_states_as_ranges) != 1 || is.na(treat_states_as_ranges))
      {
        stop("'treat_states_as_ranges' must be TRUE or FALSE (a logical vector with one value per trait is used with the maps of several traits)")
      }
    } else {
      range_trait <- resolve_range_trait(treat_states_as_ranges, trait_names = simmap_input$traits, traits_named = simmap_input$traits_named)
      treat_states_as_ranges <- FALSE # With several traits, the range trait is handled through 'range_trait'
    }
    if (!is.numeric(max_age_redraws) || length(max_age_redraws) != 1 || is.na(max_age_redraws) || max_age_redraws < 1 || max_age_redraws != round(max_age_redraws))
    {
      stop("'max_age_redraws' must be a single positive integer")
    }
    if (use_states)
    {
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

    ## Warn (once) about taxa with a single-tip target clade and 'include_target_clade_stem = FALSE'
    # The terminal branch of the tip is the only branch of a single-tip clade (it is its stem): not including the stem would lead to an error and is ignored.
    # Only values explicitly set to FALSE are flagged (not NA with a global 'include_target_clade_stem = FALSE', which is documented as ignored for single tips)
    if (!is.null(stem_requested))
    {
      is_single_tip_target <- vapply(missing_taxa_list$Target_clade, function (tips) length(tips) == 1, logical(1))
      single_tip_ignored <- which(is_single_tip_target & !is.na(stem_requested) & !as.logical(stem_requested))
      if (length(single_tip_ignored) > 0)
      {
        n_show <- min(length(single_tip_ignored), 10)
        warning(length(single_tip_ignored), " taxa have a single-tip target clade with 'include_target_clade_stem = FALSE': ",
                "this value was ignored for the whole run, and the stem was used. The terminal branch of the tip is the only branch of a single-tip clade.\n",
                paste0("  - ", missing_taxa_list$Taxon_name[single_tip_ignored[seq_len(n_show)]],
                       " (target tip: ", vapply(missing_taxa_list$Target_clade[single_tip_ignored[seq_len(n_show)]], paste, character(1), collapse = ", "), ")", collapse = "\n"),
                if (length(single_tip_ignored) > n_show) paste0("\n  ... and ", length(single_tip_ignored) - n_show, " more") else "",
                call. = FALSE)
      }
    }

    ## Warn (once) about single-tip excluded subclades with a stem that is not excluded
    # A single tip has no branch but its terminal branch, which is its stem: the exclusion has no effect
    if (use_excluded_subclades && !is.null(missing_taxa_list$Excluded_subclades))
    {
      single_tip_noop <- do.call(rbind, lapply(seq_len(n_missing), function (i)
      {
        subclades_i <- missing_taxa_list$Excluded_subclades[[i]]
        noop <- which(vapply(subclades_i, length, integer(1)) == 1 & !missing_taxa_list$exclude_subclades_stem[[i]])
        if (length(noop) == 0) return(NULL)
        data.frame(Taxon_name = missing_taxa_list$Taxon_name[i], tip = label_excluded_subclades(subclades_i[noop]), stringsAsFactors = FALSE)
      }))
      if (!is.null(single_tip_noop))
      {
        n_show <- min(nrow(single_tip_noop), 10)
        warning(nrow(single_tip_noop), " single-tip excluded subclade(s) have 'exclude_subclades_stem = FALSE': they exclude nothing, ",
                "as the terminal branch of a tip is its only branch (its stem).\n",
                paste0("  - ", single_tip_noop$Taxon_name[seq_len(n_show)], " (excluded subclade: ", single_tip_noop$tip[seq_len(n_show)], ")", collapse = "\n"),
                if (nrow(single_tip_noop) > n_show) paste0("\n  ... and ", nrow(single_tip_noop) - n_show, " more") else "",
                call. = FALSE)
      }
    }

    ## Check that taxa names are truly missing from the backbone
    if (any(missing_taxa_list$Taxon_name %in% phy$tip.label))
      stop("Taxon_name must be absent from the backbone phylogeny. Present: ",
           paste(intersect(missing_taxa_list$Taxon_name, phy$tip.label), collapse = ", "))

    ## Check that all tips used to define target clades are present in the backbone
    all_target_clade_tips <- unique(unlist(missing_taxa_list$Target_clade))
    if (!all(all_target_clade_tips %in% phy$tip.label))
    {
      stop("Some tip labels provided to locate MRCA nodes of target clades are missing from the backbone phylogeny\n",
           "Missing tips: ", paste(setdiff(all_target_clade_tips, phy$tip.label), collapse = ", "))
    }

    # Check that tips defining the excluded subclades are all found on the phylogeny
    all_excluded_tips <- unique(unlist(missing_taxa_list$Excluded_subclades))
    if (!all(all_excluded_tips %in% phy$tip.label))
    {
      stop("Some tip labels provided to locate nodes of excluded subclades are missing from the backbone phylogeny\n",
           "Missing tips: ", paste(setdiff(all_excluded_tips, phy$tip.label), collapse = ", "))
    }

    ## Check the position of the excluded subclades relative to their target clade
    # Excluded subclades are expected to be nested within their target clade. Three other cases exists and are problematic:
    #  - Excluded subclade containing the target clade: all branches of the target clade, stem included, are excluded.
    #    No branch is available, which would stop the grafting anyway: stop now with a message pointing to the problem.
    #  - Excluded subclade identical to the target clade: the whole crown is excluded. No branch is left if the stem is not available,
    #    and only the stem branch is left if the stem is available. Grafting only on a stem branch is not what this function is designed for:
    #    with 'BD_clades', the divergence age is drawn before applying the exclusions and would most likely fall within the excluded crown,
    #    where no branch is available, so the run would stop in the middle of the process.
    #    Stop now, and point to graft_missing_clades(), which is designed to graft specifically on a stem branch (of a clade or of a single tip).
    #  - Excluded subclade disjoint from the target clade: the exclusion has no effect. Harmless, so only a warning.
    if (use_excluded_subclades && !is.null(missing_taxa_list$Excluded_subclades))
    {
      excluded_position <- check_excluded_subclades_position(phy = phy, missing_taxa_list = missing_taxa_list)

      # Helper to format the list of problematic taxa (max 10 shown)
      format_excluded_position <- function (df)
      {
        n_show <- min(nrow(df), 10)
        paste0(paste0("  - ", df$Taxon_name[seq_len(n_show)], " (target clade tips: ", df$Target_clade[seq_len(n_show)],
                      "; excluded subclade: ", df$Excluded_subclade[seq_len(n_show)], ")", collapse = "\n"),
               if (nrow(df) > n_show) paste0("\n  ... and ", nrow(df) - n_show, " more") else "")
      }

      contains_target <- excluded_position[excluded_position$relation == "contains_target", ]
      if (nrow(contains_target) > 0)
      {
        stop(nrow(contains_target), " excluded subclade(s) contain the whole target clade of their taxon, so all its branches would be excluded.\n",
             format_excluded_position(contains_target),
             "\nExcluded subclades must be nested within their target clade. Check the tips defining 'Target_clade' and 'Excluded_subclades'.")
      }

      identical_clade <- excluded_position[excluded_position$relation == "identical", ]
      if (nrow(identical_clade) > 0)
      {
        stop(nrow(identical_clade), " excluded subclade(s) are identical to the target clade of their taxon, so the whole clade is excluded: ",
             "no branch would be available for grafting without the stem, and only the stem branch with the stem included.\n",
             format_excluded_position(identical_clade),
             "\nTo graft taxa specifically on the stem branch of a clade (or of a single tip), use graft_missing_clades() ",
             "with this clade as 'Sister_clade'. Otherwise, check the tips defining 'Target_clade' and 'Excluded_subclades'.")
      }

      disjoint <- excluded_position[excluded_position$relation == "disjoint", ]
      if (nrow(disjoint) > 0)
      {
        warning(nrow(disjoint), " excluded subclade(s) are outside the target clade of their taxon: these exclusions have no effect.\n",
                format_excluded_position(disjoint),
                "\nExcluded subclades are expected to be nested within their target clade. Check the tips defining 'Target_clade' and 'Excluded_subclades'.",
                call. = FALSE)
      }
    }

    ## temp_file_path
    # If provided, must be a single path ending with ".rds", in an existing and writable directory
    check_temp_file_path(temp_file_path = temp_file_path)
    # Absolute path (in the main session), so the workers of a parallel run write in the same place
    if (save_temp) { temp_file_path <- resolve_temp_file_path(temp_file_path) }

    ## Richness of the clades (unsampled species)
    # Validates 'total_richness_per_clades' if provided, and derives the number of unsampled species of each node of the backbone.
    # Always computed: it also provides the number of missing taxa to include in the sampling fraction of the BD models
    richness_info <- build_richness_info(phy = phy, missing_taxa_names = missing_taxa_list$Taxon_name,
                                         missing_clades = missing_taxa_list$Target_clade,
                                         mode = "target", total_richness_per_clades = total_richness_per_clades,
                                         missing_clade_names = missing_taxa_list[["Target_clade_name"]])

    if (!is.null(total_richness_per_clades) && method == "uniform")
    {
      warning("'total_richness_per_clades' is not used during grafting when method = 'uniform' as no diversification process is modeled. ",
              "Use extract_sampling_info() on the output to get the resulting sampling fractions of the clades.",
              call. = FALSE)
    }
  }

  ## Label the internal nodes of the backbone, to identify them in the grafted tree (node IDs change after each graft) and restore the original node labels at the end
  labelled_backbone <- label_backbone_nodes(phy)
  phy <- labelled_backbone$phy

  ### ---- Prepare the maps ---- ###

  # Everything that depends on the map (the tracker of the states, the resolution of the allowed states, and the checks of the branches available) is prepared before any simulation,
  # for each map: any problem stops the run at once, and the message states the map when several maps are used
  prepare_map <- function (simmap_k)
  {
    ## With a simmap: validate it, start following the states of the branches through grafting, and resolve the allowed states of each taxon
    # The resolved list (one character vector per taxon, no NULL) replaces the item in the list, so it follows the taxa when they are reordered
    tracker <- NULL
    if (use_states)
    {
      if (inherits(simmap_k, "phylograftR_trait_maps"))
      {
        # Maps of several traits: merged in a composite simmap (one state per combination of states), the allowed states of each taxon are resolved among the composite states
        merged <- merge_trait_simmaps(phy = phy, trait_maps = simmap_k)
        tracker <- start_state_tracking(phy = phy, simmap = merged$simmap)
        tracker$traits <- merged$traits # To split the updated simmap in one simmap per trait at the end
        if (is.null(range_trait))
        {
          missing_taxa_list$allowed_states <- resolve_trait_allowed_states(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                                           traits = merged$traits, composite_states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
        } else {
          # One trait is treated as biogeographic ranges: acceptable composite states, candidate ranges to inherit (the other traits are inherited as they are) and inherited composite states
          resolved <- resolve_trait_allowed_ranges(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states, traits = merged$traits,
                                                   composite_states = tracker$states, taxon_names = missing_taxa_list$Taxon_name, range_trait = range_trait)
          missing_taxa_list$allowed_states <- resolved$allowed_states
          missing_taxa_list$range_candidates <- resolved$range_candidates
          tracker <- extend_tracker_states(tracker = tracker, extra_states = resolved$inherited_states)
        }
      } else if (treat_states_as_ranges) {
        tracker <- start_state_tracking(phy = phy, simmap = simmap_k)
        # Allowed ranges are resolved into the acceptable states (the states containing an allowed range), used to restrict grafting, and the candidate ranges of each state
        resolved <- resolve_allowed_ranges(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                           states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
        missing_taxa_list$allowed_states <- resolved$allowed_states # Allowed states adjusted to the ones observed in the simmap
        missing_taxa_list$range_candidates <- resolved$range_candidates # List of ranges that can be inherited by the taxa, based on which allowed state it is grafted on
        tracker <- extend_tracker_states(tracker = tracker, extra_states = resolved$inherited_states) # Ranges inherited that are not states of the simmap (need to be added to update the next allowed states)
      } else {
        tracker <- start_state_tracking(phy = phy, simmap = simmap_k)
        missing_taxa_list$allowed_states <- resolve_allowed_states(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                                   states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
      }
    }

    ## Check upfront that every taxon has at least one branch available for grafting (and, with a simmap, a piece of branch in its allowed states)
    # Better to stop before grafting than in the middle of a long run
    no_branch_taxa <- find_taxa_without_allowed_branch(phy = phy, missing_taxa_list = missing_taxa_list,
                                                       include_target_clade_stem = missing_taxa_list$include_target_clade_stem,
                                                       exclude_subclades_stem = missing_taxa_list$exclude_subclades_stem,
                                                       use_excluded_subclades = use_excluded_subclades,
                                                       tracker = tracker, allowed_states = missing_taxa_list$allowed_states)
    valid_fraction <- attr(no_branch_taxa, "valid_fraction")
    excluded_taxa <- no_branch_taxa[no_branch_taxa$Cause == "excluded", , drop = FALSE]
    if (nrow(excluded_taxa) > 0)
    {
      n_show <- min(nrow(excluded_taxa), 10)
      stop(nrow(excluded_taxa), " taxa have no branch available for grafting: all branches of their target clade are excluded.\n",
           paste0("  - ", excluded_taxa$Taxon_name[seq_len(n_show)], " (target clade tips: ", excluded_taxa$Target_clade[seq_len(n_show)], ")", collapse = "\n"),
           if (nrow(excluded_taxa) > n_show) paste0("\n  ... and ", nrow(excluded_taxa) - n_show, " more") else "",
           "\nAdjust 'include_target_clade_stem' (global or in missing_taxa_list), 'missing_taxa_list$Excluded_subclades' or 'exclude_subclades_stem' for these taxa.")
    }
    no_state_taxa <- no_branch_taxa[no_branch_taxa$Cause == "allowed_states", , drop = FALSE]
    if (nrow(no_state_taxa) > 0)
    {
      n_show <- min(nrow(no_state_taxa), 10)
      stop(nrow(no_state_taxa), " taxa have no branch available for grafting in their allowed states: none of the available branches of their target clade has a piece in these states.\n",
           paste0("  - ", no_state_taxa$Taxon_name[seq_len(n_show)], " (target clade tips: ", no_state_taxa$Target_clade[seq_len(n_show)],
                  "; allowed states: ", vapply(missing_taxa_list$allowed_states[no_state_taxa$Taxon_name[seq_len(n_show)]], paste, character(1), collapse = ", "), ")", collapse = "\n"),
           if (nrow(no_state_taxa) > n_show) paste0("\n  ... and ", nrow(no_state_taxa) - n_show, " more") else "",
           "\nAdjust 'allowed_states' (global or in missing_taxa_list), or the constraints on the branches ('include_target_clade_stem', 'Excluded_subclades', 'exclude_subclades_stem') for these taxa.")
    }

    ## Warn (once) about the taxa with a very small share of valid branches in their allowed states
    # 'uniform': fraction of the length of the available branches in allowed states. 'BD_clades': fraction of the time window of the divergence ages with a branch in allowed states
    if (use_states)
    {
      share_in_states <- if (method == "uniform") valid_fraction else attr(no_branch_taxa, "time_coverage")
      low_fraction <- which(!is.na(share_in_states) & share_in_states < min_valid_fraction)
      if (length(low_fraction) > 0)
      {
        n_show <- min(length(low_fraction), 10)
        low_taxa <- names(share_in_states)[low_fraction[order(share_in_states[low_fraction])][seq_len(n_show)]]
        what_share <- if (method == "uniform") "of the length of their available branches in their allowed states" else "of the time window of their divergence ages covered by available branches in their allowed states"
        what_effect <- if (method == "uniform") "their grafting positions are concentrated on a small part of their target clade" else "their divergence ages will be concentrated on a narrow time window, or many ages will need to be redrawn"
        warning(length(low_fraction), " taxa have less than ", format(100 * min_valid_fraction), "% ", what_share, " ",
                "('min_valid_fraction' = ", min_valid_fraction, "): ", what_effect, ".\n",
                paste0("  - ", low_taxa, " (target clade tips: ", vapply(missing_taxa_list$Target_clade[low_taxa], paste, character(1), collapse = ", "),
                       "; ", sprintf("%.2f%%", 100 * share_in_states[low_taxa]), " in allowed states)", collapse = "\n"),
                if (length(low_fraction) > n_show) paste0("\n  ... and ", length(low_fraction) - n_show, " more") else "",
                call. = FALSE)
      }
    }

    list(tracker = tracker, missing_taxa_list = missing_taxa_list)
  }
  prepared_maps <- lapply(seq_along(simmap_input$maps), function (k)
  {
    with_map_context(prepare_map(simmap_input$maps[[k]]), map_name = simmap_input$names[k])
  })

  ### ---- Temporary saves ---- ###

  # Tags of the maps in the names of the temporary files (maps never share a file)
  map_tags <- temp_map_tags(simmap_input$names)

  # Warn once, before any simulation, if 'save_temp_every' is larger than the number of units to graft (taxa, or groups of taxa): nothing would be saved before the phylogeny is complete
  if (save_temp && !is.null(save_temp_every))
  {
    n_save_units <- vapply(prepared_maps, function (prepared)
    {
      if (method == "uniform") { return(n_missing) }
      length(group_and_order_clades(phy, prepared$missing_taxa_list$Target_clade, include_target_clade_stem = prepared$missing_taxa_list$include_target_clade_stem,
                                    allowed_states = if (use_states) prepared$missing_taxa_list$allowed_states else NULL)$group_clade_tips)
    }, numeric(1))
    warn_save_temp_every(save_temp_every = save_temp_every, n_units = n_save_units, unit = if (method == "uniform") "taxa" else "groups of taxa", map_names = simmap_input$names)
  }

  ### ---- Simulations ---- ###

  # Taxa already grafted on the input phylogeny, if it is the output of a previous grafting (see 'Successive graftings')
  previously_grafted <- get_previously_grafted(phy)

  # Function to produce one simulation = one grafted phylogeny, from the backbone, the tracker for states and the list of taxa of a map (prepared above), and the associated seed and index
  simulate_one <- function (tracker, missing_taxa_list, map_name, sim_seed, sim_index, verbose_sim)
  {
    phy <- phy
    verbose <- verbose_sim
    # With several simulations, the progress messages state the simulation they belong to
    # (and the map, if the maps of a 'multiSimmap' are used)
    sim_label <- paste0(if (!is.null(map_name)) paste0(" - Map ", map_name), if (nb_simulations > 1) paste0(" - Simulation n\u00b0", sim_index), collapse = "")
    temp_file_path <- temp_file_path_simulation(temp_file_path = temp_file_path, index = sim_index, multiple = nb_simulations > 1, map_name = map_name,
                                                map_tag = if (is.null(map_name)) NULL else map_tags[[map_name]])

    ## Temporary saves: signature of the inputs of this simulation (stored in its files), and resume from its file if requested
    signature <- NULL
    resumed <- NULL
    if (save_temp)
    {
      signature <- temp_signature(phy = phy, missing_taxa_list = missing_taxa_list, tracker = tracker,
                                  options = list(method = method, use_excluded_subclades = use_excluded_subclades, treat_states_as_ranges = treat_states_as_ranges, min_tips = min_tips, return_BD_model_summary = return_BD_model_summary,
                                                 max_age_redraws = max_age_redraws, ultrametric_tol = ultrametric_tol, richness = richness_info$table,
                                                 sim_seed = sim_seed, sim_index = sim_index, map_name = map_name))
      if (resume)
      {
        resumed <- read_resume_state(path = temp_file_path, signature = signature)
        # A complete phylogeny is reused as it is
        if (!is.null(resumed) && identical(resumed$status, "complete"))
        {
          if (verbose) cat(paste0(Sys.time(), sim_label, " - Complete phylogeny found in the temporary file: reused\n"))
          return(resumed$phy)
        }
        if (verbose && is.null(resumed)) cat(paste0(Sys.time(), sim_label, " - No temporary file found: grafting from scratch\n"))
      }
    }

    ## Set seed of the simulation, for the duration of the simulation only (the global RNG state is restored on exit)
    set_local_seed(sim_seed)

    ## Randomize the order of taxa so it has no influence on the pattern of grafting
    # Taxa of the same target clade are exchangeable: their order is random.
      # With method = 'uniform', taxa are then re-ordered by target clade so nested clades are processed first (the random order is kept within each clade).
      # With method = 'BD_clades', taxa are grouped and ordered by clade within the loop.
    random_order <- sample.int(n_missing)
    missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[random_order] })

    # Root age: constant during grafting (all new tips reach the present) => computed once
    root_age <- max(phytools::nodeHeights(tree = phy)[, 2])

    ## ---- For "uniform" strategy: Loop per taxon, nested target clades first ---- ##
    if (method == "uniform")
    {
      # Order taxa by target clade, nested clades first (random order kept within each clade).
      # Rationale: taxa of a larger clade are grafted on a tree that already includes the taxa of the clades nested within it,
      # so they can be grafted on the branches of those taxa. This is the same rule as for the 'BD_clades' method.
      # The order of taxa within each target clade remains random.
      groups_uniform <- group_and_order_clades(phy, missing_taxa_list$Target_clade, include_target_clade_stem = missing_taxa_list$include_target_clade_stem)
      missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[groups_uniform$taxa_order] })

      if (verbose) cat(paste0(Sys.time(), sim_label, " - Start grafting process\n"))

      # Resume: restart from the phylogeny, the states and the random number generator saved with the last progress (the taxa already grafted are skipped)
      first_taxon <- 1
      if (!is.null(resumed))
      {
        phy <- resumed$phy
        tracker <- resumed$state$tracker
        set_rng_state(resumed$state$rng)
        first_taxon <- resumed$state$progress$done + 1
        if (verbose) cat(paste0(Sys.time(), sim_label, " - Resuming from the temporary file: ", first_taxon - 1, "/", n_missing, " taxa already grafted\n"))
      }

      for (i in seq_len(n_missing)[seq_len(n_missing) >= first_taxon])
      {
        focal_taxon_name <- missing_taxa_list$Taxon_name[i]
        # Is the stem of the target clade available for this taxon?
        include_target_clade_stem_i <- missing_taxa_list$include_target_clade_stem[i]

        # Get node associated to tips
          # 1 tip => tip node
          # 2 tips => MRCA node
        MRCA_node <- get_clade_node(phy, missing_taxa_list$Target_clade[[i]])

        # Get the branches of the excluded subclades: the branches within each subclade, and its stem region if its stem is excluded
        # (single tips exclude their terminal branch only if their stem is excluded)
        excluded_nodes <- NULL
        if (use_excluded_subclades)
        {
          excluded_nodes <- get_excluded_branch_nodes(phy = phy, excluded_subclades_tips_list = missing_taxa_list$Excluded_subclades[[i]],
                                                      exclude_subclades_stem = missing_taxa_list$exclude_subclades_stem[[i]],
                                                      missing_tip_labels = missing_taxa_list$Taxon_name)
        }

        # Warn if the target clade is the whole tree: it has no stem branch, so only branches inside the clade can be used
        if (include_target_clade_stem_i && length(missing_taxa_list$Target_clade[[i]]) > 1 && is.na(get_parent_node(phy, MRCA_node)))
        {
          warning("Target clade of '", focal_taxon_name, "' defined by ", paste(missing_taxa_list$Target_clade[[i]], collapse = ", "),
                  " is the whole tree: it has no stem branch, so only branches inside the clade are used")
        }

        # With a simmap: Get ages of the branches of the current tree, to cut them in pieces in the allowed states of the taxon, and follow the states through grafting
        edge_ages_i <- if (use_states) root_age - phytools::nodeHeights(phy) else NULL

        # Select branch for grafting using uniform distribution weighted by branch length to ensure a flat density over any time point along branches
        # (no diversification model). With a simmap, only the segments in allowed states are used
        # The draw returns the selected branch (tipward node) and the position of the graft on it (distance from its tipward node)
        drawn_point <- draw_random_grafting_position_along_available_branch_segments(
           phy = phy, MRCA_node = MRCA_node, include_target_clade_stem = include_target_clade_stem_i,
           excluded_nodes = excluded_nodes,
           missing_tip_labels = missing_taxa_list$Taxon_name, seed = NULL,
           tracker = tracker, allowed_states = if (use_states) missing_taxa_list$allowed_states[[i]] else NULL, edge_ages = edge_ages_i)

        # Stop with an informative message if no branch is available (all branches excluded by the user's constraints)
        if (is.null(drawn_point))
        {
          stop("No branch available to graft '", focal_taxon_name, "'.\n",
               "All branches of its target clade (tips: ", paste(missing_taxa_list$Target_clade[[i]], collapse = ", "), ") are excluded",
               if (use_states) " or not in its allowed states" else "", ".\n",
               "Adjust 'include_target_clade_stem' (global or in missing_taxa_list), 'missing_taxa_list$Excluded_subclades' or 'exclude_subclades_stem' (global or in missing_taxa_list) for this taxon.")
        }

        selected_node_ID <- drawn_point$node
        # Position of the graft: drawn from a uniform distribution along the selected segment
        position <- drawn_point$position

        # Compute the new branch length to reach present time => make the tree ultrametric
        new_branch_length <- position + (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))
        # Graft new branch
        phy_before_graft <- phy
        phy <- phytools::bind.tip(tree = phy, tip.label = focal_taxon_name, edge.length = new_branch_length,
                                  where = selected_node_ID, position = position)
        # With a simmap: the new taxon takes the state of the point where it was grafted, and the new node is labelled to keep following the branches
        if (use_states)
        {
          tracked <- track_graft(tracker = tracker, phy_before = phy_before_graft, phy_after = phy, node = selected_node_ID, position = position,
                                 tip_label = focal_taxon_name, allowed_states = missing_taxa_list$allowed_states[[i]],
                                 range_candidates = missing_taxa_list$range_candidates[[i]], edge_ages = edge_ages_i)
          phy <- tracked$phy
          tracker <- tracked$tracker
        }
        # Enforce ultrametricity for the new terminal branch only (backbone branches are left untouched)
        phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = focal_taxon_name, root_age = root_age)

        # Save the progress if requested (the phylogeny is saved with its original node labels restored)
        if (save_temp && !is.null(save_temp_every) && (i %% save_temp_every == 0))
        {
          save_temp_progress(path = temp_file_path, phy = phy, labelled_backbone = labelled_backbone, signature = signature, sim_seed = sim_seed,
                             done = i, total = n_missing, unit = "taxa", tracker = tracker)
        }

        # Print progress if requested
        if (verbose && (i %% 100 == 0))
        {
          cat(paste0(Sys.time(), sim_label, " - Missing taxa grafted n\u00b0", i, "/", n_missing, "\n"))
        }
      }
    }

    ## ---- For "BD_clades" strategy: Loop per target clade, nested ones first ---- ##
    if (method == "BD_clades")
    {
      if (verbose) cat(paste0(Sys.time(), sim_label, " - Start identification and ordering of target clades\n"))

      # Identify groups of taxa (target clade, availability of its stem and, with a simmap, allowed states) and order them by clade size
      # Since nested clades are always smaller than their higher parents, this is equivalent to order by nesting scheme
      # Groups of the same clade are consecutive, with the taxa without access to the stem first
      groups <- group_and_order_clades(phy, missing_taxa_list$Target_clade, include_target_clade_stem = missing_taxa_list$include_target_clade_stem,
                                       allowed_states = if (use_states) missing_taxa_list$allowed_states else NULL)
      n_groups <- length(groups$group_clade_tips)

      # Cache of the BD models: the model is fitted once per target clade, before grafting its first group of taxa,
      # and reused for its other groups (same clade, so same rates and sampling fraction)
      BD_cache <- list()

      if (verbose) cat(paste0(Sys.time(), sim_label, " - Start grafting process per clades\n"))

      # Initiate number of taxa and clades grafted
      nb_taxa_grafted <- 0
      nb_clades_grafted <- 0

      # Resume: restart from the phylogeny, the states, the BD models already fitted and the random number generator saved with the last progress (the groups already grafted are skipped)
      if (!is.null(resumed))
      {
        phy <- resumed$phy
        tracker <- resumed$state$tracker
        BD_cache <- resumed$state$bd_cache
        nb_taxa_grafted <- resumed$state$counters$nb_taxa_grafted
        nb_clades_grafted <- resumed$state$counters$nb_clades_grafted
        set_rng_state(resumed$state$rng)
        if (verbose) cat(paste0(Sys.time(), sim_label, " - Resuming from the temporary file: ", nb_clades_grafted, "/", n_groups, " groups of taxa already grafted\n"))
      }

      # Loop per ordered clades
      for (i in groups$group_order[seq_along(groups$group_order) > nb_clades_grafted])
      {
        # Extract ID of taxa in group n°i
        taxa_ID_i <- which(groups$group_ID == i)
        # Extract associated data
        missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[taxa_ID_i] })
        # Extract the tips defining the target clade
        target_clade_tips_i <- groups$group_clade_tips[[i]]
        # Is the stem of the target clade available for the taxa of this group?
        include_target_clade_stem_i <- groups$group_include_target_clade_stem[i]
        # Get number of taxa to graft in this group
        ntips_missing <- length(taxa_ID_i)

        # Identify MRCA node in grafted tree
        # Node is recomputed from tips because MRCA node IDs change after each graft
        MRCA_node_i <- get_clade_node(phy, target_clade_tips_i)
        is_single_tip <- length(target_clade_tips_i) == 1

        # Identify the stem node of the target clade, i.e., its parent node in the BACKBONE phylogeny.
        # Missing taxa previously grafted on the stem created new nodes between the clade and its backbone parent.
        # Walking rootward past those nodes (they only have 'missing-only' subtrees as side branches) gives the true stem node,
        # so the stem age is always the one of the backbone, and not the age of a node created by a graft.
        stem_node_i <- get_stem_region_nodes(phy, MRCA_node_i, missing_tip_labels = missing_taxa_list$Taxon_name)$top_parent

        # Warn against problematic cases with no stem available, because the target clade is the whole tree
        if (is.na(stem_node_i))
        {
          if (include_target_clade_stem_i && !is_single_tip)
          {
            warning("Target clade defined by ", paste(target_clade_tips_i, collapse = ", "),
                    " is the whole tree: it has no stem branch, so only branches inside the clade are used")
          }
          # For case where the target clade is the full tree but no stem edge is required to be included,
          # no need to send a warning; but set the stem node to equal the MRCA node
          stem_node_i <- MRCA_node_i
        }

        # Compute stem age
        stem_age_i <- root_age - phytools::nodeheight(tree = phy, node = stem_node_i)

        ## Deal with the different cases
        if (is_single_tip)
        {
          ## Case for single tips

          # Only its stem branch is available,
          # so only the stem age is provided to TreeSim::corsim()
          divergence_times_i <- stem_age_i
          upper_age_i <- stem_age_i # Use the stem age as upper bound
          ntips_present <- 1

        } else {
          ## Case for clades

          # Extract target clade as crown clade (does not include the stem branch)
          crown_clade_tree_i <- ape::extract.clade(phy = phy, node = MRCA_node_i)
          ntips_present <- ape::Ntip(crown_clade_tree_i)

          # Get branching times within the crown clade
          divergence_times_i <- ape::branching.times(crown_clade_tree_i)
          crown_age_i <- max(divergence_times_i)
          # Use the stem age as upper bound if the stem is available for the taxa of this group
          # Otherwise, use the crown age
          upper_age_i <- if (include_target_clade_stem_i) stem_age_i else crown_age_i
        }

        ## Fit Birth-Death on the clade, or on the smallest enclosing clade with >= min_tips tips
        # Stop at the root if reached
        # Fitted only once per target clade (before its first group is grafted, so all its taxa are still pending), then reused
        clade_key_i <- as.character(groups$group_clade_ID[i])
        if (is.null(BD_cache[[clade_key_i]]))
        {
          fit_i <- fit_bd_on_clade(phy = phy, node = MRCA_node_i, min_tips = min_tips,
                                   richness_info = richness_info) # Sampling fraction accounts for taxa still to graft and unsampled species
          # Record the model for the summary (the taxa of the target clade, whatever their group, use this model)
          if (return_BD_model_summary)
          {
            taxa_clade_i <- which(groups$group_clade_ID[groups$group_ID] == groups$group_clade_ID[i])
            fit_i$record <- build_bd_record(fit = fit_i, phy = phy, initial_tips = target_clade_tips_i, taxa = missing_taxa_list$Taxon_name[taxa_clade_i],
                                            richness_info = richness_info, min_tips = min_tips,
                                            clade_name = get_clade_name_of_taxa(missing_taxa_list$Target_clade_name[taxa_clade_i]))
          }
          BD_cache[[clade_key_i]] <- fit_i
        }
        BD_i <- BD_cache[[clade_key_i]]

        ## Draw the divergence ages with TreeSim::corsim and graft the taxa of the group, from the oldest to the youngest age
        # An age can only be used if an available branch (in an allowed state with a simmap) crosses it in the current tree:
        # the taxa that could not be grafted are left pending and their ages are drawn again until they lend on a valid age (see below)
        pending <- seq_len(ntips_missing) # Taxa of the group still to graft (indices within the group)
        n_redraws <- 0
        repeat
        {
          ## Simulate missing divergence times with TreeSim::corsim (for the pending taxa only)
          new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = BD_i$b, mu = BD_i$d,
                                                    missing = length(pending), told = upper_age_i)
          # Extract the divergence times for the missing taxa
          divergence_times_missing_taxa_i <- extract_new_divergence_times(all_times = new_divergence_times_i, old_times = divergence_times_i)

          # Graft missing taxa from oldest to youngest, so younger taxa can be grafted on branches of older missing taxa
          # (extract_new_divergence_times() sorts the ages in decreasing order)
          not_grafted <- integer(0)
          for (k in seq_along(pending))
          {
            j <- pending[k]
            if (use_excluded_subclades)
            {
              excluded_subclades_tips_list_i <- missing_taxa_list_i$Excluded_subclades[[j]]
            } else {
              excluded_subclades_tips_list_i <- NULL
            }

            grafted_j <- graft_from_divergence_age(
                phy = phy,
                divergence_age = divergence_times_missing_taxa_i[k],
                tip_label = missing_taxa_list_i$Taxon_name[j],
                target_clade_tips = missing_taxa_list_i$Target_clade[[j]],
                excluded_subclades_tips_list = excluded_subclades_tips_list_i,
                exclude_subclades_stem = missing_taxa_list_i$exclude_subclades_stem[[j]],
                include_target_clade_stem = include_target_clade_stem_i,
                missing_tip_labels = missing_taxa_list$Taxon_name,
                root_age = root_age,
                seed = NULL,
                tracker = tracker,
                allowed_states = if (use_states) missing_taxa_list_i$allowed_states[[j]] else NULL,
                range_candidates = if (use_states) missing_taxa_list_i$range_candidates[[j]] else NULL,
                redraw_if_unusable = TRUE)

            # A NULL output means that no available branch (in an allowed state with a simmap) crosses the age: the age is discarded and the taxon stays pending.
            # The ages older than the oldest usable age are the ones that fail: every younger age crosses the lineage grafted at the oldest usable age
            if (is.null(grafted_j))
            {
              not_grafted <- c(not_grafted, j)
            } else if (use_states) {
              phy <- grafted_j$phy
              tracker <- grafted_j$tracker
            } else {
              phy <- grafted_j
            }
          }

          pending <- not_grafted
          if (length(pending) == 0) { break }

          # Draw again the ages of the taxa left. The ages are drawn from the same clade and not conditioned on the lineages grafted so far:
          # the ages of TreeSim::corsim() are independent draws, and the lineages grafted are accounted for by the branches that cross the ages.
          n_redraws <- n_redraws + 1
          if (n_redraws > max_age_redraws)
          {
            # Describe the time window of the divergence ages (from the present to the upper age) and the share of it crossed by available branches in the current tree, for each taxon left
            edge_ages_now <- root_age - phytools::nodeHeights(phy)
            n_show <- min(length(pending), 10)
            coverage_pending <- vapply(pending[seq_len(n_show)], function (j)
            {
              segments_j <- get_available_segments(phy = phy, target_clade_tips = missing_taxa_list_i$Target_clade[[j]],
                                                   excluded_subclades_tips_list = if (use_excluded_subclades) missing_taxa_list_i$Excluded_subclades[[j]] else NULL,
                                                   exclude_subclades_stem = missing_taxa_list_i$exclude_subclades_stem[[j]],
                                                   include_target_clade_stem = include_target_clade_stem_i, missing_tip_labels = missing_taxa_list$Taxon_name,
                                                   tracker = tracker, allowed_states = if (use_states) missing_taxa_list_i$allowed_states[[j]] else NULL, edge_ages = edge_ages_now)
              get_time_coverage(phy = phy, segments = segments_j, edge_ages = edge_ages_now, upper_age = upper_age_i)
            }, numeric(1))
            stop_simulation_failure("No usable divergence age was found for ", length(pending), " taxa of the target clade defined by ", paste(target_clade_tips_i, collapse = ", "),
                 " after ", max_age_redraws, " redraws of their divergence ages ('max_age_redraws').\n",
                 "Time window of the divergence ages: from age 0 to age ", signif(upper_age_i, 4), " (", if (include_target_clade_stem_i) "stem age" else "crown age", " of the target clade).\n",
                 "Share of this window crossed by the available branches", if (use_states) " in the allowed states" else "", " (current tree):\n",
                 paste0("  - ", missing_taxa_list_i$Taxon_name[pending[seq_len(n_show)]], ": ", sprintf("%.2f%%", 100 * coverage_pending),
                        if (use_states) paste0(" (allowed states: ", vapply(missing_taxa_list_i$allowed_states[pending[seq_len(n_show)]], paste, character(1), collapse = ", "), ")") else "", collapse = "\n"),
                 if (length(pending) > n_show) paste0("\n  ... and ", length(pending) - n_show, " more") else "", "\n",
                 "No available branch crosses the ages drawn: the part of the time window where branches are available is too narrow or too far from the ages expected from the diversification model.\n",
                 "Increase 'max_age_redraws', ", if (use_states) "extend the allowed states" else "adjust the excluded subclades or the stem availability", ", or use method = 'uniform'.")
          }
        }

        # Update counts of grafted taxa and grafted clades
        nb_taxa_grafted <- nb_taxa_grafted + ntips_missing
        nb_clades_grafted <- nb_clades_grafted + 1

        # Save the progress if requested (the phylogeny is saved with its original node labels restored)
        if (save_temp && !is.null(save_temp_every) && (nb_clades_grafted %% save_temp_every == 0))
        {
          save_temp_progress(path = temp_file_path, phy = phy, labelled_backbone = labelled_backbone, signature = signature, sim_seed = sim_seed,
                             done = nb_clades_grafted, total = n_groups, unit = "groups", tracker = tracker, bd_cache = BD_cache,
                             counters = list(nb_taxa_grafted = nb_taxa_grafted, nb_clades_grafted = nb_clades_grafted))
        }

        # Print progress if requested
        if (verbose && (nb_clades_grafted %% 10 == 0))
        {
          cat(paste0(Sys.time(), sim_label, " - Missing taxa grafted within clades n\u00b0", nb_clades_grafted, "/", n_groups,
                     " - Taxa n\u00b0", nb_taxa_grafted, "/", n_missing, "\n"))
        }
      }
    }

    # Final check of ultrametricity (can take a long time for large phylogenies)
    is_ultrametric <- ape::is.ultrametric(phy, tol = ultrametric_tol)
    if (!is_ultrametric)
    {
      if (verbose) { cat(paste0(Sys.time(), sim_label, " - Ensure the ultrametricity of the final tree\n") )}
      # Force ultrametry to avoid issue with bad rounding, but only of the grafted tips, to keep the backbone intact
      phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = missing_taxa_list$Taxon_name, root_age = root_age)
    }

    # Restore the original node labels (none if the backbone had none, and an empty string for the nodes created by grafting otherwise)
    phy_tracked <- phy # Save the tracked phylogeny to use to build the updated simmap below
    phy <- restore_node_labels(phy = phy_tracked, labelled = labelled_backbone)

    # With a simmap: the updated simmap of the grafted tree, built from the states followed during grafting (before the grafted elements are added)
    updated_simmap <- NULL
    if (use_states && return_simmap) { updated_simmap <- build_updated_simmap(tracker = tracker, phy_tracked = phy_tracked, phy = phy) }
    # Several traits: the updated composite simmap is split in one simmap per trait (named after the traits)
    if (!is.null(updated_simmap) && !is.null(tracker$traits)) { updated_simmap <- split_trait_simmap(simmap = updated_simmap, traits = tracker$traits) }

    # Add the elements identifying the grafted tips and branches
    # Previously grafted tips are kept if 'phy' was already grafted (the flags accumulate over successive graftings)
    phy <- identify_grafted(phy = phy, grafted_taxa = c(previously_grafted, missing_taxa_list$Taxon_name))
    if (!is.null(updated_simmap)) phy$simmap <- updated_simmap
    # Summary of the birth-death models fitted (one per target clade)
    if (return_BD_model_summary) { attr(phy, "BD_model_summary") <- assemble_bd_summary(lapply(BD_cache, function (x) x$record), type = "taxa") }

    # Save the complete phylogeny, which replaces the saved progress
    if (save_temp) { save_temp_complete(path = temp_file_path, phy = phy, signature = signature) }

    return(phy)
  }

  ## Run the simulations for the given nb of iterations, and forcing rule
  # With a 'multiSimmap': a list with the simulations of each map (see Details)
  phy_all_sims <- run_maps(simulate_one = simulate_one, prepared = prepared_maps, map_names = simmap_input$names, multiple = simmap_input$multiple,
                           nb_simulations = nb_simulations, force_simulations = force_simulations, seed = seed, verbose = verbose, parallel = parallel)

  # Delete all the temporary files of the run: only reached if the run was completed without crashing (a crash keeps the files, so the run can be resumed)
  if (save_temp && clean_temp_files)
  {
    removed_files <- remove_temp_files(temp_file_path = temp_file_path, map_names = simmap_input$names, map_tags = map_tags,
                                       nb_simulations = nb_simulations, force_simulations = force_simulations)
    if (verbose && length(removed_files) > 0) { cat(paste0(Sys.time(), " - ", length(removed_files), " temporary file(s) deleted\n")) }
  }

  # Return output
  return(phy_all_sims)
}


### 2/ Sub-function to draw the grafting position along the available branch segments of a clade, weighted by segment length ####

#' @title Draw the grafting position along the available branch segments of a clade
#'
#' @description Draw where a taxon is grafted, i.e., a branch and a position along it, in two combined steps:
#'   1. The available branches of the clade (target clade, stem region if available, minus excluded subclades) are described as segments
#'      (see `get_branch_segments()`), and one segment is drawn with a probability proportional to its length.
#'   2. The position of the graft is drawn uniformly along the selected segment.
#'
#'   Together, this is a uniform draw over all the available segments, weighted by their length, so every time point along any available segment
#'   has the same probability to be selected.
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param MRCA_node Node ID of the MRCA of the target clade.
#' @param include_target_clade_stem Logical. Whether the stem region of the target clade is available. Ignored for single-tip clades (their tip branch is always available). Default is `TRUE`.
#' @param excluded_nodes Integer vector. Tipward node IDs of the branches removed from the draw (branches of the excluded subclades, and their stem if excluded),
#'   see `get_excluded_branch_nodes()`. Default is `NULL`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not), used to identify the 'missing-only subtrees' of the stem region
#'   and include them in the draw if needed. Default is `character(0)`.
#' @param seed Integer. Optional seed. Default is `NULL` as the internal seed is set by the main grafting function.
#' @param tracker Optional. Output of `start_state_tracking()` (updated by `track_graft()` after each graft). If provided, only the pieces of the available branches
#'   in `allowed_states` are used as segments. Default is `NULL` (whole branches).
#' @param allowed_states Character vector. States in which grafting is allowed. Only used with a `tracker`.
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy`: `root_age - phytools::nodeHeights(phy)`. Only used with a `tracker`.
#'
#' @details The two steps drawing process ensures that this is a uniform draw over all the available segments, weighted by their length, so every time point along any available segment
#'   has the same probability to be selected. If the draw were not weighted, a time point along a short segment would have more chances to be used
#'   for grafting than a time point along a long one.
#'
#'   Both steps are done in one function because the position is drawn within the selected segment, and it is returned relative to the
#'   branch the segment belongs to (distance from the tipward node of the branch), as required to graft the new tip on this branch.
#'   This will matter when a branch is split in several available segments (e.g., when only a part of a branch is available).
#'   Returns `NULL` if no branch is available (all branches excluded).
#'
#' @importFrom ape Ntip
#'
#' @noRd
#'
#' @return `NULL` if no branch is available. Otherwise a list with 2 elements:
#'  * `node` Integer. Tipward node ID of the selected branch.
#'  * `position` Numeric. Position of the graft: distance from the tipward node of the selected branch.
#'


draw_random_grafting_position_along_available_branch_segments <- function (
    phy,
    MRCA_node, # MRCA node of the target clade where to draw the branch
    include_target_clade_stem = TRUE, # Whether to include the stem branch among the available branch for sampling
    excluded_nodes = NULL, # Tipward nodes of the branches (excluded subclades and their stem if excluded) to exclude from the drawing
    missing_tip_labels = character(0), # Labels of all taxa to graft: used to identify the missing-only subtrees of the stem region and include them in the draw if needed.
    seed = NULL, # Seed is set by the main function. Default NULL to not alter the global RNG state
    tracker = NULL, # Output of start_state_tracking(), updated by track_graft(): restricts the draw to the segments in allowed states. NULL: no restriction
    allowed_states = NULL, # States in which grafting is allowed. Only used with a tracker
    edge_ages = NULL) # Ages (rootward, tipward) of the branches of phy: root_age - phytools::nodeHeights(phy). Only used with a tracker
{
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Identify the tipward nodes of all the branches available in the target clade:
  #  - all branches within the clade
  #  - its stem region, if include_target_clade_stem = TRUE (always for single-tip clades, never for the root).
  #    The stem region includes the pieces of the stem split by previous grafts,
  #    and the 'missing-only subtrees' hanging from the stem (i.e., previously grafted taxa and their own grafts),
  #    so taxa grafted on the stem can receive further grafts, like taxa grafted within the crown. See get_stem_region_nodes().
  descendant_nodes <- get_available_branch_nodes(phy = phy, MRCA_node = MRCA_node, include_target_clade_stem = include_target_clade_stem,
                                                 missing_tip_labels = missing_tip_labels)

  # Remove branches of the excluded subclades (computed by get_excluded_branch_nodes()):
  # the branches within the subclades, and their stem region if their stem is excluded.
  # If the stem is not excluded, grafting there makes the new taxon sister to, not nested in, the excluded subclade
  if (length(excluded_nodes) > 0)
  {
    descendant_nodes <- setdiff(descendant_nodes, excluded_nodes)
  }

  # Return NULL if no valid branch remains for drawing
  # (the calling function stops with a message naming the taxon and its target clade)
  if (length(descendant_nodes) == 0)
  {
    return(NULL)
  }

  # Describe the available branches as segments (the pieces in allowed states if there is a simmap),
  # then draw a segment weighted by its length and a position within it
  if (is.null(tracker))
  {
    segments <- get_branch_segments(phy = phy, nodes = descendant_nodes)
  } else {
    segments <- get_allowed_segments_tracked(tracker = tracker, phy = phy, nodes = descendant_nodes, allowed_states = allowed_states, edge_ages = edge_ages)
    if (nrow(segments) == 0) { return(NULL) }
  }
  return(draw_point_in_segments(segments))
}


### 3/ Sub-function grafting missing taxa from divergence age #####


#' @title Graft one taxon at a given divergence age
#'
#' @description Grafts one tip on a random allowed branch of the target clade that crosses `divergence_age`.
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param divergence_age Numeric. Age at which the new branch is grafted.
#' @param tip_label Character string. Label of the new tip.
#' @param target_clade_tips Character vector. One tip, or a pair of tips, defining the target clade.
#' @param excluded_subclades_tips_list List of pairs of tips (or single tips) defining the subclades to exclude.
#' @param exclude_subclades_stem Logical vector with one value per excluded subclade, or a single value recycled. Whether the stem of the excluded subclades is also excluded. Default is `TRUE`.
#' @param include_target_clade_stem Logical. Whether the stem region of the target clade is available for grafting. Ignored for single-tip clades (always available). Default is `TRUE`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not), used to identify the 'missing-only subtrees' of the stem region.
#'   Default is `character(0)`.
#' @param root_age Numeric. Optional root age, to avoid recomputing node heights. Default is `NULL`.
#' @param seed Integer. Optional seed. Default is `NULL` as the internal seed is set by the main grafting function.
#' @param tracker Optional. Output of `start_state_tracking()` (updated by `track_graft()` after each graft). If provided, only the pieces of the available branches
#'   in `allowed_states` that cross the age are used. Default is `NULL` (whole branches).
#' @param allowed_states Character vector. States in which grafting is allowed. Only used with a `tracker`.
#' @param range_candidates `NULL`, or named list (named after the states) of character vectors: the ranges the new taxon can inherit when grafted in each state
#'   (`treat_states_as_ranges = TRUE`). Only used with a `tracker`. Default is `NULL`: the taxon takes the state of the point where it was grafted.
#' @param redraw_if_unusable Logical. If `TRUE`, the function returns `NULL` when no available branch crosses `divergence_age`, so the caller can draw another age.
#'   If `FALSE` (default), it stops with an informative message.
#'
#' @importFrom ape getMRCA
#'
#' @importFrom phytools bind.tip nodeHeights nodeheight
#'
#' @noRd
#'
#' @return Without `tracker`, the phylogeny with the new grafted tip. With a `tracker`, a list with 2 elements: the phylogeny with the new grafted tip (`phy`) and the updated `tracker`.
#'   If no available branch (in an allowed state) crosses the age, the function returns `NULL` if `redraw_if_unusable = TRUE` (the caller draws another age), and stops otherwise.
#'


graft_from_divergence_age <- function (phy, # Backbone phylogeny
                                       divergence_age, # Age at which the branch must be grafted
                                       tip_label, # Label of the tip
                                       target_clade_tips, # 1 tip or a pair of tips defining the target clade
                                       excluded_subclades_tips_list, # list of tip pairs (or single tips) defining subclades to exclude
                                       exclude_subclades_stem = TRUE, # Is the stem of the excluded subclades also excluded? (one value per subclade, or one value)
                                       include_target_clade_stem = TRUE, # Is the stem region of the target clade available?
                                       missing_tip_labels = character(0), # Labels of all taxa to graft: used to identify the missing-only subtrees of the stem region
                                       root_age = NULL, # Pass it to avoid recomputing node heights every time
                                       seed = NULL, # Seed is set by the main function. Default NULL to not alter the global RNG state
                                       tracker = NULL, # Output of start_state_tracking(), updated by track_graft(): restricts the grafting to the segments in allowed states
                                       allowed_states = NULL, # States in which grafting is allowed. Only used with a tracker
                                       range_candidates = NULL, # Ranges that can be inherited by the new taxon in each state (treat_states_as_ranges = TRUE). Only used with a tracker
                                       redraw_if_unusable = FALSE) # Return NULL (instead of stopping) if no available branch crosses the age, so the caller can draw another age
{
  # Set the seed if provided
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Compute root age if not provided
  if (is.null(root_age))
  {
    root_age <- max(phytools::nodeHeights(tree = phy)[, 2])
  }

  # Ages (rootward, tipward) of the branches, used to identify the segments crossing the target age (tipward node younger, rootward node older)
  edge_ages <- root_age - phytools::nodeHeights(tree = phy)

  # Describe the branches available for this taxon (target clade, minus excluded subclades, in allowed states if there is a simmap) as segments,
  # and keep only the segments crossing the target age
  all_segments <- get_available_segments(phy = phy, target_clade_tips = target_clade_tips, excluded_subclades_tips_list = excluded_subclades_tips_list,
                                         exclude_subclades_stem = exclude_subclades_stem, include_target_clade_stem = include_target_clade_stem,
                                         missing_tip_labels = missing_tip_labels, tracker = tracker, allowed_states = allowed_states, edge_ages = edge_ages)
  valid_segments <- get_segments_crossing_age(phy = phy, segments = all_segments, age = divergence_age, root_age = root_age, edge_ages = edge_ages)

  # No valid segment means that this age cannot be used: the caller can draw another one
  if (nrow(valid_segments) == 0 && redraw_if_unusable) { return(NULL) }

  # Send informative error if no valid segment remains
  if (nrow(valid_segments) == 0)
  {
    stop("No branch available at age ", signif(divergence_age, 4), " to graft '", tip_label,"'.\n",
         "All branches of its target clade (tips: ", paste(target_clade_tips, collapse = ", "), ") crossing that age are excluded.\n",
         "Adjust 'include_target_clade_stem' (global or in missing_taxa_list), 'missing_taxa_list$Excluded_subclades' or 'exclude_subclades_stem' (global or in missing_taxa_list) for this taxon.")
  }

  # Sample one branch (based on tipward nodes) among the branches of the valid segments, uniformly:
  # at a given age each lineage crossing it has the same chance to receive the new lineage
  selected_node_ID <- valid_segments$node[sample_one(seq_len(nrow(valid_segments)))]

  # Compute the distance from the tipward node of the selected branch up to the divergence age
  # as the position of the graft on the selected branch
  selected_branch_length <- phy$edge.length[match(selected_node_ID, phy$edge[, 2])]
  position <- divergence_age - (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))
  position <- min(max(0, position), selected_branch_length) # Ensure that the position is within the branch length

  # Graft the missing taxa on the selected branch at the target age
  phy_before_graft <- phy
  phy <- phytools::bind.tip(tree = phy, tip.label = tip_label, edge.length = divergence_age,
                            where = selected_node_ID, position = position)

  # With a simmap: the new taxon takes the state of the point where it was grafted, and the new node is labelled to keep following the branches
  if (!is.null(tracker))
  {
    tracked <- track_graft(tracker = tracker, phy_before = phy_before_graft, phy_after = phy, node = selected_node_ID, position = position,
                           tip_label = tip_label, allowed_states = allowed_states, range_candidates = range_candidates, edge_ages = edge_ages)
    phy <- tracked$phy
    tracker <- tracked$tracker
  }

  # Enforce ultrametricity for the new terminal branch only (backbone branches are left untouched)
  phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = tip_label, root_age = root_age)
  if (is.null(tracker)) { return(phy) }
  list(phy = phy, tracker = tracker)
}



### 4/ Sub-function listing the segments available for a taxon #####

#' @title Segments available to graft a taxon
#'
#' @description Describe the branches available to graft a taxon in the current tree as segments: the branches of its target clade
#'  (plus its stem region if the stem is available), minus the branches of its excluded subclades,
#'  restricted to the pieces in its allowed states if there is a simmap (`tracker`).
#'
#' @inheritParams graft_from_divergence_age
#' @param edge_ages Numeric matrix. Ages (rootward, tipward) of the branches of `phy` in the order of `phy$edge`: `root_age - phytools::nodeHeights(phy)`.
#'
#' @noRd
#'
#' @return Segments, as returned by `get_branch_segments()` (in the order of the branches of the tree). 0 rows if there is none.
#'

get_available_segments <- function (phy, target_clade_tips, excluded_subclades_tips_list, exclude_subclades_stem = TRUE, include_target_clade_stem = TRUE,
                                    missing_tip_labels = character(0), tracker = NULL, allowed_states = NULL, edge_ages)
{
  # Identify tipward nodes of the branches available in the target clade: branches within the clade + its stem region if include_target_clade_stem = TRUE.
  # The stem region contains the pieces of the stem split by previous grafts and the 'missing-only subtrees' hanging from it
  # (taxa previously grafted on the stem, and their own grafts), so taxa grafted on the stem can receive further grafts.
  MRCA_node <- get_clade_node(phy, target_clade_tips)
  stem_clade_nodes <- get_available_branch_nodes(phy = phy, MRCA_node = MRCA_node, include_target_clade_stem = include_target_clade_stem,
                                                 missing_tip_labels = missing_tip_labels)

  # Remove branches of the excluded subclades: the branches within the subclades, and their stem region if their stem is excluded
  # (single tips exclude their terminal branch only if their stem is excluded)
  available_nodes_ID <- setdiff(stem_clade_nodes, get_excluded_branch_nodes(phy = phy, excluded_subclades_tips_list = excluded_subclades_tips_list,
                                                                            exclude_subclades_stem = exclude_subclades_stem,
                                                                            missing_tip_labels = missing_tip_labels))
  # Keep the order of the branches in the tree (edge order), so the order of the candidates does not depend on how they were listed
  available_nodes_ID <- available_nodes_ID[order(match(available_nodes_ID, phy$edge[, 2]))]

  # Describe the available branches as segments (the pieces in allowed states if there is a simmap)
  if (is.null(tracker)) { return(get_branch_segments(phy = phy, nodes = available_nodes_ID)) }
  get_allowed_segments_tracked(tracker = tracker, phy = phy, nodes = available_nodes_ID, allowed_states = allowed_states, edge_ages = edge_ages)
}
