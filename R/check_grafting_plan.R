##### Function to check a grafting plan before grafting #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Check the feasibility of a grafting plan ####

#' @title Check the feasibility of a grafting plan before grafting
#'
#' @description Check, before any grafting, whether the taxa (or missing clades) of `missing_taxa_list` can be grafted** on the backbone `phy` with the constraints chosen
#'   (target clades, excluded subclades, stem availability, species attributes), and report what is available for each of them: the number of branches, the share of
#'   the available branches that is in the allowed states (if `simmap` is used as a constraint), the time window of the divergence ages and the share of it where grafting is possible,
#'   and the allowed states found along the valid branches.
#'
#'   The goal is to **help fine-tune the grafting constraints before committing to a long run**:
#'   the same checks are run at the start of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and [`phylograftR::graft_missing_clades()`][graft_missing_clades]
#'   (which stop if a taxon cannot be grafted), but here **the problems are reported for all the taxa at once**, and the table can be explored.
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny: time-calibrated (ultrametric) and binary.
#' @param missing_taxa_list A list, as given to the grafting functions (see [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] and
#'   [`phylograftR::graft_missing_clades()`][graft_missing_clades]). Items: `$Taxon_name`, and `$Target_clade` (taxa, for `graft_missing_taxa()`)
#'   or `$Sister_clade` (missing clades, for `graft_missing_clades()`), with the optional items of the grafting function.
#'   Exactly one of `$Target_clade` and `$Sister_clade` must be provided: it defines the type of the plan.
#'
#'   It can also be **the list of both types of lineages to graft**, `list(taxa = ..., clades = ...)`, as returned by [`phylograftR::build_missing_taxa_list()`][build_missing_taxa_list]
#'   when both types are provided, to be used with [`phylograftR::graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades]. Each part is then planned on the same backbone
#'   and the output is a named list (see 'Value').
#' @param min_tips Numeric. Same as in the grafting functions. With `method = "BD_clades"`, a birth-death model is fitted on the target clade (or the sister clade) of each taxon,
#'   or on the smallest enclosing clade with at least `min_tips` tips. Here, **no model is fitted: the enclosing clade is only reported.** Default is `5`.
#' @param include_target_clade_stem,use_excluded_subclades,exclude_subclades_stem Same as in [`phylograftR::graft_missing_taxa()`][graft_missing_taxa].
#'   Only used with a `$Target_clade` list.
#' @param simmap `NULL`, or stochastic character map(s) of the backbone, as in the grafting functions (a `"simmap"`, a `"multiSimmap"`, or a list of them for several traits).
#'   Used to constrain grafting based on the matching of the `allowed_states` of each taxa with the states found along branches of the stochastic character map(s).
#' @param allowed_states,treat_states_as_ranges Same as in the grafting functions. Only used with a `simmap`.
#' @param min_valid_fraction Numeric between 0 and 1. Taxa (or missing clades) with a share of available branches in allowed states below this value are flagged as `low_valid_fraction`
#'   (a warning in the grafting functions). Default is `0.1`.
#' @param ultrametric_tol Numeric. Tolerance used by `ape::is.ultrametric()` to check for ultrametricity. Default is `.Machine$double.eps^0.5`.
#' @param verbose Logical. Print a summary? Default is `TRUE`.
#'
#' @details
#' ### Taxa (`$Target_clade` list)
#'
#' For each taxon, the pool of branches available for grafting is built from its target clade, and restricted by its excluded subclades (see [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#' With a `simmap`, the pool is then restricted to the pieces of the branches in the allowed states of the taxon.
#'
#' ### Missing clades (`$Sister_clade` list)
#'
#' **Taxa with the same sister clade are grafted together as one missing clade.**
#' The missing clade is grafted on the stem branch of its sister clade (see [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]).
#' With a `simmap` the share of the stem branch in the allowed states is both the valid fraction and the share of the time window available.
#' The taxa of a missing clade must share their allowed states.
#' Results of the evaluation are reported with one row per missing clade.
#'
#' ### Enclosing clade
#'
#' The enclosing clade is the target clade (or sister clade), or the smallest clade enclosing it with at least `min_tips` tips, in the backbone.
#' With `method = "BD_clades"`, it is used to fit the Birth-Death model and estimate divergence ages of the missing taxa to graft in accordance to the fitted diversification process.
#' It is reported so users can tune `min_tips` accordingly to define more relevant enclosing clades.
#' During the grafting, taxa already grafted in a clade count in its number of tips. The chosen enclosing clade used can therefore be the reported clade or a clade nested within it.
#' It can contain more tips than the one reported here, although the reported clade is the oldest one that can be used, and the number of tips reported is its number of backbone tips.
#'
#' **The sampling fractions of the birth-death models are not reported**: they depend on the taxa already grafted in the clade, so they cannot be known before grafting.
#' Use [`phylograftR::extract_BD_model_summary()`][extract_BD_model_summary] to extract **a posteriori** the effective sampling fraction and diversification rates used during the grafting process.
#'
#' ### Several maps
#'
#' With a `"multiSimmap"`, the plan is checked for each map, as attribute-based constraints are different across maps, and the table summarizes values across all maps
#' (i.e., minimum, mean and maximum of the valid fraction and of the time coverage, and the number of maps in which the taxon cannot be grafted).
#' The table per map is stored in the attribute **`"per_map"`**.
#'
#' @importFrom ape Ntip node.depth
#' @importFrom phytools nodeHeights
#'
#' @export
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa based on target clades and excluded subclades.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades of backbone clades.
#'  * [`phylograftR::graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades] to graft missing taxa and missing clades in one run, with a single set of inputs.
#'  * [`phylograftR::check_grafting_constraints()`][check_grafting_constraints] to check a posteriori that grafted taxa respect their constraints, after grafting.
#'  * [`phylograftR::extract_grafting_info()`][extract_grafting_info] to extract a posteriori information on the location of grafted lineages.
#'  * [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract a posteriori  the richness and the sampling fraction of the clades in each grafted phylogeny.
#'  * [`phylograftR::extract_BD_model_summary()`][extract_BD_model_summary] to extract a posteriori information on the birth-death models fitted during grafting.
#'
#' @return A data.frame with one row per taxon (`$Target_clade` list) or per missing clade (`$Sister_clade` list).
#'   With a `list(taxa, clades)` covering both types of missing lineages, the output is a named list with the elements `taxa` and `clades`,
#'   each data.frame containing the output of its type (as described below), and the overarching attribute `"all_feasible"` (`TRUE` if all the missing taxa and missing clades can be grafted).
#'
#'   Columns for missing taxa:
#'   * `Taxon_name` Character string. Name of the taxon.
#'   * `Target_clade` Character string. Tips defining the target clade (first five).
#'   * `target_n_tips` Integer. Number of tips of the target clade in the backbone.
#'   * `n_branches` Integer. Number of branches available for grafting (before the restriction to the allowed states).
#'   * `window_age` Numeric. Upper end (age) of the time window of the divergence ages: the stem age of the target clade, or its crown age if the stem is not available for grafting.
#'   * `enclosing_n_tips`, `enclosing_age` Integer and numeric. Number of tips in the backbone phylogeny, and crown age of the enclosing clade (see Details).
#'
#'   Columns for missing clades:
#'   * `Missing_clade` Integer. ID of the missing clade (taxa with the same sister clade).
#'   * `Taxa` Character string. Taxa of the missing clade, separated by commas.
#'   * `Sister_clade` Character string. Tips defining the sister clade (first five).
#'   * `sister_n_tips` Integer. Number of tips of the sister clade in the backbone.
#'   * `stem_length` Numeric. Length of the stem branch of the sister clade.
#'   * `stem_age` Numeric. Age of the rootward end of the stem branch: the upper end of the divergence ages.
#'   * `enclosing_n_tips`, `enclosing_age` Integer and numeric. Number of tips in the backbone phylogeny, and crown age of the enclosing clade (see Details).
#'
#'   Columns for any type, with a `simmap`:
#'   * `allowed_states` Character string. The allowed states (resolved among the states of the simmap).
#'   * `states_available` Character string. The states found in the available branches.
#'   * `valid_fraction`, `time_coverage` Numeric. Share of the length of the available branches in allowed states, and share of the time window where at least one available branch is in an allowed state.
#'     With several maps: `valid_fraction_min`, `valid_fraction_mean`, `valid_fraction_max`, and the same for `time_coverage`.
#'   * `n_maps_infeasible` Integer. With several maps, the number of maps where the taxon (or missing clade) cannot be grafted.
#'
#'   Last columns:
#'   * `feasible` Logical. Can the taxon be grafted (in all the maps)?
#'   * `cause` Character string. Why it cannot be grafted:
#'     * `"excluded"`: no branch left by the constraints;
#'     * `"excluded_contains_target"`: an excluded subclade contains the target clade;
#'     * `"excluded_identical"`: an excluded subclade is identical to the target clade;
#'     * `"allowed_states"`: no available branch has a piece in an allowed state;
#'     * `"whole_tree"`: the sister clade is the whoel tree, thus no stem branch exists for grafting;
#'     * `"inconsistent_allowed_states"`: the taxa of a missing clade have different allowed states;
#'     * `"low_valid_fraction"`: the share of available branche pieces in allowed states is below `min_valid_fraction`;
#'     * `"none"`: no issue.
#'     With several maps, the causes found across all maps are concatenated, separated by `"|"`.
#'
#'   Attributes: `"type"` (`"taxa"` or `"clades"`), `"min_tips"`, `"min_valid_fraction"`, `"n_maps"`, and `"per_map"` (with several maps: the diagnostic table of each map, with a `map` column).
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Create a toy backbone phylogeny and associated list of missing taxa
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3"),
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5"))
#'
#' # --- Check the grafting plan for taxonomic constraints --- #
#'
#' # Check the grafting plan a priori to detect any issue at low cost
#' plan <- check_grafting_plan(phy = phy, missing_taxa_list = missing_taxa_list, min_tips = 5)
#' plan
#'
#' # --- Check the grafting plan for attribute-based constraints too --- #
#'
#' # Add a simmap to restrict grafting to some states
#' Q <- matrix(0.04, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' diag(Q) <- -0.08
#' simmap <- phytools::sim.history(phy, Q, anc = "A", message = FALSE)
#'
#' # Check the grafting plan a priori
#' plan_states <- check_grafting_plan(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                    simmap = simmap, allowed_states = "B")
#' # Inspect availability of branches in the allowed states
#' plan_states[, c("Taxon_name", "states_available", "valid_fraction", "time_coverage", "feasible")]
#'
#' # --- Check both types of missing lineages at once --- #
#'
#' # With list(taxa, clades), each part is planned on the same backbone: the output is a named list
#' both_list <- list(
#'   taxa = missing_taxa_list,
#'   clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = c("t8", "t22"), M2 = "t5")))
#' plans <- check_grafting_plan(phy = phy, missing_taxa_list = both_list, verbose = FALSE)
#' names(plans)
#' # "all_feasible" is TRUE if both types of missing lineages are feasible
#' attr(plans, "all_feasible")
#'
#'

check_grafting_plan <- function (phy, # Backbone phylogeny
                                 missing_taxa_list, # Same list as used for grafting. Items: Taxon_name, Target_clade (or Sister_clade), and the optional items
                                 min_tips = 5, # Minimum number of tips of the clade used to fit the birth-death models (the enclosing clade is reported)
                                 include_target_clade_stem = TRUE, # Same as in graft_missing_taxa()
                                 use_excluded_subclades = TRUE, # Same as in graft_missing_taxa()
                                 exclude_subclades_stem = TRUE, # Same as in graft_missing_taxa()
                                 simmap = NULL, # NULL, or stochastic character map(s) of the backbone
                                 allowed_states = NULL, # States of the simmap where grafting is allowed. Default for all taxa, overridden by missing_taxa_list$allowed_states
                                 treat_states_as_ranges = FALSE, # Are the states of the simmap biogeographic ranges? (see the grafting functions)
                                 min_valid_fraction = 0.1, # Flag the taxa with a smaller share of valid branches in allowed states
                                 ultrametric_tol = .Machine$double.eps^0.5, # Tolerance used by ape::is.ultrametric() for the ultrametricity check
                                 verbose = TRUE) # Print a summary
{
  ### ---- Check input validity ---- ###

  ## Both types of taxa, as list(taxa, clades): each part is planned on the same backbone
  typed <- split_typed_list(missing_taxa_list, arg = "missing_taxa_list")
  if (typed$typed)
  {
    if (typed$both)
    {
      plans <- lapply(names(typed$lists), function (type)
      {
        if (verbose) { print_type_title(type) }
        check_grafting_plan(phy = phy, missing_taxa_list = typed$lists[[type]], min_tips = min_tips, include_target_clade_stem = include_target_clade_stem,
                            use_excluded_subclades = use_excluded_subclades, exclude_subclades_stem = exclude_subclades_stem, simmap = simmap,
                            allowed_states = allowed_states, treat_states_as_ranges = treat_states_as_ranges, min_valid_fraction = min_valid_fraction,
                            ultrametric_tol = ultrametric_tol, verbose = verbose)
      })
      names(plans) <- names(typed$lists)
      attr(plans, "all_feasible") <- all(vapply(plans, function (x) all(x$feasible), logical(1)))
      return(plans)
    }
    missing_taxa_list <- typed$lists[[1]]
  }

  ## Type of the plan: taxa ('Target_clade', graft_missing_taxa()) or missing clades ('Sister_clade', graft_missing_clades())
  if (!is.list(missing_taxa_list)) { stop("'missing_taxa_list' must be a list") }
  has_target <- !is.null(missing_taxa_list$Target_clade)
  has_sister <- !is.null(missing_taxa_list$Sister_clade)
  if (has_target && has_sister) { stop("'missing_taxa_list' has both 'Target_clade' and 'Sister_clade': use only one, depending on the grafting function (graft_missing_taxa() or graft_missing_clades())") }
  if (!has_target && !has_sister) { stop("'missing_taxa_list' must have a 'Target_clade' item (taxa grafted by graft_missing_taxa()) or a 'Sister_clade' item (missing clades grafted by graft_missing_clades())") }
  is_clades <- has_sister
  clade_item <- if (is_clades) "Sister_clade" else "Target_clade"

  check_backbone_ultrametric(phy = phy, tol = ultrametric_tol)
  check_backbone_binary(phy = phy)
  if (!is.numeric(min_tips) || length(min_tips) != 1 || is.na(min_tips) || min_tips < 3) { stop("'min_tips' must be a single number >= 3") }
  if (!is.numeric(min_valid_fraction) || length(min_valid_fraction) != 1 || is.na(min_valid_fraction) || min_valid_fraction < 0 || min_valid_fraction > 1)
  {
    stop("'min_valid_fraction' must be a single number between 0 and 1")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) { stop("'verbose' must be TRUE or FALSE") }

  ## Structure of the list
  if (is_clades) { validate_missing_taxa_list(missing_taxa_list, clade_item = "Sister_clade", optional_items = "allowed_states") }
  else { validate_missing_taxa_list(missing_taxa_list, clade_item = "Target_clade") }
  n_missing <- length(missing_taxa_list$Taxon_name)

  if (!is_clades)
  {
    # Stem availability and stem exclusion of each taxon (global arguments overridden by the list)
    missing_taxa_list$include_target_clade_stem <- resolve_include_target_clade_stem(
      include_target_clade_stem_taxa = missing_taxa_list$include_target_clade_stem, include_target_clade_stem = include_target_clade_stem, n_taxa = n_missing,
      taxon_names = missing_taxa_list$Taxon_name)
    missing_taxa_list$exclude_subclades_stem <- resolve_exclude_subclades_stem(
      exclude_subclades_stem_taxa = missing_taxa_list$exclude_subclades_stem, exclude_subclades_stem = exclude_subclades_stem,
      excluded_subclades = missing_taxa_list$Excluded_subclades, n_taxa = n_missing, taxon_names = missing_taxa_list$Taxon_name)
  }

  ## Species attributes
  use_states <- !is.null(simmap)
  simmap_input <- check_simmap_input(simmap)
  range_trait <- NULL
  if (is.null(simmap_input$traits))
  {
    if (!is.logical(treat_states_as_ranges) || length(treat_states_as_ranges) != 1 || is.na(treat_states_as_ranges))
    {
      stop("'treat_states_as_ranges' must be TRUE or FALSE (a logical vector with one value per trait is used with the maps of several traits)")
    }
  } else {
    range_trait <- resolve_range_trait(treat_states_as_ranges, trait_names = simmap_input$traits, traits_named = simmap_input$traits_named)
    treat_states_as_ranges <- FALSE
  }
  if (!use_states)
  {
    if (treat_states_as_ranges) { stop("'treat_states_as_ranges = TRUE' requires a 'simmap', whose states are the biogeographic ranges of the branches of the backbone.") }
    states_in_list <- any(vapply(missing_taxa_list$allowed_states, function (x) length(x) > 0 && !all(is.na(x)), logical(1)))
    if (!is.null(allowed_states) || states_in_list)
    {
      stop("'allowed_states' (or 'missing_taxa_list$allowed_states') requires a 'simmap', which gives the states of the branches of the backbone.")
    }
  }

  ## Taxa absent from the backbone, tips of the clades present in the backbone
  if (any(missing_taxa_list$Taxon_name %in% phy$tip.label))
  {
    stop("Taxon_name must be absent from the backbone phylogeny. Present: ", paste(intersect(missing_taxa_list$Taxon_name, phy$tip.label), collapse = ", "))
  }
  clade_tips_used <- unique(unlist(missing_taxa_list[[clade_item]]))
  if (!all(clade_tips_used %in% phy$tip.label))
  {
    stop("Some tip labels provided to locate MRCA nodes of ", if (is_clades) "sister" else "target", " clades are missing from the backbone phylogeny\n",
         "Missing tips: ", paste(setdiff(clade_tips_used, phy$tip.label), collapse = ", "))
  }
  if (!is_clades && !is.null(missing_taxa_list$Excluded_subclades))
  {
    excluded_tips_used <- unique(unlist(missing_taxa_list$Excluded_subclades))
    if (!all(excluded_tips_used %in% phy$tip.label))
    {
      stop("Some tip labels provided to locate nodes of excluded subclades are missing from the backbone phylogeny\n",
           "Missing tips: ", paste(setdiff(excluded_tips_used, phy$tip.label), collapse = ", "))
    }
  }

  ### ---- Information independent of the maps ---- ###

  # Nodes and ages of the backbone
  labelled_backbone <- label_backbone_nodes(phy)
  phy_lab <- labelled_backbone$phy
  n_tips <- ape::Ntip(phy)
  root_node <- n_tips + 1
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  node_age <- numeric(n_tips + phy$Nnode)
  node_age[phy$edge[, 2]] <- root_age - phytools::nodeHeights(phy)[, 2]
  node_age[root_node] <- root_age
  tips_per_node <- ape::node.depth(phy, method = 1)

  # Enclosing clade: the clade, or the smallest clade enclosing it with at least 'min_tips' tips (a tip: start from its parent)
  get_enclosing <- function (node)
  {
    current <- if (node <= n_tips) get_parent_node(phy, node) else node
    while (tips_per_node[current] < min_tips)
    {
      parent <- get_parent_node(phy, current)
      if (is.na(parent)) { break }
      current <- parent
    }
    c(n_tips = tips_per_node[current], age = node_age[current])
  }
  describe_tips <- function (tips) { label_excluded_subclades(list(tips), max_tips = 5) }

  ### ---- Plan of each map ---- ###

  prepared <- lapply(seq_along(simmap_input$maps), function (k)
  {
    with_map_context(plan_prepare_map(phy = phy_lab, missing_taxa_list = missing_taxa_list, simmap_k = simmap_input$maps[[k]], use_states = use_states,
                                      allowed_states = allowed_states, treat_states_as_ranges = treat_states_as_ranges, range_trait = range_trait), map_name = simmap_input$names[k])
  })
  n_maps <- if (use_states) length(prepared) else 0

  if (!is_clades)
  {
    ## Position of the excluded subclades relative to the target clades: problems that stop the grafting functions
    position_cause <- rep(NA_character_, n_missing)
    if (!is.null(missing_taxa_list$Excluded_subclades) && use_excluded_subclades)
    {
      position <- check_excluded_subclades_position(phy = phy, missing_taxa_list = missing_taxa_list)
      for (relation in c("identical", "contains_target")) # The second overrides the first
      {
        position_cause[match(position$Taxon_name[position$relation == relation], missing_taxa_list$Taxon_name)] <- paste0("excluded_", relation)
      }
    }
    target_nodes <- vapply(missing_taxa_list$Target_clade, function (tips) get_clade_node(phy, tips), numeric(1))
    enclosing <- vapply(target_nodes, get_enclosing, numeric(2))

    per_map <- lapply(seq_along(prepared), function (k)
    {
      mtl <- prepared[[k]]$missing_taxa_list
      flagged <- find_taxa_without_allowed_branch(phy = phy_lab, missing_taxa_list = mtl, include_target_clade_stem = mtl$include_target_clade_stem,
                                                  exclude_subclades_stem = mtl$exclude_subclades_stem, use_excluded_subclades = use_excluded_subclades,
                                                  tracker = prepared[[k]]$tracker, allowed_states = mtl$allowed_states)
      cause <- rep("none", n_missing)
      cause[match(flagged$Taxon_name, mtl$Taxon_name)] <- flagged$Cause
      cause[!is.na(position_cause)] <- position_cause[!is.na(position_cause)]
      valid_fraction <- unname(attr(flagged, "valid_fraction"))
      if (is.null(valid_fraction)) { valid_fraction <- rep(NA_real_, n_missing) }
      time_coverage <- unname(attr(flagged, "time_coverage"))
      if (is.null(time_coverage)) { time_coverage <- rep(NA_real_, n_missing) }
      cause[cause == "none" & !is.na(valid_fraction) & valid_fraction < min_valid_fraction] <- "low_valid_fraction"
      states_available <- attr(flagged, "states_available")
      data.frame(
        Taxon_name = mtl$Taxon_name,
        Target_clade = unname(vapply(mtl$Target_clade, describe_tips, character(1))),
        target_n_tips = unname(tips_per_node[target_nodes]),
        n_branches = unname(attr(flagged, "n_branches")),
        window_age = unname(attr(flagged, "window_age")),
        enclosing_n_tips = unname(enclosing["n_tips", ]),
        enclosing_age = unname(enclosing["age", ]),
        allowed_states = if (use_states) vapply(mtl$allowed_states, function (x) paste(sort(unique(x)), collapse = ", "), character(1)) else NA_character_,
        states_available = if (use_states) vapply(states_available, paste, character(1), collapse = ", ") else NA_character_,
        valid_fraction = valid_fraction,
        time_coverage = time_coverage,
        cause = cause, stringsAsFactors = FALSE, row.names = NULL)
    })
    id_columns <- c("Taxon_name", "Target_clade", "target_n_tips", "n_branches", "window_age", "enclosing_n_tips", "enclosing_age")
  } else {
    ## Missing clades: taxa with the same sister clade are grafted together
    groups <- group_and_order_clades(phy, missing_taxa_list$Sister_clade)
    n_groups <- length(groups$group_clade_tips)
    sister_nodes <- vapply(groups$group_clade_tips, function (tips) get_clade_node(phy, tips), numeric(1))
    stem_row <- match(sister_nodes, phy$edge[, 2])
    stem_length <- phy$edge.length[stem_row]
    stem_age <- ifelse(is.na(stem_row), NA_real_, root_age - phytools::nodeHeights(phy)[stem_row, 1])
    enclosing <- vapply(sister_nodes, get_enclosing, numeric(2))
    group_taxa <- vapply(seq_len(n_groups), function (g) paste(missing_taxa_list$Taxon_name[groups$group_ID == g], collapse = ", "), character(1))
    backbone_edge_ages <- root_age - phytools::nodeHeights(phy)

    per_map <- lapply(seq_along(prepared), function (k)
    {
      mtl <- prepared[[k]]$missing_taxa_list
      tracker <- prepared[[k]]$tracker
      cause <- ifelse(sister_nodes == root_node, "whole_tree", "none")
      allowed_group <- rep(NA_character_, n_groups); states_group <- rep(NA_character_, n_groups); valid_fraction <- rep(NA_real_, n_groups)
      if (use_states)
      {
        describe <- function (i) # Allowed states, and the ranges inherited if any: taxa of a missing clade must share both
        {
          description <- paste(sort(unique(mtl$allowed_states[[i]])), collapse = ", ")
          if (!is.null(mtl$range_candidates))
          {
            candidates <- mtl$range_candidates[[i]]
            candidates <- candidates[order(names(candidates), method = "radix")]
            description <- paste0(description, "; ranges inherited: ", paste0(names(candidates), " -> ", vapply(candidates, function (x) paste(sort(x, method = "radix"), collapse = "|"), character(1)), collapse = ", "))
          }
          description
        }
        for (g in seq_len(n_groups))
        {
          members <- which(groups$group_ID == g)
          descriptions <- vapply(members, describe, character(1))
          allowed_g <- mtl$allowed_states[[members[1]]]
          allowed_group[g] <- paste(sort(unique(allowed_g)), collapse = ", ")
          if (cause[g] != "none") { next }
          states_group[g] <- paste(sort(unique(tracker$state_table$state[tracker$state_table$node == sister_nodes[g]])), collapse = ", ")
          if (length(unique(descriptions)) > 1) { cause[g] <- "inconsistent_allowed_states"; next }
          segments <- get_allowed_segments_tracked(tracker = tracker, phy = phy_lab, nodes = sister_nodes[g], allowed_states = allowed_g, edge_ages = backbone_edge_ages)
          valid_fraction[g] <- sum(segments$to - segments$from) / stem_length[g]
          if (valid_fraction[g] == 0) { cause[g] <- "allowed_states" } else if (valid_fraction[g] < min_valid_fraction) { cause[g] <- "low_valid_fraction" }
        }
      }
      data.frame(
        Missing_clade = seq_len(n_groups),
        Taxa = group_taxa,
        Sister_clade = unname(vapply(groups$group_clade_tips, describe_tips, character(1))),
        sister_n_tips = unname(tips_per_node[sister_nodes]),
        stem_length = stem_length,
        stem_age = stem_age,
        enclosing_n_tips = unname(enclosing["n_tips", ]),
        enclosing_age = unname(enclosing["age", ]),
        allowed_states = allowed_group,
        states_available = states_group,
        valid_fraction = valid_fraction,
        time_coverage = valid_fraction, # The stem is the only branch available: the share of its length in allowed states is the share of the time window
        cause = cause, stringsAsFactors = FALSE, row.names = NULL)
    })
    id_columns <- c("Missing_clade", "Taxa", "Sister_clade", "sister_n_tips", "stem_length", "stem_age", "enclosing_n_tips", "enclosing_age")
  }

  ### ---- Output ---- ###

  infeasible_causes <- c("excluded", "excluded_contains_target", "excluded_identical", "allowed_states", "whole_tree", "inconsistent_allowed_states")
  if (n_maps <= 1)
  {
    output <- per_map[[1]]
    if (!use_states) { output$valid_fraction <- NULL; output$time_coverage <- NULL; output$allowed_states <- NULL; output$states_available <- NULL }
    output$feasible <- !(output$cause %in% infeasible_causes)
    output <- output[, c(id_columns, if (use_states) c("allowed_states", "states_available", "valid_fraction", "time_coverage"), "feasible", "cause")]
  } else {
    names_maps <- simmap_input$names
    long <- do.call(rbind, lapply(seq_along(per_map), function (k) cbind(map = names_maps[k], per_map[[k]], feasible = !(per_map[[k]]$cause %in% infeasible_causes), stringsAsFactors = FALSE)))
    rownames(long) <- NULL
    first <- per_map[[1]]
    n_rows <- nrow(first)
    over_maps <- function (column, fun) vapply(seq_len(n_rows), function (i) { v <- vapply(per_map, function (x) x[[column]][i], numeric(1)); if (all(is.na(v))) NA_real_ else fun(v, na.rm = TRUE) }, numeric(1))
    union_of <- function (column) vapply(seq_len(n_rows), function (i) { v <- unique(unlist(strsplit(vapply(per_map, function (x) x[[column]][i], character(1)), ", ", fixed = TRUE))); paste(sort(v[!is.na(v) & nzchar(v)]), collapse = ", ") }, character(1))
    feasible_matrix <- vapply(per_map, function (x) !(x$cause %in% infeasible_causes), logical(n_rows))
    if (n_rows == 1) { feasible_matrix <- matrix(feasible_matrix, nrow = 1) }
    causes <- vapply(seq_len(n_rows), function (i)
    {
      v <- unique(vapply(per_map, function (x) x$cause[i], character(1)))
      if (length(v) > 1) { v <- setdiff(v, "none") }
      paste(v, collapse = "|")
    }, character(1))
    output <- first[, id_columns, drop = FALSE]
    output$allowed_states <- union_of("allowed_states")
    output$states_available <- union_of("states_available")
    for (column in c("valid_fraction", "time_coverage"))
    {
      output[[paste0(column, "_min")]] <- over_maps(column, min)
      output[[paste0(column, "_mean")]] <- over_maps(column, mean)
      output[[paste0(column, "_max")]] <- over_maps(column, max)
    }
    output$n_maps_infeasible <- as.integer(rowSums(!feasible_matrix))
    output$feasible <- output$n_maps_infeasible == 0
    output$cause <- causes
    attr(output, "per_map") <- long
  }
  rownames(output) <- NULL
  attr(output, "type") <- if (is_clades) "clades" else "taxa"
  attr(output, "min_tips") <- min_tips
  attr(output, "min_valid_fraction") <- min_valid_fraction
  attr(output, "n_maps") <- n_maps

  if (verbose) { print_grafting_plan_summary(output, is_clades = is_clades, use_states = use_states, n_maps = n_maps) }
  output
}


#' @title Prepare the states of one map for the check of a grafting plan
#'
#' @description Same preparation as in the grafting functions: validate the simmap, start the tracker of the states, and resolve the allowed states of each taxon among the states of the map.
#'
#' @param phy Backbone phylogeny with its nodes labelled by `label_backbone_nodes()`.
#' @param missing_taxa_list Resolved list of taxa.
#' @param simmap_k A map (`NULL` without simmap), or the maps of several traits (class `"phylograftR_trait_maps"`).
#' @param use_states,allowed_states,treat_states_as_ranges,range_trait See `check_grafting_plan()`. `range_trait` is the position of the trait treated as ranges (several traits).
#'
#' @noRd
#'
#' @return A list with the elements `tracker` (`NULL` without simmap) and `missing_taxa_list` (with `allowed_states` resolved, and `range_candidates` with ranges).
#'

plan_prepare_map <- function (phy, missing_taxa_list, simmap_k, use_states, allowed_states, treat_states_as_ranges, range_trait)
{
  tracker <- NULL
  if (use_states)
  {
    if (inherits(simmap_k, "phylograftR_trait_maps"))
    {
      merged <- merge_trait_simmaps(phy = phy, trait_maps = simmap_k)
      tracker <- start_state_tracking(phy = phy, simmap = merged$simmap)
      tracker$traits <- merged$traits
      if (is.null(range_trait))
      {
        missing_taxa_list$allowed_states <- resolve_trait_allowed_states(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                                         traits = merged$traits, composite_states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
      } else {
        resolved <- resolve_trait_allowed_ranges(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states, traits = merged$traits,
                                                 composite_states = tracker$states, taxon_names = missing_taxa_list$Taxon_name, range_trait = range_trait)
        missing_taxa_list$allowed_states <- resolved$allowed_states
        missing_taxa_list$range_candidates <- resolved$range_candidates
        tracker <- extend_tracker_states(tracker = tracker, extra_states = resolved$inherited_states)
      }
    } else if (treat_states_as_ranges) {
      tracker <- start_state_tracking(phy = phy, simmap = simmap_k)
      resolved <- resolve_allowed_ranges(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                         states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
      missing_taxa_list$allowed_states <- resolved$allowed_states
      missing_taxa_list$range_candidates <- resolved$range_candidates
      tracker <- extend_tracker_states(tracker = tracker, extra_states = resolved$inherited_states)
    } else {
      tracker <- start_state_tracking(phy = phy, simmap = simmap_k)
      missing_taxa_list$allowed_states <- resolve_allowed_states(allowed_states_taxa = missing_taxa_list$allowed_states, allowed_states = allowed_states,
                                                                 states = tracker$states, taxon_names = missing_taxa_list$Taxon_name)
    }
  }
  list(tracker = tracker, missing_taxa_list = missing_taxa_list)
}


#' @title Print the summary of a grafting plan
#'
#' @param plan Output of `check_grafting_plan()`.
#' @param is_clades,use_states Logical.
#' @param n_maps Integer. Number of maps (0 without simmap).
#'
#' @noRd
#'
#' @return `NULL`, invisibly.
#'

print_grafting_plan_summary <- function (plan, is_clades, use_states, n_maps)
{
  unit <- if (is_clades) "missing clades" else "taxa"
  label <- if (is_clades) plan$Taxa else plan$Taxon_name
  min_tips <- attr(plan, "min_tips")
  min_valid_fraction <- attr(plan, "min_valid_fraction")
  causes <- unlist(strsplit(plan$cause, "|", fixed = TRUE))
  count_cause <- function (x) sum(causes == x)
  n_infeasible <- sum(!plan$feasible)
  valid <- if (n_maps > 1) plan$valid_fraction_min else plan$valid_fraction
  coverage <- if (n_maps > 1) plan$time_coverage_min else plan$time_coverage
  target_n <- if (is_clades) plan$sister_n_tips else plan$target_n_tips
  below_min_tips <- sum(target_n < min_tips)

  format_range <- function (x) if (min(x) == max(x)) as.character(min(x)) else paste0(min(x), " to ", max(x))
  cat(paste0("Grafting plan: ", nrow(plan), " ", unit, if (n_maps > 1) paste0(", ", n_maps, " maps") else "", "\n",
             "  - Can be grafted", if (n_maps > 1) " in all the maps" else "", ": ", sum(plan$feasible), "/", nrow(plan), "\n"))
  if (n_infeasible > 0)
  {
    cat(paste0("  - Cannot be grafted (", n_infeasible, "):\n"))
    for (cause in c("excluded", "excluded_contains_target", "excluded_identical", "allowed_states", "whole_tree", "inconsistent_allowed_states"))
    {
      n <- count_cause(cause)
      if (n > 0) { cat(paste0("      ", format(cause, width = 28), n, "\n")) }
    }
    worst <- which(!plan$feasible)
    n_show <- min(length(worst), 10)
    cat(paste0("      e.g., ", paste(label[worst[seq_len(n_show)]], collapse = ", "), if (length(worst) > n_show) paste0(", ... and ", length(worst) - n_show, " more") else "", "\n"))
  }
  if (use_states)
  {
    feasible <- plan$feasible
    cat(paste0("  - Share of the available branches in allowed states (", if (n_maps > 1) "worst map, " else "", "feasible ", unit, "): ",
               if (any(feasible)) paste0("min ", sprintf("%.1f%%", 100 * min(valid[feasible], na.rm = TRUE)), ", median ", sprintf("%.1f%%", 100 * stats::median(valid[feasible], na.rm = TRUE))) else "none", "\n",
               "  - Share of the time window with grafting possible (", if (n_maps > 1) "worst map, " else "", "feasible ", unit, "): ",
               if (any(feasible)) paste0("min ", sprintf("%.1f%%", 100 * min(coverage[feasible], na.rm = TRUE)), ", median ", sprintf("%.1f%%", 100 * stats::median(coverage[feasible], na.rm = TRUE))) else "none", "\n"))
    n_low <- sum(feasible & valid < min_valid_fraction, na.rm = TRUE)
    if (n_low > 0)
    {
      low <- which(feasible & valid < min_valid_fraction)
      n_show <- min(length(low), 10)
      cat(paste0("  - Below 'min_valid_fraction' = ", min_valid_fraction, " (", n_low, "): ", paste(label[low[order(valid[low])][seq_len(n_show)]], collapse = ", "), if (length(low) > n_show) paste0(", ... and ", length(low) - n_show, " more") else "", "\n"))
    }
  }
  cat(paste0("  - ", if (is_clades) "Sister" else "Target", " clades smaller than 'min_tips' = ", min_tips, " (birth-death model fitted on an enclosing clade): ", below_min_tips, "/", nrow(plan),
             if (below_min_tips > 0) paste0("; enclosing clades of ", format_range(plan$enclosing_n_tips[target_n < min_tips]), " tips in the backbone; in practice, the clade used can be nested within them, as previously grafted taxa count towards 'min_tips'") else "", "\n",
             "  => ", if (n_infeasible == 0) "All can be grafted with these constraints\n" else "Some cannot be grafted with these constraints (see the table)\n"))
  invisible(NULL)
}
