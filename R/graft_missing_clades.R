##### Function to graft missing clades #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


#' @title Graft missing clades on a time-calibrated phylogeny
#'
#' @description This function grafts missing taxa onto an ultrametric, time-calibrated backbone phylogeny as clades sister to clades of the backbone.
#'   Missing taxa sharing the same sister clade are grafted together as ONE missing clade, on the stem branch of their sister clade.
#'   The divergence time of the missing clade is drawn from a birth-death model fitted on the sister clade,
#'   and the missing clade is grown under the same model, conditioned on its number of taxa and its stem age.
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which taxa are grafted.
#'   It must be rooted and ultrametric (i.e., be time-calibrated).
#' @param missing_taxa_list A list with the information on the taxa to graft. Items:
#'   * `$Taxa_name`: Character vector of the names of the taxa to graft. Names must be unique and absent from `phy$tip.label`.
#'   * `$Sister_clade`: List (same length as `$Taxa_name`) of character vectors. For each taxon, the tip labels defining the clade of the backbone
#'     that is sister to the missing clade. A pair of tips defines the clade as their MRCA node. A single tip defines the "clade" as the terminal branch of this tip.
#'     Taxa with the same sister clade are grafted together as one missing clade.
#' @param min_tips Integer. Minimum number of tips of the clade used to fit the birth-death model.
#'   If the sister clade has fewer tips, the smallest enclosing clade with at least `min_tips` tips is used instead. Must be at least 3. Default is `5`.
#' @param total_richness_per_clades (Optional) List with the total number of species of the clades, to account for species that are neither in the backbone
#'   nor in `missing_taxa_list` (i.e., unsampled species) when fitting the birth-death models. Items:
#'   * `$Sister_clade`: List of character vectors of tips, defining clades as in `missing_taxa_list$Sister_clade`.
#'     It must cover all the unique sister clades of `missing_taxa_list`. It can include other clades, e.g., enclosing clades.
#'   * `$total_richness`: Vector of positive integers, same length as `$Sister_clade`. Total number of species of each clade,
#'     i.e., backbone species, missing species nested within the clade, and unsampled species. Use `NA` if unknown: the unsampled species
#'     are then assumed to be those of the clades with an informed total richness nested within the clade (none if there are no such clades).
#'
#'   Missing taxa grafted on the stem of a clade are sister to this clade: they are not counted in its richness.
#'   Default is `NULL`: no unsampled species. Missing species to graft and backbone species are assumed to cover all species in the focal group.
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
#'   If `total_richness_per_clades` is provided, the tree has an attribute `"sampling_fraction"`: a data frame with one row per clade of `total_richness_per_clades`
#'   (columns `Clade_tips`, `total_richness`, `total_richness_provided`, `nb_backbone_species`, `nb_missing_species`, `nb_unsampled_species`
#'   and `sampling_fraction`). `sampling_fraction` is the fraction of the clade richness present in the grafted phylogeny:
#'   `(nb_backbone_species + nb_missing_species) / total_richness`. `total_richness` is derived when `NA` was provided.
#'   Note that functions such as `ape::drop.tip()` do not keep this attribute.
#'
#' @details
#' The function stops, before any grafting, if the backbone is not ultrametric, or if a sister clade is the whole tree (it has no stem branch to graft on).
#'
#' ### Grafting process
#'
#' Groups of missing taxa sharing the same sister clade are ordered by the size of the sister clade, so nested sister clades are processed first.
#' For each group:
#' 1. A birth-death model is fitted with `phytools::fit.bd()` on the sister clade, or on the smallest enclosing clade with at least `min_tips` tips.
#'    The sampling fraction is computed as if only one lineage was missing (the branch leading to the missing clade).
#' 2. The divergence time of the missing clade is drawn on the stem branch of the sister clade with `TreeSim::corsim()`:
#'    between the crown age of the sister clade and its stem age (or between the present and the stem age for a single-tip sister clade).
#' 3. A single missing taxon is grafted as a new tip at this divergence time. For several taxa, a missing clade is grown with `TreeSim::sim.bd.taxa.age()`
#'    under the fitted model, conditioned on the number of taxa and on its stem age (the drawn divergence time), and grafted on the stem branch of the sister clade.
#'
#' ### Enclosing clades
#'
#' When a sister clade contains fewer tips than `min_tips`, an alternative higher clade is selected to fit the birth-death model.
#' This 'enclosing clade' is identified by walking rootward the tree node by node, starting from the sister clade, until a clade meets the `min_tips` condition.
#' This ensures that the diversification dynamics modeled are based on enough tips to be relevant, while being closely associated with the sister clade.
#'
#' ### Sampling fraction and unsampled species
#'
#' Without `total_richness_per_clades`, the sampling fraction of the birth-death model is computed as if only one lineage was missing (the missing-only clade to be grafted on the stem branch).
#' With `total_richness_per_clades`, it is `n_present / (n_present + n_pending + n_unsampled)`, with `n_present` the tips of the fitted clade in the current tree
#' (backbone tips and taxa already grafted, so missing clades grafted in nested sister clades processed earlier are not missing anymore),
#' `n_pending` the taxa of `missing_taxa_list` still to graft whose sister clade is strictly nested within the fitted clade,
#' and `n_unsampled` the unsampled species derived from `total_richness_per_clades`.
#' Taxa that are sister to the fitted clade itself are outside of it, so they are not counted.
#' The divergence time of each missing clade is still drawn with `TreeSim::corsim(missing = 1)`: unsampled species only influence the fitted rates,
#' they do not take a branch.
#'
#' For a clade with a `total_richness` value provided: `nb_unsampled = total_richness - nb_backbone - nb_missing`, with `nb_missing` the number of taxa to graft
#' whose sister clade is strictly nested within the clade. For a clade with no `total_richness` value (`NA`, or enclosing a sister clade without being listed),
#' `nb_unsampled` is the sum of `nb_unsampled` over the largest clades nested within it with a `total_richness` value provided. If there is none, `nb_unsampled = 0`.
#' The rationale is that even if `nb_unsampled` cannot be derived from its `total_richness` (as it is not provided for this clade), it has to at least account for the `nb_unsampled`
#' species that have being identified in its own nested clades.
#'
#' The function stops, before any grafting, if `total_richness_per_clades` does not cover all the sister clades, if a `total_richness` value is lower than
#' the number of backbone species plus missing species of the clade, or if a clade has fewer unsampled species than the sum of the clades with
#' an informed `total_richness` nested within it.
#'
#' @importFrom ape bind.tree branching.times extract.clade is.ultrametric
#'
#' @importFrom phytools bind.tip nodeHeights nodeheight
#'
#' @importFrom TreeSim corsim sim.bd.taxa.age
#'
#' @export
#'
#' @seealso [graft_missing_taxa()] to graft taxa within a target clade with constraints (target clades and excluded subclades),
#'   and [check_grafting_constraints()] to check a posteriori that grafted taxa respect their constraints.
#'
#' @author Maël Doré
#'
#' @examples
#' # Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' # List the missing taxa to graft
#' missing_clades_list <- list(
#'   Taxa_name = c("M1", "M2", "M3", "M4"),
#'   # M1 and M2 share the same sister clade: they are grafted together as one missing clade
#'   Sister_clade = list(c("t8", "t22"), c("t8", "t22"), c("t26", "t18"), "t5"))
#'
#' # Graft the missing clades
#' grafted_phy <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
#'                                     verbose = FALSE, seed = 1234)
#' plot(grafted_phy, cex = 0.7)
#'
#' # Account for unsampled species with the total richness of the sister clades.
#' # All the unique sister clades must be listed (NA = unknown richness)
#' total_richness_per_clades <- list(
#'   Sister_clade = list(c("t8", "t22"), c("t26", "t18"), "t5"),
#'   total_richness = c(20, NA, 3))
#' grafted_phy_richness <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
#'                                              total_richness_per_clades = total_richness_per_clades,
#'                                              verbose = FALSE, seed = 1234)
#' attr(grafted_phy_richness, "sampling_fraction")
#'


graft_missing_clades <- function (phy, # Backbone phylogeny
                                  missing_taxa_list,
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

    ## Check length of items in missing_taxa_list
    n_missing <- length(missing_taxa_list$Taxa_name)

    if (n_missing == 0) stop("'missing_taxa_list$Taxa_name' is empty")

    if (length(missing_taxa_list$Sister_clade) != n_missing)
      stop("Sister_clade must have the same length as Taxa_name")

    ## Check that taxa names are truly missing from the backbone
    if (anyDuplicated(missing_taxa_list$Taxa_name) || any(missing_taxa_list$Taxa_name %in% phy$tip.label))
      stop("Taxa_name must be unique and absent from the backbone phylogeny")

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

    ## Richness of the clades (unsampled species), only if 'total_richness_per_clades' is provided
    # Validates 'total_richness_per_clades', and derives the number of unsampled species of each node of the backbone.
    # Without 'total_richness_per_clades', the sampling fraction keeps the default 'as if one lineage was missing'
    richness_info <- NULL
    if (!is.null(total_richness_per_clades))
    {
      richness_info <- build_richness_info(phy = phy, missing_taxa_names = missing_taxa_list$Taxa_name,
                                           missing_clades = missing_taxa_list$Sister_clade,
                                           mode = "sister", total_richness_per_clades = total_richness_per_clades)
    }
  }

  ## Set seed if provided, for the duration of this function only (the global RNG state is restored on exit)
  set_local_seed(seed)

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
    # Sampling fraction "as if missing only one lineage" (the branch leading to the missing clade)
    # With 'total_richness_per_clades': missing taxa still to graft and unsampled species of the fitted clade are used instead of the single missing lineage
    BD_i <- fit_bd_on_clade(phy = phy, node = MRCA_node_i, n_missing = 1, min_tips = min_tips, richness_info = richness_info)

    ## Simulate missing divergence times with TreeSim::corsim
    new_divergence_times_i <- TreeSim::corsim(
        x = divergence_times_i, lambda = BD_i$b, mu = BD_i$d,
        missing = 1, # The branch leading to the missing clade
        told = upper_age_i, # The stem age
        tyoung = lower_age_i) # The crown age (or 0 for tips)
    # Extract the divergence times for the missing clade
    divergence_time_missing_clade_i <- extract_new_divergence_times(all_times = new_divergence_times_i, old_times = divergence_times_i)

    # Compute crown age of the sister-clade Age of the sister clade's MRCA node and
    MRCA_age_i <- root_age - phytools::nodeheight(tree = phy, node = MRCA_node_i)
    # Extract stem branch length
    stem_branch_length_i <- phy$edge.length[match(MRCA_node_i, phy$edge[, 2])]

    # Compute position of the graft on the stem branch
    position_i <- min(max(0, divergence_time_missing_clade_i - MRCA_age_i), stem_branch_length_i)

    ## Graft missing clades
    if (nb_missing_taxa == 1)
    {
      ## Case where the missing 'clade' is a single tip

      # Graft the unique tip
      phy <- phytools::bind.tip(tree = phy, tip.label = missing_taxa_list_i$Taxa_name,
                                edge.length = divergence_time_missing_clade_i,
                                where = MRCA_node_i, position = position_i)
    } else {

      ## Case where the missing 'clade' is a true clade (nb_missing_taxa > 1)

      ## Grow a tree with the requested number of missing taxa under the same process as the sister-clade,
      # conditioned on BOTH the number of taxa and the stem age (= the divergence time drawn with TreeSim::corsim)
      missing_clade_i <- TreeSim::sim.bd.taxa.age(n = nb_missing_taxa, numbsim = 1,
                                                  lambda = BD_i$b, mu = BD_i$d,
                                                  frac = 1, age = divergence_time_missing_clade_i,
                                                  mrca = FALSE)[[1]] # mrca = FALSE: 'age' is the stem age (time since origin)
      missing_clade_i$tip.label <- missing_taxa_list_i$Taxa_name

      # Graft the generated subclade
      phy <- ape::bind.tree(x = phy, y = missing_clade_i, where = MRCA_node_i, position = position_i)
    }

    # Enforce ultrametricity for the new terminal branches only (backbone branches are left untouched)
    phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = missing_taxa_list_i$Taxa_name, root_age = root_age)

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
      saveRDS(object = phy, file = file.path(temp_file_path))
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
    phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = missing_taxa_list$Taxa_name, root_age = root_age)
  }

  # Attach the sampling fraction table of the clades if 'total_richness_per_clades' was provided
  if (!is.null(richness_info$table)) attr(phy, "sampling_fraction") <- richness_info$table

  return(phy)
}
