##### Functions to graft missing taxa #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Main function to graft missing taxa ####

#' @title Graft missing taxa on a time-calibrated phylogeny
#'
#' @description This function grafts missing taxa onto an ultrametric, time-calibrated backbone phylogeny.
#'   The position of each taxon is constrained by its taxonomic assignment (the target clade, defined in `Target_clade`),
#'   and optionally by exclusions based on species attributes (the excluded subclades, defined in `Excluded_subclades`).
#'   The divergence times of the grafted branches are either drawn uniformly along the branches (`method = "uniform"`),
#'   or drawn from a birth-death model fitted on the target clade (`method = "BD_clades"`, default).
#'   The backbone is never modified: only new tips are added.
#'
#' @param phy An object of class `"phylo"`. The backbone phylogeny on which taxa are grafted.
#'   It must be rooted and ultrametric (i.e., be time-calibrated).
#' @param missing_taxa_list A list with the information on the taxa to graft. Items:
#'   * `$Taxa_name`: Character vector of the names of the taxa to graft. Names must be unique and absent from `phy$tip.label`.
#'   * `$Target_clade`: List (same length as `$Taxa_name`) of character vectors. For each taxon, the tip labels defining its target clade,
#'     i.e., the clade where the taxon is known to belong (e.g., its genus). A pair of tips defines the clade as their MRCA node.
#'     A single tip defines the "clade" as the terminal branch of this tip.
#'   * `$Excluded_subclades`: (Optional) List (same length as `$Taxa_name`) of lists of character vectors. For each taxon, the subclades of its target clade
#'     where the taxon is known not to belong. Each subclade is defined by a pair of tips. Use `NULL` for taxa with no excluded subclades.
#'     Single-tip subclades exclude nothing.
#' @param include_stem Logical. Whether the stem branch of target clades is available for grafting.
#'   If `TRUE`, a taxon can be grafted on the stem branch, and becomes sister to the target clade.
#'   Taxa already grafted on the stem can receive further grafts (see 'stem region' in Details section).
#'   Ignored for single-tip target clades as their terminal branch is their only branch, and is therefore always available for grafting.
#'   Ignored if the target clade is the whole tree: it has no stem. A warning is issued if `include_stem = TRUE`.
#'   Default is `TRUE`.
#' @param use_excluded_subclades Logical. Whether to avoid grafting within the excluded subclades provided in `missing_taxa_list$Excluded_subclades`. Default is `TRUE`.
#' @param method Character string. Method used to draw the divergence times and the branches used for grafting.
#'   Must be one of `"uniform"` or `"BD_clades"` (default). See Details.
#' @param min_tips Integer. Only used for `method = "BD_clades"`. Minimum number of tips of the clade used to fit the birth-death model.
#'   If the target clade has fewer tips, the smallest enclosing clade with at least `min_tips` tips is used instead. Must be at least 3. Default is `5`.
#' @param total_richness_per_clades (Optional) List with the total number of species of the clades, to account for species that are neither in the backbone
#'   nor in `missing_taxa_list` (unsampled species). Items:
#'   * `$Target_clade`: List of character vectors of tips, defining clades as in `missing_taxa_list$Target_clade`.
#'     It must cover all the unique target clades of `missing_taxa_list`. It can include other clades, e.g., enclosing clades.
#'   * `$total_richness`: Vector of positive integers, same length as `$Target_clade`. Total number of species of each clade,
#'     i.e., backbone species, missing species nested within the clade, and unsampled species. Use `NA` if unknown: the unsampled species
#'     are then assumed to be those of the clades with an informed total richness nested within the clade (none if there are no such clades).
#'
#'   Default is `NULL`: no unsampled species. See 'Sampling fraction and unsampled species' in Details section.
#' @param verbose Logical. Whether to display progress. Messages print every 100 taxa for `method = "uniform"`, and every 10 target clades for `method = "BD_clades"`.
#'   Default is `TRUE`.
#' @param save_temp Logical. Whether to save a temporary version of the grafted phylogeny during the process, to recover from a crash on long runs. Default is `FALSE`.
#' @param save_temp_every Integer. Frequency of the temporary saves: every `save_temp_every` taxa for `method = "uniform"`,
#'   and every `save_temp_every` target clades for `method = "BD_clades"`. Only used if `save_temp = TRUE`. Default is `10`.
#' @param temp_file_path Character string. Path of the `.rds` file used for temporary saves. It must end with `.rds` and its directory must exist and be writable. Default is `NULL`, which saves in `"./temp_phy_from_grafting.rds"`.
#' @param ultrametric_tol Numeric. Tolerance used by `ape::is.ultrametric()` for the ultrametricity checks of the backbone and of the final phylogeny.
#'   Default is `.Machine$double.eps^0.5`.
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
#' The function stops, before any grafting, if the backbone is not ultrametric, or if some taxa have no branch available for grafting
#' (i.e., all branches of their target clade are excluded). The taxa concerned are listed in the error message,
#' so the constraints (`include_stem` or `Excluded_subclades`) can be adjusted.
#' With `method = "BD_clades"`, the function can also stop during grafting if no available branch crosses the divergence age drawn for a taxon.
#'
#' ### Target clades, stem regions and excluded subclades
#'
#' Each taxon is grafted within its target clade. The branches available for grafting are all the branches of the target clade
#' (descending from its MRCA node), plus its stem region if `include_stem = TRUE`. Grafting on the stem places the taxon sister to the target clade.
#'
#' The stem region is not restricted to the original stem branch. A taxon grafted on the stem splits it in two pieces, and becomes itself part of the stem region:
#' later taxa can be grafted on both pieces of the stem, and on the terminal branch of the earlier taxon (the two taxa are then sister, and together sister to the target clade).
#' More generally, the stem region includes all the branches of the 'missing-only subtrees' hanging from the stem, i.e., subtrees containing only taxa of `missing_taxa_list$Taxa_name`,
#' whatever their own target clade. This mirrors what happens within the crown group, where missing taxa can already be grafted on the branches of other missing taxa.
#'
#' Branches inside excluded subclades are removed from the pool of available branches. However, the stem branch of an excluded subclade stays available:
#' the taxon can be sister to the excluded subclade, but never nested within it.
#'
#' Excluded subclades are expected to be nested within their target clade. The function checks their position before grafting:
#' a warning is issued if an excluded subclade is outside the target clade (the exclusion then has no effect).
#' The function stops if an excluded subclade contains the whole target clade, or is identical to it: no branch (or only the stem branch) would be left for grafting.
#' To graft taxa specifically on the stem branch of a clade or of a single tip, use [graft_missing_clades()] instead.
#'
#' ### Enclosing clades
#'
#' Under the `"BD_clades"` method, when a target clade contains fewer tips than `min_tips`, an alternative higher clade is selected to fit the birth-death model.
#' This 'enclosing clade' is identified by walking rootward the tree node by node, starting from the target clade, until a clade meets the `min_tips` condition.
#' This ensures that the diversification dynamics modeled are based on enough tips to be relevant, while being closely associated with the target clade.
#'
#' ### Methods
#'
#' * `"uniform"`: Taxa are grafted one by one, per target clade from the smallest to the largest, so nested clades are processed first (see below).
#'   Within each target clade, taxa are grafted in random order. The branch used for grafting is drawn among the available branches
#'   with probabilities proportional to branch lengths, and the position on the branch is drawn from a uniform distribution.
#'   This is fast, but does not account for diversification dynamics.
#'   Note that under a birth-death process, the probability of a side lineage splitting off a branch at age t follows the exponential function e^(-rt),
#'   with r the net diversification rate (speciation minus extinction).
#'   A uniform draw implicitly assumes a flat density (r = 0), and gives too much weight to old divergence times compared to the expectation from a birth-death process.
#' * `"BD_clades"`: Taxa are grafted per target clade, from the smallest to the largest, so nested clades are processed first (see below).
#'   For each target clade, a birth-death model is fitted with `phytools::fit.bd()` on the clade, or on the smallest enclosing clade with at least `min_tips` tips.
#'   The sampling fraction accounts for the number of missing taxa provided for grafting in this clade.
#'   The divergence times of all taxa missing in the target clade are then jointly drawn with `TreeSim::corsim()`.
#'   Taxa are grafted from the oldest to the youngest divergence time. Each taxon is grafted on a branch drawn uniformly among the available branches that cross its divergence age.
#'   Younger taxa can therefore be grafted on the branches of older missing taxa. This is the recommended method, but it can be time consuming for large phylogenies.
#'
#' With both methods, nested clades are processed first, so the taxa of a larger clade are grafted on a tree that already
#' includes the taxa of the clades nested within it, and can be grafted on their branches.
#' For `"BD_clades"`, it also ensures that the divergence times of nested clades are accounted for when drawing the divergence times of the larger clade.
#' Taxa of the same target clade are grafted in random order, so they are exchangeable.
#' Clades that are not nested do not interact, so the order in which they are processed has no influence.
#'
#' ### Sampling fraction and unsampled species
#'
#' Under `method = "BD_clades"`, the birth-death model of a clade is fitted with the sampling fraction
#' `n_present / (n_present + n_pending + n_unsampled)`, with `n_present` the tips of the fitted clade in the current tree (backbone tips and taxa already grafted),
#' `n_pending` the taxa of `missing_taxa_list` still to graft in this clade, and `n_unsampled` the unsampled species derived from `total_richness_per_clades`.
#' Taxa grafted in nested clades processed earlier are already in the tree, so they are not counted as missing anymore.
#' Unsampled species are only used to fit the rates: they are not grafted, and do not compete with the grafted taxa for the branches.
#'
#' For a clade with a `total_richness` value provided: `nb_unsampled = total_richness - nb_backbone - nb_missing`, with `nb_missing` the number of taxa to graft
#' whose target clade is this clade or is nested within it. For a clade with no `total_richness` value (`NA`, or enclosing a target clade without being listed),
#' `nb_unsampled` is the sum over the largest clades nested within it with a `total_richness` value provided.
#' As such, the unsampled species of nested clades are summed across disjoint clades but never counted twice within nested clades.
#' If `total_richness_per_clades = NULL`, `nb_unsampled = 0`.
#'
#' The function stops, before any grafting, if `total_richness_per_clades` does not cover all the target clades, if a `total_richness` provided is lower than
#' the number of backbone species plus species to graft of the clade, or if a clade has fewer unsampled species than the sum of the clades with
#' a `total_richness` provided nested within it.
#' With `method = "uniform"`, no diversification model is used: `total_richness_per_clades` is only used to derive the `"sampling_fraction"` table (a warning is issued).
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
#' @seealso [graft_missing_clades()] to graft groups of missing taxa as clades sister to a clade of the backbone,
#'   and [check_grafting_constraints()] to check a posteriori that grafted taxa respect their target clades and excluded subclades.
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
#' missing_taxa_list <- list(
#'   Taxa_name = c("X1", "X2", "X3", "X4"),
#'   # Target clades defined by pairs of tips (MRCA), or by a single tip
#'   Target_clade = list(c("t26", "t18"), c("t8", "t22"), "t5", c("t26", "t18")),
#'   # Excluded subclades defined by pairs of tips. Only X1 has an excluded subclade
#'   Excluded_subclades = list(list(c("t10", "t19")), NULL, NULL, NULL))
#'
#' # Graft with the default method ("BD_clades")
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   method = "BD_clades", verbose = FALSE, seed = 1234)
#' plot(grafted_phy, cex = 0.7)
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
#' # Account for unsampled species with the total richness of the clades.
#' # All the unique target clades must be listed
#' # (NA = unknown richness, i.e., no unsampled species beyond nested clades)
#' total_richness_per_clades <- list(
#'   Target_clade = list(c("t26", "t18"), c("t8", "t22"), "t5"),
#'   total_richness = c(25, NA, 3))
#' grafted_phy_richness <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                            total_richness_per_clades = total_richness_per_clades,
#'                                            method = "BD_clades", verbose = FALSE, seed = 1234)
#' # Sampling fraction of the clades in the grafted phylogeny
#' attr(grafted_phy_richness, "sampling_fraction")
#'

graft_missing_taxa <- function (phy, # Backbone phylogeny
                                missing_taxa_list, # List with Taxa_name = c(), Target_clade = list(c(), ...), Excluded_subclades = list(list(c(), c()), NULL, ...)
                                include_stem = TRUE, # Include the stem of the target clades when grafting?
                                use_excluded_subclades = TRUE, # Avoid grafting within Excluded_subclades?
                                method = "BD_clades", # "uniform" or "BD_clades"
                                min_tips = 5, # Minimum number of tips of the clade used to fit the birth-death model for "BD_clades" (the enclosing clade is used if the focal clade is smaller)
                                total_richness_per_clades = NULL, # NULL, or list with Target_clade = list(c(), ...) and total_richness = c(...) (NA allowed). Total number of species of the clades
                                verbose = TRUE, # To print progress every 100 taxa or every 10 target clades
                                save_temp = FALSE, # Save a temporary phylogeny after each target clade (BD_clades)
                                save_temp_every = 10,
                                temp_file_path = NULL, # Default = "./temp_phy_from_grafting.rds"
                                ultrametric_tol = .Machine$double.eps^0.5, # Tolerance used by ape::is.ultrametric() for the ultrametricity checks
                                seed = 1234)
{
  ### ---- Check input validity ---- ###
  {
    ## Check 'method'
    if (!(method %in% c("uniform", "BD_clades")))
    {
      stop("'method' must be either 'uniform' or 'BD_clades'")
    }

    ## The backbone must be ultrametric (error otherwise). It is never modified: making it ultrametric is left to the user
    check_backbone_ultrametric(phy = phy, tol = ultrametric_tol)

    ## Check 'min_tips'
    if (!is.numeric(min_tips) || length(min_tips) != 1 || min_tips < 3)
    {
      stop("'min_tips' must be a single number >= 3")
    }

    ## The backbone tree must have at least 'min_tips' tips, so a BD model can be fit
    if (method == "BD_clades" && length(phy$tip.label) < min_tips)
    {
      stop("Backbone phylogeny too small to fit a birth-death model (< ", min_tips, " tips).\n",
           "Choose method = 'uniform', lower 'min_tips', or provide a larger backbone phylogeny.")
    }

    ## Check length of items in missing_taxa_list
    n_missing <- length(missing_taxa_list$Taxa_name)

    if (n_missing == 0) stop("'missing_taxa_list$Taxa_name' is empty")

    if (length(missing_taxa_list$Target_clade) != n_missing)
      stop("Target_clade must have the same length as Taxa_name")

    if (!is.null(missing_taxa_list$Excluded_subclades) && length(missing_taxa_list$Excluded_subclades) != n_missing)
      stop("Excluded_subclades must have the same length as Taxa_name")

    ## Check that taxa names are truly missing from the backbone
    if (anyDuplicated(missing_taxa_list$Taxa_name) || any(missing_taxa_list$Taxa_name %in% phy$tip.label))
      stop("Taxa_name must be unique and absent from the backbone phylogeny")

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
    #  - Excluded subclade identical to the target clade: the whole crown is excluded. No branch is left if include_stem = FALSE,
    #    and only the stem branch is left if include_stem = TRUE. Grafting only on a stem branch is not what this function is designed for:
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
        paste0(paste0("  - ", df$Taxa_name[seq_len(n_show)], " (target clade tips: ", df$Target_clade[seq_len(n_show)],
                      "; excluded subclade tips: ", df$Excluded_subclade[seq_len(n_show)], ")", collapse = "\n"),
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
             "no branch would be available for grafting with include_stem = FALSE, and only the stem branch with include_stem = TRUE.\n",
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

    ## Richness of the clades (unsampled species)
    # Validates 'total_richness_per_clades' if provided, and derives the number of unsampled species of each node of the backbone.
    # Always computed: it also provides the number of missing taxa to include in the sampling fraction of the BD models
    richness_info <- build_richness_info(phy = phy, missing_taxa_names = missing_taxa_list$Taxa_name,
                                         missing_clades = missing_taxa_list$Target_clade,
                                         mode = "target", total_richness_per_clades = total_richness_per_clades)

    if (!is.null(total_richness_per_clades) && method == "uniform")
    {
      warning("'total_richness_per_clades' is only used to derive the 'sampling_fraction' table when method = 'uniform': ",
              "it is not used during the grafting process, which does not account for the diversification process.",
              call. = FALSE)
    }
  }

  ## Check upfront that every taxon has at least one branch available for grafting
  # Better to stop before grafting than in the middle of a long run
  no_branch_taxa <- find_taxa_without_allowed_branch(phy = phy, missing_taxa_list = missing_taxa_list,
                                                     include_stem = include_stem,
                                                     use_excluded_subclades = use_excluded_subclades)
  if (nrow(no_branch_taxa) > 0)
  {
    n_show <- min(nrow(no_branch_taxa), 10)
    stop(nrow(no_branch_taxa), " taxa have no branch available for grafting: all branches of their target clade are excluded.\n",
         paste0("  - ", no_branch_taxa$Taxa_name[seq_len(n_show)], " (target clade tips: ", no_branch_taxa$Target_clade[seq_len(n_show)], ")", collapse = "\n"),
         if (nrow(no_branch_taxa) > n_show) paste0("\n  ... and ", nrow(no_branch_taxa) - n_show, " more") else "",
         "\nAdjust 'include_stem' or 'missing_taxa_list$Excluded_subclades' for these taxa.")
  }

  ## Set seed if provided, for the duration of this function only (the global RNG state is restored on exit)
  set_local_seed(seed)

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
    groups_uniform <- group_and_order_clades(phy, missing_taxa_list$Target_clade)
    missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[groups_uniform$taxa_order] })

    if (verbose) cat(paste0(Sys.time(), " - Start grafting process\n"))

    for (i in seq_len(n_missing))
    {
      focal_taxon_name <- missing_taxa_list$Taxa_name[i]

      # Get node associated to tips
        # 1 tip => tip node
        # 2 tips => MRCA node
      MRCA_node <- get_clade_node(phy, missing_taxa_list$Target_clade[[i]])

      # Get MRCA nodes of excluded subclades (single-tip or empty subclades exclude nothing)
      excluded_subclades_MRCA_nodes <- NULL
      if (use_excluded_subclades)
      {
        # Extract excluded subclades with at least 2 tips (i.e., true clades)
        excluded_pairs <- Filter(function (x) length(x) >= 2, missing_taxa_list$Excluded_subclades[[i]])
        if (length(excluded_pairs) > 0)
        {
          # Extract MRCA for all excluded subclades
          excluded_subclades_MRCA_nodes <- vapply(excluded_pairs, function (tips) ape::getMRCA(phy, tips), numeric(1))
        }
      }

      # Warn if the target clade is the whole tree: it has no stem branch, so only branches inside the clade can be used
      if (include_stem && length(missing_taxa_list$Target_clade[[i]]) > 1 && is.na(get_parent_node(phy, MRCA_node)))
      {
        warning("Target clade of '", focal_taxon_name, "' defined by ", paste(missing_taxa_list$Target_clade[[i]], collapse = ", "),
                " is the whole tree: it has no stem branch, so only branches inside the clade are used")
      }

      # Select branch for grafting using uniform distribution weighted by branch length to ensure a flat density over any time point along branches
      # (no diversification model)
      selected_node_ID <- draw_random_branch_weighted_per_length(
         phy = phy, MRCA_node = MRCA_node, include_stem = include_stem,
         excluded_subclades_MRCA_nodes = excluded_subclades_MRCA_nodes,
         missing_tip_labels = missing_taxa_list$Taxa_name, seed = NULL)

      # Stop with an informative message if no branch is available (all branches excluded by the user's constraints)
      if (is.na(selected_node_ID))
      {
        stop("No branch available to graft '", focal_taxon_name, "'.\n",
             "All branches of its target clade (tips: ", paste(missing_taxa_list$Target_clade[[i]], collapse = ", "), ") are excluded.\n",
             "Adjust 'include_stem' or 'missing_taxa_list$Excluded_subclades' for this taxon.")
      }

      # Extract branch length
      selected_branch_length <- phy$edge.length[match(selected_node_ID, phy$edge[, 2])]
      # Define position from a uniform distribution
      position <- stats::runif(n = 1, min = 0, max = selected_branch_length)

      # Compute the new branch length to reach present time => make the tree ultrametric
      new_branch_length <- position + (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))
      # Graft new branch
      phy <- phytools::bind.tip(tree = phy, tip.label = focal_taxon_name, edge.length = new_branch_length,
                                where = selected_node_ID, position = position)
      # Enforce ultrametricity for the new terminal branch only (backbone branches are left untouched)
      phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = focal_taxon_name, root_age = root_age)

      # Save temporary phylo if requested
      if (save_temp && (i %% save_temp_every == 0))
      {
        if (is.null(temp_file_path))
        {
          temp_file_path <- "./temp_phy_from_grafting.rds"
        }
        saveRDS(object = phy, file = file.path(temp_file_path))
      }

      # Print progress if requested
      if (verbose && (i %% 100 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted n\u00b0", i, "/", n_missing, "\n"))
      }

    }
  }

  ## ---- For "BD_clades" strategy: Loop per target clade, nested ones first ---- ##
  if (method == "BD_clades")
  {
    if (verbose) cat(paste0(Sys.time(), " - Start identification and ordering of target clades\n"))

    # Identify target clades and order by clade size
    # Since nested clades are always smaller than their higher parents, this is equivalent to order by nesting scheme
    groups <- group_and_order_clades(phy, missing_taxa_list$Target_clade)
    n_groups <- length(groups$group_clade_tips)

    if (verbose) cat(paste0(Sys.time(), " - Start grafting process per clades\n"))

    # Initiate number of taxa and clades grafted
    nb_taxa_grafted <- 0
    nb_clades_grafted <- 0

    # Loop per ordered clades
    for (i in groups$group_order)
    {
      # Extract ID of taxa in group n°i
      taxa_ID_i <- which(groups$group_ID == i)
      # Extract associated data
      missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[taxa_ID_i] })
      # Extract the tips defining the target clade
      target_clade_tips_i <- groups$group_clade_tips[[i]]
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
      stem_node_i <- get_stem_region_nodes(phy, MRCA_node_i, missing_tip_labels = missing_taxa_list$Taxa_name)$top_parent

      # Warn against problematic cases with no stem available, because the target clade is the whole tree
      if (is.na(stem_node_i))
      {
        if (include_stem && !is_single_tip)
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
        # Use the stem age as upper bound if include_stem = TRUE
        # Otherwise, use the crown age
        upper_age_i <- if (include_stem) stem_age_i else crown_age_i
      }

      ## Fit Birth-Death on the clade, or on the smallest enclosing clade with >= min_tips tips
      # Stop at the root if reached
      BD_i <- fit_bd_on_clade(phy = phy, node = MRCA_node_i, n_missing = ntips_missing, min_tips = min_tips,
                             richness_info = richness_info) # Sampling fraction accounts for taxa still to graft and unsampled species

      ## Simulate missing divergence times with TreeSim::corsim
      new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = BD_i$b, mu = BD_i$d,
                                                missing = ntips_missing, told = upper_age_i)
      # Extract the divergence times for the missing taxa
      divergence_times_missing_taxa_i <- extract_new_divergence_times(all_times = new_divergence_times_i, old_times = divergence_times_i)

      # Graft missing taxa from oldest to youngest, so younger taxa can be grafted on branches of older missing taxa
      divergence_times_missing_taxa_i <- sort(divergence_times_missing_taxa_i, decreasing = TRUE)
      for (j in seq_len(ntips_missing))
      {
        if (use_excluded_subclades)
        {
          excluded_subclades_tips_list_i <- missing_taxa_list_i$Excluded_subclades[[j]]
        } else {
          excluded_subclades_tips_list_i <- NULL
        }

        phy <- graft_from_divergence_age(
            phy = phy,
            divergence_age = divergence_times_missing_taxa_i[j],
            tip_label = missing_taxa_list_i$Taxa_name[j],
            target_clade_tips = missing_taxa_list_i$Target_clade[[j]],
            excluded_subclades_tips_list = excluded_subclades_tips_list_i,
            include_stem = include_stem,
            missing_tip_labels = missing_taxa_list$Taxa_name,
            root_age = root_age,
            seed = NULL)
      }

      # Update counts of grafted taxa and grafted clades
      nb_taxa_grafted <- nb_taxa_grafted + ntips_missing
      nb_clades_grafted <- nb_clades_grafted + 1

      # Save temporary phylo if requested
      if (save_temp && (nb_clades_grafted %% save_temp_every == 0))
      {
        if (is.null(temp_file_path))
        {
          temp_file_path <- "./temp_phy_from_grafting.rds"
        }
        saveRDS(object = phy, file = file.path(temp_file_path))
      }

      # Print progress if requested
      if (verbose && (nb_clades_grafted %% 10 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted within clades n\u00b0", nb_clades_grafted, "/", n_groups,
                   " - Taxa n\u00b0", nb_taxa_grafted, "/", n_missing, "\n"))
      }
    }
  }

  # Final check of ultrametricity (can take a long time for large phylogenies)
  is_ultrametric <- ape::is.ultrametric(phy, tol = ultrametric_tol)
  if (!is_ultrametric)
  {
    if (verbose) { cat(paste0(Sys.time(), " - Ensure the ultrametricity of the final tree\n") )}
    # Force ultrametry to avoid issue with bad rounding, but only of the grafted tips, to keep the backbone intact
    phy <- enforce_ultrametry_new_tips(phy = phy, tip_labels = missing_taxa_list$Taxa_name, root_age = root_age)
  }

  # Attach the sampling fraction table of the clades if 'total_richness_per_clades' was provided
  if (!is.null(richness_info$table)) attr(phy, "sampling_fraction") <- richness_info$table

  return(phy)
}


### 2/ Sub-function to draw randomly a branch within a clade weighted by its length ####

#' @title Draw a branch within a clade, weighted by branch length
#'
#' @description The draw is weighted by branch length to ensure the probability to be grafted at any point along any branch in the clade is uniform.
#'   If not weighted, a time point along shorter branch would have more chances to be used for grafting than a time point along a long branch.
#'   Returns the selected node ID, or NA if no branch is available (all branches excluded).
#'
#' @param phy Backbone phylogeny of class `"phylo"`.
#' @param MRCA_node Node ID of the MRCA of the target clade.
#' @param which_node Character string. Return the `"tipward"` or `"rootward"` node of the selected branch. Default is `"tipward"`.
#' @param include_stem Logical. Whether the stem region of the target clade is available. Ignored for single-tip clades (their tip branch is always available). Default is `TRUE`.
#' @param excluded_subclades_MRCA_nodes Integer vector. MRCA nodes of subclades removed from the draw. Default is `NULL`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not), used to identify the 'missing-only subtrees' of the stem region
#'   and include them in the draw if needed. Default is `character(0)`.
#' @param seed Integer. Optional seed. Default is `NULL` as the internal seed is set by the main grafting function.
#'
#' @importFrom ape Ntip
#'
#' @noRd
#'
#' @return The ID of the selected node, or `NA` if no branch is available.
#'


draw_random_branch_weighted_per_length <- function (
    phy,
    MRCA_node, # MRCA node of the target clade where to draw the branch
    which_node = "tipward", # Return the "tipward" or "rootward" node of the selected branch
    include_stem = TRUE, # Whether to include the stem branch among the available branch for sampling
    excluded_subclades_MRCA_nodes = NULL, # Subclades within the target clade to exclude from the drawing
    missing_tip_labels = character(0), # Labels of all taxa to graft: used to identify the missing-only subtrees of the stem region and include them in the draw if needed.
    seed = NULL) # Seed is set by the main function. Default NULL to not alter the global RNG state
{
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Identify the tipward nodes of all the branches available in the target clade:
  #  - all branches within the clade
  #  - its stem region, if include_stem = TRUE (always for single-tip clades, never for the root).
  #    The stem region includes the pieces of the stem split by previous grafts,
  #    and the 'missing-only subtrees' hanging from the stem (i.e., previously grafted taxa and their own grafts),
  #    so taxa grafted on the stem can receive further grafts, like taxa grafted within the crown. See get_stem_region_nodes().
  descendant_nodes <- get_available_branch_nodes(phy = phy, MRCA_node = MRCA_node, include_stem = include_stem,
                                                 missing_tip_labels = missing_tip_labels)

  # Remove branches inside excluded subclades
  # Do this by removing their descending nodes from the draw
  # Their own stem stays available: grafting there makes the new taxon sister to, not nested in, the excluded subclade
  if (length(excluded_subclades_MRCA_nodes) > 0)
  {
    excluded_nodes <- unlist(lapply(excluded_subclades_MRCA_nodes, function (n) get_descendant_nodes(phy, n)))
    descendant_nodes <- setdiff(descendant_nodes, excluded_nodes)
  }

  # Return NA if no valid branch remains for drawing
  # (the calling function stops with a message naming the taxon and its target clade)
  if (length(descendant_nodes) == 0)
  {
    return(NA)
  }

  # Convert node IDs into subtending branch IDs
  branches_ID <- match(descendant_nodes, phy$edge[, 2])
  # Select one branch among the remaining available branches
  selected_branch_ID <- sample_one(branches_ID, prob = phy$edge.length[branches_ID])

  # Convert back branch ID to tipward/rootward node ID
  if (which_node == "tipward")
  {
    selected_node_ID <- phy$edge[selected_branch_ID, 2]
  } else {
    selected_node_ID <- phy$edge[selected_branch_ID, 1]
  }

  # Export selected node ID
  return(selected_node_ID)
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
#' @param excluded_subclades_tips_list List of pairs of tips defining the subclades to exclude.
#' @param include_stem Logical. Whether the stem region of the target clade is available for grafting. Ignored for single-tip clades (always available). Default is `TRUE`.
#' @param missing_tip_labels Character vector. Labels of all the taxa to graft (grafted or not), used to identify the 'missing-only subtrees' of the stem region.
#'   Default is `character(0)`.
#' @param root_age Numeric. Optional root age, to avoid recomputing node heights. Default is `NULL`.
#' @param seed Integer. Optional seed. Default is `NULL` as the internal seed is set by the main grafting function.
#'
#' @importFrom ape getMRCA
#'
#' @importFrom phytools bind.tip nodeHeights nodeheight
#'
#' @noRd
#'
#' @return The phylogeny with the new grafted tip.
#'


graft_from_divergence_age <- function (phy, # Backbone phylogeny
                                       divergence_age, # Age at which the branch must be grafted
                                       tip_label, # Label of the tip
                                       target_clade_tips, # 1 tip or a pair of tips defining the target clade
                                       excluded_subclades_tips_list, # list of tip pairs defining subclades to exclude
                                       include_stem = TRUE, # Is the stem region of the target clade available?
                                       missing_tip_labels = character(0), # Labels of all taxa to graft: used to identify the missing-only subtrees of the stem region
                                       root_age = NULL, # Pass it to avoid recomputing node heights every time
                                       seed = NULL) # Seed is set by the main function. Default NULL to not alter the global RNG state
{
  # Set the seed if provided
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Compute root age if not provided
  if (is.null(root_age))
  {
    root_age <- max(phytools::nodeHeights(tree = phy)[, 2])
  }

  # Identify tipward_node_ID of branches crossing the target age (tipward node younger, rootward node older)
  edge_ages <- root_age - phytools::nodeHeights(tree = phy)
  compatible_tipward_nodes_ID <- phy$edge[(edge_ages[, 1] >= divergence_age) & (edge_ages[, 2] <= divergence_age), 2]

  # Identify tipward nodes of the branches available in the target clade: branches within the clade + its stem region if include_stem = TRUE.
  # The stem region contains the pieces of the stem split by previous grafts and the 'missing-only subtrees' hanging from it
  # (taxa previously grafted on the stem, and their own grafts), so taxa grafted on the stem can receive further grafts.
  # Only branches crossing the target age are kept below, so older pieces of the stem are naturally discarded when grafting younger taxa.
  MRCA_node <- get_clade_node(phy, target_clade_tips)
  stem_clade_nodes <- get_available_branch_nodes(phy = phy, MRCA_node = MRCA_node, include_stem = include_stem,
                                                 missing_tip_labels = missing_tip_labels)

  # Subset valid branches by keeping only the one crossing the target age
  valid_nodes_ID <- intersect(compatible_tipward_nodes_ID, stem_clade_nodes)

  # Remove branches inside excluded subclades (single-tip subclades exclude nothing)
  for (tips in Filter(function (x) length(x) >= 2, excluded_subclades_tips_list))
  {
    valid_nodes_ID <- setdiff(valid_nodes_ID, get_descendant_nodes(phy, ape::getMRCA(phy, tips)))
  }

  # Send informative error if no valid branches remain
  if (length(valid_nodes_ID) == 0)
  {
    stop("No branch available at age ", signif(divergence_age, 4), " to graft '", tip_label,"'.\n",
         "All branches of its target clade (tips: ", paste(target_clade_tips, collapse = ", "), ") crossing that age are excluded.\n",
         "Adjust 'include_stem' or 'missing_taxa_list$Excluded_subclades' for this taxon.")
  }

  # Sample one branch (based on tipward nodes) among the valid branches
  selected_node_ID <- sample_one(valid_nodes_ID)

  # Compute the distance from the tipward node of the selected branch up to the divergence age
  # as the position of the graft on the selected branch
  selected_branch_length <- phy$edge.length[match(selected_node_ID, phy$edge[, 2])]
  position <- divergence_age - (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))
  position <- min(max(0, position), selected_branch_length) # Ensure that the position is within the branch length

  # Graft the missing taxa on the selected branch at the target age
  phy <- phytools::bind.tip(tree = phy, tip.label = tip_label, edge.length = divergence_age,
                            where = selected_node_ID, position = position)

  # Enforce ultrametricity for the new terminal branch only (backbone branches are left untouched)
  enforce_ultrametry_new_tips(phy = phy, tip_labels = tip_label, root_age = root_age)
}
