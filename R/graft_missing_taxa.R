##### Functions to graft missing taxa #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# See REVIEW_NOTES.md for suggested future implementations
# See my Notes on the Word documents


## Methods:
#  "uniform"  : Graft per taxon. Branch used for grafting is drawn with weights proportional to branch length (Branches from excluded subclades are removed from the draw).
#               Position on the branch is drawn from a uniform position (which bias divergence times towards the past compare to expectation from a BD model)
#  "BD_clades": Graft per nested clade. Divergence times of all missing taxa for a given clade are drawn with TreeSim::corsim.
#               Branches used for grafting are drawn uniformly across all branches crossing those ages (Branches from taxon-specific excluded subclades are removed from the draw for each associated taxon).
#               This is the default approach recommended, but it can be time consuming for large phylogenies.


## Notes to include in Roxygen header
#  Why "uniform" is the biased option and "BD_clades" is not:
#  Under a birth-death process, the probability of a side lineage splitting off a branch at age t follows
#  the exponential function e^(-rt) with r the net diversification rate as lambda - mu.
#  A uniform position/branch draw implicitly assumes a flat density (r = 0): it ignores the diversification
#  dynamics and puts too much weight on old positions compared to expectation from a birth-death process.


### Example on a toyset ####
#
# phy <- sim.bdtree(b = 1, d = 0, n = 10)
# plot(phy)
#
# test_phy <- phy
#
# plot(test_phy)
# nodelabels()
# edgelabels()
#
# MRCA_node <- 18
# MRCA_node <- 5
# MRCA_node <- 12
# excluded_subclades_MRCA_nodes <- c(13, 16)
#
# missing_taxa_list <- list(Taxa_name = c("A", "B", "C"),
#                           Target_clade = list(c("s3", "s10"), c("s3"), c("s3", "s8")),
#                           Excluded_subclades = list(list(c("s3", "s8"), c("s9", "s10")), NULL, list(c("s7", "s3"))))
#
#
#
# grafted_phy <- graft_missing_taxa(phy = test_phy, # Backbone phylogeny
#                                      missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), Target_clade as list(c(), c(), ...), Excluded_subclades as list(list(c(), c()), list(c()), ...)
#                                      include_stem = T, # Include or not the stem of the total clade when grafting?
#                                      use_excluded_subclades = T, # Use Excluded_subclades to avoid grafting within specific subclades
#                                      method = "uniform", # Either "uniform" or "BD_clades"
#                                      verbose = T,
#                                      seed = 1234) # Set the seed
#
# plot(grafted_phy)
#
####


### 1/ Main function to graft missing taxa ####

graft_missing_taxa <- function (phy, # Backbone phylogeny
                                missing_taxa_list, # List with Taxa_name = c(), Target_clade = list(c(), ...), Excluded_subclades = list(list(c(), c()), NULL, ...)
                                include_stem = TRUE, # Include the stem of the target clades when grafting?
                                use_excluded_subclades = TRUE, # Avoid grafting within Excluded_subclades?
                                method = "BD_clades", # "uniform" or "BD_clades"
                                min_tips = 5, # Minimum number of tips of the clade used to fit the birth-death model for "BD_clades" (the enclosing clade is used if the focal clade is smaller)
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
  random_order <- sample.int(n_missing)
  missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[random_order] })

  # Root age: constant during grafting (all new tips reach the present) => computed once
  root_age <- max(phytools::nodeHeights(tree = phy)[, 2])

  ## ---- For "uniform" strategy: Loop per taxon ---- ##
  if (method == "uniform")
  {
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
         excluded_subclades_MRCA_nodes = excluded_subclades_MRCA_nodes, seed = NULL)

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
        saveRDS(object = phy, file = temp_file_path)
      }

      # Print progress if requested
      if (verbose && (i %% 100 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted n°", i, "/", n_missing, "\n"))
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
    n_groups <- length(groups$clade_tips)

    if (verbose) cat(paste0(Sys.time(), " - Start grafting process per clades\n"))

    # Initiate number of taxa and clades grafted
    nb_taxa_grafted <- 0
    nb_clades_grafted <- 0

    # Loop per ordered clades
    for (i in groups$order)
    {
      # Extract ID of taxa in group n°i
      taxa_ID_i <- which(groups$group == i)
      # Extract associated data
      missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[taxa_ID_i] })
      # Extract the tips defining the target clade
      target_clade_tips_i <- groups$clade_tips[[i]]
      # Get number of taxa to graft in this group
      ntips_missing <- length(taxa_ID_i)

      # Identify MRCA node in grafted tree
      # Node is recomputed from tips because MRCA node IDs change after each graft
      MRCA_node_i <- get_clade_node(phy, target_clade_tips_i)
      is_single_tip <- length(target_clade_tips_i) == 1

      # Identify the stem node of the target clade
      stem_node_i <- get_parent_node(phy, MRCA_node_i)

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
      BD_i <- fit_bd_on_clade(phy = phy, node = MRCA_node_i, n_missing = ntips_missing, min_tips = min_tips)

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
        saveRDS(object = phy, file = temp_file_path)
      }

      # Print progress if requested
      if (verbose && (nb_clades_grafted %% 10 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted within clades n°", nb_clades_grafted, "/", n_groups,
                   " - Taxa n°", nb_taxa_grafted, "/", n_missing, "\n"))
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

  return(phy)
}


### 2/ Sub-function to draw randomly a branch within a clade weighted by its length ####

# Draw is weighted by branch length to ensure the probability to be grafted at any point along any branch in the clade is uniform
# If not weighted, a time point along shorter branch would have more chances to be used for grafting than a time point along a long branch
# Returns the selected node ID, or NA if no branch is available (all branches excluded)

draw_random_branch_weighted_per_length <- function (
    phy,
    MRCA_node, # MRCA node of the target clade where to draw the branch
    which_node = "tipward", # Return the "tipward" or "rootward" node of the selected branch
    include_stem = TRUE, # Whether to include the stem branch among the available branch for sampling
    excluded_subclades_MRCA_nodes = NULL, # Subclades within the target clade to exclude from the drawing
    seed = NULL) # Seed is set by the main function. Default NULL to not alter the global RNG state
{
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Identify all descending nodes within the target clade
  descendant_nodes <- get_descendant_nodes(phy = phy, node = MRCA_node)

  # Deal with the special case of single-tip clades: the tip branch is the only branch available (it is the stem of the clade),
  # so include_stem is ignored (set to TRUE)
  if (MRCA_node <= ape::Ntip(phy)) { include_stem <- TRUE }

  # Deal with the special case of the root.
  # The root has no stem, so include_stem is ignored (set to FALSE)
  root_node <- setdiff(phy$edge[, 1], phy$edge[, 2])
  if (MRCA_node == root_node) { include_stem <- FALSE }

  # If requested, include the stem branch by including the MRCA node among 'descending nodes'
  if (include_stem) { descendant_nodes <- c(descendant_nodes, MRCA_node) }

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

graft_from_divergence_age <- function (phy, # Backbone phylogeny
                                       divergence_age, # Age at which the branch must be grafted
                                       tip_label, # Label of the tip
                                       target_clade_tips, # 1 tip or a pair of tips defining the target clade
                                       excluded_subclades_tips_list, # list of tip pairs defining subclades to exclude
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

  # Identify tipward nodes of the branches of the grafting clade, stem included
  MRCA_node <- get_clade_node(phy, target_clade_tips)
  stem_clade_nodes <- c(get_descendant_nodes(phy, MRCA_node), MRCA_node)

  # Subset valid branches by keping only the one crossing the target age
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
