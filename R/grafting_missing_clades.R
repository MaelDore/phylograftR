##### Function to graft missing clades #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


# See REVIEW_NOTES.md for suggested revisions

# Each group of missing taxa sharing the same sister clade (same MRCA node) is grafted as ONE clade on the stem branch of that sister clade.
# Divergence time drawn with TreeSim::corsim using a BD model fitted on the sister clade (or the smallest enclosing clade with >= min_tips tips).
# Grow a tree with the requested number of taxa and fitted diversification process
# Use the speciation rate of the sister-clade to define the crown age
# Stretch grown tree so its root age matches the crown age
# Graft missing subclade on the stem branch at the given stem age

# Changes:
#   - backbone ultrametricity checked before grafting (error; the backbone is never modified, the user fixes it beforehand)
#   - ultrametricity enforced on the NEW terminal branches only, backbone branch lengths untouched


# Input format
# List: Taxa name = c("Taxa_A", "Taxa_B")
#       MRCA_tips = list(c("Tip_AL", "Tip_AR), c("Tip_BL", "Tip_BR), ...)


grafting_missing_clades <- function (phy, # Backbone phylogeny
                                     missing_taxa_list,
                                     min_tips = 5, # Minimum number of tips of the clade used to fit the birth-death model (the enclosing clade is used if the sister clade is smaller)
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

    if (length(missing_taxa_list$MRCA_tips) != n_missing)
      stop("MRCA_tips must have the same length as Taxa_name")

    ## Check that taxa names are truly missing from the backbone
    if (anyDuplicated(missing_taxa_list$Taxa_name) || any(missing_taxa_list$Taxa_name %in% phy$tip.label))
      stop("Taxa_name must be unique and absent from the backbone phylogeny")

    ## Check that all tips used to identify MRCA_tips are present in the backbone
    all_MRCA_tips <- unique(unlist(missing_taxa_list$MRCA_tips))
    if (!all(all_MRCA_tips %in% phy$tip.label))
    {
      stop("Some tip labels provided to locate MRCA nodes of inclusive clades are missing from the backbone phylogeny\n",
           "Missing tips: ", paste(setdiff(all_MRCA_tips, phy$tip.label), collapse = ", "))
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
  groups <- group_and_order_clades(phy, missing_taxa_list$MRCA_tips)
  n_groups <- length(groups$MRCA_tips)

  # Stop before grafting if a targeted sister clade is the whole tree: it has no stem to graft on
  root_node <- setdiff(phy$edge[, 1], phy$edge[, 2])
  is_whole_tree <- vapply(groups$MRCA_tips, function (tips) get_clade_node(phy, tips) == root_node, logical(1))
  if (any(is_whole_tree))
  {
    stop("Target sister clade defined by ", paste(groups$MRCA_tips[[which(is_whole_tree)[1]]], collapse = ", "),
         " is the whole tree: it has no stem to graft on")
  }

  if (verbose) cat(paste0(Sys.time(), " - Start grafting process per missing clades\n"))

  # Initiate number of taxa and missing clades grafted
  nb_taxa_grafted <- 0
  nb_missing_clades_grafted <- 0

  # Loop per ordered clades
  for (i in groups$order)
  {
    # Extract ID of taxa in group n°i
    taxa_ID_i <- which(groups$group == i)
    # Extract associated data
    missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[taxa_ID_i] })
    # Extract MRCA tips
    MRCA_tips_i <- groups$MRCA_tips[[i]]
    # Get number of taxa to include in this missing clade
    nb_missing_taxa <- length(taxa_ID_i)

    # Identify MRCA node and stem node in grafted tree
    # Node is recomputed from tips because MRCA node IDs change after each graft
    MRCA_node_i <- get_clade_node(phy, MRCA_tips_i)
    stem_node_i <- get_parent_node(phy, MRCA_node_i)
    is_single_tip <- length(MRCA_tips_i) == 1

    # Stop for special case where the targeted sister-clade is the whoel tree
    if (is.na(stem_node_i))
    {
      stop("Target sister clade defined by ", paste(MRCA_tips_i, collapse = ", "), " is the whole tree: it has no stem to graft on")
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
    BD_i <- fit_bd_on_clade(phy = phy, node = MRCA_node_i, n_missing = 1, min_tips = min_tips)

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
      phy <- phytools::bind.tip(tree = phy, tip.label = missing_taxa_list_i$Taxa_name,
                                edge.length = divergence_time_missing_clade_i,
                                where = MRCA_node_i, position = position_i)
    } else {

      ## Case where the missing 'clade' is a true clade (nb_missing_taxa > 1)

      ## Grow a tree with the requested number of missing taxa under the same process as the sister-clade
      missing_clade_i <- TreeSim::sim.bd.taxa(n = nb_missing_taxa, numbsim = 1,
                                              lambda = BD_i$b, mu = BD_i$d,
                                              frac = 1, complete = FALSE)[[1]]
      missing_clade_i$tip.label <- missing_taxa_list_i$Taxa_name

      ## Stretch it so its stem age matches the divergence time
      current_stem_age_i <- max(phytools::nodeHeights(missing_clade_i)) + missing_clade_i$root.edge
      stretch_factor <- divergence_time_missing_clade_i / current_stem_age_i
      missing_clade_i$edge.length <- missing_clade_i$edge.length * stretch_factor
      missing_clade_i$root.edge <- missing_clade_i$root.edge * stretch_factor

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
      saveRDS(object = phy, file = temp_file_path)
    }

    if (verbose && (nb_missing_clades_grafted %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Missing taxa grafted within clades n°", nb_missing_clades_grafted, "/", n_groups,
                 " - Taxa n°", nb_taxa_grafted, "/", n_missing, "\n"))
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

  return(phy)
}
