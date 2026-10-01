##### Function to graft missing clades #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


# Fit BD model on sister-clade.
# Select position on the stem branch drawn from exponential distribution according to sister-clade speciation rate (not a uniform prior!)
# Check corsim to see if possible to ask for an age comprise between stem age and crown age
# Position also gives stem age of the missing clade
# Use it to grow a tree with the requested number of missing taxa
# Use the speciation rate of the sister-clade to define the crown age
# Stretch grown tree so its root age matches the crown age
# Graft missing subclade on the stem branch at the given stem age


## Currently, does not work to graft on single branch that are less than 2 nodes away from the root, as their is no substitution clade available to compute the BD

# Input format
# List: Taxa name c("Taxa_1", "Taxa_2")
#       MRCA_tips list(c("Tip_1", "Tip_2), c("Tip_1", "Tip_2), ...)



grafting_missing_clades <- function (phy, # Backbone phylogeny
                                     missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c("T1", "T2"), MRCA_tips as list(T1 = c(Tip_1 = "", Tip_2 = ""), T2 = c(Tip_1 = "", Tip_2 = ""), ...).
                                     force_ultrametry = FALSE, # Whether to add an extra step to ensure the final tree is ultrametric (can be time-consumming for large trees)
                                     verbose = T, # To print progress every 100 taxa
                                     save_temp = FALSE, # To save a temporary phylogeny every 10 stem clades
                                     seed = 1234) # Set the seed
{
  # Check that tips for locating MRCA are all found on the phylogeny
  all_MRCA_tips <- unique(unlist(missing_taxa_list$MRCA_tips))
  if(!all(all_MRCA_tips %in% phy$tip.label))
  {
    non_matching_tips <- all_MRCA_tips[!(all_MRCA_tips %in% phy$tip.label)]
    stop(paste0("Some tip labels provided to locate MRCA nodes of grafting clades are missing from the phylogeny\n",
                "Missing tips: ", paste(non_matching_tips, collapse = ", ")))
  }

  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Randomize the order of taxa in the missing_taxa_list so it has no influence on the pattern of grafting
  random_order <- sample(x = 1:length(missing_taxa_list$Taxa_name), size = length(missing_taxa_list$Taxa_name), replace = F)
  missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[random_order] })

  if (verbose)
  {
    cat(paste0(Sys.time()," - Start identification of sister-clades\n"))
  }

  ## Extract root age of the full tree
  full_tree_age <- max(phytools::nodeHeights(phy)[,2])

  ## Get info per missing clades

  # Identify all missing clades used for grafting
  MRCA_tips_merged <- unlist(lapply(missing_taxa_list$MRCA_tips, FUN = paste0, collapse = ","))
  # Group missing taxa per sister-clades
  MRCA_tips_groups <- as.factor(MRCA_tips_merged)
  as.numeric(MRCA_tips_groups)

  # Get info per sister-clades
  groups_ID <- 1:nlevels(MRCA_tips_groups)
  group_indices <- match(x = groups_ID, table = as.numeric(MRCA_tips_groups))

  sister_clades_list <- list(ID = groups_ID,
                             MRCA_tips = missing_taxa_list$MRCA_tips[group_indices])
  names(sister_clades_list$MRCA_tips) <- groups_ID

  # Identify MRCA nodes & descending nodes
  sister_clades_list$MRCA_nodes <- list()
  sister_clades_list$descendant_nodes <- list()
  for (i in 1:length(groups_ID))
  {
    if (length(sister_clades_list$MRCA_tips[[i]]) > 1)
    {
      # Case with pairs of tips
      sister_clades_list$MRCA_nodes[[i]] <- ape::getMRCA(phy = phy, tip = sister_clades_list$MRCA_tips[[i]])
      sister_clades_list$descendant_nodes[[i]] <- phytools::getDescendants(tree = phy, node = sister_clades_list$MRCA_nodes[[i]])
    } else {
      # Case with a single tip => MRCA node = tip node ; no descendants
      sister_clades_list$MRCA_nodes[[i]] <- which(phy$tip.label == sister_clades_list$MRCA_tips[[i]])
      sister_clades_list$descendant_nodes[i] <- list(NULL)
    }

  }
  names(sister_clades_list$MRCA_nodes) <- sister_clades_list$ID
  names(sister_clades_list$descendant_nodes) <- sister_clades_list$ID

  ## Order sister-groups such as the nested ones are selected fist

  if (verbose)
  {
    cat(paste0(Sys.time()," - Start ordering of sister-clades per nestedness\n"))
  }

  # Initiate vector of sister-clades order
  sister_clades_order <- c()
  sister_clades_list_temp <- sister_clades_list
  while (length(sister_clades_order) < length(groups_ID))
  {
    # Initiate stopping criterion as check for nestedness of the sister-clade
    # Selected sister-clade must have no other sister-clade nested within
    no_nested_subclades_within <- F
    j <- 0 # Initiate counter
    while (!no_nested_subclades_within)
    {
      j <- j + 1
      target_sister_clade_descendant_nodes <- unlist(sister_clades_list_temp$descendant_nodes[j])
      target_MCRA_other_sister_clades <- unlist(sister_clades_list_temp$MRCA_nodes)[(j+1):length(sister_clades_list_temp$MRCA_nodes)]
      no_nested_subclades_within <- !any(target_MCRA_other_sister_clades %in% target_sister_clade_descendant_nodes)
    }
    # Extract selected sister-clade ID
    selected_clade_ID <- sister_clades_list_temp$ID[j]
    # Assign sister-clade in the order
    sister_clades_order <- c(sister_clades_order, selected_clade_ID)
    # Remove it from selection
    sister_clades_list_temp <- lapply(X = sister_clades_list_temp, FUN = function (x) { x[-j] })

    # if (verbose & length(sister_clades_order) %% 100 == 0)
    # {
    #   cat(paste0(Sys.time()," - ",length(sister_clades_order), " sister-clades out of ",nlevels(MRCA_tips_groups)," ordered by nestedness\n"))
    # }

    # Stop when all sister-clades are ordered
  }

  if (verbose)
  {
    cat(paste0(Sys.time()," - All sister-clades ordered by nestedness\n"))
  }

  ## Loop per missing clades used for grafting, starting with the nested ones

  if (verbose)
  {
    cat(paste0(Sys.time()," - Start grafting process per missing clades\n"))
  }

  # Initiate number of taxa / missing clades grafted
  nb_taxa_grafted <- 0
  nb_missing_clades_grafted <- 0

  for (i in sister_clades_order)
  {
    # i <- sister_clades_order[1]

    ## Subset missing_taxa_list for the given sister-clade
    missing_taxa_ID_sister_clade_i <- which(as.numeric(MRCA_tips_groups) == i)
    missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[missing_taxa_ID_sister_clade_i] })

    # Extract MRCA tips
    MRCA_tips_i <- sister_clades_list$MRCA_tips[[i]]

    ## Draw divergence time between missing clade and sister-clade based on BD model fit in the sister-clade
    if (length(MRCA_tips_i) > 1)
    {
      ## Case with pairs of tips used to define the sister-clade

      # Extract MRCA node using tips (cannot use previous extract as the tree as been updated with the missing taxa from the previous missing group)
      MRCA_node_i <- ape::getMRCA(phy = phy, tip = MRCA_tips_i)

      # Extract the MRCA crown clade tree = the sister-clade
      sister_clade_tree_i <- ape::extract.clade(phy = phy, node = MRCA_node_i)

      # plot(sister_clade_tree_i)

      # Extract divergence times in the sister-clade tree
      crown_age_i <- max(phytools::nodeHeights(sister_clade_tree_i)[,2])
      divergence_times_i <- crown_age_i - phytools::nodeHeights(sister_clade_tree_i)[,1]

      # Set lower bound to draw divergence age as the crown age (because we want to draw a divergence age on the stem branch only)
      lower_age_i <- crown_age_i

      # Add stem branch by using the stem age as the upper limit to draw the divergence time
      stem_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i),1]
      stem_age_i <- full_tree_age - phytools::nodeheight(tree = phy, node = stem_node_i)
      upper_age_i <- stem_age_i # Set upper bound to drawn divergence age as the stem age, thus include the possibility of branching on the stem branch

      # Compute sampling fraction as if missing only one lineage (the branch leading to the missing clade)
      ntips_present <- length(sister_clade_tree_i$tip.label)
      ntips_missing <- 1
      rho <- ntips_present / (ntips_present + ntips_missing)

      ## Fit BD on crown clade, or substitute extended crown clade according to number of tips in the crown clade
      if (ntips_present > 2)
      {
        ## Case with more than 2 lineages in the crown clade: fit a BD on the crown clade

        # Estimate BD rates
        sister_clade_tree_i <- phytools::force.ultrametric(tree = sister_clade_tree_i, method = "nnls", message = F)
        sister_clade_BD_i <- phytools::fit.bd(tree = sister_clade_tree_i, rho = rho)

        ## Draw the divergence time on the stem branch
        # ?TreeSim::corsim

        new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = sister_clade_BD_i$b, mu = sister_clade_BD_i$d,
                                                  missing = ntips_missing, told = upper_age_i, tyoung = lower_age_i)
        divergence_time_missing_clade_i <- setdiff(new_divergence_times_i, divergence_times_i)

        # divergence_time_missing_clade_i

        # Save BD parameters
        BD_i <- sister_clade_BD_i

      } else {

        ## Case with only 2 lineages in the crown clade: fit a BD on a substitute clade of at least 3 tips, encompassing the crown clade

        # Retrieve MRCA node of the extended substitute clade
        substitute_clade_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i), 1]

        # Extract the substitute tree
        substitute_clade_tree_i <- ape::extract.clade(phy = phy, node = substitute_clade_node_i)

        # plot(substitute_clade_tree_i)

        # Reidentify number of terminals in the BD clade
        nb_desc_tips_substitute_clade <- length(substitute_clade_tree_i$tip.label)
        # Adjust sampling fraction
        rho_BD_clade <- nb_desc_tips_substitute_clade / (nb_desc_tips_substitute_clade + ntips_missing)

        # Estimate BD rates from the substitute tree including at least 3 tips
        substitute_clade_tree_i <- phytools::force.ultrametric(tree = substitute_clade_tree_i, method = "nnls", message = F)
        substitute_clade_BD_i <- phytools::fit.bd(tree = substitute_clade_tree_i, rho = rho_BD_clade)

        ## Draw divergence times on the stem branch

        # Draw divergence ages of missing taxa
        # ?TreeSim::corsim
        new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = substitute_clade_BD_i$b, mu = substitute_clade_BD_i$d,
                                                  missing = ntips_missing, told = upper_age_i, tyoung = lower_age_i)
        divergence_time_missing_clade_i <- setdiff(new_divergence_times_i, divergence_times_i)

        # divergence_time_missing_clade_i

        # Save BD parameters
        BD_i <- substitute_clade_BD_i

      }

    } else {

      ## Case with a single tip => MRCA node = tip node ; no descendants

      # MRCA node = unique tip node
      MRCA_node_i <- which(phy$tip.label == MRCA_tips_i)

      # Sister-clade tree is a single branch
      # Divergence times = the stem age
      stem_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i),1]
      stem_age_i <- full_tree_age - phytools::nodeheight(tree = phy, node = stem_node_i)
      divergence_times_i <- stem_age_i

      # Set upper bound to drawn divergence age as the stem age
      upper_age_i <- stem_age_i

      # Compute sampling fraction as if missing only one lineage (the branch leading to the missing clade)
      ntips_present <- 1
      ntips_missing <- 1
      rho <- ntips_present / (ntips_present + ntips_missing)

      ## Fit a BD on a substitute clade of at least 3 tips, encompassing the target grafting tip

      # Retrieve MRCA node of the extended substitute clade
      substitute_clade2_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i), 1]
      substitute_clade3_node_i <- phy$edge[(phy$edge[,2] == substitute_clade2_node_i), 1]

      # Extract the substitute tree
      substitute_clade_tree_i <- ape::extract.clade(phy = phy, node = substitute_clade3_node_i)

      # plot(substitute_clade_tree_i)

      # Reidentify number of terminals in the BD clade
      nb_desc_tips_substitute_clade <- length(substitute_clade_tree_i$tip.label)
      # Adjust sampling fraction
      rho_BD_clade <- nb_desc_tips_substitute_clade / (nb_desc_tips_substitute_clade + ntips_missing)

      # Estimate BD rates from the substitute tree including at least 3 tips
      substitute_clade_tree_i <- phytools::force.ultrametric(tree = substitute_clade_tree_i, method = "nnls", message = F)
      substitute_clade_BD_i <- phytools::fit.bd(tree = substitute_clade_tree_i, rho = rho_BD_clade)

      ## Draw divergence times

      # Draw divergence ages of missing taxa
      # ?TreeSim::corsim
      new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = substitute_clade_BD_i$b, mu = substitute_clade_BD_i$d,
                                                missing = ntips_missing, told = upper_age_i, tyoung = 0)
      divergence_time_missing_clade_i <- setdiff(new_divergence_times_i, divergence_times_i)

      # divergence_time_missing_clade_i

      # Save BD parameters
      BD_i <- substitute_clade_BD_i

    }

    # Bind missing taxa on backbone tree
    nb_missing_taxa <- length(missing_taxa_list_i$Taxa_name)

    if (nb_missing_taxa == 1)
    {
      ## Case for a single tip

      # Extract the MRCA age
      MRCA_edge_i <- which((phy$edge[, 2] == MRCA_node_i))
      root_age <- max(phytools::nodeHeights(phy))
      MRCA_age_i <- root_age - phytools::nodeHeights(phy)[MRCA_edge_i,2]

      # Compute position on the stem branch
      position_i <- divergence_time_missing_clade_i - MRCA_age_i

      # Bind missing clade on phylogeny
      phy <- phytools::bind.tip(tree = phy, tip.label = missing_taxa_list_i$Taxa_name,
                                edge.length = divergence_time_missing_clade_i,
                                where = MRCA_node_i, # MRCA is the tipward node use to locate the branch for grafting
                                position = position_i) # Distance from tipward node to root

      # plot(phy)

    } else {

      ## Case for a clade

      ## Grow a tree with the requested number of missing taxa
      missing_clade_i <- TreeSim::sim.bd.taxa(n = nb_missing_taxa, numbsim = 1,
                                              lambda = BD_i$b, mu = min(BD_i$d, 0.9*BD_i$b),
                                              frac = 1, complete = FALSE)[[1]]
      # Assign taxa names as tip labels
      missing_clade_i$tip.label <- missing_taxa_list_i$Taxa_name

      # plot(missing_clade_i)

      # ## Use the speciation rate of the sister-clade to define the crown age
      #
      # # # Plot the associated density probability function
      # # plot(x = seq(from = divergence_time_missing_clade_i, to = 0, by = -divergence_time_missing_clade_i/100),
      # #      y = dexp(x = seq(from = 0, to = divergence_time_missing_clade_i, by = divergence_time_missing_clade_i/100), rate = BD_i$b),
      # #      xlim = c(divergence_time_missing_clade_i, 0),
      # #      type = "l", xlab = "Crown age", ylab = "Probability density")
      #
      # # Draw a random position following an exponential law in accordance to the fitted speciation rate in the clade
      # # crown_age cannot be negative (cannot excess the stem age)
      # crown_age <- -Inf
      # while (crown_age <= 0)
      # {
      #   position <- rexp(n = 1, rate = BD_i$b)
      #   crown_age <- divergence_time_missing_clade_i - position
      # }

      ## Stretch the grown tree so its stem age matches the divergence age

      # Get current stem age of the grown tree
      current_crown_age_i <- max(phytools::nodeHeights(missing_clade_i))
      current_stem_age_i <- current_crown_age_i + missing_clade_i$root.edge

      # Compute stretch factor
      stretch_factor <- divergence_time_missing_clade_i/current_stem_age_i

      # Adjust internal branch lengths
      missing_clade_i$edge.length <- missing_clade_i$edge.length*stretch_factor
      # Adjust stem branch length
      missing_clade_i$root.edge <- missing_clade_i$root.edge*stretch_factor

      # # Check adjustment
      # (max(phytools::nodeHeights(missing_clade_i)) + missing_clade_i$root.edge) == divergence_time_missing_clade_i

      ## Graft missing taxa using the divergence time between the sister-clade and the missing clade

      # Extract the MRCA age
      MRCA_edge_i <- which((phy$edge[, 2] == MRCA_node_i))
      root_age <- max(phytools::nodeHeights(phy))
      MRCA_age_i <- root_age - phytools::nodeHeights(phy)[MRCA_edge_i,2]

      # Compute position on the stem branch
      position_i <- divergence_time_missing_clade_i - MRCA_age_i

      # Bind missing clade on phylogeny
      phy <- ape::bind.tree(x = phy, y = missing_clade_i,
                            where = MRCA_node_i, # MRCA is the tipward node use to locate the branch for grafting
                            position = position_i) # Distance from tipward node to root

      # plot(phy)
    }

    # Update number of missing taxa / missing clades grafted
    nb_taxa_grafted <- nb_taxa_grafted + length(missing_taxa_list_i$Taxa_name)
    nb_missing_clades_grafted <- nb_missing_clades_grafted + 1

    # Print progress every missing clades
    if (verbose & (nb_missing_clades_grafted %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Missing taxa grafted within clades n°",nb_missing_clades_grafted,"/",length(groups_ID)," - Taxa n°",nb_taxa_grafted,"/",length(missing_taxa_list$Taxa_name),"\n"))
    }

    # Save temp phy in case of crash
    if (save_temp & (nb_missing_clades_grafted %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Temporary grafted phylogeny saved in './temp_phy_from_grafting.rds'\n"))
      saveRDS(object = phy, file = "./temp_phy_from_grafting.rds")
    }
  }

  # Force ultrametry to avoid issue with bad rounding if requested
  if (force_ultrametry)
  {
    phy <- phytools::force.ultrametric(tree = phy, method = "nnls", message = F)
  }

  # Return new phylogeny
  return(phy)
}


# ### Toyset ###
#
# missing_taxa_list_test <- list(Taxa_name = c("Q1", "Q2", "D1", "D2", "D3", "A1"),
#                                MRCA_tips = list(Q1 = c(Tip_1 = "Quiberoni_bigoudensis", Tip_2 = "Quiberoni_creposucri"),
#                                                 Q2 = c(Tip_1 = "Quiberoni_bigoudensis", Tip_2 = "Quiberoni_creposucri"),
#                                                 D1 = c(Tip_1 = "Duckeanthidium_amazonense", Tip_2 = "Duckeanthidium_paraense"),
#                                                 D2 = c(Tip_1 = "Duckeanthidium_amazonense", Tip_2 = "Duckeanthidium_paraense"),
#                                                 D3 = c(Tip_1 = "Duckeanthidium_amazonense", Tip_2 = "Duckeanthidium_paraense"),
#                                                 A1 = c(Tip_1 = "Anthodioctes_mapirensis")))
#
# missing_taxa_list_test <- list(Taxa_name = c("Q1", "Q2", "D1", "D2", "D3", "A1"),
#                                MRCA_tips = list(Q1 = c(Tip_1 = "Quiberoni_bigoudensis", Tip_2 = "Quiberoni_creposucri"),
#                                                 Q2 = c(Tip_1 = "Quiberoni_bigoudensis", Tip_2 = "Quiberoni_creposucri"),
#                                                 D1 = c(Tip_1 = "Anthodioctes_mapirensis"),
#                                                 D2 = c(Tip_1 = "Anthodioctes_mapirensis"),
#                                                 D3 = c(Tip_1 = "Anthodioctes_mapirensis"),
#                                                 A1 = c(Tip_1 = "Duckeanthidium_amazonense", Tip_2 = "Duckeanthidium_paraense")))
#
# grafted_phy_test <- grafting_missing_clades(phy = test_main_phy,
#                                             missing_taxa_list = missing_taxa_list_test,
#                                             seed = 25665)
#
# plot(test_main_phy)
# plot(grafted_phy_test)
#

