##### Functions to graft missing taxa #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### 1/ Main function to graft missing taxa ####

# Arguments for MRCA, include stem or not, exclusive subclades
# Argument for the type of grafting
# Purely random = "uniform"
# branch selection probability based on branch length (excluding exclusive subclades), position on the branch drawn follows a uniform prior
# Accounting for subclade diversification dynamics = fit birth-death models per subclades
# Per taxa: "BD_taxa"
# Fit BD with one missing taxa
# branch selection probability based on branch length (excluding exclusive subclades), but position on the branch drawn from exponential distribution according to subclade speciation rate (not a uniform prior!)
# Per subclades/MRCA: "BD_clades"
# Account for the proper number of missing taxa during the BD fit
# Estimate age of missing taxa within a clade with TreeSim::corsim, then select randomly the branch to graft (Prevent branching to occur within the exclusive(s) subclade(s))
# Need to graft taxa per MRCA (if nested taxa, need to graft the nested clades first)

grafting_missing_taxa <- function (phy, # Backbone phylogeny
                                   missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
                                   include_stem = T, # Include or not the stem of the total clade when grafting?
                                   use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
                                   method = "uniform", # Either "uniform", "BD_taxa", or "BD_clades"
                                   verbose = T, # To print progress every 100 taxa
                                   save_temp = FALSE, # To save a temporary phylogeny every 10 stem clades
                                   temp_file_path = NULL, # Default = "./temp_phy_from_grafting.rds"
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

  # Check that tips for locating excluding_subclades are all found on the phylogeny
  all_MRCA_tips <- unique(unlist(missing_taxa_list$Excluding_subclades))
  if(!all(all_MRCA_tips %in% phy$tip.label))
  {
    non_matching_tips <- all_MRCA_tips[!(all_MRCA_tips %in% phy$tip.label)]
    stop(paste0("Some tip labels provided to locate nodes of excluding subclades are missing from the phylogeny\n",
                "Missing tips: ", paste(non_matching_tips, collapse = ", ")))
  }

  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Randomize the order of taxa in the missing_taxa_list so it has no influence on the pattern of grafting
  random_order <- sample(x = 1:length(missing_taxa_list$Taxa_name), size = length(missing_taxa_list$Taxa_name), replace = F)
  missing_taxa_list <- lapply(X = missing_taxa_list, FUN = function (x) { x[random_order] })

  # For "uniform" and BD_clades => Loop per taxa

  if (method %in% c("uniform", "BD_taxa"))
  {
    if (verbose)
    {
      cat(paste0(Sys.time()," - Start grafting process\n"))
    }

    for (i in 1:length(missing_taxa_list$Taxa_name))
    {
      # i <- 1

      # Extract taxon name
      focal_taxon_name <- missing_taxa_list$Taxa_name[i]

      # Extract MRCA node
      if (length(missing_taxa_list$MRCA_tips[[i]]) == 2) # If using pairs
      {
        MRCA_node <- ape::getMRCA(phy = phy, tip = missing_taxa_list$MRCA_tips[[i]])
      }
      if (length(missing_taxa_list$MRCA_tips[[i]]) == 1) # If using a single tip
      {
        MRCA_node <- which(phy$tip.label == missing_taxa_list$MRCA_tips[[i]])
      }

      # If using exclusive subclades, retrieve the MRCA nodes of forbidden subclades
      if (use_exclusive_subclades)
      {
        # Extract exclusive subclades pairs for taxon i
        exclusive_subclades_pairs <- missing_taxa_list$Exclusive_subclades[[i]]

        # Initiate vector of excluding subclades MRCA nodes
        exclusive_subclades_MRCA_nodes <- c()

        # Loop per pairs of tips defining excluding subclades
        for (j in 1:length(exclusive_subclades_pairs))
        {
          # j <- 1

          # Extract pair j
          exclusive_subclades_pair_j <- exclusive_subclades_pairs[[j]]
          # Retrieve MRCA node of subclade j
          exclusive_subclades_MRCA_node_j <- ape::getMRCA(phy = phy, tip = exclusive_subclades_pair_j)
          # Store MRCA node of subclade j
          exclusive_subclades_MRCA_nodes <- c(exclusive_subclades_MRCA_nodes, exclusive_subclades_MRCA_node_j)
        }

      } else {
        exclusive_subclades_MRCA_nodes <- NULL
      }

      # Draw node of random branch used to graft the taxa
      selected_node_ID <- drawing_random_branch(phy = phy, MRCA_node = MRCA_node, include_stem = include_stem, exclusive_subclades_MRCA_nodes = exclusive_subclades_MRCA_nodes, seed = NULL)
      # Extract branch length
      selected_branch_ID <- match(x = selected_node_ID, table = phy$edge[, 2])
      selected_branch_length <- phy$edge.length[selected_branch_ID]

      # Draw position of the new tip on the selected branch
      if (method == "uniform")
      {
        # For "uniform" method, draw a random position following a uniform PDF
        position <- runif(n = 1, min = 0, max = selected_branch_length)
      } else {
        # For "BD_taxa" method, draw a random position following an exponential law in accordance to the fitted speciation rate in the clade

        # Identify number of descending taxa in the stem clade
        nb_desc_tips <- length(geiger::tips(phy = phy, node = MRCA_node))

        ## Run BD model according to the number of descending taxa
        # Run BD model on a subclade of at least 3 tips encomapssing the stem group

        if (nb_desc_tips == 1)
        {
          ## Case with a single tip
          # Need to move backward twice to ensure having a clade of at least 3 tips

          # Retrieve MRCA node of the extended substitute clade
          BD_clade_node <- phy$edge[(phy$edge[,2] == MRCA_node), 1]
          BD_clade_node <- phy$edge[(phy$edge[,2] == BD_clade_node), 1]
        }
        if (nb_desc_tips == 2)
        {
          ## Case with two tips
          # Need to move backward once to ensure having a clade of at least 3 tips

          # Retrieve MRCA node of the extended substitute clade
          BD_clade_node <- phy$edge[(phy$edge[,2] == MRCA_node), 1]
        }
        if (nb_desc_tips > 2)
        {
          ## Case with more than two tips
          # Can use the MRCA node
          BD_clade_node <- MRCA_node
        }

        # Extract the crown group or the extended substitute tree
        BD_clade_tree <- ape::extract.clade(phy = phy, node = BD_clade_node)

        # plot(BD_clade_tree)

        # Reidentify number of terminals in the BD clade
        nb_desc_tips_BD_clade <- length(BD_clade_tree$tip.label)

        # Estimate sampling fraction accounting for one missing taxa (the one we are grafting)
        sampling_fraction <- nb_desc_tips_BD_clade / (nb_desc_tips_BD_clade + 1)
        # Fit BD model in the BD clade, accounting for one missing taxa
        BD_clade_tree <- phytools::force.ultrametric(tree = BD_clade_tree, method = "nnls", message = F)
        lambda_BD <- phytools::fit.bd(tree = BD_clade_tree, rho = sampling_fraction)$b

        # Draw a random position following an exponential law in accordance to the fitted speciation rate in the clade
        # plot(x = seq(from = 0, to = 5, by = 0.01), y = dexp(x = seq(from = 0, to = 5, by = 0.01), rate = lambda_BD), type = "l", xlab = "position", ylab = "probability density")
        # Position cannot excess the branch length
        position <- Inf
        while (position > selected_branch_length)
        {
          position <- rexp(n = 1, rate = lambda_BD)
        }
      }

      # Extract root age
      root_age <- max(phytools::nodeHeights(tree = phy)[, 2])

      # Compute branch length needed to remain ultrametric
      new_branch_length <- position + (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))

      # Bind new taxa in the tree
      phy <- phytools::bind.tip(tree = phy, tip.label = focal_taxon_name, edge.length = new_branch_length, where = selected_node_ID, position = position)

      # Check phylo
      # plot(phy)

      # Print progress every 100 tips
      if (verbose & (i %% 100 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted n°",i,"/",length(missing_taxa_list$Taxa_name),"\n"))
      }
    }
  }

  # For "BD_clades" => Loop per stem clades used for grafting, starting with the nested ones

  if (method == "BD_clades")
  {
    if (verbose)
    {
      cat(paste0(Sys.time()," - Start identification of stem clades\n"))
    }

    ## Extract root age of the full tree
    full_tree_age <- max(phytools::nodeHeights(phy)[,2])

    ## Get info per stem clades

    # Identify all stem clades used for grafting
    MRCA_tips_merged <- unlist(lapply(missing_taxa_list$MRCA_tips, FUN = paste0, collapse = ","))
    # Group missing taxa per stem clades
    MRCA_tips_groups <- as.factor(MRCA_tips_merged)
    as.numeric(MRCA_tips_groups)

    # Get info per stem clades
    groups_ID <- 1:nlevels(MRCA_tips_groups)
    group_indices <- match(x = groups_ID, table = as.numeric(MRCA_tips_groups))

    stem_clades_list <- list(ID = groups_ID,
                             MRCA_tips = missing_taxa_list$MRCA_tips[group_indices])
    names(stem_clades_list$MRCA_tips) <- groups_ID

    # Identify MRCA nodes & descending nodes
    stem_clades_list$MRCA_nodes <- list()
    stem_clades_list$descendant_nodes <- list()
    for (i in 1:length(groups_ID))
    {
      if (length(stem_clades_list$MRCA_tips[[i]]) > 1)
      {
        # Case with pairs of tips
        stem_clades_list$MRCA_nodes[[i]] <- ape::getMRCA(phy = phy, tip = stem_clades_list$MRCA_tips[[i]])
        stem_clades_list$descendant_nodes[[i]] <- phytools::getDescendants(tree = phy, node = stem_clades_list$MRCA_nodes[[i]])
      } else {
        # Case with a single tip => MRCA node = tip node ; no descendants
        stem_clades_list$MRCA_nodes[[i]] <- which(phy$tip.label == stem_clades_list$MRCA_tips[[i]])
        stem_clades_list$descendant_nodes[i] <- list(NULL)
      }

    }
    names(stem_clades_list$MRCA_nodes) <- stem_clades_list$ID
    names(stem_clades_list$descendant_nodes) <- stem_clades_list$ID

    ## Order stem groups such as the nested ones are selected fist

    if (verbose)
    {
      cat(paste0(Sys.time()," - Start ordering of stem clades per nestedness\n"))
    }

    # Initiate vector of stem clades order
    stem_clades_order <- c()
    stem_clades_list_temp <- stem_clades_list
    while (length(stem_clades_order) < length(groups_ID))
    {
      # Initiate stopping criterion as check for nestedness of the stem clade
      # Selected stem clade must have no other stem clade nested within
      no_nested_subclades_within <- F
      j <- 0 # Initiate counter
      while (!no_nested_subclades_within)
      {
        j <- j + 1
        target_stem_clade_descendant_nodes <- unlist(stem_clades_list_temp$descendant_nodes[j])
        target_MCRA_other_stem_clades <- unlist(stem_clades_list_temp$MRCA_nodes)[(j+1):length(stem_clades_list_temp$MRCA_nodes)]
        no_nested_subclades_within <- !any(target_MCRA_other_stem_clades %in% target_stem_clade_descendant_nodes)
      }
      # Extract selected stem clade ID
      selected_clade_ID <- stem_clades_list_temp$ID[j]
      # Assign stem clade in the order
      stem_clades_order <- c(stem_clades_order, selected_clade_ID)
      # Remove it from selection
      stem_clades_list_temp <- lapply(X = stem_clades_list_temp, FUN = function (x) { x[-j] })

      # if (verbose & length(stem_clades_order) %% 100 == 0)
      # {
      #   cat(paste0(Sys.time()," - ",length(stem_clades_order), " stem clades out of ",nlevels(MRCA_tips_groups)," ordered by nestedness\n"))
      # }

      # Stop when all stem clades are ordered
    }

    if (verbose)
    {
      cat(paste0(Sys.time()," - All stem clades ordered by nestedness\n"))
    }

    ## Loop per stem clades used for grafting, starting with the nested ones

    if (verbose)
    {
      cat(paste0(Sys.time()," - Start grafting process per clades\n"))
    }

    # Initiate number of taxa/stem clades grafted
    nb_taxa_grafted <- 0
    nb_stem_clades_grafted <- 0

    for (i in stem_clades_order)
    {
      # i <- stem_clades_order[1]

      ## Subset missing_taxa_list for the given stem clade
      missing_taxa_ID_stem_clade_i <- which(as.numeric(MRCA_tips_groups) == i)
      missing_taxa_list_i <- lapply(X = missing_taxa_list, FUN = function (x) { x[missing_taxa_ID_stem_clade_i] })

      # Extract MRCA tips
      MRCA_tips_i <- stem_clades_list$MRCA_tips[[i]]

      ## Draw divergence times according to the number of tips used to define the crown clade, and within the crown clade
      if (length(MRCA_tips_i) > 1)
      {
        ## Case with pairs of tips used to define the crown clade

        # Extract MRCA node using tips (cannot use previous extract as the tree as been updated with the missing taxa from the previous stem group)
        MRCA_node_i <- ape::getMRCA(phy = phy, tip = MRCA_tips_i)

        # Extract the MRCA crown clade tree
        crown_clade_tree_i <- ape::extract.clade(phy = phy, node = MRCA_node_i)

        # plot(crown_clade_tree_i)

        # Extract divergence times in the crown clade tree
        crown_age_i <- max(phytools::nodeHeights(crown_clade_tree_i)[,2])
        divergence_times_i <- crown_age_i - phytools::nodeHeights(crown_clade_tree_i)[,1]

        # Set upper bound to drawn divergence age as the crown age
        upper_age_i <- crown_age_i

        # Add stem branch by using the stem age as the upper limit to draw divergence times instead of the MRCA node age
        if (include_stem)
        {
          stem_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i),1]
          stem_age_i <- full_tree_age - phytools::nodeheight(tree = phy, node = stem_node_i)
          upper_age_i <- stem_age_i # Set upper bound to drawn divergence age as the stem age, thus include the possibility of branching on the stem branch
        }

        # Compute sampling fraction
        ntips_present <- length(crown_clade_tree_i$tip.label)
        ntips_missing <- length(missing_taxa_list_i$Taxa_name)
        rho <- ntips_present / (ntips_present + ntips_missing)

        ## Fit BD on crown clade, or substitute extended crown clade according to number of tips in the crown clade
        if (ntips_present > 2)
        {
          ## Case with more than 2 lineages in the crown clade: fit a BD on the crown clade

          # Estimate BD rates
          crown_clade_tree_i <- phytools::force.ultrametric(tree = crown_clade_tree_i, method = "nnls", message = F)
          crown_clade_BD_i <- phytools::fit.bd(tree = crown_clade_tree_i, rho = rho)

          ## Draw divergence times

          # Draw divergence ages of missing taxa
          # ?TreeSim::corsim
          new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = crown_clade_BD_i$b, mu = crown_clade_BD_i$d, missing = ntips_missing, told = upper_age_i)
          divergence_times_missing_taxa_i <- setdiff(new_divergence_times_i, divergence_times_i)

          # divergence_times_missing_taxa_i

        } else {

          ## Case with only 2 lineages in the crown clade: fit a BD on a substitute clade of at least 3 tips, encompassing the crown clade

          # Retrieve MRCA node of the extended substitute clade
          substitute_clade_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i), 1]

          # Extract the substitute tree
          substitute_clade_tree_i <- ape::extract.clade(phy = phy, node = substitute_clade_node_i)

          # plot(crown_clade_tree_i)
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
          new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = substitute_clade_BD_i$b, mu = substitute_clade_BD_i$d, missing = ntips_missing, told = upper_age_i)
          divergence_times_missing_taxa_i <- setdiff(new_divergence_times_i, divergence_times_i)

          # divergence_times_missing_taxa_i

        }

      } else {

        ## Case with a single tip => MRCA node = tip node ; no descendants

        # MRCA node = unique tip node
        MRCA_node_i <- which(phy$tip.label == MRCA_tips_i)

        # Stem clade tree is a single branch
        # Divergence times = the stem age
        stem_node_i <- phy$edge[(phy$edge[,2] == MRCA_node_i),1]
        stem_age_i <- full_tree_age - phytools::nodeheight(tree = phy, node = stem_node_i)
        divergence_times_i <- stem_age_i

        # Set upper bound to drawn divergence age as the stem age
        upper_age_i <- stem_age_i

        # Compute sampling fraction
        ntips_present <- 1
        ntips_missing <- length(missing_taxa_list_i$Taxa_name)
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
        new_divergence_times_i <- TreeSim::corsim(x = divergence_times_i, lambda = substitute_clade_BD_i$b, mu = substitute_clade_BD_i$d, missing = ntips_missing, told = upper_age_i)
        divergence_times_missing_taxa_i <- setdiff(new_divergence_times_i, divergence_times_i)

        # divergence_times_missing_taxa_i

      }

      ## Graft missing taxa using the divergence times

      # Graft new taxa from oldest to youngest to allow the possibility to graft missing taxa on missing taxa branches!
      for (j in 1:length(missing_taxa_list_i$Taxa_name))
      {
        divergence_age <- divergence_times_missing_taxa_i[j]
        tip_label <- missing_taxa_list_i$Taxa_name[j]
        MRCA_tips <- missing_taxa_list_i$MRCA_tips[[j]]
        excluding_subclades_tips_list <- missing_taxa_list_i$Excluding_subclades[[j]]

        phy <- grafting_from_divergence_age(phy = phy, # Backbone phylogeny
                                            divergence_age = divergence_age, # Divergence age
                                            tip_label = tip_label, # Missing taxa label
                                            MRCA_tips = MRCA_tips, # To provide the pair of tips defining the MRCA node of the total clade within which to graft the missing taxa
                                            excluding_subclades_tips_list = excluding_subclades_tips_list, # List of pairs defining subclades to exclude from grafting options
                                            seed = NULL) # Set seed outside of the loop
      }

      # Update number of missing taxa / stem clades grafted
      nb_taxa_grafted <- nb_taxa_grafted + length(missing_taxa_list_i$Taxa_name)
      nb_stem_clades_grafted <- nb_stem_clades_grafted + 1

      # Print progress every stem groups
      if (verbose & (nb_stem_clades_grafted %% 1 == 0))
      {
        cat(paste0(Sys.time(), " - Missing taxa grafted within clades n°",nb_stem_clades_grafted,"/",length(groups_ID)," - Taxa n°",nb_taxa_grafted,"/",length(missing_taxa_list$Taxa_name),"\n"))
      }

      # Save temp phy in case of crash
      if (save_temp & (nb_stem_clades_grafted %% 1 == 0))
      {
        if (is.null(temp_file_path))
        {
          temp_file_path <- "./temp_phy_from_grafting.rds"
        }
        cat(paste0(Sys.time(), " - Temporary grafted phylogeny saved in '",temp_file_path,"'\n"))
        saveRDS(object = phy, file = temp_file_path)
      }
    }
  }

  # Force ultrametry to avoid issue with bad rounding
  phy <- phytools::force.ultrametric(tree = phy, method = "nnls", message = F)

  # Return new phylogeny
  return(phy)
}


### 2/ Sub-function to draw randomly a branch within a clade ####

# phy <- sim.bdtree(b = 1, d = 0, n = 10)
# plot(phy)
# nodelabels()
# edgelabels()
#
# MRCA_node <- 18
# MRCA_node <- 5
# MRCA_node <- 12
# exclusive_subclades_MRCA_nodes <- c(13, 16)

drawing_random_branch <- function (phy,
                                   MRCA_node, # To provide the MRCA node of the total clades where to draw the selected branch from
                                   which_node = "tipward", # Chose to provide the rootward node or the tipward node of the selected branch
                                   include_stem = T, # Include the stem branch among the set of possible branches
                                   exclusive_subclades_MRCA_nodes = NULL, # To provide MRCA nodes of exclusive subclades where not to draw the selected branch
                                   seed = 1234) # Set the seed
{
  # Node ID: N_tips (N = N_tips), root (N = 1), N_internal (N = N_tips - 2)
  # Total = N_tips - 1

  # N_edges = N_tips - 2 because the root node does not have an edge
  # By default, edges are not ordered as internal/terminal, but rather cladewise

  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Extract descendant nodes from MRCA
  descendant_nodes <- phytools::getDescendants(tree = phy, node = MRCA_node)

  # If MRCA is the root, prevent from drawing the stem (does not exist)
  root_node <- length(phy$tip.label) + 1
  if (MRCA_node == root_node)
  {
    include_stem <- F
  }

  # Include node of the stem if needed (= MRCA node)
  if (include_stem)
  {
    descendant_nodes <- unique(c(descendant_nodes, MRCA_node))
  }
  # Extract branches ID
  branches_ID <- match(x = descendant_nodes, table = phy$edge[, 2])
  # Extract branch length
  branches_length <- setNames(object = phy$edge.length[branches_ID], nm = branches_ID)

  # Remove branch from exclusive subclades from the draw
  if (!is.null(exclusive_subclades_MRCA_nodes))
  {
    # Initiate vector of forbidden branches
    excluded_branches_ID <- c()
    # Loop per exclusive subclade
    for (i in 1:length(exclusive_subclades_MRCA_nodes))
    {
      # Extract exclusive subclade MRCA node
      exclusive_subclade_MRCA_node_i <- exclusive_subclades_MRCA_nodes[i]
      # Extract descendant nodes from exclusive subclade MRCA
      exclusive_subclade_descendant_nodes_i <- phytools::getDescendants(tree = phy, node = exclusive_subclade_MRCA_node_i)
      # Extract forbidden branches ID
      excluded_branches_ID_i <- match(x = exclusive_subclade_descendant_nodes_i, table = phy$edge[, 2])
      # Store in vector of forbidden branches
      excluded_branches_ID <- c(excluded_branches_ID, excluded_branches_ID_i)
    }

    # Remove excluded branches from the draw
    branches_ID <- setdiff(branches_ID, excluded_branches_ID)
    branches_length <- branches_length[!(names(branches_length) %in% excluded_branches_ID)]
  }

  # Draw randomly the branch
  if (length(branches_ID) == 1)
  {
    selected_branch_ID <- branches_ID
  } else {
    selected_branch_ID <- sample(x = branches_ID, size = 1, prob = branches_length)
  }

  # Provide the tipward or rootward node
  if (which_node == "tipward")
  {
    selected_node <- phy$edge[selected_branch_ID, 2]
  } else {
    selected_node <- phy$edge[selected_branch_ID, 1]
  }

  # Return selected node ID
  return(selected_node)
}


### 3/ Sub-function grafting missing taxa from divergence age #####

grafting_from_divergence_age <- function (phy, # Backbone phylogeny
                                          divergence_age, # Divergence age
                                          tip_label, # Missing taxa label
                                          MRCA_tips, # To provide the pair of tips defining the MRCA node of the total clade within which to graft the missing taxa. Can also provide a single tip used to graft the missing taxa on.
                                          excluding_subclades_tips_list, # List of pairs defining subclades to exclude from grafting options
                                          seed = 1234)

{
  # Set the seed
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Extract root age
  root_age <- max(phytools::nodeHeights(tree = phy)[, 2])

  # Subset tipward nodes of branches crossing the target divergence age
  branch_times_df <- as.data.frame(root_age - phytools::nodeHeights(tree = phy))
  branch_times_df$rootward_age_match <- branch_times_df[,1] >= divergence_age
  branch_times_df$tipward_age_match <- branch_times_df[,2] <= divergence_age
  branch_times_df$compatible <- branch_times_df$rootward_age_match & branch_times_df$tipward_age_match
  compatible_branches_ID <- which(branch_times_df$compatible)
  compatible_nodes_ID <- phy$edge[compatible_branches_ID, 2]

  # Case with two tips used to define MRCA of grafting clade
  if (length(MRCA_tips) == 2)
  {
    # Get MRCA node
    MRCA_node <- ape::getMRCA(phy = phy, tip = MRCA_tips)
    # Extract descendant nodes from MRCA
    descendant_nodes <- phytools::getDescendants(tree = phy, node = MRCA_node)
    # Include node of the stem (= MRCA node)
    stem_clade_nodes <- unique(c(descendant_nodes, MRCA_node))
  } else {
    # Case with a single grafting tip
    MRCA_node <- which(phy$tip.label == MRCA_tips)
    # Stem clade nodes = stem node of the unique grafting tip
    stem_clade_nodes <- MRCA_node
  }

  # Subset nodes of clade to use for grafting as the intersect between the nodes within the stem clades, and the nodes of branching crossing the divergence age
  valid_nodes_ID <- intersect(compatible_nodes_ID, stem_clade_nodes)

  # Exclude nodes from select options based on list of pairs defining excluding subclades
  if (length(excluding_subclades_tips_list) > 0)
  {
    # Initiate vector of forbidden nodes
    excluded_nodes_ID <- c()
    # Loop per exclusive subclade
    for (j in 1:length(excluding_subclades_tips_list))
    {
      # j <- 1

      # Extract exclusive subclade tips
      excluding_subclade_tips_j <- excluding_subclades_tips_list[[j]]
      # Retrieve exclusive subclade MRCA node
      exclusive_subclade_MRCA_node_j <- ape::getMRCA(phy = phy, tip = excluding_subclade_tips_j)
      # Extract descendant nodes from exclusive subclade MRCA
      excluded_nodes_ID_i <- phytools::getDescendants(tree = phy, node = exclusive_subclade_MRCA_node_j)
      # Store in vector of forbidden branches
      excluded_nodes_ID <- c(excluded_nodes_ID, excluded_nodes_ID_i)
    }

    # Remove excluded nodes from the draw
    valid_nodes_ID <- setdiff(valid_nodes_ID, excluded_nodes_ID)
  }

  # Select randomly the branch/node among the valid options
  selection_index <- sample(x = 1:length(valid_nodes_ID), size = 1)
  selected_node_ID <- valid_nodes_ID[selection_index]

  # Compute position needed for the tree to remain ultrametric such as position + selected node age = divergence age
  position <- divergence_age - (root_age - phytools::nodeheight(tree = phy, node = selected_node_ID))

  # Bind new taxa in the tree
  phy <- phytools::bind.tip(tree = phy, tip.label = tip_label, edge.length = divergence_age, where = selected_node_ID, position = position)

  # Force ultrametry to avoid issue with bad rounding
  phy <- phytools::force.ultrametric(tree = phy, method = "nnls", message = F)

  # Return phylogeny with the new grafted taxa
  return(phy)
}


# ### 4/ Example on a toyset ####
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
# exclusive_subclades_MRCA_nodes <- c(13, 16)
#
# missing_taxa_list <- list(Taxa_name = c("A", "B", "C"),
#                           MRCA_tips = list(c("s3", "s10"), c("s3"), c("s3", "s8")),
#                           Exclusive_subclades = list(list(c("s3", "s8"), c("s9", "s10")), NULL, list(c("s7", "s3"))))
#
#
#
# grafted_phy <- grafting_missing_taxa(phy = test_phy, # Backbone phylogeny
#                                      missing_taxa_list, # List with taxa to graft. Items: Taxa_name as c(), MRCA_tips as list(c(), c(), ...), Exclusive_subclades as list(list(c(), c()), list(c()), ...)
#                                      include_stem = T, # Include or not the stem of the total clade when grafting?
#                                      use_exclusive_subclades = T, # Use Exclusive_subclades to avoid grafting within specific subclades
#                                      method = "uniform", # Either "uniform", "BD_taxa", or "BD_clades"
#                                      verbose = T,
#                                      seed = 1234) # Set the seed
#
# plot(grafted_phy)
#
