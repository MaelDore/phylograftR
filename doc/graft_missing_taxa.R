## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----setup_data, fig.width = 6, fig.height = 5, out.width = "100%"------------
data(toy_backbone)

# Plot the backbone with its node IDs, to identify the clades used below
ape::plot.phylo(toy_backbone, cex = 0.6, no.margin = TRUE)
ape::nodelabels(cex = 0.6, frame = "rect", bg = "lightblue", col = "darkblue")

## ----helpers------------------------------------------------------------------
# Helper 1: build a missing_taxa_list of n taxa to graft in a target clade defined by tips
make_list <- function(n, target_tips, prefix = "New_sp", excluded = NULL)
{
  taxa <- paste0(prefix, seq_len(n))
  out <- list(Taxon_name = taxa,
              Target_clade = setNames(rep(list(target_tips), n), taxa))
  if (!is.null(excluded)) out$Excluded_subclades <- setNames(rep(list(excluded), n), taxa)
  return(out)
}

# Helper 2: for each taxon, is it nested within a clade of the grafted phylogeny (defined by tips)?
is_nested_in <- function(phy_grafted, clade_tips, taxa)
{
  clade_node <- ape::getMRCA(phy_grafted, clade_tips)
  taxa %in% ape::extract.clade(phy_grafted, clade_node)$tip.label
}

# Helper 3: plot two grafted phylogenies side by side
plot_two <- function(phy_1, title_1, phy_2, title_2)
{
  par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
  plot_grafted(phy_1, use_simmap = FALSE, cex = 0.55, main = title_1)
  plot_grafted(phy_2, use_simmap = FALSE, cex = 0.55, main = title_2)
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
}

## ----stem, fig.width = 10, fig.height = 5.5, out.width = "100%"---------------
# ------ Stem of the target clade: available or not ------ #

# Create list of taxa to graft
taxa_gen1 <- make_list(n = 6, target_tips = c("t13", "t27"))

# Graft with target clade's stem available (default)
grafted_stem <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_gen1,
                                   include_target_clade_stem = TRUE, seed = 2, verbose = FALSE)

# Graft with target clade's stem not available
grafted_no_stem <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_gen1,
                                      include_target_clade_stem = FALSE, seed = 2, verbose = FALSE)

# Check which taxa ended nested within the target clade
is_nested_in(grafted_stem, c("t13", "t27"), taxa_gen1$Taxon_name)
is_nested_in(grafted_no_stem, c("t13", "t27"), taxa_gen1$Taxon_name)

plot_two(grafted_stem, "Stem available", grafted_no_stem, "Stem not available")

## ----stem_per_taxon-----------------------------------------------------------
# Taxon-level values: the stem is available for the first three taxa only
taxa_gen1$include_target_clade_stem <- setNames(c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
                                                taxa_gen1$Taxon_name)
taxa_gen1$include_target_clade_stem

## ----excluded, fig.width = 10, fig.height = 11, out.width = "100%"------------
# ------ Excluded subclades: not used, excluded with stem, excluded without stem ------ #

# Define the excluded subclade using tips
excluded_tips <- c("t28", "t26")
# Check its node ID in the backbone
excluded_node <- ape::getMRCA(toy_backbone, excluded_tips)

# Define the target clade using tips, and check its node ID in the backbone
target_tips <- c("t26", "t18")
target_node <- ape::getMRCA(toy_backbone, target_tips)

taxa_A <- make_list(n = 10, target_tips = target_tips, excluded = list(excluded_tips))

# The excluded subclade is ignored
grafted_not_used <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A,
                                       use_excluded_subclades = FALSE, seed = 202, verbose = FALSE)
# The excluded subclade and its stem are excluded (default)
grafted_excl_stem <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A,
                                        exclude_subclades_stem = TRUE, seed = 202, verbose = FALSE)
# Only the branches within the excluded subclade are excluded
grafted_excl_no_stem <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A,
                                           exclude_subclades_stem = FALSE, seed = 202, verbose = FALSE)

# Summarize where the taxa were grafted in each setup, in a table
grafted_list <- list(grafted_not_used, grafted_excl_stem, grafted_excl_no_stem)
excluded_summary <- data.frame(
  Setup = c("Excluded subclade ignored", "Excluded with its stem", "Excluded without its stem"),
  # Number of taxa nested within the excluded subclade
  nb_nested = sapply(grafted_list, function(phy_grafted)
    sum(is_nested_in(phy_grafted, excluded_tips, taxa_A$Taxon_name))),
  # Number of taxa grafted on the stem of the excluded subclade (parent node = the node of the subclade)
  nb_on_stem = sapply(grafted_list, function(phy_grafted)
    sum(extract_grafting_info(phy_grafted, taxa_A, toy_backbone)$parent_node == excluded_node)))
excluded_summary

# Plot the results from the different options
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))

# Top left: the backbone with its node IDs, the target clade (node 32) in limegreen,
# and the excluded subclade (node 38) in orange
ape::plot.phylo(toy_backbone, cex = 0.55, main = "Backbone: target clade (green), excluded subclade (orange)")
ape::nodelabels(cex = 0.55, frame = "rect", bg = "lightblue", col = "darkblue")
ape::nodelabels(node = target_node, cex = 0.55, frame = "rect", bg = "limegreen", col = "black")
ape::nodelabels(node = excluded_node, cex = 0.55, frame = "rect", bg = "orange", col = "black")

# Top right: without excluded subclade
plot_grafted(grafted_not_used, use_simmap = FALSE, cex = 0.55, main = "Excluded subclade ignored")

# Bottom left: exclude subclade and its stem
plot_grafted(grafted_excl_stem, use_simmap = FALSE, cex = 0.55, main = "Excluded with its stem")

# Bottom right: exclude subclade but keep its stem available for grafting
plot_grafted(grafted_excl_no_stem, use_simmap = FALSE, cex = 0.55, main = "Excluded without its stem")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----method, fig.width = 10, fig.height = 5.5, out.width = "100%"-------------
# ------ Method: BD_clades vs uniform ------ #

# Create list of taxa to graft
taxa_A_all <- make_list(n = 8, target_tips = c("t26", "t18"))

# Graft using both methods
grafted_BD <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A_all,
                                 method = "BD_clades", seed = 1234, verbose = FALSE)
grafted_uniform <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A_all,
                                      method = "uniform", seed = 1234, verbose = FALSE)

plot_two(grafted_BD, "method = 'BD_clades'", grafted_uniform, "method = 'uniform'")

# Explore divergence ages of the taxa (the age at which they split from their sister lineage)
# across both methods
ages <- data.frame(
  BD_clades = extract_grafting_info(grafted_BD, taxa_A_all, toy_backbone)$divergence_age,
  uniform = extract_grafting_info(grafted_uniform, taxa_A_all, toy_backbone)$divergence_age)
rownames(ages) <- taxa_A_all$Taxon_name
round(ages, 2)
round(colMeans(ages), 2)

## ----min_tips-----------------------------------------------------------------
# ------ min_tips: clade used to fit the birth-death model ------ #

# Create list of taxa to graft
taxa_small <- make_list(n = 1, target_tips = c("t5", "t12"))

# Graft using different thresholds for the minimum number of tips in enclosing clades
grafted_min3 <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_small,
                                   min_tips = 3, return_BD_model_summary = TRUE,
                                   seed = 1234, verbose = FALSE)
grafted_min5 <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_small,
                                    min_tips = 5, return_BD_model_summary = TRUE,
                                    seed = 1234, verbose = FALSE)

# Summary of the BD models: initial clade vs enclosing clade used to fit the diversification model
cols <- c("initial_clade_n_tips", "enclosing", "fitted_clade_n_tips", "min_tips",
          "net_diversification", "relative_extinction", "rho")
rbind(extract_BD_model_summary(grafted_min3)[, cols], extract_BD_model_summary(grafted_min5)[, cols])

## ----max_age_redraws, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# ------ max_age_redraws: constraints that leave a narrow window for grafting ------ #

# Target clade (node 33) and excluded subclades (nodes 47 and 34), defined by tips
target_tips_B <- c("t8", "t22")
excluded_47 <- c("t5", "t22")
excluded_34 <- c("t8", "t7")

# Check the node IDs in the backbone
sapply(list(target = target_tips_B, excluded_47 = excluded_47, excluded_34 = excluded_34),
       function(tips) ape::getMRCA(toy_backbone, tips))

# Missing taxon with its constraints:
# the stem of node 47 is excluded (TRUE), the stem of node 34 is kept (FALSE)
taxon_narrow <- make_list(n = 1, target_tips = target_tips_B, excluded = list(excluded_47, excluded_34))
taxon_narrow$exclude_subclades_stem <- list(New_sp1 = c(TRUE, FALSE))

# With a low value of max_age_redraws, the function stops with an informative message
tryCatch(
  graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxon_narrow,
                     include_target_clade_stem = FALSE, max_age_redraws = 100,
                     seed = 1234, verbose = FALSE),
  error = function(e) cat("Error:", conditionMessage(e)))

## ----max_age_redraws_grid-----------------------------------------------------
# Try increasing values of max_age_redraws, with the same seed
redraws_values <- c(10, 100, 300, 500, 1000)
redraws_success <- sapply(redraws_values, function(redraws)
  !inherits(try(graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxon_narrow,
                                   include_target_clade_stem = FALSE, max_age_redraws = redraws,
                                   seed = 1234, verbose = FALSE), silent = TRUE), "try-error"))
# Summarize results
data.frame(max_age_redraws = redraws_values, grafted = redraws_success)

## ----max_age_redraws_plot, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# Graft with the default number of redraws
grafted_narrow <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxon_narrow,
                                     include_target_clade_stem = FALSE, seed = 1234, verbose = FALSE)

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
# Left: the backbone with its node IDs, the target clade (node 33) in limegreen,
# the excluded subclade with its stem excluded (node 47) in orange, and with its stem kept (node 34) in gold
ape::plot.phylo(toy_backbone, cex = 0.55, main = "Target (green), excluded (orange, gold)")
ape::nodelabels(cex = 0.55, frame = "rect", bg = "lightblue", col = "darkblue")
ape::nodelabels(node = 33, cex = 0.55, frame = "rect", bg = "limegreen", col = "black")
ape::nodelabels(node = 47, cex = 0.55, frame = "rect", bg = "orange", col = "black")
ape::nodelabels(node = 34, cex = 0.55, frame = "rect", bg = "gold", col = "black")
# Right: the grafted phylogeny, with the taxon sister to the subclade of node 34
plot_grafted(grafted_narrow, use_simmap = FALSE, cex = 0.55, main = "Grafted on the only branch available")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----richness, fig.width = 10, fig.height = 5.5, out.width = "100%"-----------
# ------ Total richness: no unsampled species vs many unsampled species ------ #

# Create list of taxa to graft 
taxa_A3 <- make_list(n = 3, target_tips = c("t26", "t18"), prefix = "Fam_A_sp")
taxa_A3$Target_clade_name <- setNames(rep("Family_A", 3), taxa_A3$Taxon_name)

# Create richness table to account for unsampled species
# Total richness of Family_A: 60 species, so 60 - 17 - 3 = 40 unsampled species
richness_A <- list(Target_clade = list(c("t26", "t18")), Target_clade_name = "Family_A", total_richness = 60)

# Graft with and without accounting for unsampled species
grafted_no_richness <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A3,
                                          return_BD_model_summary = TRUE, seed = 1234, verbose = FALSE)
grafted_richness <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A3,
                                       total_richness_per_clades = richness_A,
                                       return_BD_model_summary = TRUE, seed = 1234, verbose = FALSE)

# Inspect the effective sampling fraction and the rates of the BD models
cols <- c("Clade_name", "n_pending", "n_unsampled", "rho", "b", "d", "net_diversification", "relative_extinction")
rbind(extract_BD_model_summary(grafted_no_richness)[, cols], extract_BD_model_summary(grafted_richness)[, cols])

## The sampling fraction (`rho`) is lower with many unsampled species
## As a result, estimates of speciation (`b`) and extinction (`d`) rates are higher


## ----richness_ages, fig.width = 6, fig.height = 4.5, out.width = "100%"-------
# ------ Total richness: effect on the divergence ages (100 simulations) ------ #

# Graft 100 simulations without and with total richness
grafted_no_richness_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A3,
                                               nb_simulations = 100, seed = 1234, verbose = FALSE)
grafted_richness_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A3,
                                            total_richness_per_clades = richness_A,
                                            nb_simulations = 100, seed = 1234, verbose = FALSE)

# Extract the divergence ages of the 3 taxa in each of the 100 simulations
ages_richness <- list(
  "No total richness" = extract_grafting_info(grafted_no_richness_sims, taxa_A3, toy_backbone)$divergence_age,
  "Total richness = 60" = extract_grafting_info(grafted_richness_sims, taxa_A3, toy_backbone)$divergence_age)

# Summarize the divergence ages
data.frame(Setup = names(ages_richness),
           mean = sapply(ages_richness, mean),
           median = sapply(ages_richness, median),
           quantile_90 = sapply(ages_richness, quantile, probs = 0.9),
           row.names = NULL)

## The divergence ages are younger with many unsampled species, as expected


## ----simmap, fig.width = 6, fig.height = 5, out.width = "100%"----------------
# ------ Load simmap of color states ------ #

# Lod toy simmap
data(toy_simmap)
# Set colors
state_colors <- c(blue = "blue", green = "green3", red = "red")
# Plot simmap
phytools::plotSimmap(toy_simmap, colors = state_colors, fsize = 0.6, lwd = 2, ftype = "i")

## ----states, fig.width = 10, fig.height = 5.5, out.width = "100%"-------------
# ------ Allowed states: all states vs blue only ------ #

# Create list of taxa to graft
taxa_A_states <- make_list(n = 6, target_tips = c("t26", "t18"), prefix = "New_sp")

# A simmap with no allowed states only records the states of the grafted lineages: grafting is not restricted
grafted_all_states <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A_states,
                                         simmap = toy_simmap, allowed_states = NULL,
                                         seed = 67, verbose = FALSE)
# Graft the same taxa, but with only blue as 'allowed_state'
grafted_blue <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_A_states,
                                   simmap = toy_simmap, allowed_states = "blue",
                                   seed = 67, verbose = FALSE)

# The simmap is updated: the grafted tips have the state of the point where they were grafted
phytools::getStates(grafted_all_states$simmap, type = "tips")[taxa_A_states$Taxon_name]
phytools::getStates(grafted_blue$simmap, type = "tips")[taxa_A_states$Taxon_name]

# Helper: plot the updated simmap of a grafted phylogeny
# Grafted branches are dashed, and grafted tips are marked with a point in the color of their final state
plot_states <- function(phy_grafted, main)
{
  plot_grafted(phy_grafted, which = "edges", simmap_colors = state_colors, fsize = 0.55, main = main)
  grafted_tips <- which(phy_grafted$grafted_tips)
  tip_states <- phytools::getStates(phy_grafted$simmap, type = "tips")[phy_grafted$tip.label[grafted_tips]]
  ape::tiplabels(tip = grafted_tips, pch = 19, cex = 1.5, col = state_colors[tip_states])
}

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
# Panel A: all states allowed
plot_states(grafted_all_states, main = "A: All states allowed")
# Panel B: blue only
plot_states(grafted_blue, main = "B: Blue only")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----min_valid_fraction, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# Restrict grafting to the "green" state only
# Green covers a small part of Family_A: a warning lists the taxa with a narrow room for grafting
taxa_restricted <- make_list(n = 4, target_tips = c("t26", "t18"), prefix = "Green_sp")
grafted_restricted <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_restricted,
                                  simmap = toy_simmap, allowed_states = "green",
                                  min_valid_fraction = 0.5, seed = 67, verbose = FALSE)

# Plot the result
plot_states(grafted_restricted, main = "C: Green only")
## All grafted taxa concentrate on the "green" branch portion of
## the backbone tip 't18' in the target clade


