## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----data---------------------------------------------------------------------
# ------ Toy data ------ #

data(toy_backbone)
data(toy_tip_data)
data(toy_multiSimmap)
data(toy_multiSimmap_elevation)
data(toy_multiSimmap_ranges)

# Tip data of the three attributes
head(toy_tip_data)

# Inspect distribution of data across tips
lapply(toy_tip_data, table)

## ----tip_vectors--------------------------------------------------------------
# ------ Named vectors of tip data ------ #

# Prepare tip data as named vectors
tip_colors <- setNames(toy_tip_data$color, row.names(toy_tip_data))
tip_elevations <- setNames(toy_tip_data$elevation, row.names(toy_tip_data))
tip_ranges <- setNames(toy_tip_data$range, row.names(toy_tip_data))

# The names are the tip labels of the backbone
all(names(tip_colors) == toy_backbone$tip.label)

## ----colors-------------------------------------------------------------------
# ------ Colors of the states ------ #

# Colors of the three states of the color trait
color_states <- c(blue = "dodgerblue", green = "limegreen", red = "firebrick")

# Colors of the two states of the elevation
elevation_states <- c(high = "darkorange", low = "steelblue")

# Colors of the ranges: same color scheme as in vignette("biogeographic_data") (A, B, AB, BC),
# completed with the other ranges that can be found in the maps (C, AC, and the null range "_")
range_colors <- c(A = "gold", B = "dodgerblue", AB = "limegreen", BC = "firebrick2",
                  C = "mediumpurple", AC = "darkorange", "_" = "grey50")

## ----helpers------------------------------------------------------------------
# ------ Helpers to plot the maps before and after grafting ------ #

# Plot an initial simmap of the backbone
plot_initial_map <- function (simmap, colors, main)
{
  phytools::plotSimmap(simmap, colors = colors, fsize = 0.5, ftype = "i", lwd = 2)
  title(main, cex.main = 1, line = -1) # Inside the plot region, as plotSimmap() resets the margins
}

# Plot the updated simmap of a grafted phylogeny
# 'trait' is the name (or position) of the trait to plot when the maps of several traits were used
plot_updated_map <- function (phy_grafted, colors, main, trait = 1)
{
  # Grafted branches are dashed
  plot_grafted(phy_grafted, which = "edges", trait = trait, simmap_colors = colors,
               fsize = 0.5, lwd = 2)
  title(main, cex.main = 1, line = -1)
  # Grafted tips are marked with a point in the color of the state they inherited
  simmap <- phy_grafted$simmap
  if (!inherits(simmap, "simmap")) { simmap <- simmap[[trait]] } # Maps of several traits
  grafted_tips <- which(phy_grafted$grafted_tips)
  tip_states <- phytools::getStates(simmap, type = "tips")[phy_grafted$tip.label[grafted_tips]]
  ape::tiplabels(tip = grafted_tips, pch = 19, cex = 1.2, col = colors[tip_states])
}

## ----deepSTRAPP_color, eval = FALSE-------------------------------------------
# # ------ Not run: build the maps of the color with deepSTRAPP ------ #
# 
# color_data <- deepSTRAPP::prepare_trait_data(tip_data = tip_colors, # Named vector of tip states
#                                              trait_data_type = "categorical",
#                                              phylo = toy_backbone, # Time-calibrated phylogeny
#                                              evolutionary_models = c("ER", "SYM", "ARD"), # Models to compare. Default is "ARD"
#                                              seed = 1234, # For reproducibility
#                                              nb_simulations = 100, # Number of stochastic maps
#                                              return_simmaps = TRUE, # Needed to get the maps
#                                              return_model_selection_df = TRUE, # Needed to get the comparison of models
#                                              plot_map = FALSE, plot_overlay = FALSE)
# 
# # Comparison of the models (ER is the best model here)
# color_data$model_selection_df
# 
# # The stochastic maps: a multiSimmap of 100 maps
# my_multiSimmap <- color_data$simmaps

## ----single_simmap------------------------------------------------------------
# ------ One simmap ------ #

# Extract one map (with [[ ]]): an object of class "simmap"
my_simmap <- toy_multiSimmap[[1]]
class(my_simmap)

# Time spent in each state (sum on all branches)
round(colSums(my_simmap$mapped.edge) / sum(my_simmap$mapped.edge), 2)

# The tips in the simmap are in the states of the tip data
all(phytools::getStates(my_simmap, type = "tips")[names(tip_colors)] == tip_colors)

## ----single_graft-------------------------------------------------------------
# ------ Graft with a single simmap ------ #

# Missing taxa with their allowed colors
data(toy_missing_taxa_list_with_allowed_states)
taxa_colors <- toy_missing_taxa_list_with_allowed_states

# New_sp1 is only allowed in "red" in the current mising_taxa_list. We modify it to also allow "green", 
# so that a branch is available for it in each of the 10 maps used in Section 3 (see vignette("check_grafting"))
taxa_colors$allowed_states[[1]] <- c("red", "green")

# Graft the taxa using the map as constraint
grafted_single <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                     simmap = my_simmap, verbose = FALSE, seed = 1)

# The output includes the updated map of the grafted phylogeny
class(grafted_single$simmap)

## ----single_plot, fig.width = 10, fig.height = 5.5, out.width = "100%"--------
# ------ Compare the initial map and the updated map ------ #

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot_initial_map(my_simmap, colors = color_states, main = "A: Initial map")
plot_updated_map(grafted_single, colors = color_states, main = "B: Grafted map")
legend("bottomleft", legend = names(color_states), fill = color_states, title = "Color", bty = "n", cex = 0.8)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----multiSimmap--------------------------------------------------------------
# ------ A multiSimmap ------ #

class(toy_multiSimmap)
length(toy_multiSimmap)

# All the maps are in the same states
unique(unlist(lapply(toy_multiSimmap, function (x) colnames(x$mapped.edge))))

## ----multiSimmap_graft--------------------------------------------------------
# ------ Graft with a multiSimmap ------ #

# Graft the taxa with each of the 10 maps (one simulation per map by default)
grafted_multi <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                    simmap = toy_multiSimmap, verbose = FALSE, seed = 1)
length(grafted_multi)

## ----multiSimmap_plot, fig.width = 10, fig.height = 10, out.width = "100%"----
# ------ Initial maps (upper line) and grafted maps (lower line) ------ #

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
# Upper line: initial maps
for (i in 1:2) { plot_initial_map(toy_multiSimmap[[i]], colors = color_states, main = paste("Initial map", i)) }
# Lower line: grafted maps
for (i in 1:2) { plot_updated_map(grafted_multi[[i]], colors = color_states, main = paste("Grafted map", i)) }
legend("bottomleft", legend = names(color_states), fill = color_states, title = "Color", bty = "n", cex = 0.8)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----multiSimmap_pies, fig.width = 6, fig.height = 6, out.width = "60%"-------
# ------ Posterior probability of the states at the nodes ------ #

node_probs <- phytools::describe.simmap(toy_multiSimmap, plot = FALSE)$ace

phytools::plotTree(toy_backbone, ftype = "i", fsize = 0.6, lwd = 1)
ape::nodelabels(pie = node_probs, piecol = color_states[colnames(node_probs)], cex = 0.5)

## ----deepSTRAPP_elevation, eval = FALSE---------------------------------------
# # ------ Not run: build the maps of the elevation with deepSTRAPP ------ #
# 
# elevation_data <- deepSTRAPP::prepare_trait_data(tip_data = tip_elevations,
#                                                  trait_data_type = "categorical",
#                                                  phylo = toy_backbone,
#                                                  evolutionary_models = c("ER", "ARD"),
#                                                  seed = 2345,
#                                                  nb_simulations = 100, # The same number of maps as for the color
#                                                  return_simmaps = TRUE,
#                                                  plot_map = FALSE, plot_overlay = FALSE)
# 
# # Comparison of the models (ER is the best model here)
# elevation_data$model_selection_df
# 
# # The stochastic maps: a multiSimmap of 100 maps
# my_multiSimmap_elevation <- elevation_data$simmaps

## ----two_traits---------------------------------------------------------------
# ------ List of simmaps: one element per trait ------ #

# Combine map 1 for color and map 1 for elevation in a list of two simmaps
maps_two_traits_single <- list(color = toy_multiSimmap[[1]], elevation = toy_multiSimmap_elevation[[1]])

# States of each trait
lapply(maps_two_traits_single, function (x) colnames(x$mapped.edge))

## ----two_traits_graft---------------------------------------------------------
# ------ Graft with two traits ------ #

# Missing taxa with the allowed states of the color and of the elevation
data(toy_missing_taxa_list_with_two_traits)
taxa_two_traits <- toy_missing_taxa_list_with_two_traits

# New_sp1 is only allowed in "red" in the current mising_taxa_list. We modify it to also allow "green", 
# so that a branch is available for it in each of the 10 pairs of maps used here
taxa_two_traits$allowed_states[[1]] <- list(c("red", "green"), NULL)

# Graft the taxa using the map of each trait
grafted_two_traits <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_two_traits,
                                         simmap = maps_two_traits_single, verbose = FALSE, seed = 1)

# We obtain one grafted phylogeny, with one updated map per trait
names(grafted_two_traits$simmap)

## ----two_traits_plot, fig.width = 10, fig.height = 10, out.width = "100%"-----
# ------ Initial maps (upper line) and grafted maps (lower line) ------ #

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
# Upper line: initial maps of the color and of the elevation
plot_initial_map(maps_two_traits_single$color, colors = color_states, main = "Initial map: color")
plot_initial_map(maps_two_traits_single$elevation, colors = elevation_states, main = "Initial map: elevation")
# Lower line: grafted maps
plot_updated_map(grafted_two_traits, colors = color_states, trait = "color", main = "Grafted map: color")
plot_updated_map(grafted_two_traits, colors = elevation_states, trait = "elevation", main = "Grafted map: elevation")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----two_traits_multi---------------------------------------------------------
# ------ List of multiSimmaps: one element per trait ------ #

# Combine all maps for color with all maps for elevation in a list of two multiSimmap
maps_two_traits <- list(color = toy_multiSimmap, elevation = toy_multiSimmap_elevation)

# Number of maps per trait: must be the same
sapply(maps_two_traits, length)

# Graft the taxa using the map of each trait
grafted_two_traits_multi <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_two_traits,
                                         simmap = maps_two_traits, verbose = FALSE, seed = 1)

# We obtain 10 grafted phylogeny, one per map of each trait. 
length(grafted_two_traits_multi)
# Each grafted phylogeny contains one updated map per trait
names(grafted_two_traits_multi[[1]]$simmap)

## ----two_traits_duplicate-----------------------------------------------------
# ------ Pair a single simmap with a multiSimmap ------ #

# Duplicate the map 1 of the elevation 10 times (with [ ], the class multiSimmap is kept)
elevation_duplicated <- toy_multiSimmap_elevation[rep(1, length(toy_multiSimmap))]
sapply(list(color = toy_multiSimmap, elevation = elevation_duplicated), length)

## ----deepSTRAPP_ranges, eval = FALSE------------------------------------------
# # ------ Not run: build the maps of the ranges with deepSTRAPP and BioGeoBEARS ------ #
# 
# range_data <- deepSTRAPP::prepare_trait_data(tip_data = tip_ranges, # Named vector of ranges (or binary matrix)
#                                              trait_data_type = "biogeographic",
#                                              phylo = toy_backbone,
#                                              evolutionary_models = "DEC", # Models to compare, e.g., c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J")
#                                              max_range_size = 2, # Maximum number of areas in a range
#                                              BioGeoBEARS_directory_path = "./BioGeoBEARS_files/", # Directory for input/output files
#                                              keep_BioGeoBEARS_files = FALSE,
#                                              nb_cores = 1,
#                                              seed = 3456,
#                                              nb_simulations = 100,
#                                              return_simmaps = TRUE,
#                                              return_model_selection_df = TRUE,
#                                              return_BSM = FALSE, # TRUE to also get the tables of events
#                                              plot_map = FALSE, plot_overlay = FALSE)
# 
# # Comparison of the models (Here, we fit only a DEC model)
# range_data$model_selection_df
# 
# # The stochastic maps of ranges, converted from the BSM events: a multiSimmap of 100 maps
# my_multiSimmap_ranges <- range_data$simmaps

## ----ranges_check-------------------------------------------------------------
# ------ Maps of ranges ------ #

class(toy_multiSimmap_ranges)
length(toy_multiSimmap_ranges)

# Ranges of the tips
sort(unique(tip_ranges))

# Ranges visited in the maps, i.e., with a positive residence time on the branches
# (the columns of 'mapped.edge' list all the ranges of the model, including those never visited)
visited_ranges <- sort(unique(unlist(lapply(toy_multiSimmap_ranges, function (x) {
  colnames(x$mapped.edge)[colSums(x$mapped.edge) > 0] }))))
visited_ranges

## ----ranges_graft-------------------------------------------------------------
# ------ Graft with a multiSimmap of ranges ------ #

# Missing taxa with their allowed ranges
data(toy_missing_taxa_list_with_allowed_ranges)

# Graft the taxa using the first two maps of ranges
grafted_ranges <- graft_missing_taxa(phy = toy_backbone,
                                     missing_taxa_list = toy_missing_taxa_list_with_allowed_ranges,
                                     simmap = toy_multiSimmap_ranges[1:2],
                                     allowed_states = c("A", "B"), treat_states_as_ranges = TRUE,
                                     nb_simulations = 1, verbose = FALSE, seed = 1)
length(grafted_ranges)

## ----ranges_plot, fig.width = 10, fig.height = 10, out.width = "100%"---------
# ------ Initial maps (upper line) and grafted maps (lower line) ------ #

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
# Upper line: initial maps of ranges
for (i in 1:2) { plot_initial_map(toy_multiSimmap_ranges[[i]], colors = range_colors, main = paste("Initial map", i)) }
# Lower line: grafted maps
for (i in 1:2) { plot_updated_map(grafted_ranges[[i]], colors = range_colors, main = paste("Grafted map", i)) }
legend("bottomleft", legend = names(range_colors)[names(range_colors) %in% visited_ranges],
       fill = range_colors[names(range_colors) %in% visited_ranges], title = "Range", bty = "n", cex = 0.8)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----deepSTRAPP_continuous, eval = FALSE--------------------------------------
# # ------ Not run: build the continuous stochastic maps with deepSTRAPP ------ #
# 
# data(toy_tip_values) # Named vector of continuous values
# my_trait_data <- deepSTRAPP::prepare_trait_data(tip_data = toy_tip_values,
#                                                 trait_data_type = "continuous",
#                                                 phylo = toy_backbone,
#                                                 evolutionary_models = "BM", # Models to compare, e.g., c("BM", "OU", "EB")
#                                                 run_stochastic_maps = TRUE, # Simulate the continuous stochastic maps
#                                                 nb_simulations = 100,
#                                                 seed = 1234,
#                                                 plot_map = FALSE)
# 
# # Point estimate (interpolated ML estimates of the ancestral values)
# my_contMap <- my_trait_data$contMap
# 
# # List of 100 contMaps: the stochastic maps of the continuous trait
# my_contMaps_list <- my_trait_data$contMaps

## ----continuous_convert-------------------------------------------------------
# ------ Discretize the contMaps into a multiSimmap ------ #

data(toy_tip_values)
data(toy_contMaps)

state_names <- c("low", "medium", "high") # From the lowest to the highest values
state_colors <- c(low = "steelblue", medium = "grey75", high = "firebrick")

multiSimmap_continuous <- convert_contMap_to_simmap(contMap = toy_contMaps,
                                                    probs = c(1/3, 2/3), # Terciles of the tip values
                                                    states = state_names,
                                                    tip_values = toy_tip_values,
                                                    verbose = FALSE)
class(multiSimmap_continuous)
length(multiSimmap_continuous)

## ----continuous_graft---------------------------------------------------------
# ------ Graft with the discretized maps ------ #

data(toy_missing_taxa_list)

# Continuous values of the 7 missing taxa
missing_values <- c(New_sp1 = 0.5, New_sp2 = 2.5, New_sp3 = 1.4, New_sp4 = -3,
                    New_sp5 = -1, New_sp6 = 3.5, New_sp7 = 1.6)

# Convert into allowed states: provide the multiSimmap itself, as it holds the cut points
taxa_continuous <- toy_missing_taxa_list
taxa_continuous$allowed_states <- convert_values_to_states(values = missing_values,
                                                           simmap = multiSimmap_continuous)

# Graft the taxa using the first two maps
grafted_continuous <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_continuous,
                                         simmap = multiSimmap_continuous[1:2],
                                         nb_simulations = 1, verbose = FALSE, seed = 1)

## ----continuous_plot, fig.width = 10, fig.height = 10, out.width = "100%"-----
# ------ Initial maps (upper line) and grafted maps (lower line) ------ #

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
# Upper line: initial discretized maps
for (i in 1:2) { plot_initial_map(multiSimmap_continuous[[i]], colors = state_colors, main = paste("Initial map", i)) }
# Lower line: grafted maps
for (i in 1:2) { plot_updated_map(grafted_continuous[[i]], colors = state_colors, main = paste("Grafted map", i)) }
legend("bottomleft", legend = names(state_colors), fill = state_colors, title = "State", bty = "n", cex = 0.8)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----pitfall_error, error = TRUE----------------------------------------------
try({
# ------ Error with allowed states ------ #

# Load missing taxa list
data(toy_missing_taxa_list)

# A set of maps for a unique trait is a "multiSimmap", "multiPhylo" object
class(toy_multiSimmap)

# Convert multiSimmap into a plain list: class lost
plain_list <- lapply(toy_multiSimmap, identity)

# It becomes a list of simmaps with both classes: "phylo" and "simmap".
class(plain_list)
table(unlist(lapply(plain_list, FUN = class)))

# Check the grafting plan using a list of 10 simmaps, but allowed_states for a unique trait (color)
check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                    simmap = plain_list, allowed_states = c("blue", "green"), verbose = FALSE)
})

## ----pitfall_fix--------------------------------------------------------------
# ------ Restore the class ------ #

class(plain_list) <- c("multiSimmap", "multiPhylo")
plan_fixed <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                  simmap = plain_list, verbose = FALSE)
attr(plan_fixed, "n_maps")

## ----pitfall_n_maps, error = TRUE---------------------------------------------
try({
# ------ Different numbers of maps among traits ------ #

check_grafting_plan(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                    simmap = list(color = toy_multiSimmap, elevation = toy_multiSimmap_elevation[1:5]),
                    verbose = FALSE)
})

