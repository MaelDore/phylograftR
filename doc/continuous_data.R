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
data(toy_tip_values)
data(toy_contMap)
data(toy_contMaps)
data(toy_missing_taxa_list)

# Provide the continuous trait values of the 7 missing taxa
missing_values <- c(New_sp1 = 0.5, New_sp2 = 2.5, New_sp3 = 1.4, New_sp4 = -3,
                    New_sp5 = -1, New_sp6 = 3.5, New_sp7 = 1.6)

## ----colors-------------------------------------------------------------------
# ------ Names and colors of the states ------ #

state_names <- c("low", "medium", "high") # From the lowest to the highest values
state_colors <- c(low = "steelblue", medium = "grey75", high = "firebrick")

# ------ Build a continuous color scale matching the colors of the states ------ #

terciles <- quantile(toy_tip_values, probs = c(1/3, 2/3)) # Cut points used in Section 1 below
scale_colors <- colorRampPalette(state_colors)(101) # Continuous scale from blue to grey to red

# To match the discrete scale, we want each color to match with the median values of each group
# such as it partly adjusts to the distribution of trait values inside their range
#  - "grey" is placed at the median of the "medium" state
#  - "steelblue" is placed in the median of the "low" state
#  - "firebrick" is placed in the median of the "high" state

# Build a regular sequence of trait values within their range
trait_values <- seq(toy_contMap$lims[1], toy_contMap$lims[2], length.out = length(toy_contMap$cols))
# Extract the median trait values for each group
anchors <- c(median(toy_tip_values[toy_tip_values < terciles[1]]),
             median(toy_tip_values[toy_tip_values >= terciles[1] & toy_tip_values < terciles[2]]),
             median(toy_tip_values[toy_tip_values >= terciles[2]]))
# Map the index of the color scale to fit the distribution of trait values
mapped_indices <- round(approx(x = anchors, y = c(1, 51, 101), xout = trait_values, rule = 2)$y)

# Update color scheme in the contMap
toy_contMap_colored <- toy_contMap
toy_contMap_colored$cols[] <- scale_colors[mapped_indices]


## ----convert_single-----------------------------------------------------------
# ------ Convert one contMap into a simmap ------ #

simmap_terciles <- convert_contMap_to_simmap(contMap = toy_contMap,
                                             probs = c(1/3, 2/3), # Terciles of the tip values
                                             states = state_names, # From the lowest to the highest values
                                             tip_values = toy_tip_values, # Observed tip values
                                             verbose = FALSE)

# The cut points are stored in the output
attr(simmap_terciles, "breaks")

# Number of tips per state: is equal by design
table(factor(phytools::getStates(simmap_terciles, type = "tips"), levels = state_names))

## ----convert_single_plot, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# ------ Compare the contMap and the simmap ------ #

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
# plot(toy_contMap_colored, legend = FALSE, fsize = c(0.55, 0.8), lwd = 2, outline = FALSE, mar = c(0.5, 0.5, 2, 0.5))
plot(toy_contMap_colored, legend = FALSE, fsize = c(0.55, 0.8), lwd = 3, outline = FALSE)
title("A: contMap", line = -1, cex.main = 1)
phytools::plotSimmap(simmap_terciles, colors = state_colors, fsize = 0.55, ftype = "i", lwd = 3)
title("B: simmap with three states", line = -1, cex.main = 1)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----allowed_states-----------------------------------------------------------
# ------ Convert the values of the missing taxa into allowed states ------ #

# Use the same cut points as for the simmap (the output is a named list, as expected for the allowed states)
taxon_states <- convert_values_to_states(values = missing_values, simmap = simmap_terciles)

# Use the states as the allowed states of each taxon in the missing_taxa_list
taxa_continuous <- toy_missing_taxa_list
taxa_continuous$allowed_states <- taxon_states

# Print as a data.frame for inspection
data.frame(Taxon_name = names(missing_values), Value = missing_values, Allowed_state = unlist(taxon_states), row.names = NULL)

## ----graft_single-------------------------------------------------------------
# ------ Graft with the simmap as constraint ------ #

# Graft the missing taxa
grafted_single <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_continuous,
                                     simmap = simmap_terciles, verbose = FALSE, seed = 1)

# State at the point of grafting (parent branch) and inherited by each taxon
info_single <- extract_grafting_info(grafted_single, missing_taxa_list = taxa_continuous)
data.frame(info_single[, c("Taxon_name", "parent_state", "tip_state")], Allowed_state = unlist(taxon_states), row.names = NULL)

## ----thresholds---------------------------------------------------------------
# ------ Compare sets of thresholds ------ #

# Arguments defining the cut points of each set
threshold_sets <- list(
  Terciles = list(probs = c(1/3, 2/3)),
  Quartiles = list(probs = c(0.25, 0.75)),
  Equal_width = list(n_bins = 3),
  Terciles_up = list(breaks = attr(simmap_terciles, "breaks") + 0.5))

# Helper: discretize the contMap with one set of thresholds,
# and convert the values of the missing taxa into allowed states
convert_with_thresholds <- function (threshold_args)
{
  simmap <- do.call(convert_contMap_to_simmap,
                    c(list(contMap = toy_contMap, states = state_names, tip_values = toy_tip_values, verbose = FALSE), threshold_args))
  list(simmap = simmap, states = convert_values_to_states(values = missing_values, simmap = simmap))
}
threshold_results <- lapply(threshold_sets, convert_with_thresholds)

# Extract the cut points of each set
sapply(threshold_results, function(x) round(attr(x$simmap, "breaks"), 2))

## ----thresholds_table---------------------------------------------------------
# ------ Allowed state per taxon and set of thresholds ------ #

data.frame(Value = missing_values, sapply(threshold_results, function(x) unlist(x$states)))

## ----thresholds_plot, fig.width = 10, fig.height = 9, out.width = "100%"------
# ------ Compare the simmaps ------ #

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
for (set_name in names(threshold_results))
{
  phytools::plotSimmap(threshold_results[[set_name]]$simmap, colors = state_colors, fsize = 0.55, ftype = "i", lwd = 3)
  title(gsub("_", " ", set_name), line = -1, cex.main = 1)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----thresholds_graft, fig.width = 10, fig.height = 10, out.width = "100%"----
# ------ Graft with each set of thresholds ------ #

# Helper: set allowed states in a list of missing taxa
set_allowed_states <- function (states)
{
  taxa <- toy_missing_taxa_list
  taxa$allowed_states <- states
  taxa
}

# Graft and plot the updated simmap of the grafted phylogeny
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
for (set_name in names(threshold_results))
{
  grafted_set <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = set_allowed_states(threshold_results[[set_name]]$states),
                                    simmap = threshold_results[[set_name]]$simmap, verbose = FALSE, seed = 1)
  plot_grafted(grafted_set, simmap_colors = state_colors, color = "purple",
               fsize = 0.55, lwd = 2, main = gsub("_", " ", set_name))
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----plan_single--------------------------------------------------------------
# ------ Check the grafting plan ------ #

# Check the initial grafting plan, without tolerance in the allowed states
plan_single <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_continuous,
                                   simmap = simmap_terciles, verbose = FALSE)
plan_single[, c("Taxon_name", "allowed_states", "valid_fraction", "feasible")]

## ----tolerance----------------------------------------------------------------
# ------ Allowed states with a tolerance around the value ------ #

# Convert the values into all the states covered by [value - 0.3, value + 0.3]
allowed_with_tolerance <- convert_values_to_states(values = missing_values, simmap = simmap_terciles, tolerance = 0.3)
allowed_with_tolerance

# Update the missing_taxa_list with allowed states based on tolerance
taxa_tolerance <- toy_missing_taxa_list
taxa_tolerance$allowed_states <- allowed_with_tolerance

# Check the new grafting plan
plan_tolerance <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_tolerance,
                                      simmap = simmap_terciles, verbose = FALSE)

# Compare valid fractions between the two grafting plans
data.frame(Taxon_name = plan_single$Taxon_name,
           Allowed_state = plan_single$allowed_states, Valid_fraction = round(plan_single$valid_fraction, 2),
           Allowed_with_tolerance = plan_tolerance$allowed_states, Valid_fraction_with_tolerance = round(plan_tolerance$valid_fraction, 2))

## ----contsimmap, eval = FALSE-------------------------------------------------
# # ------ Not run: build the list of contMaps with contsimmap ------ #
# 
# # Estimate the evolutionary rate (sigma²) of a BM
# BM_fit <- phytools::anc.ML(tree = toy_backbone, x = toy_tip_values, model = "BM")
# sigma2 <- BM_fit$sig2
# 
# # Simulate 100 continuous stochastic maps conditioned on the tip data, the phylogeny and the BM
# my_contsimmap <- contsimmap::make.contsimmap(tree = toy_backbone, trait.data = toy_tip_values,
#                                              Xsig2 = list(sigma2), nsim = 100)
# 
# # Convert into a list of contMaps
# my_contMaps_list <- deepSTRAPP::convert_contsimmap_to_contMaps(contsimmap = my_contsimmap)
# 
# # Keep only the first 10 maps in the shipped dataset `toy_contMaps`, to keep the package light:
# # 10 maps are enough for illustration, but use all your maps (100 or more) in actual analyses
# my_contMaps_list <- my_contMaps_list[1:10]
# 
# # Alternative (recommended): all the steps above in one call with deepSTRAPP
# my_trait_data <- deepSTRAPP::prepare_trait_data(tip_data = toy_tip_values,
#                                                 trait_data_type = "continuous",
#                                                 phylo = toy_backbone, seed = 1234,
#                                                 evolutionary_models = "BM",
#                                                 run_stochastic_maps = TRUE,
#                                                 nb_simulations = 100, plot_map = FALSE)
# 
# my_contMap <- my_trait_data$contMap # Point estimate of the interpolated ML estimates
# my_contMaps_list <- my_trait_data$contMaps # List of contMaps from the stochastic maps

## ----multisimmap_convert------------------------------------------------------
# ------ Discretize the list of contMaps into a multiSimmap ------ #

# Use the shipped list of continuous stochastic maps loaded in Section 0 (10 maps)

# Discretize all contMaps at once using the same cut points
multiSimmap_continuous <- convert_contMap_to_simmap(contMap = toy_contMaps,
                                                    probs = c(1/3, 2/3), # Terciles of the tip values
                                                    states = state_names,
                                                    tip_values = toy_tip_values, # Observed tip values
                                                    verbose = FALSE)
class(multiSimmap_continuous)
length(multiSimmap_continuous)

## ----multisimmap_plot, fig.width = 10, fig.height = 4.5, out.width = "100%"----
# ------ Plot three histories ------ #

par(mfrow = c(1, 3), mar = c(0.5, 0.5, 2, 0.5))
for (i in 1:3)
{
  phytools::plotSimmap(multiSimmap_continuous[[i]], colors = state_colors, fsize = 0.5, ftype = "i", lwd = 2)
  title(paste("Map", i), line = -1, cex.main = 1)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----multisimmap_plan---------------------------------------------------------
# ------ Grafting plan across the maps ------ #

plan_multi <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_continuous,
                                  simmap = multiSimmap_continuous, verbose = FALSE)
plan_multi[, c("Taxon_name", "allowed_states", "valid_fraction_min", "valid_fraction_mean", "valid_fraction_max", "n_maps_infeasible")]

## ----multisimmap_low----------------------------------------------------------
# ------ Number of maps with a low valid fraction ------ #

plan_per_map <- attr(plan_multi, "per_map") # Plan of each map
tapply(plan_per_map$valid_fraction < 0.1, plan_per_map$Taxon_name, sum)

## ----multisimmap_graft--------------------------------------------------------
# ------ Graft with the multiSimmap ------ #

grafted_multi <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_continuous,
                                    simmap = multiSimmap_continuous,
                                    min_valid_fraction = 0, nb_simulations = 1,
                                    seed = 1, verbose = FALSE)

# One element per map
length(grafted_multi)

## ----multisimmap_info---------------------------------------------------------
# ------ Position of the grafts across the maps ------ #

info_multi <- extract_grafting_info(grafted_multi, missing_taxa_list = taxa_continuous, phy_backbone = toy_backbone)

# Number of backbone branches used for the grafting, and range of grafting ages across the 10 maps
data.frame(Allowed_state = unlist(taxon_states),
           N_parent_branches = tapply(info_multi$parent_node, info_multi$Taxon_name, function(x) length(unique(x))),
           Graft_age_min = round(tapply(info_multi$graft_age, info_multi$Taxon_name, min), 2),
           Graft_age_median = round(tapply(info_multi$graft_age, info_multi$Taxon_name, median), 2),
           Graft_age_max = round(tapply(info_multi$graft_age, info_multi$Taxon_name, max), 2))

## ----multisimmap_ages, fig.width = 8, fig.height = 4.5, out.width = "80%"-----
# ------ Distribution of the grafting ages across the maps ------ #

par(mar = c(4.5, 4.5, 1, 1))
boxplot(graft_age ~ Taxon_name, data = info_multi, ylab = "Grafting age", xlab = "", las = 2,
        col = state_colors[unlist(taxon_states)[levels(factor(info_multi$Taxon_name))]])
legend("topright", legend = state_names, fill = state_colors[state_names], title = "Allowed state", bty = "n")
par(mar = c(5.1, 4.1, 4.1, 2.1))

## ----transitions--------------------------------------------------------------
# ------ Number of changes of state per map ------ #

# Number of changes of state per map (number of segments - number of branches)
n_changes <- sapply(multiSimmap_continuous, function(x) sum(lengths(x$maps) - 1))
summary(n_changes)

# Compare with the contMap (based on a monotonic interpolation of the values along each branch)
sum(lengths(simmap_terciles$maps) - 1)

## ----make_simmap--------------------------------------------------------------
# ------ Discretize the tips, then build simmaps ------ #

# Discretize the observed tip values with the same cut points as in Section 1
tip_states <- unlist(convert_values_to_states(values = toy_tip_values, simmap = simmap_terciles))
tip_states <- factor(tip_states, levels = state_names) # Order of the states

# Ordered model: transitions only between neighbor states (low <-> medium <-> high) with symmetrical rates
ordered_model <- matrix(c(0, 1, 0,
                          1, 0, 2,
                          0, 2, 0), nrow = 3, byrow = TRUE, dimnames = list(state_names, state_names))

# Simulate 10 stochastic maps (enough for illustration, but use more for actual analyses)
set.seed(1)
simmaps_discrete <- phytools::make.simmap(toy_backbone, tip_states, model = ordered_model, nsim = 10, message = FALSE)

# Number of changes of state per map
summary(sapply(simmaps_discrete, function(x) sum(lengths(x$maps) - 1)))

## ----make_simmap_plot, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# ------ Compare the simmaps ------ #

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
phytools::plotSimmap(multiSimmap_continuous[[1]], colors = state_colors, fsize = 0.55, ftype = "i", lwd = 2)
title("A: continuous stochastic map, discretized", line = -1, cex.main = 1)
phytools::plotSimmap(simmaps_discrete[[1]], colors = state_colors, fsize = 0.55, ftype = "i", lwd = 2)
title("B: make.simmap on discretized tips", line = -1, cex.main = 1)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----make_simmap_nodes, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# ------ Probability of the states at the nodes ------ #

nodes <- ape::Ntip(toy_backbone) + seq_len(toy_backbone$Nnode)
node_probs_continuous <- phytools::describe.simmap(multiSimmap_continuous, plot = FALSE)$ace[as.character(nodes), state_names]
node_probs_discrete <- phytools::describe.simmap(simmaps_discrete, plot = FALSE)$ace[as.character(nodes), state_names]

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
for (panel in 1:2)
{
  phytools::plotTree(toy_backbone, ftype = "i", fsize = 0.55, lwd = 1, mar = c(0.5, 0.5, 2, 0.5))
  ape::nodelabels(pie = list(node_probs_continuous, node_probs_discrete)[[panel]], piecol = state_colors[state_names], cex = 0.6)
  title(c("A: continuous stochastic maps, discretized", "B: make.simmap on discretized tips")[panel], line = -1, cex.main = 1)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----make_simmap_compare------------------------------------------------------
# ------ Compatibility with the values of the contMap ------ #

# Ancestral values of the contMap (ML under BM), and their state
node_values <- phytools::fastAnc(toy_backbone, toy_tip_values)[as.character(nodes)]
node_states <- unlist(convert_values_to_states(values = node_values, simmap = simmap_terciles))

# Probability of the state of the reference value at each node: 
# how compatible is the reconstruction with the ML estimates of the continuous values?
prob_reference_continuous <- node_probs_continuous[cbind(seq_along(nodes), match(node_states, state_names))]
prob_reference_discrete <- node_probs_discrete[cbind(seq_along(nodes), match(node_states, state_names))]

# Probability of being two states away from the reference state (low vs high):
# how often the reconstruction strongly disagree with the reference?
prob_opposite <- function (node_probs)
{
  opposite <- c(low = "high", medium = NA, high = "low")[node_states]
  sum(node_probs[cbind(seq_along(nodes), match(opposite, state_names))], na.rm = TRUE) / length(nodes)
}

# Summarize comparison in a data.frame
data.frame(Approach = c("Cont. stochastic maps, discretized", "make.simmap on discretized tips"),
           Mean_prob_reference_state = round(c(mean(prob_reference_continuous),
                                               mean(prob_reference_discrete)), 2),
           Mean_prob_two_states_away = round(c(prob_opposite(node_probs_continuous),
                                               prob_opposite(node_probs_discrete)), 2))


