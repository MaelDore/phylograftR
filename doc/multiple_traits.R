## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----data, fig.width = 10, fig.height = 5, out.width = "100%"-----------------
# ------ Toy data ------ #

data(toy_backbone)
data(toy_simmap)
data(toy_simmap_elevation)
data(toy_missing_taxa_list_with_two_traits)

# Group the maps of the traits in a named list: the names are the names of the traits
maps_color_elevation <- list(color = toy_simmap, elevation = toy_simmap_elevation)

color_colors <- c(blue = "blue", green = "green3", red = "red")
elevation_colors <- c(low = "steelblue", high = "darkorange")

# Plot the two maps side by side
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
phytools::plotSimmap(toy_simmap, colors = color_colors, fsize = 0.55, ftype = "i", lwd = 2)
title("A: Color", line = -1, cex.main = 1)
phytools::plotSimmap(toy_simmap_elevation, colors = elevation_colors, fsize = 0.55, ftype = "i", lwd = 2)
title("B: Elevation", line = -1, cex.main = 1)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----list---------------------------------------------------------------------
# ------ Allowed states of two traits for each taxon ------ #

# Rename the loaded missing_taxa_list for simplicity
taxa_two_traits <- toy_missing_taxa_list_with_two_traits

# Helper: extract the taxon-level allowed states of one trait
# Elements of the list are either named after the traits, or in the order of the maps
get_allowed <- function(taxa, trait_name, trait_index)
{
  vapply(taxa$allowed_states, function(x)
  {
    values <- if (is.null(names(x))) x[trait_index][[1]] else x[[trait_name]]
    if (is.null(values)) "NULL" else paste(values, collapse = ", ")
  }, character(1))
}

# Build and print allowed states as a data.frame
data.frame(Taxon_name = taxa_two_traits$Taxon_name,
           Color = get_allowed(taxa_two_traits, "color", 1),
           Elevation = get_allowed(taxa_two_traits, "elevation", 2),
           row.names = NULL)

# NULL will be using the global 'allowed_states' provided to the grafting function


## ----csv----------------------------------------------------------------------
# The initial table of the toy list
csv_file <- system.file("extdata", "toy_missing_taxa_with_two_traits.csv", package = "phylograftR")
csv_table <- read.csv(csv_file)
csv_table[, c("Taxon_name", "allowed_states")]

## ----joint_list---------------------------------------------------------------
# ------ Add elevation constraints to three taxa ------ #

taxa_two_traits$allowed_states$New_sp1 <- list(color = "red", elevation = "high")
taxa_two_traits$allowed_states$New_sp3 <- list(color = NULL, elevation = "low") # NULL: the global colors are used
taxa_two_traits$allowed_states$New_sp4 <- list(color = "blue", elevation = "high")

## ----joint_plan---------------------------------------------------------------
# ------ Valid fraction with the color only, the elevation only, and both ------ #

# Helper: keep the taxon-level allowed states of one trait only
keep_trait <- function(taxa, trait, traits = c("color", "elevation"))
{
  taxa$allowed_states <- lapply(taxa$allowed_states, function(x)
  {
    if (is.null(x)) return(NULL)
    out <- setNames(vector("list", length(traits)), traits)
    for (i in seq_along(x))
    {
      name <- if (is.null(names(x))) traits[i] else names(x)[i]
      out[name] <- list(x[[i]])
    }
    out[setdiff(traits, trait)] <- list(NULL)
    if (all(vapply(out, is.null, logical(1)))) NULL else out
  })
  taxa
}

# Set the global allowed states: blue or green for the color, no restriction for the elevation
global_allowed_states <- list(c("blue", "green"), NULL)

# Check grafting plan with colors only
plan_colors <- check_grafting_plan(phy = toy_backbone, 
                                   missing_taxa_list = keep_trait(taxa_two_traits, "color"),
                                   simmap = maps_color_elevation,
                                   allowed_states = global_allowed_states, verbose = FALSE)
# Check grafting plan with elevation only
plan_elevation <- check_grafting_plan(phy = toy_backbone, 
                                      missing_taxa_list = keep_trait(taxa_two_traits, "elevation"),
                                   simmap = maps_color_elevation,
                                   allowed_states = list(NULL, NULL), verbose = FALSE)
# Check grafting plan with both
plan_both <- check_grafting_plan(phy = toy_backbone, 
                                 missing_taxa_list = taxa_two_traits,
                                 simmap = maps_color_elevation,
                                 allowed_states = global_allowed_states, verbose = FALSE)

# Extract the valid fractions for each plan
valid_fractions <- data.frame(Taxon_name = taxa_two_traits$Taxon_name,
  Color_only = plan_colors$valid_fraction,
  Elevation_only = plan_elevation$valid_fraction,
  Both = plan_both$valid_fraction)
valid_fractions[, -1] <- round(valid_fractions[, -1], 2)
valid_fractions


## ----joint_plan_flag----------------------------------------------------------
# The plan flags New_sp4 for its low valid fraction
plan_both[plan_both$Taxon_name == "New_sp4", c("Taxon_name", "valid_fraction", "feasible", "cause")]

# Inspect overall feasibility
plan_both[, c("Taxon_name", "feasible")]


## ----graft_two_traits---------------------------------------------------------
# ------ Graft with two categorical traits ------ #

grafted_two_traits <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_two_traits,
                                         simmap = maps_color_elevation,
                                         allowed_states = list(c("blue", "green"), NULL),
                                         verbose = FALSE, seed = 4)

# The output has one updated simmap per trait, named after the traits
names(grafted_two_traits$simmap)

## ----plot_two_traits, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# ------ Compare the updated maps of the two traits ------ #

# Helper: plot the updated simmap of one trait of a grafted phylogeny
plot_trait <- function(phy_grafted, trait, colors, main)
{
  plot_grafted(phy_grafted, which = "edges", trait = trait, simmap_colors = colors, fsize = 0.55, main = main)
  grafted_tips <- which(phy_grafted$grafted_tips)
  tip_states <- phytools::getStates(phy_grafted$simmap[[trait]], type = "tips")[phy_grafted$tip.label[grafted_tips]]
  ape::tiplabels(tip = grafted_tips, pch = 19, cex = 1.5, col = colors[tip_states])
}

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot_trait(grafted_two_traits, trait = "color", colors = color_colors, main = "A: Color")
plot_trait(grafted_two_traits, trait = "elevation", colors = elevation_colors, main = "B: Elevation")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----info_two_traits----------------------------------------------------------
# ------ States of the grafted taxa ------ #

info_two_traits <- extract_grafting_info(grafted_two_traits, missing_taxa_list = taxa_two_traits)

# Helper: allowed states of a trait for each taxon, taxon-level value or else global value
allowed_of_trait <- function(taxa, trait_name, trait_index, global)
{
  values <- get_allowed(taxa, trait_name, trait_index)
  ifelse(values == "NULL", if (is.null(global)) "all" else paste(global, collapse = ", "), values)
}

data.frame(Taxon_name = info_two_traits$Taxon_name,
           Allowed_color = allowed_of_trait(taxa_two_traits, "color", 1, c("blue", "green")),
           Color = info_two_traits$tip_state_color,
           Allowed_elevation = allowed_of_trait(taxa_two_traits, "elevation", 2, NULL),
           Elevation = info_two_traits$tip_state_elevation,
           row.names = NULL)

## ----check_two_traits---------------------------------------------------------
# ------ A posteriori check ------ #

constraints_two_traits <- check_grafting_constraints(phy_grafted = grafted_two_traits, phy_backbone = toy_backbone,
                                                     missing_taxa_list = taxa_two_traits,
                                                     allowed_states = list(c("blue", "green"), NULL),
                                                     use_excluded_subclades = FALSE, verbose = FALSE)
constraints_two_traits[, c("Taxon_name", "states_ok_color", "states_ok_elevation", "states_ok")]

## ----data_ranges--------------------------------------------------------------
# ------ Toy data: ranges and elevation ------ #

data(toy_simmap_ranges)
data(toy_missing_taxa_list_with_ranges_and_elevation)

# Rename the missing taxa list for consistency
taxa_ranges_elevation <- toy_missing_taxa_list_with_ranges_and_elevation
# Combine simmap for biogeographic ranges and elevation
maps_ranges_elevation <- list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation)

# Print in a data.frame the taxon-level allowed states of the two traits
data.frame(Taxon_name = taxa_ranges_elevation$Taxon_name,
           Ranges = get_allowed(taxa_ranges_elevation, "biogeo", 1),
           Elevation = get_allowed(taxa_ranges_elevation, "elevation", 2),
           row.names = NULL)

# NULL will be using the global 'allowed_states' provided to the grafting function


## ----graft_ranges_elevation---------------------------------------------------
# ------ Graft with ranges and elevation ------ #

grafted_ranges_elevation <- graft_missing_taxa(phy = toy_backbone, 
                                               missing_taxa_list = taxa_ranges_elevation,
                                               simmap = maps_ranges_elevation,
                                               allowed_states = list(c("A", "B"), NULL),
                                               treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
                                               verbose = FALSE, seed = 9)

## ----plot_ranges_elevation, fig.width = 10, fig.height = 5.5, out.width = "100%"----
# ------ Compare the updated maps of the two traits ------ #

range_colors <- c(A = "gold", B = "dodgerblue", AB = "limegreen", BC = "firebrick2")

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot_trait(grafted_ranges_elevation, trait = "biogeo", colors = range_colors, main = "A: Ranges")
plot_trait(grafted_ranges_elevation, trait = "elevation", colors = elevation_colors, main = "B: Elevation")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----info_ranges_elevation----------------------------------------------------
# ------ States at grafting and inherited states ------ #

info_ranges_elevation <- extract_grafting_info(grafted_ranges_elevation, missing_taxa_list = taxa_ranges_elevation)
info_ranges_elevation[, c("Taxon_name", "parent_state_biogeo", "tip_state_biogeo",
                          "parent_state_elevation", "tip_state_elevation")]

## ----check_ranges_elevation---------------------------------------------------
# ------ A posteriori check ------ #

constraints_ranges_elevation <- check_grafting_constraints(phy_grafted = grafted_ranges_elevation, phy_backbone = toy_backbone,
                                                           missing_taxa_list = taxa_ranges_elevation,
                                                           allowed_states = list(c("A", "B"), NULL),
                                                           treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
                                                           use_excluded_subclades = FALSE, verbose = FALSE)
constraints_ranges_elevation[, c("Taxon_name", "states_ok_biogeo", "states_ok_elevation", "ranges_ok", "ok")]

## ----multisimmap_maps---------------------------------------------------------
# ------ 10 stochastic maps for the color and for the elevation ------ #

# States of the tips of each trait
tip_colors <- phytools::getStates(toy_simmap, type = "tips")
tip_elevations <- phytools::getStates(toy_simmap_elevation, type = "tips")

# Simulate 10 stochastic maps per trait
set.seed(33)
multi_color <- phytools::make.simmap(toy_backbone, tip_colors, model = "ER", nsim = 10, message = FALSE)
multi_elevation <- phytools::make.simmap(toy_backbone, tip_elevations, model = "ER", nsim = 10, message = FALSE)

# Group them in a named list, as with a single map per trait
multi_maps <- list(color = multi_color, elevation = multi_elevation)

## ----multisimmap_plan---------------------------------------------------------
# ------ Grafting plan across the 10 pairs of maps ------ #

plan_paired <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_two_traits,
                                   simmap = multi_maps,
                                   allowed_states = list(c("blue", "green"), NULL), verbose = FALSE)
plan_paired[, c("Taxon_name", "valid_fraction_min", "valid_fraction_mean", "valid_fraction_max", "n_maps_infeasible", "feasible")]

## ----multisimmap_graft--------------------------------------------------------
# ------ Graft with paired maps ------ #

grafted_paired <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_two_traits,
                                     simmap = multi_maps,
                                     allowed_states = list(c("blue", "green"), NULL),
                                     method = "uniform", min_valid_fraction = 0,
                                     nb_simulations = 5, seed = 1, verbose = FALSE)

## ----multisimmap_output-------------------------------------------------------
# ------ Structure of the output ------ #

# One element per pair of maps
length(grafted_paired)
names(grafted_paired)

# Each element contains 5 grafted phylogenies
sapply(grafted_paired, length)

# Total number of grafted phylogenies
sum(sapply(grafted_paired, length))

## ----multisimmap_simmaps------------------------------------------------------
# ------ Updated simmaps of a grafted phylogeny ------ #

# First phylogeny of the first pair of maps
phy_1_1 <- grafted_paired$map_1[[1]]
names(phy_1_1$simmap)
class(phy_1_1$simmap$color)
class(phy_1_1$simmap$elevation)

## ----multisimmap_cross, eval = FALSE------------------------------------------
# # ------ Not run: cross all the pairs of maps ------ #
# 
# # Each map of the color repeated 10 times in a row: c1, c1, ..., c2, c2, ..., c10
# crossed_color <- multi_color[rep(1:10, each = 10)]
# # The 10 maps of the elevation repeated 10 times: e1, ..., e10, e1, ..., e10, ...
# crossed_elevation <- multi_elevation[rep(1:10, times = 10)]
# 
# # 100 pairs of maps = all the combinations of the 10 maps of each trait
# grafted_crossed <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_two_traits,
#                                       simmap = list(color = crossed_color, elevation = crossed_elevation),
#                                       allowed_states = list(c("blue", "green"), NULL),
#                                       nb_simulations = 5, seed = 1)

