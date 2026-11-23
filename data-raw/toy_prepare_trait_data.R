##### Build the toy datasets used to illustrate the preparation of stochastic maps from tip data #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Requires the R packages 'deepSTRAPP' (on CRAN) and 'BioGeoBEARS' (not on CRAN: remotes::install_github("nmatzke/BioGeoBEARS")).
# Run after 'toy_data.R' (it needs 'toy_backbone', 'toy_simmap', 'toy_simmap_elevation' and 'toy_simmap_ranges').
# Produces:
#  * data/toy_tip_data.rda: data.frame with the states of the 30 tips of 'toy_backbone' for three attributes: color, elevation and range.
#      They are the states at the tips of the toy simmaps, i.e., the 'tip data' from which stochastic maps are inferred in the vignette 'prepare_stochastic_maps'.
#  * data/toy_multiSimmap.rda: multiSimmap of 10 stochastic maps of the color, inferred from the tip data with deepSTRAPP::prepare_trait_data().
#  * data/toy_multiSimmap_elevation.rda: same, for the elevation.
#  * data/toy_multiSimmap_ranges.rda: multiSimmap of 10 stochastic maps of the ranges, inferred with BioGeoBEARS (BSM) through deepSTRAPP::prepare_trait_data().
#
# 100 maps are simulated for each attribute, as recommended for actual analyses (deepSTRAPP warns when 'nb_simulations' is lower than 100),
# and only the first 10 are kept in the package, to keep it light. 10 maps are enough for illustration.
# The tip data are always the states at the tips of the existing toy simmaps, so the new maps are consistent with the other toy datasets.


##### 1/ Tip data #####

data(toy_backbone)
data(toy_simmap)
data(toy_simmap_elevation)
data(toy_simmap_ranges)

tips <- toy_backbone$tip.label
toy_tip_data <- data.frame(color = phytools::getStates(toy_simmap, type = "tips")[tips],
                           elevation = phytools::getStates(toy_simmap_elevation, type = "tips")[tips],
                           range = phytools::getStates(toy_simmap_ranges, type = "tips")[tips],
                           row.names = tips,
                           stringsAsFactors = FALSE)
stopifnot(identical(row.names(toy_tip_data), tips), !anyNA(toy_tip_data))
print(lapply(toy_tip_data, table))


##### 2/ Helper to check and prepare the multiSimmap #####

# Keep the first 'n_keep' maps of the output of deepSTRAPP, with checks.
# '[' keeps the classes c("multiSimmap", "multiPhylo"), contrary to lapply() or unclass().
keep_first_maps <- function (maps, tip_states, n_keep = 10)
{
  stopifnot(inherits(maps, "multiSimmap"), length(maps) >= n_keep)
  maps <- maps[seq_len(n_keep)]
  names(maps) <- NULL # Unnamed maps ("map_1", "map_2", ... are used as names by the grafting functions)
  # Same tree as the backbone, and tips in the states of the tip data (maps are conditioned on the tip data)
  for (i in seq_along(maps))
  {
    stopifnot(identical(maps[[i]]$edge, toy_backbone$edge),
              isTRUE(all.equal(maps[[i]]$edge.length, toy_backbone$edge.length)),
              identical(unname(phytools::getStates(maps[[i]], type = "tips")[names(tip_states)]), unname(tip_states)))
    phylograftR:::validate_simmap(maps[[i]], toy_backbone)
  }
  stopifnot(inherits(maps, "multiSimmap"), length(maps) == n_keep)
  print(object.size(maps), units = "Kb")
  return(maps)
}


##### 3/ Color: categorical trait with three states #####

tip_colors <- setNames(toy_tip_data$color, tips)

# Compare the equal-rates (ER), symmetrical (SYM) and all-rates-different (ARD) models, and keep the best one (lowest AICc)
color_data <- deepSTRAPP::prepare_trait_data(tip_data = tip_colors,
                                             trait_data_type = "categorical",
                                             phylo = toy_backbone,
                                             evolutionary_models = c("ER", "SYM", "ARD"),
                                             seed = 1234,
                                             nb_simulations = 100,
                                             return_simmaps = TRUE,
                                             return_model_selection_df = TRUE,
                                             plot_map = FALSE, plot_overlay = FALSE,
                                             verbose = TRUE)
print(color_data$model_selection_df) # Model comparison
toy_multiSimmap <- keep_first_maps(color_data$simmaps, tip_states = tip_colors)
print(table(sapply(toy_multiSimmap, function (x) sum(lengths(x$maps) - 1)))) # Number of changes per map


##### 4/ Elevation: categorical trait with two states #####

tip_elevations <- setNames(toy_tip_data$elevation, tips)

elevation_data <- deepSTRAPP::prepare_trait_data(tip_data = tip_elevations,
                                                 trait_data_type = "categorical",
                                                 phylo = toy_backbone,
                                                 evolutionary_models = c("ER", "ARD"),
                                                 seed = 2345,
                                                 nb_simulations = 100,
                                                 return_simmaps = TRUE,
                                                 return_model_selection_df = TRUE,
                                                 plot_map = FALSE, plot_overlay = FALSE,
                                                 verbose = TRUE)
print(elevation_data$model_selection_df)
toy_multiSimmap_elevation <- keep_first_maps(elevation_data$simmaps, tip_states = tip_elevations)


##### 5/ Ranges: biogeographic ranges inferred with BioGeoBEARS #####

# Ranges are coded as in BioGeoBEARS: one capital letter per area, and multi-area ranges with letters in alphabetical order.
tip_ranges <- setNames(toy_tip_data$range, tips)

# BioGeoBEARS runs in a directory that stores its input and output files (removed afterwards with keep_BioGeoBEARS_files = FALSE)
# The best model is selected among the models compared (AICc). Here, only the DEC model is used, as in the default of prepare_trait_data()
# ('max_range_size = 2': ranges include at most two areas).
range_data <- deepSTRAPP::prepare_trait_data(tip_data = tip_ranges,
                                             trait_data_type = "biogeographic",
                                             phylo = toy_backbone,
                                             evolutionary_models = "DEC",
                                             BioGeoBEARS_directory_path = file.path(tempdir(), "BioGeoBEARS_toy_ranges"),
                                             keep_BioGeoBEARS_files = FALSE,
                                             max_range_size = 2,
                                             seed = 3456,
                                             nb_simulations = 100,
                                             return_simmaps = TRUE,
                                             return_model_selection_df = TRUE,
                                             plot_map = FALSE, plot_overlay = FALSE,
                                             verbose = TRUE)
print(range_data$model_selection_df)
toy_multiSimmap_ranges <- keep_first_maps(range_data$simmaps, tip_states = tip_ranges)

# Ranges found in the maps: the ranges of the tips, plus ranges only found at internal nodes (e.g., "C" or "AC")
print(colnames(toy_multiSimmap_ranges[[1]]$mapped.edge))
print(sort(unique(unlist(lapply(toy_multiSimmap_ranges, function (x) colnames(x$mapped.edge))))))


##### 6/ Save objects #####

dir.create("data", showWarnings = FALSE)
save(toy_tip_data, file = "data/toy_tip_data.rda", compress = "gzip", version = 2)
save(toy_multiSimmap, file = "data/toy_multiSimmap.rda", compress = "xz", version = 2)
save(toy_multiSimmap_elevation, file = "data/toy_multiSimmap_elevation.rda", compress = "xz", version = 2)
save(toy_multiSimmap_ranges, file = "data/toy_multiSimmap_ranges.rda", compress = "xz", version = 2)
print(file.size(c("data/toy_tip_data.rda", "data/toy_multiSimmap.rda", "data/toy_multiSimmap_elevation.rda", "data/toy_multiSimmap_ranges.rda")) / 1024) # Sizes in Kb
