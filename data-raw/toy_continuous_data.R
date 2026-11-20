##### Build the toy datasets of a continuous trait of phylograftR #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Run from the root of the package, after loading it (devtools::load_all() or devtools::install()).
# Requires the R packages 'contsimmap' (not on CRAN: remotes::install_github("bstaggmartin/contsimmap")) and 'deepSTRAPP' (on CRAN).
# Produces:
#  * data/toy_tip_values.rda: tip values of a continuous trait simulated on 'toy_backbone'
#  * data/toy_contMap.rda: contMap of this trait (point estimate: ML ancestral values interpolated along the branches)
#  * data/toy_multiSimmap_continuous.rda: multiSimmap of 100 histories of this trait (continuous stochastic mapping), discretized into three states
#
# The toy missing taxa list is not changed: the values of the continuous trait of the missing taxa are defined in the vignette 'continuous_data'.


##### 1/ Simulate tip values of a continuous trait on the backbone #####

data(toy_backbone)

# Seed chosen among the first 150 seeds, based on the discretization of the contMap with the terciles of the tip values (see 2/):
# the three states cover similar shares of the branches (35%, 43% and 22% for "low", "medium" and "high"),
# and, for the 7 toy missing taxa, the states allow different shares of the available branches (14 of the 21 combinations taxon x state
# have between 15% and 85% of their available branches in the state), so that the continuous constraints are informative.
set.seed(37)
toy_tip_values <- phytools::fastBM(toy_backbone,
                                   a = 0, # Ancestral value at the root
                                   sig2 = 1) # Evolutionary rate as instantaneous variance of the Brownian Motion
stopifnot(identical(names(toy_tip_values), toy_backbone$tip.label))


##### 2/ Build the contMap: point estimate #####

# ML ancestral values under BM, interpolated along the branches
toy_contMap <- phytools::contMap(tree = toy_backbone, x = toy_tip_values, plot = FALSE)

# Check the discretization of this point estimate.
# 'tip_values' provides the observed values: the contMap only approximates them (interpolation within segments),
# so the cut points (terciles) and the states of the tips are based on the observed data
toy_simmap_check <- convert_contMap_to_simmap(toy_contMap, probs = c(1/3, 2/3), states = c("low", "medium", "high"),
                                              tip_values = toy_tip_values)
print(attr(toy_simmap_check, "breaks"))
print(round(colSums(toy_simmap_check$mapped.edge) / sum(toy_simmap_check$mapped.edge), 2))


##### 3/ Build the multiSimmap: 100 continuous stochastic maps, discretized #####

# We use 'contsimmap' directly to credit the authors of the continuous stochastic mapping approach and show the steps.
# For your own data, deepSTRAPP::prepare_trait_data() is recommended: it fits the evolutionary model(s), runs the continuous stochastic mapping
# and converts the maps into contMaps in a single call, with options to fine-tune the process.
# It can also output directly both the interpolated point estimate as in phytools::contMap(),
# and the set of stochastic continuous maps simulated with contsimmap::make.contsimmap().
# See the commented alternative below.

# Estimate the evolutionary rate (sigma²) of the BM
BM_fit <- phytools::anc.ML(tree = toy_backbone, x = toy_tip_values, model = "BM")
sigma2 <- BM_fit$sig2

# Simulate 100 continuous stochastic maps (histories of the trait) conditioned on the tip data, the phylogeny and the BM model
set.seed(1234)
toy_contsimmap <- contsimmap::make.contsimmap(tree = toy_backbone,
                                              trait.data = toy_tip_values,
                                              Xsig2 = list(sigma2),
                                              nsim = 100)
dim(toy_contsimmap) # 58 edges, 1 trait, 100 simulations

# Convert into a list of contMaps
toy_contMaps_list <- deepSTRAPP::convert_contsimmap_to_contMaps(contsimmap = toy_contsimmap, verbose = TRUE)

# Alternative: all the steps above in one call with deepSTRAPP (recommended for your own data). Not run.
# toy_trait_data <- deepSTRAPP::prepare_trait_data(tip_data = toy_tip_values, trait_data_type = "continuous", phylo = toy_backbone,
#                                                  seed = 1234, evolutionary_models = "BM",
#                                                  run_stochastic_maps = TRUE, nb_simulations = 100, plot_map = FALSE)
# toy_contMaps_list <- toy_trait_data$contMaps

# Convert the contMaps into simmaps with three states, using the terciles of the observed tip values
toy_multiSimmap_continuous <- convert_contMap_to_simmap(toy_contMaps_list, probs = c(1/3, 2/3),
                                                        states = c("low", "medium", "high"),
                                                        tip_values = toy_tip_values)

# Checks
stopifnot(inherits(toy_multiSimmap_continuous, "multiSimmap"),
          length(toy_multiSimmap_continuous) == 100,
          identical(attr(toy_multiSimmap_continuous, "states"), c("low", "medium", "high")),
          isTRUE(all.equal(attr(toy_multiSimmap_continuous, "breaks"), unname(quantile(toy_tip_values, probs = c(1/3, 2/3))))))
print(attr(toy_multiSimmap_continuous, "breaks"))
print(object.size(toy_multiSimmap_continuous), units = "Mb") # Keep the data of the package small (CRAN: less than 5 Mb in total)


##### 4/ Save objects #####

dir.create("data", showWarnings = FALSE)
save(toy_tip_values, file = "data/toy_tip_values.rda", compress = "gzip", version = 2)
save(toy_contMap, file = "data/toy_contMap.rda", compress = "gzip", version = 2)
save(toy_multiSimmap_continuous, file = "data/toy_multiSimmap_continuous.rda", compress = "xz", version = 2)
print(file.size(c("data/toy_tip_values.rda", "data/toy_contMap.rda", "data/toy_multiSimmap_continuous.rda")) / 1024) # Sizes in Kb

