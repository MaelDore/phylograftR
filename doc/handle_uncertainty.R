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
data(toy_missing_taxa_list) # Taxonomic constraints only (target clades and excluded subclades)
data(toy_missing_taxa_list_with_allowed_states) # Same, with colors allowed for some taxa
data(toy_multiSimmap) # 10 stochastic maps of the color

# Missing taxa list with the colors allowed for some taxa
taxa_colors <- toy_missing_taxa_list_with_allowed_states

# New_sp1 is only allowed in "red" in the shipped list. We also allow "green", so that a branch
# is available for it in each of the 10 maps (see the end of Section 2 for what happens otherwise)
taxa_colors$allowed_states[["New_sp1"]] <- c("red", "green")

## ----simulations--------------------------------------------------------------
# ------ N simulations ------ #

# Graft the missing taxa 100 times, with taxonomic constraints only (no stochastic map)
grafted_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                   nb_simulations = 100, verbose = FALSE, seed = 1)

# The output is a multiPhylo with one phylogeny per simulation
class(grafted_sims)
length(grafted_sims)
head(names(grafted_sims))

# The seed of each simulation is stored in the attribute "seeds"
head(attr(grafted_sims, "seeds"))

## ----simulations_info, fig.width = 8, fig.height = 4.5, out.width = "80%"-----
# ------ Variability of the grafting among simulations ------ #

info_sims <- extract_grafting_info(grafted_sims, missing_taxa_list = toy_missing_taxa_list,
                                   phy_backbone = toy_backbone)

# Summarize the number of different branches used for each taxon, and their median grafting age
data.frame(N_simulations = tapply(info_sims$simulation, info_sims$Taxon_name, length),
           N_parent_branches = tapply(info_sims$parent_node, info_sims$Taxon_name, function (x) length(unique(x))),
           Graft_age_median = round(tapply(info_sims$graft_age, info_sims$Taxon_name, median), 2))

# Distribution of the grafting ages
par(mar = c(4.5, 4.5, 1, 1))
boxplot(graft_age ~ Taxon_name, data = info_sims, ylab = "Grafting age", xlab = "",
        las = 2, col = "grey85")
par(mar = c(5.1, 4.1, 4.1, 2.1))

## ----reproduce----------------------------------------------------------------
# ------ Reproduce one simulation ------ #

# Rerun the simulation 3 alone, with its seed
seed_3 <- attr(grafted_sims, "seeds")["sim_3"]
sim_3_alone <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                  nb_simulations = 1, verbose = FALSE, seed = seed_3)

# Same grafted phylogeny
all.equal(sim_3_alone$edge, grafted_sims[["sim_3"]]$edge)
all.equal(sim_3_alone$edge.length, grafted_sims[["sim_3"]]$edge.length)

## ----maps---------------------------------------------------------------------
# ------ P maps ------ #

# Graft with each of the 10 maps (one simulation per map by default)
grafted_maps <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                   simmap = toy_multiSimmap, verbose = FALSE, seed = 1)

# The output is a list with one grafted phylogeny per map, named after the maps
class(grafted_maps)
names(grafted_maps)
class(grafted_maps$map_1)

## ----maps_info----------------------------------------------------------------
# ------ Variability of the grafting among maps ------ #

info_maps <- extract_grafting_info(grafted_maps, missing_taxa_list = taxa_colors,
                                   phy_backbone = toy_backbone)

# Number of different branches used for each taxon across the 10 maps
tapply(info_maps$parent_node, info_maps$Taxon_name, function (x) length(unique(x)))

## ----maps_dropped-------------------------------------------------------------
# ------ A map dropped at random ------ #

# Only 3 attempts to draw a usable divergence age
grafted_dropped <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                      simmap = toy_multiSimmap, max_age_redraws = 3,
                                      verbose = FALSE, seed = 1)

# An output is produced, without the map(s) that failed
names(grafted_dropped)
attr(grafted_dropped, "failed_maps")$map

## ----maps_error, error = TRUE-------------------------------------------------
try({
# ------ A map where the constraints cannot be met ------ #

# Colors allowed for New_sp1 in the shipped list
toy_missing_taxa_list_with_allowed_states$allowed_states[["New_sp1"]]

# Graft with each of the 10 maps
graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
                   simmap = toy_multiSimmap, verbose = FALSE, seed = 1)
})

## ----maps_plan----------------------------------------------------------------
# ------ Check the grafting plan across the maps ------ #

plan_maps <- check_grafting_plan(phy = toy_backbone,
                                 missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
                                 simmap = toy_multiSimmap, verbose = FALSE)

# Number of maps where each taxon cannot be grafted
plan_maps[, c("Taxon_name", "allowed_states", "n_maps_infeasible")]

# Which maps?
per_map <- attr(plan_maps, "per_map")
per_map[!per_map$feasible, c("map", "Taxon_name", "cause")]

## ----crossed------------------------------------------------------------------
# ------ P x N ------ #

# Graft with each of the 10 maps, 5 simulations per map: 50 grafted phylogenies
grafted_crossed <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                      simmap = toy_multiSimmap, nb_simulations = 5,
                                      verbose = FALSE, seed = 1)

# One element per map, each a multiPhylo with the 5 simulations of this map
length(grafted_crossed)
sapply(grafted_crossed, length)

# Total number of grafted phylogenies
sum(sapply(grafted_crossed, length))

# Names of the phylogenies in the first map: "<map>_sim_<k>"
names(grafted_crossed$map_1)

## ----crossed_reproduce--------------------------------------------------------
# ------ Reproduce one phylogeny of the crossed output ------ #

# Seed of the simulation 2 of the map 3
seed_map_3_sim_2 <- attr(grafted_crossed$map_3, "seeds")["map_3_sim_2"]
seed_map_3_sim_2

# Graft again with the map 3 alone (a single simmap) and this seed
map_3_sim_2_alone <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                        simmap = toy_multiSimmap[[3]], nb_simulations = 1,
                                        verbose = FALSE, seed = seed_map_3_sim_2)

# Same grafted phylogeny
all.equal(map_3_sim_2_alone$edge, grafted_crossed$map_3[["map_3_sim_2"]]$edge)
all.equal(map_3_sim_2_alone$edge.length, grafted_crossed$map_3[["map_3_sim_2"]]$edge.length)

## ----crossed_info-------------------------------------------------------------
# ------ Variation among maps and within maps ------ #

info_crossed <- extract_grafting_info(grafted_crossed, missing_taxa_list = taxa_colors,
                                      phy_backbone = toy_backbone)

# For each taxon: standard deviation of the grafting age among simulations of the same map (within),
# and standard deviation of the mean grafting age of the maps (between)
sd_within <- tapply(info_crossed$graft_age, info_crossed$Taxon_name, function (x) NA)
for (taxon in names(sd_within))
{
  ages <- info_crossed[info_crossed$Taxon_name == taxon, ]
  sd_within[taxon] <- sqrt(mean(tapply(ages$graft_age, ages$map, var)))
}
sd_between <- tapply(info_crossed$graft_age, list(info_crossed$Taxon_name, info_crossed$map), mean)
data.frame(SD_within_maps = round(sd_within, 2), SD_between_maps = round(apply(sd_between, 1, sd), 2))

## ----inflation, eval = FALSE--------------------------------------------------
# # ------ Not run: more than 1000 grafted phylogenies ------ #
# 
# # 10 maps x 101 simulations = 1010 phylogenies
# graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
#                    simmap = toy_multiSimmap, nb_simulations = 101, seed = 1)
# #> Warning: 10 maps x 101 simulations = 1010 grafted phylogenies will be produced (more than 1000):
# #> this can take a long time and use a lot of memory.

## ----flatten------------------------------------------------------------------
# ------ Flatten the grafted phylogenies ------ #

# We use the output of Section 3, where 5 simulations per map where produced
# stored as a list of 10 multiPhylo, with 5 simulations each
length(grafted_crossed)
sapply(grafted_crossed, length)

# We flatten this output into a single multiPhylo with 50 grafted phylogenies
flat_phylos <- flatten_grafted_phylogenies(grafted_crossed)
class(flat_phylos)
length(flat_phylos)
head(names(flat_phylos))

# The map, the simulation and the seed of each phylogeny are kept in the attribute "phylogenies"
head(attr(flat_phylos, "phylogenies"))

## ----flatten_ltt, fig.width = 8, fig.height = 5, out.width = "80%"------------
# ------ Use the multiPhylo: lineages through time ------ #

# phytools::ltt() takes the flattened multiPhylo as such
ltt_grafted <- phytools::ltt(flat_phylos, plot = FALSE, gamma = FALSE)
ltt_backbone <- phytools::ltt(toy_backbone, plot = FALSE, gamma = FALSE)

# Plot the 50 grafted phylogenies, then add the backbone
plot(ltt_grafted, col = "grey70", log.lineages = FALSE,
     xlab = "Time since the root", ylab = "Number of lineages")
plot(ltt_backbone, add = TRUE, col = "firebrick", lwd = 2, log.lineages = FALSE)
legend("topleft", legend = c("Grafted phylogenies", "Backbone"),
       col = c("grey70", "firebrick"), lwd = 2, bty = "n")

## ----flatten_subset-----------------------------------------------------------
# ------ Subset and remove the simmaps ------ #

# Keep the first two simulations of two maps only, without the simmaps
subset_phylos <- flatten_grafted_phylogenies(grafted_crossed, maps = c("map_1", "map_3"),
                                             simulations = 1:2, keep_simmap = FALSE)
attr(subset_phylos, "phylogenies")

# Memory used with and without the simmaps (in Mb)
round(c(with_simmaps = as.numeric(object.size(flat_phylos)),
        without_simmaps = as.numeric(object.size(flatten_grafted_phylogenies(grafted_crossed,
                                                                             keep_simmap = FALSE)))) / 1e6, 2)

## ----backbones----------------------------------------------------------------
# ------ Alternative backbones ------ #

data(toy_backbone_taxonomic_table) # Taxonomy of the tips (the same for all the backbones)
data(toy_missing_taxa_taxonomic_table) # Taxonomic groups of the missing taxa

# For illustration, we build two alternative backbones from the toy backbone
# Backbone 2: the tips t16 and t3 are swapped, which changes the topology
backbone_2 <- toy_backbone
swapped <- match(c("t16", "t3"), backbone_2$tip.label)
backbone_2$tip.label[swapped] <- backbone_2$tip.label[rev(swapped)]

# Backbone 3: same topology, older ages (all the branches are 30% longer)
backbone_3 <- toy_backbone
backbone_3$edge.length <- backbone_3$edge.length * 1.30

backbones <- list(backbone_1 = toy_backbone, backbone_2 = backbone_2, backbone_3 = backbone_3)

# One color per backbone, used in the figures below
backbone_colors <- c(backbone_1 = "black", backbone_2 = "firebrick", backbone_3 = "dodgerblue")

## ----backbones_clade----------------------------------------------------------
# ------ The same list in two topologies ------ #

# Identify the tips of Tribe_A2a across the different backbone phylogenies
excluded_tips <- toy_missing_taxa_list$Excluded_subclades$New_sp1$Tribe_A2a
clade_tips <- lapply(backbones, FUN = function (phy) 
   { sort(ape::extract.clade(phy, phytools::findMRCA(phy, excluded_tips))$tip.label) }) 
clade_tips

# Compare to the tips of Tribe_A2a as identified in the taxonomic table used
# to describe the initial backbone (i.e., Backbone 1)
tribe_A2a <- toy_backbone_taxonomic_table$Taxa[which(toy_backbone_taxonomic_table$Tribe == "Tribe_A2a")]
sapply(clade_tips, function (x) all(x %in% tribe_A2a))

## ----backbones_lists----------------------------------------------------------
# ------ One missing_taxa_list per backbone ------ #

# Rebuild the missing_taxa_list based on the taxonomic table, across all backbone phylogenies
taxa_lists <- lapply(backbones, function (phy)
{
  build_missing_taxa_list_from_taxonomy(backbone_taxonomic_table = toy_backbone_taxonomic_table,
                                        missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
                                        phy = phy, verbose = FALSE)$taxa
})

# The excluded clade corresponding to the smallest clade that contains all tips labeled as
# belonging to Tribe_A2a has different tips in the backbone 2 (since we swapped tips on purpose)
sapply(taxa_lists, function (x) paste(x$Excluded_subclades$New_sp1$Tribe_A2a, collapse = ", "))

## ----backbones_graft----------------------------------------------------------
# ------ Graft on each backbone ------ #

# Loop across all backbone, using each dedicated missing_taxa_list
grafted_backbones <- lapply(seq_along(backbones), function (i)
{
  graft_missing_taxa(phy = backbones[[i]], missing_taxa_list = taxa_lists[[i]],
                     nb_simulations = 10, verbose = FALSE, seed = 1)
})
names(grafted_backbones) <- names(backbones)

# Since we provided no stochastic maps for constraints,
# we only produce 10 simulations per backbone.
# The output is therefore a list with one multiPhylo of 10 phylogenies per backbone
sapply(grafted_backbones, length)

## ----backbones_maps, eval = FALSE---------------------------------------------
# # ------ Not run: loop over backbones, with one multiSimmap per backbone ------ #
# 
# grafted_backbones_maps <- lapply(seq_along(backbones), function (i)
# {
#   graft_missing_taxa(phy = backbones[[i]], missing_taxa_list = taxa_lists[[i]],
#                      simmap = multiSimmaps_per_backbone[[i]],
#                      nb_simulations = 5, verbose = FALSE, seed = 1)
# })
# names(grafted_backbones_maps) <- names(backbones)

## ----backbones_info, fig.width = 9, fig.height = 4.5, out.width = "90%"-------
# ------ Compare the grafting among backbones ------ #

info_backbones <- do.call(rbind, lapply(names(backbones), function (name)
{
  info <- extract_grafting_info(grafted_backbones[[name]], 
                                missing_taxa_list = taxa_lists[[name]],
                                phy_backbone = backbones[[name]])
  cbind(backbone = name, info)
}))

# Grafting ages of the two taxa of the family A ("New_sp1", "New_sp2"), per backbone
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))
for (taxon in c("New_sp1", "New_sp2"))
{
  boxplot(graft_age ~ backbone, data = info_backbones[info_backbones$Taxon_name == taxon, ],
          main = taxon, ylab = "Grafting age", xlab = "Backbone",
          names = 1:3, col = adjustcolor(backbone_colors, alpha.f = 0.35), border = backbone_colors)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----backbones_flatten--------------------------------------------------------
# ------ Pool the phylogenies of all the backbones ------ #

# Flatten results of each backbone and name the resulting grafted phylogenies accordingly
flat_backbones <- lapply(names(grafted_backbones), function (name)
{
  flat <- flatten_grafted_phylogenies(grafted_backbones[[name]])
  names(flat) <- paste0(name, "_", names(flat))
  flat
})
# Concatenate in a shared list all flatten outputs from each backbone
all_phylos <- do.call(c, flat_backbones)

class(all_phylos)
length(all_phylos)
head(names(all_phylos))
table(sub("_sim_[0-9]+$", "", names(all_phylos)))

## ----backbones_ltt, fig.width = 8, fig.height = 5, out.width = "80%"----------
# ------ Lineages through time of the backbones and of the grafted phylogenies ------ #

backbone_of_phylo <- sub("_sim_[0-9]+$", "", names(all_phylos))

# phytools::ltt() takes the pooled multiPhylo: we split it by backbone to color the curves
ltt_pooled <- lapply(names(backbones), function (name)
{
  phytools::ltt(all_phylos[backbone_of_phylo == name], plot = FALSE, gamma = FALSE)
})
names(ltt_pooled) <- names(backbones)
ltt_backbones <- lapply(backbones, phytools::ltt, plot = FALSE, gamma = FALSE)

# The backbone 3 is the oldest: start the plot with it, and zoom on the recent part of the time axis
plot(ltt_pooled$backbone_3, col = adjustcolor(backbone_colors["backbone_3"], alpha.f = 0.25),
     log.lineages = FALSE, xlim = c(30, 65),
     xlab = "Time since the root", ylab = "Number of lineages")
for (name in c("backbone_2", "backbone_1"))
{
  plot(ltt_pooled[[name]], add = TRUE, col = adjustcolor(backbone_colors[name], alpha.f = 0.25), log.lineages = FALSE)
}

# Backbones: thick lines
for (name in names(backbones))
{
  plot(ltt_backbones[[name]], add = TRUE, col = backbone_colors[name], lwd = 2.5, log.lineages = FALSE)
}
legend("topleft", legend = names(backbones), col = backbone_colors, lwd = 2.5, bty = "n")

