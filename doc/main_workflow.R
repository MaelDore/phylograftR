## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----backbone-----------------------------------------------------------------
# ------ Step 1: Load the backbone phylogeny and its taxonomy ------ #

data(toy_backbone)
data(toy_backbone_taxonomic_table)

# Backbone phylogeny
toy_backbone

# Taxonomy of the backbone tips
head(toy_backbone_taxonomic_table)

## ----tip_data-----------------------------------------------------------------
# ------ Step 2: Load the tip data ------ #

data(toy_simmap)

# Named vector of states: one value per tip of the backbone
tip_data <- phytools::getStates(toy_simmap, type = "tips")
head(tip_data)
table(tip_data)

## ----simmap-------------------------------------------------------------------
# ------ Step 3: Build 10 stochastic maps of the trait ------ #

# Sample 10 stochastic maps with an equal-rates model of evolution
set.seed(1234)
simmaps <- phytools::make.simmap(tree = toy_backbone, x = tip_data[toy_backbone$tip.label],
                                 model = "ER", nsim = 10, message = FALSE)

# A set of maps is a 'multiSimmap'
class(simmaps)
length(simmaps)

## ----missing_taxa_table-------------------------------------------------------
# ------ Step 4: Describe the missing lineages ------ #

data(toy_missing_taxa_taxonomic_table)
missing_taxa_table <- toy_missing_taxa_taxonomic_table

# The trait states in which each taxon can be grafted, separated by commas (empty: all states allowed)
missing_taxa_table$allowed_states <- c("green, red", NA, NA, "blue, green", "green",
                                       NA, "blue, red", NA, NA, NA)
missing_taxa_table[, c("Taxa", "Type", "Target_clade", "Sister_clade", "allowed_states")]

## ----missing_taxa_list--------------------------------------------------------
# ------ Step 5: Build the missing_taxa_list ------ #

missing_taxa_lists <- build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table = toy_backbone_taxonomic_table,
  missing_taxa_taxonomic_table = missing_taxa_table,
  phy = toy_backbone,
  verbose = FALSE)

# Both types were provided: one list for the missing taxa, one for the missing clades
names(missing_taxa_lists)
missing_taxa_lists$taxa$Taxon_name
missing_taxa_lists$clades$Taxon_name

## ----total_richness-----------------------------------------------------------
# ------ Step 6: Provide the total richness of some clades ------ #

richness_table <- data.frame(Clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"),
                             total_richness = c(25, NA, 14, 8))

richness_per_clades <- build_total_richness_per_clades(data = richness_table, phy = toy_backbone,
                                                       missing_taxa_list = missing_taxa_lists,
                                                       verbose = FALSE)

## ----check_plan---------------------------------------------------------------
# ------ Step 7: Check that the grafting is feasible ------ #

grafting_plan <- check_grafting_plan(phy = toy_backbone,
                                     missing_taxa_list = missing_taxa_lists,
                                     simmap = simmaps)

## ----graft--------------------------------------------------------------------
# ------ Step 8: Graft the missing taxa and the missing clades ------ #

grafted_phylos <- graft_missing_taxa_and_clades(
  phy = toy_backbone,
  missing_taxa_list = missing_taxa_lists,
  total_richness_per_clades = richness_per_clades,
  simmap = simmaps,
  nb_simulations = 10,
  seed = 1234,
  verbose = FALSE)

# One element per stochastic map, each with 10 grafted phylogenies (one per simulation)
class(grafted_phylos)
names(grafted_phylos)
class(grafted_phylos$map_1)
length(grafted_phylos$map_1)

## ----graft_phylo--------------------------------------------------------------
# Extract the first grafted phylogeny
first_phylo <- grafted_phylos$map_1[[1]]
# Explore tips of the first grafted phylogeny
ape::Ntip(first_phylo)
# Identify the grafted tips
first_phylo$tip.label[first_phylo$grafted_tips]

## ----check_constraints--------------------------------------------------------
# ------ Step 9: Check the constraints on the grafted phylogenies ------ #

constraints <- check_grafting_constraints(phy_grafted = grafted_phylos,
                                          phy_backbone = toy_backbone,
                                          missing_taxa_list = missing_taxa_lists)

## ----extract_info-------------------------------------------------------------
# ------ Step 10: Extract the grafting information ------ #

grafting_info <- extract_grafting_info(phy_grafted = grafted_phylos,
                                       missing_taxa_list = missing_taxa_lists,
                                       phy_backbone = toy_backbone)
names(grafting_info)

# Information on the missing taxa, recorded per taxon, per grafted phylogeny
# Here, only info for the first map and first simulation are printed
head(grafting_info$taxa[, c("map", "simulation", "Taxon_name", "parent_tips", "divergence_age")])

# The divergence age of the taxa varies among the maps and simulations
# You can use this output to summarize all grafting information per taxon across grafted phylgoenies
aggregate(divergence_age ~ Taxon_name, data = grafting_info$taxa,
          FUN = function(x) round(c(mean = mean(x), median = median(x), min = min(x), max = max(x)), 1))

## ----plot_grafted, fig.width = 7, fig.height = 6------------------------------
# ------ Step 11: Plot one grafted phylogeny ------ #

# Extract the first grafted phylogeny
first_phylo <- grafted_phylos$map_1[[1]]

# Plot its simmap (default, when available)
# Grafted branches are dashed lines, grafted tips in red
plot_grafted(first_phylo, fsize = 0.7)

# Plot the plain phylogeny
# Grafted branches and tips are red
plot_grafted(first_phylo, use_simmap = FALSE, cex = 0.7)

## ----flatten------------------------------------------------------------------
# ------ Step 12: Flatten the grafted phylogenies ------ #

flat_phylos <- flatten_grafted_phylogenies(grafted_phylos)
class(flat_phylos)
length(flat_phylos)
head(names(flat_phylos))

# Information on each phylogeny: map, simulation, seed
head(attr(flat_phylos, "phylogenies"))

## ----ltt, fig.width = 6, fig.height = 4.5-------------------------------------
# ------ Step 13: Downstream analysis: LTT plots of the grafted phylogenies ------ #

# LTT of the 100 grafted phylogenies, and of the backbone
ltt_grafted <- phytools::ltt(flat_phylos, plot = FALSE, gamma = FALSE)
ltt_backbone <- phytools::ltt(toy_backbone, plot = FALSE, gamma = FALSE)

# Plot the uncertainty due to the grafting and the maps (blue), against the backbone (black)
plot(ltt_grafted, col = adjustcolor("steelblue", alpha.f = 0.3),
     main = "Lineages through time", log.lineages = FALSE)
lines(ltt_backbone$times, ltt_backbone$ltt, lwd = 2, type = "s")

