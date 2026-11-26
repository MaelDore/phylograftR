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
# Helper 1: build a missing_taxa_list with n taxa to graft as one missing clade, sister to a clade defined by tips
make_clade_list <- function(n, sister_tips, prefix = "New_sp", start = 1, sister_name = NULL)
{
  taxa <- paste0(prefix, seq(start, length.out = n))
  out <- list(Taxon_name = taxa,
              Sister_clade = setNames(rep(list(sister_tips), n), taxa))
  if (!is.null(sister_name)) out$Sister_clade_name <- setNames(rep(sister_name, n), taxa)
  return(out)
}

# Helper 2: combine the lists of several missing clades
combine_lists <- function(...)
{
  lists <- list(...)
  list(Taxon_name = unlist(lapply(lists, `[[`, "Taxon_name")),
       Sister_clade = do.call(c, lapply(lists, `[[`, "Sister_clade")))
}

# Helper 3: plot two grafted phylogenies side by side
plot_two <- function(phy_1, title_1, phy_2, title_2)
{
  par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
  plot_grafted(phy_1, use_simmap = FALSE, cex = 0.55, main = title_1)
  plot_grafted(phy_2, use_simmap = FALSE, cex = 0.55, main = title_2)
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
}

# Helper 4: summarize the grafting of the missing clades in a table
summarize_clades <- function(phy_grafted, missing_taxa_list)
{
  info <- extract_grafting_info(phy_grafted, missing_taxa_list, phy_backbone = toy_backbone)
  data.frame(Missing_clade = info$Missing_clade, Taxa = info$Taxa, n_taxa = info$n_taxa,
             sister_tips = info$parent_tips,
             graft_age = round(info$graft_age, 2), crown_age = round(info$crown_age, 2))
}

## ----two_clades, fig.width = 10, fig.height = 5.5, out.width = "100%"---------
# ------ Missing clades: taxa grouped by their sister clade ------ #

# Missing clade 1: three taxa, sister to Subfam_B1 (node 34)
clade_1 <- make_clade_list(n = 3, sister_tips = c("t8", "t7"))
# Missing clade 2: two taxa, sister to Subfam_A2 (node 38)
clade_2 <- make_clade_list(n = 2, sister_tips = c("t28", "t26"), start = 4)

# Missing taxa list with both missing clades
taxa_two_clades <- combine_lists(clade_1, clade_2)

# Graft the missing clades
grafted_two_clades <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_two_clades,
                                           seed = 2, verbose = FALSE)

# Summarize grafting info on the missing clades
summarize_clades(grafted_two_clades, taxa_two_clades)

# Plot the backbone vs the grafted phylogeny
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
# Left: the backbone with its node IDs, and the sister clades: node 34 in limegreen, node 38 in gold
ape::plot.phylo(toy_backbone, cex = 0.55, main = "Backbone: sister clades (green, gold)")
ape::nodelabels(cex = 0.55, frame = "rect", bg = "lightblue", col = "darkblue")
ape::nodelabels(node = 34, cex = 0.55, frame = "rect", bg = "limegreen", col = "black")
ape::nodelabels(node = 38, cex = 0.55, frame = "rect", bg = "gold", col = "black")
# Right: the grafted phylogeny
plot_grafted(grafted_two_clades, use_simmap = FALSE, cex = 0.55, main = "Grafted missing clades")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----min_tips, fig.width = 10, fig.height = 5.5, out.width = "100%"-----------
# ------ min_tips: clade used to fit the birth-death model ------ #

taxa_gen1 <- make_clade_list(n = 4, sister_tips = c("t13", "t27"))

grafted_min3 <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_gen1,
                                     min_tips = 3, return_BD_model_summary = TRUE,
                                     seed = 1234, verbose = FALSE)
grafted_min5 <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_gen1,
                                     min_tips = 5, return_BD_model_summary = TRUE,
                                     seed = 1234, verbose = FALSE)

# Summary of the BD models: sister clade vs enclosing clade used to fit the diversification model
cols <- c("n_taxa", "initial_clade_n_tips", "enclosing", "fitted_clade_n_tips", "min_tips",
          "b", "d", "net_diversification", "relative_extinction", "rho")
rbind(extract_BD_model_summary(grafted_min3)[, cols], extract_BD_model_summary(grafted_min5)[, cols])


## ----richness, fig.width = 10, fig.height = 5.5, out.width = "100%", warning = FALSE----
# ------ Total richness: no unsampled species vs many unsampled species ------ #

# Create list of taxa to graft as a missing clade sister to Subfam_A1
taxa_subfam_a1 <- make_clade_list(n = 4, sister_tips = c("t19", "t10"), sister_name = "Subfam_A1")


# Create richness table to account for unsampled species
# Total richness of Subfam_A1: 60 species, so 60 - 8 = 52 unsampled species
richness_subfam_a1 <- list(Sister_clade = list(c("t19", "t10")), Sister_clade_name = "Subfam_A1",
                       total_richness = 60)

# Graft with and without accounting for unsampled species
grafted_no_richness <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_subfam_a1,
                                            return_BD_model_summary = TRUE, seed = 1234, verbose = FALSE)
grafted_richness <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_subfam_a1,
                                         total_richness_per_clades = richness_subfam_a1,
                                         return_BD_model_summary = TRUE, seed = 1234, verbose = FALSE)

# Inspect the effective sampling fraction and the rates of the BD models
cols <- c("n_taxa", "n_unsampled", "rho", "b", "d", "net_diversification")
rbind(extract_BD_model_summary(grafted_no_richness)[, cols],
      extract_BD_model_summary(grafted_richness)[, cols])

##  The sampling fraction (rho) is lower with many unsampled species
##  As a result, the speciation (b) and extinction (d) rates are higher


## ----richness_ages, warning = FALSE-------------------------------------------

# Graft 100 times with and without accounting for unsampled species
sims_no_richness <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_subfam_a1,
                                         nb_simulations = 100, seed = 1234, verbose = FALSE)
sims_richness <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_subfam_a1,
                                      total_richness_per_clades = richness_subfam_a1,
                                      nb_simulations = 100, seed = 1234, verbose = FALSE)

# Extract the graft ages and crown ages of the missing clade from the simulations
info_no_richness <- extract_grafting_info(sims_no_richness, taxa_subfam_a1, phy_backbone = toy_backbone)
info_richness <- extract_grafting_info(sims_richness, taxa_subfam_a1, phy_backbone = toy_backbone)

# Summarize the ages across the simulations
ages_summary <- function(info) c(mean_graft_age = mean(info$graft_age),
                                 median_graft_age = median(info$graft_age),
                                 mean_crown_age = mean(info$crown_age),
                                 median_crown_age = median(info$crown_age))
round(rbind("No total richness" = ages_summary(info_no_richness),
            "Total richness = 60" = ages_summary(info_richness)), 2)

## As expected, the crown age of the missing clade is much younger with many unsampled species,
## while the graft age is only slightly younger.


## ----simmap_plot, fig.width = 7, fig.height = 6.5, out.width = "100%"---------
# ------ The backbone simmap ------ #

data(toy_simmap)
state_colors <- c(blue = "blue", green = "green3", red = "red")

# Plot the simmap, with the sister clade (node 33) in limegreen
par(mar = c(0.5, 0.5, 2.5, 0.5))
phytools::plotSimmap(toy_simmap, colors = state_colors, fsize = 0.6, ftype = "i", lwd = 2)
title("Backbone simmap (sister clade: node 33)", line = -1, cex.main = 1)
ape::nodelabels(node = 33, cex = 0.6, frame = "rect", bg = "limegreen", col = "black")
par(mar = c(5.1, 4.1, 4.1, 2.1))

## ----states, fig.width = 10, fig.height = 6, out.width = "100%"---------------
# ------ Allowed states: blue only vs green only ------ #

# Missing clade sister to Family_B
taxa_fam_b <- make_clade_list(n = 3, sister_tips = c("t8", "t22"))

# Graft the missing clade with only blue as 'allowed_states'
grafted_blue <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_fam_b,
                                     simmap = toy_simmap, allowed_states = "blue",
                                     seed = 48, verbose = FALSE)
# Graft the same missing clade, with only green as 'allowed_states'
grafted_green <- graft_missing_clades(phy = toy_backbone, missing_taxa_list = taxa_fam_b,
                                      simmap = toy_simmap, allowed_states = "green",
                                      seed = 48, verbose = FALSE)

# The simmap is updated: the whole missing clade has the state of the point where it was grafted
sapply(list("Blue only" = grafted_blue, "Green only" = grafted_green),
       function(phy_grafted) unique(phytools::getStates(phy_grafted$simmap, type = "tips")[taxa_fam_b$Taxon_name]))

# Helper: plot the updated simmap of a grafted phylogeny
# Grafted branches are dashed, and grafted tips are marked with a point in the color of their final state
plot_states <- function(phy_grafted, main)
{
  plot_grafted(phy_grafted, which = "edges", simmap_colors = state_colors, fsize = 0.55, main = main)
  grafted_tips <- which(phy_grafted$grafted_tips)
  tip_states <- phytools::getStates(phy_grafted$simmap, type = "tips")[phy_grafted$tip.label[grafted_tips]]
  ape::tiplabels(tip = grafted_tips, pch = 19, cex = 1.5, col = state_colors[tip_states])
}

# Plot both updated simmaps
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot_states(grafted_blue, main = "A: Blue only")
plot_states(grafted_green, main = "B: Green only")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

