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

# Plot the backbone with its node IDs, and the clades used below:  Family_B in limegreen, Subfam_B2 in gold
ape::plot.phylo(toy_backbone, cex = 0.6, no.margin = TRUE)
ape::nodelabels(cex = 0.6, frame = "rect", bg = "lightblue", col = "darkblue")
ape::nodelabels(node = 33, cex = 0.6, frame = "rect", bg = "limegreen", col = "black")
ape::nodelabels(node = 47, cex = 0.6, frame = "rect", bg = "gold", col = "black")

## ----helpers------------------------------------------------------------------
# Helper 1: build a missing_taxa_list with n taxa to graft within a target clade, defined by tips
make_taxa_list <- function(n, target_tips, prefix = "Tax_sp", start = 1)
{
  taxa <- paste0(prefix, seq(start, length.out = n))
  list(Taxon_name = taxa, Target_clade = setNames(rep(list(target_tips), n), taxa))
}

# Helper 2: build a missing_taxa_list with n taxa to graft as one missing clade,
# sister to a clade defined by tips
make_clade_list <- function(n, sister_tips, prefix = "New_sp", start = 1)
{
  taxa <- paste0(prefix, seq(start, length.out = n))
  list(Taxon_name = taxa, Sister_clade = setNames(rep(list(sister_tips), n), taxa))
}

# Helper 3: plot a grafted phylogeny, with the taxa of each grafting step in a different color
# Optionally, the branches of the missing clade (the subtree descending from the MRCA of its taxa)
# are drawn in red
plot_steps <- function(phy, steps, colors, main, clade_taxa = NULL)
{
  tip_colors <- rep("black", ape::Ntip(phy))
  for (i in seq_along(steps)) tip_colors[phy$tip.label %in% steps[[i]]] <- colors[i]
  edge_colors <- rep("black", nrow(phy$edge))
  if (!is.null(clade_taxa))
  {
    clade_node <- ape::getMRCA(phy, clade_taxa)
    clade_nodes <- c(clade_node, phytools::getDescendants(phy, clade_node))
    edge_colors[phy$edge[, 2] %in% clade_nodes] <- "red"
  }
  ape::plot.phylo(phy, tip.color = tip_colors, edge.color = edge_colors, edge.width = 1.5,
                  cex = 0.55, no.margin = FALSE, main = main)
}

## ----both_types, fig.width = 6, fig.height = 6, out.width = "100%"------------
# ------ Graft missing taxa and missing clades in one run ------ #

# Step 1: Build the lists of missing taxa and missing clades
taxa_list <- make_taxa_list(n = 5, target_tips = c("t8", "t22"))
clades_list <- make_clade_list(n = 3, sister_tips = c("t22", "t5"))
missing_taxa_list <- list(taxa = taxa_list, clades = clades_list)

# Step 2: Graft both types using the default order for grafting steps: 1/ Missing taxa; 2/ Missing clades
# The rationale for the order of the steps is discussed in the Section 2
grafted <- graft_missing_taxa_and_clades(phy = toy_backbone, missing_taxa_list = missing_taxa_list,
                                         seed = 24, verbose = FALSE)

# Step 3: Extract the grafting information: one table for each type
info <- extract_grafting_info(grafted, missing_taxa_list, phy_backbone = toy_backbone)
names(info)
info$taxa[, c("Taxon_name", "parent_node", "parent_tips", "graft_age")]
info$clades[, c("Missing_clade", "Taxa", "parent_node", "graft_age", "crown_age")]

# Step 4: Plot the grafted phylogeny, with the missing taxa in blue, and the taxa of the missing clade in red
plot_steps(grafted, steps = list(taxa_list$Taxon_name, clades_list$Taxon_name), colors = c("blue", "red"),
           main = "Missing taxa (blue) and missing clade (red)")

## ----order, fig.width = 10, fig.height = 6, out.width = "100%"----------------
# ------ Order of the steps: taxa first vs clades first ------ #

# Graft with the default order: missing taxa first, then missing clades
grafted_taxa_first <- graft_missing_taxa_and_clades(phy = toy_backbone, missing_taxa_list = missing_taxa_list,
                                                    grafting_order = c("taxa", "clades"),
                                                    seed = 24, verbose = FALSE)
# Graft with the reverse order: missing clades first, then missing taxa
grafted_clades_first <- graft_missing_taxa_and_clades(phy = toy_backbone, missing_taxa_list = missing_taxa_list,
                                                      grafting_order = c("clades", "taxa"),
                                                      seed = 24, verbose = FALSE)

# Are the taxa of the missing clade still a monophyletic group?
ape::is.monophyletic(grafted_taxa_first, clades_list$Taxon_name)
ape::is.monophyletic(grafted_clades_first, clades_list$Taxon_name)

# Plot the phylogenies: the subtree of the missing clade is in red
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot_steps(grafted_taxa_first, steps = list(taxa_list$Taxon_name, clades_list$Taxon_name),
           colors = c("blue", "red"),
           clade_taxa = clades_list$Taxon_name, main = "Taxa first, then clades")
plot_steps(grafted_clades_first, steps = list(taxa_list$Taxon_name, clades_list$Taxon_name),
           colors = c("blue", "red"),
           clade_taxa = clades_list$Taxon_name, main = "Clades first, then taxa")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----manual, fig.width = 6.5, fig.height = 7, out.width = "100%"--------------
# ------ Three successive graftings ------ #

# Clades defined by tips: Subfam_B2 and Family_B
subfam_b2 <- c("t22", "t5")
family_b <- c("t8", "t22")

# Lists of missing taxa and clades
taxa_1_list <- make_taxa_list(n = 2, target_tips = subfam_b2, prefix = "Subfam_sp")
clade_2_list <- make_clade_list(n = 3, sister_tips = subfam_b2, prefix = "Clade_sp")
taxa_3_list <- make_taxa_list(n = 5, target_tips = family_b, prefix = "Fam_sp")

# Total richness of the clades used at each step
richness_1 <- list(Target_clade = list(subfam_b2), total_richness = 12)
richness_2 <- list(Sister_clade = list(subfam_b2), total_richness = 12)
richness_3 <- list(Target_clade = list(family_b), total_richness = 40)

# Step 1: graft the taxa within Subfam_B2, on the backbone
grafted_1 <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_1_list,
                                total_richness_per_clades = richness_1, seed = 30, verbose = FALSE)
# Step 2: graft the missing clade as sister of Subfam_B2, on the output of step 1
grafted_2 <- graft_missing_clades(phy = grafted_1, missing_taxa_list = clade_2_list,
                                  total_richness_per_clades = richness_2, seed = 30, verbose = FALSE)
# Step 3: graft the taxa within Family_B, on the output of step 2
grafted_3 <- graft_missing_taxa(phy = grafted_2, missing_taxa_list = taxa_3_list,
                                total_richness_per_clades = richness_3, seed = 30, verbose = FALSE)

# Plot the final phylogeny, with the taxa of each step in a different color 
# and the subtree of the missing clade in red
plot_steps(grafted_3, steps = list(taxa_1_list$Taxon_name,
                                   clade_2_list$Taxon_name,
                                   taxa_3_list$Taxon_name),
           colors = c("blue", "red", "darkgreen"), clade_taxa = clade_2_list$Taxon_name,
           main = "Step 1 (blue), Step 2 (red), Step 3 (green)")

## ----node_ids-----------------------------------------------------------------
# Node ID of Subfam_B2 in the backbone, and after each step
c(Backbone = ape::getMRCA(toy_backbone, subfam_b2),
  Step_1 = ape::getMRCA(grafted_1, subfam_b2),
  Step_2 = ape::getMRCA(grafted_2, subfam_b2),
  Step_3 = ape::getMRCA(grafted_3, subfam_b2))

