## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----data, fig.width = 6, fig.height = 5, out.width = "100%"------------------
# ------ Toy data ------ #

data(toy_backbone)
data(toy_simmap_ranges)
data(toy_missing_taxa_list_with_allowed_ranges)

# Plot the stochastic map of ranges
range_colors <- c(A = "gold", B = "dodgerblue", AB = "limegreen", BC = "firebrick2")
phytools::plotSimmap(toy_simmap_ranges, colors = range_colors, fsize = 0.6, ftype = "i", lwd = 2)
legend("bottomleft", legend = names(range_colors), fill = range_colors, title = "Range", bty = "n")

# Time spent in each range, summed on all the branches
round(colSums(toy_simmap_ranges$mapped.edge))

## ----list---------------------------------------------------------------------
# Load missing taxa list with the allowed ranges of the 7 missing taxa
taxa_ranges <- toy_missing_taxa_list_with_allowed_ranges

# Convert to a data.frame for printing neatly the information
data.frame(Taxon_name = taxa_ranges$Taxon_name,
           Target_clade_name = taxa_ranges$Target_clade_name,
           Allowed_ranges = sapply(taxa_ranges$allowed_states,
                                   function(x) if (is.null(x)) "NULL"
                                               else paste(x, collapse = ", ")),
           row.names = NULL)

## ----matching-----------------------------------------------------------------
# ------ Which allowed ranges match each range of the map? ------ #

map_ranges <- c("A", "B", "AB", "BC")
allowed <- c("A", "B", "AB", "BC")

# A range of the map matches an allowed range if it contains all of its areas
contains <- function(map_range, allowed_range)
  all(strsplit(allowed_range, "")[[1]] %in% strsplit(map_range, "")[[1]])

data.frame(Map_range = map_ranges,
           Exact_matching = map_ranges,
           Range_matching = sapply(map_ranges, function(r) 
              paste(allowed[sapply(allowed, contains, map_range = r)], collapse = ", ")),
           row.names = NULL)

## ----plan---------------------------------------------------------------------
# ------ Check the plan with the two options ------ #

# Exact matching, with the global allowed ranges "A" and "B"
plan_exact <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_ranges,
                                  simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                                  treat_states_as_ranges = FALSE, verbose = FALSE)
plan_exact[, c("Taxon_name", "allowed_states", "states_available", "valid_fraction", "feasible", "cause")]

# Range matching, with the same global allowed ranges
plan_ranges <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_ranges,
                                   simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                                   treat_states_as_ranges = TRUE, verbose = FALSE)
plan_ranges[, c("Taxon_name", "allowed_states", "states_available", "valid_fraction", "feasible", "cause")]

## ----graft_exact_error--------------------------------------------------------
# ------ Exact matching: infeasible plan ------ #

grafted_error <- tryCatch(
  graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_ranges,
                     simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                     treat_states_as_ranges = FALSE, verbose = FALSE, seed = 1),
  error = function(e) cat(conditionMessage(e)))

## ----graft_exact, warning = FALSE---------------------------------------------
# ------ Exact matching: list all the ranges ------ #

grafted_exact <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_ranges,
                                    simmap = toy_simmap_ranges, allowed_states = c("A", "B", "BC"),
                                    treat_states_as_ranges = FALSE, verbose = FALSE, seed = 1)

## ----graft_ranges-------------------------------------------------------------
# ------ Range matching ------ #

grafted_ranges <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_ranges,
                                     simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                                     treat_states_as_ranges = TRUE, verbose = FALSE, seed = 1)

## ----plot_ranges, fig.width = 10, fig.height = 5.5, out.width = "100%"--------
# ------ Compare the grafted phylogenies ------ #

# Helper: plot the updated simmap of a grafted phylogeny
plot_ranges <- function(phy_grafted, main)
{
  plot_grafted(phy_grafted, which = "edges", simmap_colors = range_colors, fsize = 0.55, main = main)
  grafted_tips <- which(phy_grafted$grafted_tips)
  tip_ranges <- phytools::getStates(phy_grafted$simmap, type = "tips")[phy_grafted$tip.label[grafted_tips]]
  ape::tiplabels(tip = grafted_tips, pch = 19, cex = 1.5, col = range_colors[tip_ranges])
}

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot_ranges(grafted_exact, main = "A: Exact matching")
plot_ranges(grafted_ranges, main = "B: Range matching")
legend("bottomleft", legend = names(range_colors), fill = range_colors, title = "Range", bty = "n", cex = 0.8)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## ----info---------------------------------------------------------------------
# ------ Ranges at grafting and inherited ranges ------ #

# Extract grafting information for the two matching types
info_exact <- extract_grafting_info(grafted_exact, missing_taxa_list = taxa_ranges)
info_ranges <- extract_grafting_info(grafted_ranges, missing_taxa_list = taxa_ranges)

# Summarize information in a data.frame
compare_ranges <- data.frame(Taxon_name = info_exact$Taxon_name,
                             Exact_parent = info_exact$parent_state, Exact_inherited = info_exact$tip_state,
                             Ranges_parent = info_ranges$parent_state, Ranges_inherited = info_ranges$tip_state)
# Identify subset sympatric speciation event: the taxon inherits a range smaller than the range of the branch
compare_ranges$Ranges_subset_sympatry <- compare_ranges$Ranges_parent != compare_ranges$Ranges_inherited
compare_ranges

## ----check, warning = FALSE---------------------------------------------------
# ------ A posteriori check ------ #

constraints <- check_grafting_constraints(phy_grafted = grafted_ranges, phy_backbone = toy_backbone,
                                          missing_taxa_list = taxa_ranges, allowed_states = c("A", "B"),
                                          treat_states_as_ranges = TRUE, use_excluded_subclades = FALSE,
                                          verbose = FALSE)
constraints[, c("Taxon_name", "in_tree", "target_ok", "states_ok", "ranges_ok")]

