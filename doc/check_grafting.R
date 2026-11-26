## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----setup_data---------------------------------------------------------------
data(toy_backbone)
data(toy_simmap)

# One map: toy_simmap (states: blue, green, red)

# Many maps: sample 5 stochastic maps from the tip data, with an equal-rates model of evolution
tip_data <- phytools::getStates(toy_simmap, type = "tips")
set.seed(1234)
simmaps <- phytools::make.simmap(tree = toy_backbone, x = tip_data[toy_backbone$tip.label],
                                 model = "ER", nsim = 5, message = FALSE)
class(simmaps)

## ----helpers------------------------------------------------------------------
# Helper 1: build a missing_taxa_list with n taxa to graft within a target clade defined by tips
make_taxa_list <- function(n, target_tips, prefix, allowed_states = NULL, start = 1)
{
  taxa <- paste0(prefix, seq(start, length.out = n))
  out <- list(Taxon_name = taxa, Target_clade = setNames(rep(list(target_tips), n), taxa))
  if (!is.null(allowed_states)) out$allowed_states <- setNames(rep(list(allowed_states), n), taxa)
  return(out)
}

# Helper 2: combine several lists of missing taxa
combine_lists <- function(...)
{
  lists <- list(...)
  out <- list(Taxon_name = unlist(lapply(lists, `[[`, "Taxon_name")),
              Target_clade = do.call(c, lapply(lists, `[[`, "Target_clade")))
  out$allowed_states <- do.call(c, lapply(lists, function(x)
  {
    if (is.null(x$allowed_states)) setNames(rep(list(NULL), length(x$Taxon_name)), x$Taxon_name) else x$allowed_states
  }))
  return(out)
}

## ----plan_one-----------------------------------------------------------------
# ------ Check the plan with one stochastic map ------ #

# Build the list of missing taxa
taxa_blue <- make_taxa_list(n = 2, target_tips = c("t8", "t22"), prefix = "Fam_blue",
                            allowed_states = "blue")
taxa_small <- make_taxa_list(n = 1, target_tips = c("t13", "t27"), prefix = "Small_sp")
taxa_narrow <- make_taxa_list(n = 1, target_tips = c("t5", "t22"), prefix = "Narrow_sp",
                              allowed_states = "red")
taxa_mistake <- make_taxa_list(n = 1, target_tips = c("t13", "t27"), prefix = "Mistake_sp",
                               allowed_states = "red")
taxa_list <- combine_lists(taxa_blue, taxa_small, taxa_narrow, taxa_mistake)

# Check the plan with the toy stochastic map
plan_one <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_list, simmap = toy_simmap)

# The table has one row per taxon: the main columns are
plan_one[, c("Taxon_name", "target_n_tips", "allowed_states",
             "valid_fraction", "time_coverage", "feasible", "cause")]


## ----graft_error--------------------------------------------------------------
# If you try to graft anyway, the grafting will stops with an informative error
# highlighting the taxon that could not be grafted
tryCatch(graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_list, simmap = toy_simmap,
                            seed = 1, verbose = FALSE),
         error = function(e) cat(conditionMessage(e)))

## ----plan_one_fixed-----------------------------------------------------------
# Fix the allowed states
taxa_narrow_fixed <- make_taxa_list(n = 1, target_tips = c("t5", "t22"), prefix = "Narrow_sp",
                                    allowed_states = c("red", "green")) # Use red AND green as allowed states
taxa_mistake_fixed <- make_taxa_list(n = 1, target_tips = c("t13", "t27"), prefix = "Mistake_sp",
                                     allowed_states = "green") # use green instead of red
taxa_list_fixed <- combine_lists(taxa_blue, taxa_small, taxa_narrow_fixed, taxa_mistake_fixed)

# Check the fixed plan, with min_tips = 3
# By using a lower min_tips value, we ensure that the BD models are fitted strictly on our target clades
# However, as a trade-off, we may have lower accuracy in our estimates of diversification rates as we used smaller clades for fitting models
plan_fixed <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_list_fixed,
                                  simmap = toy_simmap, min_tips = 3)

## ----plan_many----------------------------------------------------------------
# ------ Check the plan with many stochastic maps ------ #

# Check grafting plan for all simmaps
plan_many <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_list, simmap = simmaps)

# One row per taxon, summarizing the maps
plan_many[, c("Taxon_name", "valid_fraction_min", "valid_fraction_mean", "valid_fraction_max", "n_maps_infeasible", "feasible", "cause")]

## ----plan_many_per_map--------------------------------------------------------

# Extract checks per map
per_map <- attr(plan_many, "per_map")

# Inspect issues with "Narrow_sp1" across all maps
per_map[per_map$Taxon_name == "Narrow_sp1", c("map", "valid_fraction", "time_coverage", "feasible", "cause")]

## ----plan_many_fixed----------------------------------------------------------
# Check grafting plan for all simmaps with the fixed list
plan_many_fixed <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_list_fixed,
                                       simmap = simmaps, min_tips = 3)

## ----check_one----------------------------------------------------------------
# ------ Check the constraints with one stochastic map ------ #

# Graft the missing taxa with the fixed list
grafted_one <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_list_fixed,
                                  simmap = toy_simmap, min_tips = 3,
                                  nb_simulations = 10, seed = 1, verbose = FALSE)

# Check the constraints
constraints_one <- check_grafting_constraints(phy_grafted = grafted_one, phy_backbone = toy_backbone,
                                              missing_taxa_list = taxa_list_fixed,
                                              use_excluded_subclades = FALSE)

## ----check_one_table----------------------------------------------------------
# Inspect table with checks per taxon per grafted phylogeny
head(constraints_one[, c("map", "simulation", "Taxon_name", "in_tree",
                         "target_ok", "excluded_ok", "states_ok", "ok")])

# Inspect overall result
attr(constraints_one, "all_ok")

# Inspect the summary table of the checks per phylogeny
head(attr(constraints_one, "phylogenies"))

## ----check_violation----------------------------------------------------------
# A stricter list: the constraints on two taxa are changed after grafting
taxa_list_strict <- taxa_list_fixed
taxa_list_strict$Target_clade[["Fam_blue1"]] <- c("t5", "t22")
taxa_list_strict$allowed_states[["Fam_blue2"]] <- "green"

# Check the constraints
constraints_strict <- check_grafting_constraints(phy_grafted = grafted_one, phy_backbone = toy_backbone,
                                                 missing_taxa_list = taxa_list_strict,
                                                 use_excluded_subclades = FALSE)

# We can use the output table to summarize the proportion of the grafted phylogenies
# in which each taxon respects each constraint
# Here, we focus on the target clade and allowed states constraints
aggregate(cbind(target_ok, states_ok, ok) ~ Taxon_name, data = constraints_strict, FUN = mean)

## ----check_many---------------------------------------------------------------
# ------ Check the constraints with many stochastic maps ------ #

# Graft with the 5 maps, and 4 simulations per map
grafted_many <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_list_fixed,
                                   simmap = simmaps, min_tips = 3,
                                   nb_simulations = 4, seed = 1, verbose = FALSE)

# Check the constraints
constraints_many <- check_grafting_constraints(phy_grafted = grafted_many, phy_backbone = toy_backbone,
                                               missing_taxa_list = taxa_list_fixed,
                                               use_excluded_subclades = FALSE)

# One table per map
names(constraints_many)
# Are all constraints respected in all the maps?
sapply(constraints_many, attr, "all_ok")

## ----info_one-----------------------------------------------------------------
# ------ Grafting information with one stochastic map ------ #

# Extract grafting information
info_one <- extract_grafting_info(phy_grafted = grafted_one, missing_taxa_list = taxa_list_fixed,
                                  phy_backbone = toy_backbone)
# Inspect the first rows describing the results for the first simulation (among 10)
head(info_one[, c("simulation", "Taxon_name", "parent_tips", "graft_age", "divergence_age", "tip_state")], 4)

# Summarize the ages of each taxon across the 10 simulations
aggregate(divergence_age ~ Taxon_name, data = info_one,
          FUN = function(x) round(c(mean = mean(x), min = min(x), max = max(x)), 2))

## ----info_states--------------------------------------------------------------
# Summarize the number of simulations with each state, for each taxon
table(info_one$Taxon_name, info_one$tip_state)

## ----info_many----------------------------------------------------------------
# ------ Grafting information with many stochastic maps ------ #

info_many <- extract_grafting_info(phy_grafted = grafted_many, missing_taxa_list = taxa_list_fixed,
                                   phy_backbone = toy_backbone)

# We grafted 5 taxa over 5 stochastic maps, and 4 simulations per map
# Therefore, the table report 5 × 5 × 4 = 100 entries
dim(info_many)

# We can then aggregate information per taxon, per map
# to inspect the effect of each map on the graftong

# Ex: Mean divergence age of each taxon, in each map (4 simulations per map)
round(tapply(info_many$divergence_age, list(info_many$Taxon_name, info_many$map), mean), 2)


## ----sampling_one-------------------------------------------------------------
# ------ Sampling fractions with one stochastic map ------ #

# # Build the list of missing taxa
taxa_sampling <- combine_lists(make_taxa_list(n = 6, target_tips = c("t8", "t22"),
                                              prefix = "Fam_sp", allowed_states = "blue"),
                               make_taxa_list(n = 2, target_tips = c("t13", "t27"), prefix = "Gen_sp"))
taxa_sampling$Target_clade_name <- c(rep("Family_B", 6), rep("Gen_1", 2))

# Build the list informing on Total richness for four clades
richness_list <- list(Target_clade = list(c("t8", "t22"), c("t8", "t7"), c("t5", "t22"), c("t13", "t27")),
                 Target_clade_name = c("Family_B", "Subfam_B1", "Subfam_B2", "Gen_1"),
                 total_richness = c(40, 20, 12, 12))

# Graft 10 simulations, and return the birth-death model summaries for the next section
grafted_sampling <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_sampling,
                                       total_richness_per_clades = richness_list,
                                       simmap = toy_simmap,
                                       # Ask explicitly to retain BD model info
                                       # (useful in Section 5)
                                       return_BD_model_summary = TRUE, 
                                       nb_simulations = 10,
                                       seed = 1, verbose = FALSE)

# Extract sampling info across all 10 grafted phylogenies
sampling_one <- extract_sampling_info(phy_grafted = grafted_sampling,
                                      total_richness_per_clades = richness_list)

# Inspect the first rows describing the results for the first simulations (among 10)
head(sampling_one[, c("simulation", "Clade_name", "total_richness",
                      "nb_backbone_species", "nb_grafted_species", "nb_unsampled_species",
                      "sampling_fraction")])

## ----sampling_one_comparison, fig.width = 7, fig.height = 4.5, out.width = "100%"----
# Compute pre-grafting sampling fraction based on backbone species only 
pre_grafting <- unique(sampling_one[, c("Clade_name", "total_richness", "nb_backbone_species")])
pre_grafting$pre_grafting <- round(pre_grafting$nb_backbone_species / pre_grafting$total_richness, 2)

# Compute mean effective sampling fractions (after grafting) across all 10 simulations
post_grafting <- aggregate(sampling_fraction ~ Clade_name, data = sampling_one,
                           FUN = function(x) round(c(mean = mean(x), min = min(x), max = max(x)), 2))
merge(pre_grafting[, c("Clade_name", "total_richness", "pre_grafting")], post_grafting, by = "Clade_name")

# Plot the distributions of the effective sampling fractions, with the pre-grafting values in red
sampling_one$Clade <- factor(sampling_one$Clade_name, levels = pre_grafting$Clade_name)
boxplot(sampling_fraction ~ Clade, data = sampling_one, ylab = "Sampling fraction", xlab = "Clade", ylim = c(0, 1))
points(seq_len(nrow(pre_grafting)), pre_grafting$pre_grafting, pch = 19, col = "red", cex = 1.5)
legend("topleft", legend = "Before grafting", pch = 19, col = "red", bty = "n")

## ----bd_one-------------------------------------------------------------------
# ------ Birth-death model summaries ------ #

# Export info on BD models
bd_one <- extract_BD_model_summary(grafted_sampling)

# Look at the first rows: one row per model, per simulation
head(bd_one[, c("simulation", "Clade_name", "initial_clade_n_tips", "fitted_clade_n_backbone_tips", "enclosing",
                "b", "d", "net_diversification", "rho")])

# Summarize the variation of each model across the 10 simulations
# Loop across 'Clade_name' to summarize the range of 'b', 'd', and 'rho'
bd_variation <- do.call(rbind, lapply(split(bd_one, bd_one$Clade_name), function(x) {
  data.frame(n_models = nrow(x),
             n_distinct_fits = nrow(unique(round(x[, c("b", "d")], 6))),
             b_range = paste(round(range(x$b), 3), collapse = " - "),
             d_range = paste(round(range(x$d), 3), collapse = " - "),
             rho_range = paste(round(range(x$rho), 2), collapse = " - "))
}))
bd_variation

