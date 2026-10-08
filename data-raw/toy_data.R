##### Build the toy datasets of phylograftR #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Run from the root of the package, after loading it (devtools::load_all() or devtools::install()).
# Produces:
#  * data/toy_backbone.rda, data/toy_missing_taxa_list.rda, data/toy_total_richness_per_clades.rda, data/toy_simmap.rda
#  * inst/extdata/toy_missing_taxa.csv, inst/extdata/toy_total_richness.csv


##### 1/ Create backbone phylogeny #####

# Backbone: 30 tips, ultrametric, root age of 50 -----------------------------------------------

set.seed(1234)
toy_backbone <- ape::rcoal(30)
toy_backbone$edge.length <- toy_backbone$edge.length / max(ape::node.depth.edgelength(toy_backbone)) * 50

# Minimum length of the internal branches -----------------------------------------------------------------

# Coalescent trees nearly always include very short internal branches (here, two of 0.012 and 0.027 for a root age of 50),
# that look like polytomies on a plot. The package does not handle polytomies, so an example must not suggest it does.
# The nodes are moved towards the root where needed, so every internal branch is at least 'min_internal_branch' long.
# The topology (tips, clades and node IDs), the root age and the ultrametry are unchanged.
min_internal_branch <- 1 # 2% of the root age
n_tips <- ape::Ntip(toy_backbone)
node_ages <- max(ape::node.depth.edgelength(toy_backbone)) - ape::node.depth.edgelength(toy_backbone) # Age of every node (tips: 0)
edges_postorder <- toy_backbone$edge[ape::postorder(toy_backbone), , drop = FALSE] # Branches of the tips towards the root
for (i in seq_len(nrow(edges_postorder)))
{
  parent <- edges_postorder[i, 1]
  child <- edges_postorder[i, 2]
  if (child > n_tips) { node_ages[parent] <- max(node_ages[parent], node_ages[child] + min_internal_branch) }
}
toy_backbone$edge.length <- node_ages[toy_backbone$edge[, 1]] - node_ages[toy_backbone$edge[, 2]]

# Check that the root age and the ultrametry are unchanged, and that no internal branch is shorter than the minimum
stopifnot(isTRUE(all.equal(max(ape::node.depth.edgelength(toy_backbone)), 50)),
          ape::is.ultrametric(toy_backbone),
          min(toy_backbone$edge.length[toy_backbone$edge[, 2] > n_tips]) >= min_internal_branch - 1e-8)

# Reference clades (node: tips) used below
# 32: t26 t2 t21 t16 t11 t4 t3 t25 t28 t10 t6 t20 t9 t30 t24 t18 t19  (clade A)
# 33: t8 t1 t15 t29 t14 t13 t17 t27 t7 t5 t23 t12 t22                 (clade B)
# 34: t8 t1 t15 t29 t14 t13 t17 t27 t7                                (clade B1, nested in B)
# 37: t10 t6 t20 t9 t30 t24 t18 t19                                   (group_1, nested in A)
# 39: t26 t2 t21 t16 t11 t4                                           (group_2, nested in A)
# 50: t13 t17 t27                                                     (group_3, nested in B1)
# 53: t5 t23 t12                                                      (clade C, nested in B)

##### 2/ Create table of the missing taxa as .csv #####

toy_missing_taxa_table <- data.frame(
  Taxon_name = c("Genus_a_sp1", "Genus_a_sp2", "Genus_a_sp3", "Genus_b_sp1", "Genus_b_sp2", "Genus_c_sp1", "Genus_c_sp2"),
  Target_clade = c("t26,t18", "t26,t18", "t26,t18", "t8,t22", "t8,t7", "t5", "t5,t12"),
  Target_clade_name = c("clade_A", "clade_A", "clade_A", "clade_B", "clade_B1", "tip_t5", "clade_C"),
  include_target_clade_stem = c(TRUE, TRUE, FALSE, NA, TRUE, NA, FALSE),
  Excluded_subclades = c("t10,t19;t26,t4", "t26,t4;t30", NA, NA, "t13,t27", NA, NA),
  Excluded_subclades_names = c("group_1;group_2", "group_2;tip_t30", NA, NA, "group_3", NA, NA),
  exclude_subclades_stem = c("TRUE;FALSE", "FALSE;TRUE", NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE)

# Export in inst/extdata/
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(toy_missing_taxa_table, file = "inst/extdata/toy_missing_taxa.csv", row.names = FALSE, na = "")

##### 3/ Convert to missing_taxa_list object #####

toy_missing_taxa_list <- phylograftR::import_missing_taxa_list(
  file = "inst/extdata/toy_missing_taxa.csv", phy = toy_backbone)

##### 4/ Create table of Total richness of the target clades as .csv #####

toy_total_richness_table <- data.frame(
  Target_clade = c("t26,t18", "t8,t22", "t8,t7", "t5", "t5,t12"),
  total_richness = c(25, NA, 14, NA, 8))

# Export in inst/extdata
utils::write.csv(toy_total_richness_table, file = "inst/extdata/toy_total_richness.csv", row.names = FALSE, na = "")

toy_total_richness_per_clades <- list(
  Target_clade = strsplit(toy_total_richness_table$Target_clade, ",", fixed = TRUE),
  total_richness = toy_total_richness_table$total_richness)

##### 5/ Create a stochastic character map (simmap) of the backbone #####

# Three states named after colors, to ease their visualization: "blue", "green" and "red" -------------------

# Equal transition rates between the three states (symmetric Q matrix), so the three states are expected in similar amounts.
# The rate (0.04 per unit of time between each pair of states, i.e., one transition every 12.5 time units on average on a lineage,
# for a root age of 50) gives a moderate number of transitions: about one branch out of five includes a change of state.
toy_simmap_states <- c("blue", "green", "red")
toy_simmap_rate <- 0.04
toy_simmap_Q <- matrix(toy_simmap_rate, nrow = 3, ncol = 3, dimnames = list(toy_simmap_states, toy_simmap_states))
diag(toy_simmap_Q) <- -2 * toy_simmap_rate

# A single history is simulated from a root in state "green" (phytools::sim.history()), but the seed is chosen so that the map is balanced:
# the first seed (among 1 to 1000) for which
#  * each state covers between 29% and 38% of the total length of the branches,
#  * between 15% and 30% of the branches include a change of state, and no branch has more than 4 pieces,
#  * each state is found in 7 to 13 of the 30 tips.
simulate_toy_simmap <- function (seed)
{
  set.seed(seed)
  phytools::sim.history(tree = toy_backbone, Q = toy_simmap_Q, anc = "green", message = FALSE)
}

is_balanced_simmap <- function (simmap)
{
  duration_per_state <- colSums(simmap$mapped.edge)
  if (!all(toy_simmap_states %in% names(duration_per_state))) { return(FALSE) }
  share_of_length <- duration_per_state[toy_simmap_states] / sum(duration_per_state)
  share_of_branches_with_change <- mean(lengths(simmap$maps) > 1)
  tips_per_state <- table(factor(simmap$states, levels = toy_simmap_states))
  all(abs(share_of_length - 1/3) < 0.04) &&
    share_of_branches_with_change >= 0.15 && share_of_branches_with_change <= 0.30 && max(lengths(simmap$maps)) <= 4 &&
    all(tips_per_state >= 7) && all(tips_per_state <= 13)
}

toy_simmap_seed <- NULL
for (seed in 1:1000)
{
  if (is_balanced_simmap(simulate_toy_simmap(seed))) { toy_simmap_seed <- seed; break }
}
stopifnot(!is.null(toy_simmap_seed))
toy_simmap <- simulate_toy_simmap(toy_simmap_seed)

# Summary of the map
toy_simmap_seed # 50
round(colSums(toy_simmap$mapped.edge) / sum(toy_simmap$mapped.edge), 2) # Share of the total branch length in each state
sum(lengths(toy_simmap$maps) - 1) # Number of changes of state along the branches
table(toy_simmap$states) # States of the tips


##### 6/ Save objects #####

dir.create("data", showWarnings = FALSE)
save(toy_backbone, file = "data/toy_backbone.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list, file = "data/toy_missing_taxa_list.rda", compress = "gzip", version = 2)
save(toy_total_richness_per_clades, file = "data/toy_total_richness_per_clades.rda", compress = "gzip", version = 2)
save(toy_simmap, file = "data/toy_simmap.rda", compress = "gzip", version = 2)
