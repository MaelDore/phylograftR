##### Build the toy datasets of phylograftR #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Run from the root of the package, after loading it (devtools::load_all() or devtools::install()).
# Produces:
#  * data/toy_backbone.rda, data/toy_missing_taxa_list.rda, data/toy_total_richness_per_clades.rda, data/toy_simmap.rda, data/toy_missing_taxa_list_with_allowed_states.rda,
#    data/toy_simmap_elevation.rda, data/toy_missing_taxa_list_with_two_traits.rda,
#    data/toy_simmap_ranges.rda, data/toy_missing_taxa_list_with_allowed_ranges.rda, data/toy_missing_taxa_list_with_ranges_and_elevation.rda
#  * inst/extdata/toy_missing_taxa.csv, inst/extdata/toy_total_richness.csv, inst/extdata/toy_missing_taxa_with_allowed_states.csv,
#    inst/extdata/toy_missing_taxa_with_two_traits.csv, inst/extdata/toy_missing_taxa_with_allowed_ranges.csv, inst/extdata/toy_missing_taxa_with_ranges_and_elevation.csv
#


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


##### 6/ Create a version of the table of missing taxa with allowed states #####

# Same table as in 2/, with the extra column 'allowed_states': the states of 'toy_simmap' where each taxon can be grafted, separated by commas.
# Blank cells (NA) use the global 'allowed_states' argument of the grafting functions.
# It is a separate dataset because a list with 'allowed_states' requires a 'simmap' to graft.
# The states give 4 taxa with a taxon-level value (one or several states) and 3 taxa using the global argument. On 'toy_simmap',
# the share of the available branch length (and of the time window of the divergence ages) in the states of each taxon is above 10%:
#   Genus_a_sp1 "red": 57% (63%) ; Genus_b_sp1 "blue, green": 96% (100%) ; Genus_b_sp2 "green": 79% (100%) ; Genus_c_sp2 "blue, red": 100% (100%)
# The other taxa have between 31% and 52% (37% to 100%) of their available branches in the states c("blue", "green"), used as global states in the examples.
toy_missing_taxa_with_allowed_states_table <- toy_missing_taxa_table
toy_missing_taxa_with_allowed_states_table$allowed_states <- c("red", NA, NA, "blue, green", "green", NA, "blue, red")

utils::write.csv(toy_missing_taxa_with_allowed_states_table, file = "inst/extdata/toy_missing_taxa_with_allowed_states.csv", row.names = FALSE, na = "")

toy_missing_taxa_list_with_allowed_states <- phylograftR::import_missing_taxa_list(
  file = "inst/extdata/toy_missing_taxa_with_allowed_states.csv", phy = toy_backbone)

# The states are valid for the toy simmap
stopifnot(all(unlist(toy_missing_taxa_list_with_allowed_states$allowed_states) %in% toy_simmap_states))


##### 6b/ Create a second trait (simmap) and a table of missing taxa with allowed states of two traits #####

# Second trait: 'toy_simmap_elevation', two states "high" and "low" simulated from a root in "low" with equal rates (0.04), on the toy backbone.
# The map is used with 'toy_simmap' (the first trait: color) as the maps of two traits: simmap = list(color = toy_simmap, elevation = toy_simmap_elevation).
# The seed is the first one (among 1 to 1000) for which
#  * the map is balanced: each state covers between 42% and 58% of the total length of the branches, 15% to 30% of the branches include a change of state,
#    no branch has more than 4 pieces, and each state is found in 10 to 20 of the 30 tips,
#  * and the grafting of the table below, with allowed_states = list(c("blue", "green"), NULL), runs without error or warning with both methods ("uniform" and "BD_clades").

toy_simmap_elevation_states <- c("high", "low")
toy_simmap_elevation_Q <- matrix(toy_simmap_rate, nrow = 2, ncol = 2, dimnames = list(toy_simmap_elevation_states, toy_simmap_elevation_states))
diag(toy_simmap_elevation_Q) <- -toy_simmap_rate

# Function to produce a candidate simmap
simulate_toy_simmap_elevation <- function (seed)
{
  set.seed(seed)
  phytools::sim.history(tree = toy_backbone, Q = toy_simmap_elevation_Q, anc = "low", message = FALSE)
}
# Function to evaluate candidate based on the balance of states along branches and tips
is_balanced_simmap_elevation <- function (simmap)
{
  duration_per_state <- colSums(simmap$mapped.edge)
  if (!all(toy_simmap_elevation_states %in% names(duration_per_state))) { return(FALSE) }
  share_of_length <- duration_per_state[toy_simmap_elevation_states] / sum(duration_per_state)
  share_of_branches_with_change <- mean(lengths(simmap$maps) > 1)
  tips_per_state <- table(factor(simmap$states, levels = toy_simmap_elevation_states))
  all(abs(share_of_length - 1/2) < 0.08) &&
    share_of_branches_with_change >= 0.15 && share_of_branches_with_change <= 0.30 && max(lengths(simmap$maps)) <= 4 &&
    all(tips_per_state >= 10) && all(tips_per_state <= 20)
}

# Table of missing taxa: same as in 2/, with the column 'allowed_states' giving the states of the two traits (color, then elevation) separated by semicolons,
# the states of a trait separated by commas. Blank states have no taxon-level value for the trait (the global argument is used). Traits can be named (e.g., "elevation=high").
# Different syntaxes are used on purpose to showcase them. Blank cells (NA) use the global 'allowed_states' for both traits.

# Rewrite toy_missing_taxa_table to adjust allowed_states for two traits
toy_missing_taxa_with_two_traits_table <- toy_missing_taxa_table
toy_missing_taxa_with_two_traits_table$allowed_states <- c("red;", NA, NA, "blue, green;high", ";low", NA, "elevation=high")
utils::write.csv(toy_missing_taxa_with_two_traits_table, file = "inst/extdata/toy_missing_taxa_with_two_traits.csv", row.names = FALSE, na = "")

# Import from the .csv
toy_missing_taxa_list_with_two_traits <- phylograftR::import_missing_taxa_list(
  file = "inst/extdata/toy_missing_taxa_with_two_traits.csv", phy = toy_backbone)

# Run grafting
runs_without_warning <- function (simmap)
{
  all(vapply(c("uniform", "BD_clades"), function (method)
  {
    result <- tryCatch(withCallingHandlers(
      phylograftR::graft_missing_taxa(toy_backbone, toy_missing_taxa_list_with_two_traits, simmap = list(color = toy_simmap, elevation = simmap),
                                      allowed_states = list(c("blue", "green"), NULL), method = method, seed = 1234, verbose = FALSE),
      warning = function (w) stop("warning")), error = function (e) NULL)
    !is.null(result)
  }, logical(1)))
}

# Look for a good candidate for the simmap of elevation
toy_simmap_elevation_seed <- NULL
for (seed in 1:1000)
{
  simmap <- simulate_toy_simmap_elevation(seed)
  if (is_balanced_simmap_elevation(simmap) && runs_without_warning(simmap)) { toy_simmap_elevation_seed <- seed; break }
}
stopifnot(!is.null(toy_simmap_elevation_seed))
toy_simmap_elevation <- simulate_toy_simmap_elevation(toy_simmap_elevation_seed)
toy_simmap_elevation_seed # 94


##### 7/ Create a stochastic character map (simmap) of the backbone with biogeographic ranges #####

# Four ranges, coded as in BioGeoBEARS: "A", "B", "AB" and "BC" ------------------------------------------------
# Areas A, B and C. Range "C" alone is not a state, so "C" is only found within "BC".

# Anagenetic changes along the branches ------------------------------------------------------------------------
# Only single-step changes of range between the states, i.e., the gain (dispersal, "d" in BioGeoBEARS) or the loss (local extinction, "e") of a single area:
#   "d": "A" -> "AB", "B" -> "AB" and "B" -> "BC"; "e": "AB" -> "A", "AB" -> "B" and "BC" -> "B".
# A change of two areas at once is not allowed ("AB" <-> "BC"), so the states are connected as a chain: A - AB - B - BC.
# All the allowed changes have the same rate (0.04 per unit of time, as for 'toy_simmap').
toy_simmap_ranges_states <- c("A", "B", "AB", "BC")
toy_simmap_ranges_Q <- matrix(0, nrow = 4, ncol = 4, dimnames = list(toy_simmap_ranges_states, toy_simmap_ranges_states))
for (change in list(c("A", "AB"), c("B", "AB"), c("B", "BC")))
{
  toy_simmap_ranges_Q[change[1], change[2]] <- toy_simmap_rate
  toy_simmap_ranges_Q[change[2], change[1]] <- toy_simmap_rate
}
diag(toy_simmap_ranges_Q) <- -rowSums(toy_simmap_ranges_Q)

# Helpers on ranges and events ------------------------------------------------------------------------------------
split_areas <- function (range) { strsplit(range, split = "", fixed = TRUE)[[1]] }
is_loss_of_one_area <- function (from, to) { length(split_areas(to)) == length(split_areas(from)) - 1 && all(split_areas(to) %in% split_areas(from)) }
n_toy_tips <- ape::Ntip(toy_backbone)
toy_nodes <- (n_toy_tips + 1):(n_toy_tips + toy_backbone$Nnode)

# Range of the lineage at the end of the branch above a node (the root range for the root)
range_above_node <- function (simmap, node)
{
  parent_edge <- match(node, simmap$edge[, 2])
  if (is.na(parent_edge)) { simmap$node.states[1, 1] } else { utils::tail(names(simmap$maps[[parent_edge]]), 1) }
}

# Cladogenetic events: move part of the changes to the nodes -------------------------------------------------------
# A simulated history only has anagenetic changes, that happen anywhere along the branches. To also represent cladogenetic events, where the change
# is linked to the speciation event, some changes are moved to the parental node: the first piece of the branch is removed and the new range starts at the node.
# Only the first change of a branch that is the loss of one area can be moved, as it gives a cladogenetic event of BioGeoBEARS:
#   * if one of the two descendant branches is moved, e.g., (AB) -> (AB), (A): subset sympatry ("s");
#   * if both are moved and the descendant ranges split the ancestral range, e.g., (AB) -> (A), (B): vicariance ("v").
# The gain of an area cannot be moved: (A) -> (A), (AB) is not an event of BioGeoBEARS. Jump dispersal ("j") is not represented.
# The nodes are taken in a random order, and changes are moved until about 'target_share' of the events (anagenetic changes and cladogenetic events) are cladogenetic.

move_changes_to_nodes <- function (simmap, seed, target_share = 0.3)
{
  force(simmap)
  set.seed(seed)
  n_anagenetic <- sum(lengths(simmap$maps) - 1)
  n_cladogenetic <- 0
  for (node in sample(toy_nodes))
  {
    if (n_cladogenetic / (n_cladogenetic + n_anagenetic) >= target_share) { break }
    node_range <- range_above_node(simmap, node)
    child_edges <- which(simmap$edge[, 1] == node)
    movable <- child_edges[vapply(child_edges, function (edge)
    {
      pieces <- names(simmap$maps[[edge]])
      length(pieces) >= 2 && pieces[1] == node_range && is_loss_of_one_area(pieces[1], pieces[2])
    }, logical(1))]
    if (length(movable) == 0) { next }
    new_ranges <- vapply(movable, function (edge) names(simmap$maps[[edge]])[2], character(1))
    vicariance <- length(movable) == 2 && new_ranges[1] != new_ranges[2] && setequal(unlist(strsplit(new_ranges, "")), split_areas(node_range)) && runif(1) < 0.5
    moved <- if (vicariance) movable else movable[sample.int(length(movable), 1)]
    for (edge in moved)
    {
      pieces <- simmap$maps[[edge]]
      simmap$maps[[edge]] <- stats::setNames(c(pieces[1] + pieces[2], pieces[-(1:2)]), names(pieces)[-1])
      simmap$node.states[edge, 1] <- names(pieces)[2]
    }
    n_anagenetic <- n_anagenetic - length(moved)
    n_cladogenetic <- n_cladogenetic + 1
  }
  # Update the total duration per range and branch
  simmap$mapped.edge[] <- 0
  for (edge in seq_along(simmap$maps))
  {
    for (range in unique(names(simmap$maps[[edge]]))) { simmap$mapped.edge[edge, range] <- sum(simmap$maps[[edge]][names(simmap$maps[[edge]]) == range]) }
  }
  simmap
}

# Count the events of a simmap of ranges, with the BioGeoBEARS syntax -----------------------------------------------
#  * anagenetic (changes within branches): "d" (gain of one area), "e" (loss of one area), "other" (anything else)
#  * cladogenetic (at the nodes, from the range above the node and the ranges at the start of the two descendant branches):
#    "y" (narrow sympatry: single-area range copied), "range_copying" (widespread range copied by both descendants, e.g., (AB) -> (AB), (AB)),
#    "s" (subset sympatry), "v" (vicariance), "other" (anything else)

count_toy_events <- function (simmap)
{
  changes <- do.call(rbind, lapply(simmap$maps, function (map)
  {
    pieces <- names(map)
    if (length(pieces) < 2) { NULL } else { data.frame(from = utils::head(pieces, -1), to = utils::tail(pieces, -1), stringsAsFactors = FALSE) }
  }))
  anagenetic <- ifelse(mapply(is_loss_of_one_area, changes$from, changes$to), "e", ifelse(mapply(is_loss_of_one_area, changes$to, changes$from), "d", "other"))
  cladogenetic <- vapply(toy_nodes, function (node)
  {
    node_range <- range_above_node(simmap, node)
    descendant_ranges <- vapply(which(simmap$edge[, 1] == node), function (edge) names(simmap$maps[[edge]])[1], character(1))
    if (all(descendant_ranges == node_range)) { return(if (nchar(node_range) == 1) "y" else "range_copying") }
    if (sum(descendant_ranges == node_range) == 1)
    {
      shifted <- descendant_ranges[descendant_ranges != node_range]
      return(if (nchar(node_range) > 1 && nchar(shifted) == 1 && shifted %in% split_areas(node_range)) "s" else "other")
    }
    areas_1 <- split_areas(descendant_ranges[1])
    areas_2 <- split_areas(descendant_ranges[2])
    if (nchar(node_range) > 1 && !any(areas_1 %in% areas_2) && setequal(c(areas_1, areas_2), split_areas(node_range))) { return("v") }
    "other"
  }, character(1))
  list(anagenetic = table(factor(anagenetic, levels = c("d", "e", "other"))),
       cladogenetic = table(factor(cladogenetic, levels = c("y", "range_copying", "s", "v", "other"))))
}

# Simulation and choice of the seed --------------------------------------------------------------------------------
# A single history is simulated from a root in range "AB" (phytools::sim.history()), then changes are moved to the nodes (see above).
# The seed (the same for the two steps) is the first among 1 to 20000 giving a balanced map:
#  * each range covers between 18% and 32% of the total length of the branches,
#  * between 15% and 30% of the branches include a change of range, and no branch has more than 4 pieces,
#  * each range is found in 5 to 10 of the 30 tips,
#  * all the events are of a type of BioGeoBEARS (no "other"), with at least 2 subset sympatry events and 1 vicariance event,
#  * between 25% and 36% of the events (anagenetic changes and cladogenetic events "s" and "v") are cladogenetic.

simulate_toy_simmap_ranges <- function (seed)
{
  set.seed(seed)
  simmap <- phytools::sim.history(tree = toy_backbone, Q = toy_simmap_ranges_Q, anc = "AB", message = FALSE)
  move_changes_to_nodes(simmap, seed = seed, target_share = 0.3)
}

is_balanced_simmap_ranges <- function (simmap)
{
  duration_per_state <- colSums(simmap$mapped.edge)
  if (!all(toy_simmap_ranges_states %in% names(duration_per_state))) { return(FALSE) }
  share_of_length <- duration_per_state[toy_simmap_ranges_states] / sum(duration_per_state)
  share_of_branches_with_change <- mean(lengths(simmap$maps) > 1)
  tips_per_state <- table(factor(simmap$states, levels = toy_simmap_ranges_states))
  events <- count_toy_events(simmap)
  n_events <- sum(events$anagenetic[c("d", "e")]) + sum(events$cladogenetic[c("s", "v")])
  share_cladogenetic <- sum(events$cladogenetic[c("s", "v")]) / n_events
  all(share_of_length >= 0.18 & share_of_length <= 0.32) &&
    share_of_branches_with_change >= 0.15 && share_of_branches_with_change <= 0.30 && max(lengths(simmap$maps)) <= 4 &&
    all(tips_per_state >= 5 & tips_per_state <= 10) &&
    events$anagenetic["other"] == 0 && events$cladogenetic["other"] == 0 && events$cladogenetic["s"] >= 2 && events$cladogenetic["v"] >= 1 &&
    share_cladogenetic >= 0.25 && share_cladogenetic <= 0.36
}

toy_simmap_ranges_seed <- NULL
for (seed in 1:20000)
{
  if (is_balanced_simmap_ranges(simulate_toy_simmap_ranges(seed))) { toy_simmap_ranges_seed <- seed; break }
}
stopifnot(!is.null(toy_simmap_ranges_seed))
toy_simmap_ranges <- simulate_toy_simmap_ranges(toy_simmap_ranges_seed)

# Summary of the map
toy_simmap_ranges_seed # 3598
round(colSums(toy_simmap_ranges$mapped.edge) / sum(toy_simmap_ranges$mapped.edge), 2) # Share of the total branch length in each range
table(toy_simmap_ranges$states) # Ranges of the tips
count_toy_events(toy_simmap_ranges) # Events: anagenetic d = 8, e = 4; cladogenetic y = 13, range_copying = 10, s = 5, v = 1


##### 8/ Create a version of the table of missing taxa with allowed ranges #####

# Same table as in 2/, with the extra column 'allowed_states': the ranges of 'toy_simmap_ranges' where each taxon can be grafted, separated by commas.
# Used with 'treat_states_as_ranges = TRUE': an allowed range matches the branch pieces in a range that contains it, and the grafted lineage inherits the allowed range.
# All the allowed ranges are states of 'toy_simmap_ranges' (so no warning is sent), and some are contained in other states:
# "A" is also found in "AB", "B" in "AB" and "BC". Grafting a lineage in "A" on a branch in "AB" is a subset sympatric speciation event, (AB) -> (AB), (A).
# Blank cells (NA) use the global 'allowed_states' argument of the grafting functions (c("A", "B") in the examples).
# The states give 4 taxa with a taxon-level value (one or several ranges) and 3 taxa using the global argument.

toy_missing_taxa_with_allowed_ranges_table <- toy_missing_taxa_table
toy_missing_taxa_with_allowed_ranges_table$allowed_states <- c("A", NA, NA, "B, BC", "B", NA, "A, BC")

utils::write.csv(toy_missing_taxa_with_allowed_ranges_table, file = "inst/extdata/toy_missing_taxa_with_allowed_ranges.csv", row.names = FALSE, na = "")

toy_missing_taxa_list_with_allowed_ranges <- phylograftR::import_missing_taxa_list(
  file = "inst/extdata/toy_missing_taxa_with_allowed_ranges.csv", phy = toy_backbone)

# The ranges are states of the toy simmap
stopifnot(all(unlist(toy_missing_taxa_list_with_allowed_ranges$allowed_states) %in% toy_simmap_ranges_states))


##### 8b/ Create a table of missing taxa with allowed ranges and allowed states of elevation #####

# Same table as in 2/, with the column 'allowed_states' giving the allowed values of two traits, separated by semicolons: the biogeographic ranges of 'toy_simmap_ranges', then the states of 'toy_simmap_elevation'.
# The two maps are used as the maps of two traits, with the first one treated as ranges:
#   simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation), treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE)
# Different syntaxes are used on purpose to showcase them (positional, named, empty trait). Blank cells (NA) use the global 'allowed_states' argument (e.g., list(c("A", "B"), NULL) in the examples).
# The allowed ranges are matched by containment ("A" is also found in "AB"), the elevation states are matched exactly.
# The cells were chosen so that the grafting runs without error or warning with both methods ("uniform" and "BD_clades") and the global argument list(c("A", "B"), NULL):
#   Genus_a_sp1 "A;high" ; Genus_b_sp1 "B, BC;" (any elevation) ; Genus_b_sp2 "elevation=low" (any range) ; Genus_c_sp2 "biogeo=A, BC; elevation=high"

toy_missing_taxa_with_ranges_and_elevation_table <- toy_missing_taxa_table
toy_missing_taxa_with_ranges_and_elevation_table$allowed_states <- c("A;high", NA, NA, "B, BC;", "elevation=low", NA, "biogeo=A, BC; elevation=high")
utils::write.csv(toy_missing_taxa_with_ranges_and_elevation_table, file = "inst/extdata/toy_missing_taxa_with_ranges_and_elevation.csv", row.names = FALSE, na = "")

toy_missing_taxa_list_with_ranges_and_elevation <- phylograftR::import_missing_taxa_list(
  file = "inst/extdata/toy_missing_taxa_with_ranges_and_elevation.csv", phy = toy_backbone)

# The grafting runs without error or warning with both methods
for (method in c("uniform", "BD_clades"))
{
  withCallingHandlers(
    phylograftR::graft_missing_taxa(toy_backbone, toy_missing_taxa_list_with_ranges_and_elevation, simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
                                    allowed_states = list(c("A", "B"), NULL), treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE), method = method, seed = 1234, verbose = FALSE),
    warning = function (w) stop("Unexpected warning: ", conditionMessage(w)))
}


##### 9/ Save objects #####

dir.create("data", showWarnings = FALSE)
save(toy_backbone, file = "data/toy_backbone.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list, file = "data/toy_missing_taxa_list.rda", compress = "gzip", version = 2)
save(toy_total_richness_per_clades, file = "data/toy_total_richness_per_clades.rda", compress = "gzip", version = 2)
save(toy_simmap, file = "data/toy_simmap.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list_with_allowed_states, file = "data/toy_missing_taxa_list_with_allowed_states.rda", compress = "gzip", version = 2)
save(toy_simmap_elevation, file = "data/toy_simmap_elevation.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list_with_two_traits, file = "data/toy_missing_taxa_list_with_two_traits.rda", compress = "gzip", version = 2)
save(toy_simmap_ranges, file = "data/toy_simmap_ranges.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list_with_allowed_ranges, file = "data/toy_missing_taxa_list_with_allowed_ranges.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list_with_ranges_and_elevation, file = "data/toy_missing_taxa_list_with_ranges_and_elevation.rda", compress = "gzip", version = 2)
