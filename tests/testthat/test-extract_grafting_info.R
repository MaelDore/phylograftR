
##### Test suite for extract_grafting_info() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Name the items of a list of missing taxa after the taxa, as required by the functions
info_mtl <- function (l)
{
  for (item in setdiff(names(l), "Taxon_name")) { if (is.null(names(l[[item]]))) { names(l[[item]]) <- l$Taxon_name } }
  l
}
info_multi_simmap <- function (maps) { structure(maps, class = c("multiSimmap", "multiPhylo")) }

## Node ID, in a tree, of the MRCA of tips (or of the tip)
info_node <- function (phy, tips) { if (length(tips) == 1) match(tips, phy$tip.label) else ape::getMRCA(phy, tips) }
## Is a pair of tips defining the same clade as the node of the backbone?
info_same_clade <- function (phy, parent_tips, node) { identical(as.numeric(info_node(phy, strsplit(parent_tips, ", ")[[1]])), as.numeric(node)) }


test_that("extract_grafting_info() gives the graft age, divergence age and parent clade of taxa grafted on the backbone",
{
  toy <- toy_setup()
  # Taxa grafted by hand on known branches, at known positions
  g <- graft_tip_on_branch(toy$phy, "X1", toy$node_B, frac = 0.5)
  g <- graft_tip_on_branch(g, "X2", match(toy$s, g$tip.label), frac = 0.25)
  mtl <- info_mtl(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$B, toy$s)))
  info <- extract_grafting_info(g, mtl, phy_backbone = toy$phy)

  expect_s3_class(info, "data.frame")
  expect_equal(info$Taxon_name, c("X1", "X2"))
  expect_equal(info$simulation, c(1, 1))
  expect_true(all(is.na(info$map)))
  expect_equal(attr(info, "type"), "taxa")
  # Graft ages: independent computation of the age of the node where each taxon is attached
  expect_equal(info$graft_age, c(graft_age(g, "X1"), graft_age(g, "X2")))
  expect_equal(info$divergence_age, info$graft_age)
  expect_true(all(is.na(info$grafted_on)))
  # Parent clades: X1 on the branch of the clade B, X2 on the terminal branch of 's'
  expect_equal(info$parent_node, c(toy$node_B, match(toy$s, toy$phy$tip.label)))
  expect_true(info_same_clade(toy$phy, info$parent_tips[1], toy$node_B))
  expect_equal(info$parent_tips[2], toy$s)
  # Without the backbone: rebuilt by removing the grafted taxa, the node IDs are not given
  info_no_backbone <- extract_grafting_info(g, mtl)
  expect_true(all(is.na(info_no_backbone$parent_node)))
  expect_equal(info_no_backbone[, setdiff(names(info), "parent_node")], info[, setdiff(names(info), "parent_node")])
  # The age is the one set by the grafting position: crown age of the branch + the position from its tipward node
  root_age <- max(phytools::nodeHeights(toy$phy)[, 2])
  branch_length <- function (node) toy$phy$edge.length[match(node, toy$phy$edge[, 2])]
  expect_equal(info$graft_age[1], root_age - phytools::nodeheight(toy$phy, toy$node_B) + 0.5 * branch_length(toy$node_B))
  expect_equal(info$graft_age[2], 0.25 * branch_length(match(toy$s, toy$phy$tip.label)))
})

test_that("extract_grafting_info() reports the graft point of taxa grafted on other grafted taxa, with the graft point of their grafted subtree",
{
  toy <- toy_setup()
  # X1 is grafted on the branch of 's', then X2 on the terminal branch of X1, then X3 on the terminal branch of X2
  g <- graft_tip_on_branch(toy$phy, "X1", match(toy$s, toy$phy$tip.label), frac = 0.5)
  g <- graft_tip_on_branch(g, "X2", match("X1", g$tip.label), frac = 0.5)
  g <- graft_tip_on_branch(g, "X3", match("X2", g$tip.label), frac = 0.5)
  mtl <- info_mtl(list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(toy$s, toy$s, toy$s)))
  info <- extract_grafting_info(g, mtl, phy_backbone = toy$phy)

  s_node <- match(toy$s, toy$phy$tip.label)
  # All the taxa belong to the same grafted subtree, attached to the terminal branch of 's' at the graft point of X1
  expect_equal(info$parent_node, rep(s_node, 3))
  expect_equal(info$parent_tips, rep(toy$s, 3))
  expect_equal(info$graft_age, rep(graft_age(g, c("X1", "X2", "X3")), 3))
  # Divergence ages: the age of the node where each taxon splits from its sister lineage (nodes get younger along the chain)
  expect_equal(info$divergence_age[3], max(phytools::nodeHeights(g)[, 2]) - phytools::nodeheight(g, ape::getMRCA(g, c("X2", "X3"))))
  expect_true(all(info$divergence_age <= info$graft_age + 1e-9))
  # The sister lineage of X3 is X2, and the sister lineage of X2 is X3 (the cherry)
  expect_equal(info$divergence_age[2], info$divergence_age[3])
  expect_equal(info$grafted_on[2], "X3")
  expect_equal(info$grafted_on[3], "X2")
  # X1 is attached to the cherry (X2, X3): its sister lineage is made of grafted taxa only, and it splits from it at the node below the cherry
  expect_equal(info$grafted_on[1], "X2, X3")
  expect_equal(info$divergence_age[1], graft_age(g, c("X2", "X3")))
  expect_gt(info$graft_age[1], info$divergence_age[1])
  expect_gt(info$divergence_age[1], info$divergence_age[2])
})

test_that("extract_grafting_info() extracts the states at the graft point and inherited by the taxa from the output simmap",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- info_mtl(list(Taxon_name = c("X1", "X2", "X3", "X4"), Target_clade = list(toy$A, toy$B, toy$C, toy$s)))
  g <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = simmap)
  info <- extract_grafting_info(g, mtl, toy$phy)
  states <- graft_states_in_simmap(g, mtl$Taxon_name)

  expect_equal(info$tip_state, states$state)
  expect_equal(info$parent_state, states$parent_state)
  # Independent computation of the state at the graft point from the backbone simmap
  expected <- vapply(mtl$Taxon_name, function (x) state_at_graft(g, simmap, x, mtl$Taxon_name), character(1))
  ok <- !is.na(expected)
  expect_equal(unname(info$parent_state[ok]), unname(expected[ok]))
  # Without states: no state columns
  g_no <- drop_simmap(g)
  expect_false(any(c("parent_state", "tip_state") %in% names(extract_grafting_info(g_no, mtl, toy$phy))))
  # An output simmap that does not match the phylogeny gives a warning and NA
  g_bad <- ape::drop.tip(g, "X1"); g_bad$simmap <- g$simmap
  expect_warning(bad <- extract_grafting_info(g_bad, info_mtl(list(Taxon_name = c("X2", "X3", "X4"), Target_clade = mtl$Target_clade[2:4]))), "does not match")
  expect_true(all(is.na(bad$parent_state)))
})

test_that("extract_grafting_info() handles ranges (inherited state differs from the state of the graft point) and several traits",
{
  toy <- toy_setup()
  ranges <- make_toy_range_simmap(toy$phy)
  mtl <- info_mtl(list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(toy$A, toy$B, toy$C), allowed_states = list("A", "A", "A")))
  g <- suppressWarnings(quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = ranges, treat_states_as_ranges = TRUE))
  info <- extract_grafting_info(g, mtl, toy$phy)
  states <- graft_states_in_simmap(g, mtl$Taxon_name)
  expect_equal(info$tip_state, states$state)
  expect_equal(info$parent_state, states$parent_state)
  # The allowed range is "A": the inherited range is "A", the state of the graft point contains it
  expect_true(all(info$tip_state == "A"))
  expect_true(all(range_in_states_pairwise(info$tip_state, info$parent_state)))

  # Several traits: one column per trait
  traits <- list(biogeo = ranges, elev = make_toy_simmap(toy$phy, seed = 2))
  mtl_t <- info_mtl(list(Taxon_name = c("X1", "X2"), Target_clade = list(toy$A, toy$B)))
  g_t <- suppressWarnings(quiet_graft(graft_missing_taxa, toy$phy, mtl_t, simmap = traits))
  info_t <- extract_grafting_info(g_t, mtl_t, toy$phy)
  expect_true(all(c("parent_state_biogeo", "tip_state_biogeo", "parent_state_elev", "tip_state_elev") %in% names(info_t)))
  expect_false("parent_state" %in% names(info_t))
  expect_equal(info_t$tip_state_elev, vapply(mtl_t$Taxon_name, function (x) names(g_t$simmap$elev$maps[[match(match(x, g_t$tip.label), g_t$edge[, 2])]])[1], character(1), USE.NAMES = FALSE))
})

test_that("extract_grafting_info() gives a long table over simulations, maps and a list of maps",
{
  toy <- toy_setup()
  mtl <- info_mtl(list(Taxon_name = c("X1", "X2", "X3"), Target_clade = list(toy$A, toy$B, toy$s)))
  sims <- quiet_graft(graft_missing_taxa, toy$phy, mtl, nb_simulations = 3)
  info <- extract_grafting_info(sims, mtl, toy$phy)
  expect_equal(nrow(info), 9)
  expect_equal(info$simulation, rep(1:3, each = 3))
  expect_equal(info$seed, rep(unname(vapply(sims, function (p) attr(p, "seed"), numeric(1))), each = 3))
  # Each phylogeny gives the same rows as the extraction of the phylogeny alone
  alone <- extract_grafting_info(sims[[2]], mtl, toy$phy)
  expect_equal(info$graft_age[4:6], alone$graft_age)
  expect_equal(info$parent_tips[4:6], alone$parent_tips)
  expect_equal(info$Taxon_name, rep(mtl$Taxon_name, 3))
  # X2 is grafted in the clade B: its graft age cannot be older than the stem age of B
  stem_age_B <- max(phytools::nodeHeights(toy$phy)[, 2]) - phytools::nodeHeights(toy$phy)[match(toy$node_B, toy$phy$edge[, 2]), 1]
  expect_true(all(info$graft_age[info$Taxon_name == "X2"] <= stem_age_B + 1e-9))

  # List of maps (multiSimmap): one long table with the map name
  maps <- lapply(1:2, function (i) make_toy_simmap(toy$phy, seed = i))
  g_maps <- quiet_graft(graft_missing_taxa, toy$phy, mtl, simmap = info_multi_simmap(maps), nb_simulations = 2)
  info_maps <- extract_grafting_info(g_maps, mtl, toy$phy)
  expect_equal(nrow(info_maps), 12)
  expect_equal(unique(info_maps$map), c("map_1", "map_2"))
  expect_equal(info_maps$simulation, rep(rep(1:2, each = 3), 2))
  expect_true(all(c("parent_state", "tip_state") %in% names(info_maps)))
})

test_that("extract_grafting_info() with a 'Sister_clade' list gives one row per missing clade",
{
  toy <- toy_setup()
  simmap <- make_toy_simmap(toy$phy)
  mtl <- info_mtl(list(Taxon_name = c("M1", "M2", "M3"), Sister_clade = list(toy$B, toy$B, toy$s)))
  g <- quiet_graft(graft_missing_clades, toy$phy, mtl, simmap = simmap)
  info <- extract_grafting_info(g, mtl, toy$phy)

  expect_equal(nrow(info), 2)
  expect_equal(attr(info, "type"), "clades")
  expect_equal(info$Taxa, c("M1, M2", "M3"))
  expect_equal(info$n_taxa, c(2, 1))
  # The parent clade is the sister clade
  expect_equal(info$parent_node, c(toy$node_B, match(toy$s, toy$phy$tip.label)))
  expect_true(info_same_clade(toy$phy, info$parent_tips[1], toy$node_B))
  expect_equal(info$parent_tips[2], toy$s)
  # Graft age and crown age
  expect_equal(info$graft_age, c(graft_age(g, c("M1", "M2")), graft_age(g, "M3")))
  expect_equal(info$crown_age[1], max(phytools::nodeHeights(g)[, 2]) - phytools::nodeheight(g, ape::getMRCA(g, c("M1", "M2"))))
  expect_true(is.na(info$crown_age[2]))
  expect_true(info$crown_age[1] < info$graft_age[1])
  # States
  states <- graft_states_in_simmap(g, c("M1", "M3"))
  expect_equal(info$clade_state[2], states$state[2])
  expect_equal(info$parent_state, states$parent_state)
  expect_false("tip_state" %in% names(info))
  # Several simulations
  sims <- quiet_graft(graft_missing_clades, toy$phy, mtl, nb_simulations = 3)
  expect_equal(nrow(extract_grafting_info(sims, mtl, toy$phy)), 6)

  # The taxa of a missing clade that do not form a clade: NA with a warning
  broken <- g
  broken <- ape::drop.tip(broken, "M2"); broken <- graft_tip_on_branch(broken, "M2", match("t1", broken$tip.label), frac = 0.5)
  broken$simmap <- NULL
  expect_warning(info_broken <- extract_grafting_info(broken, mtl, toy$phy), "not monophyletic")
  expect_true(is.na(info_broken$graft_age[1]))
  expect_false(is.na(info_broken$graft_age[2]))
})

test_that("extract_grafting_info() stops or warns on invalid inputs",
{
  toy <- toy_setup()
  mtl <- info_mtl(list(Taxon_name = "X1", Target_clade = list(toy$B)))
  g <- graft_tip_on_branch(toy$phy, "X1", toy$node_B)
  expect_error(extract_grafting_info(g, list(Taxon_name = "X1")), "'Target_clade' item")
  expect_error(extract_grafting_info(g, list(Taxon_name = "X1", Target_clade = list(X1 = toy$B), Sister_clade = list(X1 = toy$B))), "both")
  expect_error(extract_grafting_info(g, "a"), "must be a list")
  expect_error(extract_grafting_info("a", mtl), "phylogeny")
  expect_error(extract_grafting_info(g, mtl, phy_backbone = "a"), "phy_backbone")
  # A backbone with other tips
  other_backbone <- toy$phy; other_backbone$tip.label[1] <- "other"
  expect_error(extract_grafting_info(g, mtl, phy_backbone = other_backbone), "differ")
  # A taxon of the list that is not in the phylogeny: NA with a warning
  expect_warning(info <- extract_grafting_info(g, info_mtl(list(Taxon_name = c("X1", "Y"), Target_clade = list(toy$B, toy$B))), toy$phy), "not tips")
  expect_true(is.na(info$graft_age[2]))
  expect_false(is.na(info$graft_age[1]))
})
