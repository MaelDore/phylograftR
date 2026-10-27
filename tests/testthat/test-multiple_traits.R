
##### Test suite for the use of several traits (list of simmaps) to constrain grafting in graft_missing_taxa() and graft_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


## Toy lists of missing taxa / missing clades
traits_missing_taxa <- function (toy)
{
  list(Taxon_name = c("X1", "X2", "X3", "X4", "X5"),
       Target_clade = list(X1 = toy$A, X2 = toy$A, X3 = toy$B, X4 = toy$C, X5 = toy$s),
       Excluded_subclades = list(X1 = list(toy$B), X2 = NULL, X3 = NULL, X4 = NULL, X5 = NULL))
}
traits_missing_clades <- function (toy)
{
  list(Taxon_name = c("M1", "M2", "M3", "M4"),
       Sister_clade = list(M1 = toy$B, M2 = toy$B, M3 = toy$C, M4 = toy$s))
}
as_multi <- function (maps, names = NULL)
{
  names(maps) <- names
  structure(maps, class = c("multiSimmap", "multiPhylo"))
}

## Second trait: the toy states relabelled ("low", "mid", "high")
make_trait2_simmap <- function (phy, seed = 2) { relabel_simmap(make_toy_simmap(phy, seed = seed), c(A = "low", B = "mid", C = "high")) }

## Relabel the states of the maps only (simmaps without 'mapped.edge', e.g., from make_young_state_simmap())
relabel_maps <- function (simmap, mapping)
{
  simmap$maps <- lapply(simmap$maps, function (x) { names(x) <- unname(mapping[names(x)]); x })
  simmap
}

## Toy traits with states in all the lineages: "color" is "A" before the age 'age_h' and "B" after, "elevation" is "low" before the age 'age_e' and "high" after
# (the combination "B|high" is found on every terminal branch, and the combinations of the three periods on most lineages)
make_period_traits <- function (phy, age_h = 20, age_e = 10)
{
  list(color = make_young_state_simmap(phy, age_h), elevation = relabel_maps(make_young_state_simmap(phy, age_e), c(A = "low", B = "high")))
}

## State of each taxon in the simmap of a trait (the grafted branch is entirely in the state inherited)
taxa_states <- function (simmap, taxa)
{
  vapply(taxa, function (taxon) names(simmap$maps[[match(match(taxon, simmap$tip.label), simmap$edge[, 2])]])[1], character(1))
}

## Same elements except '$simmap' and the attribute 'map'
same_tree <- function (x, y)
{
  x$simmap <- NULL; y$simmap <- NULL
  attr(x, "map") <- NULL; attr(y, "map") <- NULL
  isTRUE(all.equal(x, y, tolerance = 0))
}


## ---- Parsing and checking of the inputs ---- ##

test_that("check_simmap_input() lists the maps of several traits, pairs them, and names the traits and the maps",
{
  toy <- toy_setup()
  one <- make_toy_simmap(toy$phy, seed = 1)
  two <- make_trait2_simmap(toy$phy)

  # Single maps: a single run, traits named trait_1, trait_2 by default
  out <- check_simmap_input(list(one, two))
  expect_equal(out$traits, c("trait_1", "trait_2"))
  expect_false(out$multiple)
  expect_null(out$names)
  expect_length(out$maps, 1)
  expect_s3_class(out$maps[[1]], "phylograftR_trait_maps")
  expect_equal(names(out$maps[[1]]), c("trait_1", "trait_2"))
  expect_identical(out$maps[[1]]$trait_2, two)
  expect_equal(check_simmap_input(list(color = one, elevation = two))$traits, c("color", "elevation"))

  # multiSimmap: the maps are paired (map k of every trait together), names map_1... by default
  m1 <- as_multi(list(one, one, one))
  m2 <- as_multi(list(two, make_trait2_simmap(toy$phy, seed = 3), make_trait2_simmap(toy$phy, seed = 4)))
  out <- check_simmap_input(list(a = m1, b = m2))
  expect_true(out$multiple)
  expect_equal(out$names, c("map_1", "map_2", "map_3"))
  expect_length(out$maps, 3)
  expect_identical(out$maps[[2]]$b, m2[[2]])
  expect_identical(out$maps[[3]]$a, m1[[3]])
  # A single simmap and a multiSimmap with a single map: one map, the output is still a list
  out <- check_simmap_input(list(one, as_multi(list(two))))
  expect_true(out$multiple)
  expect_equal(out$names, "map_1")

  # Names of the maps: the same for every multiSimmap, or none
  out <- check_simmap_input(list(a = as_multi(list(one, one), c("h1", "h2")), b = as_multi(list(two, two), c("h1", "h2"))))
  expect_equal(out$names, c("h1", "h2"))
  expect_error(check_simmap_input(list(a = as_multi(list(one, one), c("h1", "h2")), b = as_multi(list(two, two), c("h1", "h3")))), "named differently across the traits")
  expect_error(check_simmap_input(list(a = as_multi(list(one, one), c("h1", "h2")), b = as_multi(list(two, two)))), "named differently across the traits")
})

test_that("check_simmap_input() stops for the lists of maps that are not valid",
{
  toy <- toy_setup()
  one <- make_toy_simmap(toy$phy, seed = 1)
  two <- make_trait2_simmap(toy$phy)
  expect_error(check_simmap_input(list(one)), "list of 1 element.*several traits.*not a list")
  expect_error(check_simmap_input(list(a = one, one)), "named and others are not.*2")
  expect_error(check_simmap_input(list(a = one, a = two)), "different names.*Duplicated names: a")
  expect_error(check_simmap_input(list(`a;b` = one, c = two)), "cannot contain")
  expect_error(check_simmap_input(list(`a=b` = one, c = two)), "cannot contain")
  expect_error(check_simmap_input(list(one, "map")), "Each element.*'trait_2'")
  expect_error(check_simmap_input(list(one, toy$phy)), "Each element.*'trait_2'")
  # Same number of maps for every trait
  expect_error(check_simmap_input(list(a = as_multi(list(one, one, one)), b = as_multi(list(two, two)))), "same number of maps.*'a' = 3, 'b' = 2")
  expect_error(check_simmap_input(list(a = as_multi(list(one, one)), b = two)), "same number of maps.*'a' = 2, 'b' = 1")
  # Maps in a multiSimmap named partially or with duplicates
  expect_error(check_simmap_input(list(a = as_multi(list(one, one), c("x", "x")), b = two)), "trait 'a' must have different names.*x")
})

test_that("parse_trait_states_cell() splits the traits on ';' and the states on ',', and keeps empty fields",
{
  expect_equal(parse_trait_states_cell("A1,A2;B1"), list(c("A1", "A2"), "B1"))
  expect_equal(parse_trait_states_cell(" A1 , A2 ; B1 "), list(c("A1", "A2"), "B1"))
  expect_equal(parse_trait_states_cell("A1,A2;;C3"), list(c("A1", "A2"), NULL, "C3"))
  expect_equal(parse_trait_states_cell(";B1"), list(NULL, "B1"))
  expect_equal(parse_trait_states_cell("A1;"), list("A1", NULL))
  expect_equal(parse_trait_states_cell("A1;;"), list("A1", NULL, NULL))
  # Named traits
  expect_equal(parse_trait_states_cell("color=A1,A2; elevation=C3"), list(color = c("A1", "A2"), elevation = "C3"))
  expect_equal(parse_trait_states_cell("color=A1; elevation="), list(color = "A1", elevation = NULL))
  # Errors
  expect_error(parse_trait_states_cell("color=A1; C3"), "mixes named and unnamed traits")
  expect_error(parse_trait_states_cell("=A1; elevation=C3"), "without name")
  expect_error(parse_trait_states_cell("a=A1; a=C3"), "repeated traits.*a")
  expect_error(parse_trait_states_cell("A1,,A2;B1"), "empty state")
  expect_error(parse_trait_states_cell("A1,A1;B1"), "repeated states.*A1")
  expect_error(parse_trait_states_cell("A1,A2;B1", who = "taxon X1"), NA)
  expect_error(parse_trait_states_cell("A1,,A2;B1", who = "taxon X1"), "^taxon X1 has an empty state")
})

test_that("normalize_trait_states() gives one element per trait, from a list (positional or named), a text, or nothing",
{
  traits <- c("color", "elevation", "diet")
  empty <- stats::setNames(vector("list", 3), traits)
  expect_equal(normalize_trait_states(NULL, traits, "x"), empty)
  expect_equal(normalize_trait_states(NA, traits, "x"), empty)
  expect_equal(normalize_trait_states(list(), traits, "x"), empty)
  # Positional: NULL, empty vectors and NA have no states
  expect_equal(normalize_trait_states(list(c("A1", "A2"), NULL, "C3"), traits, "x"), list(color = c("A1", "A2"), elevation = NULL, diet = "C3"))
  expect_equal(normalize_trait_states(list(c("A1", "A1"), character(0), NA), traits, "x"), list(color = "A1", elevation = NULL, diet = NULL))
  # Named: the traits absent have no states, any order
  expect_equal(normalize_trait_states(list(diet = "C3", color = "A1"), traits, "x"), list(color = "A1", elevation = NULL, diet = "C3"))
  # Text
  expect_equal(normalize_trait_states("A1,A2;;C3", traits, "x"), list(color = c("A1", "A2"), elevation = NULL, diet = "C3"))
  expect_equal(normalize_trait_states("diet=C3; color=A1", traits, "x"), list(color = "A1", elevation = NULL, diet = "C3"))
  # Errors
  expect_error(normalize_trait_states(c("A1", "A2"), traits, "taxon X1"), "taxon X1 must be given for each trait.*color, elevation, diet")
  expect_error(normalize_trait_states("A1", traits, "x"), "must be given for each trait")
  expect_error(normalize_trait_states(list("A1", "B1"), traits, "taxon X1"), "taxon X1 have 2 element\\(s\\) but there are 3 traits")
  expect_error(normalize_trait_states(list(color = "A1", size = "B"), traits, "x"), "traits that are not in 'simmap': 'size'")
  expect_error(normalize_trait_states(list(color = "A1", "B"), traits, "x"), "name some traits but not all")
  expect_error(normalize_trait_states(list(color = "A1", color = "B"), traits, "x"), "repeat traits.*color")
  expect_error(normalize_trait_states(list("A1", 2, NULL), traits, "x"), "character vectors")
  expect_error(normalize_trait_states(1:3, traits, "x"), "must be given for each trait|must be a list")
})


## ---- Merge and split ---- ##

test_that("merge_trait_simmaps() cuts the branches at the change points of all the traits and names the pieces with the state of each trait",
{
  toy <- toy_setup()
  phy <- label_backbone_nodes(toy$phy)$phy
  one <- make_toy_simmap(toy$phy, seed = 1)
  two <- make_trait2_simmap(toy$phy, seed = 2)
  merged <- merge_trait_simmaps(phy = phy, trait_maps = list(color = one, elevation = two))
  simmap <- merged$simmap

  expect_s3_class(simmap, "simmap")
  expect_equal(merged$traits$names, c("color", "elevation"))
  expect_equal(merged$traits$states, list(color = c("A", "B", "C"), elevation = c("high", "low", "mid")))
  expect_equal(merged$traits$sep, "|")
  expect_true(all(grepl("^[ABC]\\|(high|low|mid)$", unlist(lapply(simmap$maps, names)))))
  # The branches keep their length and the pieces sum to it
  expect_equal(vapply(simmap$maps, sum, numeric(1)), simmap$edge.length, tolerance = 1e-9)
  expect_equal(unname(rowSums(simmap$mapped.edge)), simmap$edge.length, tolerance = 1e-9)

  # At any point of any branch, the state of each trait is the one of its own map: compare the total time in each state
  for (t in c("color", "elevation"))
  {
    original <- list(color = one, elevation = two)[[t]]
    split_maps <- split_trait_simmap(simmap, merged$traits)[[t]]
    expect_equal(split_maps$mapped.edge[, colnames(original$mapped.edge)], original$mapped.edge, ignore_attr = TRUE, tolerance = 1e-8)
  }

  # More pieces than in each map alone (a branch is cut at the change points of both traits), never more than the sum
  n_pieces <- lengths(simmap$maps)
  expect_true(all(n_pieces <= lengths(one$maps) + lengths(two$maps) - 1))
  expect_true(any(n_pieces > lengths(one$maps)))
})

test_that("merge_trait_simmaps() on a single branch: the pieces are the union of the change points, identical change points are merged",
{
  phy <- ape::read.tree(text = "((a:2,b:2):3,c:5);")
  phy <- phy_labelled <- label_backbone_nodes(phy)$phy
  make_simmap <- function (map_list)
  {
    simmap <- phy
    simmap$maps <- map_list
    simmap$mapped.edge <- NULL
    class(simmap) <- c("simmap", "phylo")
    simmap
  }
  # Edges (rows of phy$edge): find the branch to the tip 'a'
  edge_a <- match(match("a", phy$tip.label), phy$edge[, 2])
  base <- lapply(phy$edge.length, function (l) c(X = l))
  map1 <- base; map1[[edge_a]] <- c(X = 0.5, Y = 1.5)       # change at 0.5 from the rootward end
  map2 <- base; map2[[edge_a]] <- c(P = 1, Q = 1)           # change at 1
  merged <- merge_trait_simmaps(phy, list(t1 = make_simmap(map1), t2 = make_simmap(map2)))
  expect_equal(merged$simmap$maps[[edge_a]], c(`X|P` = 0.5, `Y|P` = 0.5, `Y|Q` = 1))
  # The other branches: one piece each, with the state of each trait
  expect_equal(unname(unlist(lapply(merged$simmap$maps[-edge_a], names))), rep("X|X", length(phy$edge.length) - 1))

  # Change points at the same place (up to the tolerance): no tiny piece
  map2[[edge_a]] <- c(P = 0.5 + 1e-10, Q = 1.5 - 1e-10)
  merged <- merge_trait_simmaps(phy, list(t1 = make_simmap(map1), t2 = make_simmap(map2)))
  expect_equal(merged$simmap$maps[[edge_a]], c(`X|P` = 0.5, `Y|Q` = 1.5), tolerance = 1e-8)

  # Adjacent pieces with the same composite state are merged: t2 changes where t1 does not matter... a change in a trait always makes a new composite state
  map3 <- base; map3[[edge_a]] <- c(X = 1, X = 1)           # two pieces in the same state in a trait
  merged <- merge_trait_simmaps(phy, list(t1 = make_simmap(map3), t2 = make_simmap(base)))
  expect_equal(merged$simmap$maps[[edge_a]], c(`X|X` = 2))
})

test_that("merge_trait_simmaps() stops for a trait whose simmap does not match the backbone or whose states contain the separator",
{
  toy <- toy_setup()
  phy <- label_backbone_nodes(toy$phy)$phy
  one <- make_toy_simmap(toy$phy, seed = 1)
  two <- make_trait2_simmap(toy$phy)
  other <- make_toy_simmap(make_toy_phy(n_tips = 20, seed = 99))
  expect_error(merge_trait_simmaps(phy, list(color = one, elevation = other)), "^Trait 'elevation': ")
  expect_error(merge_trait_simmaps(phy, list(color = other, elevation = two)), "^Trait 'color': ")
  bad <- relabel_simmap(two, c(low = "l|ow", mid = "mid", high = "high"))
  expect_error(merge_trait_simmaps(phy, list(color = one, elevation = bad)), "Trait 'elevation'.*cannot contain '\\|'.*'l\\|ow'")
})

test_that("split_trait_simmap() gives back the simmap of each trait, merging the adjacent pieces in the same state",
{
  toy <- toy_setup()
  phy <- label_backbone_nodes(toy$phy)$phy
  one <- make_toy_simmap(toy$phy, seed = 1)
  two <- make_trait2_simmap(toy$phy)
  merged <- merge_trait_simmaps(phy, list(color = one, elevation = two))
  split <- split_trait_simmap(merged$simmap, merged$traits)
  expect_equal(names(split), c("color", "elevation"))
  for (t in names(split))
  {
    expect_s3_class(split[[t]], "simmap")
    expect_equal(split[[t]]$edge, phy$edge)
    expect_equal(split[[t]]$edge.length, phy$edge.length)
    # No two adjacent pieces in the same state
    expect_false(any(vapply(split[[t]]$maps, function (m) any(names(m)[-1] == names(m)[-length(m)]), logical(1))))
    expect_equal(vapply(split[[t]]$maps, sum, numeric(1)), phy$edge.length, tolerance = 1e-9)
  }
  # The maps are those of the original simmaps (pieces of the same state merged)
  original_maps <- function (simmap) lapply(simmap$maps, function (m) { r <- rle(names(m)); stats::setNames(as.numeric(tapply(m, rep(seq_along(r$lengths), r$lengths), sum)), r$values) })
  expect_equal(split$color$maps, original_maps(one), tolerance = 1e-8)
  expect_equal(split$elevation$maps, original_maps(two), tolerance = 1e-8)
  # The states of the columns of mapped.edge are the states of the trait
  expect_equal(colnames(split$color$mapped.edge), c("A", "B", "C"))
  expect_equal(colnames(split$elevation$mapped.edge), c("high", "low", "mid"))
})


## ---- Resolution of the allowed states ---- ##

test_that("resolve_trait_allowed_states() gives the composite states with an allowed state for every trait",
{
  traits <- list(names = c("color", "elevation"), states = list(color = c("A", "B"), elevation = c("high", "low")), sep = "|")
  composite <- c("A|high", "A|low", "B|high")   # The combination B|low is not found on the maps
  resolve <- function (...) resolve_trait_allowed_states(traits = traits, composite_states = composite, taxon_names = c("X1", "X2", "X3"), ...)

  # Without any states: everything is allowed
  expect_equal(resolve(), list(X1 = composite, X2 = composite, X3 = composite))
  # Global: a NULL element does not restrict the trait
  expect_equal(resolve(allowed_states = list("A", NULL)), list(X1 = c("A|high", "A|low"), X2 = c("A|high", "A|low"), X3 = c("A|high", "A|low")))
  expect_equal(resolve(allowed_states = list(NULL, "high")), list(X1 = c("A|high", "B|high"), X2 = c("A|high", "B|high"), X3 = c("A|high", "B|high")))
  expect_equal(resolve(allowed_states = list(color = "B")), list(X1 = "B|high", X2 = "B|high", X3 = "B|high"))
  expect_equal(resolve(allowed_states = list("A", "low"))$X1, "A|low")
  # Taxon-level: replaces the global value of the trait it gives, the trait with NULL (or absent) uses the global value
  taxa <- list(X1 = list("B", NULL), X2 = list(NULL, "low"), X3 = NULL)
  expect_equal(resolve(allowed_states_taxa = taxa, allowed_states = list("A", "high")),
               list(X1 = "B|high", X2 = "A|low", X3 = "A|high"))
  # A text, and named traits
  taxa <- list(X1 = "B;", X2 = "elevation=low", X3 = ";low")
  expect_equal(resolve(allowed_states_taxa = taxa, allowed_states = list("A", NULL)),
               list(X1 = "B|high", X2 = c("A|low"), X3 = "A|low"))

  # Errors: unknown state (naming the trait), no combination found, wrong structure
  expect_error(resolve(allowed_states = list("C", NULL)), "global 'allowed_states' argument \\(trait 'color'\\).*'C'.*'A', 'B'")
  expect_error(resolve(allowed_states_taxa = list(X1 = NULL, X2 = list(NULL, "mid"), X3 = NULL)), "taxon 'X2'.*\\(trait 'elevation'\\).*'mid'")
  expect_error(resolve(allowed_states = list("B", "low")), "No combination of the allowed states of taxon 'X1'.*color = B; elevation = low.*A\\|high")
  expect_error(resolve(allowed_states = c("A", "high")), "must be given for each trait")
  expect_error(resolve(allowed_states_taxa = list(X1 = c("A", "high"), X2 = NULL, X3 = NULL)), "taxon 'X1'.*must be given for each trait")
  expect_error(resolve(allowed_states_taxa = list(X1 = NULL)), "same length as 'Taxon_name'")
  expect_error(resolve(allowed_states_taxa = c("A", "B", "C")), "must be a list")
})


## ---- Grafting with several traits ---- ##

for (name in c("taxa", "clades"))
{
  graft_traits <- function (toy, ..., missing = NULL)
  {
    if (name == "taxa") graft_missing_taxa(toy$phy, if (is.null(missing)) traits_missing_taxa(toy) else missing, method = "uniform", verbose = FALSE, ...)
    else graft_missing_clades(toy$phy, if (is.null(missing)) traits_missing_clades(toy) else missing, verbose = FALSE, ...)
  }
  names_all <- if (name == "taxa") c("X1", "X2", "X3", "X4", "X5") else c("M1", "M2", "M3", "M4")

  test_that(paste0("graft_missing_", name, "(): several traits, the grafting is restricted to the points with an allowed state in every trait, and $simmap is a list of simmaps named after the traits"),
  {
    toy <- toy_setup()
    color <- make_toy_simmap(toy$phy, seed = 1)
    elevation <- make_trait2_simmap(toy$phy, seed = 2)
    taxa <- names_all
    for (seed in 1:4)
    {
      out <- graft_traits(toy, simmap = list(color = color, elevation = elevation), allowed_states = list(c("A", "B"), c("low", "mid")), seed = seed)
      expect_s3_class(out, "phylo")
      expect_type(out$simmap, "list")
      expect_equal(names(out$simmap), c("color", "elevation"))
      expect_s3_class(out$simmap$color, "simmap"); expect_s3_class(out$simmap$elevation, "simmap")
      expect_true(all(taxa_states(out$simmap$color, taxa) %in% c("A", "B")))
      expect_true(all(taxa_states(out$simmap$elevation, taxa) %in% c("low", "mid")))
      # Each simmap is a map of the grafted tree
      for (t in names(out$simmap))
      {
        expect_equal(out$simmap[[t]]$edge, out$edge)
        expect_equal(out$simmap[[t]]$edge.length, out$edge.length)
        expect_equal(vapply(out$simmap[[t]]$maps, sum, numeric(1)), out$edge.length, tolerance = 1e-8)
        expect_true(all(names(unlist(out$simmap[[t]]$maps)) != ""))
      }
      expect_null(attr(out, "map")) # A single run: no map name
    }
  })

  test_that(paste0("graft_missing_", name, "(): several traits, the backbone keeps the states of each trait and a taxon inherits the state of each trait at the graft point"),
  {
    toy <- toy_setup()
    color <- make_toy_simmap(toy$phy, seed = 1)
    elevation <- make_trait2_simmap(toy$phy, seed = 2)
    out <- graft_traits(toy, simmap = list(color = color, elevation = elevation), allowed_states = list(c("A", "B", "C"), NULL), seed = 3)
    backbone <- toy$phy$tip.label
    # The states of the backbone tips are unchanged
    expect_equal(taxa_states_tip_end <- vapply(backbone, function (tip) { m <- out$simmap$color$maps[[match(match(tip, out$tip.label), out$edge[, 2])]]; names(m)[length(m)] }, character(1)),
                 vapply(backbone, function (tip) { m <- color$maps[[match(match(tip, color$tip.label), color$edge[, 2])]]; names(m)[length(m)] }, character(1)))
    # The grafted taxa inherit the state of the lineage at the graft point, for each trait (independent computation)
    for (taxon in names_all)
    {
      if (name == "taxa") { in_missing <- traits_missing_taxa(toy) } else { in_missing <- traits_missing_clades(toy) }
      state_h <- state_at_graft(out, color, taxon, in_missing$Taxon_name)
      state_e <- state_at_graft(out, elevation, taxon, in_missing$Taxon_name)
      if (!is.na(state_h)) expect_equal(unname(taxa_states(out$simmap$color, taxon)), state_h)
      if (!is.na(state_e)) expect_equal(unname(taxa_states(out$simmap$elevation, taxon)), state_e)
    }
  })

  test_that(paste0("graft_missing_", name, "(): a trait without restriction gives the same grafting as the other trait alone, and the same simmap for this trait"),
  {
    toy <- toy_setup()
    color <- make_toy_simmap(toy$phy, seed = 1)
    constant <- relabel_simmap(make_toy_simmap(toy$phy, seed = 1), c(A = "z", B = "z", C = "z")) # A trait with a single state everywhere
    constant$maps <- lapply(constant$maps, function (m) stats::setNames(sum(m), "z"))
    constant$mapped.edge <- matrix(constant$edge.length, ncol = 1, dimnames = list(rownames(constant$mapped.edge), "z"))
    for (seed in 1:3)
    {
      single <- graft_traits(toy, simmap = color, allowed_states = c("A", "B"), seed = seed)
      several <- graft_traits(toy, simmap = list(color = color, other = constant), allowed_states = list(c("A", "B"), NULL), seed = seed)
      expect_true(same_tree(single, several))
      expect_equal(several$simmap$color$maps, single$simmap$maps, tolerance = 1e-8)
      expect_equal(unname(unlist(lapply(several$simmap$other$maps, names))), rep("z", nrow(out_edges <- several$edge)))
    }
  })

  test_that(paste0("graft_missing_", name, "(): several traits, taxon-level allowed states (list, named list or text) replace the global value trait by trait"),
  {
    toy <- toy_setup()
    maps <- make_period_traits(toy$phy)
    missing <- if (name == "taxa") traits_missing_taxa(toy) else traits_missing_clades(toy)
    n <- length(missing$Taxon_name)
    spec <- stats::setNames(vector("list", n), missing$Taxon_name)
    # Taxa grafted together as a missing clade share their states: the first two taxa have the same value
    spec[[1]] <- list("B", NULL)               # color B, elevation: global
    spec[[2]] <- if (name == "taxa") list(elevation = "high") else list("B", NULL)  # elevation high (taxa), same as the first (clades)
    spec[[3]] <- "B;high"                      # color B and elevation high
    missing$allowed_states <- spec
    for (seed in 1:3)
    {
      out <- graft_traits(toy, simmap = maps, allowed_states = list(c("A", "B"), c("low", "high")), seed = seed, missing = missing)
      taxa <- missing$Taxon_name
      expect_equal(unname(taxa_states(out$simmap$color, taxa[1])), "B")
      if (name == "taxa") expect_equal(unname(taxa_states(out$simmap$elevation, taxa[2])), "high")
      expect_equal(unname(taxa_states(out$simmap$color, taxa[3])), "B")
      expect_equal(unname(taxa_states(out$simmap$elevation, taxa[3])), "high")
    }
  })

  test_that(paste0("graft_missing_", name, "(): paired maps of several traits give one grafted phylogeny per pair, as when each pair is given alone"),
  {
    toy <- toy_setup()
    pairs <- lapply(list(c(20, 10), c(15, 8), c(25, 12)), function (ages) make_period_traits(toy$phy, ages[1], ages[2]))
    colors <- lapply(pairs, function (x) x$color)
    elevations <- lapply(pairs, function (x) x$elevation)
    allowed <- if (name == "taxa") list(c("A", "B"), "high") else list(c("A", "B"), c("low", "high")) # (the stem branches of sister clades can be old)
    out <- graft_traits(toy, simmap = list(color = as_multi(colors), elevation = as_multi(elevations)), allowed_states = allowed, seed = 11)
    # One element per pair (3 and not 9), named map_1, map_2, map_3
    expect_equal(names(out), c("map_1", "map_2", "map_3"))
    expect_equal(nrow(attr(out, "failed_maps")), 0)
    for (k in 1:3)
    {
      single <- graft_traits(toy, simmap = list(color = colors[[k]], elevation = elevations[[k]]), allowed_states = allowed, seed = 11)
      expect_true(same_tree(out[[k]], single))
      expect_equal(out[[k]]$simmap$color$maps, single$simmap$color$maps)
      expect_equal(out[[k]]$simmap$elevation$maps, single$simmap$elevation$maps)
      expect_equal(attr(out[[k]], "map"), paste0("map_", k))
      if (name == "taxa") expect_true(all(taxa_states(out[[k]]$simmap$elevation, names_all) == "high"))
    }
    # With several simulations: a multiPhylo per pair, named <map>_sim_k
    several <- graft_traits(toy, simmap = list(color = as_multi(colors, c("a", "b", "c")), elevation = as_multi(elevations, c("a", "b", "c"))),
                            allowed_states = allowed, nb_simulations = 2, seed = 11)
    expect_equal(names(several), c("a", "b", "c"))
    expect_s3_class(several$b, "multiPhylo")
    expect_equal(names(several$b), c("b_sim_1", "b_sim_2"))
    expect_equal(names(several$b[[2]]$simmap), c("color", "elevation"))
    expect_equal(attr(several$c[[1]], "map"), "c")
    # A single simmap for each trait: one phylogeny, not a list
    expect_s3_class(graft_traits(toy, simmap = list(colors[[1]], elevations[[1]]), seed = 1), "phylo")
    # A multiSimmap with a single map: a list of one element
    one_pair <- graft_traits(toy, simmap = list(as_multi(colors[1]), elevations[[1]]), seed = 1)
    expect_equal(names(one_pair), "map_1")
  })

  test_that(paste0("graft_missing_", name, "(): several traits stop before any simulation for the inputs that are not valid"),
  {
    toy <- toy_setup()
    maps <- make_period_traits(toy$phy)
    color <- make_toy_simmap(toy$phy, seed = 1)
    elevation <- make_trait2_simmap(toy$phy, seed = 2)
    calls <- 0
    original <- run_simulations
    local_mocked_bindings(run_simulations = function (...) { calls <<- calls + 1; original(...) })

    # Ranges with several traits: not yet supported by graft_missing_clades(), a single TRUE is ambiguous in graft_missing_taxa()
    if (name == "taxa") { expect_error(graft_traits(toy, simmap = maps, treat_states_as_ranges = TRUE), "ambiguous with the maps of several traits") }
    else { expect_error(graft_traits(toy, simmap = maps, treat_states_as_ranges = TRUE), "not supported with the maps of several traits") }
    # Allowed states of a single trait, or of the wrong structure
    expect_error(graft_traits(toy, simmap = maps, allowed_states = c("A", "low")), "must be given for each trait")
    expect_error(graft_traits(toy, simmap = maps, allowed_states = list("A")), "1 element\\(s\\) but there are 2 traits")
    expect_error(graft_traits(toy, simmap = maps, allowed_states = list(c("A", "low"), NULL)), "\\(trait 'color'\\).*'low'")
    expect_error(graft_traits(toy, simmap = maps, allowed_states = list("A", "C")), "\\(trait 'elevation'\\).*'C'")
    expect_error(graft_traits(toy, simmap = maps, allowed_states = list(size = "A")), "traits that are not in 'simmap': 'size'")
    # No combination of states on the maps: the states exist in each trait, but never together
    # (color "A" before the age 20, elevation "high" after the age 10)
    expect_error(graft_traits(toy, simmap = maps, allowed_states = list("A", "high")), "No combination of the allowed states.*color = A; elevation = high")
    # A map that does not match the backbone, naming the trait
    other <- make_toy_simmap(make_toy_phy(n_tips = 20, seed = 99))
    expect_error(graft_traits(toy, simmap = list(color = color, elevation = other)), "Trait 'elevation'")
    # A map of a pair that does not match: the message names the map and the trait
    expect_error(graft_traits(toy, simmap = list(color = as_multi(list(color, color)), elevation = as_multi(list(elevation, other)))), "Map 'map_2': Trait 'elevation'")
    expect_equal(calls, 0)
    # Unequal numbers of maps
    expect_error(graft_traits(toy, simmap = list(color = as_multi(list(color, color)), elevation = as_multi(list(elevation)))), "same number of maps")
    # A list of states with the map of a single trait
    expect_error(graft_traits(toy, simmap = color, allowed_states = list("A", "B")), "list \\(one element per trait\\), but 'simmap' is the map of a single trait")
    expect_equal(calls, 0)
  })

  test_that(paste0("graft_missing_", name, "(): several traits with 'return_simmap = FALSE' do not return any simmap, and a list of maps without any allowed states warns"),
  {
    toy <- toy_setup()
    maps <- make_period_traits(toy$phy)
    out <- graft_traits(toy, simmap = maps, allowed_states = list("B", NULL), return_simmap = FALSE, seed = 1)
    expect_null(out$simmap)
    expect_warning(graft_traits(toy, simmap = maps, return_simmap = FALSE, seed = 1), "has no effect")
    # Allowed states that are all NULL are not restrictions either
    expect_warning(graft_traits(toy, simmap = maps, allowed_states = list(NULL, NULL), return_simmap = FALSE, seed = 1), "has no effect")
    # All the states of every trait: the same tree as without maps, and the simmaps are returned
    unrestricted <- graft_traits(toy, simmap = maps, seed = 5)
    plain <- graft_traits(toy, seed = 5)
    expect_true(same_tree(unrestricted, plain))
    expect_equal(names(unrestricted$simmap), c("color", "elevation"))
  })
}


## ---- Other functions ---- ##

test_that("a list of maps of several traits is usable with identify_grafted() and flatten_grafted_phylogenies(), and check_grafting_constraints()",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  maps <- make_period_traits(toy$phy)
  out <- graft_missing_taxa(toy$phy, mtl, simmap = maps, allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1)
  # The grafted tips and branches are identified
  expect_equal(sum(out$grafted_tips), length(mtl$Taxon_name))
  expect_equal(out$tip.label[out$grafted_tips], out$tip.label[out$tip.label %in% mtl$Taxon_name])
  # check_grafting_constraints(): works with the simmaps (see the dedicated tests below), and without them
  checked <- check_grafting_constraints(out, toy$phy, mtl, allowed_states = list("B", NULL), verbose = FALSE)
  expect_true(all(checked$ok))
  checked <- check_grafting_constraints(drop_simmap(out), toy$phy, mtl, verbose = FALSE)
  expect_true(all(checked$ok))
  expect_true(all(is.na(checked$states_ok)))
  # Flatten: the simmaps can be dropped
  flat <- flatten_grafted_phylogenies(out, keep_simmap = FALSE)
  expect_null(flat$simmap)
  expect_s3_class(flatten_grafted_phylogenies(out), "phylo")
  # With paired multiSimmaps
  pairs <- lapply(list(c(20, 10), c(15, 8)), function (ages) make_period_traits(toy$phy, ages[1], ages[2]))
  multi <- graft_missing_taxa(toy$phy, mtl, simmap = list(color = as_multi(lapply(pairs, function (x) x$color)), elevation = as_multi(lapply(pairs, function (x) x$elevation))),
                              allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1, nb_simulations = 2)
  flat <- flatten_grafted_phylogenies(multi, keep_simmap = FALSE)
  expect_s3_class(flat, "multiPhylo")
  expect_length(flat, 4)
  expect_true(all(vapply(flat, function (x) is.null(x$simmap), logical(1))))
})

test_that("import_missing_taxa_list() reads the allowed states of several traits (positional or named) separated by ';'",
{
  toy <- toy_setup()
  write_table <- function (allowed_states)
  {
    table <- data.frame(Taxon_name = c("X1", "X2", "X3", "X4"),
                        Target_clade = vapply(c("A", "A", "B", "C"), function (x) paste(toy[[x]], collapse = ", "), character(1)),
                        allowed_states = allowed_states, stringsAsFactors = FALSE)
    path <- tempfile(fileext = ".csv")
    utils::write.csv(table, path, row.names = FALSE, na = "")
    path
  }
  import <- function (allowed_states) suppressMessages(import_missing_taxa_list(write_table(allowed_states), phy = toy$phy, verbose = FALSE))
  imported <- import(c("A, B;;high", "color=B; elevation=low, high", NA, "C"))
  expect_equal(imported$allowed_states$X1, list(c("A", "B"), NULL, "high"))
  expect_equal(imported$allowed_states$X2, list(color = "B", elevation = c("low", "high")))
  expect_null(imported$allowed_states$X3)
  expect_equal(imported$allowed_states$X4, "C") # One trait: unchanged
  expect_error(import(c("A,,B;low", NA, NA, NA)), "empty state")
  expect_error(import(c("color=A; low", NA, NA, NA)), "mixes named and unnamed")

  # Used with the maps of two traits
  maps <- make_period_traits(toy$phy)
  imported <- import(c("A, B;high", "color=B; elevation=low, high", NA, "B;"))
  imported$Excluded_subclades <- NULL
  out <- graft_missing_taxa(toy$phy, imported, simmap = maps, method = "uniform", verbose = FALSE, seed = 1)
  expect_true(unname(taxa_states(out$simmap$color, "X1")) %in% c("A", "B"))
  expect_equal(unname(taxa_states(out$simmap$elevation, "X1")), "high")
  expect_equal(unname(taxa_states(out$simmap$color, "X2")), "B")
  expect_equal(unname(taxa_states(out$simmap$color, "X4")), "B")
  # A taxon with 3 fields with two traits
  imported$allowed_states$X1 <- list("A", NULL, "high")
  expect_error(graft_missing_taxa(toy$phy, imported, simmap = maps, method = "uniform", verbose = FALSE), "taxon 'X1'.*3 element\\(s\\) but there are 2 traits")
  # With one trait (a character vector of states): the usual error with a list
  imported$allowed_states$X1 <- list("A", "high")
  expect_error(graft_missing_taxa(toy$phy, imported, simmap = maps$color, method = "uniform", verbose = FALSE), "given as a list \\(one element per trait\\)")
})


## ---- check_grafting_constraints() with several traits ---- ##

test_that("check_grafting_constraints() checks the states of each trait and of all the traits together",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  maps <- make_period_traits(toy$phy)
  grafted <- graft_missing_taxa(toy$phy, mtl, simmap = maps, allowed_states = list("B", "high"), method = "uniform", verbose = FALSE, seed = 1)
  taxa <- mtl$Taxon_name

  # Same constraints as for the grafting: all respected, with one column per trait after the joint 'states_ok'
  res <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("B", "high"), verbose = FALSE)
  expect_true(all(res$ok))
  expect_true(attr(res, "all_ok"))
  expect_true(all(res$states_ok))
  expect_true(all(res$states_ok_color))
  expect_true(all(res$states_ok_elevation))
  cols <- names(res)
  expect_equal(cols[match("states_ok", cols) + 0:2], c("states_ok", "states_ok_color", "states_ok_elevation"))
  expect_equal(attr(res, "trait_names"), c("color", "elevation"))

  # Positional or named allowed states are equivalent; NULL = all the states of the trait
  expect_equal(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list(elevation = "high", color = "B"), verbose = FALSE), res)
  expect_true(all(check_grafting_constraints(grafted, toy$phy, mtl, verbose = FALSE)$states_ok))
  expect_true(all(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list(NULL, "high"), verbose = FALSE)$ok))

  # A violation on one trait only: the joint column and the trait column are FALSE, the other trait is TRUE
  res_bad <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("A", "high"), verbose = FALSE)
  expect_false(any(res_bad$states_ok))
  expect_false(any(res_bad$states_ok_color))
  expect_true(all(res_bad$states_ok_elevation))
  expect_false(any(res_bad$ok))
  expect_false(attr(res_bad, "all_ok"))

  # The violating trait is named in the summary (counts per trait) and in the description of the failed taxa
  txt <- paste(capture.output(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("A", "high"), verbose = TRUE)), collapse = "\n")
  expect_match(txt, "Allowed states violated:\\s+5 \\(color: 5, elevation: 0\\)")
  expect_match(txt, "color: B \\(allowed: A\\)")
  expect_false(grepl("elevation: high", txt, fixed = TRUE))
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("A", "high"), verbose = FALSE, stop_on_failure = TRUE), "color: B \\(allowed: A\\)")
  txt_ok <- paste(capture.output(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("B", "high"), verbose = TRUE)), collapse = "\n")
  expect_match(txt_ok, "Allowed states violated:\\s+0 \\(color: 0, elevation: 0\\)")

  # Violations on both traits
  res_both <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("A", "low"), verbose = FALSE)
  expect_false(any(res_both$states_ok_color)); expect_false(any(res_both$states_ok_elevation))
})

test_that("check_grafting_constraints() resolves the taxon-level allowed states of each trait as the grafting functions do",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  maps <- make_period_traits(toy$phy)
  taxa <- mtl$Taxon_name
  # X1: color A or B, elevation high; X2: all colors, elevation low or high; the others use the global value (color B)
  allowed <- stats::setNames(vector("list", length(taxa)), taxa)
  allowed$X1 <- list(c("A", "B"), "high")
  allowed$X2 <- list(NULL, c("low", "high"))
  mtl$allowed_states <- allowed
  grafted <- graft_missing_taxa(toy$phy, mtl, simmap = maps, allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1)
  res <- check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("B", NULL), verbose = FALSE)
  expect_true(all(res$ok))
  # Same with the text syntax of the tables
  mtl_text <- mtl
  mtl_text$allowed_states <- stats::setNames(list("A,B;high", ";low,high", NA, NA, NA), taxa)
  expect_equal(check_grafting_constraints(grafted, toy$phy, mtl_text, allowed_states = list("B", NULL), verbose = FALSE), res)
  # The taxon-level value replaces the global one (it is not intersected): X1 can be in "A" even if the global value is "B"
  expect_equal(taxa_states(grafted$simmap$color, "X1") %in% c("A", "B"), TRUE)
  # A taxon-level value that the grafted state does not respect is a violation of this taxon and trait only
  mtl_strict <- mtl
  mtl_strict$allowed_states$X3 <- list("A", NULL)
  res_strict <- check_grafting_constraints(grafted, toy$phy, mtl_strict, allowed_states = list("B", NULL), verbose = FALSE)
  expect_equal(res_strict$Taxon_name[!res_strict$ok], "X3")
  expect_false(res_strict$states_ok_color[res_strict$Taxon_name == "X3"])
  expect_true(all(res_strict$states_ok_elevation))
})

test_that("check_grafting_constraints() uses the names of the traits, or trait_1, trait_2... if the traits are not named, and detects the states modified in a map",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  maps <- make_period_traits(toy$phy)
  unnamed <- graft_missing_taxa(toy$phy, mtl, simmap = unname(maps), allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1)
  res <- check_grafting_constraints(unnamed, toy$phy, mtl, allowed_states = list("B", NULL), verbose = FALSE)
  expect_true(all(c("states_ok_trait_1", "states_ok_trait_2") %in% names(res)))
  expect_true(all(res$ok))
  # Simmap without names in the grafted phylogeny (names removed afterwards): same default
  grafted <- graft_missing_taxa(toy$phy, mtl, simmap = maps, allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1)
  renamed <- grafted; names(renamed$simmap) <- NULL
  expect_true(all(c("states_ok_trait_1", "states_ok_trait_2") %in% names(check_grafting_constraints(renamed, toy$phy, mtl, verbose = FALSE))))
  # A modified state on the terminal branch of a taxon is detected for its trait only
  bad <- grafted
  tip_edge <- match(match("X2", bad$tip.label), bad$edge[, 2])
  names(bad$simmap$elevation$maps[[tip_edge]]) <- "low"
  res_bad <- check_grafting_constraints(bad, toy$phy, mtl, allowed_states = list("B", "high"), verbose = FALSE)
  expect_equal(res_bad$Taxon_name[!res_bad$ok], "X2")
  expect_true(res_bad$states_ok_color[res_bad$Taxon_name == "X2"])
  expect_false(res_bad$states_ok_elevation[res_bad$Taxon_name == "X2"])
  # A simmap with the wrong number of branches
  bad$simmap$color$maps <- bad$simmap$color$maps[-1]
  expect_error(check_grafting_constraints(bad, toy$phy, mtl, verbose = FALSE), "simmap of trait 'color'.*one map per branch")
})

test_that("check_grafting_constraints() stops with several traits and ranges, or with allowed states that are not given per trait or not in the maps",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  maps <- make_period_traits(toy$phy)
  grafted <- graft_missing_taxa(toy$phy, mtl, simmap = maps, allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1)
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, treat_states_as_ranges = TRUE, verbose = FALSE), "not supported with the maps of several traits")
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = "B", verbose = FALSE), "for each trait")
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list("B", "Z"), verbose = FALSE), "Z")
  expect_error(check_grafting_constraints(grafted, toy$phy, mtl, allowed_states = list(habitat = "B"), verbose = FALSE), "not in 'simmap'")
  expect_error(check_grafting_constraints(drop_simmap(grafted), toy$phy, mtl, allowed_states = list("B", NULL), verbose = FALSE), "can only be used if")
})

test_that("check_grafting_constraints() checks several traits with several simulations and with paired multiSimmaps, and one trait is unchanged",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  maps <- make_period_traits(toy$phy)
  # Several simulations of one run
  sims <- graft_missing_taxa(toy$phy, mtl, simmap = maps, allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1, nb_simulations = 3)
  res <- check_grafting_constraints(sims, toy$phy, mtl, allowed_states = list("B", NULL), verbose = FALSE)
  expect_equal(nrow(res), 3 * length(mtl$Taxon_name))
  expect_true(all(res$ok))
  expect_true(all(c("states_ok_color", "states_ok_elevation") %in% names(res)))
  txt <- paste(capture.output(check_grafting_constraints(sims, toy$phy, mtl, allowed_states = list("A", NULL), verbose = TRUE)), collapse = "\n")
  expect_match(txt, "\\(color: 15, elevation: 0\\)")
  # Paired multiSimmaps: one data.frame per map
  pairs <- lapply(list(c(20, 10), c(15, 8)), function (ages) make_period_traits(toy$phy, ages[1], ages[2]))
  multi <- graft_missing_taxa(toy$phy, mtl, simmap = list(color = as_multi(lapply(pairs, function (x) x$color)), elevation = as_multi(lapply(pairs, function (x) x$elevation))),
                              allowed_states = list("B", NULL), method = "uniform", verbose = FALSE, seed = 1, nb_simulations = 2)
  res_multi <- check_grafting_constraints(multi, toy$phy, mtl, allowed_states = list("B", NULL), verbose = FALSE)
  expect_type(res_multi, "list")
  expect_length(res_multi, 2)
  expect_true(all(vapply(res_multi, function (r) all(r$ok) && "states_ok_elevation" %in% names(r), logical(1))))
  # A single trait: no per-trait columns
  one <- graft_missing_taxa(toy$phy, mtl, simmap = maps$color, allowed_states = "B", method = "uniform", verbose = FALSE, seed = 1)
  res_one <- check_grafting_constraints(one, toy$phy, mtl, allowed_states = "B", verbose = FALSE)
  expect_false(any(grepl("^states_ok_", names(res_one))))
  expect_null(attr(res_one, "trait_names"))
  expect_false(grepl("color:", paste(capture.output(check_grafting_constraints(one, toy$phy, mtl, allowed_states = "B", verbose = TRUE)), collapse = "\n"), fixed = TRUE))
})


## ---- Ranges: one trait treated as biogeographic ranges, with other traits (graft_missing_taxa) ---- ##

## Toy traits: ranges "A" (before the age 'age_r') then "AB", and an elevation "low" (before 'age_e') then "high" (state combinations: A|low, AB|low, AB|high)
make_range_traits <- function (phy, age_r = 20, age_e = 10)
{
  list(biogeo = relabel_maps(make_young_state_simmap(phy, age_r), c(A = "A", B = "AB")), elevation = relabel_maps(make_young_state_simmap(phy, age_e), c(A = "low", B = "high")))
}
graft_ranges <- function (toy, maps, ..., mtl = NULL)
{
  suppressWarnings(graft_missing_taxa(toy$phy, if (is.null(mtl)) traits_missing_taxa(toy) else mtl, simmap = maps, method = "uniform", verbose = FALSE, ...))
}
state_of <- function (grafted, trait, taxa) taxa_states(grafted$simmap[[trait]], taxa)
parent_state_of <- function (simmap, taxa)
{
  vapply(taxa, function (taxon)
  {
    tip_edge <- match(match(taxon, simmap$tip.label), simmap$edge[, 2])
    parent_edge <- match(simmap$edge[tip_edge, 1], simmap$edge[, 2])
    x <- names(simmap$maps[[parent_edge]]); x[length(x)]
  }, character(1))
}

test_that("resolve_range_trait() validates the logical vector of treat_states_as_ranges: positional or named, one TRUE at most",
{
  traits <- c("biogeo", "color", "elevation")
  # No range trait
  expect_null(resolve_range_trait(FALSE, traits))
  expect_null(resolve_range_trait(c(FALSE, FALSE, FALSE), traits))
  expect_null(resolve_range_trait(c(biogeo = FALSE, color = FALSE, elevation = FALSE), traits))
  # Positional, named in any order
  expect_equal(resolve_range_trait(c(FALSE, TRUE, FALSE), traits), 2L)
  expect_equal(resolve_range_trait(c(color = FALSE, elevation = FALSE, biogeo = TRUE), traits), 1L)
  expect_equal(resolve_range_trait(c(FALSE, FALSE, TRUE), traits, traits_named = FALSE), 3L)
  # A single TRUE is ambiguous and lists the traits
  expect_error(resolve_range_trait(TRUE, traits), "ambiguous.*'biogeo', 'color', 'elevation'.*c\\(biogeo = TRUE, color = FALSE, elevation = FALSE\\)")
  expect_error(resolve_range_trait(TRUE, c("trait_1", "trait_2"), traits_named = FALSE), "c\\(TRUE, FALSE\\)")
  # Invalid vectors
  expect_error(resolve_range_trait(c(TRUE, FALSE), traits), "one value per trait \\(3 traits.*but has 2")
  expect_error(resolve_range_trait(c(TRUE, NA, FALSE), traits), "without NA")
  expect_error(resolve_range_trait("A", traits), "must be TRUE or FALSE")
  expect_error(resolve_range_trait(c(TRUE, TRUE, FALSE), traits), "Only one trait.*'biogeo', 'color'")
  expect_error(resolve_range_trait(c(biogeo = TRUE, color = FALSE, habitat = FALSE), traits), "exactly the names of the traits")
  expect_error(resolve_range_trait(c(biogeo = TRUE, color = FALSE, color = FALSE), traits), "exactly the names of the traits")
  expect_error(resolve_range_trait(stats::setNames(c(TRUE, FALSE, FALSE), c("biogeo", "", "elevation")), traits), "named and others are not")
  expect_error(resolve_range_trait(c(trait_1 = TRUE, trait_2 = FALSE), c("trait_1", "trait_2"), traits_named = FALSE), "named but the traits of 'simmap' are not")
})

test_that("check_simmap_input() tells whether the traits were named",
{
  toy <- toy_setup()
  maps <- make_range_traits(toy$phy)
  expect_true(check_simmap_input(maps)$traits_named)
  expect_false(check_simmap_input(unname(maps))$traits_named)
})

test_that("resolve_trait_allowed_ranges() lifts the ranges of the range trait to the composite states, keeping the states of the other traits",
{
  traits <- list(names = c("biogeo", "elevation"), states = list(biogeo = c("A", "AB"), elevation = c("low", "high")), sep = "|")
  composite <- c("A|low", "AB|low", "AB|high")
  taxa <- c("X1", "X2", "X3")
  # X1: range A and elevation high, X2: range A only (any elevation), X3: nothing (no restriction)
  allowed <- list(X1 = list("A", "high"), X2 = list("A", NULL), X3 = NULL)
  res <- resolve_trait_allowed_ranges(allowed_states_taxa = allowed, allowed_states = NULL, traits = traits, composite_states = composite, taxon_names = taxa, range_trait = 1)
  expect_equal(res$inherited_states, c("A|high"))
  expect_equal(res$allowed_states$X1, c("AB|high", "A|high"))
  expect_equal(res$range_candidates$X1, list(`AB|high` = "A|high", `A|high` = "A|high"))
  expect_equal(res$allowed_states$X2, c("A|low", "AB|low", "AB|high", "A|high"))
  expect_equal(res$range_candidates$X2, list(`A|low` = "A|low", `AB|low` = "A|low", `AB|high` = "A|high", `A|high` = "A|high"))
  # No restriction: all the composite states acceptable, each its own candidate (the lineage inherits the whole state of the branch)
  # (the inherited state "A|high" is also acceptable: the lineages grafted in this state exist for the following taxa)
  expect_equal(res$allowed_states$X3, c(composite, "A|high"))
  expect_equal(res$range_candidates$X3, list(`A|low` = "A|low", `AB|low` = "AB|low", `AB|high` = "AB|high", `A|high` = "A|high"))
  # Global values, replaced by the taxon-level values
  res_g <- resolve_trait_allowed_ranges(allowed_states_taxa = list(X1 = NULL, X2 = list(NULL, "low"), X3 = NULL), allowed_states = list("A", "high"), traits = traits,
                                        composite_states = composite, taxon_names = taxa, range_trait = 1)
  expect_equal(res_g$allowed_states$X1, c("AB|high", "A|high"))
  expect_equal(res_g$allowed_states$X2, c("A|low", "AB|low"))
  expect_equal(res_g$range_candidates$X2, list(`A|low` = "A|low", `AB|low` = "A|low"))
  # The range trait in another position
  traits_2 <- list(names = c("elevation", "biogeo"), states = list(elevation = c("low", "high"), biogeo = c("A", "AB")), sep = "|")
  res_2 <- resolve_trait_allowed_ranges(allowed_states_taxa = list(X1 = list("high", "A")), allowed_states = NULL, traits = traits_2,
                                        composite_states = c("low|A", "low|AB", "high|AB"), taxon_names = "X1", range_trait = 2)
  expect_equal(res_2$allowed_states$X1, c("high|AB", "high|A"))
  expect_equal(res_2$inherited_states, "high|A")
  # No observed combination: informative stop
  expect_error(resolve_trait_allowed_ranges(allowed_states_taxa = list(X1 = list("A", "high")), allowed_states = NULL, traits = traits, composite_states = c("A|low", "AB|low"),
                                            taxon_names = "X1", range_trait = 1), "No combination of the allowed states of taxon 'X1'.*biogeo = A; elevation = high")
  # A range contained in no state, or an area that is not in the map, stop naming the range trait
  expect_error(resolve_trait_allowed_ranges(allowed_states_taxa = list(X1 = list("C", NULL)), allowed_states = NULL, traits = traits, composite_states = composite,
                                            taxon_names = "X1", range_trait = 1), "areas.*\\(trait 'biogeo'\\)")
  expect_error(resolve_trait_allowed_ranges(allowed_states_taxa = list(X1 = list("B", NULL)), allowed_states = NULL, traits = traits, composite_states = c("A|low"),
                                            taxon_names = "X1", range_trait = 1), "Range trait 'biogeo':.*areas")
  expect_warning(resolve_trait_allowed_ranges(allowed_states_taxa = list(X1 = list(c("A", "B"), NULL)), allowed_states = NULL, traits = list(names = traits$names, states = list(biogeo = c("A", "AB", "B"), elevation = traits$states$elevation), sep = "|"),
                                              composite_states = c("A|low", "AB|low"), taxon_names = "X1", range_trait = 1), "Range trait 'biogeo':.*not states of the simmap")
})

test_that("graft_missing_taxa() treats one trait as ranges and the other traits as exact states: ranges inherited, the other traits kept from the graft point",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  taxa <- mtl$Taxon_name
  maps <- make_range_traits(toy$phy)
  ranges <- c(biogeo = TRUE, elevation = FALSE)

  # Allowed range "A" and elevation "high": taxa are grafted on AB|high, and inherit "A" and "high" ("A|high" is not observed on the backbone)
  grafted <- graft_ranges(toy, maps, treat_states_as_ranges = ranges, allowed_states = list("A", "high"), seed = 1)
  expect_s3_class(grafted, "phylo")
  expect_equal(names(grafted$simmap), c("biogeo", "elevation"))
  expect_true(all(state_of(grafted, "biogeo", taxa) == "A"))
  expect_true(all(state_of(grafted, "elevation", taxa) == "high"))
  expect_true(all(range_in_states_pairwise(state_of(grafted, "biogeo", taxa), parent_state_of(grafted$simmap$biogeo, taxa))))
  expect_true(all(parent_state_of(grafted$simmap$elevation, taxa) == "high"))
  # Both simmaps are valid on the grafted tree
  for (trait in names(grafted$simmap)) { expect_equal(length(grafted$simmap[[trait]]$maps), nrow(grafted$edge)); expect_s3_class(grafted$simmap[[trait]], "simmap") }
  expect_equal(grafted$simmap$biogeo$edge.length, grafted$edge.length)

  # Positional or named in any order: the same result
  same <- function (a, b) { expect_equal(a$edge, b$edge); expect_equal(a$simmap$biogeo$maps, b$simmap$biogeo$maps); expect_equal(a$simmap$elevation$maps, b$simmap$elevation$maps) }
  same(grafted, graft_ranges(toy, maps, treat_states_as_ranges = c(TRUE, FALSE), allowed_states = list("A", "high"), seed = 1))
  same(grafted, graft_ranges(toy, maps, treat_states_as_ranges = c(elevation = FALSE, biogeo = TRUE), allowed_states = list(elevation = "high", biogeo = "A"), seed = 1))
  same(grafted, graft_ranges(toy, rev(maps), treat_states_as_ranges = c(elevation = FALSE, biogeo = TRUE), allowed_states = list(elevation = "high", biogeo = "A"), seed = 1))
  # Unnamed traits: positional only
  unnamed <- graft_ranges(toy, unname(maps), treat_states_as_ranges = c(TRUE, FALSE), allowed_states = list("A", "high"), seed = 1)
  expect_equal(names(unnamed$simmap), c("trait_1", "trait_2"))
  expect_equal(unnamed$simmap$trait_1$maps, grafted$simmap$biogeo$maps)

  # Without ranges, the states are matched exactly: "A" cannot be grafted on "AB"
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = FALSE, allowed_states = list("A", "high"), seed = 1), "No combination")

  # Taxon-level values (list or text) replace the global ones: same results
  mtl_list <- mtl; mtl_list$allowed_states <- stats::setNames(list(list("A", "high"), NULL, NULL, NULL, NULL), taxa)
  mtl_text <- mtl; mtl_text$allowed_states <- stats::setNames(list("A;high", NA, NA, NA, NA), taxa)
  with_list <- graft_ranges(toy, maps, mtl = mtl_list, treat_states_as_ranges = ranges, allowed_states = list(NULL, "high"), seed = 3)
  with_text <- graft_ranges(toy, maps, mtl = mtl_text, treat_states_as_ranges = ranges, allowed_states = list(NULL, "high"), seed = 3)
  same(with_list, with_text)
  expect_equal(state_of(with_list, "biogeo", "X1"), c(X1 = "A"))
  expect_true(all(state_of(with_list, "elevation", taxa) == "high"))

  # No allowed ranges: the lineage inherits the whole state, the result is the one without ranges (same seed)
  free <- graft_ranges(toy, maps, treat_states_as_ranges = ranges, allowed_states = list(NULL, "high"), seed = 5)
  plain <- graft_ranges(toy, maps, allowed_states = list("AB", "high"), seed = 5)
  expect_true(same_tree(free, plain))
  expect_equal(free$simmap$biogeo$maps, plain$simmap$biogeo$maps)
  expect_true(all(state_of(free, "biogeo", taxa) == "AB"))
  # FALSE alone is the default, and a vector of FALSE too
  expect_true(same_tree(plain, graft_ranges(toy, maps, treat_states_as_ranges = c(FALSE, FALSE), allowed_states = list("AB", "high"), seed = 5)))
})

test_that("graft_missing_taxa() with ranges and traits: an allowed range that is not a state of the map is inherited and added to the simmap of the range trait, with a warning naming the trait",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  taxa <- mtl$Taxon_name
  maps <- make_range_traits(toy$phy)
  # "B" is not a state of the range map, but is contained in "AB"
  maps$biogeo <- relabel_maps(make_young_state_simmap(toy$phy, 20), c(A = "A", B = "AB"))
  mtl_b <- mtl
  expect_warning(grafted <- graft_missing_taxa(toy$phy, mtl_b, simmap = maps, treat_states_as_ranges = c(TRUE, FALSE), allowed_states = list("B", "high"), method = "uniform", verbose = FALSE, seed = 1),
                 "Range trait 'biogeo' ?:|trait_1|Range trait", fixed = FALSE)
  expect_true(all(state_of(grafted, "biogeo", taxa) == "B"))
  expect_true("B" %in% colnames(grafted$simmap$biogeo$mapped.edge))
  expect_false("B" %in% colnames(maps$biogeo$mapped.edge))
  expect_true(all(state_of(grafted, "elevation", taxa) == "high"))
  # No state of the elevation allows the taxa: informative stop
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = c(TRUE, FALSE), allowed_states = list("A", "Z"), seed = 1), "\\(trait 'elevation'\\).*'Z'")
  # An area that is not in the map
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = c(TRUE, FALSE), allowed_states = list("D", NULL), seed = 1), "areas.*\\(trait 'biogeo'\\)")
  # A trait with states that are not coded as ranges cannot be the range trait
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = c(FALSE, TRUE), seed = 1), "Range trait 'elevation'")
})

test_that("graft_missing_taxa() with ranges and traits: several simulations, paired multiSimmaps, BD_clades, and invalid arguments",
{
  toy <- toy_setup()
  mtl <- traits_missing_taxa(toy)
  taxa <- mtl$Taxon_name
  maps <- make_range_traits(toy$phy)
  ranges <- c(biogeo = TRUE, elevation = FALSE)
  # Several simulations
  sims <- graft_ranges(toy, maps, treat_states_as_ranges = ranges, allowed_states = list("A", "high"), nb_simulations = 3, seed = 2)
  expect_s3_class(sims, "multiPhylo")
  for (g in sims) { expect_true(all(state_of(g, "biogeo", taxa) == "A")); expect_true(all(state_of(g, "elevation", taxa) == "high")) }
  # Paired multiSimmaps: one result per pair
  pairs <- lapply(list(c(20, 10), c(15, 8)), function (ages) make_range_traits(toy$phy, ages[1], ages[2]))
  multi <- graft_ranges(toy, list(biogeo = as_multi(lapply(pairs, function (x) x$biogeo), c("p", "q")), elevation = as_multi(lapply(pairs, function (x) x$elevation), c("p", "q"))),
                        treat_states_as_ranges = ranges, allowed_states = list("A", "high"), seed = 2)
  expect_equal(names(multi), c("p", "q"))
  for (k in c("p", "q")) { expect_true(all(state_of(multi[[k]], "biogeo", taxa) == "A")); expect_true(all(state_of(multi[[k]], "elevation", taxa) == "high")) }
  # BD_clades
  bd <- suppressWarnings(graft_missing_taxa(toy$phy, mtl, simmap = maps, treat_states_as_ranges = ranges, allowed_states = list("A", "high"), method = "BD_clades", verbose = FALSE, seed = 4))
  expect_true(all(state_of(bd, "biogeo", taxa) == "A"))
  expect_true(all(state_of(bd, "elevation", taxa) == "high"))
  # Taxa grafted on the lineages of taxa grafted before: the inherited ("A|high") branches are available for the next grafts (no restriction of the elevation, any range)
  # Invalid arguments stop before any simulation
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = c(TRUE, TRUE)), "Only one trait")
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = c(TRUE, FALSE, FALSE)), "one value per trait")
  expect_error(graft_ranges(toy, maps, treat_states_as_ranges = c(color = TRUE, elevation = FALSE)), "exactly the names")
  expect_error(graft_ranges(toy, unname(maps), treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE)), "named but the traits of 'simmap' are not")
  expect_error(graft_ranges(toy, maps$biogeo, treat_states_as_ranges = c(TRUE, FALSE)), "must be TRUE or FALSE")
  expect_error(graft_missing_taxa(toy$phy, mtl, treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE), method = "uniform", verbose = FALSE), "must be TRUE or FALSE")
})
