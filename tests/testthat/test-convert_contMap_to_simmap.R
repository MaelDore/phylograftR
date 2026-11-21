
##### Test suite for convert_contMap_to_simmap() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


# Local helpers: build a contMap of a BM trait on a toy phylogeny (the RNG of the session is restored)
make_toy_contMap <- function (phy, seed = 1, ...)
{
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  trait <- phytools::fastBM(phy)
  phytools::contMap(phy, trait, plot = FALSE, ...)
}

toy_cm_setup <- function (n_tips = 25, n_maps = 1)
{
  phy <- make_toy_phy(n_tips = n_tips, seed = 11, root_age = 50)
  list(phy = phy, cms = lapply(seq_len(n_maps), function (i) make_toy_contMap(phy, seed = i)))
}

states3 <- c("low", "medium", "high")


# --- Discretisation of one branch ---

test_that("discretise_edge merges equal neighbours and short segments", {
  # Equal neighbours merged, no minimum
  r <- discretise_edge(c(1, 1, 2, 1), c(1L, 1L, 2L, 2L), min_len = 0, terminal = FALSE)
  expect_equal(r$idx, c(1L, 2L))
  expect_equal(r$lens, c(2, 3))

  # A short segment is merged into its longest neighbour, then equal neighbours are merged
  r <- discretise_edge(c(5, 0.5, 4), c(1L, 2L, 1L), min_len = 1, terminal = FALSE)
  expect_equal(r$idx, 1L)
  expect_equal(r$lens, 9.5)

  # Shortest first, into the longest neighbour
  r <- discretise_edge(c(3, 0.2, 0.5, 6), c(1L, 2L, 3L, 1L), min_len = 1, terminal = FALSE)
  expect_equal(sum(r$lens), 9.7)
  expect_true(all(r$lens >= 1))

  # A branch shorter than the minimum has one state: the longest in total
  r <- discretise_edge(c(0.1, 0.3, 0.2), c(1L, 2L, 1L), min_len = 1, terminal = FALSE)
  expect_equal(r$idx, 1L)
  expect_equal(r$lens, 0.6)
  r <- discretise_edge(c(0.1, 0.3, 0.1), c(1L, 2L, 1L), min_len = 1, terminal = FALSE)
  expect_equal(r$idx, 2L)
})

test_that("discretise_edge keeps the state of the tip on terminal branches", {
  # Short terminal branch: state of the tip
  r <- discretise_edge(c(0.4, 0.3, 0.1), c(1L, 1L, 2L), min_len = 1, terminal = TRUE)
  expect_equal(r$idx, 2L)
  expect_equal(r$lens, 0.8)

  # The tip segment is never merged away, even if shorter than the minimum
  r <- discretise_edge(c(5, 0.2), c(1L, 2L), min_len = 1, terminal = TRUE)
  expect_equal(r$idx, c(1L, 2L))
  expect_equal(r$lens, c(5, 0.2))

  # Same segments on an internal branch: the short one is merged
  r <- discretise_edge(c(5, 0.2), c(1L, 2L), min_len = 1, terminal = FALSE)
  expect_equal(r$idx, 1L)
  expect_equal(r$lens, 5.2)
})


# --- Single contMap ---

test_that("A single contMap gives a valid simmap with the requested states", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  sm <- convert_contMap_to_simmap(cm, probs = c(1/3, 2/3), states = states3)

  expect_s3_class(sm, "simmap")
  expect_s3_class(sm, "phylo")
  expect_false(inherits(sm, "multiSimmap"))
  expect_equal(colnames(sm$mapped.edge), states3)
  expect_equal(length(attr(sm, "breaks")), 2)
  expect_equal(attr(sm, "states"), states3)

  # Same tree, durations sum to the branch lengths, states from the list
  expect_equal(sm$edge, cm$tree$edge)
  expect_equal(unname(vapply(sm$maps, sum, numeric(1))), sm$edge.length, tolerance = 1e-8)
  expect_true(all(unlist(lapply(sm$maps, names)) %in% states3))
  expect_equal(unname(rowSums(sm$mapped.edge)), sm$edge.length, tolerance = 1e-8)

  # Accepted by the package validation
  expect_no_error(validate_simmap(sm, toy$phy))
})

test_that("States follow the breaks in increasing order of the trait values", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  sm <- convert_contMap_to_simmap(cm, breaks = mean(cm$lims) + c(-0.5, 0.5), states = states3, min_segment_proportion = 0)
  vals <- extract_contMap_values(cm)$values
  for (e in seq_along(sm$maps))
  {
    expected <- states3[findInterval(vals[[e]], mean(cm$lims) + c(-0.5, 0.5)) + 1]
    got <- rep(names(sm$maps[[e]]), times = 1)
    # Each original segment lies in a new segment of its state: compare the sequence of distinct states
    expect_equal(got, rle(expected)$values)
  }
})

test_that("Tips keep the state of the tip data", {
  toy <- toy_cm_setup(n_tips = 30)
  cm <- toy$cms[[1]]
  sm <- convert_contMap_to_simmap(cm, probs = c(1/3, 2/3), states = states3, min_segment_proportion = 0.05)
  vals <- extract_contMap_values(cm)$values
  tip_edges <- match(seq_len(ape::Ntip(cm$tree)), cm$tree$edge[, 2])
  last_state_new <- vapply(sm$maps[tip_edges], function (m) names(m)[length(m)], character(1))
  last_state_orig <- vapply(vals[tip_edges], function (v) states3[findInterval(v[length(v)], attr(sm, "breaks")) + 1], character(1))
  expect_equal(last_state_new, last_state_orig)
})

test_that("min_segment_proportion removes the short segments and 0 keeps every crossing", {
  toy <- toy_cm_setup(n_tips = 40)
  cm <- toy$cms[[1]]
  sm_0 <- convert_contMap_to_simmap(cm, n_bins = 3, states = states3, min_segment_proportion = 0)
  sm_1 <- convert_contMap_to_simmap(cm, n_bins = 3, states = states3, min_segment_proportion = 0.01)
  sm_5 <- convert_contMap_to_simmap(cm, n_bins = 3, states = states3, min_segment_proportion = 0.2)
  n_seg <- function (sm) sum(vapply(sm$maps, length, integer(1)))
  expect_gte(n_seg(sm_0), n_seg(sm_1))
  expect_gte(n_seg(sm_1), n_seg(sm_5))
  expect_gt(n_seg(sm_0), n_seg(sm_5))

  # Internal branches: segments not shorter than the minimum, unless the branch is shorter than the minimum (single state)
  depth <- max(ape::node.depth.edgelength(cm$tree))
  min_len <- 0.01 * depth
  internal <- cm$tree$edge[, 2] > ape::Ntip(cm$tree)
  ok <- vapply(which(internal), function (e) length(sm_1$maps[[e]]) == 1 || all(sm_1$maps[[e]] >= min_len - 1e-10), logical(1))
  expect_true(all(ok))
})

test_that("Binning arguments produce the expected cut points", {
  toy <- toy_cm_setup(n_tips = 40)
  cm <- toy$cms[[1]]
  tips <- extract_contMap_values(cm)$tips

  sm_p <- convert_contMap_to_simmap(cm, probs = c(0.25, 0.75), states = states3)
  expect_equal(attr(sm_p, "breaks"), unname(quantile(tips, c(0.25, 0.75))))

  sm_n <- convert_contMap_to_simmap(cm, n_bins = 4, states = c("a", "b", "c", "d"))
  expect_equal(attr(sm_n, "breaks"), seq(min(tips), max(tips), length.out = 5)[2:4])

  mid <- mean(cm$lims)
  sm_b <- suppressWarnings(convert_contMap_to_simmap(cm, breaks = c(mid + 1, mid - 1), states = states3))
  expect_equal(attr(sm_b, "breaks"), c(mid - 1, mid + 1))
})


# --- Input checks ---

test_that("Exactly one binning argument is required", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  expect_error(convert_contMap_to_simmap(cm, states = states3), "Exactly one of")
  expect_error(convert_contMap_to_simmap(cm, breaks = c(-1, 1), n_bins = 3, states = states3), "Exactly one of")
  expect_error(convert_contMap_to_simmap(cm, probs = 0.5, n_bins = 3, states = c("a", "b")), "Exactly one of")
})

test_that("states must match the number of bins", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  expect_error(convert_contMap_to_simmap(cm, n_bins = 3, states = c("a", "b")), "one name per bin")
  expect_error(convert_contMap_to_simmap(cm, breaks = 0, states = states3), "one name per bin")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 3), "'states'")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 2, states = c("a", "a")), "unique")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 2, states = c("a", "")), "non-empty")
})

test_that("breaks are checked against the range of the trait values", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  lims <- cm$lims
  # Outside the range
  expect_error(convert_contMap_to_simmap(cm, breaks = c(mean(lims), lims[2] + 1), states = states3), "within the range")
  expect_error(convert_contMap_to_simmap(cm, breaks = c(lims[1] - 1, mean(lims)), states = states3), "within the range")

  # On a boundary: ignored with a warning
  expect_warning(sm <- convert_contMap_to_simmap(cm, breaks = c(lims[1], mean(lims), lims[2]), states = c("low", "high")), "ignored")
  expect_equal(attr(sm, "breaks"), mean(lims))

  # Still too many states after the removal
  expect_error(suppressWarnings(convert_contMap_to_simmap(cm, breaks = c(lims[1], mean(lims)), states = states3)), "one name per bin")
  # Only boundaries
  expect_error(suppressWarnings(convert_contMap_to_simmap(cm, breaks = lims, states = c("a", "b"))), "No interior cut point")

  # Other invalid breaks
  expect_error(convert_contMap_to_simmap(cm, breaks = c(0, NA), states = states3), "finite")
  expect_error(convert_contMap_to_simmap(cm, breaks = "a", states = states3), "numeric")
})

test_that("n_bins, probs, min_segment_proportion and verbose are checked", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  expect_error(convert_contMap_to_simmap(cm, n_bins = 1, states = c("a", "b")), "n_bins")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 2.5, states = c("a", "b", "c")), "n_bins")
  expect_error(convert_contMap_to_simmap(cm, probs = c(0, 0.5), states = states3), "probs")
  expect_error(convert_contMap_to_simmap(cm, probs = c(0.5, 1), states = states3), "probs")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 2, states = c("a", "b"), min_segment_proportion = 1), "min_segment_proportion")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 2, states = c("a", "b"), min_segment_proportion = -1), "min_segment_proportion")
  expect_error(convert_contMap_to_simmap(cm, n_bins = 2, states = c("a", "b"), verbose = NA), "verbose")
})

test_that("Invalid contMap inputs are rejected", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  expect_error(convert_contMap_to_simmap(NULL, n_bins = 2, states = c("a", "b")), "contMap")
  expect_error(convert_contMap_to_simmap(toy$phy, n_bins = 2, states = c("a", "b")), "contMap")
  expect_error(convert_contMap_to_simmap(list(cm, 1), n_bins = 2, states = c("a", "b")), "contMap")
  mcm <- structure(list(cm), class = "multiContMap")
  expect_error(convert_contMap_to_simmap(mcm, n_bins = 2, states = c("a", "b")), "multiContMap")
  bad <- cm; bad$lims <- c(1, 1)
  expect_error(convert_contMap_to_simmap(bad, n_bins = 2, states = c("a", "b")), "lims")
  bad <- cm; bad$cols <- NULL
  expect_error(convert_contMap_to_simmap(bad, n_bins = 2, states = c("a", "b")), "'\\$cols'")
})


# --- List of contMaps ---

test_that("A list of contMaps gives a named multiSimmap with shared cut points", {
  toy <- toy_cm_setup(n_maps = 3)
  sms <- convert_contMap_to_simmap(toy$cms, probs = c(1/3, 2/3), states = states3, verbose = FALSE)
  expect_s3_class(sms, "multiSimmap")
  expect_s3_class(sms, "multiPhylo")
  expect_length(sms, 3)
  expect_equal(names(sms), c("map_1", "map_2", "map_3"))
  expect_true(all(vapply(sms, inherits, logical(1), what = "simmap")))
  expect_length(attr(sms, "breaks"), 2)
  for (sm in sms) { expect_no_error(validate_simmap(sm, toy$phy)) }

  # Named, and partially named lists
  named <- setNames(toy$cms, c("a", "b", "c"))
  expect_equal(names(convert_contMap_to_simmap(named, n_bins = 2, states = c("x", "y"), verbose = FALSE)), c("a", "b", "c"))
  partial <- setNames(toy$cms, c("a", "", "c"))
  expect_equal(names(convert_contMap_to_simmap(partial, n_bins = 2, states = c("x", "y"), verbose = FALSE)), c("a", "map_2", "c"))

  # A list of one contMap gives a multiSimmap of one map
  one <- convert_contMap_to_simmap(toy$cms[1], n_bins = 2, states = c("x", "y"), verbose = FALSE)
  expect_s3_class(one, "multiSimmap")
  expect_length(one, 1)
})

test_that("The tip values used for the cut points are averaged across the list", {
  toy <- toy_cm_setup(n_maps = 2)
  # Same tips values in both: same cut points as for the single contMap
  same <- list(toy$cms[[1]], toy$cms[[1]])
  s1 <- convert_contMap_to_simmap(toy$cms[[1]], probs = 0.5, states = c("a", "b"))
  s2 <- convert_contMap_to_simmap(same, probs = 0.5, states = c("a", "b"), verbose = FALSE)
  expect_equal(attr(s1, "breaks"), attr(s2, "breaks"))
  expect_equal(s2[[1]]$maps, s1$maps)
})

test_that("The contMaps of a list must share the same tree", {
  toy <- toy_cm_setup(n_maps = 2)
  other_phy <- make_toy_phy(n_tips = 25, seed = 99, root_age = 50)
  other_phy$tip.label <- toy$phy$tip.label
  other <- make_toy_contMap(other_phy, seed = 5)
  expect_error(convert_contMap_to_simmap(list(toy$cms[[1]], other), n_bins = 2, states = c("a", "b")), "topology|branch lengths|tips")

  # Different tips
  phy2 <- make_toy_phy(n_tips = 20, seed = 11, root_age = 50)
  cm2 <- make_toy_contMap(phy2, seed = 5)
  expect_error(convert_contMap_to_simmap(list(toy$cms[[1]], cm2), n_bins = 2, states = c("a", "b")), "tips")

  # Different branch lengths
  scaled <- toy$cms[[2]]
  scaled$tree$edge.length <- scaled$tree$edge.length * 1.5
  expect_error(convert_contMap_to_simmap(list(toy$cms[[1]], scaled), n_bins = 2, states = c("a", "b")), "branch lengths")
})

test_that("The progress is printed every 10% for lists of several contMaps", {
  toy <- toy_cm_setup(n_maps = 1)
  cms <- rep(toy$cms, 10)
  msgs <- character(0)
  withCallingHandlers(
    convert_contMap_to_simmap(cms, n_bins = 2, states = c("a", "b"), verbose = TRUE),
    message = function (m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_length(msgs, 10)
  expect_match(msgs[1], "10%")
  expect_match(msgs[10], "100%")

  # Nothing when verbose = FALSE, or for a single contMap
  expect_no_message(convert_contMap_to_simmap(cms, n_bins = 2, states = c("a", "b"), verbose = FALSE))
  expect_no_message(convert_contMap_to_simmap(toy$cms[[1]], n_bins = 2, states = c("a", "b"), verbose = TRUE))
})


# --- Use in the grafting functions ---

test_that("The simmap and the multiSimmap from contMaps constrain the grafting", {
  data(toy_backbone, envir = environment())
  data(toy_missing_taxa_list, envir = environment())
  data(toy_simmap, envir = environment())
  # Deterministic traits built from the states of the tips of the toy simmap (1 = blue, 2 = green, 3 = red) plus a small variation
  tip_states <- phytools::getStates(toy_simmap, type = "tips")
  trait <- setNames(match(tip_states, c("blue", "green", "red")) + seq_along(tip_states) / 100, names(tip_states))
  trait <- trait[toy_backbone$tip.label]
  cms <- list(phytools::contMap(toy_backbone, trait, plot = FALSE),
              phytools::contMap(toy_backbone, trait + rev(seq_along(trait)) / 100, plot = FALSE))
  states <- c("blue", "green", "red")

  sm <- convert_contMap_to_simmap(cms[[1]], breaks = c(1.5, 2.5), states = states)
  grafted <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                simmap = sm, allowed_states = c("blue", "green"),
                                method = "BD_clades", verbose = FALSE, seed = 1234)
  expect_s3_class(grafted, "phylo")
  expect_true(all(toy_missing_taxa_list$Taxon_name %in% grafted$tip.label))

  sms <- convert_contMap_to_simmap(cms, breaks = c(1.5, 2.5), states = states, verbose = FALSE)
  grafted_multi <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                      simmap = sms, allowed_states = c("blue", "green"),
                                      method = "BD_clades", verbose = FALSE, seed = 1234)
  expect_equal(names(grafted_multi), c("map_1", "map_2"))
})


# --- Increasing order and different values of the cut points ---

test_that("Unsorted breaks are sorted with a warning", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  mid <- mean(cm$lims)
  expect_warning(sm <- convert_contMap_to_simmap(cm, breaks = c(mid + 1, mid - 1), states = states3), "increasing order")
  expect_equal(attr(sm, "breaks"), c(mid - 1, mid + 1))
  expect_no_warning(convert_contMap_to_simmap(cm, breaks = c(mid - 1, mid + 1), states = states3))
})

test_that("Duplicated breaks are removed with a warning", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  mid <- mean(cm$lims)
  expect_warning(sm <- convert_contMap_to_simmap(cm, breaks = c(mid, mid), states = c("a", "b")), "different values")
  expect_equal(attr(sm, "breaks"), mid)
  # Fewer bins than states: informative error after the warning
  expect_error(suppressWarnings(convert_contMap_to_simmap(cm, breaks = c(mid, mid), states = states3)), "one name per bin")
})

test_that("Unsorted or duplicated probs are corrected with a warning", {
  toy <- toy_cm_setup()
  cm <- toy$cms[[1]]
  expect_warning(sm <- convert_contMap_to_simmap(cm, probs = c(0.75, 0.25), states = states3), "increasing order")
  expect_equal(attr(sm, "breaks"), attr(convert_contMap_to_simmap(cm, probs = c(0.25, 0.75), states = states3), "breaks"))
  expect_warning(sm <- convert_contMap_to_simmap(cm, probs = c(0.5, 0.5), states = c("a", "b")), "different values")
  expect_length(attr(sm, "breaks"), 1)
})

test_that("Identical quantiles (tied tip values) are removed with a warning", {
  tied <- c(1, 1, 1, 1, 2, 3)
  expect_warning(b <- resolve_contMap_breaks(probs = c(0.25, 0.5), tip_values = tied, range_values = c(0, 4)), "different values.*share the same value")
  expect_equal(b, 1)
})


# --- Observed tip values (tip_values) ---

# Local helper: the trait simulated by make_toy_contMap() with the same seed (the RNG of the session is restored)
make_toy_trait <- function (phy, seed = 1)
{
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  phytools::fastBM(phy)
}

test_that("tip_values give cut points computed on the observed values and tips in the state of their observed value", {
  toy <- toy_cm_setup()
  x <- make_toy_trait(toy$phy, seed = 1)
  sm <- convert_contMap_to_simmap(toy$cms[[1]], probs = c(1/3, 2/3), states = states3, tip_values = x, verbose = FALSE)
  # Cut points: quantiles of the observed values, not of the values read from the contMap
  expect_equal(attr(sm, "breaks"), unname(stats::quantile(x, probs = c(1/3, 2/3))))
  # Every tip is in the state of its observed value
  expected <- unlist(convert_values_to_states(x, breaks = attr(sm, "breaks"), states = states3))
  expect_identical(phytools::getStates(sm, type = "tips")[names(x)], expected[names(x)])
  expect_no_error(validate_simmap(sm, toy$phy))
  # Bins of equal width on the observed range
  sm_n <- convert_contMap_to_simmap(toy$cms[[1]], n_bins = 3, states = states3, tip_values = x, verbose = FALSE)
  expect_equal(attr(sm_n, "breaks"), seq(min(x), max(x), length.out = 4)[2:3])
  # User cut points: the tips are in the state of their observed value too
  sm_b <- convert_contMap_to_simmap(toy$cms[[1]], breaks = c(-1, 1), states = states3, tip_values = x, verbose = FALSE)
  expect_identical(phytools::getStates(sm_b, type = "tips")[names(x)], unlist(convert_values_to_states(x, breaks = c(-1, 1), states = states3))[names(x)])
})

test_that("tip_values can be unnamed (in the order of the tips) or named in another order", {
  toy <- toy_cm_setup()
  x <- make_toy_trait(toy$phy, seed = 1)
  ref <- convert_contMap_to_simmap(toy$cms[[1]], probs = c(1/3, 2/3), states = states3, tip_values = x, verbose = FALSE)
  expect_equal(convert_contMap_to_simmap(toy$cms[[1]], probs = c(1/3, 2/3), states = states3, tip_values = unname(x), verbose = FALSE), ref)
  expect_equal(convert_contMap_to_simmap(toy$cms[[1]], probs = c(1/3, 2/3), states = states3, tip_values = rev(x), verbose = FALSE), ref)
})

test_that("tip_values are shared by all the contMaps of a list", {
  toy <- toy_cm_setup(n_maps = 3)
  x <- make_toy_trait(toy$phy, seed = 1)
  ms <- convert_contMap_to_simmap(toy$cms, probs = c(1/3, 2/3), states = states3, tip_values = x, verbose = FALSE)
  expect_s3_class(ms, "multiSimmap")
  expect_equal(attr(ms, "breaks"), unname(stats::quantile(x, probs = c(1/3, 2/3))))
  expected <- unlist(convert_values_to_states(x, simmap = ms))
  for (i in seq_along(ms)) { expect_identical(phytools::getStates(ms[[i]], type = "tips")[names(x)], expected[names(x)]) }
})

test_that("Without tip_values, the tip values are read from the contMap (default behavior unchanged)", {
  toy <- toy_cm_setup()
  expect_equal(convert_contMap_to_simmap(toy$cms[[1]], probs = c(1/3, 2/3), states = states3, tip_values = NULL, verbose = FALSE),
               convert_contMap_to_simmap(toy$cms[[1]], probs = c(1/3, 2/3), states = states3, verbose = FALSE))
})

test_that("tip_values are checked", {
  toy <- toy_cm_setup()
  x <- make_toy_trait(toy$phy, seed = 1)
  cm <- toy$cms[[1]]
  convert <- function (values) convert_contMap_to_simmap(cm, probs = c(1/3, 2/3), states = states3, tip_values = values, verbose = FALSE)
  expect_error(convert(as.character(x)), "numeric vector of finite values")
  expect_error(convert(replace(x, 1, NA)), "numeric vector of finite values")
  expect_error(convert(replace(x, 1, Inf)), "numeric vector of finite values")
  expect_error(convert(x[-1]), "one value per tip")
  expect_error(convert(setNames(x, c(names(x)[-1], "unknown_tip"))), "must be the tip labels")
  expect_error(convert(setNames(x, c(names(x)[-1], names(x)[2]))), "must be the tip labels") # Duplicated name
  expect_error(convert(x + 1000), "outside the range")
})
