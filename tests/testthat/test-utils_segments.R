
##### Test suite for helper functions handling branch sgements #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


test_that("get_branch_segments() returns one segment per branch covering the whole branch, in the order of the nodes",
{
  phy <- make_toy_phy()
  nodes <- c(5, 30, 12, 1)
  seg <- get_branch_segments(phy, nodes)
  expect_equal(nrow(seg), 4)
  expect_identical(seg$node, nodes)
  expect_identical(seg$from, rep(0, 4))
  expect_identical(seg$to, phy$edge.length[match(nodes, phy$edge[, 2])])
  expect_equal(nrow(get_branch_segments(phy, integer(0))), 0)
})

test_that("get_segments_crossing_age() keeps exactly the branches crossing the age, in the input order",
{
  phy <- make_toy_phy()
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  edge_ages <- root_age - phytools::nodeHeights(phy)
  nodes <- rev(phy$edge[, 2]) # not in edge order on purpose
  seg <- get_branch_segments(phy, nodes)
  for (age in c(0.5, 5, 20, 33.3, 49))
  {
    crossing <- get_segments_crossing_age(phy, seg, age, root_age = root_age)
    expected <- phy$edge[(edge_ages[, 1] >= age) & (edge_ages[, 2] <= age), 2]
    expect_setequal(crossing$node, expected)
    expect_identical(crossing$node, nodes[nodes %in% expected]) # input order kept
    # Same result when node heights are not provided or when edge ages are passed
    expect_identical(crossing, get_segments_crossing_age(phy, seg, age))
    expect_identical(crossing, get_segments_crossing_age(phy, seg, age, edge_ages = edge_ages))
  }
  expect_equal(nrow(get_segments_crossing_age(phy, seg, root_age + 1, root_age = root_age)), 0)
})

test_that("get_segments_crossing_age() handles segments shorter than their branch",
{
  phy <- make_toy_phy()
  root_age <- max(phytools::nodeHeights(phy)[, 2])
  node <- phy$edge[1, 2]
  len <- phy$edge.length[1]
  tip_age <- root_age - phytools::nodeheight(phy, node)
  seg <- data.frame(node = node, from = len * 0.25, to = len * 0.5)
  expect_equal(nrow(get_segments_crossing_age(phy, seg, tip_age + len * 0.4, root_age = root_age)), 1)
  expect_equal(nrow(get_segments_crossing_age(phy, seg, tip_age + len * 0.1, root_age = root_age)), 0)
  expect_equal(nrow(get_segments_crossing_age(phy, seg, tip_age + len * 0.9, root_age = root_age)), 0)
})

test_that("draw_point_in_segments() returns NULL without segments, and a node and a position within the selected segment otherwise",
{
  expect_null(draw_point_in_segments(data.frame(node = integer(0), from = numeric(0), to = numeric(0))))
  seg <- data.frame(node = c(3L, 7L, 9L), from = c(0, 2, 1), to = c(1, 5, 1.5))
  set.seed(1)
  for (i in 1:50)
  {
    d <- draw_point_in_segments(seg)
    k <- match(d$node, seg$node)
    expect_true(d$position >= seg$from[k] && d$position <= seg$to[k])
  }
  # A single segment is always selected
  expect_identical(draw_point_in_segments(seg[2, ])$node, 7L)
})

test_that("draw_point_in_segments() draws segments proportionally to their length (flat density over time)",
{
  seg <- data.frame(node = c(3L, 7L), from = c(0, 2), to = c(1, 5)) # lengths 1 and 3
  set.seed(42)
  nodes <- vapply(seq_len(4000), function (i) draw_point_in_segments(seg)$node, integer(1))
  expect_equal(mean(nodes == 7L), 0.75, tolerance = 0.05)
})
