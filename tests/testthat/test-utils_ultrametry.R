test_that("check_backbone_ultrametric() accepts an ultrametric backbone and stops otherwise", {
  phy <- make_toy_phy()
  expect_true(check_backbone_ultrametric(phy))

  # No branch lengths
  no_length <- phy
  no_length$edge.length <- NULL
  expect_error(check_backbone_ultrametric(no_length), "no branch lengths")

  # Not ultrametric
  not_ultra <- phy
  not_ultra$edge.length[1] <- not_ultra$edge.length[1] * 1.5
  expect_error(check_backbone_ultrametric(not_ultra), "not ultrametric")
})

test_that("enforce_ultrametry_new_tips() only resets the terminal branches of the new tips", {
  toy <- toy_setup()
  node_B <- ape::getMRCA(toy$phy, toy$B)
  grafted <- graft_tip_on_branch(toy$phy, "X1", node_B, frac = 0.5)
  # Break the ultrametry of the new tip
  edge_X1 <- match(match("X1", grafted$tip.label), grafted$edge[, 2])
  broken <- grafted
  broken$edge.length[edge_X1] <- broken$edge.length[edge_X1] * 1.2
  expect_false(ape::is.ultrametric(broken))

  fixed <- enforce_ultrametry_new_tips(broken, tip_labels = "X1", root_age = toy$root_age)
  expect_true(ape::is.ultrametric(fixed))
  # All the other branches are untouched
  expect_equal(fixed$edge.length[-edge_X1], broken$edge.length[-edge_X1])
})

test_that("ensure_ultrametry() leaves ultrametric trees untouched and fixes small drifts", {
  phy <- make_toy_phy()
  expect_identical(ensure_ultrametry(phy), phy)

  drifted <- phy
  tip_edge <- match(1, drifted$edge[, 2])
  drifted$edge.length[tip_edge] <- drifted$edge.length[tip_edge] * (1 + 1e-4)
  expect_false(ape::is.ultrametric(drifted))
  expect_true(ape::is.ultrametric(ensure_ultrametry(drifted)))
})
