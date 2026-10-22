##### Test suite for helper functions handling biogeographic ranges #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


# ---- Coding of ranges ----

## Test that split_range() splits ranges in one area per letter, in the order of the ranges
test_that("split_range() splits ranges in one area per letter",
{
  expect_identical(split_range(c("A", "AB", "BCD")), list("A", c("A", "B"), c("B", "C", "D")))
  expect_identical(split_range(character(0)), list())
})

## Test that check_range_coding() accepts the BioGeoBEARS coding
test_that("check_range_coding() accepts single-area and sorted multi-area ranges of capital letters",
{
  checked <- check_range_coding(c("A", "AB", "BC", "ABCD", "Z"))
  expect_identical(checked$range, c("A", "AB", "BC", "ABCD", "Z"))
  expect_true(all(checked$valid))
  expect_true(all(checked$reason == ""))
  expect_identical(colnames(checked), c("range", "valid", "reason"))
})

## Test that check_range_coding() reports one specific reason per kind of problem
test_that("check_range_coding() gives the reason why a range is incorrect",
{
  checked <- check_range_coding(c("A", "ab", "BA", "AA", "A1", "red", "A_B", "", NA, "AC", "ACB"))
  reason <- stats::setNames(checked$reason, checked$range)
  expect_identical(checked$valid, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_match(reason[["ab"]], "lowercase letters \\(a, b\\)")
  expect_match(reason[["BA"]], "not in alphabetical order \\(expected 'AB'\\)")
  expect_match(reason[["AA"]], "repeated areas \\(A\\)")
  expect_match(reason[["A1"]], "characters that are not letters \\('1'\\)")
  expect_match(reason[["red"]], "lowercase letters")
  expect_match(reason[["A_B"]], "characters that are not letters \\('_'\\)")
  expect_match(checked$reason[checked$range %in% ""], "empty or NA")
  expect_match(checked$reason[is.na(checked$range)], "empty or NA")
  expect_match(reason[["ACB"]], "expected 'ABC'")
})

## Test that check_range_coding() reports all the reasons of a range and the order is not a reason when the letters are not capitals
test_that("check_range_coding() reports several reasons for the same range, and unsorted lowercase letters are only reported as lowercase",
{
  reason <- check_range_coding("BAB")$reason
  expect_match(reason, "repeated areas \\(B\\)")
  expect_match(reason, "not in alphabetical order \\(expected 'AB'\\)")
  expect_match(check_range_coding("aB")$reason, "lowercase letters \\(a\\)")
  expect_no_match(check_range_coding("ba")$reason, "alphabetical")
  # Capital letters are sorted independently of the locale
  expect_true(check_range_coding("ABZ")$valid)
  expect_false(check_range_coding("ZAB")$valid)
})

## Test that check_range_coding() checks each unique range once
test_that("check_range_coding() returns one row per unique range",
{
  checked <- check_range_coding(c("A", "ab", "A", "ab"))
  expect_identical(checked$range, c("A", "ab"))
})

## Test that validate_range_coding() stops and lists every incorrect state with its reason
test_that("validate_range_coding() returns the ranges invisibly if correct, and otherwise lists the incorrect states and why",
{
  expect_invisible(validate_range_coding(c("A", "AB", "B")))
  expect_identical(validate_range_coding(c("A", "AB", "B")), c("A", "AB", "B"))
  expect_error(validate_range_coding(c("A", "ba", "BA", "C"), who = "the simmap"),
               "With 'treat_states_as_ranges = TRUE', the states in the simmap must be biogeographic ranges coded as in BioGeoBEARS.*Incorrect states: 'ba' \\(lowercase letters \\(b, a\\)\\); 'BA' \\(areas not in alphabetical order \\(expected 'AB'\\)\\)")
  expect_error(validate_range_coding("red", who = "the global 'allowed_states' argument"), "the global 'allowed_states' argument.*'red'")
  expect_error(validate_range_coding(character(0)), "non-empty character vector")
  expect_error(validate_range_coding(1:2), "non-empty character vector")
  expect_error(validate_range_coding(NA_character_), "empty or NA")
})


# ---- Areas ----

## Test that get_range_areas() returns the sorted unique areas
test_that("get_range_areas() returns the unique areas in alphabetical order",
{
  expect_identical(get_range_areas(c("BC", "A", "AB")), c("A", "B", "C"))
  expect_identical(get_range_areas("A"), "A")
})

## Test that check_allowed_areas() stops with the areas available if an area is absent from the simmap
test_that("check_allowed_areas() stops when an allowed area is absent from the simmap and lists the areas available",
{
  states <- c("A", "B", "AB", "BC")
  expect_identical(check_allowed_areas(c("A", "BC", "AB"), states), c("A", "BC", "AB"))
  expect_error(check_allowed_areas(c("A", "DE"), states, who = "taxon 'X1'"), "taxon 'X1'.*'D', 'E'.*Areas available in the simmap: 'A', 'B', 'C'")
  # A multi-area range made of areas that are all in the simmap passes, even if no state contains it (handled when the ranges are expanded)
  expect_no_error(check_allowed_areas("AC", states))
})


# ---- Matching of ranges ----

## Test that range_in_states() tests if the range is contained in each state
test_that("range_in_states() is TRUE for the states containing all the areas of the range",
{
  states <- c("A", "B", "AB", "BC", "ABC", "C")
  expect_identical(range_in_states("A", states), c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_identical(range_in_states("AB", states), c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_identical(range_in_states("BC", states), c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_identical(range_in_states("ABC", states), c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  # A range is contained in itself, and not in a smaller state
  expect_true(range_in_states("AB", "AB"))
  expect_false(range_in_states("AB", "A"))
  expect_identical(range_in_states("A", character(0)), logical(0))
})

## Test that expand_allowed_ranges() returns the acceptable states and the candidate ranges of each
test_that("expand_allowed_ranges() gives the acceptable states and the allowed ranges contained in each",
{
  states <- c("A", "B", "AB", "BC")
  # A single area matches the states containing it
  expanded <- expand_allowed_ranges("A", states)
  expect_identical(expanded$acceptable_states, c("A", "AB"))
  expect_identical(expanded$candidates, list(A = "A", AB = "A"))
  expect_identical(expanded$inherited_ranges, "A")
  expect_identical(expanded$unmatched_ranges, character(0))

  # Several ranges: a state can contain several of them (candidates kept in the order of the allowed ranges)
  expanded <- expand_allowed_ranges(c("B", "A"), states)
  expect_identical(expanded$acceptable_states, c("A", "B", "AB", "BC"))
  expect_identical(expanded$candidates, list(A = "A", B = "B", AB = c("B", "A"), BC = "B"))
  expect_identical(expanded$inherited_ranges, c("B", "A"))

  # A multi-area range alone only matches the states containing the whole range, and it is inherited as such
  expanded <- expand_allowed_ranges("AB", c("A", "B", "AB", "ABC", "BC"))
  expect_identical(expanded$acceptable_states, c("AB", "ABC"))
  expect_identical(expanded$candidates, list(AB = "AB", ABC = "AB"))
  expect_identical(expanded$inherited_ranges, "AB")
})

## Test that expand_allowed_ranges() reports the ranges that match no state, and handles duplicates
test_that("expand_allowed_ranges() reports the allowed ranges contained in no state",
{
  states <- c("A", "B", "AB", "BC")
  # "AC" is made of areas of the simmap but no state contains it
  expanded <- expand_allowed_ranges(c("A", "AC"), states)
  expect_identical(expanded$acceptable_states, c("A", "AB"))
  expect_identical(expanded$unmatched_ranges, "AC")
  expect_identical(expanded$inherited_ranges, "A")

  # Nothing matches: no acceptable state
  expanded <- expand_allowed_ranges("AC", states)
  expect_identical(expanded$acceptable_states, character(0))
  expect_identical(expanded$candidates, stats::setNames(list(), character(0)))
  expect_identical(expanded$inherited_ranges, character(0))
  expect_identical(expanded$unmatched_ranges, "AC")

  # Duplicated ranges or states are used once
  expanded <- expand_allowed_ranges(c("A", "A"), c("A", "AB", "A"))
  expect_identical(expanded$acceptable_states, c("A", "AB"))
  expect_identical(expanded$candidates, list(A = "A", AB = "A"))
})

## Test that with single-area states and single-area ranges the expansion is the identity (same states as without areas)
test_that("expand_allowed_ranges() with single-area states keeps exactly the allowed states, each inherited as such",
{
  states <- c("A", "B", "C")
  expanded <- expand_allowed_ranges(c("C", "A"), states)
  expect_identical(expanded$acceptable_states, c("A", "C"))
  expect_identical(expanded$candidates, list(A = "A", C = "C"))
  expect_true(all(lengths(expanded$candidates) == 1))
})
