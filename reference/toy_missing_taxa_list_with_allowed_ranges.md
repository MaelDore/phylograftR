# Toy list of missing taxa, with allowed ranges

Same as
[`phylograftR::toy_missing_taxa_list_with_allowed_states`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_states.md),
with biogeographic ranges as allowed states: the `$allowed_states` item
holds ranges of
[`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md),
coded as in BioGeoBEARS, in which each taxon can be grafted. It is used
to showcase the grafting with biogeographic ranges as states, with
`treat_states_as_ranges = TRUE` (see the Details section "Biogeographic
ranges as states" of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).

It is the output of
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
applied to the table
`system.file("extdata", "toy_missing_taxa_with_allowed_ranges.csv", package = "phylograftR")`.
Four taxa have a taxon-level value (one or several ranges) for
`$allowed_states` and three have none (`NULL`), so they use the global
`allowed_states` argument of the grafting function. All the allowed
ranges are states of the toy simmap, so grafting does not raise
warnings. Some are contained in other states of the simmap (`"A"` in
`"AB"`, `"B"` in `"AB"` and `"BC"`): a lineage grafted in a branch in
`"AB"` with the allowed range `"A"` inherits `"A"`, which is a subset
sympatric speciation event, `(AB) -> (AB), (A)`. **Such a list requires
a `simmap` to graft.**

## Usage

``` r
data(toy_missing_taxa_list_with_allowed_ranges)
```

## Format

A named list with the same 7 taxa and items as
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md),
and the extra item:

- `$allowed_states` List (one element per taxon, named after the taxa)
  of character vectors of ranges of
  [`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
  in which the taxon can be grafted, or `NULL` to use the global
  `allowed_states` argument. `"New_sp1"`: `"A"`; `"New_sp4"`: `"B"` and
  `"BC"`; `"New_sp5"`: `"B"`; `"New_sp7"`: `"A"` and `"BC"`; the others:
  `NULL`.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_missing_taxa_list_with_allowed_states`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_states.md)
  for the same list, with allowed states of colors.

- [`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
  for the stochastic character map of biogeographic ranges.

- [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  to import such a list from an external table.

- [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  to build such a list from a data.frame.

## Examples

``` r
# Import the list
data(toy_missing_taxa_list_with_allowed_ranges)
toy_missing_taxa_list_with_allowed_ranges$allowed_states
#> $New_sp1
#> [1] "A"
#> 
#> $New_sp2
#> NULL
#> 
#> $New_sp3
#> NULL
#> 
#> $New_sp4
#> [1] "B"  "BC"
#> 
#> $New_sp5
#> [1] "B"
#> 
#> $New_sp6
#> NULL
#> 
#> $New_sp7
#> [1] "A"  "BC"
#> 

# Compare to the initial '.csv' file
system.file("extdata", "toy_missing_taxa_with_allowed_ranges.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_missing_taxa_with_allowed_ranges.csv"

# Graft with the toy simmap of ranges, treating the states as ranges
# Taxa with taxon-level ranges follow their specific constraints.
# Taxa without taxon-level ranges follow the global argument
# and can only be grafted on branches in a range containing "A" or "B".
data(toy_backbone)
data(toy_simmap_ranges)

grafted_phy <- graft_missing_taxa(phy = toy_backbone,
                                  missing_taxa_list = toy_missing_taxa_list_with_allowed_ranges,
                                  simmap = toy_simmap_ranges, allowed_states = c("A", "B"),
                                  treat_states_as_ranges = TRUE,
                                  verbose = FALSE, seed = 1234)

# Plot the updated simmap
plot_grafted(grafted_phy,
             simmap_colors = c("A" = "gold", "B" = "dodgerblue",
                               "AB" = "limegreen", "BC" = "firebrick2"))


# Inspect range inherited by each grafted taxon (state of its terminal branch in the updated simmap)
taxon_names <- toy_missing_taxa_list_with_allowed_ranges$Taxon_name
tip_edges <- match(match(taxon_names, grafted_phy$tip.label), grafted_phy$edge[, 2])
taxon_ranges <- setNames(object = vapply(grafted_phy$simmap$maps[tip_edges],
                                         function (x) names(x)[1], character(1)),
                         nm = taxon_names)
# Should match with the allowed ranges. NULL default to "A" or "B"
# based on the 'allowed_ranges' argument.
taxon_ranges
#> New_sp1 New_sp2 New_sp3 New_sp4 New_sp5 New_sp6 New_sp7 
#>     "A"     "A"     "A"     "B"     "B"     "B"    "BC" 
toy_missing_taxa_list_with_allowed_ranges$allowed_states
#> $New_sp1
#> [1] "A"
#> 
#> $New_sp2
#> NULL
#> 
#> $New_sp3
#> NULL
#> 
#> $New_sp4
#> [1] "B"  "BC"
#> 
#> $New_sp5
#> [1] "B"
#> 
#> $New_sp6
#> NULL
#> 
#> $New_sp7
#> [1] "A"  "BC"
#> 
```
