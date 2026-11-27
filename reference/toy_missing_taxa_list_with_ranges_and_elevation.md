# Toy list of missing taxa, with allowed states of two traits: biogeographic ranges and elevation

Same as
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md),
with the extra item `$allowed_states` covering the states of two traits,
the biogeographic ranges of
[`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
(`"A"`, `"B"`, `"AB"`, `"BC"`) and the elevation of
[`phylograftR::toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
(`"high"`, `"low"`), in which each taxon can be grafted. It is used to
showcase the grafting constrained by a range trait and a 'normal' trait
at once, with
`treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE)` (see the
Details section "Biogeographic ranges as states" of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).

It is the output of
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
applied to the table
`system.file("extdata", "toy_missing_taxa_with_ranges_and_elevation.csv", package = "phylograftR")`.
In the column `allowed_states`, **the states of the traits are separated
by semicolons** (in the order of the maps provided in `simmap`: the
ranges, then the elevation), **and the states of a trait by commas**.
Open that file to see the different syntaxes:

- `"A;high"`: the range is `"A"`, and the elevation is `"high"`.

- `"B, BC;"`: the range is `"B"` or `"BC"`, and the elevation has no
  taxon-level value (the global value is used).

- `"elevation=low"`: traits named, so the ranges can be omitted (they
  have no taxon-level value).

- `"biogeo=A, BC; elevation=high"`: all the traits named, in any order.

- Blank cell: no taxon-level value for any trait, so the global
  `allowed_states` argument is used for both.

**The names of the traits (`biogeo` and `elevation`) must match the
names of the maps provided in `simmap`** to use the syntax with names.

## Usage

``` r
data(toy_missing_taxa_list_with_ranges_and_elevation)
```

## Format

A named list with the same 7 taxa and items as
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md),
and the extra item:

- `$allowed_states` List (one element per taxon, named after the taxa)
  of lists with the states of each trait (one element per trait, `NULL`
  for a trait without taxon-level value), or `NULL` to use the global
  `allowed_states` argument for all the traits. `"New_sp1"`:
  `list("A", "high")`; `"New_sp4"`: `list(c("B", "BC"), NULL)`;
  `"New_sp5"`: `list(elevation = "low")`; `"New_sp7"`:
  `list(biogeo = c("A", "BC"), elevation = "high")`; the others: `NULL`.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_ranges.md)
  for the same missing taxa list, with allowed biogeographic ranges.

- [`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
  for the stochastic character map of biogeographic ranges.

- [`phylograftR::toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
  for the stochastic character map of elevation.

- [`phylograftR::toy_missing_taxa_list_with_two_traits`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_two_traits.md)
  for the list with allowed states of colors and elevation.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa within a target clade based on ranges and elevation,
  using this dataset.

## Examples

``` r
# Import the list
data(toy_missing_taxa_list_with_ranges_and_elevation)
toy_missing_taxa_list_with_ranges_and_elevation$allowed_states
#> $New_sp1
#> $New_sp1[[1]]
#> [1] "A"
#> 
#> $New_sp1[[2]]
#> [1] "high"
#> 
#> 
#> $New_sp2
#> NULL
#> 
#> $New_sp3
#> NULL
#> 
#> $New_sp4
#> $New_sp4[[1]]
#> [1] "B"  "BC"
#> 
#> $New_sp4[[2]]
#> NULL
#> 
#> 
#> $New_sp5
#> $New_sp5$elevation
#> [1] "low"
#> 
#> 
#> $New_sp6
#> NULL
#> 
#> $New_sp7
#> $New_sp7$biogeo
#> [1] "A"  "BC"
#> 
#> $New_sp7$elevation
#> [1] "high"
#> 
#> 

# Compare to the initial '.csv' file
system.file("extdata", "toy_missing_taxa_with_ranges_and_elevation.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_missing_taxa_with_ranges_and_elevation.csv"

# Graft with the maps of the two traits. The first trait is treated as ranges,
# the second as 'normal', non-biogeographic states.
# Taxa without taxon-level states can be grafted on branches
# in a range containing "A" or "B" (global value),
# whatever the elevation (NULL).

# Load datasets
data(toy_backbone)
data(toy_simmap_ranges)
data(toy_simmap_elevation)

# Graft with constraints from elevation states and biogeographic ranges
grafted_phy <- graft_missing_taxa(phy = toy_backbone,
    missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
    simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
    allowed_states = list(c("A", "B"), NULL),
    treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
    verbose = FALSE, seed = 1234)

# Result stores one updated simmap per trait
names(grafted_phy$simmap)
#> [1] "biogeo"    "elevation"

# State of each grafted taxon for each trait: should match the allowed states
taxon_names <- toy_missing_taxa_list_with_ranges_and_elevation$Taxon_name
taxon_states <- sapply(grafted_phy$simmap, function (simmap)
{
  tip_edges <- match(match(taxon_names, simmap$tip.label), simmap$edge[, 2])
  setNames(vapply(simmap$maps[tip_edges], function (x) names(x)[1], character(1)), taxon_names)
})
taxon_states
#>         biogeo elevation
#> New_sp1 "A"    "high"   
#> New_sp2 "A"    "high"   
#> New_sp3 "A"    "high"   
#> New_sp4 "B"    "low"    
#> New_sp5 "B"    "low"    
#> New_sp6 "B"    "high"   
#> New_sp7 "BC"   "high"   
toy_missing_taxa_list_with_ranges_and_elevation$allowed_states
#> $New_sp1
#> $New_sp1[[1]]
#> [1] "A"
#> 
#> $New_sp1[[2]]
#> [1] "high"
#> 
#> 
#> $New_sp2
#> NULL
#> 
#> $New_sp3
#> NULL
#> 
#> $New_sp4
#> $New_sp4[[1]]
#> [1] "B"  "BC"
#> 
#> $New_sp4[[2]]
#> NULL
#> 
#> 
#> $New_sp5
#> $New_sp5$elevation
#> [1] "low"
#> 
#> 
#> $New_sp6
#> NULL
#> 
#> $New_sp7
#> $New_sp7$biogeo
#> [1] "A"  "BC"
#> 
#> $New_sp7$elevation
#> [1] "high"
#> 
#> 

# Check that all the grafted taxa respect their constraints, trait by trait
check_grafting_constraints(
    phy_grafted = grafted_phy,
    phy_backbone = toy_backbone,
    missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
    allowed_states = list(c("A", "B"), NULL),
    treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE))
#> 7 taxa checked
#>   - Missing from the grafted tree: 0
#>   - Target clade violated:          0
#>   - Excluded subclade violated:     0
#>   - Allowed states violated:        0 (biogeo: 0, elevation: 0)
#>   - Range not in graft-point state: 0
#>   - Backbone: topology preserved, branch lengths shifted by max 4.53e-09 (9.06e-11 of backbone depth; tolerance = 1e-06)
#>   => All constraints respected
#>           map simulation seed Taxon_name in_tree include_target_clade_stem
#> New_sp1 map_1          1 1234    New_sp1    TRUE                      TRUE
#> New_sp2 map_1          1 1234    New_sp2    TRUE                      TRUE
#> New_sp3 map_1          1 1234    New_sp3    TRUE                     FALSE
#> New_sp4 map_1          1 1234    New_sp4    TRUE                      TRUE
#> New_sp5 map_1          1 1234    New_sp5    TRUE                      TRUE
#> New_sp6 map_1          1 1234    New_sp6    TRUE                      TRUE
#> New_sp7 map_1          1 1234    New_sp7    TRUE                     FALSE
#>         target_position target_ok excluded_ok n_excluded_violated
#> New_sp1           crown      TRUE        TRUE                   0
#> New_sp2           crown      TRUE        TRUE                   0
#> New_sp3           crown      TRUE        TRUE                   0
#> New_sp4           crown      TRUE        TRUE                   0
#> New_sp5           crown      TRUE        TRUE                   0
#> New_sp6      tip_branch      TRUE        TRUE                   0
#> New_sp7           crown      TRUE        TRUE                   0
#>         excluded_violated states_ok states_ok_biogeo states_ok_elevation
#> New_sp1                NA      TRUE             TRUE                TRUE
#> New_sp2                NA      TRUE             TRUE                TRUE
#> New_sp3                NA      TRUE             TRUE                TRUE
#> New_sp4                NA      TRUE             TRUE                TRUE
#> New_sp5                NA      TRUE             TRUE                TRUE
#> New_sp6                NA      TRUE             TRUE                TRUE
#> New_sp7                NA      TRUE             TRUE                TRUE
#>         ranges_ok   ok
#> New_sp1      TRUE TRUE
#> New_sp2      TRUE TRUE
#> New_sp3      TRUE TRUE
#> New_sp4      TRUE TRUE
#> New_sp5      TRUE TRUE
#> New_sp6      TRUE TRUE
#> New_sp7      TRUE TRUE
```
