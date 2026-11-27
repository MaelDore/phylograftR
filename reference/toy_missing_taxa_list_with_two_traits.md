# Toy list of missing taxa, with allowed states of two traits: colors and elevation

Same as
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md),
with the extra item `$allowed_states` covering the states of two traits,
[`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
(`"blue"`, `"green"`, `"red"`) and
[`phylograftR::toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
(`"high"`, `"low"`), in which each taxon can be grafted. It is used to
showcase the grafting constrained by several traits at once (see Several
traits in Details section of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).

It is the output of
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
applied to the table
`system.file("extdata", "toy_missing_taxa_with_two_traits.csv", package = "phylograftR")`.
In the column `allowed_states`, **the states of the traits are separated
by semicolons** (in the order of the maps provided in `simmap`: the
color, then the elevation), **and the states of a trait by commas**.
Open that file to see the different syntaxes:

- `"red;"`: the color is `"red"`, and the elevation has no taxon-level
  value (the global value is used).

- `"blue, green;high"`: color `"blue"` or `"green"`, and elevation
  `"high"`.

- `";low"`: no taxon-level value for the color, elevation `"low"`.

- `"elevation=high"`: traits named, so the color can be omitted (it has
  no taxon-level value).

- Blank cell: no taxon-level value for any trait, so the global
  `allowed_states` argument is used for both.

## Usage

``` r
data(toy_missing_taxa_list_with_two_traits)
```

## Format

A named list with the same 7 taxa and items as
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md),
and the extra item:

- `$allowed_states` List (one element per taxon, named after the taxa)
  of lists with the states of each trait (one element per trait, `NULL`
  for a trait without taxon-level value), or `NULL` to use the global
  `allowed_states` argument for all the traits. `"New_sp1"`:
  `list("red", NULL)`; `"New_sp4"`: `list(c("blue", "green"), "high")`;
  `"New_sp5"`: `list(NULL, "low")`; `"New_sp7"`:
  `list(elevation = "high")`; the others: `NULL`.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_missing_taxa_list_with_allowed_states`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_states.md)
  for the same list, with allowed states of colors.

- [`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
  for the stochastic character map of colors.

- [`phylograftR::toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
  for the stochastic character map of elevation.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa within a target clade based on colors and elevation,
  using this dataset.

## Examples

``` r
# Import the list
data(toy_missing_taxa_list_with_two_traits)
toy_missing_taxa_list_with_two_traits$allowed_states
#> $New_sp1
#> $New_sp1[[1]]
#> [1] "red"
#> 
#> $New_sp1[[2]]
#> NULL
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
#> [1] "blue"  "green"
#> 
#> $New_sp4[[2]]
#> [1] "high"
#> 
#> 
#> $New_sp5
#> $New_sp5[[1]]
#> NULL
#> 
#> $New_sp5[[2]]
#> [1] "low"
#> 
#> 
#> $New_sp6
#> NULL
#> 
#> $New_sp7
#> $New_sp7$elevation
#> [1] "high"
#> 
#> 

# Compare to the initial '.csv' file
system.file("extdata", "toy_missing_taxa_with_two_traits.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_missing_taxa_with_two_traits.csv"

# Graft with the maps of the two traits. The global 'allowed_states' has one element per trait:
# grafting on "blue" or "green" branches for the color, and no restriction for the elevation (NULL).
data(toy_backbone)
data(toy_simmap)
data(toy_simmap_elevation)

grafted_phy <- graft_missing_taxa(phy = toy_backbone,
    missing_taxa_list = toy_missing_taxa_list_with_two_traits,
    simmap = list(color = toy_simmap, elevation = toy_simmap_elevation),
    allowed_states = list(c("blue", "green"), NULL),
    verbose = FALSE, seed = 1234)

# One updated simmap per trait
names(grafted_phy$simmap)
#> [1] "color"     "elevation"

# State of each grafted taxon for each trait: should match the allowed states
taxon_names <- toy_missing_taxa_list_with_two_traits$Taxon_name
taxon_states <- sapply(grafted_phy$simmap, function (simmap)
{
  tip_edges <- match(match(taxon_names, simmap$tip.label), simmap$edge[, 2])
  setNames(vapply(simmap$maps[tip_edges], function (x) names(x)[1], character(1)), taxon_names)
})
taxon_states
#>         color  elevation
#> New_sp1 "red"  "low"    
#> New_sp2 "blue" "low"    
#> New_sp3 "blue" "high"   
#> New_sp4 "blue" "high"   
#> New_sp5 "blue" "low"    
#> New_sp6 "blue" "high"   
#> New_sp7 "blue" "high"   
toy_missing_taxa_list_with_two_traits$allowed_states
#> $New_sp1
#> $New_sp1[[1]]
#> [1] "red"
#> 
#> $New_sp1[[2]]
#> NULL
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
#> [1] "blue"  "green"
#> 
#> $New_sp4[[2]]
#> [1] "high"
#> 
#> 
#> $New_sp5
#> $New_sp5[[1]]
#> NULL
#> 
#> $New_sp5[[2]]
#> [1] "low"
#> 
#> 
#> $New_sp6
#> NULL
#> 
#> $New_sp7
#> $New_sp7$elevation
#> [1] "high"
#> 
#> 
```
