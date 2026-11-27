# Toy multiSimmap of a binary trait: stochastic maps of the elevation

Set of 10 stochastic character maps (`multiSimmap`) of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
for the **elevation** of
[`phylograftR::toy_tip_data`](https://maeldore.github.io/phylograftR/reference/toy_tip_data.md)
(`"high"` and `"low"`), inferred like
[`phylograftR::toy_multiSimmap`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap.md)
(models `"ER"` and `"ARD"` compared with their AICc; the **best model**
was the equal-rates model `"ER"`; seed 2345; 100 maps simulated, the
first 10 kept).

It has the same number of maps as
[`phylograftR::toy_multiSimmap`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap.md),
so that the two can be **combined as the maps of two traits**: the maps
of the traits are paired (map 1 of each trait, then map 2 of each trait,
and so on).

## Usage

``` r
data(toy_multiSimmap_elevation)
```

## Format

An object of classes `"multiSimmap"` and `"multiPhylo"`: a list of 10
objects of class `"simmap"`, with the same tips, topology and branch
lengths as
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md).

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_tip_data`](https://maeldore.github.io/phylograftR/reference/toy_tip_data.md)
  for the tip data used to build the maps.

- [`phylograftR::toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
  for a single map of the same states.

- [`phylograftR::toy_multiSimmap`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap.md)
  for the maps of the colors.

## Examples

``` r
# Import the data
data(toy_multiSimmap)
data(toy_multiSimmap_elevation)

# Maps of two traits: a list with one element per trait,
# each a multiSimmap with the same number of maps
maps_two_traits <- list(color = toy_multiSimmap, elevation = toy_multiSimmap_elevation)
sapply(maps_two_traits, length)
#>     color elevation 
#>        10        10 
```
