# Toy multiSimmap of a categorical trait: stochastic maps of the colors

Set of 10 stochastic character maps (`multiSimmap`) of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
for the **colors** of
[`phylograftR::toy_tip_data`](https://maeldore.github.io/phylograftR/reference/toy_tip_data.md)
(`"blue"`, `"green"` and `"red"`). Each map is **one possible history**
of the color along the branches, so that together they account for the
**uncertainty** in the reconstruction of the ancestral states, unlike a
single `simmap`. All the maps have the **same states at the tips** (the
observed ones) and differ in the internal branches.

They were inferred from the tip data with
`deepSTRAPP::prepare_trait_data()` (see
[`vignette("prepare_stochastic_maps")`](https://maeldore.github.io/phylograftR/articles/prepare_stochastic_maps.md)):
the models of equal rates (`"ER"`), symmetrical rates (`"SYM"`) and all
rates different (`"ARD"`) were compared with their AICc. The **best
model** was the equal-rates model (`"ER"`), under which the maps were
simulated with
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html)
(seed 1234). 100 maps were simulated, and only the first 10 are kept to
keep the package light. 10 maps are enough for illustration, but
**actual analyses should use more maps (100 or more)**.

## Usage

``` r
data(toy_multiSimmap)
```

## Format

An object of classes `"multiSimmap"` and `"multiPhylo"`: a list of 10
objects of class `"simmap"` (see
[`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
for their structure), with the same tips, topology and branch lengths as
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md).
The class `"multiSimmap"` is **what tells the grafting functions that
the maps are alternative histories of the same trait** (see
[`vignette("prepare_stochastic_maps")`](https://maeldore.github.io/phylograftR/articles/prepare_stochastic_maps.md)).

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_tip_data`](https://maeldore.github.io/phylograftR/reference/toy_tip_data.md)
  for the tip data used to build the maps.

- [`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
  for a single map of the same states.

- [`phylograftR::toy_multiSimmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_elevation.md)
  to combine with these maps as the maps of two traits.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa with the maps as constraints: each map is used for a
  different simulation.

## Examples

``` r
# Import the data
data(toy_multiSimmap)
class(toy_multiSimmap)
#> [1] "multiSimmap" "multiPhylo" 
length(toy_multiSimmap)
#> [1] 10

# Plot the first two maps: two possible histories of the same trait
state_colors <- c(blue = "dodgerblue", green = "limegreen", red = "firebrick")
par(mfrow = c(1, 2))
phytools::plotSimmap(toy_multiSimmap[[1]], colors = state_colors, fsize = 0.6, ftype = "i")
phytools::plotSimmap(toy_multiSimmap[[2]], colors = state_colors, fsize = 0.6, ftype = "i")

par(mfrow = c(1, 1))
```
