# Toy contMap of a continuous trait

A `contMap` of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
for the continuous trait of
[`phylograftR::toy_tip_values`](https://maeldore.github.io/phylograftR/reference/toy_tip_values.md),
built with
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html).
It is a **point estimate**: the ancestral values are estimated by
maximum likelihood under a Brownian Motion and **linearly interpolated
along the branches**, so it does not cover the uncertainty of the
reconstruction.

It is an example of input of
[`phylograftR::convert_contMap_to_simmap()`](https://maeldore.github.io/phylograftR/reference/convert_contMap_to_simmap.md),
which discretizes the trait values into states (e.g., `"low"`,
`"medium"`, `"high"`) to obtain a `simmap` that can **constrain the
grafting**. The resulting `simmap` is a **constraint, not a
reconstruction of a discrete trait**.

## Usage

``` r
data(toy_contMap)
```

## Format

An object of class `"contMap"` (see
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html))
with the elements `$tree` (the backbone as a `"simmap"`, in which the
states are the color levels of the trait values), `$cols` (the color
palette) and `$lims` (range of the trait values, about -7.5 to 4.2).

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_tip_values`](https://maeldore.github.io/phylograftR/reference/toy_tip_values.md)
  for the data used to build the contMap.

- [`phylograftR::toy_contMaps`](https://maeldore.github.io/phylograftR/reference/toy_contMaps.md)
  for 10 stochastic histories of the trait (list of contMaps).

- [`phylograftR::convert_contMap_to_simmap()`](https://maeldore.github.io/phylograftR/reference/convert_contMap_to_simmap.md)
  to discretize a contMap into a simmap.

## Examples

``` r
# Import the data
data(toy_contMap)
data(toy_tip_values)

# Plot the contMap
plot(toy_contMap, legend = FALSE, fsize = 0.7)


# Discretize into a simmap of three states using the terciles of the tip values.
# Provide the observed tip values: the contMap only approximates them
# (interpolation within segments)
my_simmap <- convert_contMap_to_simmap(toy_contMap, probs = c(1/3, 2/3),
                                       states = c("low", "medium", "high"),
                                       tip_values = toy_tip_values)
attr(my_simmap, "breaks")
#> [1] -0.06642759  1.49246937
plot(my_simmap, colors = c(low = "blue", medium = "grey", high = "red"),
     fsize = 0.7, ftype = "i")

```
