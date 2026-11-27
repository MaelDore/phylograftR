# Toy tip values of a continuous trait

Values of a continuous trait for the 30 tips of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md).
They are the data used to build
[`phylograftR::toy_contMap`](https://maeldore.github.io/phylograftR/reference/toy_contMap.md)
and
[`phylograftR::toy_contMaps`](https://maeldore.github.io/phylograftR/reference/toy_contMaps.md),
examples of **continuous traits discretized into states** to constrain
the grafting (see
[`phylograftR::convert_contMap_to_simmap()`](https://maeldore.github.io/phylograftR/reference/convert_contMap_to_simmap.md)).

They were simulated with
[`phytools::fastBM()`](https://rdrr.io/pkg/phytools/man/fastBM.html)
under a Brownian Motion (ancestral value of 0 at the root, evolutionary
rate `sig2 = 1`), with seed 37. Among the first 150 seeds, this is one
of the seeds for which the discretization of the contMap into three
states using the terciles of the tip values gives states that cover
similar shares of the branches (35% `"low"`, 43% `"medium"` and 22%
`"high"`), and constraints that differ among the 7 toy missing taxa of
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md).

## Usage

``` r
data(toy_tip_values)
```

## Format

A named numeric vector of length 30, named after the tips of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
and in the same order. The values range from about -7.5 to 4.2.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_contMap`](https://maeldore.github.io/phylograftR/reference/toy_contMap.md)
  for the contMap of this trait (point estimate).

- [`phylograftR::toy_contMaps`](https://maeldore.github.io/phylograftR/reference/toy_contMaps.md)
  for 10 stochastic histories of this trait (list of contMaps).

- [`phylograftR::convert_contMap_to_simmap()`](https://maeldore.github.io/phylograftR/reference/convert_contMap_to_simmap.md)
  to discretize a contMap into a simmap, and
  [`phylograftR::convert_values_to_states()`](https://maeldore.github.io/phylograftR/reference/convert_values_to_states.md)
  to convert continuous values into the same states.

## Examples

``` r
# Import the data
data(toy_tip_values)

# Distribution of the values
summary(toy_tip_values)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -7.5157 -0.2817  0.8178  0.5176  1.7173  4.2193 

# Terciles of the tip values: cut points used to define three states
quantile(toy_tip_values, probs = c(1/3, 2/3))
#>   33.33333%   66.66667% 
#> -0.06642759  1.49246937 
```
