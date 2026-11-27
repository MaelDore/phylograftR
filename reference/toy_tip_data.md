# Toy tip data of three species attributes

States of three **species attributes** for the 30 tips of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md):
a color, an elevation and a biogeographic range. They are the **tip
data** from which stochastic maps are inferred in
[`vignette("prepare_stochastic_maps")`](https://maeldore.github.io/phylograftR/articles/prepare_stochastic_maps.md)
(see
[`phylograftR::toy_multiSimmap`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap.md),
[`phylograftR::toy_multiSimmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_elevation.md)
and
[`phylograftR::toy_multiSimmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_ranges.md)).

They are the states found at the tips of
[`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md),
[`phylograftR::toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
and
[`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md),
so the datasets are consistent with all the other toy datasets.

## Usage

``` r
data(toy_tip_data)
```

## Format

A data.frame with 30 rows (the tips of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md),
as row names and in the same order) and three character columns:

- `$color`: `"blue"`, `"green"` or `"red"`.

- `$elevation`: `"high"` or `"low"`.

- `$range`: biogeographic ranges coded as in `BioGeoBEARS`, i.e., `"A"`,
  `"B"`, `"AB"` or `"BC"`.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_multiSimmap`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap.md),
  [`phylograftR::toy_multiSimmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_elevation.md)
  and
  [`phylograftR::toy_multiSimmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_ranges.md)
  for the stochastic maps inferred from these data.

- [`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
  for the backbone phylogeny.

## Examples

``` r
# Import the data
data(toy_tip_data)
head(toy_tip_data)
#>     color elevation range
#> t26  blue       low     B
#> t2   blue      high    AB
#> t21  blue      high    AB
#> t16 green      high    AB
#> t11  blue       low     A
#> t4   blue       low     A

# Named vector of tip states of one attribute, as required by deepSTRAPP::prepare_trait_data()
tip_colors <- setNames(toy_tip_data$color, row.names(toy_tip_data))
table(tip_colors)
#> tip_colors
#>  blue green   red 
#>    11    11     8 
```
