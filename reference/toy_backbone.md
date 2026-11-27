# Toy backbone phylogeny

Simulated ultrametric phylogeny with 30 tips (`t1` to `t30`) and a root
age of 50, used in the examples and vignettes.

Generated with [`ape::rcoal()`](https://rdrr.io/pkg/ape/man/rtree.html)
(seed 1234) and rescaled. Internal branches shorter than 1 (2% of the
root age) were lengthened by moving their rootward node towards the
root, so that no node looks like a polytomy on a plot (the topology, the
root age and the ultrametry are unchanged).

## Usage

``` r
data(toy_backbone)
```

## Format

An object of class `"phylo"` with 30 tips and 29 internal nodes.

## See also

Other datasets from phylograftR:

- [`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md)
  for the list of missing taxa to graft on the backbone.

- [`phylograftR::toy_taxonomic_tables`](https://maeldore.github.io/phylograftR/reference/toy_taxonomic_tables.md)
  for the taxonomic tables of the backbone and of the missing taxa.

- [`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
  for the stochastic character map of colors.

## Examples

``` r
# Import the phylogeny
data(toy_backbone)

# Plot the phylogeny
ape::plot.phylo(toy_backbone, cex = 0.7)

```
