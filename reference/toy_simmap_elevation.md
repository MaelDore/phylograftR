# Toy stochastic character map of elevation states

Stochastic character map (simmap) of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
with **two elevation states**, `"high"` and `"low"`, to combine with
[`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
(color) or
[`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
(biogeographic ranges) as examples of attribute-based constraints
applied on two traits at once.

It was simulated with
[`phytools::sim.history()`](https://rdrr.io/pkg/phytools/man/sim.history.html)
from a root in state `"low"`, with equal transition rates between the
two states (0.04 per unit of time). Among the first 1000 seeds, the one
used (94) is the first giving a map that is balanced and compatible with
the constraints of
[`phylograftR::toy_missing_taxa_list_with_two_traits`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_two_traits.md)
used with
[`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
(at least one branch in the allowed states of each taxon):

- The two states cover similar shares of the total branch length (54%
  `"low"` and 46% `"high"`).

- A moderate number of transitions: 12 changes of state, along 9 of the
  58 branches (no branch has more than 4 pieces).

- `"high"` is found in 18 of the 30 tips and `"low"` in 12.

## Usage

``` r
data(toy_simmap_elevation)
```

## Format

An object of classes `"simmap"` and `"phylo"` (see
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html)),
with the same tips, topology and branch lengths as
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md),
and the additional elements `$maps` (for each branch, the duration spent
in each state, from the rootward node to the tipward node),
`$mapped.edge` (total duration per state and branch), `$node.states` and
`$states` (states of the nodes and tips).

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
  for the stochastic character map of colors.

- [`phylograftR::toy_missing_taxa_list_with_two_traits`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_two_traits.md)
  for the list of missing taxa to graft, with allowed states of colors
  and elevation.

- [`phylograftR::toy_missing_taxa_list_with_ranges_and_elevation`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_ranges_and_elevation.md)
  for the list of missing taxa to graft, with allowed ranges and
  elevation.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa within a target clade based on elevation, using this
  dataset.

## Examples

``` r
# Import the data
data(toy_simmap_elevation)

# Plot the map
plot(toy_simmap_elevation, colors = c(high = "darkorange", low = "steelblue"),
     fsize = 0.7, ftype = "i")


# Duration spent in each state (sum on all the branches)
colSums(toy_simmap_elevation$mapped.edge)
#>      low     high 
#> 131.1994 110.8468 
```
