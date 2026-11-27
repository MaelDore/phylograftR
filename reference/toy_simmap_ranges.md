# Toy stochastic character map of biogeographic ranges of the toy backbone phylogeny

Stochastic character map (simmap) of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
with **four states that are biogeographic ranges** coded as in
BioGeoBEARS: `"A"`, `"B"`, `"AB"` and `"BC"` (areas A, B and C; the area
C is only found in `"BC"`). It can be used as `simmap` with
`treat_states_as_ranges = TRUE` to graft taxa on branches in ranges that
contain their allowed ranges (see the Details section "Biogeographic
ranges as states" of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
and the examples of
[`phylograftR::toy_missing_taxa_list_with_allowed_ranges`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_ranges.md)).

It was simulated with
[`phytools::sim.history()`](https://rdrr.io/pkg/phytools/man/sim.history.html)
from a root in range `"AB"`, with equal rates (0.04 per unit of time)
for the anagenetic changes of range that gain or lose a single area
(`"A"` \<-\> `"AB"`, `"B"` \<-\> `"AB"` and `"B"` \<-\> `"BC"`). A
change of two areas at once (`"AB"` \<-\> `"BC"`) is not allowed. As a
simulated history only has anagenetic changes, some of them were then
moved to the parental node of their branch, to also represent
cladogenetic events (see the Details section). Among the first 20000
seeds, the one used (3598) is the first giving a balanced map:

- The four ranges cover similar shares of the total branch length (about
  30% `"B"`, 26% `"AB"`, 25% `"A"` and 20% `"BC"`).

- A moderate number of changes of range: 10 of the 58 branches include a
  change of range (no branch has more than 3 pieces).

- Each range is found in 5 to 10 of the 30 tips.

## Usage

``` r
data(toy_simmap_ranges)
```

## Format

An object of classes `"simmap"` and `"phylo"` (see
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html)),
with the same tips, topology and branch lengths as
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md),
and the additional elements `$maps` (for each branch, the duration spent
in each range, from the rootward node to the tipward node),
`$mapped.edge` (total duration per range and branch), `$node.states` and
`$states` (ranges of the nodes and tips).

## Details

Events represented in the map, with the syntax of BioGeoBEARS.

**Anagenetic changes** (along the branches, 12 in total):

- 8 `"d"` (range expansion, gain of one area): 5 `"A"` -\> `"AB"` and 3
  `"B"` -\> `"BC"`.

- 4 `"e"` (range contraction, loss of one area): 3 `"AB"` -\> `"A"` and
  1 `"AB"` -\> `"B"`.

**Cladogenetic events** (at the 29 nodes, from the range of the lineage
above the node and the ranges at the start of the two descendant
branches):

- 5 `"s"` (subset sympatry), such as `(AB) -> (AB), (A)`: one descendant
  keeps the widespread range, the other starts in a single area within
  it. They are the events that are produced by the grafting functions
  when a taxon inherits an allowed range smaller than the range of the
  lineage it is grafted on (subset sympatric speciation, see
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).

- 1 `"v"` (vicariance), `(AB) -> (A), (B)`: the widespread range is
  split between the two descendants.

- 13 `"y"` (narrow sympatry): both descendants inherit the single-area
  range of the lineage, such as `(A) -> (A), (A)`.

- 10 nodes where both descendants inherit the same widespread range
  (widespread sympatry): `(AB) -> (AB), (AB)` or `(BC) -> (BC), (BC)`
  (allowed in BAYAREALIKE-type models, but not in DEC).

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_ranges.md)
  for the list of missing taxa to graft, with associated allowed
  biogeographic ranges.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa within a target clade based on geographic ranges, using
  this dataset.

## Examples

``` r
# Import the data
data(toy_simmap_ranges)

# Plot the map
range_colors <- c(A = "gold", B = "dodgerblue", AB = "limegreen", BC = "firebrick2")
phytools::plotSimmap(toy_simmap_ranges, colors = range_colors, fsize = 0.7, ftype = "i")


# Duration spent in each range (sum on all the branches)
colSums(toy_simmap_ranges$mapped.edge)
#>       AB        B        A       BC 
#> 62.15758 71.59366 59.53341 48.76147 
```
