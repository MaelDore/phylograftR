# Toy multiSimmap of biogeographic ranges

Set of 10 stochastic character maps (`multiSimmap`) of
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
for the **biogeographic ranges** of
[`phylograftR::toy_tip_data`](https://maeldore.github.io/phylograftR/reference/toy_tip_data.md),
coded as in `BioGeoBEARS` (areas `"A"`, `"B"` and `"C"`; ranges of at
most two areas). Each map is one possible biogeographic history, with
anagenetic changes of range along the branches and cladogenetic events
at the nodes.

They were inferred from the tip data with the **DEC** model of
`BioGeoBEARS` (the only model fit), through
`deepSTRAPP::prepare_trait_data()` (see
[`vignette("prepare_stochastic_maps")`](https://maeldore.github.io/phylograftR/articles/prepare_stochastic_maps.md)):
biogeographic stochastic mapping (BSM) was run, and its output converted
into stochastic maps (seed 3456). 100 maps were simulated, and only the
first 10 are kept to keep the package light. 10 maps are enough for
illustration, but **actual analyses should use more maps (100 or
more)**.

The ranges of the maps include the ranges of the tips and may include
**other ranges only found at the internal nodes** (e.g., `"C"` or
`"AC"`). Use `treat_states_as_ranges = TRUE` to use the areas as allowed
states in the grafting functions (see
[`vignette("biogeographic_data")`](https://maeldore.github.io/phylograftR/articles/biogeographic_data.md)).

## Usage

``` r
data(toy_multiSimmap_ranges)
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

- [`phylograftR::toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
  for a single map of ranges.

- [`phylograftR::toy_missing_taxa_list_with_allowed_ranges`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_ranges.md)
  for missing taxa with allowed ranges.

## Examples

``` r
# Import the data
data(toy_multiSimmap_ranges)

# Ranges found in the first map
colnames(toy_multiSimmap_ranges[[1]]$mapped.edge)
#> [1] "_"  "A"  "B"  "C"  "AB" "AC" "BC"
```
