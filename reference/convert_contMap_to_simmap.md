# Convert contMap(s) of a continuous trait into simmap(s) with discrete states

Converts the continuous trait history(ies) mapped on a phylogeny into a
stochastic character map of class `"simmap"` with discrete states
obtained by cutting the trait values into bins (e.g., `"low"`,
`"medium"`, `"high"`).

The output can be used as a `simmap` in
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
and
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
to constrain the grafting of missing taxa by the continuous trait (e.g.,
a taxon known to have a high body size can only be grafted on branches
in the `"high"` state).

The function accepts a single contMap, typically built with
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html),
or a list of contMaps typically obtained from continuous stochastic
mapping approach such as the one featured in the `contsimmap` R package.
**A single `contMap` gives a single `simmap`.** A list of `contMap`s,
typically covering the uncertainty of the reconstruction (e.g., one
`contMap` per stochastic history of a continuous trait), gives a
**`multiSimmap`**, with **one `simmap` per `contMap`**.

## Usage

``` r
convert_contMap_to_simmap(
  contMap,
  breaks = NULL,
  n_bins = NULL,
  probs = NULL,
  states,
  tip_values = NULL,
  min_segment_proportion = 0.01,
  verbose = TRUE
)
```

## Arguments

- contMap:

  An object of class `"contMap"`, or a list of such objects. **All the
  `contMap`s of a list must be built on the same phylogeny** (same tips,
  topology, branch lengths and order of the edges). A list is named
  after the names of the `contMap`s, with `map_1`, `map_2`, ... if no
  names are provided.

- breaks:

  Numeric vector. The interior cut points between the bins, in the units
  of the trait. Values lower than the first cut point belong to the
  first state, and so on. They must be provided in increasing order, be
  different values, and lie within the range of the trait values (see
  Details). If not, a warning is sent and the cut points are sorted
  and/or the duplicated values removed. Exactly one of `breaks`,
  `n_bins` and `probs` arguments must be provided.

- n_bins:

  Integer, at least 2. The number of bins of equal width over the range
  of the trait values at the tips. Exactly one of `breaks`, `n_bins` and
  `probs` arguments must be provided.

- probs:

  Numeric vector with values strictly between 0 and 1. The quantile
  levels, computed on the trait values observed at the tips, at which to
  cut the bins. They must be provided in increasing order and be
  different values, otherwise a warning is sent and they are sorted
  and/or the duplicated values removed (e.g., `c(1/3, 2/3)` for three
  bins of equal numbers of tips). Exactly one of `breaks`, `n_bins` and
  `probs` arguments must be provided.

- states:

  Character vector. The names of the states, from the bin with the
  lowest values to the bin with the highest values. Its length must
  match the number of bins: **the number of interior cut points + 1**.

- tip_values:

  Numeric vector (optional). The **observed trait values at the tips**,
  named after the tips of the phylogeny (or provided in the order of the
  tips if unnamed). If provided, they are used instead of the values
  read from the `contMap`(s) to **compute the cut points** (`n_bins` and
  `probs`) and to **set the state of each tip** (the last segment of
  each terminal branch), so that the tips are in the state of their
  observed value. The same values are used for all the `contMap`s of a
  list. They must lie within the range of the trait values of the
  `contMap`(s). If `NULL`, the tip values are read from the
  `contMap`(s), where they are only approximated (see Details).

- min_segment_proportion:

  Numeric between 0 (included) and 1 (excluded). The minimum duration of
  a segment of a branch, as a proportion of the depth of the tree (root
  age). Default is `0.01`. Shorter segments are merged into their
  longest neighbour on the branch (see Details for why this may matter).
  Use `0` to keep every crossing of a cut point.

- verbose:

  Logical. Whether to print the progress of the conversion (every 10% of
  the `contMap`s) when a list of several `contMap`s is provided. Default
  is `TRUE`.

## Value

Depending on the nature of the input:

- If `contMap` is a single `contMap`: an object of class `"simmap"`.

- If `contMap` is a list of `contMap`s: an object of class
  `"multiSimmap"` in which each element is a `"simmap"` ordered to match
  the initial `contMap`s.

The interior cut points used to discretize the trait values are stored
in the attribute **`"breaks"`** of the output, and the names of the
states (from the lowest to the highest values) in its attribute
**`"states"`**. They allow to convert other values, e.g., the trait
values of the missing taxa, into the same states with
[`phylograftR::convert_values_to_states()`](https://maeldore.github.io/phylograftR/reference/convert_values_to_states.md).

## Details

### Principle

A `contMap` stores the trait history as a series of short segments along
each branch, with ancestral node values typically estimated from a
maximum likelihood reconstruction under a Brownian motion (BM) model and
linearly interpolated along the branches.

During discretization, each segment is assigned to the bin its value
falls in, and consecutive segments of the same bin are merged. **A value
equal to a cut point belongs to the bin above it.** The result is a
`simmap` (as in
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html))
where the states are the bins.

### Cutting the trait values

**Provide exactly one of the three arguments:**

- `breaks`: the cut points are chosen by the user.

- `n_bins`: uses bins of equal width over the range of the tip values.

- `probs`: uses bins defined by quantiles of the tip values.

`n_bins` and `probs` are computed on the tip values (averaged across the
`contMap`s of a list). By default, they are read from the `contMap`s,
which store the trait values on a finite scale of colors (from 0 to 1000
in
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)),
so they approximate the observed data. **Provide the observed values
with `tip_values` to use them instead** (see below).

### Tip values

A `contMap` stores the trait history as a series of short segments along
each branch, with the value of each segment interpolated from the
estimated values at the nodes. The value associated to a segment is a
single value for the whole segment, on a finite scale of colors. As a
result, the value of the last segment of a terminal branch is **nearly,
but not identical to the observed value at the tip**. With the default
settings, the function extrapolates the tip value from the last
segments, which is only an approximate. The difference is usually small
compared to the range of the trait, but it can move the tip of a taxon
with a value close to a cut point into the neighboring state, or shift
the cut points obtained from `n_bins` or `probs`.

To avoid this, provide the **observed tip values** with `tip_values`:
they are then used to compute the cut points, and **each tip is set in
the state of its observed value**, so that the tips of the `simmap`s are
in the same states as the observed data converted with
[`phylograftR::convert_values_to_states()`](https://maeldore.github.io/phylograftR/reference/convert_values_to_states.md).
This is recommended whenever the tip data used to build the `contMap`(s)
are available. When the observed state of a tip differs from the state
of the interpolated segments just before it, the terminal branch ends
with a short segment in the state of the observed value (the segment at
the tip is never merged away, see `min_segment_proportion`).

The values of the internal branches are reconstructed and are expected
to shrink towards the intermediate values under a BM process (or towards
one specific direction if a directional model was used to fit the
`contMap`s), because ancestral values are averaged over their
descendants. As a consequence, bins defined from the tips will contain a
smaller share of the internal branches for the extreme states than for
the intermediate states.

`breaks` must be interior cut points: they must lie within the range of
the trait values, i.e., between the minimum and the maximum of the
values of all the `contMap`s. **A cut point beyond this range is an
error.** A cut point equal to the minimum or the maximum is ignored with
a warning, because it would define an empty bin.

The cut points must be in increasing order and different values. This is
checked for the three arguments, including the cut points derived from
`n_bins` and `probs` (e.g., many tips sharing the same value can give
identical quantiles). Unsorted cut points are sorted, and duplicated
values are removed, with an informative warning.

**The number of `states` must be equal to the number of bins** (the
number of cut points + 1, after the cut points equal to the boundaries
were ignored). Removing cut points reduces the number of bins, so
`states` must match the number of bins left.

### Short segments

When using `contMap`s derived from continuous stochastic mapping
processes such as the one featured in the `contsimmap` R package, trait
values can wander around a cut point: the trait evolutionary history can
crosses the threshold several times over short distances, resulting in
tiny segments alternating states. **Segments shorter than
`min_segment_proportion` times the depth of the tree are merged**, from
the shortest to the longest, into their longest neighbour on the same
branch. A branch shorter than the minimum length gets a single state:
the one with the longest total duration, or the state of the tip for a
terminal branch. On a terminal branch, the segment at the tip is never
merged away to ensure tip states matched with observed tip data. Setting
`min_segment_proportion = 0` keeps all the crossings.

This is not an issue for `contMap` featuring linear interpolation of
ancestral trait values as in
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html),
since trait values can only vary monotonously along branches (i.e.,
decrease, increase, or remain stable, but not a combination of two),
thus cannot cross each cut point multiple times.

### Cautions for interpretation

- The resulting `simmap`s is/are a discretization of the reconstructed
  trait to constrain the grafting of missing taxa. **It is not a
  stochastic map of a discrete trait**: the number of transitions
  depends on the resolution of the `contMap` and on
  `min_segment_proportion` and should not be interpreted as a rate of
  evolution of discrete states.

- Nodes with reconstructed values close to a cut point have uncertain
  states. Providing a list of `contMap`s covering the uncertainty (see
  below) is the way to account for it, as each grafted phylogeny can
  then use a different `simmap`.

- The states depend on the cut points. A taxon whose allowed states are
  defined by a threshold can have different eligible branches if the
  threshold is moved. Ideally, you would want to compare several sets of
  cut points before relying on the results (see
  [`vignette("continuous_data")`](https://maeldore.github.io/phylograftR/articles/continuous_data.md)).

- As for the use of `simmap` derived from discrete states, constraining
  the grafting by a trait that will then be analyzed (e.g., the
  association of the trait with diversification) biases the results
  towards the constraint. Interpret such analyses accordingly.

### Accounting for uncertainty in trait estimates: continuous stochastic mapping

A single `contMap` is a point estimate. To cover the uncertainty of the
reconstruction, provide a list of `contMap`s, each one built from one
possible history of the trait. For example, the `contsimmap` R package
samples histories of a continuous trait on a phylogeny, and the function
`deepSTRAPP::convert_contsimmap_to_contMaps()` converts them into a list
of `contMap`s (see Examples). The output of this function is then a
**`multiSimmap`**, that can be used in the grafting functions like any
`multiSimmap`, to account for uncertainty in the ancestral trait
reconstructions.

The function `deepSTRAPP::prepare_trait_data()` is the **recommended
way** to obtain such a list of `contMap`s: with
`run_stochastic_maps = TRUE`, it fits the evolutionary model(s)
(selected by AICc), runs the continuous stochastic mapping with
`contsimmap::make.contsimmap()` and converts the maps into `contMap`s in
a single call, with additional options to fine-tune the process (see
Examples). The example below uses `contsimmap` directly to credit the
authors of the continuous stochastic mapping approach and to show the
steps usde internally by `deepSTRAPP::prepare_trait_data()`.

## See also

Other functions from phylograftR:

- [`phylograftR::convert_values_to_states()`](https://maeldore.github.io/phylograftR/reference/convert_values_to_states.md)
  to convert continuous values of the missing taxa into the states of
  the `simmap`.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades, using the
  states to constrain the grafting.

- [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft groups of missing taxa as sister clades, using the states to
  constrain the grafting.

## Author

Maël Doré

## Examples
