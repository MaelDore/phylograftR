# Graft missing clades on a time-calibrated phylogeny

This function grafts missing taxa onto an ultrametric, time-calibrated
backbone phylogeny as clades sister to clades of the backbone. Missing
taxa sharing the same sister clade are **grafted together as one missing
clade**, on the stem branch of their sister clade.

The divergence time of the missing clade is drawn from a birth-death
model fitted on the sister clade, and the missing clade is grown under
the same model, conditioned on its number of taxa and its stem age.

## Usage

``` r
graft_missing_clades(
  phy,
  missing_taxa_list,
  min_tips = 5,
  return_BD_model_summary = FALSE,
  total_richness_per_clades = NULL,
  simmap = NULL,
  allowed_states = NULL,
  treat_states_as_ranges = FALSE,
  return_simmap = TRUE,
  min_valid_fraction = 0.1,
  max_age_redraws = 1000,
  seed = 1234,
  nb_simulations = 1,
  force_simulations = FALSE,
  parallel = NULL,
  save_temp = FALSE,
  save_temp_every = NULL,
  temp_file_path = NULL,
  clean_temp_files = TRUE,
  resume = FALSE,
  ultrametric_tol = .Machine$double.eps^0.5,
  verbose = TRUE
)
```

## Arguments

- phy:

  An object of class `"phylo"`. The backbone phylogeny on which taxa are
  grafted. It must be rooted, binary (no polytomies) and ultrametric
  (i.e., be time-calibrated). It can be the output of a previous
  grafting: the grafted tips then accumulate (see 'Successive
  graftings').

- missing_taxa_list:

  A list with the information on the taxa to graft. **Every item
  provided must be named after `$Taxon_name`.** Items:

  - `$Taxon_name`: Character vector of the names of the taxa to graft.
    **Names must be unique and absent from `phy$tip.label`.**

  - `$Sister_clade`: List (named after `$Taxon_name`) of character
    vectors. For each taxon, the tip labels defining the **sister
    clade**, i.e., the clade of the backbone that is sister to the
    missing clade. **Taxa with the same sister clade are grafted
    together as one missing clade.** A pair of tips defines the clade as
    their MRCA node. A single tip defines the "clade" as the terminal
    branch of this tip. **Clades must be defined by tips: node IDs are
    not accepted** (the function stops, and lists some of them). Use
    [`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md)
    to convert node IDs into pairs of tips.

  - `$Sister_clade_name`: (Optional) Character vector (named after
    `$Taxon_name`) of the names of the sister clades (`NA` if a clade
    has no name). **Taxa sharing the same sister clade must have the
    same name.** The names are not used for the grafting, but, if
    `total_richness_per_clades` is provided, they are checked against
    its `$Sister_clade_name` (a clade must have the same name in both
    objects).

  - `$allowed_states`: (Optional) List (named after `$Taxon_name`) of
    character vectors. For each taxon, the states of the `simmap` in
    which its missing clade can be grafted on the stem branch of its
    sister clade. It **overrides the global argument** `allowed_states`
    for this taxon. `NULL` or `NA` use the global argument. Only used
    with a `simmap`. Taxa with the same sister clade are grafted
    together as one missing clade, so **they must have the same allowed
    states** (the function stops otherwise, and lists the clades
    concerned). With the maps of several traits, each element is a list
    with the states of each trait (`NULL`: no taxon-level value for the
    trait), e.g.,
    `list(M1 = list(c("blue", "green"), "high"), M2 = list(NULL, "low"))`
    (see 'Several traits' in Details section).

- min_tips:

  Integer. Minimum number of tips of the clade used to fit the
  birth-death model. If the sister clade has fewer tips, the smallest
  enclosing clade with at least `min_tips` tips is used instead. **Must
  be at least 3.** Default is `5`.

- return_BD_model_summary:

  Logical. Whether to report the birth-death models fitted during the
  grafting. If `TRUE`, each grafted phylogeny has the attribute
  `"BD_model_summary"`: a data.frame with one row per model (i.e., per
  missing clade), summarizing diversification rates and sampling
  fraction. It can be extracted a posteriori with
  [`phylograftR::extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md).
  Default is `FALSE`.

- total_richness_per_clades:

  (Optional) List with the total number of species of the clades, to
  account for species that are neither in the backbone nor in
  `missing_taxa_list` (i.e., unsampled species) when fitting the
  birth-death models. Items:

  - `$Sister_clade`: List of character vectors of tips, defining clades
    as in `missing_taxa_list$Sister_clade`. **Node IDs are not
    accepted**: convert them with
    [`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md)
    in the same call as `missing_taxa_list`, so the tips of the clades
    match between both objects. The easiest way to obtain this list is
    to build it from a table of clades with
    [`phylograftR::build_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/build_total_richness_per_clades.md)
    or
    [`phylograftR::import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md):
    the clades can be written with names, tips or node IDs, they are
    matched with the clades of `missing_taxa_list`, and the richness is
    checked against all the taxa to graft. **It must cover all the
    unique sister clades of `missing_taxa_list`.** It can include other
    clades, e.g., enclosing clades.

  - `$Sister_clade_name`: (Optional) Character vector of the names of
    the clades, same length as `$Sister_clade` (`NA` if a clade has no
    name). The names must be unique. If
    `missing_taxa_list$Sister_clade_name` is provided, **a clade must
    have the same name in both objects** (the function stops otherwise).
    They are not used for the grafting, but they are reported by
    [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md).

  - `$total_richness`: Vector of positive integers, same length as
    `$Sister_clade`. Total number of species of each clade, i.e.,
    backbone species, missing species nested within the clade, and
    unsampled species. Use `NA` if unknown: the unsampled species are
    then assumed to be those of the clades with an informed total
    richness nested within the clade (none if there are no such clades).

  With `NULL` (default), missing species to graft and backbone species
  are assumed to cover all species in the focal group (no unsampled
  species). If the list is provided with only `NA` in `$total_richness`,
  the fits are identical (no unsampled species in any clade,
  `total_richness` = backbone + missing species).

  To obtain a posteriori the species richness and the effective sampling
  fraction of the clades in each grafted phylogeny, use
  [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  on the output, with the same `total_richness_per_clades`.

  See 'Sampling fraction and unsampled species' in Details section.

- simmap:

  (Optional) A stochastic character map of the backbone, of class
  `"simmap"` (e.g., from
  [`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html)
  or
  [`phytools::sim.history()`](https://rdrr.io/pkg/phytools/man/sim.history.html)),
  or a set of maps of the same trait, of class `"multiSimmap"` (e.g.,
  from `phytools::make.simmap(nsim = 10)`), used as alternative
  hypotheses of the evolution of a species attribute (e.g., a trait or a
  range) with discrete states along the branches of `phy`.

  Alternatively, you can provide a list of `"simmap"` and/or
  `"multiSimmap"` objects, one per trait, that contains the maps of
  several traits used together as constraints All maps must have the
  same tips, topology and branch lengths as `phy`. With a `simmap`, the
  missing clades can only be grafted on the parts of the stem branches
  of their sister clades that are in their allowed states (see
  `allowed_states`), and each missing clade inherits the state(s) of the
  point of the stem branch where it is grafted, along all its branches,
  so the states can be followed alongside the phylogeny.

  See 'Species attributes', 'Multiple stochastic maps' and 'Several
  traits' in Details section. Default is `NULL` (no restriction for
  grafting based on states).

- allowed_states:

  Character vector. The states of `simmap` in which grafting is allowed.
  This is the default for all taxa: taxa with a value in
  `missing_taxa_list$allowed_states` use this specific value instead.
  States are matched exactly with the states of the `simmap`, unless
  `treat_states_as_ranges = TRUE`. Taxa with no allowed states (neither
  a global nor a taxon-level value) can be grafted in all the states of
  the `simmap`. The function stops if some states are not in the
  `simmap`, or if `allowed_states` is provided without a `simmap`. With
  `treat_states_as_ranges = TRUE`, a range that is not a state of the
  `simmap` is accepted if it is contained in some of them (e.g., `"A"`
  can be allowed if `"AB"` is found on the `simmap`).

  With the maps of several traits (a list in `simmap`), you must provide
  a list with one element per trait, in the order of the traits or named
  after them. In each element, the states allowed for each trait are
  named, with `NULL` for no restriction on a trait, e.g.,
  `list(c("blue", "green"), NULL, "high")` (see 'Several traits' in
  Details section). Default is `NULL` (all states allowed).

- treat_states_as_ranges:

  Logical. Only used with a `simmap`. Whether the states of the `simmap`
  are biogeographic ranges to be treated as sets of areas, rather than
  as unrelated states. Ranges must be coded as in `BioGeoBEARS`: one
  capital letter per area, and the letters of a multi-area range in
  alphabetical order (e.g., `"A"`, `"AB"`, `"ABC"`). If `TRUE`, an
  allowed range enables grafting on branch pieces in any range that
  contains it: `"A"` can be grafted on `"A"`, `"AB"` and `"ABC"`,
  whereas `"AB"` can be grafted on `"AB"` and `"ABC"` but not `"A"`.

  If stochastic maps of several traits are provided (a list in
  `simmap`), you must provide a logical vector with one value per trait,
  with at most one `TRUE` to define the trait treated as biogeographic
  ranges, the other traits being matched exactly. The vector is
  positional (e.g., `c(FALSE, TRUE, FALSE)`: the second trait is the
  range trait) or named after the traits, in any order (e.g.,
  `c(biogeo = TRUE, color = FALSE, elevation = FALSE)`). A single
  `FALSE` is accepted with several traits (no range trait), but a single
  `TRUE` stops.

  See 'Biogeographic ranges as states' and 'Several traits' in Details
  section. Default is `FALSE` (states are matched exactly).

- return_simmap:

  Logical. Only used with a `simmap`, and silently ignored without it.
  Whether to add the updated `simmap` to the output, in `$simmap` (see
  'Value'). Default is `TRUE`.

- min_valid_fraction:

  Numeric between 0 and 1. Only used with a `simmap`. Threshold of the
  fraction of the time window of the divergence time of a missing clade
  (the stem branch of its sister clade) that is in its allowed states. A
  warning lists the missing clades below the threshold for inspection,
  as their divergence times are concentrated on a small part of the stem
  branch. The fraction is computed on the backbone, before grafting.
  Default is `0.1`.

- max_age_redraws:

  Integer. Only used with a `simmap`. Maximum number of times the
  divergence time of a missing clade can be drawn again when it falls in
  a part of the stem branch that is not in an allowed state. The
  function stops with an informative message if this limit is reached.
  Default is `1000`.

- seed:

  Integer. Seed used for reproducibility. The seed only applies within
  the function: the state of the random number generator of the R
  session is restored on exit. Set to `NULL` to use the random number
  generator of the session. Default is `1234`. With several simulations,
  it is the master seed: the first simulation uses `seed`, and the seeds
  of the other simulations are derived from it (see 'Multiple
  simulations').

- nb_simulations:

  Integer. Number of simulations, i.e., of grafted phylogenies to
  produce in the same run, from the same backbone and
  `missing_taxa_list`. The grafting is stochastic, so ideally downstream
  analyses should be run on a set of grafted phylogenies rather than on
  a single one. Default is `1`. See 'Multiple simulations' in Details
  section for handling of the output, the seeds and the failures.

- force_simulations:

  Logical. Whether to try new seeds until `nb_simulations` grafted
  phylogenies are obtained, up to `10 * nb_simulations` attempts (see
  'Multiple simulations' in Details section). If `FALSE`, the
  simulations that fail are skipped. Default is `FALSE`.

- parallel:

  Logical or `NULL`. Whether to run the simulations in parallel, with
  the packages `future` and `future.apply` (**they must be installed
  separately from phylograftR**). Default is `NULL`: the simulations are
  run in parallel if the user set a `future` plan with several workers
  beforehand (e.g., `future::plan(future::multisession, workers = 4)`),
  and sequentially otherwise. `FALSE` forces a sequential run, even if a
  plan is set. `TRUE` execute the `future` plan, but stops with an
  informative error if the packages are not installed or no plan is set.
  See 'Parallel computation' in Details section for recommendations on
  whether parallelization is useful in your case.

- save_temp:

  Logical. Whether to save the grafted phylogenies during the process,
  to recover from a crash or an interruption of long runs (see
  'Temporary saves and resuming' in the Details section). Each phylogeny
  (one per simulation and per stochastic map) has its own temporary
  file. Its progress is saved every `save_temp_every` missing clades,
  and the phylogeny is always saved once fully grafted. Default is
  `FALSE`.

- save_temp_every:

  Integer or `NULL`. Frequency of the saves of the progress, in number
  of missing clades grafted (groups of taxa with the same sister clade).
  Only used if `save_temp = TRUE`. Default is `NULL`: save each
  phylogeny only once fully grafted.

- temp_file_path:

  Character string. Path of the `.rds` file used for the temporary
  saves. It must end with `.rds` and its directory must exist and be
  writable. Default is `NULL`, which saves in
  `"temp_phy_from_grafting.rds"` in the working directory. With several
  simulations and/or several stochastic maps, the path is a template:
  the suffixes `"_map_<name>"` (several maps; `"map_"` is not repeated
  if the name of the map already starts with it, e.g.,
  `"temp_map_1.rds"` for a map named `"map_1"`) and `"_sim_<k>"`
  (several simulations) are added before `".rds"`, in this order (e.g.,
  `"temp_map_A_sim_2.rds"`). See 'Temporary saves and resuming' in the
  Details section.

- clean_temp_files:

  Logical. Whether to delete all the temporary files of the run once it
  is completed without stopping/crashing. Only used if
  `save_temp = TRUE`. If the run crashes or is interrupted, nothing is
  deleted so it can be resumed. Use `FALSE` to keep the files. Default
  is `TRUE`.

- resume:

  Logical. Whether to resume an interrupted run from its temporary
  files, which **must come from a run with the exact same arguments**
  (see 'Temporary saves and resuming'). Needs `save_temp = TRUE`.
  Default is `FALSE`.

- ultrametric_tol:

  Numeric. Tolerance used by
  [`ape::is.ultrametric()`](https://rdrr.io/pkg/ape/man/is.ultrametric.html)
  for the ultrametricity checks of the backbone and of the final
  phylogeny. Default is `.Machine$double.eps^0.5`.

- verbose:

  Logical. Whether to display progress. Messages print every 10 missing
  clades. If computation is run in parallel, only a message is printed
  when each simulation is done. Default is `TRUE`.

## Value

The function returns an object of class `"phylo"`: the backbone
phylogeny with the grafted taxa added as new tips. Backbone branch
lengths are untouched, and the grafted phylogeny is ultrametric.

In addition to the standard elements of a `"phylo"` object, the
phylogeny has the following elements:

- `grafted_tips` Named logical vector with one value per tip, in the
  order of `$tip.label`. `TRUE` for grafted tips (including the tips
  grafted by a previous grafting if `phy` was already grafted). Names
  are the tip labels.

- `grafted_edges` Named logical vector with one value per branch, in the
  order of the rows of `$edge`. `TRUE` for grafted branches, i.e.,
  branches with only grafted taxa descending from them. Names are the
  IDs of the tipward nodes of the branches (`$edge[, 2]`). Use them to
  color the phylogeny on plots (see examples).

**These elements are not updated by functions that modify the
phylogeny** (e.g.,
[`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html), which
removes tips but keeps the other elements as they are, and
[`ape::write.tree()`](https://rdrr.io/pkg/ape/man/write.tree.html),
which removes them). Use
[`phylograftR::identify_grafted()`](https://maeldore.github.io/phylograftR/reference/identify_grafted.md)
to recompute `grafted_tips` and `grafted_edges` on a manipulated
phylogeny.

Optional elements can be added to the output based on what optional
argument was provided:

- `"BD_model_summary"` (attribute). Only if
  `return_BD_model_summary = TRUE`. Data.frame with one row per
  birth-death model fitted (one per missing clade), summarizing
  information about the clade used for model fitting, diversification
  rates, and sampling fraction. It can be extracted a posteriori with
  [`phylograftR::extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md).

- `simmap` Only if a `simmap` is provided and `return_simmap = TRUE`.
  The updated stochastic character map of the grafted phylogeny (class
  `"simmap"`): the pieces of the branches of the backbone keep their
  states (cut where taxa were grafted), and every grafted clades is
  entirely in the state of the point where it was grafted. With
  `treat_states_as_ranges = TRUE`, each grafted clade is entirely in the
  range it inherited (see 'Biogeographic ranges as states' in Details
  section). This element can be plotted as any `simmap`, e.g.,
  `plot(grafted_phy$simmap)`, or used to extract the states of the tips
  with `phytools::getStates(grafted_phy$simmap, type = "tips")`. It has
  the same exact topology and branch lengths as the grafted phylogeny.
  With the maps of several traits, `$simmap` is a list of such simmaps,
  one per trait, named after the traits (see 'Several traits' in Details
  section).

With `nb_simulations > 1`, the function returns an object of class
`"multiPhylo"` with one grafted phylogeny per simulation, named
`"sim_1"`, `"sim_2"`, ... Each phylogeny has the elements described
above (including `$simmap` if a `simmap` was provided and
`return_simmap = TRUE`). The `"multiPhylo"` object has the attributes
`"seeds"` (named vector of the seeds of the phylogenies, to reproduce
each of them alone with `nb_simulations = 1`; `NULL` if `seed = NULL`)
and `"failed"` (data.frame with the seed and the reason of the
simulations that failed, 0 row if none).

With multiple stochastic maps (i.e., a `"multiSimmap"`) provided as
`simmap`, the function returns a list with one element per map (named
after the maps): a phylogeny if `nb_simulations = 1`, a `"multiPhylo"`
otherwise (see 'Multiple stochastic maps'). The list has the attribute
`"failed_maps"` (data.frame with the `map` and the `reason` of the maps
without any simulation, dropped from the output). The phylogenies have
the attributes `"seed"` (their seed, unless `seed = NULL`) and, with
several maps, `"map"`, to ensure reproducibility of each simulation.

## Details

### Stop conditions

The function **stops, before any grafting,** if the backbone is not
ultrametric, or if a sister clade is the whole tree (it has no stem
branch to graft on).

### Grafting process

Groups of missing taxa sharing the same sister clade are ordered by the
size of the sister clade, so **nested sister clades are processed
first**. For each group:

1.  A birth-death model is fitted with
    [`phytools::fit.bd()`](https://rdrr.io/pkg/phytools/man/fit.bd.html)
    on the sister clade, or on the smallest enclosing clade with at
    least `min_tips` tips. The model is fitted on the crown group of
    this clade, with the sampling fraction of this crown group (see
    Sampling fraction and unsampled species).

2.  The divergence time of the missing clade is drawn on the stem branch
    of the sister clade with
    [`TreeSim::corsim()`](https://rdrr.io/pkg/TreeSim/man/corsim.html):
    between the crown age of the sister clade and its stem age (or
    between the present and the stem age for a single-tip sister clade).

3.  A single missing taxon is grafted as a new tip at this divergence
    time. For several taxa, a missing clade is grown with
    [`TreeSim::sim.bd.taxa.age()`](https://rdrr.io/pkg/TreeSim/man/sim.bd.taxa.age.html)
    under the fitted model, conditioned on the number of taxa and on its
    stem age (the drawn divergence time), and grafted on the stem branch
    of the sister clade.

Unlike
[`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md),
this function has no `"uniform"` method: **a birth-death model is always
required**. It is needed to **generate the subtree of the missing
clade** under a diversification process. The model is fitted on the
sister clade (or its enclosing clade), and used to produce the subtree
to be grafted, conditioned on the number of taxa and its divergence
time. Additionally, the fitted model is used to draw the divergence time
of the missing clade on the stem branch of its sister clade.

### Enclosing clades

When a sister clade contains fewer tips than `min_tips`, an alternative
higher clade is selected to fit the birth-death model. This 'enclosing
clade' is identified by walking rootward the tree node by node, starting
from the sister clade, until a clade meets the `min_tips` condition.
This ensures that the diversification dynamics modeled are based on
enough tips to be relevant, while being closely associated with the
sister clade.

### Sampling fraction

The birth-death model is fitted on the crown group of the sister clade
or its enclosing clade, so the sampling fraction is:
`n_present / (n_present + n_pending + n_unsampled)`, with `n_present`
the tips of the fitted clade in the current tree (backbone tips and taxa
already grafted, so missing clades grafted in nested sister clades
processed earlier are not missing anymore), `n_pending` the taxa of
`missing_taxa_list` still to graft whose sister clade is strictly nested
within the fitted clade, and `n_unsampled` the unsampled species derived
from `total_richness_per_clades` (0 without it).

Taxa that are sister to the fitted clade itself are outside of its crown
group, so they are not counted: when the sister clade itself is fitted
(i.e. have at least `min_tips` tips) and nothing is pending nor
unsampled within it, the sampling fraction is 1. When an enclosing clade
is fitted, the taxa to graft on the sister clade (nested in the
enclosing clade) are counted. The divergence time of each missing clade
is still drawn with `TreeSim::corsim(missing = 1)`: unsampled species
only influence the fitted rates, they do not take a branch.

### Unsampled species

For a focal clade with a `total_richness` value provided:
`nb_unsampled = total_richness - nb_backbone - nb_missing`, with
`nb_missing` the number of taxa to graft as missing clades within the
focal clade (because their sister clade is nested within the focal
clade), therefore to be included in this focal clade during grafting.
For a clade with no `total_richness` value (`NA`, or enclosing a sister
clade without being listed), `nb_unsampled` is the sum of `nb_unsampled`
over the largest clades nested within it with a `total_richness` value
provided. If there is none, `nb_unsampled = 0`. The rationale is that
even if `nb_unsampled` cannot be derived from its `total_richness` (as
it is not provided for this clade), it has to at least account for the
`nb_unsampled` species that have being identified in its own nested
clades.

If `total_richness_per_clades = NULL`, `nb_unsampled = 0`: the clades
have no unsampled species (the same as providing only `NA` in
`$total_richness`).

**Each missing clade is grafted on the stem branch of its sister clade
in all the grafted phylogenies: only the divergence times vary among
simulations and maps.** The species present in each listed sister clade,
and the resulting sampling fractions for downstream analyses, are
therefore the same in all the grafted phylogenies. Use
[`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
to extract them a posteriori.

The function stops, before any grafting, if `total_richness_per_clades`
does not cover all the sister clades, if a `total_richness` provided is
lower than the number of backbone species plus missing species of the
clade, or if a clade has fewer unsampled species than the sum of the
clades with an informed `total_richness` nested within it.

### Species attributes: allowed states of a simmap

A `simmap` describes the states of a discrete attribute (e.g., a trait
or a range) along the branches of the backbone: each branch is cut in
pieces, each piece being in one state. Missing clades are grafted on the
stem branch of a clade of the backbone, and the function assumes that
the user considers this stem branch a valid place for them.
**Restricting the grafting to some states of the stem branch is
therefore only useful when the user has a strong hypothesis about a
change of state along the stem branch that should affect where the
missing clade can be grafted. For most uses, the `simmap` is not
needed.**

For example, in biogeography: a stem fossil is found in range A, and the
sister clade is in range A, while the missing clade is in range B. If
the shift from range B to range A on the stem branch of the sister clade
is estimated to have occurred before the age of the fossil, the missing
clade can be restricted to the part of the stem branch before the shift,
with `allowed_states = "B"`.

With a `simmap`, **a missing clade is only grafted on the pieces of the
stem branch of its sister clade that are in its allowed states**
(`missing_taxa_list$allowed_states` of its taxa if provided, otherwise
the global argument `allowed_states`). The divergence time of the
missing clade is drawn with
[`TreeSim::corsim()`](https://rdrr.io/pkg/TreeSim/man/corsim.html) as
without a `simmap`: if it falls in a part of the stem branch that is not
in an allowed state, it is drawn again, up to `max_age_redraws` times.

The `simmap` can also be provided so that the states are updated
alongside the phylogeny: **the whole missing clade (all its taxa and
branches) inherits the state of the point of the stem branch where it is
grafted.** The updated `simmap` is returned in `$simmap`, if
`return_simmap = TRUE`). The allowed states are not required: without
any allowed states, all the states are allowed, grafting is not
restricted by the states, but the simmap will be updated, assuming
grafted clades inherit the state of the point where they are grafted.

The function stops, before any grafting, if the taxa of a same missing
clade have different allowed states, or if the stem branch of a sister
clade has no piece in the allowed states (the clades are listed), and
warns about the missing clades with a small share of their stem branch
in allowed states (see `min_valid_fraction`).

**Cautions for interpretation: risks of circularity.**

Constraining the grafting by a trait that will then be analyzed (e.g.,
the association of the trait with diversification) may bias the results
towards the constraint. The grafted taxa are placed only where the trait
allows it, so they are in the allowed states by construction. This can
create or inflate an association between the trait and the
diversification dynamics, affect the rates of transitions between
states, or inflate the phylogenetic signal of the trait. It also reduces
the uncertainty about the state of the grafted clades since they always
inherit a unique state. It is recommended to acknowledge these potential
biases when interpreting such analyses.

To limit the risk of circularity, multiple approaches are possible:

- constrain the grafting by traits that are not the object of downstream
  analysis;

- compare the results with and without the constraint (or on the
  backbone only);

- keep the constraint as loose as the information allows (e.g., allowed
  states that include neighboring states);

- use several maps and simulations to propagate the uncertainty in
  ancestral state estimates.

### Biogeographic ranges as states

With `treat_states_as_ranges = TRUE`, the states of the `simmap` are
biogeographic ranges, coded as in `BioGeoBEARS` (one capital letter per
area, and the letters of a multi-area range in alphabetical order, e.g.,
`"AB"`), and the allowed states are the ranges a missing clade can have.

**An allowed range matches a state of the `simmap` if it is contained in
it**: a missing clade allowed in `"A"` can be grafted on the parts of
the stem branch in `"A"`, `"AB"` or `"ABC"`, whereas a missing clade
allowed in `"AB"` only can be grafted on the parts of the stem branch in
`"AB"` or `"ABC"`.

When multiple allowed states are provided, the inherited range is drawn
randomly from the allowed ranges enclosed within the range of the stem
branch piece where it is grafted. When grafting on an edge with a larger
range than the inherited range (e.g., `"A"` on `"AB"`), this results in
a subset sympatric speciation event such as (AB) -\> (AB), (A). In the
returned `simmap`, the branches of the grafted clades are in their
inherited range, whereas the stem branch it is grafted on continues in
its own range.

### Multiple simulations

With `nb_simulations > 1`, the grafting is simulated `nb_simulations`
times from the same backbone and the same `missing_taxa_list`. Each
simulation starts again from the backbone: the birth-death models are
fitted on the tree of the simulation (they depend on the taxa already
grafted), and the simulations do not share anything but the inputs. The
phylogenies are returned in a `"multiPhylo"` object (a single `"phylo"`
if `nb_simulations = 1`), and the memory used grows with the number of
simulations.

**The first simulation uses `seed`**, so it is the same phylogeny as a
run with `nb_simulations = 1` and the same `seed`. The seeds of the
other simulations are drawn from a stream started by `seed`, and
returned in the attribute `"seeds"`: a simulation can be reproduced
alone with its seed, and the first simulations of a longer run are the
same as the ones of a shorter run. With `seed = NULL`, the simulations
follow the random number generator of the session.

A failure that is specific to a simulation, because of the random draws
(no usable divergence age or time after `max_age_redraws` redraws),
skips this simulation with a warning: the name of a phylogeny keeps the
index of its simulation (e.g., `"sim_3"` is missing if the third
simulation failed). The function reports the number of phylogenies
obtained if it differs from `nb_simulations`, and stops if no simulation
succeeded. With `force_simulations = TRUE`, new seeds are tried until
`nb_simulations` phylogenies are obtained (up to `10 * nb_simulations`
attempts): the phylogenies are then named `"sim_1"` to `"sim_N"` without
gaps, and the seeds of the failed attempts are in the attribute
`"failed"`. Any other error (invalid inputs, no branch available for a
taxon, ...) concerns all the simulations and stops the whole run.

With `save_temp = TRUE`, each simulation saves its own temporary file
(see 'Temporary saves and resuming').

### Multiple stochastic maps

When multiple stochastic maps of the same trait are provided as a
`"multiSimmap"`, the grafting is run with each of its maps. The maps are
alternative hypotheses of the evolution of the same trait (or range),
used as alternative constraints. **All the maps and the simulations are
crossed**: each map is used for `nb_simulations` simulations, so `P`
maps give `P * nb_simulations` grafted phylogenies. This can lead to a
quick inflation of computation time, thus a warning is sent if the total
number of grafted phylogenies (`P * nb_simulations`) is more than
**1000**.

The output is a list with one element per stochastic map, named after
the maps (the names of the `"multiSimmap"` if it has some, otherwise
`"map_1"`, ..., `"map_P"`). With `nb_simulations = 1`, each element is a
grafted phylogeny. With several simulations, each element is a
`"multiPhylo"` (see 'Multiple simulations') with phylogenies named
`"<map>_sim_1"`, `"<map>_sim_2"`, ... and the attributes `"seeds"` and
`"failed"` of the map. Each phylogeny can be reproduced alone with its
map (as a single `simmap`) and its seed.

A map for which no simulation could be obtained is dropped from the
output with a warning, and listed in the attribute `"failed_maps"` of
the list (data.frame with the columns `map` and `reason`). The function
stops if no map gives any phylogeny.

With `save_temp = TRUE`, each map and each simulation saves its own
temporary file (see 'Temporary saves and resuming').

### Several traits

The maps of several traits can be used at once to constrain grafting.
`simmap` is then a list with one element per trait, each a `"simmap"` or
a `"multiSimmap"` (e.g.,
`simmap = list(color = color_simmap, elevation = elevation_simmap)`).
The elements can be named after the traits (with different names that
cannot contain `;`, `,`, `=` or `|`), otherwise the traits are named
`"trait_1"`, `"trait_2"`, ... All the maps must have the same tips,
topology and branch lengths as the backbone `phy`.

**A taxon can only be grafted on the parts of the branches that are in
one of its allowed states, for every trait.** The allowed states are
provided for each trait, in the order of the list of maps (or named
after the traits): the global argument `allowed_states` is a list with
one element per trait, e.g., `list(c("blue", "green"), NULL, "high")`
for three traits, and `missing_taxa_list$allowed_states` has such a list
for each taxon. A trait with `NULL` (or an empty vector, or `NA`) has no
taxon-level value, so the global `allowed_states` of the trait is used
for the taxon (a taxon-level value replaces the global one, trait by
trait, as with a single trait). A trait with no states, neither at the
taxon level nor globally, has all its states allowed (no restriction).
Other traits can still be restricted: `list(NULL, "high")` only
restricts the second trait.

With `return_simmap = TRUE`, `$simmap` is a list of updated simmaps, one
per trait, named after the traits. In each simmap, a grafted branch is
entirely in the state of each trait at the point where it was grafted
(each trait is inherited from the lineage it was grafted on).

With a list of `"multiSimmap"` provided as `simmap`, the maps of the
traits are paired: the grafting is run with the map 1 of every trait,
then the map 2 of every trait, and so on. **All the traits must
therefore have the same number of maps** to be associated in a given
grafting process.

One of the traits can be treated as biogeographic ranges combined with
any number of other traits. Use a logical vector in
`treat_states_as_ranges`to define which of the trait is to be treated as
biogeographic ranges for inheritance. The range logic is the one of
'Biogeographic ranges as states' for the range trait: its allowed values
for grafting are any ranges that encompass the allowed range. For
instance, a missing clade allowed in `"A;high"` can be grafted on stem
branch pieces in `"A;high"`, `"AB;high"` or `"ABC;high"`, whereas a
missing clade allowed in `"AB;high"` only can be grafted on stem branch
pieces in `"AB;high"` or `"ABC;high"`.

The other traits are matched exactly, and are inherited as they are from
the point where the clade is grafted: with `"AB"` and `"high"` at the
graft point and the allowed range `"A"`, the whole missing clade
inherits `"A"` and `"high"`. **Only one trait can be treated as ranges**
to prevent inflation of the number of combined allowed states to
consider.

### Parallel computation

The simulations based on the different stochastic maps and seeds are
independent, so they can be run in parallel with the package `future`
and `future.apply`. **Both packages must be installed independently from
phylograftR**.

The grafting function does not manage the workers: **you must set a plan
beforehand**, e.g., `future::plan(future::multisession, workers = 4)`,
and reset it afterwards with `future::plan(future::sequential)`. With
the default `parallel = NULL`, the function uses the plan if there is
one, `parallel = FALSE` ignores it. All the `P * nb_simulations`
simulations (with `P` stochastic maps per trait) are run as one set of
tasks, so the workers are used evenly.

The output is exactly the same as without parallelization. The messages,
warnings and printed output of the simulations are relayed to the
session when the tasks are collected (so they appear after the
simulations, not live), and with `verbose = TRUE` one line is printed
when each simulation is done.

Setting up the workers and sending the data to them has a time cost, so
**parallelization only pays off with large backbone phylogenies, or many
simulations or maps**. As a rule of thumb, use parallel computation when
the sequential run would last more than a minute or so, i.e., for large
backbone phylogenies (500 tips or more) and several grafted phylogenies
(`nb_maps * nb_simulations`), or for 50 or more grafted phylogenies on
small backbones.

### Temporary saves and resuming

With `save_temp = TRUE`, each grafted phylogeny (one per simulation and
per stochastic map) is saved in its own `.rds` file, built from
`temp_file_path` by adding the suffixes `"_map_<name>"` when a
`"multiSimmap"` is used, then `"_sim_<k>"` when `nb_simulations > 1`.
For instance, with `temp_file_path = "dir/temp.rds"`:

- a single phylogeny uses `"dir/temp.rds"`;

- several simulations use `"dir/temp_sim_1.rds"`,
  `"dir/temp_sim_2.rds"`, ...;

- several maps use `"dir/temp_map_A.rds"`, `"dir/temp_map_B.rds"`, ...;

- several simulations across several maps use both suffixes as
  `"dir/temp_map_A_sim_1.rds"`;

- the prefix `"map_"` is not repeated if the name of the map already
  starts with it: maps named `"map_1"`, `"map_2"`, ... (default names of
  a `"multiSimmap"` without names) use `"dir/temp_map_1_sim_1.rds"`,
  `"dir/temp_map_2_sim_1.rds"`, ...

With `force_simulations = TRUE`, `<k>` is the number of the attempt, not
necessarily the ID to which the grafted phylogeny will be attributed to
in the final output. This is required to identify the appropriate seed
used to produce each attempt; attempts that failed leave a file too for
tractability.

Temporary files can be read with
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) and used to plot the
grafted phylogeny in its current state. The object is a regular
`"phylo"` with the attribute `"phylograftR_temp"` that stores what is
needed to resume a run. Once the phylogeny is fully grafted, the file is
replaced by the complete output of the simulation (with `$simmap`,
`$grafted_tips`, etc.).

With `clean_temp_files = TRUE` (default), all the temporary files of the
run are deleted once the function returns without crashing. **A run that
crashes or is interrupted never deletes its files**, so it can be
resumed with `resume = TRUE`. Use `clean_temp_files = FALSE` to keep the
files after a successful run.

With `resume = TRUE`, the function must be called again with **exactly
the same arguments** as the interrupted run. For each phylogeny, its
temporary file is searched: a complete phylogeny is reused as it is, an
unfinished one is continued from its last save, and a phylogeny without
file is grafted from scratch.

**Resuming is not equivalent to grafting the taxa left on the temporary
phylogeny used as a new backbone**, with a `missing_taxa_list`
restricted to these taxa: the taxa already grafted would be treated as
backbone species, so the sampling fractions of the birth-death models,
the stem of the target clades, the identification of the grafted tips
and the updated `simmap` would be wrong, and a `simmap` of the initial
backbone would not match the temporary phylogeny. Use `resume = TRUE`
instead.

### Successive graftings

To graft the missing taxa and the missing clades on the same backbone,
use
[`graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md),
which chains both steps with a single set of inputs and a single seed,
including with several simulations and/or stochastic maps.

The output of this function can also be used as the `phy` of another
grafting (e.g., to chain the steps by hand for a single grafted
phylogeny).

**The grafted tips and branches accumulate**: if `phy` has the
`grafted_tips` element of a previous grafting, its flagged tips are kept
as grafted (they are still treated as backbone species by the grafting
itself, e.g., to count the species used for the sampling fractions), and
the new taxa are added to them. `grafted_tips` and `grafted_edges` of
the output thus identify **all** the grafted taxa, not only those of the
last step. A warning is sent if `grafted_tips` does not match the tips
of `phy` (e.g., after
[`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html)): use
[`identify_grafted()`](https://maeldore.github.io/phylograftR/reference/identify_grafted.md)
to recompute them first.

If a `simmap` is used, provide the updated `simmap` of the previous step
(`phy$simmap`, obtained with `return_simmap = TRUE`) to keep the same
trait history.

## See also

Other functions from phylograftR:

- [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  to build a `missing_taxa_list` from a data.frame loaded in R, with
  clades defined by tips or node IDs.

- [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
  to build a `missing_taxa_list` from taxonomic tables, with clades
  defined by taxonomic names.

- [`phylograftR::build_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/build_total_richness_per_clades.md)
  and
  [`phylograftR::import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md)
  to prepare `total_richness_per_clades` for the unsampled species.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades.

- [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)
  to graft missing taxa and missing clades in one run, with a single set
  of inputs.

- [`phylograftR::check_grafting_plan()`](https://maeldore.github.io/phylograftR/reference/check_grafting_plan.md)
  to check a priori the feasibility of a grafting plan, before grafting.

- [`phylograftR::check_grafting_constraints()`](https://maeldore.github.io/phylograftR/reference/check_grafting_constraints.md)
  to check a posteriori that grafted taxa respect their constraints.

- [`phylograftR::extract_grafting_info()`](https://maeldore.github.io/phylograftR/reference/extract_grafting_info.md)
  to extract a posteriori information on the location of grafted
  lineages.

- [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  to extract a posteriori the richness and the sampling fraction of the
  clades in each grafted phylogeny.

- [`phylograftR::extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md)
  to extract a posteriori information on the birth-death models fitted
  during grafting (with `return_BD_model_summary = TRUE`).

- [`phylograftR::plot_grafted()`](https://maeldore.github.io/phylograftR/reference/plot_grafted.md)
  to plot grafted phylogenies.

## Author

Maël Doré

## Examples

``` r
# --- Prepare data --- #

# Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
set.seed(1234)
phy <- ape::rcoal(n = 30)
phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50

# List the missing taxa to graft
missing_clades_list <- list(
  Taxon_name = c("M1", "M2", "M3", "M4"),
  # Sister_clade is named after Taxon_name (same names, same order)
  # M1 and M2 share the same sister clade: they are grafted together as one missing clade
  Sister_clade = list(M1 = c("t8", "t22"), M2 = c("t8", "t22"), M3 = c("t26", "t18"), M4 = "t5"))

# --- Perform grafting --- #

# Graft the missing clades
grafted_phy <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
                                    verbose = FALSE, seed = 1234)
plot(grafted_phy, cex = 0.7)


# Highlight the grafted tips and branches in red
plot_grafted(grafted_phy, cex = 0.7)



# --- Produce several grafted phylogenies --- #

# The grafting is stochastic: several simulations can be run at once from the same backbone and list
# (a "multiPhylo" named sim_1, sim_2, ...; sim_1 is the same as a single run with the same seed)
grafted_phys <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
                                      nb_simulations = 3, verbose = FALSE, seed = 1234)
names(grafted_phys)
#> [1] "sim_1" "sim_2" "sim_3"
# Seeds of the simulations, to reproduce each phylogeny alone
attr(grafted_phys, "seeds")
#>     sim_1     sim_2     sim_3 
#>      1234 488349520 469344150 

# Plot phylos
old_par <- par()["mfrow"]
par(mfrow = c(1, 3))
plot_grafted(grafted_phys, index = 1, cex = 0.7, main = "Sim #1") # sim_1
plot_grafted(grafted_phys, index = 2, cex = 0.7, main = "Sim #2") # sim_2
plot_grafted(grafted_phys, index = 3, cex = 0.7, main = "Sim #3") # sim_3

par(old_par)


# --- Account for unsampled species --- #

# Account for unsampled species with the total richness of the sister clades.
# All the unique sister clades must be listed (NA = unknown richness)
total_richness_per_clades <- list(
  Sister_clade = list(c("t8", "t22"), c("t26", "t18"), "t5"),
  total_richness = c(20, NA, 3))
grafted_phy_richness <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
                                             total_richness_per_clades = total_richness_per_clades,
                                             verbose = FALSE, seed = 1234)

# Inspect the resulting clade richness and sampling fractions
extract_sampling_info(phy_grafted = grafted_phy_richness,
                      total_richness_per_clades = total_richness_per_clades)
#>    map simulation seed Clade_tips total_richness total_richness_provided
#> 1 <NA>          1 1234    t8, t22             20                    TRUE
#> 2 <NA>          1 1234   t26, t18             17                   FALSE
#> 3 <NA>          1 1234         t5              3                    TRUE
#>   nb_backbone_species nb_grafted_species nb_unsampled_species sampling_fraction
#> 1                  13                  1                    6         0.7000000
#> 2                  17                  0                    0         1.0000000
#> 3                   1                  0                    2         0.3333333


# --- Apply attribute-based exclusions --- #

# Restrict the grafting to some states of a stochastic character map (simmap) of the backbone
# The toy simmap has three states, "blue", "green" and "red"
data(toy_backbone)
data(toy_simmap)

# List the missing taxa to graft
missing_clades_list <- list(
   Taxon_name = c("M1", "M2", "M3", "M4", "M5"),
   Sister_clade = list(M1 = c("t26", "t18"), M2 = c("t26", "t18"), M3 = c("t26", "t18"),
                       M4 = "t14", M5 = "t14"))

# The simmap shows that the targeted stem branches recorded state shifts.
plot(toy_simmap, col = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))


# We can use 'allowed_states' to limit grafting on the "green" portions of the stem branches
   # The stem branch from the "bottom" clade shifts from "green" to "red"
   # The terminal (stem) branch of "t14" shift from "green" to "red" too
grafted_phy_states <- graft_missing_clades(phy = toy_backbone,
                                           missing_taxa_list = missing_clades_list,
                                           simmap = toy_simmap, allowed_states = "green",
                                           verbose = FALSE, seed = 1234)

# Plot output phylogeny with colored grafts
plot_grafted(grafted_phy_states, use_simmap = FALSE, cex = 0.7)


# The updated simmap is in the output: the new grafted clades branch from the "green" segments
# of the stem branches of their sister clades, and each newly grafted clades are entirely "green"
# Grafted branches are dashed lines, and grafted tips are labelled in red.
plot_grafted(grafted_phy_states,
             simmap_colors = c("blue" = "dodgerblue", "green" = "limegreen", "red" = "firebrick"))



# --- Use several stochastic maps (alternative hypotheses of trait evolution) as constraints --- #

# A 'multiSimmap' gives several stochastic maps of the same trait:
# the grafting is run with each of them.
# Here, 3 maps are simulated from the states of the tips of the toy simmap.
tip_states <- phytools::getStates(toy_simmap, type = "tips")

# (the random number generator of the session is restored afterwards)
had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
set.seed(118)
toy_multiSimmap <- phytools::make.simmap(toy_backbone, tip_states,
                                         model = "ER", nsim = 3, message = FALSE)
if (had_seed)
{ assign(".Random.seed", old_seed, envir = globalenv())
} else { rm(".Random.seed", envir = globalenv()) }

# Graft using constraints from all 3 simmaps X 2 simulations per simmap
grafted_phys_maps <- graft_missing_clades(phy = toy_backbone,
                                         missing_taxa_list = missing_clades_list,
                                         simmap = toy_multiSimmap, allowed_states = "green",
                                         nb_simulations = 2, verbose = FALSE, seed = 1234)

# Inspect output:
# A list with one element per simmap, each a multiPhylo with the simulations of the map
names(grafted_phys_maps)
#> [1] "map_1" "map_2" "map_3"
names(grafted_phys_maps$map_2)
#> [1] "map_2_sim_1" "map_2_sim_2"
# The map and the seed of a phylogeny are included as attributes for reproducibility
attr(grafted_phys_maps$map_2[[1]], "map")
#> [1] "map_2"
attr(grafted_phys_maps$map_2[[1]], "seed")
#> [1] 1234


# --- Constraints from several traits at once --- #

# Load the simmap for the second trait: 'toy_simmap_elevation' has two states: "high" and "low".
data(toy_simmap_elevation)

# Provide the maps of the traits as a list
# The global 'allowed_states' has one element per trait (here, "green" color and "low" elevation).
grafted_phy_traits <- graft_missing_clades(phy = toy_backbone,
                                           missing_taxa_list = missing_clades_list,
                                           simmap = list(color = toy_simmap,
                                                         elevation = toy_simmap_elevation),
                                           allowed_states = list("green", "low"),
                                           verbose = FALSE, seed = 1234)

# '$simmap' is a list of updated simmaps, one per trait
names(grafted_phy_traits$simmap)
#> [1] "color"     "elevation"

# Newly grafted clades can only be grafted on and harbor the allowed state in each trait
plot_grafted(grafted_phy_traits, trait = 1) # Only 'green'

plot_grafted(grafted_phy_traits, trait = 2) # Only 'low'



# --- Constraints from biogeographic ranges and another trait at once --- #

# One of the traits can be treated as biogeographic ranges:
# the allowed ranges do not have to match exactly the range of target stem branch pieces,
# any allowed range comprised in a target range is valid (e.g., "A" can be grafted on "AB").
# Here, the first trait holds ranges.

# Load the range simmap
data(toy_simmap_ranges)

# The clades are grafted on stem branch pieces in ranges
# containing "A" or "B" and in elevation "low"
grafted_phy_ranges_traits <- graft_missing_clades(phy = toy_backbone,
    missing_taxa_list = missing_clades_list,
    simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
    allowed_states = list(c("A", "B"), "low"),
    treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
    verbose = FALSE, seed = 1234)

# '$simmap' is a list of updated simmaps, one per trait: ranges, and elevation
names(grafted_phy_ranges_traits$simmap)
#> [1] "biogeo"    "elevation"

# Newly grafted clades can only harbor the allowed state/range in each trait
# They can be grafted on range that encompasses their allowed range (here, "AB")

# Range: Grafted on "A", "B", or "AB", but inherit only "A" or "B".
plot_grafted(grafted_phy_traits, trait = 1)

# Elevation: Only 'low'
plot_grafted(grafted_phy_traits, trait = 2)


```
