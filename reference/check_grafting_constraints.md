# Check that grafted taxa and missing clades respect their grafting constraints

This function checks, **a posteriori**, that each taxon grafted on a
phylogeny sits in a valid position. For missing taxa grafted with
[`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md):
within its target clade (or on its stem branch if allowed), and outside
its excluded subclades. For missing clades grafted with
[`graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md):
each missing clade (the taxa sharing the same sister clade) must form
one clade, grafted on the stem branch of its sister clade.

The type of grafting is detected from `missing_taxa_list`: a
`$Target_clade` item for taxa, a `$Sister_clade` item for missing clades
(see 'Missing clades' in Details section). With a `simmap`, the validity
of the states of the taxa or of the missing clades are checked too.

The function can also **check that the backbone phylogeny is preserved**
by the grafting once the grafted taxa are pruned.

## Usage

``` r
check_grafting_constraints(
  phy_grafted,
  phy_backbone,
  missing_taxa_list,
  include_target_clade_stem = TRUE,
  use_excluded_subclades = TRUE,
  exclude_subclades_stem = TRUE,
  allowed_states = NULL,
  treat_states_as_ranges = FALSE,
  check_backbone = TRUE,
  tolerance = 1e-06,
  stop_on_failure = FALSE,
  verbose = TRUE
)
```

## Arguments

- phy_grafted:

  An object of class `"phylo"`. The phylogeny after grafting. It can
  also be a set of grafted phylogenies (class `"multiPhylo"`, from
  `nb_simulations > 1`), or a list with one element per stochastic map
  (a `"phylo"` or a `"multiPhylo"`, from the grafting functions used
  with a `"multiSimmap"` to account for multiple alternative trait
  evolutionary histories). See 'Several grafted phylogenies' in Details
  section.

- phy_backbone:

  An object of class `"phylo"`. The backbone phylogeny used for grafting
  (i.e., before grafting).

- missing_taxa_list:

  A list with the information on the grafted taxa. Typically the same
  list as used for grafting. Items: `$Taxon_name`, `$Target_clade` (taxa
  grafted with
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md))
  or `$Sister_clade` (missing clades grafted with
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)),
  and optionally `$include_target_clade_stem`, `$Excluded_subclades` and
  `$exclude_subclades_stem` (taxa only), and `$allowed_states`. Exactly
  one of `$Target_clade` and `$Sister_clade` must be provided: it
  defines the type of the check. Every item provided must be named after
  `$Taxon_name` (same names, same order), as in the grafting functions,
  see these functions for details.

  It can also be **the list of both types of lineages to graft**,
  `list(taxa = ..., clades = ...)`, as returned by
  [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  when both types are provided, to check the output of
  [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md).
  Each part is then checked on the same phylogenies and the output is a
  named list (see 'Value'). **In a missing clade, tips grafted as
  missing taxa, during the other grafting step, are tolerated**: when
  the missing clades are grafted first, taxa can be grafted within them
  afterwards.

- include_target_clade_stem:

  Logical. Same as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md):
  whether grafting on the stem branch of target clades is accepted. Taxa
  with a value in `missing_taxa_list$include_target_clade_stem` use it
  instead (`NA` or no item: this global value is used). Ignored for
  missing clades. Default is `TRUE`.

- use_excluded_subclades:

  Logical. Same as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md):
  whether to check the excluded subclades. If no `$Excluded_subclades`
  item is found in `missing_taxa_list`, a warning is issued and excluded
  subclades are not checked. Ignored for missing clades. Default is
  `TRUE`.

- exclude_subclades_stem:

  Logical. Same as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md):
  whether the stem of the excluded subclades is also excluded (a taxon
  sister to an excluded subclade is then a violation). Excluded
  subclades with a value in `missing_taxa_list$exclude_subclades_stem`
  use it instead (`NA` or no item: this global value is used). Ignored
  for missing clades. Default is `TRUE`.

- allowed_states:

  Character vector. Same as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md):
  the states in which grafting is allowed. This is the default for all
  taxa: taxa with a value in `missing_taxa_list$allowed_states` use this
  specific value instead. Only used if `phy_grafted` has a `$simmap`
  item (the updated simmap returned by
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  with a `simmap`). If the `$simmap` item is a list of the simmaps of
  several traits (grafting constrained with several traits),
  `allowed_states` is a list with one element per trait (in the order of
  the traits, or named after the traits; `NULL` = all the states of the
  trait, or any range for the trait treated as ranges), and the
  taxon-level values of `missing_taxa_list$allowed_states` are given per
  trait as in the grafting functions (see 'Several traits' in Details).
  The function stops if some states are not in the `simmap`. Default is
  `NULL` (all states allowed for the taxa without a taxon-level value).

- treat_states_as_ranges:

  Logical. Same as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md):
  whether the states of the `simmap` and the allowed states are treated
  as biogeographic ranges: an allowed range enables grafting on branch
  pieces in any range that contains it: `"A"` can be grafted on `"A"`,
  `"AB"` and `"ABC"`, whereas `"AB"` can be grafted on `"AB"` and
  `"ABC"` but not `"A"`. Only used if `phy_grafted` has a `$simmap`
  item, **otherwise the function stops if `TRUE`**. With a `$simmap`
  item that is a list of the simmaps of several traits, takes a logical
  vector with one value per trait, with at most one `TRUE` to identify
  the trait treated as ranges (positional, or named after the traits in
  any order, as in the grafting functions; a single `FALSE` is accepted,
  a single `TRUE` stops). States / ranges are coded as in BioGeoBEARS
  (see the Details section "Biogeographic ranges as states"). The
  function stops if the states of the `simmap` or the allowed states are
  not correctly coded, or if an area of the allowed states is not in the
  `simmap`. Default is `FALSE` (states are matched exactly).

- check_backbone:

  Logical. Whether to check that the topology and branch lengths of the
  backbone are preserved once the grafted taxa are pruned. Default is
  `TRUE`.

- tolerance:

  Numeric. Tolerance for the maximum deviation in backbone branch
  lengths, relative to the depth (root age) of the backbone. Default is
  `1e-6`.

- stop_on_failure:

  Logical. Whether to stop with an error if any constraint is violated.
  Default is `FALSE`.

- verbose:

  Logical. Whether to print a summary of the checks. Default is `TRUE`.

## Value

The structure of the output depends on the type of elements grafted, and
their associated grafting function used in the process.

**Missing taxa.** With a `$Target_clade` list (for missing taxa), the
function assumed that
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
was used for the grafting. It returns a data.frame with one row per
taxon in `missing_taxa_list` (in the original order) for each phylogeny
checked, and the following columns:

- `$map`: Name of the stochastic map used to constrain grafting based on
  species attributes. `NA` if none was used.

- `$simulation`: Integer. Index of the simulation (the index in the name
  `"sim_k"` of the phylogeny; 1 for a single phylogeny).

- `$seed`: Integer. Seed of the simulation (attribute `"seed"` of the
  phylogeny, `NA` if unknown).

- `$Taxon_name`: Name of the taxon.

- `$in_tree`: Logical. Whether the taxon is found in the grafted
  phylogeny.

- `$include_target_clade_stem`: Logical. Whether the stem of the target
  clade was accepted for this taxon
  (`missing_taxa_list$include_target_clade_stem` of the taxon if
  provided, otherwise the global `include_target_clade_stem` argument).

- `$target_position`: Character string. Position of the taxon relative
  to its target clade: `"crown"` (within the clade), `"stem"` (sister to
  the clade, or to other missing taxa grafted on the stem),
  `"tip_branch"` (on the branch of the tip, for single-tip target
  clades), or `"outside"` (neither within the clade nor on its stem = an
  error).

- `$target_ok`: Logical. Whether the position respects the target clade
  (`"crown"` or, if the stem is accepted for the taxon, `"stem"`;
  `"tip_branch"` for single-tip clades).

- `$excluded_ok`: Logical. Whether the taxon is outside all its excluded
  subclades (and not on their stem, for the subclades with an excluded
  stem).

- `$n_excluded_violated`: Integer. Number of excluded subclades
  violated.

- `$excluded_violated`: Character string. Names of the excluded
  subclades violated (or their indices if the subclades are not named),
  separated by commas (`NA` if none).

- `$states_ok`: Logical. Only checked if `phy_grafted` has a `$simmap`
  item: whether the taxon was grafted in one of its allowed states
  (`missing_taxa_list$allowed_states` of the taxon if provided,
  otherwise the global `allowed_states` argument; all the states if
  neither is provided). With `treat_states_as_ranges = TRUE`, the state
  of the terminal branch of the taxon (the range it inherited) must be
  one of its allowed ranges (any range if the taxon has no allowed
  states). `NA` if `phy_grafted` has no `$simmap`. With the maps of
  several traits, `TRUE` if the taxon was grafted in an allowed state
  for all the traits.

- `$states_ok_<trait>`: Logical. Only with the maps of several traits
  (one column per trait, named after the traits): whether the taxon was
  grafted in one of its allowed states for this trait. Used to find
  which trait is violated when `$states_ok` is `FALSE`.

- `$ranges_ok`: Logical. Only checked with
  `treat_states_as_ranges = TRUE` (or a range trait with several traits)
  and a `$simmap` item: whether the range inherited by the taxon is
  contained in the state of the lineage at the graft point (see the
  Details section "Biogeographic ranges as states"). `NA` if not
  checked, or if the taxon is grafted at the root (no parent lineage).
  This differs from `$states_ok` that checks only the validity of the
  terminal tip state. Here, the validity of the inherited range is
  assessed against the range of the backbone lineage it was grafted on.

- **`$ok`**: Logical. Whether the taxon respects its target clade, its
  excluded subclades, and its allowed states and ranges (if checked).

**Missing clades.** With a `$Sister_clade` list (for missing clades),
the function assumed that
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
was used for the grafting. It returns a data.frame with one row per
missing clade, i.e., per group of taxa sharing the same sister clade,
numbered in the order of their first appearance in `missing_taxa_list`,
for each phylogeny checked.

The columns `$map`, `$simulation` and `$seed` are the same as above,
followed by:

- `$missing_clade`: ID or label of the missing clade: the label of its
  sister clade if `missing_taxa_list` has a `$Sister_clade_name` item
  (character string), otherwise the ID of the clade (integer).

- `$Taxon_names`: Character string. The taxa of the missing clade,
  separated by commas.

- `$n_taxa`: Integer. Number of taxa of the missing clade.

- `$Sister_clade_tips`: Character string. The tips that define the
  sister clade, separated by commas.

- `$in_tree`: Logical. Whether all the taxa of the missing clade are
  found in the grafted phylogeny. `$n_in_tree` is the number of taxa
  found.

- `$monophyletic`: Logical. Whether the taxa form one clade of the
  grafted phylogeny with no other tip (backbone tip or taxon of another
  missing clade).

- `$sister_position`: Character string. Position of the missing clade
  relative to its sister clade: `"stem"` (sister to the clade, or to
  other missing taxa grafted on its stem, with no other backbone tip in
  between), `"tip_branch"` (on the branch of the tip, for a single-tip
  sister clade), `"crown"` (nested within the sister clade), or
  `"outside"`.

- `$sister_ok`: Logical. Whether the missing clade is on the stem of its
  sister clade (`"stem"` or `"tip_branch"`).

- `$states_ok`: Logical. Only checked if `phy_grafted` has a `$simmap`
  item: whether the whole missing clade (its stem branch and all its
  descendant branches) is in one single state, and this state is allowed
  (`missing_taxa_list$allowed_states` if provided, otherwise the global
  `allowed_states` argument; all the states if neither is provided).
  `NA` if `phy_grafted` has no `$simmap`. With the maps of several
  traits, `TRUE` if this is the case for all the traits, and
  `$states_ok_<trait>` gives the result of each trait. With
  `treat_states_as_ranges = TRUE`, the single state is the range
  inherited by the missing clade, which must be one of its allowed
  ranges (any range if none).

- `$ranges_ok`: Logical. Only checked with
  `treat_states_as_ranges = TRUE` (or a range trait): whether the range
  inherited by the missing clade is contained in the state of the
  lineage at the graft point.

- **`$ok`**: Logical. Whether the missing clade is in the tree,
  monophyletic, on the stem of its sister clade, and in an allowed state
  and range (if checked).

**With a `list(taxa, clades)` covering both types of missing lineages**,
the output is a named list with the elements `taxa` and `clades`\*\*,
each element contains the output of its type (as described above), with
the overarching attribute `"all_ok"` (`TRUE` if all the constraints are
respected for both types of grafted lineages).

**Multiple maps.** With a list of maps as `phy_grafted`, the function
returns a list of such data.frames (one per map, named after the maps,
with the attribute `"all_ok"` for all the maps), since the maps have
their own constraints: in each data.frame, the rows of all the
simulations of the map are merged.

The data.frame has the following attributes (for several phylogenies:
`"topology_preserved"` and `"backbone_preserved"` are `TRUE` if they are
`TRUE` for all the phylogenies, the deviations are the maximum over the
phylogenies, and `"all_ok"` is `TRUE` if it is `TRUE` for all of them):

- `"topology_preserved"`: Logical. Whether the backbone topology is
  preserved once the grafted taxa are pruned (`NA` if
  `check_backbone = FALSE`).

- `"max_branch_length_deviation"`: Numeric. Maximum absolute deviation
  in backbone branch lengths (`NA` if the topology is not preserved).

- `"max_branch_length_deviation_relative"`: Numeric. Same, relative to
  the depth of the backbone.

- `"backbone_preserved"`: Logical. Whether both topology and branch
  lengths are preserved, within `tolerance` (`NA` if
  `check_backbone = FALSE`).

- **`"all_ok"`**: Logical. Whether all taxa respect their constraints
  and the backbone is preserved (if checked).

- `"trait_names"`: Character vector. Names of the traits checked (only
  with the maps of several traits, otherwise `NULL`).

- `"range_trait"`: Character string. Name of the trait treated as ranges
  (only with the maps of several traits, otherwise `NULL`).

- `"phylogenies"`: data.frame with one row per phylogeny checked
  (columns `map`, `simulation`, `seed`, `n_failed` (number of taxa, or
  of missing clades, that failed), `backbone_preserved` and `all_ok`).

Attributes are accessed with
`attr(x = output, which = "attribute_name")`.

## Details

### Missing taxa

**A taxon respects its target clade if it is found within the clade, or
on its stem branch if the stem is accepted for this taxon**: global
`include_target_clade_stem = TRUE`, or
`missing_taxa_list$include_target_clade_stem` of this taxon (the taxon
is then sister to the target clade, or to other missing taxa sister to
the target clade, with no other backbone tip in between). For a
single-tip target clade, the taxon must be on the branch of the tip,
with no other backbone tip in between.

The 'stem' position is not restricted to a taxon being directly sister
to the target clade: **taxa grafted on the stem can receive further
grafts**, so several missing taxa can form together a group sister to
the target clade (i.e, a **missing-only subtree** grafted on the the
stem). Only backbone tips are considered when checking the position of a
taxon (missing taxa are ignored), so such configurations are valid.

If `phy_grafted` has a `$simmap` item (the updated `simmap` returned by
the grafting functions when a `simmap` is provided), the states are also
checked: a grafted taxon is in the state of the point where it was
grafted, along its whole terminal branch, so the states of its terminal
branch in the `simmap` must be allowed. The `simmap` item is used as it
is: it is not compared to the backbone.

**A taxon respects its excluded subclades if it is not nested within any
of them.** If the stem of an excluded subclade is excluded
(`exclude_subclades_stem = TRUE`, the default), the taxon must also not
be on its stem, i.e., sister to the subclade (or to other missing taxa
sister to the subclade), with no other backbone tip in between. For a
single-tip excluded subclade, the taxon must not be on the branch of the
tip. As in
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md),
if the stem of the subclade is not excluded, the taxon can be sister to
the subclade, not within it.

Branches of the backbone are matched between the pruned grafted
phylogeny and the backbone by the set of tips they subtend. Branch
lengths are considered preserved if the maximum deviation, relative to
the depth of the backbone, is lower than `tolerance`.

### Missing clades

[`graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
grafts the taxa sharing the same sister clade together, as **one missing
clade**, on the stem branch of their sister clade. With a
`$Sister_clade` list, this function checks each missing clade (the taxa
whose sister clades have the same backbone tips, whatever the pair of
tips used to define them) for:

- its presence (`$in_tree`): all its taxa must be in the grafted
  phylogeny ;

- its monophyly (`$monophyletic`): its taxa must form one clade of the
  grafted phylogeny, free of any other tip ;

- its position (`$sister_position` and `$sister_ok`): the missing clade
  must be sister to its sister clade, or to other missing taxa grafted
  on the same stem, with no other backbone tip in between (`"stem"`), or
  on the terminal branch of a single-tip sister clade (`"tip_branch"`).
  A missing clade nested within its sister clade (`"crown"`) is a
  violation, unlike a taxon grafted by
  [`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  within its target clade;

- with a `$simmap` item, its states (`$states_ok`): the whole missing
  clade is in the state of the point where it was grafted, so its stem
  branch and all its descendant branches must have one single state,
  which must be allowed. The taxa of a missing clade are grafted
  together, so they must share their allowed states (the function stops
  otherwise, as the grafting function does). With
  `treat_states_as_ranges = TRUE`, the whole missing clade inherits one
  range, which must be contained in the state of the lineage at the
  graft point (see 'Biogeographic ranges as states').

\*e arguments `include_target_clade_stem`, `use_excluded_subclades` and
`exclude_subclades_stem` only apply to taxa and are ignored. The
backbone is checked as for taxa. Several phylogenies, maps and traits
are handled as described in the corresponding sections. With
`verbose = TRUE`, the summary counts missing clades instead of taxa.

### Several traits

If the grafting used several traits as constraints, the grafting
functions return a `$simmap` item that is a list of simmaps, one per
trait, named after the traits, (or `"trait_1"`, `"trait_2"`, ... if the
input was not named). Then this function checks the validity of states
trait by trait: the states of the terminal branch of the taxon in the
simmap of each trait must be allowed for this trait. The allowed states
of each taxon and trait are resolved as in the grafting functions: the
value of `missing_taxa_list$allowed_states` for this trait if provided,
otherwise the global `allowed_states` for this trait, otherwise all the
states of the trait (see
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).
`$states_ok_<trait>` gives the result for each trait, and **`$states_ok`
is `TRUE` only if all traits are ok**. With `verbose = TRUE`, the number
of violations is also given for each trait.

With a trait treated as ranges (`treat_states_as_ranges`, one `TRUE` in
the vector), the checks of 'Biogeographic ranges as states' are
performed on this trait: `$states_ok_<trait>` is `TRUE` if the range
inherited is one of the allowed ranges (any range if the taxon has no
allowed ranges), and `$ranges_ok` is `TRUE` if the range inherited is
contained in the state of the lineage of this trait at the graft point.
The other traits are checked as exact matches. The attribute
`"range_trait"` gives the name of the range trait.

### Several grafted phylogenies

The phylogenies are checked one by one against the same `phy_backbone`
and `missing_taxa_list`. When `allowed_states` or
`treat_states_as_ranges` are used, each phylogeny is checked with its
own `$simmap` item (the maps of a `"multiSimmap"` can differ). The
warnings about the inputs are sent once. With several phylogenies,
`verbose = TRUE` prints one summary for all of them (number of
phylogenies respecting all the constraints, number of violations of each
type, and the first phylogenies with violations), and
`stop_on_failure = TRUE` stops if any phylogeny violates a constraint.

### Biogeographic ranges as states

With `treat_states_as_ranges = TRUE`, the states are biogeographic
ranges (one capital letter per area, e.g., `"A"`, `"AB"`, `"BC"`), and a
range is a set of areas. A grafted lineage inherits one of its allowed
ranges, which can be a subset of the state of the lineage it is grafted
on. Example: the range `"A"` can be inherited from a lineage in the
state `"AB"`. This then describes a subset sympatric speciation event:
`(AB) -> (AB), (A)`.

Two checks are then performed:

- `$states_ok`: the state of the terminal branch of the taxon is exactly
  one of its allowed ranges (the range inherited).

- `$ranges_ok`: the range inherited is contained in the state of the
  lineage at the graft point (e.g., `"A"` is contained in `"AB"`, but
  `"BC"` is not). The check is skipped (`NA`) for taxa grafted at the
  root.

## See also

Other functions from phylograftR:

- [`phylograftR::check_grafting_plan()`](https://maeldore.github.io/phylograftR/reference/check_grafting_plan.md)
  to check a priori the feasibility of a grafting plan, before grafting.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades.

- [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft groups of missing taxa as sister clades of backbone clades.

- [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)
  to graft missing taxa and missing clades in one run, with a single set
  of inputs.

- [`phylograftR::extract_grafting_info()`](https://maeldore.github.io/phylograftR/reference/extract_grafting_info.md)
  to extract a posteriori information on the location of grafted
  lineages.

- [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  to extract a posteriori the richness and the sampling fraction of the
  clades in each grafted phylogeny.

- [`phylograftR::extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md)
  to extract a posteriori information on the birth-death models fitted
  during grafting.

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
missing_taxa_list <- list(
  Taxon_name = c("X1", "X2", "X3", "X4"),
  # Items are named after Taxon_name (same names, same order)
  Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")),
  Excluded_subclades = list(X1 = list(c("t10", "t19")), X2 = NULL, X3 = NULL, X4 = NULL))

# --- Perform grafting of missing taxa--- #

# Graft the missing taxa with graft_missing_taxa()
grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
                                  verbose = FALSE, seed = 1234)

# --- Check one grafted phylogenies --- #

# Check that the grafted taxa respect their constraints, and that the backbone is preserved
check_results <- check_grafting_constraints(phy_grafted = grafted_phy, phy_backbone = phy,
                                            missing_taxa_list = missing_taxa_list)
#> 4 taxa checked
#>   - Missing from the grafted tree: 0
#>   - Target clade violated:          0
#>   - Excluded subclade violated:     0
#>   - Backbone: topology preserved, branch lengths shifted by max 4.82e-09 (9.65e-11 of backbone depth; tolerance = 1e-06)
#>   => All constraints respected
# Inspect results
check_results
#>     map simulation seed Taxon_name in_tree include_target_clade_stem
#> X1 <NA>          1 1234         X1    TRUE                      TRUE
#> X2 <NA>          1 1234         X2    TRUE                      TRUE
#> X3 <NA>          1 1234         X3    TRUE                      TRUE
#> X4 <NA>          1 1234         X4    TRUE                      TRUE
#>    target_position target_ok excluded_ok n_excluded_violated excluded_violated
#> X1           crown      TRUE        TRUE                   0                NA
#> X2           crown      TRUE        TRUE                   0                NA
#> X3      tip_branch      TRUE        TRUE                   0                NA
#> X4           crown      TRUE        TRUE                   0                NA
#>    states_ok ranges_ok   ok
#> X1        NA        NA TRUE
#> X2        NA        NA TRUE
#> X3        NA        NA TRUE
#> X4        NA        NA TRUE
# Check overall results
attr(check_results, "all_ok")
#> [1] TRUE

# --- Check several grafted phylogenies --- #

# Build 3 grafted phylognies (i.e., simulations) with the same backbone and constraints
grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
                                   nb_simulations = 3, verbose = FALSE, seed = 1234)

# Check all resulting grafted phylogenies at once
check_results_multi <- check_grafting_constraints(phy_grafted = grafted_phys, phy_backbone = phy,
                                                  missing_taxa_list = missing_taxa_list)
#> 3 phylogenies checked (no map), 4 taxa each
#>   - Phylogenies respecting all the constraints: 3/3
#>   Violations, over all the phylogenies (number of taxa x phylogenies):
#>   - Missing from the grafted tree: 0
#>   - Target clade violated:          0
#>   - Excluded subclade violated:     0
#>   - Backbone: preserved in 3/3 phylogenies, branch lengths shifted by max 9.65e-11 of backbone depth (tolerance = 1e-06)
#>   => All constraints respected

# With several simulations: the data.frame has the columns 'map', 'simulation' and 'seed',
# and the summary is printed once for all the phylogenies
head(check_results_multi[, c("map", "simulation", "seed", "Taxon_name", "ok")])
#>    map simulation      seed Taxon_name   ok
#> 1 <NA>          1      1234         X1 TRUE
#> 2 <NA>          1      1234         X2 TRUE
#> 3 <NA>          1      1234         X3 TRUE
#> 4 <NA>          1      1234         X4 TRUE
#> 5 <NA>          2 488349520         X1 TRUE
#> 6 <NA>          2 488349520         X2 TRUE

# An additional attribute provide a summary per phylogeny:
# number of taxa that failed, backbone preserved, all constraints respected
attr(check_results_multi, "phylogenies")
#>    map simulation      seed n_failed backbone_preserved all_ok
#> 1 <NA>          1      1234        0               TRUE   TRUE
#> 2 <NA>          2 488349520        0               TRUE   TRUE
#> 3 <NA>          3 469344150        0               TRUE   TRUE


# --- Check missing clades grafted with graft_missing_clades() --- #

# With a 'Sister_clade' item, the missing clades are checked (one row per missing clade)

# Prepare a missing_taxa_list for missing clades
missing_clades_list <- list(
  Taxon_name = c("M1", "M2", "M3"),
  Sister_clade = list(M1 = c("t26", "t18"), M2 = c("t26", "t18"), M3 = "t5"))

# Graft the missing clades with graft_missing_clades()
grafted_clades <- graft_missing_clades(phy = phy, missing_taxa_list = missing_clades_list,
                                       verbose = FALSE, seed = 1234)

# Check results
check_clades <- check_grafting_constraints(phy_grafted = grafted_clades, phy_backbone = phy,
                                           missing_taxa_list = missing_clades_list)
#> 2 missing clades (3 taxa) checked
#>   - Missing from the grafted tree:        0
#>   - Not forming one clade:                0
#>   - Not sister to their sister clade:     0
#>   - Backbone: topology preserved, branch lengths shifted by max 4.82e-09 (9.65e-11 of backbone depth; tolerance = 1e-06)
#>   => All constraints respected
# Inspect results
check_clades[, c("missing_clade", "Taxon_names", "monophyletic", "sister_position", "ok")]
#>   missing_clade Taxon_names monophyletic sister_position   ok
#> 1             1      M1, M2         TRUE            stem TRUE
#> 2             2          M3         TRUE      tip_branch TRUE


# --- Check the constraints from several traits, with biogeographic ranges --- #

# With the stochastic maps of several traits used to assign grafting constraints,
# 'allowed_states' has one element per trait,
# and 'treat_states_as_ranges' is a logical vector with one value per trait
# used to identify which trait was treated as biogeographic ranges.

# Load datasets
data(toy_backbone)
data(toy_simmap_ranges)
data(toy_simmap_elevation)
data(toy_missing_taxa_list_with_ranges_and_elevation)

# Graft phylogeny using elevation and biogeographic ranges as constraints
grafted_phy_traits <- graft_missing_taxa(phy = toy_backbone,
    missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
    simmap = list(biogeo = toy_simmap_ranges, elevation = toy_simmap_elevation),
    allowed_states = list(c("A", "B"), NULL),
    treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE),
    verbose = FALSE, seed = 1234)

# The result has the joint columns 'states_ok' and 'ranges_ok' (ranges of the 'range trait'),
# and the column 'states_ok_<trait>' for each trait. The summary details the failures per trait.
check_traits <- check_grafting_constraints(
    phy_grafted = grafted_phy_traits,
    phy_backbone = toy_backbone,
    missing_taxa_list = toy_missing_taxa_list_with_ranges_and_elevation,
    allowed_states = list(c("A", "B"), NULL),
    treat_states_as_ranges = c(biogeo = TRUE, elevation = FALSE))
#> 7 taxa checked
#>   - Missing from the grafted tree: 0
#>   - Target clade violated:          0
#>   - Excluded subclade violated:     0
#>   - Allowed states violated:        0 (biogeo: 0, elevation: 0)
#>   - Range not in graft-point state: 0
#>   - Backbone: topology preserved, branch lengths shifted by max 4.53e-09 (9.06e-11 of backbone depth; tolerance = 1e-06)
#>   => All constraints respected
check_traits[, c("Taxon_name", "states_ok", "states_ok_biogeo", "states_ok_elevation",
                 "ranges_ok", "ok")]
#>         Taxon_name states_ok states_ok_biogeo states_ok_elevation ranges_ok
#> New_sp1    New_sp1      TRUE             TRUE                TRUE      TRUE
#> New_sp2    New_sp2      TRUE             TRUE                TRUE      TRUE
#> New_sp3    New_sp3      TRUE             TRUE                TRUE      TRUE
#> New_sp4    New_sp4      TRUE             TRUE                TRUE      TRUE
#> New_sp5    New_sp5      TRUE             TRUE                TRUE      TRUE
#> New_sp6    New_sp6      TRUE             TRUE                TRUE      TRUE
#> New_sp7    New_sp7      TRUE             TRUE                TRUE      TRUE
#>           ok
#> New_sp1 TRUE
#> New_sp2 TRUE
#> New_sp3 TRUE
#> New_sp4 TRUE
#> New_sp5 TRUE
#> New_sp6 TRUE
#> New_sp7 TRUE

# --- Check both types of missing lineages grafted in one run --- #

# With list(taxa, clades), both parts are checked on the same phylogenies
both_list <- list(
  taxa = missing_taxa_list,
  clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = c("t8", "t22"), M2 = "t5")))

grafted_both <- graft_missing_taxa_and_clades(phy = phy, missing_taxa_list = both_list,
                                              verbose = FALSE, seed = 1234)

checks_both <- check_grafting_constraints(phy_grafted = grafted_both, phy_backbone = phy,
                                          missing_taxa_list = both_list, verbose = FALSE)
names(checks_both)
#> [1] "taxa"   "clades"
attr(checks_both, "all_ok") # "all_ok" is TRUE if both types of missing lineages are ok
#> [1] TRUE

```
