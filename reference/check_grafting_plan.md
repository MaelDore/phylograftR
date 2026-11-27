# Check the feasibility of a grafting plan before grafting

Check, before any grafting, whether the taxa (or missing clades) of
`missing_taxa_list` can be grafted\*\* on the backbone `phy` with the
constraints chosen (target clades, excluded subclades, stem
availability, species attributes), and report what is available for each
of them: the number of branches, the share of the available branches
that is in the allowed states (if `simmap` is used as a constraint), the
time window of the divergence ages and the share of it where grafting is
possible, and the allowed states found along the valid branches.

The goal is to **help fine-tune the grafting constraints before
committing to a long run**: the same checks are run at the start of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
and
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
(which stop if a taxon cannot be grafted), but here **the problems are
reported for all the taxa at once**, and the table can be explored.

## Usage

``` r
check_grafting_plan(
  phy,
  missing_taxa_list,
  min_tips = 5,
  include_target_clade_stem = TRUE,
  use_excluded_subclades = TRUE,
  exclude_subclades_stem = TRUE,
  simmap = NULL,
  allowed_states = NULL,
  treat_states_as_ranges = FALSE,
  min_valid_fraction = 0.1,
  ultrametric_tol = .Machine$double.eps^0.5,
  verbose = TRUE
)
```

## Arguments

- phy:

  An object of class `"phylo"`. The backbone phylogeny: time-calibrated
  (ultrametric) and binary.

- missing_taxa_list:

  A list, as given to the grafting functions (see
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  and
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).
  Items: `$Taxon_name`, and `$Target_clade` (taxa, for
  [`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md))
  or `$Sister_clade` (missing clades, for
  [`graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)),
  with the optional items of the grafting function. Exactly one of
  `$Target_clade` and `$Sister_clade` must be provided: it defines the
  type of the plan.

  It can also be **the list of both types of lineages to graft**,
  `list(taxa = ..., clades = ...)`, as returned by
  [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  when both types are provided, to be used with
  [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md).
  Each part is then planned on the same backbone and the output is a
  named list (see 'Value').

- min_tips:

  Numeric. Same as in the grafting functions. With
  `method = "BD_clades"`, a birth-death model is fitted on the target
  clade (or the sister clade) of each taxon, or on the smallest
  enclosing clade with at least `min_tips` tips. Here, **no model is
  fitted: the enclosing clade is only reported.** Default is `5`.

- include_target_clade_stem, use_excluded_subclades,
  exclude_subclades_stem:

  Same as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md).
  Only used with a `$Target_clade` list.

- simmap:

  `NULL`, or stochastic character map(s) of the backbone, as in the
  grafting functions (a `"simmap"`, a `"multiSimmap"`, or a list of them
  for several traits). Used to constrain grafting based on the matching
  of the `allowed_states` of each taxa with the states found along
  branches of the stochastic character map(s).

- allowed_states, treat_states_as_ranges:

  Same as in the grafting functions. Only used with a `simmap`.

- min_valid_fraction:

  Numeric between 0 and 1. Taxa (or missing clades) with a share of
  available branches in allowed states below this value are flagged as
  `low_valid_fraction` (a warning in the grafting functions). Default is
  `0.1`.

- ultrametric_tol:

  Numeric. Tolerance used by
  [`ape::is.ultrametric()`](https://rdrr.io/pkg/ape/man/is.ultrametric.html)
  to check for ultrametricity. Default is `.Machine$double.eps^0.5`.

- verbose:

  Logical. Print a summary? Default is `TRUE`.

## Value

A data.frame with one row per taxon (`$Target_clade` list) or per
missing clade (`$Sister_clade` list). With a `list(taxa, clades)`
covering both types of missing lineages, the output is a named list with
the elements `taxa` and `clades`, each data.frame containing the output
of its type (as described below), and the overarching attribute
`"all_feasible"` (`TRUE` if all the missing taxa and missing clades can
be grafted).

Columns for missing taxa:

- `Taxon_name` Character string. Name of the taxon.

- `Target_clade` Character string. Tips defining the target clade (first
  five).

- `target_n_tips` Integer. Number of tips of the target clade in the
  backbone.

- `n_branches` Integer. Number of branches available for grafting
  (before the restriction to the allowed states).

- `window_age` Numeric. Upper end (age) of the time window of the
  divergence ages: the stem age of the target clade, or its crown age if
  the stem is not available for grafting.

- `enclosing_n_tips`, `enclosing_age` Integer and numeric. Number of
  tips in the backbone phylogeny, and crown age of the enclosing clade
  (see Details).

Columns for missing clades:

- `Missing_clade` Integer. ID of the missing clade (taxa with the same
  sister clade).

- `Taxa` Character string. Taxa of the missing clade, separated by
  commas.

- `Sister_clade` Character string. Tips defining the sister clade (first
  five).

- `sister_n_tips` Integer. Number of tips of the sister clade in the
  backbone.

- `stem_length` Numeric. Length of the stem branch of the sister clade.

- `stem_age` Numeric. Age of the rootward end of the stem branch: the
  upper end of the divergence ages.

- `enclosing_n_tips`, `enclosing_age` Integer and numeric. Number of
  tips in the backbone phylogeny, and crown age of the enclosing clade
  (see Details).

Columns for any type, with a `simmap`:

- `allowed_states` Character string. The allowed states (resolved among
  the states of the simmap).

- `states_available` Character string. The states found in the available
  branches.

- `valid_fraction`, `time_coverage` Numeric. Share of the length of the
  available branches in allowed states, and share of the time window
  where at least one available branch is in an allowed state. With
  several maps: `valid_fraction_min`, `valid_fraction_mean`,
  `valid_fraction_max`, and the same for `time_coverage`.

- `n_maps_infeasible` Integer. With several maps, the number of maps
  where the taxon (or missing clade) cannot be grafted.

Last columns:

- `feasible` Logical. Can the taxon be grafted (in all the maps)?

- `cause` Character string. Why it cannot be grafted:

  - `"excluded"`: no branch left by the constraints;

  - `"excluded_contains_target"`: an excluded subclade contains the
    target clade;

  - `"excluded_identical"`: an excluded subclade is identical to the
    target clade;

  - `"allowed_states"`: no available branch has a piece in an allowed
    state;

  - `"whole_tree"`: the sister clade is the whoel tree, thus no stem
    branch exists for grafting;

  - `"inconsistent_allowed_states"`: the taxa of a missing clade have
    different allowed states;

  - `"low_valid_fraction"`: the share of available branche pieces in
    allowed states is below `min_valid_fraction`;

  - `"none"`: no issue. With several maps, the causes found across all
    maps are concatenated, separated by `"|"`.

Attributes: `"type"` (`"taxa"` or `"clades"`), `"min_tips"`,
`"min_valid_fraction"`, `"n_maps"`, and `"per_map"` (with several maps:
the diagnostic table of each map, with a `map` column).

## Details

### Taxa (`$Target_clade` list)

For each taxon, the pool of branches available for grafting is built
from its target clade, and restricted by its excluded subclades (see
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).
With a `simmap`, the pool is then restricted to the pieces of the
branches in the allowed states of the taxon.

### Missing clades (`$Sister_clade` list)

**Taxa with the same sister clade are grafted together as one missing
clade.** The missing clade is grafted on the stem branch of its sister
clade (see
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)).
With a `simmap` the share of the stem branch in the allowed states is
both the valid fraction and the share of the time window available. The
taxa of a missing clade must share their allowed states. Results of the
evaluation are reported with one row per missing clade.

### Enclosing clade

The enclosing clade is the target clade (or sister clade), or the
smallest clade enclosing it with at least `min_tips` tips, in the
backbone. With `method = "BD_clades"`, it is used to fit the Birth-Death
model and estimate divergence ages of the missing taxa to graft in
accordance to the fitted diversification process. It is reported so
users can tune `min_tips` accordingly to define more relevant enclosing
clades. During the grafting, taxa already grafted in a clade count in
its number of tips. The chosen enclosing clade used can therefore be the
reported clade or a clade nested within it. It can contain more tips
than the one reported here, although the reported clade is the oldest
one that can be used, and the number of tips reported is its number of
backbone tips.

**The sampling fractions of the birth-death models are not reported**:
they depend on the taxa already grafted in the clade, so they cannot be
known before grafting. Use
[`phylograftR::extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md)
to extract **a posteriori** the effective sampling fraction and
diversification rates used during the grafting process.

### Several maps

With a `"multiSimmap"`, the plan is checked for each map, as
attribute-based constraints are different across maps, and the table
summarizes values across all maps (i.e., minimum, mean and maximum of
the valid fraction and of the time coverage, and the number of maps in
which the taxon cannot be grafted). The table per map is stored in the
attribute **`"per_map"`**.

## See also

Other functions from phylograftR:

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades.

- [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft groups of missing taxa as sister clades of backbone clades.

- [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)
  to graft missing taxa and missing clades in one run, with a single set
  of inputs.

- [`phylograftR::check_grafting_constraints()`](https://maeldore.github.io/phylograftR/reference/check_grafting_constraints.md)
  to check a posteriori that grafted taxa respect their constraints,
  after grafting.

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

# Create a toy backbone phylogeny and associated list of missing taxa
set.seed(1234)
phy <- ape::rcoal(n = 30)
phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50

missing_taxa_list <- list(
  Taxon_name = c("X1", "X2", "X3"),
  Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5"))

# --- Check the grafting plan for taxonomic constraints --- #

# Check the grafting plan a priori to detect any issue at low cost
plan <- check_grafting_plan(phy = phy, missing_taxa_list = missing_taxa_list, min_tips = 5)
#> Grafting plan: 3 taxa
#>   - Can be grafted: 3/3
#>   - Target clades smaller than 'min_tips' = 5 (birth-death model fitted on an enclosing clade): 1/3; enclosing clades of 13 tips in the backbone; in practice, the clade used can be nested within them, as previously grafted taxa count towards 'min_tips'
#>   => All can be grafted with these constraints
plan
#>   Taxon_name Target_clade target_n_tips n_branches window_age enclosing_n_tips
#> 1         X1     t26, t18            17         33 50.0000000               17
#> 2         X2      t8, t22            13         25 50.0000000               13
#> 3         X3           t5             1          1  0.6368624               13
#>   enclosing_age feasible cause
#> 1      20.03438     TRUE  none
#> 2      15.28106     TRUE  none
#> 3      15.28106     TRUE  none

# --- Check the grafting plan for attribute-based constraints too --- #

# Add a simmap to restrict grafting to some states
Q <- matrix(0.04, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
diag(Q) <- -0.08
simmap <- phytools::sim.history(phy, Q, anc = "A", message = FALSE)

# Check the grafting plan a priori
plan_states <- check_grafting_plan(phy = phy, missing_taxa_list = missing_taxa_list,
                                   simmap = simmap, allowed_states = "B")
#> Grafting plan: 3 taxa
#>   - Can be grafted: 2/3
#>   - Cannot be grafted (1):
#>       allowed_states              1
#>       e.g., X3
#>   - Share of the available branches in allowed states (feasible taxa): min 8.1%, median 13.4%
#>   - Share of the time window with grafting possible (feasible taxa): min 8.6%, median 19.1%
#>   - Below 'min_valid_fraction' = 0.1 (1): X1
#>   - Target clades smaller than 'min_tips' = 5 (birth-death model fitted on an enclosing clade): 1/3; enclosing clades of 13 tips in the backbone; in practice, the clade used can be nested within them, as previously grafted taxa count towards 'min_tips'
#>   => Some cannot be grafted with these constraints (see the table)
# Inspect availability of branches in the allowed states
plan_states[, c("Taxon_name", "states_available", "valid_fraction", "time_coverage", "feasible")]
#>   Taxon_name states_available valid_fraction time_coverage feasible
#> 1         X1          A, B, C     0.08115666    0.08617336     TRUE
#> 2         X2          A, B, C     0.18601114    0.29636500     TRUE
#> 3         X3                C     0.00000000    0.00000000    FALSE

# --- Check both types of missing lineages at once --- #

# With list(taxa, clades), each part is planned on the same backbone: the output is a named list
both_list <- list(
  taxa = missing_taxa_list,
  clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = c("t8", "t22"), M2 = "t5")))
plans <- check_grafting_plan(phy = phy, missing_taxa_list = both_list, verbose = FALSE)
names(plans)
#> [1] "taxa"   "clades"
# "all_feasible" is TRUE if both types of missing lineages are feasible
attr(plans, "all_feasible")
#> [1] TRUE

```
