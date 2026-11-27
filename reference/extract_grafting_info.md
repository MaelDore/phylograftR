# Extract the grafting information from grafted phylogenies

**Summarize, in a table, where and when each grafted taxon (or missing
clade) was grafted on the backbone**: the age of the graft, the backbone
branch it was grafted on (the parent clade), and, with a `simmap`, the
state at the graft point and the state inherited by the grafted lineage.
**The function works a posteriori**, on the output of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
or
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
(a phylogeny, a set of simulations, or the list of maps obtained with a
`"multiSimmap"`).

The table is meant to check the plausibility of the grafting against the
backbone, and, over a set of phylogenies, to describe the uncertainty of
the grafting (e.g., the distribution of the graft ages and graft states
of each taxon).

## Usage

``` r
extract_grafting_info(phy_grafted, missing_taxa_list, phy_backbone = NULL)
```

## Arguments

- phy_grafted:

  A grafted phylogeny (class `"phylo"`), a set of grafted phylogenies
  (class `"multiPhylo"`, e.g., several simulations), or a list of them
  with one element per map (output of the grafting functions with a
  `"multiSimmap"`).

- missing_taxa_list:

  The list of missing taxa used for the grafting. Items: `$Taxon_name`,
  and `$Target_clade` (taxa, grafted with
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md))
  or `$Sister_clade` (missing clades, grafted with
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).
  **Exactly one of them must be provided per list**.

  Alternatively, you can provide **the list of both types of missing
  lineages taxa** as `list(taxa = ..., clades = ...)`, as returned by
  [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  when both types are provided (e.g., to use with
  [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)).
  The output is then a named list (see 'Value').

- phy_backbone:

  Optional. The backbone phylogeny used for the grafting (class
  `"phylo"`). If `NULL` (default), the backbone is rebuilt by removing
  the grafted taxa from the grafted phylogeny (same topology and branch
  lengths): the taxa of `missing_taxa_list` (both types if a
  `list(taxa, clades)` is provided), and the tips flagged in
  `$grafted_tips`. **The node IDs of the backbone (`parent_node`) can
  only be reported if `phy_backbone` is provided.**

## Value

A data.frame, with one row per taxon (or per missing clade) and
phylogeny.

With a `list(taxa, clades)` covering both types of missing lineages, the
output is a named list with the elements `taxa` and `clades`, each
data.frame containing the output of its type (as described below).

The columns are:

- `map` Character string. Name of the stochastic map (`NA` if no map was
  used).

- `simulation` Integer. Index of the simulation.

- `seed` Integer. Seed of the simulation (`NA` if unknown).

For **missing taxa**:

- `Taxon_name` Character string. Name of the taxon.

- `parent_node` Integer. Node ID in `phy_backbone` of the tipward end of
  the backbone branch of the graft point (`NA` without `phy_backbone`).

- `parent_tips` Character string. Tips defining the parent clade.

- `graft_age` Numeric. Age of the graft point.

- `divergence_age` Numeric. Age of the node where the taxon splits from
  its sister lineage.

- `grafted_on` Character string. Taxa grafted on the sister lineage
  (first five), `NA` if the taxon is attached directly to the backbone.

For **missing clades**:

- `Missing_clade` Integer. Clade ID.

- `Taxa` Character string. Taxa of the missing clade, separated by
  commas.

- `n_taxa` Integer. Number of taxa in the grafted clade.

- `parent_node` Integer. Node ID in `phy_backbone` of the tipward end of
  the backbone branch of the graft point (`NA` without `phy_backbone`).

- `parent_tips` Character string. Tips defining the parent clade.

- `graft_age` Numeric. Age of the graft point = stem age of the grafted
  clade.

- `crown_age` Numeric. Crown age of the grafted clade.

With a `simmap`: `parent_state` and `tip_state` (taxa) or `clade_state`
(missing clades), or one column per trait.

Attributes: `"type"` (`"taxa"` or `"clades"`).

## Details

### What is extracted depends on the type of grafting

**Missing taxa** (`$Target_clade` list): **one row per taxon**. A taxon
can be grafted on a branch of the backbone, or on a branch of a taxon
grafted before it (and so on): the grafted taxa that hang together form
a *grafted subtree* attached to the backbone at a *graft point*. The
summary table reports:

- `graft_age`: the age of the graft point, where the grafted subtree is
  attached to a branch of the backbone;

- `divergence_age`: the age of the node where the taxon splits from its
  sister lineage. It is equal to `graft_age` for a taxon attached
  directly to the backbone, and younger for a taxon attached to another
  grafted taxon (`grafted_on` then lists the grafted taxa of its sister
  lineage).

- The parent clade is the clade of the backbone subtended by the branch
  of the graft point: `parent_tips` gives two tips defining it (one if
  the branch is terminal), and `parent_node` its node ID in
  `phy_backbone`.

**Missing clades** (`$Sister_clade` list): **one row per missing clade**
(the taxa with the same sister clade, grafted together). The missing
clade is attached to the stem branch of its sister clade at a *graft
point*. The summary table reports:

- `graft_age`: the age of the graft point;

- `crown_age`: the age of the most recent common ancestor of the taxa of
  the clade (`NA` for a single taxon);

- The parent clade is the backbone clade subtended by the branch it was
  attached to (the sister clade, if the constraints were respected, see
  [`phylograftR::check_grafting_constraints()`](https://maeldore.github.io/phylograftR/reference/check_grafting_constraints.md)).
  **If the taxa of a missing clade do not form a clade in the phylogeny,
  the values of the row are `NA`, with a warning.**

### States

If the phylogeny has a `$simmap` (output of the grafting with a `simmap`
and `return_simmap = TRUE`), the states are extracted: `parent_state` is
the state of the branch of the backbone at the graft point, and
`tip_state` (taxa) or `clade_state` (missing clades) the state inherited
by the grafted lineage.

**They are different if `treat_states_as_ranges = TRUE`**: a taxon
grafted on a branch in the range "AB" can inherit the range "A". With
the maps of several traits (`$simmap` is a list), there is one column
per trait (`parent_state_<trait>`, `tip_state_<trait>`). **The states
are `NA` for phylogenies without `$simmap`**, or whose `$simmap` does
not match the phylogeny (e.g., after removing tips).

## See also

Other functions from phylograftR:

- [`phylograftR::check_grafting_plan()`](https://maeldore.github.io/phylograftR/reference/check_grafting_plan.md)
  to check a priori the feasibility of a grafting plan, before grafting.

- [`phylograftR::check_grafting_constraints()`](https://maeldore.github.io/phylograftR/reference/check_grafting_constraints.md)
  to check a posteriori that grafted taxa respect their constraints.

- [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  to extract the richness and the sampling fraction of the clades in
  each grafted phylogeny.

- [`phylograftR::extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md)
  to extract a posteriori information on the birth-death models fitted
  during grafting.

- [`phylograftR::identify_grafted()`](https://maeldore.github.io/phylograftR/reference/identify_grafted.md)
  to identify the grafted tips and branches of a phylogeny.

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
  Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")))

# --- Graft the taxa (several simulations) --- #

grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
                                  nb_simulations = 5, verbose = FALSE, seed = 1234)

# --- Extract the grafting information --- #

info <- extract_grafting_info(phy_grafted = grafted_phy, missing_taxa_list = missing_taxa_list,
                              phy_backbone = phy)
head(info)
#>    map simulation      seed Taxon_name parent_node parent_tips  graft_age
#> 1 <NA>          1      1234         X1          52     t6, t20 0.74762887
#> 2 <NA>          1      1234         X2          30         t22 1.85178529
#> 3 <NA>          1      1234         X3          27          t5 0.00441298
#> 4 <NA>          1      1234         X4           9         t28 2.95133160
#> 5 <NA>          2 488349520         X1          15         t24 0.53229148
#> 6 <NA>          2 488349520         X2          59    t13, t17 0.38557333
#>   divergence_age grafted_on
#> 1     0.74762887       <NA>
#> 2     1.85178529       <NA>
#> 3     0.00441298       <NA>
#> 4     2.95133160       <NA>
#> 5     0.53229148       <NA>
#> 6     0.38557333       <NA>

# Inspect uncertainty of the graft age of each taxon over the simulations
# Average minimum and maximum graft age across simulations
aggregate(graft_age ~ Taxon_name, data = info, FUN = function (x) round(range(x), 1))
#>   Taxon_name graft_age.1 graft_age.2
#> 1         X1         0.5         4.3
#> 2         X2         0.1         2.4
#> 3         X3         0.0         0.6
#> 4         X4         0.4         3.0

# --- Both types of missing lineages, grafted in one run --- #

# With list(taxa, clades), the information is extracted for each type, with the same backbone
both_list <- list(
  taxa = missing_taxa_list,
  clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = c("t8", "t22"), M2 = "t5")))
grafted_both <- graft_missing_taxa_and_clades(phy = phy, missing_taxa_list = both_list,
                                              nb_simulations = 3, verbose = FALSE, seed = 1234)
info_both <- extract_grafting_info(phy_grafted = grafted_both, missing_taxa_list = both_list,
                                   phy_backbone = phy)
names(info_both)
#> [1] "taxa"   "clades"
head(info_both$clades)
#>    map simulation      seed Missing_clade Taxa n_taxa parent_node parent_tips
#> 1 <NA>          1      1234             1   M1      1          33      t8, t5
#> 2 <NA>          1      1234             2   M2      1          27          t5
#> 3 <NA>          2 488349520             1   M1      1          33      t8, t5
#> 4 <NA>          2 488349520             2   M2      1          27          t5
#> 5 <NA>          3 469344150             1   M1      1          33      t8, t5
#> 6 <NA>          3 469344150             2   M2      1          27          t5
#>      graft_age crown_age
#> 1 27.886880636        NA
#> 2  0.002746001        NA
#> 3 19.135357456        NA
#> 4  0.075010203        NA
#> 5 38.930312846        NA
#> 6  0.099970264        NA

```
