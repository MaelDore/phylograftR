# Extract the richness and sampling fraction of clades from grafted phylogenies

Summarize, in each grafted phylogeny, the species of each clade listed
in `total_richness_per_clades`: species of the backbone, grafted
species, and unsampled species (species of the clade that are neither in
the backbone nor grafted). It then derives the **sampling fraction** of
each clade that can be used in downstream analyses to account for the
clade-specific sampling fractions (e.g., diversification models).

The function works **a posteriori**, on the output of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
or
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
(a phylogeny, a set of simulations, or the list of maps obtained with a
`"multiSimmap"`). The counts are taken in each grafted phylogeny. With
[`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md),
**they can differ among simulations and maps**, as a taxon can end up in
different clades. With
[`graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md),
**they are the same in all the phylogenies**, as each missing clade is
always grafted on the stem of its sister clade.

## Usage

``` r
extract_sampling_info(phy_grafted, total_richness_per_clades)
```

## Arguments

- phy_grafted:

  A grafted phylogeny (class `"phylo"`), a set of grafted phylogenies
  (class `"multiPhylo"`, e.g., several simulations), or a list of them
  with one element per map (output of the grafting functions with a
  `"multiSimmap"`). The phylogenies must have the elements
  `$grafted_tips` and `$grafted_edges` (see
  [`phylograftR::identify_grafted()`](https://maeldore.github.io/phylograftR/reference/identify_grafted.md)).

- total_richness_per_clades:

  A list with the total number of species of the clades, as used for the
  grafting. Items:

  - `$Target_clade` (clades defined for the taxa grafted with
    [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md))
    or `$Sister_clade` (sister clades of the missing clades grafted with
    [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).
    **Exactly one of them must be provided.** List of character vectors
    of tips of the backbone: a pair of tips defines the clade as their
    most recent common ancestor (MRCA), and a single tip defines a
    terminal branch. **Node IDs are not accepted**: convert them with
    [`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md).

  - `$Target_clade_name` or `$Sister_clade_name`: (Optional) Character
    vector of the names of the clades, same length as the clades (`NA`
    if a clade has no name). The names must be unique. They are reported
    in the `Clade_name` column of the output.

  - `$total_richness`: Vector of positive integers, same length as the
    clades. Total number of species of each clade, i.e., backbone
    species, grafted species, and unsampled species. Use `NA` if
    unknown.

  It can also be **the list of both types of missing lineages** as
  `list(taxa = ..., clades = ...)`, as returned by
  [`phylograftR::build_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/build_total_richness_per_clades.md)
  and
  [`phylograftR::import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md)
  when both types are provided (e.g., to use with
  [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)).

## Value

A data.frame, with one row per clade and phylogeny.

With a `list(taxa, clades)` covering both types of missing lineages, the
output is a named list with the elements `taxa` and `clades`, each
data.frame containing the following columns:

- `map` Character string. Name of the stochastic map (`NA` if no map was
  used).

- `simulation` Integer. Index of the simulation.

- `seed` Integer. Seed of the simulation (`NA` if unknown).

- `Clade_name` Character string. Name of the clade (`NA` if it has no
  name). Only if the names of the clades are provided in
  `total_richness_per_clades`.

- `Clade_tips` Character string. Tips used to identify the clade.

- `total_richness` Numeric. Total number of species of the clade:
  provided, or derived if `NA` was provided.

- `total_richness_provided` Logical. Whether the total richness of the
  clade was provided.

- `nb_backbone_species` Integer. Number of species of the clade that are
  in the backbone.

- `nb_grafted_species` Integer. Number of species of the clade that were
  grafted.

- `nb_unsampled_species` Numeric. Number of species of the clade that
  are neither in the backbone nor grafted. Negative (conflictual) if the
  clade has more species than expected.

- `sampling_fraction` Numeric. Fraction of the clade richness present in
  the grafted phylogeny,
  `(nb_backbone_species + nb_grafted_species) / total_richness`.

Attributes: `"type"` (`"taxa"` if the clades are `$Target_clade`,
`"clades"` if they are `$Sister_clade`).

## Details

### Species of a clade

**The species of a clade are the tips that descend from its MRCA in the
grafted phylogeny**: the backbone species of the clade, and the grafted
taxa that were grafted within it. It depends on the phylogeny, as a
taxon grafted in a target clade can fall in a clade nested within it,
and as the taxa grafted on the stem branch of a clade are sister to the
clade, so they are not part of it. The exception is the clade of a
single tip of `$Target_clade`: its terminal branch is its only branch,
so the taxa grafted on it are part of the clade. For `$Sister_clade`,
the missing clades are grafted on the stem branch of the sister clade
and are never part of it.

### Unsampled species and sampling fraction

For a clade with a total richness provided,
`nb_unsampled_species = total_richness - nb_backbone_species - nb_grafted_species`.
For a clade without a total richness (`NA`), the unsampled species are
those identified in the largest clades nested within it that have a
total richness (none if there is none), and its total richness is
derived as
`nb_backbone_species + nb_grafted_species + nb_unsampled_species` (it
can differ among phylogenies). The **sampling fraction** is
`(nb_backbone_species + nb_grafted_species) / total_richness`.

In problematic cases, the number of species of a clade can be higher
than its total richness: for example, a taxon known to belong to a genus
can be grafted in a subclade of this genus whose provided total richness
did not include it. The number of unsampled species is then negative and
the sampling fraction is higher than 1. These values are returned as
they are, with a warning listing the clades and the number of
phylogenies concerned, so the clades and the total richness provided can
be reviewed. When a clade without a total richness is derived from
nested clades, the negative values of the nested clades are counted as
`0`.

## See also

Other functions from phylograftR:

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades.

- [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft groups of missing taxa as sister clades on backbone clades.

- [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)
  to graft missing taxa and missing clades in one run, with a single set
  of inputs.

- [`phylograftR::extract_grafting_info()`](https://maeldore.github.io/phylograftR/reference/extract_grafting_info.md)
  to extract a posteriori information on the location of grafted
  lineages.

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

# Two nested clades: B is nested within A
clade_A <- c("t26", "t18") # 17 backbone species
clade_B <- c("t26", "t3")  # 9 backbone species, nested within A

# List the missing taxa to graft: 3 taxa in clade A and 1 taxon in clade B
missing_taxa_list <- list(
  Taxon_name = c("X1", "X2", "X3", "X4"),
  Target_clade = list(X1 = clade_A, X2 = clade_A, X3 = clade_A, X4 = clade_B),
  # Optional names of the clades (one per taxon): they must match those of total_richness_per_clades
  Target_clade_name = c("clade_A", "clade_A", "clade_A", "clade_B"))

# Total richness of the target clades
# NA = unknown richness (i.e., no unsampled species beyond nested clades)
total_richness_per_clades <- list(
  Target_clade = list(clade_A, clade_B),
  Target_clade_name = c("clade_A", "clade_B"),
  total_richness = c(25, 12))

# --- Graft the taxa (several simulations) --- #

grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
                                   total_richness_per_clades = total_richness_per_clades,
                                   nb_simulations = 5, verbose = FALSE, seed = 1234)

# --- Extract the sampling information --- #

sampling_info <- extract_sampling_info(phy_grafted = grafted_phys,
                                       total_richness_per_clades = total_richness_per_clades)
sampling_info[, c("simulation", "Clade_name", "nb_backbone_species", "nb_grafted_species",
                  "nb_unsampled_species", "sampling_fraction")]
#>    simulation Clade_name nb_backbone_species nb_grafted_species
#> 1           1    clade_A                  17                  4
#> 2           1    clade_B                   9                  3
#> 3           2    clade_A                  17                  4
#> 4           2    clade_B                   9                  1
#> 5           3    clade_A                  17                  4
#> 6           3    clade_B                   9                  2
#> 7           4    clade_A                  17                  4
#> 8           4    clade_B                   9                  1
#> 9           5    clade_A                  17                  4
#> 10          5    clade_B                   9                  3
#>    nb_unsampled_species sampling_fraction
#> 1                     4         0.8400000
#> 2                     0         1.0000000
#> 3                     4         0.8400000
#> 4                     2         0.8333333
#> 5                     4         0.8400000
#> 6                     1         0.9166667
#> 7                     4         0.8400000
#> 8                     2         0.8333333
#> 9                     4         0.8400000
#> 10                    0         1.0000000

# The sampling fraction of clade B varies among simulations:
# the taxa assigned to clade A can be grafted within clade B (nested in A, so they count for B)
# or in A outside of B. As such, the richness of B differs between grafted phylogenies,
# while all the taxa are always in A (constant sampling fraction).
tapply(X = sampling_info$sampling_fraction,
       INDEX = list(sampling_info$Clade_name, sampling_info$simulation),
       FUN = mean)
#>            1         2         3         4    5
#> clade_A 0.84 0.8400000 0.8400000 0.8400000 0.84
#> clade_B 1.00 0.8333333 0.9166667 0.8333333 1.00

# --- Both types of missing lineages, grafted in one run --- #

# With graft_missing_taxa_and_clades(), the totals of both types are provided as list(taxa, clades)
both_list <- list(
  taxa = missing_taxa_list,
  clades = list(Taxon_name = c("M1", "M2"), Sister_clade = list(M1 = clade_B, M2 = "t5")))
both_richness <- list(
  taxa = total_richness_per_clades,
  clades = list(Sister_clade = list(clade_B, "t5"), total_richness = c(12, 4)))

grafted_both <- graft_missing_taxa_and_clades(phy = phy, missing_taxa_list = both_list,
                                              total_richness_per_clades = both_richness,
                                              nb_simulations = 3, verbose = FALSE, seed = 1234)

# Output: a named list with the table of each type
sampling_both <- extract_sampling_info(phy_grafted = grafted_both,
                                       total_richness_per_clades = both_richness)
names(sampling_both)
#> [1] "taxa"   "clades"
sampling_both$clades[, c("simulation", "Clade_tips", "nb_grafted_species", "sampling_fraction")]
#>   simulation Clade_tips nb_grafted_species sampling_fraction
#> 1          1    t26, t3                  3         1.0000000
#> 2          1         t5                  0         0.2500000
#> 3          2    t26, t3                  1         0.8333333
#> 4          2         t5                  0         0.2500000
#> 5          3    t26, t3                  2         0.9166667
#> 6          3         t5                  0         0.2500000

```
