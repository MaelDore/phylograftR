# Toy taxonomic tables

Tables to build a `missing_taxa_list` from taxonomic names with
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md),
consistent with
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md).

`toy_backbone_taxonomic_table` gives the taxonomy of the 30 tips of the
backbone, with four levels (`Family`, `Subfamily`, `Tribe`, `Genus`).
**All the groups are monophyletic.** Some tips have no group at the
finest levels (`NA`), as often in real data.

`toy_missing_taxa_taxonomic_table` describes 10 taxa to graft. The first
seven (**`Type = "taxa"`**) are the taxa of
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md),
written with taxonomic names: the clades, excluded subclades and stem
settings are the same. The last three (**`Type = "clades"`**) are taxa
to graft as sister clades with
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md):
`Genus_a_sp1` and `Genus_a_sp2` belong to a new genus sister to the
genus `Gen_1`, so they are grafted together as a missing clade, and
`Genus_b_sp1` is sister to its own tribe. The target and sister clades
are written with the name of a group (e.g., `"Family_B"`) or of a level
(e.g., `"Family"`, meaning the family of the taxon itself). Both tables
are also provided as `.csv` files in
`system.file("extdata", package = "phylograftR")`
(`toy_backbone_taxonomic_table.csv` and
`toy_missing_taxa_taxonomic_table.csv`).

**The clades are named after their taxonomic groups** (e.g.,
`"Family_A"`), as in
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md)
and
[`phylograftR::toy_total_richness_per_clades`](https://maeldore.github.io/phylograftR/reference/toy_total_richness_per_clades.md),
so the lists built from the taxonomy can be used with this richness
list.

## Usage

``` r
data(toy_backbone_taxonomic_table)
data(toy_missing_taxa_taxonomic_table)

toy_backbone_taxonomic_table

toy_missing_taxa_taxonomic_table
```

## Format

`toy_backbone_taxonomic_table` is a data frame with 30 rows (one per tip
of the backbone) and the columns `Taxa` (tip labels), `Family`,
`Subfamily`, `Tribe` and `Genus`.

`toy_missing_taxa_taxonomic_table` is a data frame with 10 rows (one per
taxon to graft) and the columns `Taxa`, `Type`, the taxonomic levels
`Family`, `Subfamily`, `Tribe` and `Genus` of the taxa (`NA` if
unknown), `Target_clade`, `Sister_clade`, `Excluded_subclades`,
`include_target_clade_stem` and `exclude_subclades_stem`.

An object of class `data.frame` with 30 rows and 5 columns.

An object of class `data.frame` with 10 rows and 11 columns.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
  for the backbone phylogeny.

- [`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md)
  for the list of missing taxa to graft on the backbone.

- [`phylograftR::import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md)
  to read such tables from files.

- [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
  to build `missing_taxa_list` objects from these tables.

## Examples

``` r
data(toy_backbone)
data(toy_backbone_taxonomic_table)
data(toy_missing_taxa_taxonomic_table)

head(toy_backbone_taxonomic_table)
#>   Taxa   Family Subfamily     Tribe Genus
#> 1  t26 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 2   t2 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 3  t21 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 4  t16 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 5  t11 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 6   t4 Family_A Subfam_A2 Tribe_A2a  <NA>
toy_missing_taxa_taxonomic_table
#>           Taxa   Type   Family Subfamily     Tribe   Genus Target_clade
#> 1      New_sp1   taxa Family_A      <NA>      <NA>    <NA>       Family
#> 2      New_sp2   taxa Family_A      <NA>      <NA>    <NA>       Family
#> 3      New_sp3   taxa Family_A      <NA>      <NA>    <NA>       Family
#> 4      New_sp4   taxa Family_B      <NA>      <NA>    <NA>     Family_B
#> 5      New_sp5   taxa Family_B Subfam_B1      <NA>    <NA>    Subfamily
#> 6      New_sp6   taxa Family_B Subfam_B2  Tribe_B2    <NA>        Gen_2
#> 7      New_sp7   taxa Family_B Subfam_B2  Tribe_B2    <NA>        Tribe
#> 8  Genus_a_sp1 clades Family_B Subfam_B1      <NA> Genus_a         <NA>
#> 9  Genus_a_sp2 clades Family_B Subfam_B1      <NA> Genus_a         <NA>
#> 10 Genus_b_sp1 clades Family_A Subfam_A1 Tribe_A1b Genus_b         <NA>
#>    Sister_clade  Excluded_subclades include_target_clade_stem
#> 1          <NA> Subfam_A1;Tribe_A2a                      TRUE
#> 2          <NA>     Tribe_A2a;Gen_3                      TRUE
#> 3          <NA>                <NA>                     FALSE
#> 4          <NA>                <NA>                        NA
#> 5          <NA>               Gen_1                      TRUE
#> 6          <NA>                <NA>                        NA
#> 7          <NA>                <NA>                     FALSE
#> 8         Gen_1                <NA>                        NA
#> 9         Gen_1                <NA>                        NA
#> 10        Tribe                <NA>                        NA
#>    exclude_subclades_stem
#> 1              TRUE;FALSE
#> 2              FALSE;TRUE
#> 3                    <NA>
#> 4                    <NA>
#> 5                    <NA>
#> 6                    <NA>
#> 7                    <NA>
#> 8                    <NA>
#> 9                    <NA>
#> 10                   <NA>

# Build the lists: one for graft_missing_taxa() and one for graft_missing_clades()
toy_lists <- build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table = toy_backbone_taxonomic_table,
  missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
  phy = toy_backbone)
#> 7 taxa for graft_missing_taxa() and 3 for graft_missing_clades() built from 10 taxonomic group(s), 0 not monophyletic. See attr(, "taxonomic_groups").
str(toy_lists$taxa, max.level = 1)
#> List of 6
#>  $ Taxon_name               : chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Target_clade             :List of 7
#>  $ Target_clade_name        : Named chr [1:7] "Family_A" "Family_A" "Family_A" "Family_B" ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ include_target_clade_stem: Named logi [1:7] TRUE TRUE FALSE NA TRUE NA ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Excluded_subclades       :List of 7
#>  $ exclude_subclades_stem   :List of 7
str(toy_lists$clades, max.level = 1)
#> List of 3
#>  $ Taxon_name       : chr [1:3] "Genus_a_sp1" "Genus_a_sp2" "Genus_b_sp1"
#>  $ Sister_clade     :List of 3
#>  $ Sister_clade_name: Named chr [1:3] "Gen_1" "Gen_1" "Tribe_A1b"
#>   ..- attr(*, "names")= chr [1:3] "Genus_a_sp1" "Genus_a_sp2" "Genus_b_sp1"

# Compare to the initial '.csv' files
system.file("extdata", "toy_backbone_taxonomic_table.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_backbone_taxonomic_table.csv"
system.file("extdata", "toy_missing_taxa_taxonomic_table.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_missing_taxa_taxonomic_table.csv"
```
