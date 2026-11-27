# Toy list of missing taxa

Valid `missing_taxa_list` for
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md),
with 7 taxa to graft by
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md).

It is the output of
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
applied to the table
`system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")`:
open that file to see how to fill the table. It shows target clades
defined by one or two tips, excluded subclades (named or not, with one
tip or several), and the stem settings. Clades and excluded subclades
are named after the taxonomic groups of
[`phylograftR::toy_backbone_taxonomic_table`](https://maeldore.github.io/phylograftR/reference/toy_taxonomic_tables.md)
(e.g., `"Family_A"`).

## Usage

``` r
data(toy_missing_taxa_list)
```

## Format

A named list with 7 taxa and the following items. **Every item is named
after the taxa (`$Taxon_name`), in the same order.**

- `$Taxon_name` Character vector. Names of the taxa to graft.

- `$Target_clade` List of character vectors of tips defining the target
  clade of each taxon (one tip: the terminal branch of this tip).

- `$Target_clade_name` Character vector. Labels of the target clades.
  Not used for the grafting, but checked against the names of
  `total_richness_per_clades`, if provided.

- `$include_target_clade_stem` Logical vector. Whether the stem of the
  target clade of each taxon is available for grafting (`NA` uses the
  global argument).

- `$Excluded_subclades` List (one element per taxon) of named lists of
  character vectors of tips defining the excluded subclades, or `NULL`.

- `$exclude_subclades_stem` List (one element per taxon) of logical
  vectors (one value per excluded subclade), or `NULL`.

See
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
for the use of the items in the grafting.

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
  for the backbone phylogeny.

- [`phylograftR::toy_total_richness_per_clades`](https://maeldore.github.io/phylograftR/reference/toy_total_richness_per_clades.md)
  for the total richness of the target clades.

- [`phylograftR::toy_taxonomic_tables`](https://maeldore.github.io/phylograftR/reference/toy_taxonomic_tables.md)
  for the taxonomic tables of the backbone and of the missing taxa.

- [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  to import such a list from an external table.

- [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  to build such a list from a data.frame.

## Examples

``` r
# Import the list from loaded datasets
data(toy_missing_taxa_list)
str(toy_missing_taxa_list, max.level = 1)
#> List of 6
#>  $ Taxon_name               : chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Target_clade             :List of 7
#>  $ Target_clade_name        : Named chr [1:7] "Family_A" "Family_A" "Family_A" "Family_B" ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ include_target_clade_stem: Named logi [1:7] TRUE TRUE FALSE NA TRUE NA ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Excluded_subclades       :List of 7
#>  $ exclude_subclades_stem   :List of 7

# Compare to the initial '.csv' file
system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_missing_taxa.csv"
```
