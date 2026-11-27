# Toy total richness of the target clades

\*List usable as `total_richness_per_clades` with
[`phylograftR::toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
and
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md).
It is built from the table
`system.file("extdata", "toy_total_richness.csv", package = "phylograftR")`,
with the tips of each clade separated by commas (see the examples).

## Usage

``` r
data(toy_total_richness_per_clades)
```

## Format

A list with the items `$Target_clade` (list of character vectors of
tips), `$Target_clade_name` (character vector of the names of the
clades, the same as in
[`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md))
and `$total_richness` (numeric vector, `NA` if unknown).

## See also

Other datasets and functions from phylograftR:

- [`phylograftR::toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md)
  for the list of missing taxa to graft on the backbone.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa within a target clade, using this dataset.

- [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  to extract the richness and the sampling fraction of the clades of
  grafted phylogenies.

## Examples

``` r
# Import the data
data(toy_total_richness_per_clades)
str(toy_total_richness_per_clades)
#> List of 3
#>  $ Target_clade     :List of 5
#>   ..$ : chr [1:2] "t26" "t18"
#>   ..$ : chr [1:2] "t8" "t22"
#>   ..$ : chr [1:2] "t8" "t7"
#>   ..$ : chr "t5"
#>   ..$ : chr [1:2] "t5" "t12"
#>  $ Target_clade_name: chr [1:5] "Family_A" "Family_B" "Subfam_B1" "Gen_2" ...
#>  $ total_richness   : int [1:5] 25 NA 14 NA 8

# Compare to the initial '.csv' file
system.file("extdata", "toy_total_richness.csv", package = "phylograftR")
#> [1] "/home/runner/.cache/R/renv/library/phylograftR-802cedf8/linux-ubuntu-resolute/R-4.6/x86_64-pc-linux-gnu/phylograftR/extdata/toy_total_richness.csv"

# Rebuild the data from the '.csv' table
tab <- read.csv(system.file("extdata", "toy_total_richness.csv", package = "phylograftR"))
total_richness_per_clades <- list(
  Target_clade = strsplit(tab$Target_clade, ",", fixed = TRUE),
  Target_clade_name = tab$Target_clade_name,
  total_richness = tab$total_richness)
```
