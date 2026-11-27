# Build the `total_richness_per_clades` object from a table of clades

Build the `total_richness_per_clades` object used by
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
and
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
to account for the species that are neither in the backbone nor in the
`missing_taxa_list` (i.e., unsampled species), from **a single table**
with one row per clade (e.g., read from a .csv file) and the total
number of species of the clade.

The clades can be written with tips, node IDs, or **names** (see
Details), and are matched with the clades of the `missing_taxa_list`
object built to be used in the same call to a grafting function: this
ensures that **the tips and the names of a clade are the same in all the
objects**. The function checks that the richness is **consistent with
the full sequence of grafting**, before any grafting.

The output follows the structure of `missing_taxa_list`: if it contains
both types of taxa (as returned by
[`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
and
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)),
a list for each grafting function is returned.

Use
[`phylograftR::import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md)
to read the table directly from a file.

## Usage

``` r
build_total_richness_per_clades(
  data,
  phy,
  missing_taxa_list,
  return_richness_table = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  A data.frame, with one row per clade. Columns are:

  - `Clade`: a monophyletic group from the backbone phylogeny defined by
    tips, a node ID, or a name (see Details), from which the total
    richness is known.

  - `total_richness`: Total number of species of the clade. It includes
    the backbone species, **all the taxa to graft that belong to the
    clade**, and the unsampled species that are known to exist but will
    not be grafted. Use `NA` if unknown.

  - `Clade_name` (Optional): Name of the clade, used to check matching
    with other inputs, and to report warnings in downstream functions.
    If the clade is written with a name in `Clade`, it is its name.

  The columns `Target_clade` and `Target_clade_name`, or `Sister_clade`
  and `Sister_clade_name`, are accepted instead of `Clade` and
  `Clade_name` (e.g., for tables prepared for a single grafting
  function). All the columns are read as text. Blank cells and `NA` are
  missing values.

- phy:

  An object of class `"phylo"`. The backbone phylogeny on which the taxa
  will be grafted. It is used to check the tips and clades, and to
  convert node IDs into tips.

- missing_taxa_list:

  A `missing_taxa_list`, or the list of two `missing_taxa_list` (`$taxa`
  and `$clades`) returned by
  [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md),
  [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  or
  [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md).
  It provides the tips and the names used for the clades when they are
  only named, and is used to check the validity of the richness values
  provided.

- return_richness_table:

  Logical. Whether to return the expected richness and sampling fraction
  of each clade in the attribute `"richness_table"` of the output. It is
  computed before any grafting, from the backbone, the
  `missing_taxa_list` and the total richness. It is a tool to **examine
  a priori the sampling fraction expected in each clade**. The actual
  **effective sampling fractions of the grafted phylogenies**,
  accessible with
  [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md),
  may differ (see 'Richness table' in Details). Default is `TRUE`.

- verbose:

  Logical. Whether to print the conversions of node IDs into tips (if
  any), and a summary of the table. Default is `TRUE`.

## Value

If `data` contains one type of taxa, a `total_richness_per_clades` list;
if it contains both types (`$taxa` and `$clades`), a list with two
`total_richness_per_clades` named `$taxa` (for
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md))
and `$clades` (for
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).

Each `total_richness_per_clades` list contains all the clades of the
table, and has the following items:

- `$Target_clade` (or `$Sister_clade`, for the list of missing clades)
  List of character vectors of tips defining the clades.

- `$Target_clade_name` (or `$Sister_clade_name`) Only if some clades
  have a name. Character vector (`NA` if a clade has no name).

- `$total_richness` Integer vector of the total richness of the clades
  (`NA` if unknown).

If `return_richness_table = TRUE`, the output has the attribute
`"richness_table"`, a data.frame with one row per clade, and the
columns:

- `Clade_name`: Name of the clade (`NA` if none).

- `Clade_tips`: Tips defining the clade, as in the output.

- `total_richness_provided`: Total richness of the table (`NA` if
  unknown).

- `total_richness`: Total richness as used for the grafting, i.e., the
  value provided, or, if unknown, the species known in the clade
  (backbone and taxa to graft) plus the unsampled species derived from
  the clades nested within it (none if there are no such clades).

- `nb_backbone_species`: Number of tips of the backbone in the clade.

- `nb_missing_species`: Number of taxa to graft expected in the clade,
  of both types (see 'Richness table' in Details).

- `nb_unsampled_species`: Number of unsampled species, i.e., species
  that are neither in the backbone nor to graft (`total_richness` -
  `nb_backbone_species` - `nb_missing_species`).

- `sampling_fraction`: Expected sampling fraction of the clade in the
  grafted phylogenies,
  `(nb_backbone_species + nb_missing_species) / total_richness`. It is
  computed **a priori**, and the effective sampling fraction may differ
  after grafting: use
  [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  to extract it **a posteriori**.

## Details

### Identifying clades

Clades can be defined by tips, a node ID, or a name:

- A pair of tips defines the clade as their MRCA node (e.g.,
  `"t26, t18"`), separated by commas. A single tip defines the "clade"
  as the terminal branch of this tip. **Tip labels must include at least
  one other character than digits.**

- A node ID of `phy` (digits only, e.g., `"34"`) also defines the clade,
  as in
  [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md).

- A single value that is not a tip of `phy` is a **clade name** (e.g.,
  `"Family_A"`). Names are matched with the names of the clades of
  `missing_taxa_list` (`$Target_clade_name`, `$Sister_clade_name`, and
  the names of the excluded subclades) and, for lists built from
  taxonomic tables, with all the taxonomic groups of the backbone
  taxonomy. The taxonomic groups are stored in the attribute
  `"taxonomic_groups"` of the object returned by
  [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md):
  they are only available if `missing_taxa_list` is this object itself,
  not one of its elements (`$taxa` or `$clades`) alone.

**Clades are matched by their MRCA node in `phy`**: a clade can only be
listed once, and the tips written in the output are those already used
in `missing_taxa_list` for the same clade, so the tips match in all the
objects. The grafting functions require **all the clades used in
`missing_taxa_list` to be listed**: the clades that are not in the table
are added with an unknown richness (`NA`), and their number is reported.

### Unsampled species

For a focal clade with a `total_richness` value provided:
`total_richness = nb_backbone + nb_missing + nb_unsampled`, with
`nb_missing` the number of lineages to graft within the focal clade, and
`nb_unsampled` the number of species in this clade that are not
represented in the backbone, nor will be grafted, but must be accounted
in the sampling fraction for when fitting diversification models.

For a clade with no `total_richness` value (`NA` or `NULL`),
`nb_unsampled` are estimated as the sum of `nb_unsampled` over the
largest clades nested within it with a `total_richness` value provided.
If there is none, `nb_unsampled = 0`. The rationale is that even if
`nb_unsampled` cannot be derived from its `total_richness` (as it is not
provided for this clade), it has to at least account for the
`nb_unsampled` species that have being identified in its own nested
clades.

### Full sequence of grafting

The same table is valid for the grafting of the taxa within target
clades and for the grafting of the missing clades as sister to backbone
clades, in either order. With the recommended sequence (taxa first, then
missing clades on the grafted phylogeny):

- When the taxa are grafted, the taxa to graft as missing clades are not
  in the backbone nor in the list, so they are counted among the
  unsampled species of the clades.

- When the missing clades are grafted, the taxa grafted before are tips
  of the backbone, and the unsampled species are the species still
  unaccounted for.

The function checks the richness against all the taxa to graft at once:
**the function stops if the total richness of a clade is lower than the
number of species it contains** (backbone species and all the taxa to
graft within it, or within clades nested within it), or if a clade has
fewer unsampled species than the clades with a total richness nested
within it.

### Richness table

With `return_richness_table = TRUE`, the output has the attribute
`"richness_table"`: a data.frame that reports, **before any grafting**,
the richness and the sampling fraction expected in each clade (see Value
for details on its columns). It can be used to **examine a priori the
sampling fraction of the clades** and to spot errors in the total
richness: e.g., an unlikely low sampling fraction, or clades with an
unknown richness (`NA`) for which unsampled species were derived from
nested clades, or assumed to be none. It can also be used to compare the
sampling of different clades, or to decide which clades need a better
estimate of their richness.

**The sampling fractions of the table are expectations, and the
effective sampling fractions of the grafted phylogenies may differ.**
The table assumes that all the taxa to graft end up in the clade where
they belong: a taxon counts for all the clades that are, or enclose, its
target clade, and a missing clade counts for the clades that enclose its
sister clade (it is grafted on the stem of the sister clade, so not
within the sister clade itself). However, where taxa are actually
grafted varies among the grafted phylogenies (and among the stochastic
maps, if any): a taxon can be grafted on the stem of its target clade,
which is outside of the clade, or within a clade nested in its target
clade, which changes the richness of these clades. The actual effective
richness and sampling fraction of each clade, in each grafted phylogeny,
are extracted **a posteriori** with
[`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
on the output of the grafting functions, using the same
`total_richness_per_clades`.

## See also

Other functions from phylograftR:

- [`phylograftR::import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md)
  to read the table directly from an external file.

- [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  and
  [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  to prepare the `missing_taxa_list`.

- [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
  to prepare the `missing_taxa_list` from taxonomic tables.

- [`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md)
  to convert the node IDs of a `total_richness_per_clades` written 'by
  hand'.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  and
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft the taxa using the richness.

- [`phylograftR::extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  to extract the richness and sampling fraction of the clades after
  grafting.

## Author

Maël Doré

## Examples

``` r
# Load the backbone, the taxonomic tables, and the toy lists of taxa and missing clades
data(toy_backbone)
data(toy_backbone_taxonomic_table)
data(toy_missing_taxa_taxonomic_table)

# Build the missing taxa list for both types of missing lineages (taxa and clades)
lists <- build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table = toy_backbone_taxonomic_table,
  missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
  phy = toy_backbone, verbose = FALSE)
names(lists)
#> [1] "taxa"   "clades"

# Inform the total richness of some clades, written with their taxonomic names.
# Use NA if the total richness is unknown
richness_table <- data.frame(Clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"),
                             total_richness = c(25, NA, 14, 8))

# Use the the function to validate this richness table
# against the backbone tree and missing taxon lists
richness_per_clades <- build_total_richness_per_clades(data = richness_table, phy = toy_backbone,
                                                       missing_taxa_list = lists)
#> 7 clade(s) in total richness: 3 with a value, 4 with an unknown richness (including 3 clade(s) of 'missing_taxa_list' added to the 4 of the table). See attr(, "richness_table") for the expected sampling fraction of the clades.

# Output = one list for each grafting function,
# with the clades of both lists and the same tips as in 'lists'
str(richness_per_clades, 2)
#> List of 2
#>  $ taxa  :List of 3
#>   ..$ Target_clade     :List of 7
#>   ..$ Target_clade_name: chr [1:7] "Family_A" "Family_B" "Subfam_B1" "Tribe_B2" ...
#>   ..$ total_richness   : int [1:7] 25 NA 14 8 NA NA NA
#>  $ clades:List of 3
#>   ..$ Sister_clade     :List of 7
#>   ..$ Sister_clade_name: chr [1:7] "Family_A" "Family_B" "Subfam_B1" "Tribe_B2" ...
#>   ..$ total_richness   : int [1:7] 25 NA 14 8 NA NA NA
#>  - attr(*, "richness_table")='data.frame':   7 obs. of  8 variables:
#>   ..$ Clade_name             : chr [1:7] "Family_A" "Family_B" "Subfam_B1" "Tribe_B2" ...
#>   ..$ Clade_tips             : chr [1:7] "t26, t19" "t8, t22" "t8, t7" "t5, t12" ...
#>   ..$ total_richness_provided: num [1:7] 25 NA 14 8 NA NA NA
#>   ..$ total_richness         : num [1:7] 25 24 14 8 2 3 2
#>   ..$ nb_backbone_species    : num [1:7] 17 13 9 3 1 3 2
#>   ..$ nb_missing_species     : num [1:7] 4 6 3 2 1 0 0
#>   ..$ nb_unsampled_species   : num [1:7] 4 5 2 3 0 0 0
#>   ..$ sampling_fraction      : num [1:7] 0.84 0.792 0.857 0.625 1 ...

# Expected richness and sampling fraction of each clade, examined before any grafting
attr(richness_per_clades, "richness_table")
#>   Clade_name Clade_tips total_richness_provided total_richness
#> 1   Family_A   t26, t19                      25             25
#> 2   Family_B    t8, t22                      NA             24
#> 3  Subfam_B1     t8, t7                      14             14
#> 4   Tribe_B2    t5, t12                       8              8
#> 5      Gen_2         t5                      NA              2
#> 6      Gen_1   t13, t27                      NA              3
#> 7  Tribe_A1b   t18, t19                      NA              2
#>   nb_backbone_species nb_missing_species nb_unsampled_species sampling_fraction
#> 1                  17                  4                    4         0.8400000
#> 2                  13                  6                    5         0.7916667
#> 3                   9                  3                    2         0.8571429
#> 4                   3                  2                    3         0.6250000
#> 5                   1                  1                    0         1.0000000
#> 6                   3                  0                    0         1.0000000
#> 7                   2                  0                    0         1.0000000

# Use the lists to build the fully grafted phylogenies in the recommended sequence:
# the missing taxa first, and then the missing clades
grafted_taxa <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = lists$taxa,
                                   total_richness_per_clades = richness_per_clades$taxa,
                                   seed = 1234, verbose = FALSE)
grafted_phy <- graft_missing_clades(phy = grafted_taxa, missing_taxa_list = lists$clades,
                                    total_richness_per_clades = richness_per_clades$clades,
                                    seed = 1234, verbose = FALSE)
```
