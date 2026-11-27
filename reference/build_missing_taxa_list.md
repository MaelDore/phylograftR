# Build a `missing_taxa_list` from an internal data.frame of missing taxa

Build a valid `missing_taxa_list` from a data.frame describing the taxa
to graft (one row per taxon, e.g., read from a .csv file), for
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
or
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md).
Lists of tips (or node IDs, see Details) are written as values of the
data.frame, and the function parses them and names every item after the
taxa (as required by the grafting functions), and runs all the checks of
validity against the backbone phylogeny.

The same data.frame can describe **both types of taxa**: taxa to graft
within a target clade
([`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)),
and taxa to graft as missing clades sister to a clade
([`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).
The function then returns **one list per type**, built from the same
backbone, as done by
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
(see 'Taxa and missing clades in the same data.frame' in Details).

Use
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
to read the table directly from an external file.

Clades are defined by tips or node IDs, which allows for fine-tuning the
constraints of each taxon (see 'Fine-tuning the clades' in Details). For
large phylogenies with neatly defined taxonomic levels, taxonomic tables
are likely easier to prepare: see
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
and
[`phylograftR::import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md).

## Usage

``` r
build_missing_taxa_list(data, phy, verbose = TRUE)
```

## Arguments

- data:

  A data frame, with one row per taxon to graft. Its columns are the
  items of `missing_taxa_list` as in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  or
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md).
  See Details for a full description. All the columns are read as text.
  Blank cells and `NA` are missing values. A column `Type` (`"taxa"` or
  `"clades"`), or both columns `Target_clade` and `Sister_clade`, allow
  for describing both types of taxa in the same data.frame.

- phy:

  An object of class `"phylo"`. The backbone phylogeny on which the taxa
  will be grafted. It is used to check the tips and clades, and to
  convert node IDs into tips.

- verbose:

  Logical. Whether to print the conversions of node IDs into tips (if
  any). Default is `TRUE`.

## Value

A `missing_taxa_list`: a list, ready to be used in the grafting
functions, with the following items. If the table describes both types
of taxa, a list with two `missing_taxa_list`: `$taxa` (with
`$Target_clade`, for
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md))
and `$clades` (with `$Sister_clade`, for
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).

- `$Taxon_name` Character vector.

- `$Target_clade` (or `$Sister_clade`) List of character vectors of
  tips.

- `$Target_clade_name` (or `$Sister_clade_name`) Only if the column was
  provided. Character vector of labels (`NA` if none). Not used for the
  grafting, but checked against the names of
  `total_richness_per_clades`, if provided (see the grafting functions).

- `$include_target_clade_stem` Only if the column was provided. Logical
  vector (`NA` uses the global argument).

- `$Excluded_subclades` Only if the column was provided. List (one
  element per taxon) of lists of character vectors of tips. `NULL` for
  taxa with no excluded subclades. The subclades are named after their
  labels if `Excluded_subclades_names` was provided.

- `$exclude_subclades_stem` Only if the column was provided. List (one
  element per taxon) of logical vectors (one value per excluded
  subclade), or `NULL`.

- `$allowed_states` Only if the column was provided. List (one element
  per taxon) of character vectors of states, or `NULL` for taxa without
  a taxon-level value (the global `allowed_states` argument is used).

Every item is named after the taxa (`$Taxon_name`), in the same order.

## Details

The columns of `data` are the items of `missing_taxa_list` as in
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
or
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md).
Column names must be written exactly as below. Tips are separated by
commas (`,`), and the subclades by semicolons (`;`). Spaces around
separators are ignored. Values are `TRUE` or `FALSE` (also `T`, `F`,
`yes`, `no`, `1`, `0`), or `NA`.

- `Taxon_name` (Required): Names of the taxa to graft. **Names must be
  unique and absent from `phy$tip.label`.**

- `Type`: (Optional) `"taxa"` for taxa to graft within a target clade,
  or `"clades"` for taxa to graft as missing clades sister to a clade
  (case is ignored). Only needed to describe both types in the same
  table with a single clade column (see 'Taxa and missing clades in the
  same table' below).

- `Target_clade` and/or `Sister_clade` (Required): Tip labels defining a
  **target clade** per taxon, i.e., the clade where the taxon is known
  to belong (e.g., its genus); or a **sister clade**, i.e., the clade on
  whose stem branch the taxon is grafted, as its sister. A pair of tips
  defines the clade as their MRCA node (e.g., `"t26, t18"`). A single
  tip defines the "clade" as the terminal branch of this tip (e.g.,
  `"t15"`). A node ID of `phy` (digits only, e.g., `"34"`) also defines
  the clade, see 'Node IDs' below. Use `Sister_clade` to prepare a list
  for
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md):
  the columns specific to
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  below are then ignored, with a warning.

- `include_target_clade_stem`: (Optional) Logical values. For each
  taxon, whether the stem branch of its target clade is available for
  grafting. `NA` values and blank cells use the global
  `include_target_clade_stem` argument of
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md).

- `Excluded_subclades`: (Optional) For each taxon, tip labels defining
  the **excluded subclades** of its target clade where the taxon is
  known not to belong. Each subclade is defined by a pair of tips, a
  single tip, or a node ID. Ex: `"t10, t19; 50; t7"` for three
  subclades, the two first being clades with multiple species defined by
  a pair of tips (`"t10, t19"`), and an internal node ID (`"50"`). The
  last one is a single tip define by its tip label (`"t7"`). Use blank
  cells for taxa without excluded subclades.

- `exclude_subclades_stem`: (Optional) Logical values. For each taxon,
  whether the stem of each of its excluded subclades is also excluded,
  with one value per subclade of `Excluded_subclades` separated by
  semicolons, or a single value used for all of them. Ex:
  `"FALSE; TRUE; TRUE"`. `NA` values and blank cells use the global
  `exclude_subclades_stem` argument of
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md).

- `allowed_states`: (Optional) For each taxon, the states of a `simmap`
  in which it can be grafted (see the `simmap` argument of
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  or
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)),
  separated by commas. It is accepted with both `Target_clade` and
  `Sister_clade`. Ex: `"blue, green"`. Blank cells and `NA` have no
  taxon-level value, so the global `allowed_states` argument of the
  grafting function is used for the taxon. **States are matched exactly
  with the states of the `simmap`** (e.g., case matters) when grafting,
  and so cannot contain commas. If the `simmap` argument is a list of
  the maps of several traits, give the states of each trait separated by
  semicolons, in the order of the traits (e.g., `"A1, A2;; C3"` for
  three traits: no states for the second trait, so the global
  `allowed_states` of the trait, or all its states, are used), or name
  the traits (e.g., `"color=A1,A2; elevation=C3"`). The names of the
  states are not checked here as the `simmap` is not provided: the
  grafting function checks them.

Optional columns that do not match items in `missing_taxa_list` can be
used to provide names of clades:

- `Target_clade_name` or `Sister_clade_name` (Optional): Label of the
  target (or sister) clade. It is only used to check that the clades are
  defined consistently (see 'Labels of the clades' section below).

- `Excluded_subclades_names` (Optional): Labels of the excluded
  subclades, separated by semicolons, in the same order as in
  `Excluded_subclades`. Ex: `"group_1; group_2; tip_t7"`. **There must
  be one label per subclade.** If provided, they are used as names of
  the subclades in the `Excluded_subclades` item of the list, and in the
  messages of the grafting functions.

Other columns are ignored, with a warning.

### Taxa and missing clades in the same data.frame

**A data.frame can describe both types of taxa**, so a single call
prepares the lists for
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
and
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md),
on the same backbone and with the same clades written with the same
tips. The type of each taxon is given by:

- The column `Type`, with `"taxa"` (target clade) or `"clades"` (sister
  clade): the clade is then read from `Target_clade` or `Sister_clade`,
  respectively.

- Otherwise, without `Type`, by the clade column that is filled in its
  row, when the data.frame has both `Target_clade` and `Sister_clade`
  columns. **Each row must fill exactly one of them** (the function
  stops if both or none are filled).

The columns specific to a type are only read for the taxa of this type
(e.g., `Excluded_subclades` for taxa with a target clade), and are
ignored for the others, with a warning if they are filled.
`allowed_states` is read for both types. A clade written in the column
of the other type than `Type` is ignored with a warning.

The function returns:

- **A single `missing_taxa_list`** if all the taxa have the same type,
  whether `Type` is provided or not.

- **A list of two `missing_taxa_list`**, named `$taxa` and `$clades`, if
  both types are present. Each one is a standard `missing_taxa_list`,
  and the taxa of both lists are named differently (names must be unique
  across the table). Use `$taxa` in
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  and `$clades` in
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md).

Clade labels are also checked across both lists (the same clade has the
same label as a target clade and as a sister clade), so that they can be
matched with the names of `total_richness_per_clades`.

**Recommended use:** provide the whole list (`$taxa` and `$clades`) to
[`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md),
which grafts the missing taxa first and then the missing clades, with a
single set of inputs. The two steps can also be chained by hand: graft
the missing taxa with
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md),
then use the grafted phylogeny as the backbone to graft the missing
clades with
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md).
The clades are defined by tips of the backbone, so they stay valid even
when the grafted phylogeny is used for the second step.

### Node IDs

Clades (`Target_clade`, `Sister_clade`, and each subclade of
`Excluded_subclades`) are typically defined by (pairs of) tips. However,
the function also accepts clades defined by the ID of a node of `phy`
instead of tips, in the same table, and even in the same cell (e.g.,
`"50; t5, t12; t7"` for three excluded subclades defined respectively by
an internal node, a pair of tips, a terminal tip label). Node IDs are
integers: `1` to N, with N being the number of tips (`ape::Ntip(phy)`)
correspond to terminal nodes, and thus define single-tip clades.
Following IDs are the internal nodes, as shown by
[`ape::nodelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html), and
correspond to clades with multiple descending species.

A value made only of digits is a node ID, so **tip labels must include
at least one other character than digits** to be distinguished. A clade
is defined by tips or by one node ID, not both, and one node ID defines
the clade alone. Node IDs are converted into tips before any other
check, as done by
[`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md),
and the conversions are printed if `verbose = TRUE`.

The warnings and errors that follow, and the grafting functions,
therefore use the tips that are traced back to your node IDs in the
messages. If a clade is already defined by tips in the table, a node ID
of the same clade is converted into the same tips. Convert
`total_richness_per_clades` (if you use it) at the same time as the
list, with
[`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md),
so the tips of the clades match between both objects. Alternatively,
build it from a table of clades with
[`phylograftR::build_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/build_total_richness_per_clades.md),
which matches the clades with those of the `missing_taxa_list`.

### Fine-tuning the clades

Defining clades with tips or node IDs gives full control on the
constraints of each taxon. **Any node of `phy` can be used** as a
`Target_clade`, as a `Sister_clade`, or as an excluded subclade, simply
by writing its node ID (e.g., `"34"`), or two tips whose MRCA is the
node (e.g., `"t8, t7"`). The node does not need to match a taxonomic
group: it can be a clade without a name, a clade that groups several
taxonomic units (e.g., two neighboring genera), a part of a taxonomic
group only, or a clade you identified from the phylogeny itself or from
another source (e.g., a molecular study). Clades can be nested at any
depth, and each taxon has its own set of excluded subclades, with or
without their stem.

This is the best route to fine-tune the grafting of a limited number of
taxa, but the clades have to be written for each taxon. When the
backbone is large, and taxonomic levels (e.g., Family, Subfamily, Tribe,
Genus) are defined for its tips and for the taxa to graft, **taxonomic
tables are the go-to option**: clades are retrieved from their taxonomic
names, which ease the attribution of clades across many taxa. See
[`phylograftR::import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md)
and
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md).

Both routes return the same `missing_taxa_list`, so they can be
combined: build a list from taxonomic tables, and then adjust the clades
of a few taxa with node IDs or tips.

### Labels of the clades

Labels (`Target_clade_name`, `Sister_clade_name`, and
`Excluded_subclades_names`) are optional but, when used, they identify
the clades across all the rows of the table: **the same clade must
always have the same label, and different clades must have different
labels.** A clade is identified by the MRCA node of its tips in `phy`,
so the same clade can be defined by different pairs of tips, as long as
they share this MRCA node. The function stops if a label is used for
different clades, or if a clade is given different labels, and lists the
labels and tips in conflict. This catches typos in the tips or in the
labels. Subclades without a label are not checked.

### Validity checks

The function **stops with an informative message** if:

- Required columns are missing, or one `Taxon_name` is missing,
  duplicated, or already in `phy`.

- A clade has no tips, or tips are repeated, or are absent from `phy`.

- The number of values in `exclude_subclades_stem`, or of labels in
  `Excluded_subclades_names`, does not match the number of excluded
  subclades.

- Labels are inconsistent (see above).

- A cell of `allowed_states` has an empty state (check the commas) or
  repeated states.

- An excluded subclade is identical to, or contains, the target clade of
  its taxon (no branch would remain). A warning is issued if it is
  outside of the target clade (the exclusion has no effect).

The messages give the names of the problematic taxa and the numbers of
the rows in the table.

## See also

Other functions from phylograftR:

- [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  to read the table directly from an external file.

- [`phylograftR::import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md)
  to read the taxonomic tables from external files.

- [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
  to build a `missing_taxa_list` from taxonomic tables, with clades
  defined by taxonomic names instead of tips or node IDs.

- [`phylograftR::convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md)
  to convert the node IDs of a list built 'by hand' into pairs of tips
  defining their MRCA.

- [`phylograftR::build_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/build_total_richness_per_clades.md)
  and
  [`phylograftR::import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md)
  to prepare `total_richness_per_clades` for the unsampled species.

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades using the
  list.

- [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft groups of missing taxa as sister clades using the list.

- [`phylograftR::graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)
  to graft missing taxa and missing clades in one run, with a single set
  of inputs.

## Author

Maël Doré

## Examples

``` r
# Load a table of missing taxa as a data.frame in R
file_path <- system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")
missing_taxa_table <- read.csv(file_path)

# Convert it into a 'missing_taxa_list', checking it against the backbone phylogeny
missing_taxa_list <- build_missing_taxa_list(data = missing_taxa_table, phy = toy_backbone)
# Items are named after the taxa, tips are split, and excluded subclades are named
str(missing_taxa_list, 2)
#> List of 6
#>  $ Taxon_name               : chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Target_clade             :List of 7
#>   ..$ New_sp1: chr [1:2] "t26" "t18"
#>   ..$ New_sp2: chr [1:2] "t26" "t18"
#>   ..$ New_sp3: chr [1:2] "t26" "t18"
#>   ..$ New_sp4: chr [1:2] "t8" "t22"
#>   ..$ New_sp5: chr [1:2] "t8" "t7"
#>   ..$ New_sp6: chr "t5"
#>   ..$ New_sp7: chr [1:2] "t5" "t12"
#>  $ Target_clade_name        : Named chr [1:7] "Family_A" "Family_A" "Family_A" "Family_B" ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ include_target_clade_stem: Named logi [1:7] TRUE TRUE FALSE NA TRUE NA ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Excluded_subclades       :List of 7
#>   ..$ New_sp1:List of 2
#>   ..$ New_sp2:List of 2
#>   ..$ New_sp3: NULL
#>   ..$ New_sp4: NULL
#>   ..$ New_sp5:List of 1
#>   ..$ New_sp6: NULL
#>   ..$ New_sp7: NULL
#>  $ exclude_subclades_stem   :List of 7
#>   ..$ New_sp1: Named logi [1:2] TRUE FALSE
#>   .. ..- attr(*, "names")= chr [1:2] "Subfam_A1" "Tribe_A2a"
#>   ..$ New_sp2: Named logi [1:2] FALSE TRUE
#>   .. ..- attr(*, "names")= chr [1:2] "Tribe_A2a" "Gen_3"
#>   ..$ New_sp3: NULL
#>   ..$ New_sp4: NULL
#>   ..$ New_sp5: NULL
#>   ..$ New_sp6: NULL
#>   ..$ New_sp7: NULL

# A table can also be built in R and then converted to a valid 'missing_taxa_list' object
my_table <- data.frame(Taxon_name = c("M1", "M2"),
                       Target_clade = c("t26, t18", "t5"),
                       Excluded_subclades = c("t10, t19; t26, t4", NA),
                       Excluded_subclades_names = c("group_1; group_2", NA))
my_missing_taxa_list <- build_missing_taxa_list(data = my_table, phy = toy_backbone)
str(my_missing_taxa_list, 2)
#> List of 3
#>  $ Taxon_name        : chr [1:2] "M1" "M2"
#>  $ Target_clade      :List of 2
#>   ..$ M1: chr [1:2] "t26" "t18"
#>   ..$ M2: chr "t5"
#>  $ Excluded_subclades:List of 2
#>   ..$ M1:List of 2
#>   ..$ M2: NULL

# Clades can also be defined by node IDs (digits only), alone or mixed with tips.
# Node IDs are converted into tips, and the conversions are printed.
ape::plot.phylo(toy_backbone, cex = 0.7, no.margin = TRUE)
ape::nodelabels(cex = 0.6, frame = "none", col = "blue")

my_table_nodes <- data.frame(Taxon_name = c("M1", "M2"),
                             Target_clade = c("34", "t5, t12"),
                             Excluded_subclades = c("50", NA))
str(build_missing_taxa_list(data = my_table_nodes, phy = toy_backbone), 2)
#> 2 node ID(s) converted into tips (a tip node gives its tip, an internal node gives one pair of tips with this MRCA node):
#>   - node 34 (internal node, 9 tips) -> t8, t7   [Target_clade: M1]
#>   - node 50 (internal node, 3 tips) -> t13, t27   [Excluded_subclades: M1]
#> List of 3
#>  $ Taxon_name        : chr [1:2] "M1" "M2"
#>  $ Target_clade      :List of 2
#>   ..$ M1: chr [1:2] "t8" "t7"
#>   ..$ M2: chr [1:2] "t5" "t12"
#>  $ Excluded_subclades:List of 2
#>   ..$ M1:List of 1
#>   ..$ M2: NULL

# A table can describe both types of taxa: here, with a column 'Type'.
# Taxa to graft within a target clade and taxa to graft as missing clades
# are returned as two lists, for the two grafting functions.
my_table_both <- data.frame(Taxon_name = c("M1", "M2", "M3"),
                            Type = c("taxa", "taxa", "clades"),
                            Target_clade = c("t26, t18", "t5, t12", NA),
                            Sister_clade = c(NA, NA, "t10, t19"),
                            Excluded_subclades = c("t26, t4", NA, NA))
lists <- build_missing_taxa_list(data = my_table_both, phy = toy_backbone)
str(lists, 2)
#> List of 2
#>  $ taxa  :List of 3
#>   ..$ Taxon_name        : chr [1:2] "M1" "M2"
#>   ..$ Target_clade      :List of 2
#>   ..$ Excluded_subclades:List of 2
#>  $ clades:List of 2
#>   ..$ Taxon_name  : chr "M3"
#>   ..$ Sister_clade:List of 1

# Without 'Type', each row fills the clade column of its type
my_table_both$Type <- NULL
identical(build_missing_taxa_list(data = my_table_both, phy = toy_backbone), lists)
#> [1] TRUE
```
