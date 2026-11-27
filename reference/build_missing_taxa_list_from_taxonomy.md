# Build a `missing_taxa_list` from taxonomic tables

Build a ready-to-be-used `missing_taxa_list` object based on taxonomic
information from the backbone phylogeny and the taxa to graft. The
`missing_taxa_list` can then be used as input for
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
or
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md).

Target clades, sister clades and excluded subclades are identified with
taxonomic names (e.g., `"Homo"`) or with taxonomic levels (e.g.,
`"Genus"`) instead of tips or node IDs. The function retrieves the MRCA
node of each taxonomic group in the backbone phylogeny, checks that the
group is monophyletic, and writes each clade as a pair of tips, as
expected by
[`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md).

Taxonomic tables are the go-to option for large phylogenies with neatly
defined taxonomic levels (see 'When to use taxonomic tables' in
Details). Use
[`phylograftR::import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md)
to read the tables directly from files.

For a finer level of control on the clades to use as target clades,
sister clades, and/or excluded subclades, you can use alternative
functions that define those clades based on any node IDs or pairs of
tips used to identify MRCA: see
[`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
and
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md).

## Usage

``` r
build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table,
  missing_taxa_taxonomic_table,
  phy,
  non_monophyletic = c("exclude", "mrca", "error"),
  verbose = TRUE
)
```

## Arguments

- backbone_taxonomic_table:

  Data frame. Taxonomy of the tips of the backbone phylogeny, with one
  row per tip. The first column must be named `Taxa` and holds the tip
  labels of `phy`. All the other columns are **taxonomic levels** (e.g.,
  `Family`, `Subfamily`, `Tribe`, `Genus`, but they can be anything),
  usually from the broadest to the narrowest. The values of a level are
  the names of the **taxonomic groups** (e.g., `"Homo"`). Levels do not
  have to be informed for all the tips: blank cells and `NA` are missing
  values. See Details for the constraints on the names of the levels and
  groups.

- missing_taxa_taxonomic_table:

  Data frame. Description of the taxa to graft, with one row per taxon.
  Columns `Taxa` (names of the taxa to graft) and `Type` (`"taxa"` or
  `"clades"`) are required, as well as the column of the clade to use
  for grafting (`Target_clade` or `Sister_clade`). Entries can be either
  the name of the taxonomic group (e.g., `"Homo"`) or the taxonomic
  level of the taxa to use (e.g., `Genus`). Next columns provide the
  taxonomic levels of the taxa to graft using the same names as in
  `backbone_taxonomic_table` (e.g., `Family`, `Subfamily`, `Tribe`,
  `Genus`, but they can be anything). Other optional columns provide
  additional information on constraints to use for grafting:
  `Excluded_subclades`, `include_target_clade_stem`,
  `exclude_subclades_stem`, and `allowed_states`, as defined in
  `missing_taxa_list` argument of the grafting functions
  ([`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  or
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)).
  See 'Missing taxa taxonomic table' in Details section for an extended
  description.

- phy:

  An object of class `"phylo"`. The backbone phylogeny on which the taxa
  will be grafted.

- non_monophyletic:

  Character string. What to do when a taxonomic group is not
  monophyletic in `phy`, i.e., when the clade defined by the MRCA of its
  tips also holds tips of other groups of the same level. Must be one of
  `"exclude"` (default), `"mrca"` or `"error"`. See 'Non-monophyletic
  groups' in Details.

- verbose:

  Logical. Whether to print a summary of the conversion. Default is
  `TRUE`.

## Value

A `missing_taxa_list` object as returned by
[`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md):
a list ready-to-be-used as input of
[`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
if all the taxa have `Type = "taxa"`, or of
[`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
if they all have `Type = "clades"`. If both types are present, a list
with two elements, `taxa` and `clades`, each holding the
`missing_taxa_list` of one type.

The attribute `"taxonomic_groups"` is a data frame with one row per
taxonomic group used, with the columns `Group`, `Level`, `n_members`
(number of tips of the group), `clade_node` (MRCA node ID in `phy`),
`clade_tips` (pair of tips defining the clade), `n_tips_clade` (number
of tips of the clade), `n_foreign` (tips of other groups in the clade),
`n_unannotated` (tips without value at this level in the clade), and
`monophyletic`.

## Details

**When to use taxonomic tables**

**Taxonomic tables are the go-to option for large phylogenies with
neatly defined taxonomic levels** (e.g., Family, Subfamily, Tribe,
Genus). The taxonomy of the backbone is provided once, and each taxon to
graft only needs its own taxonomic levels: the constraints (target
clade, sister clade, excluded subclades) are written with the names of
the groups or of the taxonomic levels to use, and the clades are
retrieved from the phylogeny. The same group is therefore always the
same clade, whatever the number of taxa, which avoids the typos and
inconsistencies of clades written taxon by taxon, and the tables can be
taken as such from a taxonomic database.

The limit is that only the taxonomic groups of the backbone table can be
used as clades. To fine-tune the grafting of some taxa with any other
node of the phylogeny (e.g., a clade without a name, a part of a group,
or a group of several genera), **define their clades by tips or node
IDs** with
[`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
or
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
(see 'Fine-tuning the clades' in their Details sections).

Both routes return the same `missing_taxa_list`, so a list built from
taxonomic tables can be adjusted afterwards.

**Backbone taxonomic table**

Names of taxonomic levels (the column names; e.g., `"Genus"`) and names
of taxonomic groups (their values; e.g., `"Homo"`) must be unambiguous:

- A level name cannot be the name of a tip, of a taxon to graft, or of a
  group.

- A group name belongs to a single level: the same name cannot be found
  in two columns (e.g., a Family cannot be named the same as a Genus).

- Level names and group names cannot contain `";"`. Level names cannot
  be `"Taxa"`, `"Type"` or the names of the other columns of
  `missing_taxa_taxonomic_table`.

- All the tips of `backbone_taxonomic_table` must be in `phy`. Tips of
  `phy` absent from the table are allowed (with a warning).

**Missing taxa taxonomic table**

- `Taxa` (Required): Names of the taxa to graft. **Names must be unique
  and absent from `phy$tip.label`.**

- `Type` (Required): `"taxa"` for taxa grafted with
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  within a target clade, or `"clades"` for taxa grafted with
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  as a sister of a clade.

- Taxonomic levels (Optional): Same columns as in
  `backbone_taxonomic_table`. They are used to retrieve the taxonomic
  group of the taxon when `Target_clade`, `Sister_clade`, and/or
  `Excluded_subclades` (see below) use a taxonomic level (e.g.,
  `"Genus"`) rather than directly a taxonomic group (e.g., `"Homo"`) as
  identification. Columns that are not levels of the backbone, nor one
  of the columns listed here, are ignored with a warning.

- `Target_clade` (Required for `Type = "taxa"`) or `Sister_clade`
  (Required for `Type = "clades"`): A single taxonomic group, written
  either as the name of the group (e.g., `"Homo"`), or as the name of a
  taxonomic level (e.g., `"Genus"`). If a taxonomic level is provided,
  the group will be identified as recorded in the column of the same
  name in `missing_taxa_taxonomic_table` (e.g., `"Genus"` with
  `Genus = "Camponotus"` gives the clade of the genus `Camponotus`).
  **It is an error if the group is absent from the backbone, or if the
  taxon has no value at the requested level.**

- `Excluded_subclades` (Optional, `Type = "taxa"`): Subclades of the
  target clade where the taxon is known not to belong, separated by
  semicolons (e.g., `"Pheidole; Camponotus"`). Each one is a group or a
  level name, as in `Target_clade`.

- `$include_target_clade_stem` (Optional, `Type = "taxa"`): Logical.
  Whether the stem branch of the target clade is a valid placement to
  graft this taxa (`NA` uses the global argument set in the grafting
  function).

- `$exclude_subclades_stem` (Optional, `Type = "taxa"`): Logical.
  Whether the stem branch of the excluded clade(s) is/are a excluded
  from grafting too. Separated by semicolons if multiple
  `Excluded_subclades` are provided (`NA` uses the global argument set
  in the grafting function).

- `$allowed_states` (Optional): For each taxon, the states of a `simmap`
  in which it can be grafted (see the `simmap` argument of
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  or
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)),
  separated by commas. If multiple traits are used, each set of states
  must be separated by semicolons (e.g., `"blue, red; low"`). It is
  accepted with both `Target_clade` and `Sister_clade`. (`NA` uses the
  global argument set in the grafting function).

**Non-monophyletic groups**

What to do when a taxonomic group is not monophyletic in `phy`, i.e.,
when the clade defined by the MRCA of its tips also holds tips of other
groups of the same level?

- `non_monophyletic = "exclude"` (default): For target clades,
  **intruders** (i.e., tips from other groups descending from the MRCA)
  are added to the `Excluded_subclades` of the taxon (with their stem,
  `exclude_subclades_stem = TRUE`), so the taxon is never grafted among
  intruders. Intruders nested within an excluded subclade provided by
  the user are not added. Intruders cannot be removed from a sister
  clade or from a group of `Excluded_subclades`: the clade of the MRCA
  is used with a warning (as for `"mrca"`). For `Excluded_subclades`,
  all the branches descending from the MRCA are then excluded, including
  the ones of intruders.

- `non_monophyletic = "mrca"`: The clade of the MRCA is used as such,
  with a warning, including its intruders.

- `non_monophyletic = "error"`: Stop if one group used is not
  monophyletic.

## See also

Other functions from phylograftR:

- [`phylograftR::import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md)
  to read the taxonomic tables from external files.

- [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  to read a table of missing taxa from a file, with clades defined by
  tips or node IDs.

- [`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  to build a `missing_taxa_list` from a data.frame loaded in R, with
  clades defined by tips or node IDs.

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
## Load the toy backbone phylogeny and its taxonomy
data(toy_backbone)
data(toy_backbone_taxonomic_table)
head(toy_backbone_taxonomic_table)
#>   Taxa   Family Subfamily     Tribe Genus
#> 1  t26 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 2   t2 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 3  t21 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 4  t16 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 5  t11 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 6   t4 Family_A Subfam_A2 Tribe_A2a  <NA>

## Load the taxonomic table of the taxa to graft
data(toy_missing_taxa_taxonomic_table)
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

# The taxa to graft with graft_missing_taxa() (Type = "taxa") use different types of constraints:
# * New_sp1: grafted in its family (Target_clade = "Family" is the family of the taxon,
#   i.e., "Family_A"), but not in the subfamily "Subfam_A1" nor in the tribe "Tribe_A2a"
#   (Excluded_subclades).
#   The stem of the family is available (include_target_clade_stem),
#   the stem of "Subfam_A1" is excluded too,
#   but not the one of "Tribe_A2a" (exclude_subclades_stem).
# * New_sp2: same target clade, but excluded subclades are a tribe
#   and a genus with a single species ("Gen_3"), and the stem is excluded for the genus only.
# * New_sp3: same target clade, without excluded subclades and
#   without access to the stem of the family.
# * New_sp4: grafted in the family "Family_B", written directly with the name of the group.
# * New_sp5: grafted in its subfamily (Target_clade = "Subfamily", i.e., "Subfam_B1"),
#   except in the genus "Gen_1".
# * New_sp6: grafted on the terminal branch of the only species
#   of the genus "Gen_2" (a single tip).
# * New_sp7: grafted in its tribe (Target_clade = "Tribe", i.e., "Tribe_B2"),
#   without access to the stem.

# The taxa to graft with graft_missing_clades() (Type = "clades") use sister clades:
# * Genus_a_sp1 and Genus_a_sp2: belong to a new genus, and have the same sister clade ("Gen_1"),
#   so they are grafted together as a missing clade.
# * Genus_b_sp1: grafted as the sister of its own tribe (Sister_clade = "Tribe", i.e., "Tribe_A1b").

## Build the missing taxa lists associated to these taxonomic tables
missing_taxa_lists <- build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table = toy_backbone_taxonomic_table,
  missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
  phy = toy_backbone)
#> 7 taxa for graft_missing_taxa() and 3 for graft_missing_clades() built from 10 taxonomic group(s), 0 not monophyletic. See attr(, "taxonomic_groups").

# Since both types of taxa to graft were provided,
# the output contains one list per grafting function
names(missing_taxa_lists)
#> [1] "taxa"   "clades"
str(missing_taxa_lists$taxa, max.level = 1) # Missing_taxa_list for graft_missing_taxa()
#> List of 6
#>  $ Taxon_name               : chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Target_clade             :List of 7
#>  $ Target_clade_name        : Named chr [1:7] "Family_A" "Family_A" "Family_A" "Family_B" ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ include_target_clade_stem: Named logi [1:7] TRUE TRUE FALSE NA TRUE NA ...
#>   ..- attr(*, "names")= chr [1:7] "New_sp1" "New_sp2" "New_sp3" "New_sp4" ...
#>  $ Excluded_subclades       :List of 7
#>  $ exclude_subclades_stem   :List of 7
str(missing_taxa_lists$clades, max.level = 1) # Missing_taxa_list for graft_missing_clades()
#> List of 3
#>  $ Taxon_name       : chr [1:3] "Genus_a_sp1" "Genus_a_sp2" "Genus_b_sp1"
#>  $ Sister_clade     :List of 3
#>  $ Sister_clade_name: Named chr [1:3] "Gen_1" "Gen_1" "Tribe_A1b"
#>   ..- attr(*, "names")= chr [1:3] "Genus_a_sp1" "Genus_a_sp2" "Genus_b_sp1"

# Clades are written as pairs of tips, and named after their taxonomic groups
missing_taxa_lists$taxa$Target_clade
#> $New_sp1
#> [1] "t26" "t19"
#> 
#> $New_sp2
#> [1] "t26" "t19"
#> 
#> $New_sp3
#> [1] "t26" "t19"
#> 
#> $New_sp4
#> [1] "t8"  "t22"
#> 
#> $New_sp5
#> [1] "t8" "t7"
#> 
#> $New_sp6
#> [1] "t5"
#> 
#> $New_sp7
#> [1] "t5"  "t12"
#> 
missing_taxa_lists$taxa$Excluded_subclades$New_sp1
#> $Subfam_A1
#> [1] "t10" "t19"
#> 
#> $Tribe_A2a
#> [1] "t26" "t4" 
#> 
missing_taxa_lists$taxa$exclude_subclades_stem$New_sp1
#> Subfam_A1 Tribe_A2a 
#>      TRUE     FALSE 
missing_taxa_lists$clades$Sister_clade
#> $Genus_a_sp1
#> [1] "t13" "t27"
#> 
#> $Genus_a_sp2
#> [1] "t13" "t27"
#> 
#> $Genus_b_sp1
#> [1] "t18" "t19"
#> 

# The taxa to graft as clades share their sister clade, so they will be grafted together
missing_taxa_lists$clades$Sister_clade_name
#> Genus_a_sp1 Genus_a_sp2 Genus_b_sp1 
#>     "Gen_1"     "Gen_1" "Tribe_A1b" 

# Information on the groups used
attr(missing_taxa_lists, "taxonomic_groups")
#>        Group     Level n_members clade_node clade_tips n_tips_clade n_foreign
#> 1   Family_A    Family        17         32    t26,t19           17         0
#> 2   Family_B    Family        13         33     t8,t22           13         0
#> 3  Subfam_B1 Subfamily         9         34      t8,t7            9         0
#> 4      Gen_2     Genus         1         27         t5            1         0
#> 5   Tribe_B2     Tribe         3         53     t5,t12            3         0
#> 6      Gen_1     Genus         3         50    t13,t27            3         0
#> 7  Tribe_A1b     Tribe         2         40    t18,t19            2         0
#> 8  Subfam_A1 Subfamily         8         37    t10,t19            8         0
#> 9  Tribe_A2a     Tribe         6         39     t26,t4            6         0
#> 10     Gen_3     Genus         1         14        t30            1         0
#>    n_unannotated monophyletic
#> 1              0         TRUE
#> 2              0         TRUE
#> 3              0         TRUE
#> 4              0         TRUE
#> 5              0         TRUE
#> 6              0         TRUE
#> 7              0         TRUE
#> 8              0         TRUE
#> 9              0         TRUE
#> 10             0         TRUE

# The lists can be used for grafting
grafted_taxa <- graft_missing_taxa(toy_backbone, missing_taxa_lists$taxa, verbose = FALSE)
grafted_clades <- graft_missing_clades(toy_backbone, missing_taxa_lists$clades, verbose = FALSE)
```
