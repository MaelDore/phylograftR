# Import an external taxonomic table as a data frame

Read a taxonomic table from a file (e.g., a .csv file, or an Excel
sheet), check its format, and return it as a data frame ready-to-be-used
as input of
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md).

It is the mirror of
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
for taxonomic tables. Two tables are used to build a `missing_taxa_list`
from taxonomy: the taxonomy of the tips of the backbone phylogeny
(`type = "backbone"`), and the description of the taxa to graft
(`type = "missing_taxa"`). Each one is read with a call of the function.

**Taxonomic tables are the go-to option for large phylogenies with
neatly defined taxonomic levels** (see Details).

## Usage

``` r
import_taxonomic_table(
  file,
  type,
  phy = NULL,
  sep = NULL,
  sheet = 1,
  reader = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- file:

  Character string. Path to the file with the table. The format is
  chosen from the extension of the file: `".csv"` (comma-separated),
  `".tsv"`, `".tab"` or `".txt"` (tab-separated), `".xlsx"` or `".xls"`
  (Excel, requires the package `readxl`). Use `sep` or `reader` for
  other formats. All the columns are read as text. Blank cells and `NA`
  are missing values.

- type:

  Character string. The type of table to read:

  - `"backbone"` for the taxonomy of the tips of the backbone phylogeny
    (one row per tip),

  - `"missing_taxa"` for the description of the taxa to graft (one row
    per taxon). See Details for the columns of each one.

- phy:

  An object of class `"phylo"`, or `NULL`. Optional backbone phylogeny,
  used to check the table against its tips: **all the tips of a
  `"backbone"` table must be tips of `phy`, and no taxon of a
  `"missing_taxa"` table can already be a tip of `phy`.** Default is
  `NULL`, with no check against the phylogeny (it is anyway done by
  [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)).

- sep:

  Character string. Field separator of a text file. Default is `NULL`,
  to choose it from the extension of the file: `","` for `.csv`, and a
  tab for `.tsv`, `.tab` and `.txt`. Semicolons are also used inside the
  cells of the table of the missing taxa (to separate the excluded
  subclades, and their stems), so do not use them as field separator,
  and **quote the cells that contain commas** (e.g., `allowed_states`).

- sheet:

  Character string or integer. Only used for Excel files. Name or
  position of the sheet to read. Default is `1`. Both taxonomic tables
  can be stored in two sheets of the same Excel file, and read with two
  calls.

- reader:

  Function. Optional function used to read `file`, to import any other
  format (e.g., `.ods`, or an export of an online spreadsheet). It must
  take the path of the file as first argument and return a data frame.
  If provided, it is used whatever the extension of the file. Default is
  `NULL`.

- verbose:

  Logical. Whether to print a summary of the table. Default is `TRUE`.

- ...:

  Further arguments passed to the function reading the file
  ([`utils::read.table()`](https://rdrr.io/r/utils/read.table.html) for
  text files, `readxl::read_excel()` for Excel files, or `reader`).

## Value

A data frame with only character columns, ready to be used as
`backbone_taxonomic_table` (if `type = "backbone"`) or
`missing_taxa_taxonomic_table` (if `type = "missing_taxa"`) in
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md).
Blank cells are `NA`.

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

**The limit is that only the taxonomic groups of the backbone table can
be used as clades.** To fine-tune the grafting of some taxa with any
other node of the phylogeny (e.g., a clade without a name, a part of a
group, or a group of several genera), **define their clades by tips or
node IDs** with
[`phylograftR::build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
or
[`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
(see 'Fine-tuning the clades' in their Details sections).

**Both routes return the same `missing_taxa_list`**, so a list built
from taxonomic tables can be adjusted afterwards.

**Backbone table** (`type = "backbone"`)

One row per tip of the backbone phylogeny. The column `Taxa` (Required)
holds the tip labels, and **must have unique values**. All the other
columns are **taxonomic levels** (e.g., `Family`, `Subfamily`, `Tribe`,
`Genus`, but they can be anything), usually from the broadest to the
narrowest. **At least one level is required.** The values of a level are
the names of the **taxonomic groups** (e.g., `"Homo"` for a `"Genus"`).
Levels do not have to be informed for all the tips: blank cells and `NA`
are missing values.

**Missing taxa table** (`type = "missing_taxa"`)

One row per taxon to graft. The columns `Taxa` (Required; unique names
of the taxa to graft) and `Type` (Required; `"taxa"` or `"clades"`) are
required. The next columns are the taxonomic levels of the taxa to
graft, with the same names as in the backbone table, and the addtional
constraints:

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
  graft this taxa (`NA` or empty cell uses the global argument set in
  the grafting function).

- `$exclude_subclades_stem` (Optional, `Type = "taxa"`): Logical.
  Whether the stem branch of the excluded clade(s) is/are a excluded
  from grafting too. Separated by semicolons if multiple
  `Excluded_subclades` are provided (`NA` or empty cell uses the global
  argument set in the grafting function).

- `$allowed_states` (Optional): For each taxon, the states of a `simmap`
  in which it can be grafted (see the `simmap` argument of
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  or
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)),
  separated by commas. If multiple traits are used, each set of states
  must be separated by semicolons (e.g., `"blue, red; low"`). It is
  accepted with both `Target_clade` and `Sister_clade`. (`NA` or empty
  cell uses the global argument set in the grafting function).

**Checks**

The function checks that the required columns are present, that names of
columns are unique, that `Taxa` is never empty and always unique, and,
for the table of the missing taxa, that `Type` is valid. Names of
columns and values are stripped from the spaces around them. Rows with
only blank cells (e.g., the trailing rows of a spreadsheet) are dropped.

The other checks (e.g., names of levels and groups, consistency with the
backbone, and the groups used as clades) are done when the clades are
retrieved, by
[`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md).

## See also

Other functions from phylograftR:

- [`phylograftR::import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  to read a table of missing taxa from a file, with clades defined by
  tips or node IDs (for fine-tuning).

- [`phylograftR::build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
  to build the `missing_taxa_list` from the taxonomic tables.

## Author

Maël Doré

## Examples

``` r
# Paths to example tables, provided with the package
backbone_path <- system.file("extdata", "toy_backbone_taxonomic_table.csv", package = "phylograftR")
missing_taxa_path <- system.file("extdata", "toy_missing_taxa_taxonomic_table.csv",
                                 package = "phylograftR")

# Load the associated backbone phylogeny
data(toy_backbone)

# Import the taxonomy of the backbone: one row per tip, one column per taxonomic level
backbone_taxonomic_table <- import_taxonomic_table(file = backbone_path, type = "backbone",
                                                   phy = toy_backbone)
#> Backbone taxonomic table: 30 tips, with 4 taxonomic level(s): Family, Subfamily, Tribe, Genus
head(backbone_taxonomic_table)
#>   Taxa   Family Subfamily     Tribe Genus
#> 1  t26 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 2   t2 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 3  t21 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 4  t16 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 5  t11 Family_A Subfam_A2 Tribe_A2a  <NA>
#> 6   t4 Family_A Subfam_A2 Tribe_A2a  <NA>

# Import the description of the taxa to graft: one row per taxon
missing_taxa_taxonomic_table <- import_taxonomic_table(file = missing_taxa_path,
                                                       type = "missing_taxa",
                                                       phy = toy_backbone)
#> Missing taxa taxonomic table: 10 taxa to graft (7 as 'taxa', 3 as 'clades'), with 9 other column(s): Family, Subfamily, Tribe, Genus, Target_clade, Sister_clade, Excluded_subclades, include_target_clade_stem, exclude_subclades_stem
missing_taxa_taxonomic_table
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
#> 4          <NA>                <NA>                      <NA>
#> 5          <NA>               Gen_1                      TRUE
#> 6          <NA>                <NA>                      <NA>
#> 7          <NA>                <NA>                     FALSE
#> 8         Gen_1                <NA>                      <NA>
#> 9         Gen_1                <NA>                      <NA>
#> 10        Tribe                <NA>                      <NA>
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

# The tables are the datasets provided with the package
data(toy_backbone_taxonomic_table)
data(toy_missing_taxa_taxonomic_table)
identical(backbone_taxonomic_table, toy_backbone_taxonomic_table)
#> [1] TRUE
# (the logical column include_target_clade_stem of the dataset is read as text from the file)
all(mapply(function (x, y) identical(as.character(x), y),
           toy_missing_taxa_taxonomic_table, missing_taxa_taxonomic_table))
#> [1] TRUE

# Use them to build the lists of missing taxa
missing_taxa_lists <- build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table = backbone_taxonomic_table,
  missing_taxa_taxonomic_table = missing_taxa_taxonomic_table,
  phy = toy_backbone, verbose = FALSE)
names(missing_taxa_lists)
#> [1] "taxa"   "clades"

# Use them to graft the taxa
grafted_phy <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = missing_taxa_lists$taxa,
                                  verbose = FALSE, seed = 1234)
```
