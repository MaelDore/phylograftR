# Flatten grafted phylogenies into a single set of phylogenies

The grafting functions return a list with one element per stochastic map
when they are used with several maps (a `"multiSimmap"`), each element
being a phylogeny or a set of phylogenies (one per simulation). This
function flattens such an output into a single object of class
`"multiPhylo"`, which can be used directly by the functions that expect
a set of phylogenies under this class (e.g., from the packages `ape` or
`phytools`). The associated stochastic map ID, the simulation ID and the
seed of each phylogeny are kept, so the phylogenies can still be traced
back and selected.

The arguments `maps` and `simulations` can be used to subset the object
and retain only a specific sets of grafted phylogenies based on their
associated stochastic map, and/or their simulation index.

Additionally, the function can **trim out the updated stochastic maps**
stored in the `$simmap` element of each grafted function to save memory
space if those updated grafted stochastic maps are not required for
downstream analyses.

## Usage

``` r
flatten_grafted_phylogenies(
  phy_grafted,
  maps = NULL,
  simulations = NULL,
  keep_simmap = TRUE
)
```

## Arguments

- phy_grafted:

  The output of
  [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  or
  [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md):
  a list with one element per stochastic map (each a `"phylo"` or a
  `"multiPhylo"`), a `"multiPhylo"`, or a `"phylo"`. A `"multiPhylo"` or
  a `"phylo"` object is already flat, but the function can still be used
  to select a subset of phylogenies and remove their associated simmaps
  if requested.

- maps:

  Character vector. Names of the maps to keep (the names of the elements
  of the list, as in the column `map` of the attribute `"phylogenies"`).
  Default is `NULL`, which keeps all the maps. **The function stops if
  some names are not found.**

- simulations:

  Integer vector. Indices of the simulations to keep (the numbers in the
  names `"sim_k"` of the phylogenies, which are the same in each map:
  `1` for the first simulation of each map, etc.). Default is `NULL`,
  which keeps all the simulations.

- keep_simmap:

  Logical. Whether to keep the element `$simmap` (the updated simmap) of
  the phylogenies, if they have one. Set to `FALSE` **to save memory**
  when the stochastic maps are not needed in downstream analyses.
  Default is `TRUE`.

## Value

If `phy_grafted` is an object of class `"phylo"`, the output is also a
phylogeny, with its `$simmap` removed if `keep_simmap = FALSE`
(otherwise unchanged). Otherwise, the output is an object of class
`"multiPhylo"` with all the phylogenies selected, and the following
attributes:

- `"phylogenies"`: data.frame with one row per phylogeny, in the order
  of the object, with the columns `name`, `map` (`NA` if no map was
  used), `simulation` (index of the simulation) and `seed` (`NA` if
  unknown).

- `"seeds"`: named integer vector with the seed of each phylogeny
  (`NULL` if the seeds are unknown, i.e., if the grafting was run with
  `seed = NULL`).

- `"failed"`: data.frame with the simulations that failed in the maps
  kept, with the columns `map`, `seed` and `reason` (0 row if none). Not
  filtered by `simulations`.

- `"failed_maps"`: data.frame with the maps dropped from the output of
  the grafting function because none of their simulations could be
  obtained, with the columns `map` and `reason` (0 row if none).

## Details

The phylogenies are returned in the order of the input (map by map, then
simulation by simulation), with their names as in the input: the names
of the maps for a list of phylogenies (one simulation per map), and
`"<map>_sim_k"` for a list of sets of phylogenies (several simulations
per map). **If two phylogenies have the same name, the names are made
unique with a warning.**

Each phylogeny keeps its attributes `"map"` and `"seed"`. When `maps`
and `simulations` are used, **the selection is applied to all the
phylogenies together**: the phylogenies kept are the ones from the maps
in `maps` and with the simulation index in `simulations`.

## See also

Other functions from phylograftR:

- [`phylograftR::graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  to graft taxa based on target clades and excluded subclades using the
  list.

- [`phylograftR::graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  to graft groups of missing taxa as sister clades using the list.

## Author

Maël Doré

## Examples

``` r
# --- Prepare data --- #

# Load backbone phylogeny, associated simmap, and missing taxa list
data(toy_backbone)
data(toy_simmap)
data(toy_missing_taxa_list_with_allowed_states)

# Create a toy backbone phylogeny, with a set of 3 stochastic maps of a trait
tip_states <- phytools::getStates(toy_simmap, type = "tips")
had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
set.seed(118)
toy_multiSimmap <- phytools::make.simmap(toy_backbone, tip_states,
                                         model = "ER", nsim = 3, message = FALSE)
if (had_seed)
{ assign(".Random.seed", old_seed, envir = globalenv())
} else { rm(".Random.seed", envir = globalenv()) }

# --- Graft missing taxa --- #

# Graft with the 3 maps and 2 simulations per map as constraints:
# the output is a list with one element per map hosting 2 simulations each.
grafted_list <- graft_missing_taxa(phy = toy_backbone,
                              missing_taxa_list = toy_missing_taxa_list_with_allowed_states,
                              simmap = toy_multiSimmap, allowed_states = c("blue", "green"),
                              method = "BD_clades", nb_simulations = 2,
                              verbose = FALSE, seed = 1234)
names(grafted_list)
#> [1] "map_1" "map_2" "map_3"
str(grafted_list, 2)
#> List of 3
#>  $ map_1:Class "multiPhylo"
#> List of 2
#>   ..$ map_1_sim_1:List of 7
#>   .. ..- attr(*, "class")= chr "phylo"
#>   .. ..- attr(*, "order")= chr "cladewise"
#>   .. ..- attr(*, "seed")= int 1234
#>   .. ..- attr(*, "map")= chr "map_1"
#>   ..$ map_1_sim_2:List of 7
#>   .. ..- attr(*, "class")= chr "phylo"
#>   .. ..- attr(*, "order")= chr "cladewise"
#>   .. ..- attr(*, "seed")= int 488349520
#>   .. ..- attr(*, "map")= chr "map_1"
#>   ..- attr(*, "seeds")= Named num [1:2] 1.23e+03 4.88e+08
#>   .. ..- attr(*, "names")= chr [1:2] "map_1_sim_1" "map_1_sim_2"
#>   ..- attr(*, "failed")='data.frame':    0 obs. of  2 variables:
#>  $ map_2:Class "multiPhylo"
#> List of 2
#>   ..$ map_2_sim_1:List of 7
#>   .. ..- attr(*, "class")= chr "phylo"
#>   .. ..- attr(*, "order")= chr "cladewise"
#>   .. ..- attr(*, "seed")= int 1234
#>   .. ..- attr(*, "map")= chr "map_2"
#>   ..$ map_2_sim_2:List of 7
#>   .. ..- attr(*, "class")= chr "phylo"
#>   .. ..- attr(*, "order")= chr "cladewise"
#>   .. ..- attr(*, "seed")= int 488349520
#>   .. ..- attr(*, "map")= chr "map_2"
#>   ..- attr(*, "seeds")= Named num [1:2] 1.23e+03 4.88e+08
#>   .. ..- attr(*, "names")= chr [1:2] "map_2_sim_1" "map_2_sim_2"
#>   ..- attr(*, "failed")='data.frame':    0 obs. of  2 variables:
#>  $ map_3:Class "multiPhylo"
#> List of 2
#>   ..$ map_3_sim_1:List of 7
#>   .. ..- attr(*, "class")= chr "phylo"
#>   .. ..- attr(*, "order")= chr "cladewise"
#>   .. ..- attr(*, "seed")= int 1234
#>   .. ..- attr(*, "map")= chr "map_3"
#>   ..$ map_3_sim_2:List of 7
#>   .. ..- attr(*, "class")= chr "phylo"
#>   .. ..- attr(*, "order")= chr "cladewise"
#>   .. ..- attr(*, "seed")= int 488349520
#>   .. ..- attr(*, "map")= chr "map_3"
#>   ..- attr(*, "seeds")= Named num [1:2] 1.23e+03 4.88e+08
#>   .. ..- attr(*, "names")= chr [1:2] "map_3_sim_1" "map_3_sim_2"
#>   ..- attr(*, "failed")='data.frame':    0 obs. of  2 variables:
#>  - attr(*, "failed_maps")='data.frame':  0 obs. of  2 variables:
#>   ..$ map   : chr(0) 
#>   ..$ reason: chr(0) 

# --- Flatten the output --- #

# Flatten the output into a single multiPhylo with the 6 grafted phylogenies
grafted_phys <- flatten_grafted_phylogenies(grafted_list)
names(grafted_phys)
#> [1] "map_1_sim_1" "map_1_sim_2" "map_2_sim_1" "map_2_sim_2" "map_3_sim_1"
#> [6] "map_3_sim_2"
# The map, simulation and seed of each phylogeny are stored as an attribute
attr(grafted_phys, "phylogenies")
#>          name   map simulation      seed
#> 1 map_1_sim_1 map_1          1      1234
#> 2 map_1_sim_2 map_1          2 488349520
#> 3 map_2_sim_1 map_2          1      1234
#> 4 map_2_sim_2 map_2          2 488349520
#> 5 map_3_sim_1 map_3          1      1234
#> 6 map_3_sim_2 map_3          2 488349520

# --- Flatten the output, subset, and remove simmaps all-at-once --- #

# Keep the first simulation of two maps only, without the simmaps (to save memory)
grafted_subset <- flatten_grafted_phylogenies(grafted_list, maps = c("map_1", "map_3"),
                                              simulations = 1, keep_simmap = FALSE)
attr(grafted_subset, "phylogenies")
#>          name   map simulation seed
#> 1 map_1_sim_1 map_1          1 1234
#> 2 map_3_sim_1 map_3          1 1234

# Use the set of phylogenies with functions that expect a multiPhylo,
# e.g., to compute the tree heights.
sapply(grafted_phys, function (x) max(ape::node.depth.edgelength(x)))
#> map_1_sim_1 map_1_sim_2 map_2_sim_1 map_2_sim_2 map_3_sim_1 map_3_sim_2 
#>          50          50          50          50          50          50 

```
