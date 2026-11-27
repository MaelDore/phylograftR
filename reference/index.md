# Package index

## Format inputs

Build or import the list of missing taxa and the total richness of
clades.

- [`import_taxonomic_table()`](https://maeldore.github.io/phylograftR/reference/import_taxonomic_table.md)
  : Import an external taxonomic table as a data frame

- [`build_missing_taxa_list_from_taxonomy()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list_from_taxonomy.md)
  :

  Build a `missing_taxa_list` from taxonomic tables

- [`build_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/build_missing_taxa_list.md)
  :

  Build a `missing_taxa_list` from an internal data.frame of missing
  taxa

- [`import_missing_taxa_list()`](https://maeldore.github.io/phylograftR/reference/import_missing_taxa_list.md)
  :

  Import an external table of missing taxa as a `missing_taxa_list`
  object

- [`convert_nodes_to_tips()`](https://maeldore.github.io/phylograftR/reference/convert_nodes_to_tips.md)
  :

  Convert the MRCA node IDs of a `missing_taxa_list` into tips

- [`build_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/build_total_richness_per_clades.md)
  :

  Build the `total_richness_per_clades` object from a table of clades

- [`import_total_richness_per_clades()`](https://maeldore.github.io/phylograftR/reference/import_total_richness_per_clades.md)
  :

  Import an external table of clades as `total_richness_per_clades`
  objects

## Prepare species attributes

Prepare stochastic maps (simmaps) used to constrain grafting.

- [`convert_contMap_to_simmap()`](https://maeldore.github.io/phylograftR/reference/convert_contMap_to_simmap.md)
  : Convert contMap(s) of a continuous trait into simmap(s) with
  discrete states
- [`convert_values_to_states()`](https://maeldore.github.io/phylograftR/reference/convert_values_to_states.md)
  : Convert continuous trait values into the states of a discretized
  trait

## Graft missing lineages

Graft missing taxa and/or missing clades on a backbone phylogeny.

- [`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md)
  : Graft missing taxa on a time-calibrated phylogeny
- [`graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
  : Graft missing clades on a time-calibrated phylogeny
- [`graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md)
  : Graft missing taxa and missing clades on a time-calibrated phylogeny

## Check grafting and extract information

Check the grafting plan before grafting and the constraints after
grafting, and extract information on the grafted phylogenies.

- [`check_grafting_plan()`](https://maeldore.github.io/phylograftR/reference/check_grafting_plan.md)
  : Check the feasibility of a grafting plan before grafting
- [`check_grafting_constraints()`](https://maeldore.github.io/phylograftR/reference/check_grafting_constraints.md)
  : Check that grafted taxa and missing clades respect their grafting
  constraints
- [`extract_grafting_info()`](https://maeldore.github.io/phylograftR/reference/extract_grafting_info.md)
  : Extract the grafting information from grafted phylogenies
- [`extract_sampling_info()`](https://maeldore.github.io/phylograftR/reference/extract_sampling_info.md)
  : Extract the richness and sampling fraction of clades from grafted
  phylogenies
- [`extract_BD_model_summary()`](https://maeldore.github.io/phylograftR/reference/extract_BD_model_summary.md)
  : Extract the birth-death models fitted during grafting

## Handle grafted phylogenies

Identify grafted elements, flatten sets of grafted phylogenies, and plot
them.

- [`identify_grafted()`](https://maeldore.github.io/phylograftR/reference/identify_grafted.md)
  : Identify the grafted tips and branches of a phylogeny
- [`flatten_grafted_phylogenies()`](https://maeldore.github.io/phylograftR/reference/flatten_grafted_phylogenies.md)
  : Flatten grafted phylogenies into a single set of phylogenies
- [`plot_grafted()`](https://maeldore.github.io/phylograftR/reference/plot_grafted.md)
  : Plot a grafted phylogeny with the grafted taxa highlighted

## Datasets

Toy datasets used to run examples and tutorials.

- [`toy_backbone`](https://maeldore.github.io/phylograftR/reference/toy_backbone.md)
  : Toy backbone phylogeny
- [`toy_contMap`](https://maeldore.github.io/phylograftR/reference/toy_contMap.md)
  : Toy contMap of a continuous trait
- [`toy_contMaps`](https://maeldore.github.io/phylograftR/reference/toy_contMaps.md)
  : Toy list of contMaps: stochastic histories of a continuous trait
- [`toy_missing_taxa_list`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list.md)
  : Toy list of missing taxa
- [`toy_missing_taxa_list_with_allowed_ranges`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_ranges.md)
  : Toy list of missing taxa, with allowed ranges
- [`toy_missing_taxa_list_with_allowed_states`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_allowed_states.md)
  : Toy list of missing taxa, with allowed states
- [`toy_missing_taxa_list_with_ranges_and_elevation`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_ranges_and_elevation.md)
  : Toy list of missing taxa, with allowed states of two traits:
  biogeographic ranges and elevation
- [`toy_missing_taxa_list_with_two_traits`](https://maeldore.github.io/phylograftR/reference/toy_missing_taxa_list_with_two_traits.md)
  : Toy list of missing taxa, with allowed states of two traits: colors
  and elevation
- [`toy_multiSimmap`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap.md)
  : Toy multiSimmap of a categorical trait: stochastic maps of the
  colors
- [`toy_multiSimmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_elevation.md)
  : Toy multiSimmap of a binary trait: stochastic maps of the elevation
- [`toy_multiSimmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_multiSimmap_ranges.md)
  : Toy multiSimmap of biogeographic ranges
- [`toy_simmap`](https://maeldore.github.io/phylograftR/reference/toy_simmap.md)
  : Toy stochastic character map of the toy backbone phylogeny
- [`toy_simmap_elevation`](https://maeldore.github.io/phylograftR/reference/toy_simmap_elevation.md)
  : Toy stochastic character map of elevation states
- [`toy_simmap_ranges`](https://maeldore.github.io/phylograftR/reference/toy_simmap_ranges.md)
  : Toy stochastic character map of biogeographic ranges of the toy
  backbone phylogeny
- [`toy_backbone_taxonomic_table`](https://maeldore.github.io/phylograftR/reference/toy_taxonomic_tables.md)
  [`toy_missing_taxa_taxonomic_table`](https://maeldore.github.io/phylograftR/reference/toy_taxonomic_tables.md)
  : Toy taxonomic tables
- [`toy_tip_data`](https://maeldore.github.io/phylograftR/reference/toy_tip_data.md)
  : Toy tip data of three species attributes
- [`toy_tip_values`](https://maeldore.github.io/phylograftR/reference/toy_tip_values.md)
  : Toy tip values of a continuous trait
- [`toy_total_richness_per_clades`](https://maeldore.github.io/phylograftR/reference/toy_total_richness_per_clades.md)
  : Toy total richness of the target clades
