## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, fig.height = 5   # Set device size at rendering time (when plots are generated)
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----plot_backbone, fig.height = 6.5, out.width = "100%"----------------------
data(toy_backbone)

# Plot the backbone with the tip labels in green, and the node IDs in label boxes
ape::plot.phylo(toy_backbone, cex = 0.7, no.margin = TRUE, tip.color = "green3")
ape::nodelabels(cex = 0.7, frame = "rect", bg = "lightblue", col = "darkblue")

## ----taxonomy_look------------------------------------------------------------
# ------ Step 1: Look at the taxonomic tables ------ #

# Taxonomy of the backbone: the first column holds the tip labels (`Taxa`)
data(toy_backbone_taxonomic_table)
head(toy_backbone_taxonomic_table)

# Description of the taxa to graft: one row per taxon
data(toy_missing_taxa_taxonomic_table)
toy_missing_taxa_taxonomic_table

## ----taxonomy_import----------------------------------------------------------
# ------ Step 2: Import the tables from external files (optional) ------ #

backbone_path <- system.file("extdata", "toy_backbone_taxonomic_table.csv", package = "phylograftR")
missing_taxa_path <- system.file("extdata", "toy_missing_taxa_taxonomic_table.csv",
                                 package = "phylograftR")

# Import the taxonomy of the backbone, checked against the tips of the phylogeny
backbone_taxonomic_table <- import_taxonomic_table(file = backbone_path, type = "backbone",
                                                   phy = toy_backbone, verbose = FALSE)

# Import the description of the taxa to graft
missing_taxa_taxonomic_table <- import_taxonomic_table(file = missing_taxa_path, type = "missing_taxa",
                                                       phy = toy_backbone, verbose = FALSE)

## ----taxonomy_build-----------------------------------------------------------
# ------ Step 3: Build the missing_taxa_list ------ #

# Taxonomic names are converted into clades of the backbone phylogeny
missing_taxa_lists <- build_missing_taxa_list_from_taxonomy(
  backbone_taxonomic_table = toy_backbone_taxonomic_table,
  missing_taxa_taxonomic_table = toy_missing_taxa_taxonomic_table,
  phy = toy_backbone,
  verbose = FALSE)

# Since both types of taxa were provided, the output has one list per grafting function
names(missing_taxa_lists)
str(missing_taxa_lists$taxa, max.level = 1)

# Clades are written as pairs of tips, and named after their taxonomic groups
missing_taxa_lists$taxa$Target_clade$New_sp1
missing_taxa_lists$taxa$Target_clade_name

# Taxa to graft as a clade have the same sister clade: they will be grafted together
missing_taxa_lists$clades$Sister_clade_name

# Information on the groups used (clade node, number of tips, monophyly...) can be inspected
# in the attribute "taxonomic_groups"
attr(missing_taxa_lists, "taxonomic_groups")[, c("Group", "Level", "clade_node",
                                                 "n_tips_clade", "monophyletic")]

## ----build_df-----------------------------------------------------------------
# ------ Step 1: Build the missing_taxa_list from a data.frame ------ #

# Write the table in R, with a mix of tips and node IDs (see the plot of the backbone)
my_table <- data.frame(Taxon_name = c("M1", "M2", "M3"),
                       Target_clade = c("34", "t5, t12", "t26, t18"),
                       Excluded_subclades = c("50", NA, "t10, t19; t26, t4"),
                       Excluded_subclades_names = c(NA, NA, "group_1; group_2"),
                       include_target_clade_stem = c(TRUE, FALSE, NA))
my_table

# The table is checked against the backbone phylogeny.
# Node IDs are converted into tips, and the conversions are printed
my_missing_taxa_list <- build_missing_taxa_list(data = my_table, phy = toy_backbone)
str(my_missing_taxa_list, max.level = 1)

# Node 34 is now defined by a pair of tips whose MRCA is node 34
my_missing_taxa_list$Target_clade$M1

## ----import_csv---------------------------------------------------------------
# ------ Step 2 (optional): Import the table from an external file ------ #

# Path to a table provided with the package
file_path <- system.file("extdata", "toy_missing_taxa.csv", package = "phylograftR")

# One row per taxon, one column per item of the list
read.csv(file_path)

# Import it as a missing_taxa_list, checked against the backbone phylogeny
missing_taxa_list <- import_missing_taxa_list(file = file_path, phy = toy_backbone)

# Items are named after the taxa, tips are split, and excluded subclades are named
str(missing_taxa_list, max.level = 1)

# Example of the clades of one taxon
missing_taxa_list$Target_clade$New_sp1
missing_taxa_list$Excluded_subclades$New_sp1
missing_taxa_list$exclude_subclades_stem$New_sp1

# The result is the same as reading the file and building the list from the data.frame
identical(missing_taxa_list, build_missing_taxa_list(data = read.csv(file_path), phy = toy_backbone))

# The same list is provided as a dataset
data(toy_missing_taxa_list)
identical(missing_taxa_list, toy_missing_taxa_list)

## ----import_xlsx, eval = FALSE------------------------------------------------
# # Import an Excel file
# missing_taxa_list <- import_missing_taxa_list(file = "my_missing_taxa.xlsx",
#                                               phy = toy_backbone, sheet = "missing_taxa")
# 
# # Import any other format with a reader function (e.g., '.ods')
# missing_taxa_list <- import_missing_taxa_list(file = "my_missing_taxa.ods", phy = toy_backbone,
#                                               reader = readODS::read_ods)

## ----convert_nodes------------------------------------------------------------
# A list with clades defined by node IDs, tips, or both
list_with_nodes <- list(
  Taxon_name = c("X1", "X2", "X3"),
  Target_clade = list(X1 = 34, X2 = c("t5", "t12"), X3 = "t26"),
  Excluded_subclades = list(X1 = list(group_1 = 50), X2 = NULL, X3 = NULL))

# Convert the node IDs into tips: the conversions are printed
list_with_tips <- convert_nodes_to_tips(list_with_nodes, phy = toy_backbone)
list_with_tips$Target_clade$X1

## ----allowed_states-----------------------------------------------------------
file_path <- system.file("extdata", "toy_missing_taxa_with_allowed_states.csv", package = "phylograftR")
read.csv(file_path)[, c("Taxon_name", "Target_clade_name", "allowed_states")]

# In this example, the allowed states define color (i.e., "red", "blue", "green") that evolve along the branches of a `simmap`
# and are used to contraint the grafting of each taxon depending on its allowed states.

missing_taxa_list_states <- import_missing_taxa_list(file = file_path, phy = toy_backbone, verbose = FALSE)
missing_taxa_list_states$allowed_states

## ----both_types---------------------------------------------------------------
# ------ Step 1: A table with both types ------ #

both_table <- data.frame(Taxon_name = c("M1", "M2", "M3", "M4"),
                         Type = c("taxa", "taxa", "clades", "clades"),
                         Target_clade = c("t26, t18", "t5, t12", NA, NA),
                         Sister_clade = c(NA, NA, "t10, t19", "t10, t19"),
                         Excluded_subclades = c("t26, t4", NA, NA, NA))

# ------ Step 2: Build the lists ------ #

lists <- build_missing_taxa_list(data = both_table, phy = toy_backbone)
names(lists)
lists$taxa$Taxon_name      # Taxa to graft with graft_missing_taxa()
lists$clades$Taxon_name    # Taxa to graft with graft_missing_clades()
lists$clades$Sister_clade  # M3 and M4 have the same sister clade: they form one missing clade

## ----manual_edit--------------------------------------------------------------
# ------ Step 1: Edit the list ------ #

# Start from the list of taxa built from the taxonomy (Section 3)
taxa_list <- missing_taxa_lists$taxa

# Allow New_sp3 to be grafted on the stem branch of its target clade
taxa_list$include_target_clade_stem[["New_sp3"]] <- TRUE

# Add an excluded subclade (and exclude its stem) for New_sp4
taxa_list$Excluded_subclades[["New_sp4"]] <- list(Gen_1 = c("t13", "t27"))
taxa_list$exclude_subclades_stem[["New_sp4"]] <- c(Gen_1 = TRUE)

# Remove a taxon: subset all the items with the same selection
keep <- taxa_list$Taxon_name != "New_sp7"
taxa_list <- lapply(taxa_list, function (item) item[keep])

# ------ Step 2: Check the edited list ------ #

# The structure of the list and the availability of branches for grafting
# can be checked with check_grafting_plan() (see `vignette("check_grafting")`)
plan <- check_grafting_plan(phy = toy_backbone, missing_taxa_list = taxa_list, verbose = FALSE)
plan[, c("Taxon_name", "n_branches", "feasible")]

## ----total_richness-----------------------------------------------------------
# ------ Step 1: A table with the clades and their total richness ------ #

# Clades are written with their taxonomic names, as in the list built from the taxonomy.
# Use NA if the total richness is unknown
richness_table <- data.frame(Clade = c("Family_A", "Family_B", "Subfam_B1", "Tribe_B2"),
                             total_richness = c(25, NA, 14, 8))

# ------ Step 2: Build the total_richness_per_clades ------ #

# Check against the backbone and the lists of both types
richness_per_clades <- build_total_richness_per_clades(data = richness_table, phy = toy_backbone,
                                                       missing_taxa_list = missing_taxa_lists,
                                                       verbose = FALSE)

# With both types, the output has one element per type. Clades are defined by the same tips as in the lists
names(richness_per_clades)
str(richness_per_clades, max.level = 2)

# ------ Step 3: Check the richness before grafting ------ #

# Expected richness and sampling fraction of each clade
attr(richness_per_clades, "richness_table")[, c("Clade_name", "total_richness", "nb_backbone_species",
                                                "nb_missing_species", "nb_unsampled_species",
                                                "sampling_fraction")]

## ----total_richness_import----------------------------------------------------
file_path <- system.file("extdata", "toy_total_richness.csv", package = "phylograftR")
read.csv(file_path)

# The richness table of the toy dataset is defined with tips, and goes with the toy missing_taxa_list
richness <- import_total_richness_per_clades(file = file_path, phy = toy_backbone,
                                             missing_taxa_list = toy_missing_taxa_list,
                                             verbose = FALSE)
str(richness, max.level = 1)

## ----total_richness_convert---------------------------------------------------
# A missing_taxa_list written by hand with node IDs and tips (see Section 5)
list_with_nodes <- list(
  Taxon_name = c("X1", "X2"),
  Target_clade = list(X1 = 34, X2 = c("t5", "t12")))

# A total_richness_per_clades written by hand: node 34 is the clade of X1
richness_with_nodes <- list(Target_clade = list(34, c("t5", "t12")),
                            total_richness = c(14, 8))

# Convert both objects at the same time
converted <- convert_nodes_to_tips(missing_taxa_list = list_with_nodes, phy = toy_backbone,
                                   total_richness_per_clades = richness_with_nodes)
names(converted)

# Node 34 is defined by the same tips in the list and in the total richness
converted$missing_taxa_list$Target_clade$X1
converted$total_richness_per_clades$Target_clade[[1]]

## ----common_errors------------------------------------------------------------
# Helper to display the message of an error
show_error <- function (expr) cat("Error:", conditionMessage(tryCatch(expr, error = function (e) e)), "\n")

# A tip that is not in the backbone (e.g., a typo, or a wrong separator)
show_error(build_missing_taxa_list(data.frame(Taxon_name = c("M1", "M2"), Target_clade = c("t1, t999", "t5")),
                                   phy = toy_backbone, verbose = FALSE))

# Tips and a node ID mixed in the same clade
show_error(build_missing_taxa_list(data.frame(Taxon_name = c("M1", "M2"), Target_clade = c("t1, 34", "t5")),
                                   phy = toy_backbone, verbose = FALSE))

# A node ID that does not exist
show_error(build_missing_taxa_list(data.frame(Taxon_name = c("M1", "M2"), Target_clade = c("t8, t7", "900")),
                                   phy = toy_backbone, verbose = FALSE))

# Duplicated names of taxa
show_error(build_missing_taxa_list(data.frame(Taxon_name = c("M1", "M1"), Target_clade = c("t8, t7", "t5")),
                                   phy = toy_backbone, verbose = FALSE))

# A taxon to graft that is already a tip of the backbone
show_error(build_missing_taxa_list(data.frame(Taxon_name = c("t5", "M2"), Target_clade = c("t8, t7", "t5")),
                                   phy = toy_backbone, verbose = FALSE))

# A taxon without clade
show_error(build_missing_taxa_list(data.frame(Taxon_name = c("M1", "M2"), Target_clade = c("t8, t7", NA)),
                                   phy = toy_backbone, verbose = FALSE))

