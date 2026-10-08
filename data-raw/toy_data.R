##### Build the toy datasets of phylograftR #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Run from the root of the package, after loading it (devtools::load_all() or devtools::install()).
# Produces:
#  * data/toy_backbone.rda, data/toy_missing_taxa_list.rda, data/toy_total_richness_per_clades.rda
#  * inst/extdata/toy_missing_taxa.csv, inst/extdata/toy_total_richness.csv


##### 1/ Create backbone phylogeny #####

# Backbone: 30 tips, ultrametric, root age of 50 -----------------------------------------------

set.seed(1234)
toy_backbone <- ape::rcoal(30)
toy_backbone$edge.length <- toy_backbone$edge.length / max(ape::node.depth.edgelength(toy_backbone)) * 50

# Reference clades (node: tips) used below
# 32: t26 t2 t21 t16 t11 t4 t3 t25 t28 t10 t6 t20 t9 t30 t24 t18 t19  (clade A)
# 33: t8 t1 t15 t29 t14 t13 t17 t27 t7 t5 t23 t12 t22                 (clade B)
# 34: t8 t1 t15 t29 t14 t13 t17 t27 t7                                (clade B1, nested in B)
# 37: t10 t6 t20 t9 t30 t24 t18 t19                                   (group_1, nested in A)
# 39: t26 t2 t21 t16 t11 t4                                           (group_2, nested in A)
# 50: t13 t17 t27                                                     (group_3, nested in B1)
# 53: t5 t23 t12                                                      (clade C, nested in B)

##### 2/ Create table of the missing taxa as .csv #####

toy_missing_taxa_table <- data.frame(
  Taxon_name = c("Genus_a_sp1", "Genus_a_sp2", "Genus_a_sp3", "Genus_b_sp1", "Genus_b_sp2", "Genus_c_sp1", "Genus_c_sp2"),
  Target_clade = c("t26,t18", "t26,t18", "t26,t18", "t8,t22", "t8,t7", "t5", "t5,t12"),
  Target_clade_name = c("clade_A", "clade_A", "clade_A", "clade_B", "clade_B1", "tip_t5", "clade_C"),
  include_target_clade_stem = c(TRUE, TRUE, FALSE, NA, TRUE, NA, FALSE),
  Excluded_subclades = c("t10,t19;t26,t4", "t26,t4;t30", NA, NA, "t13,t27", NA, NA),
  Excluded_subclades_names = c("group_1;group_2", "group_2;tip_t30", NA, NA, "group_3", NA, NA),
  exclude_subclades_stem = c("TRUE;FALSE", "FALSE;TRUE", NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE)

# Export in inst/extdata/
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(toy_missing_taxa_table, file = "inst/extdata/toy_missing_taxa.csv", row.names = FALSE, na = "")

##### 3/ Convert to missing_taxa_list object #####

toy_missing_taxa_list <- phylograftR::import_missing_taxa_list(
  file = "inst/extdata/toy_missing_taxa.csv", phy = toy_backbone)

##### 4/ Create table of Total richness of the target clades as .csv #####

toy_total_richness_table <- data.frame(
  Target_clade = c("t26,t18", "t8,t22", "t8,t7", "t5", "t5,t12"),
  total_richness = c(25, NA, 14, NA, 8))

# Export in inst/extdata
utils::write.csv(toy_total_richness_table, file = "inst/extdata/toy_total_richness.csv", row.names = FALSE, na = "")

toy_total_richness_per_clades <- list(
  Target_clade = strsplit(toy_total_richness_table$Target_clade, ",", fixed = TRUE),
  total_richness = toy_total_richness_table$total_richness)

##### 5/ Save objects #####

dir.create("data", showWarnings = FALSE)
save(toy_backbone, file = "data/toy_backbone.rda", compress = "gzip", version = 2)
save(toy_missing_taxa_list, file = "data/toy_missing_taxa_list.rda", compress = "gzip", version = 2)
save(toy_total_richness_per_clades, file = "data/toy_total_richness_per_clades.rda", compress = "gzip", version = 2)
