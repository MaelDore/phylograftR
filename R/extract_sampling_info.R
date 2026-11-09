##### Function to extract the sampling fractions of clades from grafted phylogenies #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Extract the sampling information ####

#' @title Extract the richness and sampling fraction of clades from grafted phylogenies
#'
#' @description Summarize, in each grafted phylogeny, the species of each clade listed in `total_richness_per_clades`:
#'   species of the backbone, grafted species, and unsampled species (species of the clade that are neither in the backbone nor grafted).
#'   It then derives the **sampling fraction** of each clade that can be used in downstream analyses
#'   to account for the clade-specific sampling fractions (e.g., diversification models).
#'
#'   The function works **a posteriori**, on the output of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades]
#'   (a phylogeny, a set of simulations, or the list of maps obtained with a `"multiSimmap"`).
#'   The counts are taken in each grafted phylogeny. With `graft_missing_taxa()`, **they can differ among simulations and maps**, as a taxon can end up in different clades.
#'   With `graft_missing_clades()`, **they are the same in all the phylogenies**, as each missing clade is always grafted on the stem of its sister clade.
#'
#' @param phy_grafted A grafted phylogeny (class `"phylo"`), a set of grafted phylogenies (class `"multiPhylo"`, e.g., several simulations),
#'   or a list of them with one element per map (output of the grafting functions with a `"multiSimmap"`).
#'   The phylogenies must have the elements `$grafted_tips` and `$grafted_edges` (see [`phylograftR::identify_grafted()`][identify_grafted]).
#' @param total_richness_per_clades A list with the total number of species of the clades, as used for the grafting. Items:
#'   * `$Target_clade` (clades defined for the taxa grafted with [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]) or `$Sister_clade`
#'     (sister clades of the missing clades grafted with [`phylograftR::graft_missing_clades()`][graft_missing_clades]). **Exactly one of them must be provided.**
#'     List of character vectors of tips of the backbone: a pair of tips defines the clade as their most recent common ancestor (MRCA), and a single tip defines a terminal branch.
#'     **Node IDs are not accepted**: convert them with [`phylograftR::convert_nodes_to_tips()`][convert_nodes_to_tips].
#'   * `$Target_clade_name` or `$Sister_clade_name`: (Optional) Character vector of the names of the clades, same length as the clades (`NA` if a clade has no name).
#'     The names must be unique. They are reported in the `Clade_name` column of the output.
#'   * `$total_richness`: Vector of positive integers, same length as the clades. Total number of species of each clade,
#'     i.e., backbone species, grafted species, and unsampled species. Use `NA` if unknown.
#'
#' @details
#' ### Species of a clade
#'
#' **The species of a clade are the tips that descend from its MRCA in the grafted phylogeny**: the backbone species of the clade, and the grafted taxa that were grafted within it.
#' It depends on the phylogeny, as a taxon grafted in a target clade can fall in a clade nested within it, and as the taxa grafted on the stem branch of a clade
#' are sister to the clade, so they are not part of it. The exception is the clade of a single tip of `$Target_clade`: its terminal branch is its only branch,
#' so the taxa grafted on it are part of the clade. For `$Sister_clade`, the missing clades are grafted on the stem branch of the sister clade and are never part of it.
#'
#' ### Unsampled species and sampling fraction
#'
#' For a clade with a total richness provided, `nb_unsampled_species = total_richness - nb_backbone_species - nb_grafted_species`.
#' For a clade without a total richness (`NA`), the unsampled species are those identified in the largest clades nested within it that have a total richness (none if there is none),
#' and its total richness is derived as `nb_backbone_species + nb_grafted_species + nb_unsampled_species` (it can differ among phylogenies).
#' The **sampling fraction** is `(nb_backbone_species + nb_grafted_species) / total_richness`.
#'
#' In problematic cases, the number of species of a clade can be higher than its total richness: for example, a taxon known to belong to a genus can be grafted in a subclade
#' of this genus whose provided total richness did not include it. The number of unsampled species is then negative and the sampling fraction is higher than 1.
#' These values are returned as they are, with a warning listing the clades and the number of phylogenies concerned, so the clades and the total richness provided can be reviewed.
#' When a clade without a total richness is derived from nested clades, the negative values of the nested clades are counted as `0`.
#'
#' @importFrom ape Ntip drop.tip getMRCA extract.clade
#'
#' @export
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa based on target clades and excluded subclades.
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as sister clades on backbone clades.
#'  * [`phylograftR::extract_grafting_info()`][extract_grafting_info] to extract a posteriori information on the location of grafted lineages.
#'  * [`phylograftR::extract_BD_model_summary()`][extract_BD_model_summary] to extract a posteriori information on the birth-death models fitted during grafting.
#'
#' @return A data.frame, with one row per clade and phylogeny. The columns are:
#'   * `map` Character string. Name of the stochastic map (`NA` if no map was used).
#'   * `simulation` Integer. Index of the simulation.
#'   * `seed` Integer. Seed of the simulation (`NA` if unknown).
#'   * `Clade_name` Character string. Name of the clade (`NA` if it has no name). Only if the names of the clades are provided in `total_richness_per_clades`.
#'   * `Clade_tips` Character string. Tips used to identify the clade.
#'   * `total_richness` Numeric. Total number of species of the clade: provided, or derived if `NA` was provided.
#'   * `total_richness_provided` Logical. Whether the total richness of the clade was provided.
#'   * `nb_backbone_species` Integer. Number of species of the clade that are in the backbone.
#'   * `nb_grafted_species` Integer. Number of species of the clade that were grafted.
#'   * `nb_unsampled_species` Numeric. Number of species of the clade that are neither in the backbone nor grafted. Negative (conflictual) if the clade has more species than expected.
#'   * `sampling_fraction` Numeric. Fraction of the clade richness present in the grafted phylogeny, `(nb_backbone_species + nb_grafted_species) / total_richness`.
#'
#'   Attributes: `"type"` (`"taxa"` if the clades are `$Target_clade`, `"clades"` if they are `$Sister_clade`).
#'
#' @author Maël Doré
#'
#' @examples
#' # --- Prepare data --- #
#'
#' # Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' # Two nested clades: B is nested within A
#' clade_A <- c("t26", "t18") # 17 backbone species
#' clade_B <- c("t26", "t3")  # 9 backbone species, nested within A
#'
#' # List the missing taxa to graft: 3 taxa in clade A and 1 taxon in clade B
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   Target_clade = list(X1 = clade_A, X2 = clade_A, X3 = clade_A, X4 = clade_B),
#'   # Optional names of the clades (one per taxon): they must match those of total_richness_per_clades
#'   Target_clade_name = c("clade_A", "clade_A", "clade_A", "clade_B"))
#'
#' # Total richness of the target clades
#' # NA = unknown richness (i.e., no unsampled species beyond nested clades)
#' total_richness_per_clades <- list(
#'   Target_clade = list(clade_A, clade_B),
#'   Target_clade_name = c("clade_A", "clade_B"),
#'   total_richness = c(25, 12))
#'
#' # --- Graft the taxa (several simulations) --- #
#'
#' grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                    total_richness_per_clades = total_richness_per_clades,
#'                                    nb_simulations = 5, verbose = FALSE, seed = 1234)
#'
#' # --- Extract the sampling information --- #
#'
#' sampling_info <- extract_sampling_info(phy_grafted = grafted_phys,
#'                                        total_richness_per_clades = total_richness_per_clades)
#' sampling_info[, c("simulation", "Clade_name", "nb_backbone_species", "nb_grafted_species",
#'                   "nb_unsampled_species", "sampling_fraction")]
#'
#' # The sampling fraction of clade B varies among simulations:
#' # the taxa assigned to clade A can be grafted within clade B (nested in A, so they count for B)
#' # or in A outside of B. As such, the richness of B differs between grafted phylogenies,
#' # while all the taxa are always in A (constant sampling fraction).
#' tapply(X = sampling_info$sampling_fraction,
#'        INDEX = list(sampling_info$Clade_name, sampling_info$simulation),
#'        FUN = mean)
#'
#'

extract_sampling_info <- function (phy_grafted, # Output of the grafting functions: phylo, multiPhylo, or list of maps
                                   total_richness_per_clades) # Same list as used for grafting. Items: Target_clade (or Sister_clade), total_richness
{
  ### ---- Check input validity ---- ###

  if (!is.list(total_richness_per_clades)) { stop("'total_richness_per_clades' must be a list with the items 'Target_clade' (or 'Sister_clade') and 'total_richness'") }
  has_target <- !is.null(total_richness_per_clades[["Target_clade"]])
  has_sister <- !is.null(total_richness_per_clades[["Sister_clade"]])
  if (has_target && has_sister) { stop("'total_richness_per_clades' has both 'Target_clade' and 'Sister_clade': use only one, depending on the grafting function (graft_missing_taxa() or graft_missing_clades())") }
  if (!has_target && !has_sister) { stop("'total_richness_per_clades' must have a 'Target_clade' item (taxa grafted by graft_missing_taxa()) or a 'Sister_clade' item (missing clades grafted by graft_missing_clades())") }
  is_sister <- has_sister
  clade_element <- if (is_sister) "Sister_clade" else "Target_clade"
  clades <- total_richness_per_clades[[clade_element]]
  values <- total_richness_per_clades[["total_richness"]]
  if (is.null(values)) { stop("'total_richness_per_clades' must have the item 'total_richness'") }

  node_clades <- which(vapply(clades, is_node_id_clade, logical(1)))
  if (is.list(clades) && length(node_clades) > 0)
  {
    stop(length(node_clades), " clade(s) of 'total_richness_per_clades$", clade_element, "' are defined by node IDs, but tips are required. ",
         "Convert them with convert_nodes_to_tips() (see its argument 'total_richness_per_clades'), and use the converted list.", call. = FALSE)
  }
  if (!is.list(clades) || length(clades) == 0 || !all(vapply(clades, function (x) is.character(x) && length(x) >= 1, logical(1))))
  {
    stop("'total_richness_per_clades$", clade_element, "' must be a non-empty list of character vectors of tips")
  }
  if (length(values) != length(clades)) { stop("'total_richness_per_clades$total_richness' must have the same length as 'total_richness_per_clades$", clade_element, "'") }
  if (is.logical(values) && all(is.na(values))) { values <- as.numeric(values) } # A vector with only NA is logical: accept it
  if (!is.numeric(values) || any(!is.na(values) & (values < 0 | values != round(values))))
  {
    stop("'total_richness_per_clades$total_richness' must be a numeric vector of non-negative whole numbers (NA allowed)")
  }

  ## List the phylogenies, and check that they identify their grafted tips
  listed <- list_grafted_phylogenies(phy_grafted, arg = "phy_grafted")
  for (phy in listed$phylos)
  {
    if (is.null(phy$grafted_tips) || length(phy$grafted_tips) != ape::Ntip(phy))
    {
      stop("The phylogenies must have a '$grafted_tips' element matching the tips (as returned by the grafting functions). Use identify_grafted() to add or update it.")
    }
  }

  ## Backbone: the first phylogeny without its grafted taxa (all the phylogenies share the backbone)
  first <- listed$phylos[[1]]
  backbone <- ape::drop.tip(first, first$tip.label[first$grafted_tips])
  all_tips <- unique(unlist(clades))
  if (!all(all_tips %in% backbone$tip.label))
  {
    stop("Some tip labels provided in 'total_richness_per_clades$", clade_element, "' are not tips of the backbone of the grafted phylogenies\n",
         "Missing tips: ", paste(setdiff(all_tips, backbone$tip.label), collapse = ", "))
  }

  ## Clades are matched by their MRCA node in the backbone: they cannot be listed twice
  clade_nodes <- vapply(clades, function (tips) get_clade_node(backbone, tips), numeric(1))
  clade_labels <- vapply(clades, function (tips) paste(tips, collapse = ", "), character(1))
  if (anyDuplicated(clade_nodes))
  {
    dup <- unique(clade_nodes[duplicated(clade_nodes)])
    stop("Some clades are listed more than once in 'total_richness_per_clades$", clade_element, "' (clades are matched by their MRCA node):\n",
         paste0("  - ", vapply(dup, function (nd) paste(clade_labels[clade_nodes == nd], collapse = "  =  "), character(1)), collapse = "\n"))
  }

  ## Names of the clades (optional): valid and unique
  name_item <- paste0(clade_element, "_name")
  has_names <- !is.null(total_richness_per_clades[[name_item]])
  check_clade_names(clade_names = total_richness_per_clades[[name_item]], clade_nodes = clade_nodes, clade_labels = clade_labels,
                    item = paste0("total_richness_per_clades$", name_item), n_clades = length(clades))
  clade_names <- if (has_names) { x <- unname(as.character(total_richness_per_clades[[name_item]])); x[!is.na(x) & !nzchar(x)] <- NA_character_; x } else NULL

  ## Nesting of the clades: largest clades nested within each clade (backbone tips of the clades)
  backbone_tips <- lapply(clade_nodes, function (nd) if (nd <= ape::Ntip(backbone)) backbone$tip.label[nd] else ape::extract.clade(backbone, nd)$tip.label)
  nested_in <- lapply(seq_along(clades), function (k)
  {
    inside <- which(vapply(seq_along(clades), function (j) j != k && all(backbone_tips[[j]] %in% backbone_tips[[k]]), logical(1)))
    # Keep the largest clades: those not nested within another clade of 'inside'
    inside[!vapply(inside, function (j) any(vapply(setdiff(inside, j), function (i) all(backbone_tips[[j]] %in% backbone_tips[[i]]), logical(1))), logical(1))]
  })

  ### ---- Extract the counts of each phylogeny ---- ###

  results <- lapply(listed$phylos, function (phy)
  {
    grafted <- stats::setNames(as.logical(phy$grafted_tips), phy$tip.label)
    members <- lapply(clades, function (tips) clade_members(phy = phy, tips = tips, grafted = grafted, climb = !is_sister))
    nb_grafted <- vapply(members, function (m) sum(grafted[m]), numeric(1))
    nb_backbone <- vapply(members, function (m) sum(!grafted[m]), numeric(1))
    nb_tips <- nb_backbone + nb_grafted

    provided <- !is.na(values)
    nb_unsampled <- ifelse(provided, values - nb_tips, NA_real_)
    # Clades without a total richness: unsampled species of the largest nested clades with a total richness (negative values count as 0)
    # Computed from the smallest to the largest clade (in the backbone), as an unknown clade can be nested in another unknown clade
    for (k in order(lengths(backbone_tips)))
    {
      if (!provided[k])
      {
        nb_unsampled[k] <- sum(vapply(nested_in[[k]], function (j) max(nb_unsampled[j], 0), numeric(1)))
      }
    }
    total <- ifelse(provided, values, nb_tips + nb_unsampled)
    data.frame(Clade_tips = clade_labels, total_richness = total, total_richness_provided = provided,
               nb_backbone_species = as.integer(nb_backbone), nb_grafted_species = as.integer(nb_grafted), nb_unsampled_species = nb_unsampled,
               sampling_fraction = nb_tips / total, stringsAsFactors = FALSE)
  })

  ### ---- Long table ---- ###

  n_rows <- vapply(results, nrow, integer(1))
  output <- do.call(rbind, lapply(seq_along(results), function (k)
  {
    cbind(map = rep(listed$map[k], n_rows[k]), simulation = rep(listed$simulation[k], n_rows[k]), seed = rep(listed$seed[k], n_rows[k]), results[[k]], stringsAsFactors = FALSE)
  }))
  if (has_names) { output <- cbind(output[, c("map", "simulation", "seed")], Clade_name = rep(clade_names, length(results)), output[, -(1:3)], stringsAsFactors = FALSE) }
  rownames(output) <- NULL
  attr(output, "type") <- if (is_sister) "clades" else "taxa"

  ## Warning for the clades with more species than expected
  negative <- output$nb_unsampled_species < 0
  if (any(negative))
  {
    concerned <- unique(if (has_names) ifelse(is.na(output$Clade_name), output$Clade_tips, paste0(output$Clade_name, " (", output$Clade_tips, ")"))[negative] else output$Clade_tips[negative])
    n_show <- min(length(concerned), 5)
    warning(sum(negative), " row(s) have a negative number of unsampled species (more species in the clade than its total richness): ",
            length(concerned), " clade(s) in ", length(unique(paste(output$map[negative], output$simulation[negative]))), " phylogeny(ies).\n",
            paste0("  - ", concerned[seq_len(n_show)], collapse = "\n"), if (length(concerned) > n_show) paste0("\n  ... and ", length(concerned) - n_show, " more") else "",
            "\nTaxa grafted in a clade can come from a larger clade: review the total richness provided, or the clades listed.", call. = FALSE)
  }

  return(output)
}


#' @title Tips of a clade in a grafted phylogeny
#'
#' @description Tips of the grafted phylogeny that descend from the MRCA of the tips of the clade (backbone species and grafted taxa within the clade).
#'   For a clade defined by a single tip with `climb = TRUE`, the clade includes the taxa grafted on the terminal branch of the tip: the node is moved rootward
#'   as long as the other descendant lineage of its parent is made only of grafted taxa.
#'
#' @param phy Grafted phylogeny (class `"phylo"`).
#' @param tips Character vector. Tips of the backbone defining the clade.
#' @param grafted Named logical vector. Is the tip grafted? Named after the tip labels of `phy`.
#' @param climb Logical. Whether the taxa grafted on the terminal branch of a single tip are part of its clade.
#'
#' @importFrom ape getMRCA extract.clade
#'
#' @noRd
#'
#' @return A character vector of tip labels.
#'

clade_members <- function (phy, tips, grafted, climb)
{
  n_tips <- ape::Ntip(phy)
  if (length(tips) > 1)
  {
    node <- ape::getMRCA(phy, tips)
    return(ape::extract.clade(phy, node)$tip.label)
  }
  node <- match(tips, phy$tip.label)
  if (climb)
  {
    repeat
    {
      parent_row <- match(node, phy$edge[, 2])
      if (is.na(parent_row)) { break }
      parent <- phy$edge[parent_row, 1]
      sibling <- setdiff(phy$edge[phy$edge[, 1] == parent, 2], node)
      sibling_tips <- if (sibling <= n_tips) phy$tip.label[sibling] else ape::extract.clade(phy, sibling)$tip.label
      if (!all(grafted[sibling_tips])) { break }
      node <- parent
    }
  }
  if (node <= n_tips) phy$tip.label[node] else ape::extract.clade(phy, node)$tip.label
}
