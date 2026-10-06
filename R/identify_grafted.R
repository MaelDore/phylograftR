##### Function to identify grafted tips and branches #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Identify the grafted tips and branches of a phylogeny ####

#' @title Identify the grafted tips and branches of a phylogeny
#'
#' @description Identify which tips and which branches of a phylogeny are grafted, as opposed to backbone.
#'   The information is added to the phylogeny as two elements, `grafted_tips` and `grafted_edges`, which can be used to color the phylogeny on plots (see examples).
#'   The function is also useful to update these elements after the phylogeny has been manipulated
#'   (e.g., with `ape::drop.tip()`, which does not update the elements of the phylogeny it does not know).
#'
#' @param phy An object of class `"phylo"`. The phylogeny, including the grafted taxa.
#' @param grafted_taxa Character vector. Labels of the grafted taxa (e.g., `missing_taxa_list$Taxon_name`).
#'   Labels that are not tips of `phy` are ignored, with a warning.
#'
#' @details A tip is grafted if its label is in `grafted_taxa`.
#'   A branch is grafted if all the tips descending from it are grafted taxa (this includes the terminal branches of grafted taxa,
#'   and the internal branches of clades made only of grafted taxa).
#'   A branch with at least one backbone tip below it is a backbone branch, even if a grafted taxon was added on it
#'   (the backbone branch is then split in two pieces, which are both backbone branches).
#'
#'   Node IDs are renumbered by functions that remove tips, such as `ape::drop.tip()`:
#'   use this function again to update the `grafted_tips` and `grafted_edges` elements in that case.
#'
#' @importFrom ape Ntip reorder.phylo
#'
#' @export
#'
#' @seealso Grafting functions, which already add internally the two elements to the grafted phylogeny with this function:
#'   * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade with constraints (target clades and excluded subclades),
#'   * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as clades sister to a clade of the backbone.
#'
#' @return The phylogeny `phy` (class `"phylo"`), with two elements added, or replaced if `phy` already had them:
#'  * `grafted_tips` Named logical vector with one value per tip of `phy`, in the order of `phy$tip.label`.
#'    `TRUE` for grafted tips. Names are the tip labels.
#'  * `grafted_edges` Named logical vector with one value per branch (edge) of `phy`, in the order of the rows of `phy$edge`.
#'    `TRUE` for grafted branches. Names are the IDs of the tipward nodes of the branches (i.e., `phy$edge[, 2]`).
#'
#' @author Maël Doré
#'
#' @examples
#' # Create a toy backbone phylogeny, and graft taxa
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")))
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   verbose = FALSE, seed = 1234)
#'
#' # Identify the grafted tips and branches:
#' # the function returns the phylogeny with two additional elements
#' grafted_phy <- identify_grafted(phy = grafted_phy, grafted_taxa = missing_taxa_list$Taxon_name)
#' table(grafted_phy$grafted_tips)
#' table(grafted_phy$grafted_edges)
#'
#' # Color the grafted tips and branches in red
#' plot(grafted_phy, cex = 0.7,
#'      edge.color = c("black", "red")[grafted_phy$grafted_edges + 1],
#'      tip.color = c("black", "red")[grafted_phy$grafted_tips + 1])
#'
#' # Prune the phylogeny: drop a grafted tip and a backbone tip
#' pruned_phy <- ape::drop.tip(phy = grafted_phy, tip = c("X2", "t1"))
#' # The elements 'grafted_tips' and 'grafted_edges' are NOT updated by ape::drop.tip():
#' # they still have their previous lengths, and are not valid anymore
#' length(pruned_phy$grafted_tips) == ape::Ntip(pruned_phy)
#' length(pruned_phy$grafted_edges) == nrow(pruned_phy$edge)
#'
#' # Update them with identify_grafted(), using the grafted taxa still in the pruned phylogeny
#' # (grafted taxa that are not tips are ignored with a warning)
#' pruned_phy <- identify_grafted(
#'    phy = pruned_phy,
#'    grafted_taxa = intersect(missing_taxa_list$Taxon_name, pruned_phy$tip.label))
#'
#' # Plot the pruned phylogeny with the updated colors
#' plot(pruned_phy, cex = 0.7,
#'      edge.color = c("black", "red")[pruned_phy$grafted_edges + 1],
#'      tip.color = c("black", "red")[pruned_phy$grafted_tips + 1])
#'

identify_grafted <- function (phy, grafted_taxa)
{
  ## Check input
  if (!inherits(phy, "phylo")) stop("'phy' must be an object of class 'phylo'")
  if (!is.character(grafted_taxa)) stop("'grafted_taxa' must be a character vector of tip labels")

  # Get total number of tips
  n_tips <- ape::Ntip(phy)

  absent <- setdiff(grafted_taxa, phy$tip.label)
  if (length(absent) > 0)
  {
    warning(length(absent), " taxa of 'grafted_taxa' are not tips of 'phy' and are ignored: ",
            paste(utils::head(absent, 10), collapse = ", "), if (length(absent) > 10) ", ..." else "", call. = FALSE)
  }

  ## Identify grafted tips by matching with tip.label
  is_grafted_tip <- phy$tip.label %in% grafted_taxa

  ## Count the number of backbone tips below each node (tips included): 1 for backbone tips, 0 for grafted tips
  # Get total number of nodes (tips + internal)
  n_backbone <- numeric(n_tips + phy$Nnode)
  # Inform tip status
  n_backbone[seq_len(n_tips)] <- as.numeric(!is_grafted_tip)
  # Order phylo in 'postorder' fashion
  # In 'postorder', the branches of a subtree come before the branch leading to its root: children are transverse before being added to their parent
  postorder_edges <- ape::reorder.phylo(phy, order = "postorder")$edge
  # Loop per edge in postorder fashion
  for (i in seq_len(nrow(postorder_edges)))
  {
    # Add counts of grafted childrens iteratively
    n_backbone[postorder_edges[i, 1]] <- n_backbone[postorder_edges[i, 1]] + n_backbone[postorder_edges[i, 2]]
  }

  ## Identify grafted branches with no backbone tip below their tipward node
  is_grafted_edge <- n_backbone[phy$edge[, 2]] == 0

  ## Name outputs
  # Use tip.label for tips
  names(is_grafted_tip) <- phy$tip.label
  # Use tipward node ID for branches/edges
  names(is_grafted_edge) <- as.character(phy$edge[, 2])

  ## Add the elements to the phylogeny (replaced if already present) and return it
  phy$grafted_tips <- is_grafted_tip
  phy$grafted_edges <- is_grafted_edge
  return(phy)
}
