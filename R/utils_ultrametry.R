##### Internal helpers for ultrametricity, shared by grafting_missing_taxa() and grafting_missing_clades() #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

## Utility functions used to enforce phylogeny ultrametry and displaying proper warnings depending on context


##### 1/ Check ultrametricity of the backbone tree #####

# The backbone must be ultrametric BEFORE grafting.
# It is checked, never modified by the grafting functions: making it ultrametric (e.g. when the deviation is a rounding error) is left to the user.
# Stop if the backbone has no branch lengths or is not ultrametric

check_backbone_ultrametric <- function (phy, # Backbone phylogeny
                                        tol = .Machine$double.eps^0.5) # Tolerance used by ape::is.ultrametric()
{
  if (is.null(phy$edge.length))
  {
    stop("The backbone phylogeny has no branch lengths. A time-calibrated (ultrametric) phylogeny is required.")
  }
  if (!ape::is.ultrametric(phy, tol = tol))
  {
    stop(paste0("The backbone phylogeny is not ultrametric (tolerance = ", format(tol), ").\n",
                "Grafting requires an ultrametric backbone. If you think this is a rounding error and the tree should be considered ultrametric, ",
                "you can use phytools::force.ultrametric(tree = phy, method = \"nnls\") (or method = \"extend\") ",
                "on the backbone BEFORE grafting and provide this ultrametric backbone as input."))
  }
  invisible(TRUE)
}

##### 2/ Enforce ultrametricity of new tips #####

# Reset the length of the terminal branch of each new tip such as the tip is exactly at root_age from the root.
# Only the new terminal branches are modified, so the branch lengths of the backbone remain untouched
# (contrary to phytools::force.ultrametric() applied on the whole tree).

enforce_ultrametry_new_tips <- function (phy, # Phylogeny including the new tips
                                         tip_labels, # Labels of the new tips
                                         root_age) # Age of the root = expected distance from the root to all tips
{
  for (tip_label in tip_labels)
  {
    tip_node <- match(x = tip_label, table = phy$tip.label)
    tip_edge_ID <- match(x = tip_node, table = phy$edge[, 2])
    parent_node <- phy$edge[tip_edge_ID, 1]
    phy$edge.length[tip_edge_ID] <- max(0, root_age - phytools::nodeheight(tree = phy, node = parent_node))
  }
  return(phy)
}


##### 3/ Check before enforcing overall ultrametry #####

# (nnls method is very costly: dense ntips x nedges matrix)
# Favor using extend by default

ensure_ultrametry <- function (phy, # Phylogeny
                               method = "extend",
                               tol = .Machine$double.eps^0.5) # Tolerance used by ape::is.ultrametric()
{
  is_ultrametric <- ape::is.ultrametric(phy, tol = tol)
  if (!is_ultrametric)
  {
    # Force ultrametry to avoid issue with bad rounding
    phy <- phytools::force.ultrametric(tree = phy, method = method, message = FALSE)
  }

  return(phy)
}


