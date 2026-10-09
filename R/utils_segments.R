
##### Internal helpers used to handle branch segments: the unit used to draw where a taxon is grafted #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

# Rationale
#   Previously, the branches available for a taxon were handled as whole branches (a set of tipward node IDs).
#   To be able to restrict grafting to a part of a branch (e.g., the section of the branch mapped in an allowed state),
#   the available branches are now handled as SEGMENTS: portions of branches described by their tipward node and by two distances
#   measured from this tipward node, up the branch (0 = tipward node, branch length = rootward node).
#   For the moment each available branch is a single segment covering the whole branch, so the draws are identical to the previous ones.
#   Any restriction of the available time along branches only has to shorten or split these segments: the draw functions below do not change.


#' @title Get the segments of a set of branches
#'
#' @description Describe each branch (identified by its tipward node) as a segment.
#'   Each branch is currently a single segment covering the whole branch (`from = 0`, `to = branch length`).
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param nodes Integer vector. Tipward node IDs of the branches. The order is kept in the output.
#'
#' @noRd
#'
#' @return A data.frame with one row per segment, in the order of `nodes`, and 3 columns:
#'  * `node` Integer. Tipward node ID of the branch the segment belongs to.
#'  * `from` Numeric. Distance from the tipward node to the tipward end of the segment.
#'  * `to` Numeric. Distance from the tipward node to the rootward end of the segment. The length of the segment is `to - from`.
#'

get_branch_segments <- function (phy, nodes)
{
  edge_ID <- match(nodes, phy$edge[, 2])
  data.frame(node = nodes, from = rep(0, length(nodes)), to = phy$edge.length[edge_ID])
}


#' @title Get the segments crossing a given age
#'
#' @description Keep the segments that contain a given age (age measured from the present, i.e., root age minus node height).
#'   The segments are kept in their input order.
#'
#' @param phy Phylogeny of class `"phylo"`.
#' @param segments Segments, as returned by `get_branch_segments()`.
#' @param age Numeric. Age (time before present).
#' @param root_age Numeric. Optional root age, to avoid recomputing node heights. Default is `NULL`.
#' @param edge_ages Numeric matrix. Optional ages (rootward, tipward) of the branches of `phy` in the order of `phy$edge`, i.e.,
#'   `root_age - phytools::nodeHeights(phy)`, to avoid recomputing them. Default is `NULL`.
#'
#' @importFrom phytools nodeHeights
#'
#' @noRd
#'
#' @return The rows of `segments` that cross `age` (same columns).
#'

get_segments_crossing_age <- function (phy, segments, age, root_age = NULL, edge_ages = NULL)
{
  if (is.null(edge_ages))
  {
    node_heights <- phytools::nodeHeights(phy)
    if (is.null(root_age)) { root_age <- max(node_heights[, 2]) }
    edge_ages <- root_age - node_heights
  }
  edge_ID <- match(segments$node, phy$edge[, 2])
  branch_length <- phy$edge.length[edge_ID]

  # Ages of the tipward and rootward ends of each segment. Written so that a segment covering the whole branch
  # gets exactly the ages of the branch (adding 0 and subtracting 0 leaves the ages unchanged)
  age_tipward <- edge_ages[edge_ID, 2] + segments$from
  age_rootward <- edge_ages[edge_ID, 1] - (branch_length - segments$to)

  segments[(age_rootward >= age) & (age_tipward <= age), , drop = FALSE]
}


#' @title Draw a point uniformly along segments
#'
#' @description Draw a segment with a probability proportional to its length, then a position uniformly within it.
#'   Therefore every time point along the segments has the same chance to be drawn, whatever the length of the segment it belongs to.
#'
#' @param segments Segments, as returned by `get_branch_segments()`.
#'
#' @noRd
#'
#' @return `NULL` if there is no segment. Otherwise a list with 2 elements:
#'  * `node` Integer. Tipward node ID of the branch of the selected segment.
#'  * `position` Numeric. Distance from the tipward node of this branch to the selected point (between `from` and `to` of the segment).
#'

draw_point_in_segments <- function (segments)
{
  if (nrow(segments) == 0) { return(NULL) }
  selected <- sample_one(seq_len(nrow(segments)), prob = segments$to - segments$from)
  list(node = segments$node[selected],
       position = stats::runif(n = 1, min = segments$from[selected], max = segments$to[selected]))
}


#' @title Length of the union of intervals
#'
#' @description Total length covered by a set of intervals, counting the parts covered by several intervals only once.
#'   This is used to evaluate the proportion of a time-window that is available for grafting based on
#'   the time interval covered by the union of the valid segments.
#'
#' @param from Numeric vector. Start of the intervals.
#' @param to Numeric vector. End of the intervals (same length as `from`).
#'
#' @noRd
#'
#' @return Numeric. The length of the union of the intervals (0 if there is none).
#'

get_union_length <- function (from, to)
{
  if (length(from) == 0) { return(0) }
  ord <- order(from)
  from <- from[ord]; to <- to[ord]
  # Each interval adds the part of it that is beyond the end of all the intervals starting before it
  covered_until <- c(-Inf, cummax(to)[-length(to)])
  sum(pmax(0, to - pmax(from, covered_until)))
}
