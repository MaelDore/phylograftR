##### Function to check grafting constraints #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Main function to check that grafted taxa respect their clade assignments ####

#' @title Check that grafted taxa respect their grafting constraints
#'
#' @description This function checks, a posteriori, that each taxon grafted with `graft_missing_taxa()` sits at a valid position in the grafted phylogeny:
#' within its target clade (or on its stem branch if allowed), and outside its excluded subclades.
#' It can also check that the backbone phylogeny is preserved by the grafting once the grafted taxa are pruned.
#'
#' @param phy_grafted An object of class `"phylo"`. The phylogeny after grafting.
#' @param phy_backbone An object of class `"phylo"`. The backbone phylogeny used for grafting (i.e., before grafting).
#' @param missing_taxa_list A list with the information on the grafted taxa. Typically the same list as used for grafting.
#'   Items: `$Taxon_name`, `$Target_clade`, and optionally `$include_target_clade_stem`, `$Excluded_subclades` and `$exclude_subclades_stem`.
#'   Every item provided must be named after `$Taxon_name` (same names, same order), as in [graft_missing_taxa()], see this function for details.
#' @param include_target_clade_stem Logical. Same as in [graft_missing_taxa()]: whether grafting on the stem branch of target clades is accepted.
#'   Taxa with a value in `missing_taxa_list$include_target_clade_stem` use it instead (`NA` or no item: this global value is used). Default is `TRUE`.
#' @param use_excluded_subclades Logical. Same as in [graft_missing_taxa()]: whether to check the excluded subclades.
#'   If no `$Excluded_subclades` item is found in `missing_taxa_list`, a warning is issued and excluded subclades are not checked. Default is `TRUE`.
#' @param exclude_subclades_stem Logical. Same as in [graft_missing_taxa()]: whether the stem of the excluded subclades is also excluded
#'   (a taxon sister to an excluded subclade is then a violation). Excluded subclades with a value in `missing_taxa_list$exclude_subclades_stem` use it instead
#'   (`NA` or no item: this global value is used). Default is `TRUE`.
#' @param check_backbone Logical. Whether to check that the topology and branch lengths of the backbone are preserved once the grafted taxa are pruned. Default is `TRUE`.
#' @param tolerance Numeric. Tolerance for the maximum deviation in backbone branch lengths, relative to the depth (root age) of the backbone.
#'   Default is `1e-6`.
#' @param stop_on_failure Logical. Whether to stop with an error if any constraint is violated. Default is `FALSE`.
#' @param verbose Logical. Whether to print a summary of the checks. Default is `TRUE`.
#'
#' @return The function returns a data.frame with one row per taxon in `missing_taxa_list` (in the original order), and the following columns:
#'
#' * `$Taxon_name`: Name of the taxon.
#' * `$in_tree`: Logical. Whether the taxon is found in the grafted phylogeny.
#' * `$include_target_clade_stem`: Logical. Whether the stem of the target clade was accepted for this taxon (`missing_taxa_list$include_target_clade_stem` of the taxon if provided, otherwise the global `include_target_clade_stem` argument).
#' * `$target_position`: Character string. Position of the taxon relative to its target clade: `"crown"` (within the clade), `"stem"` (sister to the clade, or to other missing taxa grafted on the stem),
#'      `"tip_branch"` (on the branch of the tip, for single-tip target clades), or `"outside"` (neither within the clade nor on its stem = an error).
#' * `$target_ok`: Logical. Whether the position respects the target clade (`"crown"` or, if the stem is accepted for the taxon, `"stem"`; `"tip_branch"` for single-tip clades).
#' * `$excluded_ok`: Logical. Whether the taxon is outside all its excluded subclades (and not on their stem, for the subclades with an excluded stem).
#' * `$n_excluded_violated`: Integer. Number of excluded subclades violated.
#' * `$excluded_violated`: Character string. Names of the excluded subclades violated (or their indices if the subclades are not named), separated by commas (`NA` if none).
#' * `$ok`: Logical. Whether the taxon respects both its target clade and its excluded subclades.
#'
#' The data.frame has the following attributes:
#'
#' * `"topology_preserved"`: Logical. Whether the backbone topology is preserved once the grafted taxa are pruned (`NA` if `check_backbone = FALSE`).
#' * `"max_branch_length_deviation"`: Numeric. Maximum absolute deviation in backbone branch lengths (`NA` if the topology is not preserved).
#' * `"max_branch_length_deviation_relative"`: Numeric. Same, relative to the depth of the backbone.
#' * `"backbone_preserved"`: Logical. Whether both topology and branch lengths are preserved, within `tolerance` (`NA` if `check_backbone = FALSE`).
#' * `"all_ok"`: Logical. Whether all taxa respect their constraints and the backbone is preserved (if checked).
#'
#' Attributes are accessed with `attr(x = output, which = "attribute_name")`.
#'
#'
#' @details
#' A taxon respects its target clade if it is found within the clade, or on its stem branch if the stem is accepted for this taxon:
#' global `include_target_clade_stem = TRUE`, or `missing_taxa_list$include_target_clade_stem` of this taxon
#' (the taxon is then sister to the target clade, or to other missing taxa sister to the target clade, with no other backbone tip in between).
#' For a single-tip target clade, the taxon must be on the branch of the tip, with no other backbone tip in between.
#'
#' The 'stem' position is not restricted to a taxon being directly sister to the target clade: taxa grafted on the stem can receive further grafts,
#' so several missing taxa can form together a group sister to the target clade (i.e, a missing-only subtree grafted on the the stem).
#' Only backbone tips are considered when checking the position of a taxon (missing taxa are ignored), so such configurations are valid.
#'
#' A taxon respects its excluded subclades if it is not nested within any of them.
#' If the stem of an excluded subclade is excluded (`exclude_subclades_stem = TRUE`, the default), the taxon must also not be on its stem, i.e., sister to the subclade
#' (or to other missing taxa sister to the subclade), with no other backbone tip in between. For a single-tip excluded subclade, the taxon must not be on the branch of the tip.
#' As in [graft_missing_taxa()], if the stem of the subclade is not excluded, the taxon can be sister to the subclade, not within it.
#'
#' Branches of the backbone are matched between the pruned grafted phylogeny and the backbone by the set of tips they subtend.
#' Branch lengths are considered preserved if the maximum deviation, relative to the depth of the backbone, is lower than `tolerance`.
#'
#' @importFrom ape all.equal.phylo getMRCA keep.tip node.depth.edgelength prop.part
#'
#' @importFrom stats setNames
#'
#' @export
#'
#' @seealso Grafting functions:
#'   * [graft_missing_taxa()] to graft taxa within a target clade with constraints (target clades and excluded subclades),
#'   * [graft_missing_clades()] to graft groups of missing taxa as clades sister to a clade of the backbone.
#'
#' @author Maël Doré
#'
#' @examples
#' # Create a toy backbone phylogeny: ultrametric, with a root age of 50 time units
#' set.seed(1234)
#' phy <- ape::rcoal(n = 30)
#' phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * 50
#'
#' # List the missing taxa to graft
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   # Items are named after Taxon_name (same names, same order)
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")),
#'   Excluded_subclades = list(X1 = list(c("t10", "t19")), X2 = NULL, X3 = NULL, X4 = NULL))
#'
#' # Graft the missing taxa
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   verbose = FALSE, seed = 1234)
#'
#' # Check that the grafted taxa respect their constraints, and that the backbone is preserved
#' check_results <- check_grafting_constraints(phy_grafted = grafted_phy, phy_backbone = phy,
#'                                             missing_taxa_list = missing_taxa_list)
#' check_results
#' attr(check_results, "all_ok")
#'


check_grafting_constraints <- function (phy_grafted, # Phylogeny after grafting
                                        phy_backbone, # Backbone phylogeny used for the grafting (before grafting)
                                        missing_taxa_list, # Same list as used for grafting. Items: Taxon_name, Target_clade, include_target_clade_stem, Excluded_subclades, exclude_subclades_stem
                                        include_target_clade_stem = TRUE, # Same as in graft_missing_taxa(): accept grafting on the stem branch of the target clade?
                                        use_excluded_subclades = TRUE, # Same as in graft_missing_taxa(): check the excluded subclades?
                                        exclude_subclades_stem = TRUE, # Same as in graft_missing_taxa(): is the stem of the excluded subclades also excluded?
                                        check_backbone = TRUE, # Check that the backbone topology and branch lengths are preserved after pruning the grafted taxa
                                        tolerance = 1e-6, # Tolerance for the max deviation in backbone branch lengths, relative to the depth (root age) of the backbone
                                        stop_on_failure = FALSE, # Stop with an error if any constraint is violated
                                        verbose = TRUE) # Print a summary
{
  ## 1.1/ Prepare inputs ####

  # Check the structure of missing_taxa_list (Taxon_name, and every item provided with the same length and names as Taxon_name), and extract taxa
  taxa_names <- validate_missing_taxa_list(missing_taxa_list, clade_item = "Target_clade")
  n_taxa <- length(taxa_names)

  # Stem availability of each taxon: 'missing_taxa_list$include_target_clade_stem' of the taxon if provided, otherwise the global 'include_target_clade_stem'
  include_target_clade_stem_taxa <- resolve_include_target_clade_stem(include_target_clade_stem_taxa = missing_taxa_list$include_target_clade_stem,
                                                      include_target_clade_stem = include_target_clade_stem, n_taxa = n_taxa,
                                                      taxon_names = taxa_names)

  # Check presence of excluded subclades
  excluded_list <- missing_taxa_list$Excluded_subclades
  if (use_excluded_subclades & is.null(excluded_list))
  {
    warning("No 'Excluded_subclades' item found in missing_taxa_list: excluded subclades are not checked")
    use_excluded_subclades <- FALSE
  }

  # Stem exclusion of each excluded subclade: 'missing_taxa_list$exclude_subclades_stem' of the subclade if provided, otherwise the global 'exclude_subclades_stem'
  exclude_stem_list <- resolve_exclude_subclades_stem(exclude_subclades_stem_taxa = missing_taxa_list$exclude_subclades_stem,
                                                      exclude_subclades_stem = exclude_subclades_stem,
                                                      excluded_subclades = if (use_excluded_subclades) excluded_list else NULL, n_taxa = n_taxa,
                                                      taxon_names = taxa_names)

  # Check tips used to define the clades are in the backbone
  defining_tips <- unique(unlist(c(missing_taxa_list$Target_clade, if (use_excluded_subclades) excluded_list)))
  if (!all(defining_tips %in% phy_backbone$tip.label))
  {
    stop(paste0("Some tip labels used to define clades are missing from the backbone phylogeny\n",
                "Missing tips: ", paste(setdiff(defining_tips, phy_backbone$tip.label), collapse = ", ")))
  }

  # Check the backbone tips are all in the grafted phylogeny
  if (!all(phy_backbone$tip.label %in% phy_grafted$tip.label))
  {
    stop("Some backbone tips are missing from the grafted phylogeny")
  }

  # Precompute tips descending from each internal node in both trees
  n_tips_bb <- length(phy_backbone$tip.label)
  n_tips_g <- length(phy_grafted$tip.label)
  parts_bb <- ape::prop.part(phy_backbone, check.labels = FALSE)
  parts_g <- ape::prop.part(phy_grafted, check.labels = FALSE)

  # Functions to extract descending tips for a given node in both the backbone and grafted trees
  tips_of_node_bb <- function (node) { phy_backbone$tip.label[parts_bb[[node - n_tips_bb]]] }
  tips_of_node_g <- function (node) { parts_g[[node - n_tips_g]] }

  # Backbone tips found in the grafted tree
  is_backbone_tip_g <- phy_grafted$tip.label %in% phy_backbone$tip.label

  # Cache for the tips of the backbone clades
  clade_cache <- new.env()
  # Function to retrieve tips of a given node in the backbone trees, but using a cache system
  # to not have to recompute it for each taxa if it was already done for its associated clade
  get_backbone_clade_tips <- function (tips)
  {
    tips <- unique(tips)
    if (length(tips) == 1) { return(tips) }
    key <- paste(sort(tips), collapse = "|")
    if (is.null(clade_cache[[key]]))
    {
      clade_cache[[key]] <- tips_of_node_bb(ape::getMRCA(phy = phy_backbone, tip = tips))
    }
    return(clade_cache[[key]])
  }

  ## 1.2/ Check each taxon ####

  # Initiate resulting data.frame
  res <- data.frame(Taxon_name = taxa_names,
                    in_tree = taxa_names %in% phy_grafted$tip.label,
                    include_target_clade_stem = include_target_clade_stem_taxa,
                    target_position = NA, # "crown", "stem", "outside", or "tip_branch" (single-tip clade)
                    target_ok = NA,
                    excluded_ok = NA,
                    n_excluded_violated = NA,
                    excluded_violated = NA, # Indices of the violated excluded subclades
                    ok = FALSE,
                    stringsAsFactors = FALSE)

  ## Loop per taxa
  for (i in seq_len(n_taxa))
  {
    # i <- 1
    if (!res$in_tree[i]) { next }

    taxon_idx <- match(taxa_names[i], phy_grafted$tip.label)
    target_clade_tips_i <- unique(unname(missing_taxa_list$Target_clade[[i]]))

    ## Target clade
    C_tips <- get_backbone_clade_tips(target_clade_tips_i) # Backbone tips of the target clade

    # Smallest clade of the grafted tree containing the target clade and the taxon
    Y_node <- ape::getMRCA(phy = phy_grafted, tip = c(target_clade_tips_i, taxa_names[i]))
    Y_tips_idx <- tips_of_node_g(Y_node)
    Y_backbone_tips <- phy_grafted$tip.label[Y_tips_idx[is_backbone_tip_g[Y_tips_idx]]]
    # Check that no foreign tip (outside of the backbone target clade) are present in the associated clade in the grafted tree
    no_foreign_backbone_tip <- all(Y_backbone_tips %in% C_tips)

    ## Deal with the different cases
    if (length(target_clade_tips_i) == 1)
    {
      ## Case with a single tip

      # Single-tip clade: taxon must be on the tip branch, with no other backbone tip in between (no crown possible)
      res$target_position[i] <- ifelse(no_foreign_backbone_tip, "tip_branch", "outside")
      res$target_ok[i] <- no_foreign_backbone_tip
    } else {
      ## Case with a clade (defined by two tips)

      # Get MRCA of the target clade in the grafted tree
      X_node <- ape::getMRCA(phy = phy_grafted, tip = target_clade_tips_i)
      # Check that the taxa is found within the target clade on the grafted phylogeny
      in_crown <- taxon_idx %in% tips_of_node_g(X_node)
      # If not, use 'no_foreign_backbone_tip' to distinguish between taxa grafted on the stem branch, and taxa outside of the target clade
      res$target_position[i] <- ifelse(in_crown, "crown", ifelse(no_foreign_backbone_tip, "stem", "outside"))
      # Target is okay if in_crown or on the stem when the stem is accepted for this taxon
      res$target_ok[i] <- in_crown | (include_target_clade_stem_taxa[i] & no_foreign_backbone_tip)
    }

    ## Excluded subclades
    if (use_excluded_subclades)
    {
      excluded_i <- excluded_list[[i]]
      stem_excluded_i <- exclude_stem_list[[i]] # Is the stem of each excluded subclade also excluded?
      violated <- c()
      for (j in seq_along(excluded_i))
      {
        E_tips <- unique(unname(excluded_i[[j]]))
        if (length(E_tips) == 0) { next }

        # Case with only a single tip as Excluded subclade: it has no descendants, only its terminal branch (its stem) can be excluded
        if (length(E_tips) == 1)
        {
          in_E_crown <- FALSE
        } else {
          # Case with a clade (pairs of tips) as Excluded subclade(s)
          E_node <- ape::getMRCA(phy = phy_grafted, tip = E_tips)
          # Check that the clade containing the target taxon on the grafted tree does not include tips from the excluded subclade(s)
          in_E_crown <- taxon_idx %in% tips_of_node_g(E_node)
        }
        if (in_E_crown) { violated <- c(violated, j); next }

        # If the stem of the excluded subclade is excluded: the taxon must not be on its stem, i.e., sister to the subclade
        # (or to other missing taxa sister to it), with no other backbone tip in between
        if (stem_excluded_i[j])
        {
          Z_node <- ape::getMRCA(phy = phy_grafted, tip = c(E_tips, taxa_names[i]))
          Z_tips_idx <- tips_of_node_g(Z_node)
          Z_backbone_tips <- phy_grafted$tip.label[Z_tips_idx[is_backbone_tip_g[Z_tips_idx]]]
          if (all(Z_backbone_tips %in% get_backbone_clade_tips(E_tips))) { violated <- c(violated, j) }
        }
      }
      # Record the number of excluded subclades violated
      res$n_excluded_violated[i] <- length(violated)
      # Record the excluded subclades violated: their names if they are named, otherwise their indices
      subclade_labels <- names(excluded_i)
      violated_labels <- if (is.null(subclade_labels)) as.character(violated) else ifelse(is.na(subclade_labels[violated]) | subclade_labels[violated] == "", as.character(violated), subclade_labels[violated])
      res$excluded_violated[i] <- ifelse(length(violated) == 0, NA, paste(violated_labels, collapse = ","))
      # Taxon is ok if no excluded subclades are violated
      res$excluded_ok[i] <- length(violated) == 0
    } else {
      res$excluded_ok[i] <- TRUE
      res$n_excluded_violated[i] <- 0L
    }

    # Overall, a taxon is okay if both rules for target clade and excluded subclades are respected
    res$ok[i] <- res$target_ok[i] & res$excluded_ok[i]
  }

  ## 1.3/ Check the backbone is preserved ####

  # Two separate results:
  #   - topology_preserved: same topology (and tip labels) once the grafted taxa are pruned
  #   - max_branch_length_deviation: maximum absolute difference in branch length between matched backbone branches,
  #     also given relative to the depth (root age) of the backbone. NA if the topology is not preserved.
  # Branch lengths are considered preserved if the relative deviation <= tolerance.

  # Initiate variables
  topology_preserved <- NA
  max_dev_abs <- NA
  max_dev_rel <- NA
  lengths_preserved <- NA

  # Check backbone only if requested
  if (check_backbone)
  {
    # Prune the grafted tree from the grafted tips
    pruned <- ape::keep.tip(phy = phy_grafted, tip = phy_backbone$tip.label)
    # Compare the topology of the pruned tree to the initial backbone tree
    topology_preserved <- isTRUE(ape::all.equal.phylo(target = phy_backbone, current = pruned,
                                                      use.edge.length = FALSE, use.tip.label = TRUE))

    # if the topology match, check branch lengths
    if (topology_preserved & !is.null(phy_backbone$edge.length) & !is.null(pruned$edge.length))
    {
      # Function to identify each branch by the set of backbone tips it subtends, to match branches between the two trees
      branch_keys <- function (tree)
      {
        n <- length(tree$tip.label)
        parts <- ape::prop.part(tree, check.labels = FALSE)
        node_keys <- character(n + tree$Nnode)
        node_keys[1:n] <- as.character(match(tree$tip.label, phy_backbone$tip.label))
        for (k in seq_along(parts))
        {
          node_keys[n + k] <- paste(sort(match(tree$tip.label[parts[[k]]], phy_backbone$tip.label)), collapse = ",")
        }
        return(node_keys[tree$edge[, 2]])
      }
      # Extract matching branch lengths
      lengths_bb <- stats::setNames(phy_backbone$edge.length, branch_keys(phy_backbone))
      lengths_pruned <- stats::setNames(pruned$edge.length, branch_keys(pruned))
      # Compute maximum deviation in length across all matching branches
      max_dev_abs <- max(abs(lengths_pruned[names(lengths_bb)] - lengths_bb))
      # Adjust deviation to the depth of the tree = relative deviation
      backbone_depth <- max(ape::node.depth.edgelength(phy_backbone))
      max_dev_rel <- max_dev_abs / backbone_depth
      # Compare maximum relative deviation to tolerance as the check for preserved branch lengths
      lengths_preserved <- max_dev_rel <= tolerance
    }
  }
  # backbone is fully preserved if both topology and branch lengths are matching (within tolerance)
  backbone_preserved <- if (check_backbone) isTRUE(topology_preserved) & isTRUE(lengths_preserved) else NA

  ## 1.4/ Summary ####

  # Summarize check results as attributes
  all_ok <- all(res$ok) & !isFALSE(backbone_preserved)
  attr(res, "topology_preserved") <- topology_preserved
  attr(res, "max_branch_length_deviation") <- max_dev_abs
  attr(res, "max_branch_length_deviation_relative") <- max_dev_rel
  attr(res, "backbone_preserved") <- backbone_preserved
  attr(res, "all_ok") <- all_ok

  # Describe the taxa that failed (max 10), for the summary and the error message
  describe_failures <- function ()
  {
    failed <- which(!res$ok | !res$in_tree)
    if (length(failed) == 0) { return("") }
    n_show <- min(length(failed), 10)
    reason <- vapply(failed[seq_len(n_show)], function (i)
    {
      if (!res$in_tree[i]) { return("missing from the grafted tree") }
      paste(c(if (isFALSE(res$target_ok[i])) paste0("outside its target clade (position: ", res$target_position[i], ")"),
              if (isFALSE(res$excluded_ok[i])) paste0("within or on the stem of excluded subclade(s): ", res$excluded_violated[i])), collapse = "; ")
    }, character(1))
    paste0("\n  Failed taxa:\n", paste0("    - ", taxa_names[failed[seq_len(n_show)]], ": ", reason, collapse = "\n"),
           if (length(failed) > n_show) paste0("\n    ... and ", length(failed) - n_show, " more") else "")
  }

  # Display check results if requested
  if (verbose)
  {
    backbone_msg <- if (!check_backbone) { "not checked" } else if (!topology_preserved) {
      "topology NOT preserved"
    } else {
      paste0("topology preserved, branch lengths shifted by max ", format(max_dev_abs, digits = 3),
             " (", format(max_dev_rel, digits = 3), " of backbone depth; tolerance = ", format(tolerance), ")")
    }
    cat(paste0(n_taxa, " taxa checked\n",
               "  - Missing from the grafted tree: ", sum(!res$in_tree), "\n",
               "  - Target clade violated:          ", sum(res$target_ok %in% FALSE), "\n",
               "  - Excluded subclade violated:     ", sum(res$excluded_ok %in% FALSE), "\n",
               "  - Backbone: ", backbone_msg, if (!all(res$ok)) paste0(describe_failures(), "\n") else "", "\n",
               "  => ", ifelse(all_ok, "All constraints respected\n", "Some constraints NOT respected (see returned data.frame and attributes)\n")))
  }

  # Display error message
  if (stop_on_failure & !all_ok)
  {
    stop("Some grafting constraints are not respected: ", sum(!res$ok), " taxa failed; backbone preserved = ", backbone_preserved, describe_failures())
  }

  return(res)
}
