##### Function to check grafting constraints #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### 1/ Main function to check that grafted taxa respect their clade assignments ####

# Double-checks, a posteriori, that each grafted taxon sits in the right place in the grafted phylogeny:
#   - Target clade = clade defined by Target_clade (a pair of tips, or a single tip). The taxon must fall within this clade
#                       (crown), or on its stem branch if include_stem = TRUE (taxon then sister to the clade, with no other
#                       backbone tip outside the clade in between).
#                       With a single tip, the clade is the tip branch: the taxon must be sister to this tip
#                       (the same behaviour as graft_missing_taxa() / graft_from_divergence_age()).
#   - Excluded subclades = subclades (defined by pairs of tips) nested in the target clade where the taxon must NOT fall.
#                       As in the grafting functions, the stem branch of an excluded subclade is allowed
#                       (the taxon is then sister to the subclade, not within it).
# Also (optionally) checks that the backbone phylogeny is untouched by the grafting once the grafted taxa are pruned:
# topology and branch lengths are reported separately (branch lengths: max deviation, absolute and relative to backbone depth).
#
# Returns a data.frame with one row per taxon in missing_taxa_list (original order).
# The data.frame has attributes: "topology_preserved", "max_branch_length_deviation", "max_branch_length_deviation_relative",
# "backbone_preserved" (TRUE/FALSE/NA) and "all_ok" (TRUE/FALSE)

check_grafting_constraints <- function (phy_grafted, # Phylogeny after grafting
                                        phy_backbone, # Backbone phylogeny used for the grafting (before grafting)
                                        missing_taxa_list, # Same list as used for grafting. Items: Taxa_name, Target_clade, Excluded_subclades
                                        include_stem = TRUE, # Same as in graft_missing_taxa(): accept grafting on the stem branch of the target clade?
                                        use_excluded_subclades = TRUE, # Same as in graft_missing_taxa(): check the excluded subclades?
                                        check_backbone = TRUE, # Check that the backbone topology and branch lengths are preserved after pruning the grafted taxa
                                        tolerance = 1e-6, # Tolerance for the max deviation in backbone branch lengths, relative to the depth (root age) of the backbone
                                        stop_on_failure = FALSE, # Stop with an error if any constraint is violated
                                        verbose = TRUE) # Print a summary
{
  ## 0/ Prepare inputs ####

  taxa_names <- as.character(missing_taxa_list$Taxa_name)
  n_taxa <- length(taxa_names)

  excluded_list <- missing_taxa_list$Excluded_subclades
  if (use_excluded_subclades & is.null(excluded_list))
  {
    warning("No 'Excluded_subclades' item found in missing_taxa_list: excluded subclades are not checked")
    use_excluded_subclades <- FALSE
  }

  # Check list consistency
  if (length(missing_taxa_list$Target_clade) != n_taxa)
  {
    stop("'Target_clade' must have the same length as 'Taxa_name'")
  }
  if (use_excluded_subclades & length(excluded_list) != n_taxa)
  {
    stop("'Excluded_subclades' must have the same length as 'Taxa_name' (use NULL for taxa without excluded subclades)")
  }

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
  tips_of_node_bb <- function (node) { phy_backbone$tip.label[parts_bb[[node - n_tips_bb]]] }
  tips_of_node_g <- function (node) { parts_g[[node - n_tips_g]] } # Tip indices

  # Backbone tips found in the grafted tree
  is_backbone_tip_g <- phy_grafted$tip.label %in% phy_backbone$tip.label

  # Cache for the tips of the backbone clades
  clade_cache <- new.env()
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

  ## 1/ Check each taxon ####

  res <- data.frame(Taxon = taxa_names,
                    in_tree = taxa_names %in% phy_grafted$tip.label,
                    target_position = NA_character_, # "crown", "stem", "outside", or "tip_branch" (single-tip clade)
                    target_ok = NA,
                    excluded_ok = NA,
                    n_excluded_violated = NA_integer_,
                    excluded_violated = NA_character_, # Indices of the violated excluded subclades
                    ok = FALSE,
                    stringsAsFactors = FALSE)

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
    no_foreign_backbone_tip <- all(Y_backbone_tips %in% C_tips)

    if (length(target_clade_tips_i) == 1)
    {
      # Single-tip clade: taxon must be sister to the tip (no crown possible)
      res$target_position[i] <- ifelse(no_foreign_backbone_tip, "tip_branch", "outside")
      res$target_ok[i] <- no_foreign_backbone_tip
    } else {
      X_node <- ape::getMRCA(phy = phy_grafted, tip = target_clade_tips_i) # MRCA of the target clade in the grafted tree
      in_crown <- taxon_idx %in% tips_of_node_g(X_node)
      res$target_position[i] <- ifelse(in_crown, "crown", ifelse(no_foreign_backbone_tip, "stem", "outside"))
      res$target_ok[i] <- in_crown | (include_stem & no_foreign_backbone_tip)
    }

    ## Excluded subclades
    if (use_excluded_subclades)
    {
      excluded_i <- excluded_list[[i]]
      violated <- c()
      for (j in seq_along(excluded_i))
      {
        E_tips <- unique(unname(excluded_i[[j]]))
        if (length(E_tips) < 2) { next } # A single tip has no descendants: nothing to violate
        E_node <- ape::getMRCA(phy = phy_grafted, tip = E_tips)
        if (taxon_idx %in% tips_of_node_g(E_node)) { violated <- c(violated, j) }
      }
      res$n_excluded_violated[i] <- length(violated)
      res$excluded_violated[i] <- ifelse(length(violated) == 0, NA_character_, paste(violated, collapse = ","))
      res$excluded_ok[i] <- length(violated) == 0
    } else {
      res$excluded_ok[i] <- TRUE
      res$n_excluded_violated[i] <- 0L
    }

    res$ok[i] <- res$target_ok[i] & res$excluded_ok[i]
  }

  ## 2/ Check the backbone is preserved ####

  # Two separate results:
  #   - topology_preserved: same topology (and tip labels) once the grafted taxa are pruned
  #   - max_branch_length_deviation: maximum absolute difference in branch length between matched backbone branches,
  #     also given relative to the depth (root age) of the backbone. NA if the topology is not preserved.
  # Branch lengths are considered preserved if the relative deviation <= tolerance.
  topology_preserved <- NA
  max_dev_abs <- NA_real_
  max_dev_rel <- NA_real_
  lengths_preserved <- NA
  if (check_backbone)
  {
    pruned <- ape::keep.tip(phy = phy_grafted, tip = phy_backbone$tip.label)
    topology_preserved <- isTRUE(ape::all.equal.phylo(target = phy_backbone, current = pruned,
                                                      use.edge.length = FALSE, use.tip.label = TRUE))
    if (topology_preserved & !is.null(phy_backbone$edge.length) & !is.null(pruned$edge.length))
    {
      # Identify each branch by the set of backbone tips it subtends, to match branches between the two trees
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
      lengths_bb <- setNames(phy_backbone$edge.length, branch_keys(phy_backbone))
      lengths_pruned <- setNames(pruned$edge.length, branch_keys(pruned))
      max_dev_abs <- max(abs(lengths_pruned[names(lengths_bb)] - lengths_bb))
      backbone_depth <- max(ape::node.depth.edgelength(phy_backbone))
      max_dev_rel <- max_dev_abs / backbone_depth
      lengths_preserved <- max_dev_rel <= tolerance
    }
  }
  backbone_preserved <- if (check_backbone) isTRUE(topology_preserved) & isTRUE(lengths_preserved) else NA

  ## 3/ Summary ####

  all_ok <- all(res$ok) & !isFALSE(backbone_preserved)
  attr(res, "topology_preserved") <- topology_preserved
  attr(res, "max_branch_length_deviation") <- max_dev_abs
  attr(res, "max_branch_length_deviation_relative") <- max_dev_rel
  attr(res, "backbone_preserved") <- backbone_preserved
  attr(res, "all_ok") <- all_ok

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
               "  - Backbone: ", backbone_msg, "\n",
               "  => ", ifelse(all_ok, "All constraints respected\n", "Some constraints NOT respected (see returned data.frame and attributes)\n")))
  }
  if (stop_on_failure & !all_ok)
  {
    stop("Some grafting constraints are not respected: ", sum(!res$ok), " taxa failed; backbone preserved = ", backbone_preserved)
  }

  return(res)
}
