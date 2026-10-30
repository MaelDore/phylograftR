##### Function to extract the grafting information from grafted phylogenies #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Extract the grafting information ####

#' @title Extract the grafting information from grafted phylogenies
#'
#' @description Summarize, in a table, where and when each grafted taxon (or missing clade) was grafted on the backbone:
#'   the age of the graft, the backbone branch it was grafted on (the parent clade), and, with a `simmap`, the state at the graft point and the state inherited by the grafted lineage.
#'   The function works a posteriori, on the output of [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] or [`phylograftR::graft_missing_clades()`][graft_missing_clades]
#'   (a phylogeny, a set of simulations, or the list of maps obtained with a `"multiSimmap"`).
#'
#'   The table is meant to check the plausibility of the grafting against the backbone, and, over a set of phylogenies, to describe the uncertainty of the grafting
#'   (e.g., the distribution of the graft ages and graft states of each taxon).
#'
#' @param phy_grafted A grafted phylogeny (class `"phylo"`), a set of grafted phylogenies (class `"multiPhylo"`, e.g., several simulations),
#'   or a list of them with one element per map (output of the grafting functions with a `"multiSimmap"`).
#' @param missing_taxa_list The list of missing taxa used for the grafting. Items: `$Taxon_name`, and `$Target_clade` (taxa, grafted with
#'   [`phylograftR::graft_missing_taxa()`][graft_missing_taxa]) or `$Sister_clade` (missing clades, grafted with [`phylograftR::graft_missing_clades()`][graft_missing_clades]).
#'   Exactly one of them must be provided: it defines the type of the table.
#' @param phy_backbone Optional. The backbone phylogeny used for the grafting (class `"phylo"`). If `NULL` (default), the backbone is rebuilt by removing the grafted taxa from
#'   the grafted phylogeny (same topology and branch lengths). The node IDs of the backbone (`parent_node`) can only be given if `phy_backbone` is provided.
#'
#' @details
#' ### What is extracted depends on the type of grafting
#'
#' **Taxa** (`$Target_clade` list): one row per taxon.
#' A taxon can be grafted on a branch of the backbone, or on a branch of a taxon grafted before it (and so on):
#' the grafted taxa that hang together form a *grafted subtree* attached to the backbone at a *graft point*.
#' The summary table reports:
#' * `graft_age`: the age of the graft point, where the grafted subtree is attached to a branch of the backbone;
#' * `divergence_age`: the age of the node where the taxon splits from its sister lineage. It is equal to `graft_age` for a taxon attached directly to the backbone,
#'   and younger for a taxon attached to another grafted taxon (`grafted_on` then lists the grafted taxa of its sister lineage).
#' * The parent clade is the clade of the backbone subtended by the branch of the graft point: `parent_tips` gives two tips defining it (one if the branch is terminal),
#'   and `parent_node` its node ID in `phy_backbone`.
#'
#' **Missing clades** (`$Sister_clade` list): one row per missing clade (the taxa with the same sister clade, grafted together).
#' The missing clade is attached to the stem branch of its sister clade at a *graft point*.
#' The summary table reports:
#' * `graft_age`: the age of the graft point;
#' * `crown_age`: the age of the most recent common ancestor of the taxa of the clade (`NA` for a single taxon);
#' * The parent clade is the backbone clade subtended by the branch it was attached to (the sister clade, if the constraints were respected, see [`phylograftR::check_grafting_constraints()`][check_grafting_constraints]).
#' If the taxa of a missing clade do not form a clade in the phylogeny, the values of the row are `NA`, with a warning.
#'
#' ### States
#'
#' If the phylogeny has a `$simmap` (output of the grafting with a `simmap` and `return_simmap = TRUE`), the states are extracted:
#' `parent_state` is the state of the branch of the backbone at the graft point, and `tip_state` (taxa) or `clade_state` (missing clades) the state inherited by the grafted lineage.
#'
#' They are different if `treat_states_as_ranges = TRUE`: a taxon grafted on a branch in the range "AB" can inherit the range "A".
#' With the maps of several traits (`$simmap` is a list), there is one column per trait (`parent_state_<trait>`, `tip_state_<trait>`).
#' The states are `NA` for phylogenies without `$simmap`, or whose `$simmap` does not match the phylogeny (e.g., after removing tips).
#'
#' @importFrom ape Ntip reorder.phylo node.depth.edgelength drop.tip getMRCA
#'
#' @export
#'
#' @seealso [`phylograftR::check_grafting_plan()`][check_grafting_plan] to check the plan before grafting,
#'   [`phylograftR::check_grafting_constraints()`][check_grafting_constraints] to check the constraints after grafting,
#'   and [`phylograftR::identify_grafted()`][identify_grafted] to identify the grafted tips and branches.
#'
#' @return A data.frame, with one row per taxon (or per missing clade) and phylogeny.
#'   The columns are:
#'   * `map` Character string. Name of the stochastic map (`NA` if no map was used).
#'   * `simulation` Integer. Index of the simulation.
#'   * `seed` Integer. Seed of the simulation (`NA` if unknown).
#'
#'   For taxa:
#'   * `Taxon_name` Character string. Name of the taxon.
#'   * `parent_node` Integer. Node ID in `phy_backbone` of the tipward end of the backbone branch of the graft point (`NA` without `phy_backbone`).
#'   * `parent_tips` Character string. Tips defining the parent clade.
#'   * `graft_age` Numeric. Age of the graft point.
#'   * `divergence_age` Numeric. Age of the node where the taxon splits from its sister lineage.
#'   * `grafted_on` Character string. Taxa grafted on the sister lineage (first five), `NA` if the taxon is attached directly to the backbone.
#'
#'   For missing clades:
#'   * `Missing_clade` Integer. Clade ID.
#'   * `Taxa` Character string. Taxa of the missing clade, separated by commas.
#'   * `n_taxa` Integer. Number of taxa in the grafted clade.
#'   * `parent_node` Integer. Node ID in `phy_backbone` of the tipward end of the backbone branch of the graft point (`NA` without `phy_backbone`).
#'   * `parent_tips` Character string. Tips defining the parent clade.
#'   * `graft_age` Numeric. Age of the graft point = stem age of the grafted clade.
#'   * `crown_age` Numeric. Crown age of the grafted clade.
#'
#'   With a `simmap`: `parent_state` and `tip_state` (taxa) or `clade_state` (missing clades), or one column per trait.
#'
#'   Attributes: `"type"` (`"taxa"` or `"clades"`).
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
#' # List the missing taxa to graft
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3", "X4"),
#'   Target_clade = list(X1 = c("t26", "t18"), X2 = c("t8", "t22"), X3 = "t5", X4 = c("t26", "t18")))
#'
#' # --- Graft the taxa (several simulations) --- #
#'
#' grafted_phy <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                   nb_simulations = 5, verbose = FALSE, seed = 1234)
#'
#' # --- Extract the grafting information --- #
#'
#' info <- extract_grafting_info(phy_grafted = grafted_phy, missing_taxa_list = missing_taxa_list,
#'                               phy_backbone = phy)
#' head(info)
#'
#' # Inspect uncertainty of the graft age of each taxon over the simulations
#' aggregate(graft_age ~ Taxon_name, data = info, FUN = function (x) round(range(x), 1))
#'
#'

extract_grafting_info <- function (phy_grafted, # Output of the grafting functions: phylo, multiPhylo, or list of maps
                                   missing_taxa_list, # Same list as used for grafting. Items: Taxon_name, Target_clade (or Sister_clade)
                                   phy_backbone = NULL) # Backbone phylogeny used for the grafting. NULL: rebuilt by removing the grafted taxa
{
  ## Type of the table: taxa ('Target_clade') or missing clades ('Sister_clade')
  if (!is.list(missing_taxa_list)) { stop("'missing_taxa_list' must be a list") }
  has_target <- !is.null(missing_taxa_list$Target_clade)
  has_sister <- !is.null(missing_taxa_list$Sister_clade)
  if (has_target && has_sister) { stop("'missing_taxa_list' has both 'Target_clade' and 'Sister_clade': use only one, depending on the grafting function (graft_missing_taxa() or graft_missing_clades())") }
  if (!has_target && !has_sister) { stop("'missing_taxa_list' must have a 'Target_clade' item (taxa grafted by graft_missing_taxa()) or a 'Sister_clade' item (missing clades grafted by graft_missing_clades())") }
  is_clades <- has_sister
  validate_missing_taxa_list(missing_taxa_list, clade_item = if (is_clades) "Sister_clade" else "Target_clade", optional_items = character(0))
  taxa <- missing_taxa_list$Taxon_name

  if (!is.null(phy_backbone) && !inherits(phy_backbone, "phylo")) { stop("'phy_backbone' must be NULL or an object of class 'phylo'") }

  ## List the phylogenies
  listed <- list_grafted_phylogenies(phy_grafted, arg = "phy_grafted")

  ## Backbone: provided, or rebuilt from the first phylogeny (all the phylogenies share the topology of the backbone)
  rebuild_backbone <- function (phy) { ape::drop.tip(phy, intersect(taxa, phy$tip.label)) }
  if (!is.null(phy_backbone))
  {
    expected <- setdiff(listed$phylos[[1]]$tip.label, taxa)
    if (!setequal(expected, phy_backbone$tip.label))
    {
      stop("The tips of 'phy_backbone' differ from the tips of 'phy_grafted' without the grafted taxa (",
           length(setdiff(expected, phy_backbone$tip.label)), " tips of the grafted phylogeny are absent from the backbone, ",
           length(setdiff(phy_backbone$tip.label, expected)), " tips of the backbone are absent from the grafted phylogeny).")
    }
  }
  backbone <- if (is.null(phy_backbone)) rebuild_backbone(listed$phylos[[1]]) else phy_backbone

  ## Missing clades: groups of taxa with the same sister clade
  groups <- NULL
  if (is_clades)
  {
    absent_tips <- setdiff(unique(unlist(missing_taxa_list$Sister_clade)), backbone$tip.label)
    if (length(absent_tips) > 0)
    {
      stop("Some tip labels provided to locate the sister clades are missing from the backbone phylogeny\nMissing tips: ", paste(absent_tips, collapse = ", "))
    }
    groups <- group_and_order_clades(backbone, missing_taxa_list$Sister_clade)
  }

  ## Extract the information of each phylogeny (the same warnings are sent once)
  results <- with_unique_warnings(lapply(listed$phylos, function (phy)
  {
    extract_grafting_info_one(phy = phy, taxa = taxa, backbone = backbone, is_clades = is_clades, groups = groups, node_ids = !is.null(phy_backbone))
  }))

  ## Long table
  n_rows <- vapply(results, nrow, integer(1))
  output <- do.call(rbind, lapply(seq_along(results), function (k)
  {
    cbind(map = rep(listed$map[k], n_rows[k]), simulation = rep(listed$simulation[k], n_rows[k]), seed = rep(listed$seed[k], n_rows[k]), results[[k]], stringsAsFactors = FALSE)
  }))
  rownames(output) <- NULL
  attr(output, "type") <- if (is_clades) "clades" else "taxa"
  output
}


#' @title Extract the grafting information of one grafted phylogeny
#'
#' @param phy Grafted phylogeny (class `"phylo"`), possibly with a `$simmap`.
#' @param taxa Character vector. Names of the grafted taxa (`missing_taxa_list$Taxon_name`).
#' @param backbone Backbone phylogeny (provided by the user, or rebuilt).
#' @param is_clades Logical. Missing clades (rows per missing clade) or taxa (rows per taxon).
#' @param groups Output of `group_and_order_clades()` for the sister clades (missing clades only).
#' @param node_ids Logical. Is `backbone` provided by the user? (Its node IDs are returned.)
#'
#' @noRd
#'
#' @return A data.frame with the columns of `extract_grafting_info()` except `map`, `simulation` and `seed`.
#'

extract_grafting_info_one <- function (phy, taxa, backbone, is_clades, groups, node_ids)
{
  n_tips <- ape::Ntip(phy)
  n_nodes <- n_tips + phy$Nnode
  edge <- phy$edge

  ## Structure of the tree
  parent <- integer(n_nodes); parent[edge[, 2]] <- edge[, 1]
  children <- split(edge[, 2], factor(edge[, 1], levels = seq_len(n_nodes)))
  age <- { depth <- ape::node.depth.edgelength(phy); max(depth[seq_len(n_tips)]) - depth }
  is_grafted <- phy$tip.label %in% taxa
  # Number of tips and of backbone tips (not grafted) below each node (tips included)
  n_below <- c(rep(1, n_tips), rep(0, phy$Nnode))
  n_backbone <- c(as.numeric(!is_grafted), rep(0, phy$Nnode))
  post <- ape::reorder.phylo(phy, order = "postorder")$edge
  for (i in seq_len(nrow(post)))
  {
    n_below[post[i, 1]] <- n_below[post[i, 1]] + n_below[post[i, 2]]
    n_backbone[post[i, 1]] <- n_backbone[post[i, 1]] + n_backbone[post[i, 2]]
  }
  tips_by_node <- NULL # Tips below each internal node, computed if needed
  tips_below <- function (node)
  {
    if (node <= n_tips) { return(node) }
    if (is.null(tips_by_node)) { tips_by_node <<- ape::prop.part(phy) }
    tips_by_node[[node - n_tips]]
  }
  other_child <- function (node, child) { setdiff(children[[node]], child)[1] }

  ## Missing from the phylogeny
  absent <- setdiff(taxa, phy$tip.label)
  if (length(absent) > 0)
  {
    warning(length(absent), " taxa of 'missing_taxa_list' are not tips of the grafted phylogeny: their values are NA. ",
            paste(utils::head(absent, 10), collapse = ", "), if (length(absent) > 10) ", ..." else "", call. = FALSE)
  }

  ## States of the branches: the simmap of the phylogeny (a simmap, or a list of simmaps, one per trait), if it matches the phylogeny
  simmaps <- NULL
  if (!is.null(phy$simmap))
  {
    candidates <- if (inherits(phy$simmap, "simmap")) list(phy$simmap) else phy$simmap
    if (!is.null(names(phy$simmap)) && !inherits(phy$simmap, "simmap")) { names(candidates) <- names(phy$simmap) }
    matching <- vapply(candidates, function (s) inherits(s, "simmap") && identical(s$edge, edge), logical(1))
    if (all(matching)) { simmaps <- candidates } else { warning("The '$simmap' of a phylogeny does not match the phylogeny (edges differ): the states are NA.", call. = FALSE) }
  }
  edge_row <- integer(n_nodes); edge_row[edge[, 2]] <- seq_len(nrow(edge))
  state_of_edge <- function (smap, node) # State at the tipward end of the branch leading to the node
  {
    if (is.na(node) || edge_row[node] == 0) { return(NA_character_) }
    map <- smap$maps[[edge_row[node]]]
    names(map)[length(map)]
  }

  ## Parent clade: the clade of the backbone that the lineage 'c' of the graft point continues (its backbone tips define the clade)
  cache <- new.env()
  parent_clade <- function (lineage)
  {
    key <- as.character(lineage)
    if (is.null(cache[[key]]))
    {
      backbone_tips <- phy$tip.label[tips_below(lineage)]
      backbone_tips <- backbone_tips[!(backbone_tips %in% taxa)]
      node <- if (length(backbone_tips) == 1) match(backbone_tips, backbone$tip.label) else ape::getMRCA(backbone, backbone_tips)
      cache[[key]] <- list(node = if (node_ids) as.integer(node) else NA_integer_, tips = paste(node_to_tips(backbone, node), collapse = ", "))
    }
    cache[[key]]
  }
  describe_taxa <- function (nodes_tips) { label_excluded_subclades(list(phy$tip.label[nodes_tips]), max_tips = 5) }

  ## One record per taxon (or per missing clade): the node of the grafted lineage (the tip, or the root of the missing clade)
  if (!is_clades)
  {
    tip_of <- match(taxa, phy$tip.label)
    records <- lapply(tip_of, function (tip)
    {
      if (is.na(tip)) { return(NULL) }
      # Maximal grafted subtree containing the taxon (without backbone tips): attached to the backbone at its parent
      root <- tip
      while (parent[root] != 0 && n_backbone[parent[root]] == 0) { root <- parent[root] }
      graft_node <- parent[root]
      lineage <- other_child(graft_node, root)
      sister <- other_child(parent[tip], tip)
      clade <- parent_clade(lineage)
      list(tip = as.integer(tip), root = as.integer(root), graft_node = as.integer(graft_node), parent_node = clade$node, parent_tips = clade$tips,
           graft_age = age[graft_node], divergence_age = age[parent[tip]],
           grafted_on = if (n_backbone[sister] == 0) describe_taxa(tips_below(sister)) else NA_character_)
    })
    column <- function (name, type) vapply(records, function (r) if (is.null(r)) na_of(type) else r[[name]], vector(type, 1))
    output <- data.frame(Taxon_name = taxa,
                         parent_node = column("parent_node", "integer"), parent_tips = column("parent_tips", "character"),
                         graft_age = column("graft_age", "numeric"), divergence_age = column("divergence_age", "numeric"),
                         grafted_on = column("grafted_on", "character"), stringsAsFactors = FALSE)
    lineage_node <- column("tip", "integer") # Node whose state is the state inherited by the grafted lineage
  } else {
    n_groups <- length(groups$group_clade_tips)
    group_taxa <- lapply(seq_len(n_groups), function (g) taxa[groups$group_ID == g])
    records <- lapply(seq_len(n_groups), function (g)
    {
      members <- match(group_taxa[[g]], phy$tip.label)
      if (anyNA(members)) { return(NULL) }
      root <- if (length(members) == 1) members else ape::getMRCA(phy, phy$tip.label[members])
      if (n_below[root] != length(members) || parent[root] == 0)
      {
        return(list(monophyletic = FALSE))
      }
      graft_node <- parent[root]
      clade <- parent_clade(other_child(graft_node, root))
      list(monophyletic = TRUE, root = as.integer(root), graft_node = as.integer(graft_node), parent_node = clade$node, parent_tips = clade$tips,
           graft_age = age[graft_node], crown_age = if (length(members) == 1) NA_real_ else age[root])
    })
    broken <- which(vapply(records, function (r) !is.null(r) && !isTRUE(r$monophyletic), logical(1)))
    if (length(broken) > 0)
    {
      warning(length(broken), " missing clade(s) are not monophyletic in a phylogeny (the taxa of the clade do not form a clade): their values are NA. Missing clade ID(s): ",
              paste(utils::head(broken, 10), collapse = ", "), if (length(broken) > 10) ", ..." else "", call. = FALSE)
    }
    valid <- function (r) !is.null(r) && isTRUE(r$monophyletic)
    column <- function (name, type) vapply(records, function (r) if (valid(r)) r[[name]] else na_of(type), vector(type, 1))
    output <- data.frame(Missing_clade = seq_len(n_groups),
                         Taxa = vapply(group_taxa, paste, character(1), collapse = ", "), n_taxa = lengths(group_taxa),
                         parent_node = column("parent_node", "integer"), parent_tips = column("parent_tips", "character"),
                         graft_age = column("graft_age", "numeric"), crown_age = column("crown_age", "numeric"), stringsAsFactors = FALSE)
    lineage_node <- column("root", "integer")
    graft_nodes <- column("graft_node", "integer")
  }
  if (!is_clades) { graft_nodes <- column("graft_node", "integer") }

  ## States at the graft point and of the grafted lineage
  if (!is.null(simmaps))
  {
    lineage_name <- if (is_clades) "clade_state" else "tip_state"
    for (trait in seq_along(simmaps))
    {
      suffix <- if (length(simmaps) > 1) paste0("_", if (is.null(names(simmaps)) || !nzchar(names(simmaps)[trait])) trait else names(simmaps)[trait]) else ""
      output[[paste0("parent_state", suffix)]] <- vapply(graft_nodes, function (n) state_of_edge(simmaps[[trait]], n), character(1))
      output[[paste0(lineage_name, suffix)]] <- vapply(lineage_node, function (n) state_of_edge(simmaps[[trait]], n), character(1))
    }
  }
  output
}


#' @title Missing value of a type
#'
#' @param type Character string: `"integer"`, `"numeric"` or `"character"`.
#'
#' @noRd
#'
#' @return The missing value of the type.
#'

na_of <- function (type) { switch(type, integer = NA_integer_, numeric = NA_real_, character = NA_character_) }
