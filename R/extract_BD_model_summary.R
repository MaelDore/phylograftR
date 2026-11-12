##### Function to extract the birth-death models fitted during grafting #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################


### 1/ Extract the summary of the birth-death models ####

#' @title Extract the birth-death models fitted during grafting
#'
#' @description Extract the summary of the birth-death (BD) models fitted to simulate the divergence times of the grafted taxa.
#'   Returns a data.frame with one row per model, per grafted phylogeny (simulation and stochastic map).
#'   It is similar to [`phylograftR::extract_grafting_info()`][extract_grafting_info], but the summary is not computed a posteriori:
#'   the grafting functions record them if `return_BD_model_summary = TRUE`, and this function extracts the records.
#'
#' @param phy_grafted A grafted phylogeny (class `"phylo"`), a set of grafted phylogenies (class `"multiPhylo"`, e.g., several simulations),
#'   or a list of them with one element per map (output of the grafting functions with a `"multiSimmap"`).
#'   The phylogenies must have the attribute `"BD_model_summary"`, so they must be the output of
#'   [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] (with `method = "BD_clades"`) or [`phylograftR::graft_missing_clades()`][graft_missing_clades]
#'   run with `return_BD_model_summary = TRUE`.
#'
#' @details
#' With [`phylograftR::graft_missing_taxa()`][graft_missing_taxa], **a model is fitted per target clade**, and used by all the taxa of the target clade,
#' whatever their stem availability or excluded subclades.
#' With [`phylograftR::graft_missing_clades()`][graft_missing_clades], **a model is fitted per missing clade**, i.e., per sister clade.
#'
#' The initial clade is the target clade (or the sister clade). The model is fitted on the initial clade if it has at least `min_tips` tips,
#' otherwise on the **smallest enclosing clade** that has `min_tips` tips (`enclosing = TRUE`).
#' The number of tips of the fitted clade counts the taxa already grafted when the model is fitted.
#'
#' Clades are identified by their node in the backbone phylogeny provided to the grafting function, and by a pair of tips that defines the clade
#' as their most recent common ancestor (a single tip for a clade made of one tip).
#'
#' The death rate used to simulate the divergence times (`d`) can differ from the one fitted by `phytools::fit.bd()` (`d_fitted`):
#' it is set to `0` if it is not finite or negative (`d_adjustment = "set_to_0"`), and capped at 99% of the birth rate (`d_adjustment = "capped"`),
#' so the net diversification rate is always positive, which is needed to simulate divergence ages.
#'
#' @export
#'
#' @seealso Other functions from phylograftR:
#'  * [`phylograftR::graft_missing_taxa()`][graft_missing_taxa] to graft taxa within a target clade, with constraints (target clades and excluded subclades).
#'  * [`phylograftR::graft_missing_clades()`][graft_missing_clades] to graft groups of missing taxa as clades sister to clades of the backbone.
#'  * [`phylograftR::graft_missing_taxa_and_clades()`][graft_missing_taxa_and_clades] to graft missing taxa and missing clades in one run, with a single set of inputs.
#'  * [`phylograftR::extract_grafting_info()`][extract_grafting_info] to extract a posteriori information on the location of grafted lineages.
#'  * [`phylograftR::extract_sampling_info()`][extract_sampling_info] to extract the richness and the sampling fraction of the clades in each grafted phylogeny.
#'
#' @return A data.frame, with one row per model and phylogeny. The columns are:
#'   * `map` Character string. Name of the stochastic map (`NA` if no map was used).
#'   * `simulation` Integer. Index of the simulation.
#'   * `seed` Integer. Seed of the simulation (`NA` if unknown).
#'   * `Clade_name` Character string. Name of the initial clade (`NA` if it has no name). Only if `Target_clade_name` (or `Sister_clade_name`) was provided in `missing_taxa_list`.
#'   * `Taxa` Character string. Taxa that were grafted under this model, separated by commas.
#'   * `n_taxa` Integer. Number of taxa that were grafted under this model.
#'   * `initial_clade_node` Integer. Node ID in the backbone of the MRCA of the initial clade (target clade or sister clade). The ID of the tip for a single tip.
#'   * `initial_clade_tips` Character string. Pair of tips that defines the initial clade.
#'   * `initial_clade_n_tips` Integer. Number of tips of the initial clade in the backbone.
#'   * `initial_clade_crown_age` Numeric. Crown age of the initial clade (`NA` for a single tip).
#'   * `enclosing` Logical. Whether the model was fitted on an enclosing clade, because the initial clade has fewer than `min_tips` tips.
#'   * `fitted_clade_node` Integer. Node ID in the backbone of the MRCA of the clade used to fit the model (the same as `initial_clade_node` if `enclosing = FALSE`).
#'   * `fitted_clade_tips` Character string. Pair of tips that defines the fitted clade.
#'   * `fitted_clade_n_tips` Integer. Number of tips of the fitted clade when the model was fitted, including the taxa already grafted.
#'   * `fitted_clade_n_backbone_tips` Integer. Number of tips of the fitted clade in the backbone.
#'   * `fitted_clade_crown_age` Numeric. Crown age of the fitted clade.
#'   * `min_tips` Integer. Minimum number of tips required to fit the model.
#'   * `b` Numeric. Birth rate (speciation rate), in events per time unit per lineage.
#'   * `d` Numeric. Death rate (extinction rate) used to simulate the divergence times.
#'   * `d_fitted` Numeric. Death rate fitted by `phytools::fit.bd()`, before the adjustments.
#'   * `d_adjustment` Character string. `"none"`, `"set_to_0"` or `"capped"` (see Details).
#'   * `net_diversification` Numeric. `b - d`.
#'   * `relative_extinction` Numeric. `d / b`.
#'   * `rho` Numeric. Sampling fraction of the fitted clade: tips present divided by tips present, taxa still to graft in the clade and unsampled species.
#'   * `n_pending` Integer. Number of taxa still to graft in the fitted clade when the model was fitted.
#'   * `n_unsampled` Numeric. Number of unsampled species in the fitted clade (see `total_richness_per_clades` in the grafting functions).
#'
#'   Attributes: `"type"` (`"taxa"` for models of `graft_missing_taxa()`, `"clades"` for models of `graft_missing_clades()`, `"taxa_and_clades"` for models of `graft_missing_taxa_and_clades()`, which have a column `step`).
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
#' # List the missing taxa to graft: 2 taxa in a clade of 9 tips, 1 taxon in a clade of 3 tips
#' missing_taxa_list <- list(
#'   Taxon_name = c("X1", "X2", "X3"),
#'   Target_clade = list(X1 = c("t26", "t3"), X2 = c("t26", "t3"), X3 = c("t5", "t12")))
#'
#' # --- Graft the taxa, and report the BD models fitted --- #
#'
#' grafted_phys <- graft_missing_taxa(phy = phy, missing_taxa_list = missing_taxa_list,
#'                                    return_BD_model_summary = TRUE,
#'                                    nb_simulations = 3, verbose = FALSE, seed = 1234)
#'
#' # --- Extract the BD models --- #
#'
#' BD_summary <- extract_BD_model_summary(phy_grafted = grafted_phys)
#' BD_summary[, c("simulation", "Taxa", "initial_clade_n_tips",
#'                "enclosing", "fitted_clade_n_tips", "b", "d", "rho")]
#'
#' # The clade of 3 tips is smaller than min_tips = 5: the model is fitted on its enclosing clade
#' BD_summary[BD_summary$enclosing, c("Taxa", "initial_clade_tips", "fitted_clade_tips")]
#'

extract_BD_model_summary <- function (phy_grafted) # Output of the grafting functions run with return_BD_model_summary = TRUE: phylo, multiPhylo, or list of maps
{
  ### ---- Check input validity ---- ###

  listed <- list_grafted_phylogenies(phy_grafted, arg = "phy_grafted")
  tables <- lapply(listed$phylos, function (phy) attr(phy, "BD_model_summary"))
  missing_table <- which(vapply(tables, is.null, logical(1)))
  if (length(missing_table) > 0)
  {
    stop(length(missing_table), " of the ", length(tables), " phylogenies have no 'BD_model_summary' attribute.\n",
         "Run the grafting function with 'return_BD_model_summary = TRUE' (it requires method = 'BD_clades' with graft_missing_taxa()). ",
         "The models are fitted during the grafting: they cannot be recovered a posteriori.")
  }

  ### ---- Long table ---- ###

  n_rows <- vapply(tables, nrow, integer(1))
  output <- do.call(rbind, lapply(seq_along(tables), function (k)
  {
    cbind(map = rep(listed$map[k], n_rows[k]), simulation = rep(listed$simulation[k], n_rows[k]), seed = rep(listed$seed[k], n_rows[k]), tables[[k]], stringsAsFactors = FALSE)
  }))
  rownames(output) <- NULL
  attr(output, "type") <- attr(tables[[1]], "type")

  return(output)
}
