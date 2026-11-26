## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----install, eval = FALSE----------------------------------------------------
# install.packages(c("future", "future.apply"))

## ----plan, eval = FALSE-------------------------------------------------------
# # How many workers? You can use up to the number of physical cores of your computer, minus one or two to keep it responsive
# parallel::detectCores(logical = FALSE)
# 
# # Set the plan: Here, 4 workers
# future::plan(future::multisession, workers = 4)
# future::nbrOfWorkers()
# #> [1] 4
# 
# # Run the grafting as usual
# # The argument 'parallel = NULL' (default) uses the plan if one is set
# grafted_phylo <- graft_missing_taxa(phy = toy_backbone,
#                                     missing_taxa_list = toy_missing_taxa_list,
#                                     nb_simulations = 100,
#                                     seed = 1)
# 
# # Reset the plan when you are done, to free the workers
# future::plan(future::sequential)

## ----parallel_true, eval = FALSE----------------------------------------------
# # No plan set: error
# graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
#                    nb_simulations = 4, parallel = TRUE)
# #> Error: 'parallel = TRUE' needs a 'future' plan with several workers, set before the call,
# #> e.g., future::plan(future::multisession, workers = 4). Use 'parallel = FALSE' to run sequentially.

## ----tasks, eval = FALSE------------------------------------------------------
# # 10 stochastic maps x 5 simulations = 50 tasks, shared among the workers
# future::plan(future::multisession, workers = 4)
# grafted_phylo_par <- graft_missing_taxa(phy = toy_backbone,
#                                         missing_taxa_list = taxa_colors, # See `vignette("handle_uncertainty")`
#                                         simmap = toy_multiSimmap,
#                                         nb_simulations = 5,
#                                         seed = 1)
# future::plan(future::sequential)

## ----identical, eval = FALSE--------------------------------------------------
# 
# # Run the same example without a parallel plan = sequentially
# grafted_phylo_seq <- graft_missing_taxa(phy = toy_backbone,
#                                         missing_taxa_list = taxa_colors,
#                                         simmap = toy_multiSimmap,
#                                         nb_simulations = 5,
#                                         seed = 1,
#                                         parallel = FALSE) # Ignore the plan
# 
# # Compare outputs
# identical(grafted_phylo_seq, grafted_phylo_par)
# #> [1] TRUE

## ----model_plot, echo = FALSE, fig.width = 7, fig.height = 5------------------
# Parameters estimated on the benchmark machine
backbones <- data.frame(tips = c(100, 500, 1000),
                        time_one = c(0.49, 10.9, 40.6),        # Seconds per grafted phylogeny (sequential)
                        overhead = c(0.17, 0.09, 0.01))        # Relative overhead of the parallel run
fixed_cost <- c("2" = 13.5, "4" = 18.2, "10" = 32.3)           # Seconds to set up the workers

speed_up <- function (N, time_one, workers, overhead)
{
  overhead <- if (workers == 2) 0 else overhead
  time_par <- fixed_cost[as.character(workers)] + (N / workers) * time_one * (1 + overhead)
  return(N * time_one / time_par)
}

N_grid <- exp(seq(log(4), log(1000), length.out = 150))
cols <- c("100" = "#1B9E77", "500" = "#D95F02", "1000" = "dodgerblue")
ltys <- c("2" = 3, "4" = 2, "10" = 1)

plot(NA, xlim = c(4, 1000), ylim = c(0, 10.5), log = "x",
     xlab = "Number of grafted phylogenies (P x N)", ylab = "Speed-up (sequential time / parallel time)")
abline(h = 1, lty = 1, lwd = 2, col = "red")
for (i in seq_len(nrow(backbones)))
{
  b <- backbones[i, ]
  for (W in c(2, 4, 10))
  {
    lines(N_grid, sapply(N_grid, speed_up, time_one = b$time_one, workers = W, overhead = b$overhead),
          col = cols[as.character(b$tips)], lty = ltys[as.character(W)], lwd = 2)
  }
}

legend("topleft", legend = c("100 tips", "500 tips", "1000 tips", "2 workers", "4 workers", "10 workers"),
       col = c(cols, rep("grey30", 3)), lty = c(1, 1, 1, 3, 2, 1), lwd = 2, bty = "n", cex = 0.8)

## ----time_one, eval = FALSE---------------------------------------------------
# # Time of one grafted phylogeny (sequential), including the preparation of the grafting
# time_one <- system.time(graft_missing_taxa(phy = my_backbone,
#                                            missing_taxa_list = my_missing_taxa_list,
#                                            simmap = my_simmap[[1]], # Optional: a single map
#                                            nb_simulations = 1, seed = 1,
#                                            parallel = FALSE))["elapsed"]
# 
# # Expected sequential time of the full run, in minutes
# time_one * P * N / 60

