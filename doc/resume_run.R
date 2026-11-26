## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----setup, include = FALSE---------------------------------------------------
library(phylograftR)

## ----data---------------------------------------------------------------------
# ------ Toy data ------ #

data(toy_backbone)
data(toy_missing_taxa_list) # Taxonomic constraints only
data(toy_missing_taxa_list_with_allowed_states) # Same, with colors allowed for some taxa
data(toy_multiSimmap) # 10 stochastic maps of the color

# Allow "green" for New_sp1 in addition to "red", so that a branch is available in each map
taxa_colors <- toy_missing_taxa_list_with_allowed_states
taxa_colors$allowed_states[["New_sp1"]] <- c("red", "green")

# Folder for the temporary files of this tutorial
temp_dir <- file.path(tempdir(), "temp_grafting")
dir.create(temp_dir)

## ----save_temp----------------------------------------------------------------
# Path of the temporary file
temp_file <- file.path(temp_dir, "temp_phy.rds")

# One grafted phylogeny, with a save after each group of taxa grafted
grafted_phylo <- graft_missing_taxa(phy = toy_backbone,
                                    missing_taxa_list = toy_missing_taxa_list,
                                    seed = 1, verbose = FALSE,
                                    save_temp = TRUE, save_temp_every = 1,
                                    temp_file_path = temp_file)

# The temporary file is deleted once the run is completed (default clean_temp_files = TRUE)
list.files(temp_dir)

## ----temp_names---------------------------------------------------------------
# 3 simulations
grafted_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                   nb_simulations = 3, seed = 1, verbose = FALSE,
                                   save_temp = TRUE, temp_file_path = temp_file,
                                   clean_temp_files = FALSE) # Keep the files
list.files(temp_dir)

# 2 maps x 2 simulations
grafted_maps_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                        simmap = toy_multiSimmap[1:2], nb_simulations = 2,
                                        seed = 1, verbose = FALSE,
                                        save_temp = TRUE, temp_file_path = temp_file,
                                        clean_temp_files = FALSE)
list.files(temp_dir)

## ----read_temp----------------------------------------------------------------
temp_phy <- readRDS(file.path(temp_dir, "temp_phy_sim_1.rds"))
class(temp_phy)
attr(temp_phy, "phylograftR_temp")$status

## ----resume-------------------------------------------------------------------
# The run of 3 simulations of the previous section, whose files were kept
reference <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                nb_simulations = 3, seed = 1, verbose = FALSE)

# Mimic an interruption before the third simulation by deleting the temp save for simulation 3
file.remove(file.path(temp_dir, "temp_phy_sim_3.rds"))

# Same arguments + resume = TRUE
resumed <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                              nb_simulations = 3, seed = 1, verbose = FALSE,
                              save_temp = TRUE, temp_file_path = temp_file,
                              resume = TRUE)

# Same output as a run never interrupted
identical(resumed, reference)

# Since we did not ask to retain the save files (i.e., clean_temp_files = TRUE by default)
# the files of this run were deleted at the end of the resumed run (the files of the other 2x2 run remain)
list.files(temp_dir)

## ----clean_dir, include = FALSE-----------------------------------------------
unlink(temp_dir, recursive = TRUE)

## ----long_run_script, eval = FALSE--------------------------------------------
# # Content of the script 'long_run.R'
# library(phylograftR)
# 
# data(toy_backbone)
# data(toy_multiSimmap)
# data(toy_missing_taxa_list_with_allowed_states)
# taxa_colors <- toy_missing_taxa_list_with_allowed_states
# taxa_colors$allowed_states[["New_sp1"]] <- c("red", "green")
# 
# grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
#                                     simmap = toy_multiSimmap, nb_simulations = 100,
#                                     seed = 1, verbose = TRUE,
#                                     save_temp = TRUE, save_temp_every = 1,
#                                     temp_file_path = "temp_files/temp_grafting.rds")
# 
# saveRDS(grafted_phylo, "grafted_phylo.rds")

## ----kill, eval = FALSE-------------------------------------------------------
# dir.create("temp_files")
# 
# job <- callr::r_bg(function () source("long_run.R"), stdout = "run_log.txt", stderr = "run_errors.txt")
# Sys.sleep(20)
# job$kill() # Kill the run
# 
# # The temporary files left by the run
# temp_files <- list.files("temp_files", full.names = TRUE)
# length(temp_files)
# #> [1] 275
# head(basename(temp_files), 3)
# #> [1] "temp_grafting_map_1_sim_1.rds"   "temp_grafting_map_1_sim_10.rds"
# #> [3] "temp_grafting_map_1_sim_100.rds"

## ----inspect, eval = FALSE----------------------------------------------------
# # Status of each temporary file
# status <- sapply(temp_files, function (f) {
#   state <- attr(readRDS(f), "phylograftR_temp")
#   if (is.null(state)) "complete" else state$status
# })
# table(status)
# #> status
# #> complete  partial
# #>      274        1
# 
# # The unfinished phylogeny: map 3, simulation 75
# partial <- readRDS("temp_files/temp_grafting_map_3_sim_75.rds")
# attr(partial, "phylograftR_temp")$progress
# #> $done
# #> [1] 4
# #>
# #> $total
# #> [1] 7
# #>
# #> $unit
# #> [1] "groups"
# 
# # The phylogeny grafted so far can be plotted as any phylogeny
# plot(partial)

## ----resume_real, eval = FALSE------------------------------------------------
# grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
#                                     simmap = toy_multiSimmap, nb_simulations = 100,
#                                     seed = 1, verbose = TRUE,
#                                     save_temp = TRUE, save_temp_every = 1,
#                                     temp_file_path = "temp_files/temp_grafting.rds",
#                                     resume = TRUE) # Only change
# #> 2026-11-26 12:37:52 - Map map_1 - Simulation n°1 - Complete phylogeny found in the temporary file: reused
# #> 2026-11-26 12:37:52 - Map map_1 - Simulation n°2 - Complete phylogeny found in the temporary file: reused
# #> [...]
# #> 2026-11-26 12:37:52 - Map map_3 - Simulation n°75 - Resuming from the temporary file: 4/7 groups of taxa already grafted
# #> 2026-11-26 12:37:53 - Map map_3 - Simulation n°76 - No temporary file found: grafting from scratch
# #> [...]
# #> 2026-11-26 12:38:39 - 1000 temporary file(s) deleted

## ----check_identical, eval = FALSE--------------------------------------------
# # The output is the same as the output of a run that was never interrupted
# grafted_phylo_ref <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
#                                         simmap = toy_multiSimmap, nb_simulations = 100,
#                                         seed = 1, verbose = FALSE)
# identical(grafted_phylo, grafted_phylo_ref)
# #> [1] TRUE

## ----parallel, eval = FALSE---------------------------------------------------
# future::plan(future::multisession, workers = 2)
# 
# grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
#                                     simmap = toy_multiSimmap, nb_simulations = 100,
#                                     seed = 1, save_temp = TRUE, save_temp_every = 1,
#                                     temp_file_path = "temp_files/temp_grafting.rds")
# # ---> Interrupted after 20 seconds: 421 temporary files, including 2 unfinished ones (1 per worker)
# 
# # Resume: same plan, same arguments, + resume = TRUE
# grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
#                                     simmap = toy_multiSimmap, nb_simulations = 100,
#                                     seed = 1, save_temp = TRUE, save_temp_every = 1,
#                                     temp_file_path = "temp_files/temp_grafting.rds",
#                                     resume = TRUE)
# 
# future::plan(future::sequential)

## ----signature, eval = FALSE--------------------------------------------------
# graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
#                    nb_simulations = 3, seed = 2, # Not the seed of the interrupted run
#                    save_temp = TRUE, temp_file_path = "temp_files/temp_phy.rds", resume = TRUE)
# #> Error: Cannot resume from the temporary file '/path/to/temp_files/temp_phy_sim_1.rds': it was produced with different inputs (other arguments (method, seed, number of simulations or maps, ...)).
# #> Differences between the interrupted run (previous) and the current call (new):
# #>   - Seed of the simulation (derived from 'seed'): previous = 1, new = 2
# #> Use the same inputs as in the interrupted run, or run with 'resume = FALSE' to overwrite the temporary files.

