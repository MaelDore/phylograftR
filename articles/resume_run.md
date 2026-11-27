# Resume a crashed or interrupted run

Grafting **thousands of taxa** on a large backbone phylogeny, with
**many simulations and stochastic maps**, can last **hours or days**
(see
[`vignette("parallel_grafting")`](https://maeldore.github.io/phylograftR/articles/parallel_grafting.md)).
Since grafting functions work with a single call, a power cut, a closed
session, a memory crash or a cluster time limit would then **destroy all
the work done so far**.

To prevent such issue, the grafting functions can **save their progress
in temporary files**, and **resume** from them after an interruption.
The resumed run gives **exactly the same output** as a run that would
never have been interrupted.

Here, we cover:

1.  How to **save** the progress with `save_temp`, `save_temp_every` and
    `temp_file_path`.
2.  How the **temporary files** are named and what they contain.
3.  How to **resume** a run with `resume = TRUE`, and when the files are
    **deleted** (`clean_temp_files`).
4.  A **real example**: a run is killed, then resumed.
5.  Resume a run **in parallel**.
6.  The **pitfalls** to avoid.
7.  How to choose **`save_temp_every`**.

Everything below applies in the same way to the three grafting
functions:
[`graft_missing_taxa()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa.md),
[`graft_missing_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_clades.md)
and
[`graft_missing_taxa_and_clades()`](https://maeldore.github.io/phylograftR/reference/graft_missing_taxa_and_clades.md).

## 0/ Data

We use the toy backbone, the toy list of **7 missing taxa**, and the toy
**`multiSimmap`** of 10 stochastic maps of a color with three states
(see
[`vignette("handle_uncertainty")`](https://maeldore.github.io/phylograftR/articles/handle_uncertainty.md)).

``` r

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
```

The toy runs are very fast. To show how the saves work, we first use
them on small runs, and then on a **real interrupted run** (*Section
4*).

## 1/ Save the progress

Three arguments control the temporary saves:

- `save_temp = TRUE` activates the temporary save system. Default is
  `FALSE`, so you need to activate the system by setting
  `save_temp = TRUE`.
- `save_temp_every`: the **frequency of the saves of the progress**, in
  number of taxa grafted (`method = "uniform"`), or groups of taxa
  grafted together (`method = "BD_clades"`, taxa that share the same
  target clade, stem and allowed states). With `NULL` (default), a
  phylogeny is only saved **once it is fully grafted**.
- `temp_file_path`: the path of the `.rds` file used for the saves,
  which must end with `.rds`. Its folder must exist and be writable.
  Default is `"temp_phy_from_grafting.rds"` in the working directory.

``` r

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
#> character(0)
```

By default (`clean_temp_files = TRUE`), the files are **deleted when the
function returns without crashing**. **A run that crashes or is
interrupted never deletes its files**, so that it can be resumed. Use
`clean_temp_files = FALSE` to **keep the files after a successful run**.

## 2/ The temporary files

Each **grafted phylogeny** (one per simulation and per stochastic map)
has **its own file**. The path is a **template**: the suffixes
`"_map_<name>"` (for a `multiSimmap`) and `"_sim_<k>"` (for several
simulations) are added before `".rds"`, in this order. The prefix
`"map_"` is **not repeated** if the name of the map already starts with
it, as in the default names `map_1`, `map_2`, … of a `multiSimmap`
(e.g., `"temp_phy_map_1_sim_1.rds"`).

``` r

# 3 simulations
grafted_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                   nb_simulations = 3, seed = 1, verbose = FALSE,
                                   save_temp = TRUE, temp_file_path = temp_file,
                                   clean_temp_files = FALSE) # Keep the files
list.files(temp_dir)
#> [1] "temp_phy_sim_1.rds" "temp_phy_sim_2.rds" "temp_phy_sim_3.rds"

# 2 maps x 2 simulations
grafted_maps_sims <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                        simmap = toy_multiSimmap[1:2], nb_simulations = 2,
                                        seed = 1, verbose = FALSE,
                                        save_temp = TRUE, temp_file_path = temp_file,
                                        clean_temp_files = FALSE)
list.files(temp_dir)
#> [1] "temp_phy_map_1_sim_1.rds" "temp_phy_map_1_sim_2.rds"
#> [3] "temp_phy_map_2_sim_1.rds" "temp_phy_map_2_sim_2.rds"
#> [5] "temp_phy_sim_1.rds"       "temp_phy_sim_2.rds"      
#> [7] "temp_phy_sim_3.rds"
```

A temporary file is a **regular `"phylo"` object**, so it can be read
with [`readRDS()`](https://rdrr.io/r/base/readRDS.html) and plotted.
Once a phylogeny is fully grafted, its file holds the **complete output
of the simulation** (`$grafted_tips`, `$simmap`, etc.). For an
**unfinished phylogeny**, the file holds the phylogeny **grafted so
far**, and the attribute `"phylograftR_temp"` stores what is needed to
resume: the state of the grafting, the random numbers generator, and a
signature of the inputs used to ensure the run is resumed with identical
inputs (see Pitfalls in *Section 6*).

``` r

temp_phy <- readRDS(file.path(temp_dir, "temp_phy_sim_1.rds"))
class(temp_phy)
#> [1] "phylo"
attr(temp_phy, "phylograftR_temp")$status
#> [1] "complete"
```

With `force_simulations = TRUE`, `<k>` is the **number of the attempt**,
not necessarily the ID of the simulation in the output, as failed
attempts leave a file too.

## 3/ Resume a run

To resume, **call the function again with exactly the same arguments,
and add `resume = TRUE`** (it needs `save_temp = TRUE`). For each
grafted phylogeny, the function searches its temporary file, and:

- a **complete** phylogeny is **reused** as it is;
- an **unfinished** phylogeny is **continued from its last save**;
- a phylogeny **without file** is **grafted from scratch**.

Here, we mimic an interruption after the second simulation by **deleting
the file of the third simulation**, and resume the run.

``` r

# The run of 3 simulations of the previous section, whose files were kept
reference <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                                nb_simulations = 3, seed = 1, verbose = FALSE)

# Mimic an interruption before the third simulation by deleting the temp save for simulation 3
file.remove(file.path(temp_dir, "temp_phy_sim_3.rds"))
#> [1] TRUE

# Same arguments + resume = TRUE
resumed <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                              nb_simulations = 3, seed = 1, verbose = FALSE,
                              save_temp = TRUE, temp_file_path = temp_file,
                              resume = TRUE)

# Same output as a run never interrupted
identical(resumed, reference)
#> [1] TRUE

# Since we did not ask to retain the save files (i.e., clean_temp_files = TRUE by default)
# the files of this run were deleted at the end of the resumed run (the files of the other 2x2 run remain)
list.files(temp_dir)
#> [1] "temp_phy_map_1_sim_1.rds" "temp_phy_map_1_sim_2.rds"
#> [3] "temp_phy_map_2_sim_1.rds" "temp_phy_map_2_sim_2.rds"
```

The output is **identical** to the output of a run that was never
interrupted: the state of the **random numbers generator** is saved with
the progress, so the grafting continues with the **same draws** as it
would have without interruption. The resumed run keeps saving its
progress in the same files, so it can be **interrupted and resumed
several times**.

The files deleted at the end of a successful run are **only the files of
this run**. Files of other runs in the same folder, as the files of the
2 maps x 2 simulations above, are not touched.

## 4/ A real interrupted run

We now **kill a real run** and resume it. This example is **not run in
the tutorial** as it takes about a minute and relies on the timing of
the interruption, but you can run it as it is. The run grafts the 7 toy
taxa on the 10 maps, with 100 simulations per map = **1000 grafted
phylogenies** in total. We save the progress after each group of taxa
(i.e., sharing their target clade, stem availability, and allowed
states).

First, we write the script of the run in a file, `long_run.R`. The
folder `temp_files` must exist.

``` r

# Content of the script 'long_run.R'
library(phylograftR)

data(toy_backbone)
data(toy_multiSimmap)
data(toy_missing_taxa_list_with_allowed_states)
taxa_colors <- toy_missing_taxa_list_with_allowed_states
taxa_colors$allowed_states[["New_sp1"]] <- c("red", "green")

grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                    simmap = toy_multiSimmap, nb_simulations = 100,
                                    seed = 1, verbose = TRUE,
                                    save_temp = TRUE, save_temp_every = 1,
                                    temp_file_path = "temp_files/temp_grafting.rds")

saveRDS(grafted_phylo, "grafted_phylo.rds")
```

Then we start the script in a **background R session** with the package
`callr` (`install.packages("callr")`), and **kill it after 20 seconds**.
Any other way to stop the run works the same: the stop button of your R
console, a closed session, a power cut, or the time limit of a cluster
job.

``` r

dir.create("temp_files")

job <- callr::r_bg(function () source("long_run.R"), stdout = "run_log.txt", stderr = "run_errors.txt")
Sys.sleep(20)
job$kill() # Kill the run

# The temporary files left by the run
temp_files <- list.files("temp_files", full.names = TRUE)
length(temp_files)
#> [1] 275
head(basename(temp_files), 3)
#> [1] "temp_grafting_map_1_sim_1.rds"   "temp_grafting_map_1_sim_10.rds"
#> [3] "temp_grafting_map_1_sim_100.rds"
```

The run was stopped after **275 grafted phylogenies had a file**. All
are complete, except the one that was being grafted when the run was
killed: it was saved after 4 of the 7 groups of taxa.

``` r

# Status of each temporary file
status <- sapply(temp_files, function (f) {
  state <- attr(readRDS(f), "phylograftR_temp")
  if (is.null(state)) "complete" else state$status
})
table(status)
#> status
#> complete  partial
#>      274        1

# The unfinished phylogeny: map 3, simulation 75
partial <- readRDS("temp_files/temp_grafting_map_3_sim_75.rds")
attr(partial, "phylograftR_temp")$progress
#> $done
#> [1] 4
#>
#> $total
#> [1] 7
#>
#> $unit
#> [1] "groups"

# The phylogeny grafted so far can be plotted as any phylogeny
plot(partial)
```

To resume, we **run the same script again**, changing only
`resume = TRUE`. In a real situation, you can simply edit your script,
or re-run your own code, **keeping all the other arguments exactly the
same**.

``` r

grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                    simmap = toy_multiSimmap, nb_simulations = 100,
                                    seed = 1, verbose = TRUE,
                                    save_temp = TRUE, save_temp_every = 1,
                                    temp_file_path = "temp_files/temp_grafting.rds",
                                    resume = TRUE) # Only change
#> 2026-11-26 12:37:52 - Map map_1 - Simulation n°1 - Complete phylogeny found in the temporary file: reused
#> 2026-11-26 12:37:52 - Map map_1 - Simulation n°2 - Complete phylogeny found in the temporary file: reused
#> [...]
#> 2026-11-26 12:37:52 - Map map_3 - Simulation n°75 - Resuming from the temporary file: 4/7 groups of taxa already grafted
#> 2026-11-26 12:37:53 - Map map_3 - Simulation n°76 - No temporary file found: grafting from scratch
#> [...]
#> 2026-11-26 12:38:39 - 1000 temporary file(s) deleted
```

The 274 complete phylogenies were **reused**, the 75th simulation of the
map 3 **resumed from its 4th group of taxa**, and the others **grafted
from scratch**. The run was completed in 48 seconds, instead of 65
seconds for the full run, and the temporary files were deleted at the
end. With the real, long runs, the gain is larger: what was done before
the interruption is never repeated.

``` r

# The output is the same as the output of a run that was never interrupted
grafted_phylo_ref <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                        simmap = toy_multiSimmap, nb_simulations = 100,
                                        seed = 1, verbose = FALSE)
identical(grafted_phylo, grafted_phylo_ref)
#> [1] TRUE
```

## 5/ Resume a parallel run

The temporary saves work with the **parallel computation** (see
[`vignette("parallel_grafting")`](https://maeldore.github.io/phylograftR/articles/parallel_grafting.md))
with nothing more to do: each worker saves the files of the phylogenies
it grafts, in the same folder. If the main session is killed, **the
workers stop with it**. To resume, **set the plan again** and call the
function again with `resume = TRUE`.

``` r

future::plan(future::multisession, workers = 2)

grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                    simmap = toy_multiSimmap, nb_simulations = 100,
                                    seed = 1, save_temp = TRUE, save_temp_every = 1,
                                    temp_file_path = "temp_files/temp_grafting.rds")
# ---> Interrupted after 20 seconds: 421 temporary files, including 2 unfinished ones (1 per worker)

# Resume: same plan, same arguments, + resume = TRUE
grafted_phylo <- graft_missing_taxa(phy = toy_backbone, missing_taxa_list = taxa_colors,
                                    simmap = toy_multiSimmap, nb_simulations = 100,
                                    seed = 1, save_temp = TRUE, save_temp_every = 1,
                                    temp_file_path = "temp_files/temp_grafting.rds",
                                    resume = TRUE)

future::plan(future::sequential)
```

The resumed parallel run gave the **same output** as a sequential run
that was never interrupted. The number of workers can differ between the
interrupted and the resumed runs: each grafted phylogeny has its own
file and its own seed, whatever the worker that grafts it.

## 6/ Pitfalls

- **Use exactly the same arguments** to resume (`phy`,
  `missing_taxa_list`, `simmap`, `method`, `seed`, `nb_simulations`,
  etc.). A **signature of the inputs** is stored in each file, and the
  function **stops** if it does not match. The error message lists
  **what differs** (backbone, missing taxa, simmap, or other arguments),
  with the **previous value** (interrupted run) and the **new value**
  (current call):

``` r

graft_missing_taxa(phy = toy_backbone, missing_taxa_list = toy_missing_taxa_list,
                   nb_simulations = 3, seed = 2, # Not the seed of the interrupted run
                   save_temp = TRUE, temp_file_path = "temp_files/temp_phy.rds", resume = TRUE)
#> Error: Cannot resume from the temporary file '/path/to/temp_files/temp_phy_sim_1.rds': it was produced with different inputs (other arguments (method, seed, number of simulations or maps, ...)).
#> Differences between the interrupted run (previous) and the current call (new):
#>   - Seed of the simulation (derived from 'seed'): previous = 1, new = 2
#> Use the same inputs as in the interrupted run, or run with 'resume = FALSE' to overwrite the temporary files.
```

- **Do not forget `resume = TRUE`**: a call with `save_temp = TRUE` but
  without `resume = TRUE` starts again from scratch and **overwrites the
  temporary files**, and so the progress already done.
- **Do not use the temporary phylogeny as a new backbone**: grafting the
  taxa left on the phylogeny saved so far is **not equivalent** to
  resuming. The taxa already grafted would be treated as backbone
  species, so the sampling fractions of the birth-death models, the
  stems of the target clades, the identification of the grafted tips and
  the updated stochastic maps would be wrong. Use `resume = TRUE`
  instead.
- **Keep the same path and the same working directory**: a relative
  `temp_file_path` is resolved when the function is called. Use the same
  folder, or an absolute path, to find the files again.
- **Do not edit or move the files** of the run between the interruption
  and the resumption.
- **Temporary files are not a way to store results**: once a run is
  completed, they are deleted (by default), and the final output
  contains all relevant information for your analyses.

## 7/ How often should the progress be saved?

There are **two levels of protection**:

- **Without `save_temp_every`** (`NULL`), each grafted phylogeny is
  saved **once complete**. After an interruption, **only the phylogeny
  being grafted is lost** per worker. It is enough when you have **many
  phylogenies, each grafted quickly** (many simulations and maps on
  small or medium backbones).
- **With `save_temp_every`**, the progress **within** a grafted
  phylogeny is saved too. It matters when **a single grafted phylogeny
  is long**: a backbone with thousands of tips and hundreds or thousands
  of taxa to graft, where one phylogeny can take several minutes to more
  than ten minutes (see
  [`vignette("parallel_grafting")`](https://maeldore.github.io/phylograftR/articles/parallel_grafting.md)).

As a rule of thumb, set `save_temp_every` to **a few percent of the taxa
to graft** (e.g., `save_temp_every = 5` for 100 groups of taxa to
graft).

**Saving is cheap**: in a test on a backbone of 1000 tips with 300 taxa
to graft, saving after **every** taxon increased the time of the run by
about 20%, and saving every 10 taxa or less often had no measurable
cost. If `save_temp_every` is larger than the number of taxa (or of
groups of taxa) to graft, no progress is saved before the phylogeny is
complete (as if `save_temp_every = NULL`), and a warning is raised.

## Next steps

- Check how to **combine** missing taxa and missing clades in one run,
  with the temporary files of both steps: see
  [`vignette("combine_missing_taxa_and_clades")`](https://maeldore.github.io/phylograftR/articles/combine_missing_taxa_and_clades.md).
  Each step has its own file (suffixes `"_step_taxa"` and
  `"_step_clades"`), and a run interrupted during the second step
  resumes from it, reusing the first step.
- Run **many grafted phylogenies at the same time**: see
  [`vignette("parallel_grafting")`](https://maeldore.github.io/phylograftR/articles/parallel_grafting.md).
- **Check** the resulting grafted phylogenies against the constraints to
  ensure the (resumed) run was successful: see
  [`vignette("check_grafting")`](https://maeldore.github.io/phylograftR/articles/check_grafting.md).
