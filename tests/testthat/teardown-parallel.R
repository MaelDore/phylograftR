## Runs once after all the tests: make sure no parallel worker is left alive, and restore the temporary directory variables changed by setup-parallel.R

if (requireNamespace("future", quietly = TRUE))
{
  # Switching to a sequential plan stops the workers of the multisession plans used in the tests
  future::plan(future::sequential)
}
invisible(gc())

# Leave a moment to the worker processes to exit and release their temporary files
Sys.sleep(1)
old_tmp <- getOption("phylograftR.test_old_tmp")
if (!is.null(old_tmp))
{
  for (v in names(old_tmp))
  {
    if (is.na(old_tmp[[v]])) { Sys.unsetenv(v) } else { do.call(Sys.setenv, stats::setNames(list(old_tmp[[v]]), v)) }
  }
  options(phylograftR.test_old_tmp = NULL)
}
