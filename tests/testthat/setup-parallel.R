## Runs once before all the tests.

# Rationale: On Windows, every background R process (the workers of the 'future' plans used in the parallel tests) creates a file 'Rscript<hex>' in its temporary directory
# at startup and does not delete it (a known issue of R for Windows). Workers inherit the temporary directory of the check, so these files are reported by R CMD check
# ('checking for detritus in the temp directory'). The workers are given their own temporary directory instead, inside the temporary directory of the test session,
# which is deleted at the end of the session. Only the processes started by the tests are concerned (the temporary directory of the session itself is already set).

.tmp_vars <- c("TMPDIR", "TMP", "TEMP")
options(phylograftR.test_old_tmp = Sys.getenv(.tmp_vars, unset = NA))
.workers_tmp <- tempfile("workers_tmp_")
dir.create(.workers_tmp, showWarnings = FALSE)
do.call(Sys.setenv, stats::setNames(as.list(rep(.workers_tmp, length(.tmp_vars))), .tmp_vars))
rm(.tmp_vars, .workers_tmp)
