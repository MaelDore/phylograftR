
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phylograftR

<!-- badges: start -->

<!-- 
usethis::use_cran_badge() reports the current version of your package on CRAN.
usethis::use_coverage() reports test coverage.
use_github_actions()  reports the R CMD check status of your development package. 
&#10;# Check badges in mFD package.
&#10;# Examples to adjust
&#10;[![CRAN status](https://www.r-pkg.org/badges/version/deepSTRAPP)](https://cran.r-project.org/package=deepSTRAPP) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads - Total](http://cranlogs.r-pkg.org/badges/grand-total/deepSTRAPP?color=blue)](https://cran.r-project.org/package=deepSTRAPP) # Total (looks better)
[![Downloads - Per month](https://cranlogs.r-pkg.org/badges/deepSTRAPP)](https://r-pkg.org:443/pkg/deepSTRAPP) # Per month
[![DOI](https://zenodo.org/badge/TBA.svg)](https://zenodo.org/doi/TBA)
&#10;-->

<!-- badges: end -->

## :dart: Summary

The goal of **phylograftR** is to graft missing taxa on time-calibrated
phylogenies combining taxonomic assigments and species attributes to
define inclusive/exclusive clades for grafting with diversification
modeling for inferring divergence times of imputed branches.

## :envelope_with_arrow: Installation

<br> phylograftR works on R version 4.6 or more. Be sure to have an R
version that is compatible. <br> See <https://CRAN.R-project.org/>.

From **CRAN**, for the latest release:

``` r
install.packages("phylograftR")
```

From **GitHub**, for the current development version:

``` r
library(devtools)

# If you want to have access to the vignettes/tutorials locally, run:
remotes::install_github(repo = "MaelDore/phylograftR")

# Altough, this is time-consuming, you can also opt for this light installation:
remotes::install_github(repo = "MaelDore/phylograftR", build_vignettes = FALSE)
# You will not have access to the vignettes/tutorials within R, but can still acess them through this website.
# See the dedicated Sections below.
```

## :desktop_computer: Website

<br> A companion website is available to browse interactively the
different tutorials and functions of **phylograftR** at this URL:
<https://maeldore.github.io/phylograftR/>.

An overview of all functions and datasets is available
[here](https://maeldore.github.io/phylograftR/reference/index.html).
<br><br>

## :joystick: Quick-to-run example

<br> A **simple use-case** that shows how phylograftR can be used to
**graft missing species based on taxonomic assignments** is available
[here](https://maeldore.github.io/phylograftR/articles/main_tutorial.html)
and within R: `vignette("main_tutorial")`.

## :bug: Found a bug?

<br> Thank you for finding it! Head over to the [GitHub Issues
tab](https://github.com/MaelDore/phylograftR/issues) and let me know
about it. <br> You can also [send me an
e-mail](mailto:mael.dore@gmail.com). <br><br>

## :black_nib: How to cite phylograftR

<br> For any use of **phylograftR**:

> Doré, M. phylograftR: Grafting missing taxa on time-calibrated
> phylogenies combining diversification modeling, taxonomic assignments,
> and species attributes. (DOI TBA)

**phylograftR** extends the **TACT framework** to enable the
identification of exclusive subclades where NOT to graft missing species
based on taxonomic assignments and/or species attributes. Please cite
the TACT paper as a foundation work for the development of
**phylograftR**:

> Chang, J., Rabosky, D. L., & Alfaro, M. E. (2019). Estimating
> diversification rates on incompletely-sampled phylogenies: theoretical
> concerns and practical solutions. *Systematic Biology*.
> <https://doi.org/10.1093/sysbio/syz081>

**TACT** is available as Python-based software at:
<https://tact.jonathanchang.org/>

As **phylograftR** relies strongly on functions designed for the
[phytools R package](https://blog.phytools.org/) and the [TreeSim R
package](https://CRAN.R-project.org/package=TreeSim), it is good
practice to also cite these packages:

> Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for
> phylogenetic comparative methods (and other things). *PeerJ*, 12,
> e16505. <https://doi.org/10.7717/peerj.16505>.

> Stadler T (2019). TreeSim: Simulating Phylogenetic Trees. R package
> version 2.4, <https://CRAN.R-project.org/package=TreeSim>
