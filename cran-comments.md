## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* Note 1:

> Suggests or Enhances not in mainstream repositories:
    contsimmap
  Availability using Additional_repositories specification:
    contsimmap   yes   https://maeldore.github.io/drat
    
  One feature of phylograftR relies on the contsimmap package, which is not currently hosted on CRAN 
  but is an actively maintained R package used in macroevolutionary research.
  This package is indicated as Suggests and conditions are implemented to check for its presence when running examples that rely on contsimmap,
  ensuring the package is functional even if contsimmap is not installed.
