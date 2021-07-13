# zalpha
  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/chorscroft/zalpha/branch/master/graph/badge.svg)](https://codecov.io/gh/chorscroft/zalpha?branch=master)
  [![Travis build status](https://travis-ci.org/chorscroft/zalpha.svg?branch=master)](https://travis-ci.org/chorscroft/zalpha)
  [![](https://www.r-pkg.org/badges/version/zalpha?color=green)](https://cran.r-project.org/package=zalpha)
  [![CRAN checks](https://cranchecks.info/badges/summary/zalpha)](https://cran.r-project.org/web/checks/check_results_zalpha.html)
  [![DOI](https://joss.theoj.org/papers/10.21105/joss.02638/status.svg)](https://doi.org/10.21105/joss.02638)
  <!-- badges: end -->
R package for the zalpha suite of selection statistics. This package allows the user to run a suite of selection statistics over genomic data. The aim is to identify regions of the genome which have experienced a selective sweep.

## Installing

Version 0.2.0 is on CRAN.

**To install this package from CRAN:**

```
install.packages("zalpha")
```

**To install the development version of this package from github:**

1) Install the devtools package

```
install.packages("devtools")
```

2) Load the devtools package

```
library(devtools)
```

3) Install the zalpha package direct from Github using the install_github function

```
install_github("chorscroft/zalpha")
```

## How to use the zalpha package

For a full explanation of how to use the zalpha package, please see the "zalpha" vignette.

```
browseVignettes("zalpha")
# or
vignette("zalpha")
```
If this does not work and you installed the packages directly from Github, try installing using the following code:

```
install_github("chorscroft/zalpha", build_vignettes=TRUE)  # might also have to include the option force=TRUE if the package is already installed
```

## Dependencies

This package requires R version 2.10 or higher.

Optional packages:

* fitdistrplus

## Authors

* **Clare Horscroft**

## Acknowledgments

Guy Jacobs for the original development of the Zalpha statistics

## Community Guidelines

I am very open to feedback, bug reports and requests!

* If you have a bug to report or a new feature request, please create a new [issue](https://github.com/chorscroft/zalpha/issues) on GitHub.

* If you would like to contribute, please feel free to fork the project and create a [pull request](https://github.com/chorscroft/zalpha/pulls) on GitHub. A guide to how this works can be found [here](https://guides.github.com/activities/forking/).

* If you have any other feedback or need support, please email the author.
