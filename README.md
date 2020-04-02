# zalpha

R package for the zalpha suite of selection statistics. This package allows the user to run a suite of selection statistics over genomic data. The aim is to identify regions of the genome which have experienced a selective sweep.

## Installing

Version 0.1.0 is on CRAN.

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
vignette("zalpha")
```
If this does not work and you installed the packages directly from Github, try installing using the following code:

```
install_github("chorscroft/zalpha", build_vignettes=TRUE)  # might also have to include the option force=TRUE if the package is already installed
```
## Authors

* **Clare Horscroft**

## Acknowledgments

Guy Jacobs for the original development of the Zalpha statistics
