# zalpha

R package for the zalpha suite of selection statistics. This package allows the user to run a suite of selection statistics over genomic data. The aim is to identify regions of the genome which have experienced a selective sweep.

## Installing

To install this package:

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
install_github("chorscroft/zalpha", build_vignettes=TRUE)  # might also have to include the option force=TRUE if the package is already installed
vignette("zalpha")
```

## Authors

* **Clare Horscroft**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

Guy Jacobs for the original development of the Zalpha statistics
