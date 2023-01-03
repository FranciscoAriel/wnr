# wnr

This repository contains source code for the R-package `wnr`.

It estimates parameters of a weibull distribution using a Newton-Raphson Algorithm.

This package is experimental and was done with educational purposes.

R-tools is required in order to compile the package in windows. You can download them from [here](https://cran.r-project.org/bin/windows/Rtools/).

You can install this package from github using:

````r
library(devtools)
install_github(repo = "FranciscoAriel/wnr")
````

This package was written and built in Windows 11, some tests were done in Ubuntu using [WSL](https://learn.microsoft.com/en-us/windows/wsl/about).

## History versions

* 0.3 Stable version ([issue #1](https://github.com/FranciscoAriel/wnr/issues/1#issue-1516963574) solved)
* 0.2 Some bugs solved
* 0.1 First release