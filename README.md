# metaquant: Estimating Means, Standard Deviations, and Visualizing Distributions Using Quantiles

The `metaquant` package implements a novel density-based approach for estimating unknown means, visualizing distributions, and performing meta-analyses of quantiles, as introduced by De Livera et al.

## Overview

The `metaquant` package provides functions for estimating means, standard deviations, and visualizing distributions using quantile summary data. It is designed with a focus on meta-analyses of continuous outcomes, especially when studies report only quantile-based information (e.g., medians, quartiles, or extremes). By leveraging flexible quantile-based distributions, `metaquant` enables researchers to integrate quantile summary measures into a comprehensive meta-analysis framework.

## Installation

You can install  `metaquant` directly from CRAN as follows:
```R
# Install from CRAN
install.packages("metaquant")
library(metaquant)

Alternatively, the development version of the package is available on GitHub. To install this version, the user needs to ensure that Rtools has been installed and integrated beforehand.

```R
# Install from GitHub
# install.packages("devtools")
library(devtools)
devtools::install_github("metaanalysisR/metaquant")
