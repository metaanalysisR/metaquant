# metaquant: Estimating Means, Standard Deviations, and Visualizing Distributions Using Quantiles

The `metaquant` package implements a novel density-based approach for estimating unknown means, visualizing distributions, and meta-analyses of quantiles, as introduced by De Livera et al. The package offers functions for estimating means and standard deviations using both novel and existing methods, as well as visualizing sample distributions using quantile summary data. It is designed with a focus on meta-analyses of continuous outcomes, especially when studies report only quantile-based information (e.g., medians, quartiles, or extremes). Using flexible quantile-based distributions, `metaquant` enables researchers to integrate quantile summary measures into a comprehensive meta-analysis framework.

## Installation

You can install  `metaquant` directly from CRAN as follows:
```R
# Install from CRAN
install.packages("metaquant")
library(metaquant)
```
Alternatively, the development version of the package is available on GitHub. To install this version, the user needs to ensure that Rtools has been installed and integrated beforehand.

```R
# install.packages("devtools")
library(devtools)
devtools::install_github("metaanalysisR/metaquant")
```
## Usage

A detailed vignette with example datasets and code to prepare data and analyses are available at <https://bookdown.org/a2delivera/metaquant/>

## Reference

Alysha De Livera, Luke Prendergast, and Udara Kumaranathunga. A novel density-based approach for estimating unknown means, distribution visualisations, and meta-analyses of quantiles. *Submitted for Review*, 2024, pre-print available here: https://arxiv.org/abs/2411.10971 

