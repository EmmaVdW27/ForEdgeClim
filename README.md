# ForEdgeClim

ForEdgeClim is an R package for modelling microclimates in forests, including radiation processes and heat transfer.
The model is initially written to simulate microclimate gradients along transect lines from a forest’s core towards its edge.
A more detailed overview with the applied physical equations can be found in the attached PDF file 'Formulae ForEdgeClim.pdf'.

## Installation
Install the package directly from GitHub:
```r
# Using devtools
install.packages("devtools")
devtools::install_github("EmmaVdW27/ForEdgeClim")
```

## Example functions
```r
library(ForEdgeClim)
run_foredgeclim()
```
