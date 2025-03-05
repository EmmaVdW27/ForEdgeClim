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

#########
# INPUT #
#########

create_input_data(summer_day = T)

#######################
# Creation of 3D grid #
#######################

voxel_TLS = generate_DTM_grid_TLS(las_file = 'Data/June_2023_Gontrode_forest.las', voxel_size = 1)
saveRDS(voxel_TLS, 'Output/DTM_grid_TLS_DATA.rds')

#############
# Run model #
#############

res = run_foredgeclim(voxel_TLS$grid)
saveRDS(res, 'Output/model_results.rds')

############
# Plotting #
############

# Digital Terrain Model & structural grid plot
plots_dtm_struct(dtm = voxel_TLS$dtm, grid = voxel_TLS$grid)

# Shortwave radiation plots
plots_sw(sw_rad_2D = res$sw_rad_2D)

# Longwave radiation plots
plots_lw(lw_rad_2D = res$lw_rad_2D)

# Temperature & flux plots
plots_temp_flux(res$micro_grid, res$air_temperature, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux)

```
