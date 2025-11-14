#' Generate a voxelized point density grid from a LAS file
#'
#' This function reads a LAS file and creates a 3D voxel grid where each voxel represents
#' a 1x1x1m^3 volume (or a user-defined size). The density of each voxel is computed
#' based on the number of points in the point cloud. The output is a data.table containing
#' rescaled X, Y, and Z coordinates and a normalized density value.
#'
#' Additionally, a 3D plot of the voxel grid and a digital terrain model (DTM) are generated and saved as PNG files in
#' the "Output" directory, i.e. 'output_path' as defined by the user (see script 'testing_foredgeclim.R')
#'
#' @param las_file Path to the LAS file
#' @param voxel_size Numeric value defining the voxel size (default = 1m).
#' @return A data.table with columns: X, Y, Z, and density.
#' @importFrom dplyr group_by summarise
#' @import ggplot2
#' @importFrom lidR readLAS classify_ground rasterize_terrain filter_poi voxel_metrics csf knnidw normalize_height
#' @importFrom data.table as.data.table setnames CJ uniqueN
#' @export
generate_DTM_grid_TLS <- function(las_file, voxel_size = 1) {

  start = Sys.time()

  # Read the LAS file
  # Voxel density (voxel point cloud value) was determined via octree filtering in RiSCAN:
  # Laser beam divergence of RIEGL VZ-400i = 0.02°, at 38m distance (~height canopy), this translates
  # into a distance error of 0.0256m: tan(0.02°) = (x/2)/38m => x = 0.0256m
  # Because of this, a 5cm downsampling should do to correct for distance to scan positions.
  # Remark that this does not correct for occlusion (nor for the very dense ground layer).
  las <- readLAS(las_file)

  #######
  # DTM #
  #######

  # Use a filter to retain only the ground points (create DTM)
  dtm <- classify_ground(las, csf())
  # Create a raster model of the ground (DTM)
  dtm_raster <- rasterize_terrain(dtm, res = 1, algorithm = knnidw(k=10L))
  # Convert the raster to a dataframe
  dtm_df <- as.data.frame(dtm_raster, xy = TRUE)
  colnames(dtm_df) <- c("X", "Y", "Z")

  # Normalize height of the LAS file (everywhere set ground level to Z = 0)
  las <- normalize_height(las, dtm_raster)
  # Remove the points below ground level (for points that were e.g. wrongly not classified as ground points)
  las <- filter_poi(las, Z >= 0)

  ###################
  # STRUCTURAL GRID #
  ###################

  # Compute the 3D voxel grid with point count per voxel
  voxels <- voxel_metrics(las, ~length(Z), res = voxel_size)

  # Convert to a data.table and rename columns
  voxel_grid <- as.data.table(voxels)
  setnames(voxel_grid, c("X", "Y", "Z", "n"))

  # Remove voxels on min and max X & Y boundaries
  voxel_grid <- voxel_grid[voxel_grid$X > min(voxel_grid$X) & voxel_grid$X < max(voxel_grid$X), ]
  voxel_grid <- voxel_grid[voxel_grid$Y > min(voxel_grid$Y) & voxel_grid$Y < max(voxel_grid$Y), ]

  # Rescale coordinates so that the minimum value for each axis is 1
  voxel_grid$X <- voxel_grid$X - min(voxel_grid$X, na.rm = TRUE) + 1
  voxel_grid$Y <- voxel_grid$Y - min(voxel_grid$Y, na.rm = TRUE) + 1
  voxel_grid$Z <- voxel_grid$Z - min(voxel_grid$Z, na.rm = TRUE) + 1

  # Fill missing voxels with density = 0
  full_grid <- CJ(
    X = seq(min(voxel_grid$X), max(voxel_grid$X), by = voxel_size),
    Y = seq(min(voxel_grid$Y), max(voxel_grid$Y), by = voxel_size),
    Z = seq(min(voxel_grid$Z), max(voxel_grid$Z), by = voxel_size)
  )
  full_grid <- merge(full_grid, voxel_grid, by = c("X", "Y", "Z"), all.x = TRUE)
  full_grid[is.na(n), n := 0]  # Replace NA with 0 for voxels without points

  # Normalize the density values
  # Log transformation: Reduces dominance of high-density areas (mainly ground level) (necessary even after downsampling*) -> make distribution less extreme
  # *downsampling is done with the RiSCAN octree filter previous to the import of the las file
  # Exponential scaling: Boosts contribution of higher voxels even with lower density -> sort of correction for occlusion
  # The factor 0.5 determines how strong higher layers get more weight. 0.5 seems reasonable for a temperate forest.
  # For this application (a structural proxy in a microclimate model), this transformation is likely sufficient,
  # especially since the goal is not to quantify each voxel precisely, but rather to capture a spatial pattern.
  scaling <- log1p(full_grid$n) * exp(0.5*full_grid$Z/max(full_grid$Z))
  full_grid$density <- scaling/max(scaling)

  finish = Sys.time()
  print(paste0('Time to create DTM and structural grid from TLS data = ', round(as.numeric(finish - start, units = "secs"), 2), ' s'))

  # Return the voxel grid
  return(list(dtm = dtm_df,
              grid = as.data.frame(full_grid)))
}
