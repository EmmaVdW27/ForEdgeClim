#' Function to calculate ground heat flux G
#'
#' @param net_rad_ground Net radiation at ground surface
#' @return G Ground heat flux
#' @export
#'
# G is positive if E is entering the soil
# More dense voxels lose less heat to ground
calculate_G <- function(net_rad_ground){
  G = p_ground*(1-den)*net_rad_ground
  return(G)
}

#' Function to simulate air to air diffusion via Fourier's law of heat conduction. This is done in 3D.
#' In this function, air temperature is updated and ready to use for the sensible heat flux calculation.
#'
#' @param temp_air Air temperature
#' @param temp_soil Ground surface temperature
#' @return total_diffusion Diffusion flux per voxel
#' @export
#'
calculate_D <- function(temp_air, temp_soil){
  area = voxel_length^2  # Voxel 1 side area (m2)
  V <<- voxel_length^3     # Voxel volume (m3)

  nx <- x_dim
  ny <- y_dim
  nz <- z_dim
  dx = voxel_length
  dy = voxel_length
  dz = voxel_length

  # Convert vector to 3D-array
  T_air_grid = array(temp_air, dim = c(nx, ny, nz))
  T_soil_grid = array(temp_soil, dim = c(nx, ny, nz))

  # Initialize air heat diffusion flux 3D array
  D_x <- array(0, dim = c(nx, ny, nz))
  D_y <- array(0, dim = c(nx, ny, nz))
  D_z <- array(0, dim = c(nx, ny, nz))

  # Simulate air diffusion via Fourier's law
  # Zero-flux boundary condition at 'core boundaries', ie, min x, min y & max y.
  # At these boundaries the forest is assumed to continue indefinitely.

  # Boundary conditions for x-direction (including eastern and western forest edge)
  D_x[nx, , ] <-  h * area * (T_air_grid[nx, , ] - macro_temp) # Eastern edge (max x)
  D_x[1, , ] <- 0 # Western edge (min x)

  D_x[2:nx, , ] <- D_x[2:nx, , ] + h * area * (T_air_grid[2:nx, , ] - T_air_grid[1:(nx-1), , ])
  D_x[1:(nx-1), , ] <- D_x[1:(nx-1), , ] + h * area * (T_air_grid[1:(nx-1), , ] - T_air_grid[2:nx, , ])

  # Boundary conditions for y-direction (including northern and southern forest edge)
  D_y[, ny, ] <- 0 # Northern edge (max y)
  D_y[, 1, ] <- 0 # Southern edge (min y)

  D_y[, 1:(ny-1), ] <- D_y[, 1:(ny-1), ] + h * area * (T_air_grid[, 1:(ny-1), ] - T_air_grid[, 2:ny, ])
  D_y[, 2:ny, ] <- D_y[, 2:ny, ] + h * area * (T_air_grid[, 2:ny, ] - T_air_grid[, 1:(ny-1), ])

  # Boundary conditions for z-direction (including upper canopy and ground)
  D_z[, , nz] <- h * area * (T_air_grid[, , nz] - macro_temp)  # Upper canopy (max z)
  D_z[, , 1]  <- h * area * (T_air_grid[, , 1] - T_soil_grid[ , , 1])  # Ground (min z)

  D_z[, , 2:nz] <- D_z[, , 2:nz] + h * area * (T_air_grid[, , 2:nz] - T_air_grid[, , 1:(nz-1)])
  D_z[, , 1:(nz-1)] <- D_z[, , 1:(nz-1)] + h * area * (T_air_grid[, , 1:(nz-1)] - T_air_grid[, , 2:nz])

  # Total air diffusion flux per voxel
  # This value is positive if heat energy leaves the voxel
  total_diffusion <- D_x + D_y + D_z

  return(as.vector(total_diffusion))
}

#' Function to calculate sensible heat transfer between forest surface and air
#'
#' @param T_surf Forest surface temperature
#' @param T_air Air temperature
#' @return H Sensible heat flux
#' @export
#'
calculate_H <- function(T_surf,T_air){
  # H is simulated via 'convective conductance' where g functions as an effective conductance that accounts for
  # both conduction within the boundary layer and convection around it.
  # H is positive when energy is lost from the surface to the air.
  H = den * g_forest * (T_surf - T_air)

  return(H)

}

#' Function for saturated vapor pressure; this is the empirical equation of Tetens
#'
#' @param temp Forest surface temperature
#' @return Saturated vapor pressure
#' @export
#'
saturated_vapor_pressure <- function(temp) {
  0.6108 * exp(17.27 * temp / (temp + 237.3)) # temp in °C, result in kPa
}

#' Function to calculate latent heat flux LE (~ evapotranspiration, ET: LE = L_v . ET with L_v = latent heat of vaporization of water (J/kg))
#'
#' @param temperature Forest surface temperature
#' @param net_rad Net radiation
#' @return LE Latent heat flux
#' @export
#'
calculate_LE <- function(temperature, net_rad) {
  # slope of the saturation pressure curve; temp in °C; slope in kPa/K = kPa/°C
  slope = 4098 * saturated_vapor_pressure(temperature) / (temperature + 237.3)^2
  # Latent heat flux by the empirical method of Priestly-Taylor
  LE = den * alpha_PT * (net_rad) * slope / (slope + gamma_psy)
  # LE cannot be negative
  LE[LE<0] = 0

  return(LE) # LE in W/m²
}
