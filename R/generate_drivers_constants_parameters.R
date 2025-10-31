#' Generate input driver variables
#'
#' @importFrom readxl read_excel
#' @return input driver variables
#' @export
create_input_drivers <- function() {

  # input files
  RMI_input_file <<- read.csv("Data/RMI_Melle.csv")
  RMI_radiation_input_file <<- read.csv("Data/RMI_radiation.csv", sep = ";")
  PE_input_file <<- "Data/Macro_temp_plant_eco.txt"
  pyr_input_file <<-  read.csv("Data/pyranometer_tower.dat", skip = 4, header = FALSE, stringsAsFactors = FALSE)
  TOMST_input_file <<- read_excel("Data/TOMST_hourly.xlsx")

  # plot parameters
  req_height <<- 1

  # spacetime specifics
  lat <<- 50.980
  lon <<- 3.816
  length_transect <<- 135
  height_canopy <<- 38

  # parameters for energy balance convergence
  energy_balance_tolerance <<- 2  # Maximum energy balance closure error between successive iteration steps to reach convergence (W/m2)

}

#' Generate (physical) constants
#'
#' @return (physical) constants
#' @export
create_physical_constants <- function() {

  # constants for longwave RTM
  sigma_SB <<- 5.67e-8  # Stefan-Boltzmann constant (W/m²/K⁴)

  # constants for air to air heat convection (C)
  voxel_length <<- 1    # Voxel edge length (m)
  Cp <<- 1000           # Specific heat air (J/kg/K)
  rho <<- 1.225         # Density air (kg/m3)

  # constants for ground heat flux (G) and soil temperature
  stable_soil_depth <<- 0.08 # Depth of measured soil temperature (assumed to be stable over 6-12h)

  # constants (empirical) for latent heat flux (LE)
  alpha_PT <<- 1.26         # Replaces aerodynamic terms from Penman-Monteith equation (unitless)
                            # 1.26 is the value for non-water-stressed conditions (arid regions lower, wetlands higher)
  gamma_psy <<- 0.066       # Psychrometric constant (kPa/K = kPa/°C)

}


#' Generate ForEdgeClim model parameters (to be calibrated)
#'
#' @return model parameters
#' @export
create_model_parameters <- function() {

  # parameters for shortwave RTM
  betad <<- 0.325          # Fraction of scattered diffuse radiation in backward direction
  beta0 <<- 0.325         # Fraction of scattered direct beam radiation in backward direction
  omega <<- 0.52        # Scattering coefficient
  # -> Vertical RTM
  Kd_v <<- 0.775          # Diffuse extinction coefficient for vertical radiation, per unit density
  Kb_v <<- 1.25         # Direct beam extinciton coefficient for vertical radiation, per unit density
  omega_g_v <<- 0.13    # Ground scattering
  # -> Horizontal RTM
  Kd_h <<- 0.725          # Diffuse extinciton coefficient for lateral radiation, per unit density
  Kb_h <<- 1.15          # Direct beam extinction coefficient for lateral radiation, per unit density
  omega_g_h <<- 0.15       # This is not ground scattering, but scattering by the inner forest. This is a balanced
  # value between reflection by leaves and transmission by open spaces.

  # parameters for longwave RTM
  e_forest <<- 0.965     # Emissivity forest (leaves and wood)
  beta_lw <<- 0.325      # Fraction of scattered longwave radiation in backward direction
  omega_lw <<- 0.035      # Scattering coefficient for longwave radiation
  # -> Vertical RTM
  Kd_lw_v <<- 0.3       # Longwave extinction coefficient for vertical radiation, per unit density
  omega_g_lw_v <<- 0.055 # Ground scattering for longwave radiation
  # -> Horizontal RTM
  Kd_lw_h <<- 0.3       # Longwave extinciton coefficient for lateral radiation, per unit density
  omega_g_lw_h <<- 0.035    # This is not ground scattering, but scattering by the inner forest. This is a balanced
  # value between reflection by leaves and transmission by open spaces.

  # parameters for air to air heat convection (C)
  h <<- 10              # Convection coefficient of air (W/m2/K), higher value means more wind/more turbulence

  # parameters for update linearisation air temperature
  g_macro <<- 25          # Convective heat transfer coefficient between (voxel) air and (macro) air (W/m2/K)
  infl_macro <<- 32.5       # Distance over which the influence of macro temp on air temp is reduced by 50% (m)
  infl_soil <<- 5         # Distance over which the influence of soil temp on air temp is reduced by 50% (m)
  infl_forest <<- 5       # Distance over which the influence of forest temp on air temp is reduced by 50% (m)

  # parameters for sensible heat flux (H)
  g_forest <<- 12.5   # Combined conductive & convective heat transfer coefficient between (voxel) air and structure (leaf) (W/m2/K)

  # parameters for ground heat flux (G) and soil temperature
  p_ground <<- 0.225          # Fraction of net ground radiation to define ground flux
  g_soil <<- 10              # Convective heat transfer coefficient between (voxel) air of air layer just above the ground and ground surface (W/m2/K)
  k_soil <<- 1.225             # Thermal conductance soil (W/m/K)

}

