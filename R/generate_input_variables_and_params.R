#' Generate standard input data for ForEdgeClim
#'
#' @param summer_day Logical value: TRUE for day, FALSE for night
#' @return global variables and parameters
#' @export
create_input_data <- function(summer_day = TRUE) {

  if (summer_day) {
    datetime <- as.POSIXct("2023-06-23 12:00:00", tz = "CET")
    F_sky_dir_init <- 600         # Direct solar beam radiation (W/m2); summer 400, autumn 200, winter 100, spring 275
    F_sky_diff_init <- 200        # Diffuse radiation (W/m2); summer 75, autumn 75, winter 55, spring 80
    F_sky_lw <- 400               # Longwave radiation (W/m2); summer D 375 - N 275
    macro_temp <- 25 + 273.15     # macrotemperature (K); summer D 25 - N 17, autumn D 13 - N 7, winter D 5 - N -2, spring D 17 - N 11
  } else {
    datetime <- as.POSIXct("2023-06-23 02:00:00", tz = "CET")
    F_sky_dir_init <- 0
    F_sky_diff_init <- 0
    F_sky_lw <- 275
    macro_temp <- 17 + 273.15
  }


  # Save al input variables as global variables (<<-)
  datetime <<- datetime
  F_sky_dir_init <<- F_sky_dir_init
  F_sky_diff_init <<- F_sky_diff_init
  F_sky_lw <<- F_sky_lw
  macro_temp <<- macro_temp
  # spacetime specifics
  #datetime = as.POSIXct("2023-06-23 12:00:00", tz = "CET")
  lat <<- 50.980
  lon <<- 3.816

  # macro environment variables
  # F_sky_dir_init <- 400          # Direct solar beam radiation (W/m2); summer 400, autumn 200, winter 100, spring 275
  # F_sky_diff_init <- 75          # Diffuse radiation (W/m2); summer 75, autumn 75, winter 55, spring 80
  # macro_temp <- 25 + 273.15      # macrotemperature (K); summer D 25 - N 17, autumn D 13 - N 7, winter D 5 - N -2, spring D 17 - N 11
  macro_temp_max <<- 26 + 273.15   # maximum daily or seasonal macro temperature (for lagged air temp just above ground)
  macro_temp_min <<- 16 + 273.15   # minimum daily or seasonal macro temperature (for lagged air temp just above ground)
  t_max <<- 13                     # Time point (hour) of max macro temp

  # parameters for shortwave RTM
  betad <<- 0.3         # Fraction of scattered diffuse radiation in backward direction
  beta0 <<- 0.2         # Fraction of scattered direct beam radiation in backward direction
  omega <<- 0.85        # Scattering coefficient
  # -> Vertical RTM
  Kd_v <<- 0.3          # 0.3 Diffuse extinction coefficient for vertical radiation, per unit density
  Kb_v <<- 0.75         # 0.75 Direct beam extinciton coefficient for vertical radiation, per unit density
  omega_g_v <<- 0.2     # 0.2 Ground scattering
  # -> Horizontal RTM
  Kd_h <<- 0.2          # 0.2 Diffuse extinciton coefficient for lateral radiation, per unit density
  Kb_h <<- 0.3          # 0.3 Direct beam extinction coefficient for lateral radiation, per unit density
  omega_g_h <<- 0       # This is not ground scattering, but scattering by the inner forest. This is a balanced
                        # value between reflection by leaves and transmission by open spaces.

  # parameters for longwave RTM
  sigma_SB <<- 5.67e-8  # Stefan-Boltzmann constant (W/m²/K⁴)
  e_atm <<- 0.2         # Emissivity atmosphere (0.9 on cloudy days)
  e_forest <<- 0.95     # Emissivity forest (leaves and wood)
  beta_lw <<- 0.03      # Fraction of scattered diffuse radiation in backward direction
  omega_lw <<- 0.1      # Scattering coefficient
  # -> Vertical RTM
  Kd_lw_v <<- 0.2       # 0.2 Diffuse extinction coefficient for vertical radiation, per unit density
  omega_g_lw_v <<- 0.02 # Ground scattering
  # -> Horizontal RTM
  Kd_lw_h <<- 0.1       # 0.1 Diffuse extinciton coefficient for lateral radiation, per unit density
  omega_g_lw_h <<- 0    # This is not ground scattering, but scattering by the inner forest. This is a balanced
                        # value between reflection by leaves and transmission by open spaces.

  # parameters for air to air heat convection (C)
  h <<- 15              # Convection coefficient of air (W/m2/K), higher value means more wind/more turbulence
  voxel_length <<- 1    # Voxel edge length (m)
  Cp <<- 1000           # Specific heat air (J/kg/K)
  rho <<- 1.225         # Density air (kg/m3)

  # parameters for update linearisation air temperature
  g_macro <<- 10          # Convective heat transfer coefficient between (voxel) air and (macro) air (W/m2/K)
  g_soil <<- 10           # Combined conductive & convective heat transfer coefficient between air and ground (W/m2/K)
  infl_macro <<- 30       # Distance over which the influence of macro temp on air temp is reduced by 50% (m)
  infl_soil <<- 5         # Distance over which the influence of soil temp on air temp is reduced by 50% (m)
  infl_forest <<- 6       # Distance over which the influence of forest temp on air temp is reduced by 50% (m)

  # parameters for sensible heat flux (H)
  g_forest <<- 10   # Combined conductive & convective heat transfer coefficient between air and structure (leaf) (W/m2/K)

  # parameters for ground heat flux (G)
  p_ground <<- 0.10   # Fraction of net ground radiation to define ground flux

  # parameters for latent heat flux (LE)
  alpha_PT <<- 1.26         # Replaces aerodynamic terms from Penman-Monteith equation (unitless)
                            # 1.26 is the value for non-water-stressed conditions (arid regions lower, wetlands higher)
  gamma_psy <<- 0.066       # Psychrometric constant (kPa/K = kPa/°C)

  # parametes for energy balance convergence
  energy_balance_tolerance <<- 10  # Maximum energy balance closure error between successive iteration steps to reach convergence (W/m2)

}
