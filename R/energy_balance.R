#' Function to calculate the derivative of the energy balance closure error to surface temperature
#'
#' @param temp Surface temperature
#' @return Derivative of energy balance wrt surface temperature
#' @export
#'
de_dT <- function(temp, net_rad){
  # derivative of net radiation to SURF temp:
  dR_dT = -4*e_forest*sigma_SB*temp^3
  # derivative of sensible heat transfer to SURF temp:
  dH_dT = den*g_forest
  # derivative of latent heat flux to SURF temp:
  T_c = temp - 273.15 # temperature in degrees Celsius
  s = 4098 * saturated_vapor_pressure(T_c) / (T_c + 237.3)^2 # slope in Priestly-Taylor equation
  des_dT = 0.6108 * exp(17.27*T_c/(T_c+237.3)) * ((T_c+237.3)*17.27 - 17.27*T_c)/(T_c+237.3)^2 # derivative of saturated vapor pressure to temperature
  ds_dT = (4098*des_dT*(T_c+237.3)^2 - 2*(T_c+237.3)*4098*saturated_vapor_pressure(T_c))/(T_c+237.3)^4 # derivative of slope to temperature
  df_dT = (ds_dT*(s+gamma_psy)-s*ds_dT)/(s+gamma_psy)^2 # derivative of s/(s+gamma_psy) to temperature
  dLE_dT = den*alpha_PT * (dR_dT)* s / (s + gamma_psy) + den*alpha_PT * (net_rad) * df_dT
  # print('de_dT components')
  # print(summary(dR_dT))
  # print(summary(dH_dT))
  # print(summary(dLE_dT))

  return(dR_dT - dH_dT - dLE_dT)

}
