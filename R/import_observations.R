#' Function to import hourly DTS data observations
#'
#' @return Filtered CSV file
#' @importFrom dplyr mutate
#' @importFrom lubridate ymd_hms with_tz
#' @export
import_DTS_observations <- function(){

  DTS_data <- DTS_input_file

  # DTS data is collected in CET
  DTS_data$datetime_brussels <- ymd_hms(paste(DTS_data$date, DTS_data$hour, "00", "00"), tz = "CET")

  # Convert to UTC
  DTS_data$datetime_utc <- with_tz(DTS_data$datetime_brussels, "UTC")

  # Filter distances between 82 and 219 m (the transect line length, at 1m height from west to east)
  DTS_data_filtered <- DTS_data[DTS_data$distance >= 82 & DTS_data$distance <= 219, ]

  # Transform distance to 0 - 135 m scale
  DTS_data_filtered$distance <- (DTS_data_filtered$distance - 82) / (219 - 82) * 135

  # Reverse scaled distance from east to west
  DTS_data_filtered$distance <- 135 - DTS_data_filtered$distance

  # # Between 49 (52) and 62 (58) m there is a wooden construction that heats up the DTS fiber, we will remove these datapoints
  # DTS_data_filtered <- DTS_data_filtered |>
  #   mutate(temp = ifelse(distance >= 49 & distance <= 62, NA, temp))

  # Extract datetime of interest and output columns
  DTS_data_output <- DTS_data_filtered[DTS_data_filtered$datetime_utc == datetime, c("distance", "temp")]

  # Save dataframe as CSV
  write.csv(DTS_data_output, "Data/DTS_filtered_distance_temp.csv", row.names = FALSE)
}


#' Function to import RMI national weather station observation (macro temperature)
#'
#' @return macro temperature and downward longwave radiation
#' @importFrom lubridate ymd_hms
#' @export
import_RMI_observations <- function(){

  RMI_data <- RMI_input_file

  # Set timestamp in RMI_data as a POSIXct (RMI data always is in UTC)
  RMI_data$timestamp <- ymd_hms(RMI_data$timestamp, tz = "UTC")

  # Extract hour of interest
  RMI_hour <- RMI_data[RMI_data$timestamp == datetime, ]

  # macro_temp <<- RMI_hour$temp + 273.15 # in Kelvin
  #
  # F_sky_lw <<- e_atm*sigma_SB*macro_temp^4 # in W/m2

}

#' Function to import Delta-T pyranometer observations (direct and diffuse light)
#'
#' @return Direct solar beam and diffuse radiation
#' @importFrom lubridate ymd_hms hours hour
#' @export
import_pyr_observations <- function(){

  pyr_data <- pyr_input_file

  # Add column names
  colnames(pyr_data) <- c("datetime", "record", "battery", "temp", "total", "diffuse")

  # Set datetime as a POSIXct (without forcing a time zone)
  pyr_data$datetime <- ymd_hms(pyr_data$datetime)

  # pyr_data is always in summer time (CEST = UTC + 2) => reduce by 2h to set into UTC
  pyr_data$datetime_utc <- pyr_data$datetime - hours(2)

  # Get date and hour from pyr_data$datetime_utc
  pyr_data$date <- as.Date(pyr_data$datetime_utc)
  pyr_data$hour <- hour(pyr_data$datetime_utc)

  # Get date and hour from datetime input parameter
  datetime_date <- as.Date(datetime)
  datetime_hour <- hour(datetime)

  # Get light values for corresponding date and hour
  # pyr_data has measurements for every 15 min, so 4 measurements every hour => we will use the hourly-mean value (see below)
  pyr_filtered <- pyr_data[pyr_data$date == datetime_date & pyr_data$hour == datetime_hour, ]

  # Conversion from mV to W/m2 by factor 0.5
  F_sky_dir_init <<- mean(pyr_filtered$total - pyr_filtered$diffuse)/2    # Direct solar beam radiation (W/m2)
  F_sky_diff_init <<- mean(pyr_filtered$diffuse)/2                        # Diffuse radiation (W/m2)
  macro_temp <<-  mean(pyr_filtered$temp) + 273.15 # in Kelvin
  F_sky_lw <<- e_atm*sigma_SB*macro_temp^4 # in W/m2

}

#' Function to import TOMST observations (soil temperature at tower)
#'
#' @return T_soil_deep The soil temperature at 6cm deep at the position of the tower
#' @export
import_soil_temperature <- function(){

  TOMST_hourly <- TOMST_input_file

  # Set 'datehour' as POSIXct-object
  TOMST_hourly$datehour_posix <- as.POSIXct(TOMST_hourly$datehour, format = "%Y.%m.%d %H", tz = "UTC")

  # Filter on datetime and name = "C75"
  filtered_data <- subset(TOMST_hourly, datehour_posix == datetime & name == "C75")

  T_soil_deep <<- filtered_data$Tsoi + 273.15

}
