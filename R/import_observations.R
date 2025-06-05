#' Function to import hourly DTS data observations
#'
#' @return Filtered CSV file
#' @importFrom dplyr mutate
#' @importFrom lubridate ymd_hms with_tz
#' @export
import_DTS_observations <- function(){

  DTS_data <- DTS_input_file

  # DTS data is collected in CET = UTC+1
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

  macro_temp <<- RMI_hour$temp + 273.15 # in Kelvin

  e_atm <<- 9.37e-6*macro_temp^2 # formula of Swinbank, mainly for clear sky (when only having macro_temp as input, no cloudiness...)

  F_sky_lw <<- e_atm*sigma_SB*macro_temp^4 # in W/m2

}

#' Function to import macro temperature from Plant Ecology lab
#'
#' @return macro temperature and downward longwave radiation
#' @importFrom dplyr mutate filter summarise
#' @importFrom lubridate dmy_hm with_tz hours
#' @importFrom readr read_delim col_character locale cols col_double
#' @export
import_PE_observations <- function(){

  # read txt file
  df <- read_delim(
    PE_input_file,
    delim     = "\t",
    skip      = 1,
    col_names = c("Date", "Time", "Temperature"),
    locale    = locale(decimal_mark = ","),
    col_types = cols(
      Date        = col_character(),
      Time        = col_character(),
      Temperature = col_double()
    )
  )

  # merge Date and Time to posixct object in Brussel's time zone
  df <- df |>
    mutate(
      timestamp_local = dmy_hm(paste(Date, Time), tz = "Europe/Brussels"),
      timestamp_utc = with_tz(timestamp_local, tzone = "UTC") # convert to UTC
    )

  target_utc <- datetime

  # filter on 1h interval in UTC
  result <- df |>
    filter(timestamp_utc >= target_utc,
           timestamp_utc <  target_utc + hours(1)) |>
    summarise(
      MeanTemperature = mean(Temperature, na.rm = TRUE)
    )

  macro_temp <<- result$MeanTemperature + 273.15 # in Kelvin

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
  # macro_temp <<-  mean(pyr_filtered$temp) + 273.15 # in Kelvin


}

#' Function to import TOMST observations (soil temperature at tower)
#'
#' @return T_soil_deep The soil temperature at 6cm deep at the position of the tower
#' @export
import_soil_temperature <- function(){

  TOMST_hourly <- TOMST_input_file

  # Set 'datehour' as POSIXct-object
  TOMST_hourly$datehour_posix <- as.POSIXct(TOMST_hourly$datehour, format = "%Y.%m.%d %H", tz = "UTC")

  # Filter on datetime and bodemhoogte
  filtered_data <- subset(TOMST_hourly, datehour_posix == datetime & height == 0)
  # Filter on datetime
  filtered_data_vertical <- subset(TOMST_hourly, datehour_posix == datetime)

  T_soil_deep <<- mean(filtered_data$Tsoi) + 273.15

  # Extract TOMST sensors of interest and output columns for horizontal measurements
  TOMST_air_output <- filtered_data[grepl("^C", filtered_data$name),  c("D_edge", "Tair", 'Tsoi')]

  # Extract TOMST sensors of interest and output columns for vertical measurements
  TOMST_air_output_vertical <- filtered_data_vertical[grepl("^V", filtered_data_vertical$name) | filtered_data_vertical$name == "C75",
                                                      c("Tair", "height")]



  # To define multiple T_soil_deep along the transect
  # x_new = 1:135
  # Tsoi_interp <- data.frame(
  #   x = x_new,
  #   Tsoi = approx(x = 135 - TOMST_air_output$D_edge, y = TOMST_air_output$Tsoi, xout = x_new)$y
  # )
  #
  # # Repeat the Tsoi value closest to the edge 15 times to cover the soil between forest edge and max X (over the street)
  # T_soil_deep <<- c(Tsoi_interp$Tsoi, rep(Tsoi_interp$Tsoi[length(Tsoi_interp$Tsoi)], 15)) + 273.15

  # Save dataframes as CSV
  write.csv(TOMST_air_output, "Data/TOMST_filtered_distance_temp.csv", row.names = FALSE)
  write.csv(TOMST_air_output_vertical, "Data/TOMST_filtered_height_temp.csv", row.names = FALSE)

}

#' Function to import Odyssey PAR observations (PAR radiation)
#'
#' @return PAR radiation
#' @export
import_PAR_observations <- function(){

  PAR_data <- PAR_input_file

  # Set datehour as a POSIXct, time zone of PAR is local, Brussels time (switching between CET and CEST)
  PAR_data$datehour <- as.POSIXct(PAR_data$datehour, format = "%Y-%m-%d %H", tz = "Europe/Brussels")

  # Convert to UTC
  PAR_data$timestamp <- format(PAR_data$datehour, tz = "UTC", usetz = TRUE)
  PAR_data$timestamp <- as.POSIXct(PAR_data$timestamp, tz = "UTC")

  # Filter on datetime and only sensors at ground level
  filtered_data <- subset(PAR_data, timestamp == datetime & height == 0 & grepl("^C", name))

  # Reverse distance from east to west
  filtered_data$distance <- 135 - filtered_data$D_edge

  # Rescale PAR value in micromol/m2/s to total solar spectrum value in W/m2
  # unit conversion = 0.217 (average energy value of photons in PAR range)
  # spectrum range conversion = 0.45
  filtered_data$rad <- filtered_data$PAR * 0.217 / 0.45


  # Extract datetime of interest and output columns
  PAR_output <- filtered_data[,  c("distance", "rad")]

  # Save dataframe as CSV
  write.csv(PAR_output, "Data/PAR_filtered_distance_rad.csv", row.names = FALSE)


}

