#' Function to import RMI national weather station observations (macro
#' temperature and incoming longwave radiation)
#'
#' @return macro temperature and downward from sky longwave radiation
#' @importFrom lubridate ymd_hms ymd_hm floor_date
#' @importFrom dplyr mutate group_by summarise
#' @export
import_RMI_observations <- function(datetime){

  # MACRO TEMPERATURE
  #------------------
  RMI_data <- RMI_input_file

  # Set timestamp in RMI_data as a POSIXct (RMI data always is in UTC)
  RMI_data$timestamp <- ymd_hms(RMI_data$timestamp, tz = "UTC")

  # Extract hour of interest
  RMI_hour <- RMI_data[RMI_data$timestamp == datetime, ]

  macro_temp <<- RMI_hour$temp + 273.15 # in Kelvin

  # LONGWAVE RADIATION FROM SKY
  #----------------------------
  RMI_radiation <- RMI_radiation_input_file

  # Set timestamp in RMI_data as a POSIXct (RMI data always is in UTC)
  RMI_radiation <- RMI_radiation |>
    mutate(DATE = ymd_hm(DATE, tz = "UTC")) |>
    mutate(DATE_hour = floor_date(DATE, unit = "hour")) |>
    group_by(DATE_hour) |>
    summarise(IR_mean = mean(IR_FROM_SKY_AVG..W.m2.min., na.rm = TRUE))

  # Extract hour of interest
  RMI_radiation <- RMI_radiation[RMI_radiation$DATE_hour == datetime, ]

  F_sky_lw <<- RMI_radiation$IR_mean # in W/m2

}

#' Function to import macro temperature from Plant Ecology lab
#'
#' @return macro temperature and downward longwave radiation
#' @importFrom dplyr mutate filter summarise
#' @importFrom lubridate dmy_hm with_tz hours
#' @importFrom readr read_delim col_character locale cols col_double
#' @export
import_PE_observations <- function(datetime){

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
import_pyr_observations <- function(datetime){

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

#' Function to import TOMST observations (soil temperature at tower and
#' observations along the transect)
#'
#' @return T_soil_deep The soil temperature at 8cm deep at the position of the tower
#' @export
import_soil_temperature <- function(datetime){

  TOMST_hourly <- TOMST_input_file

  # Set 'datehour' as POSIXct-object
  TOMST_hourly$datehour_posix <- as.POSIXct(TOMST_hourly$datehour, format = "%Y.%m.%d %H", tz = "UTC")

  # Filter on datetime and height
  filtered_data <- subset(TOMST_hourly, datehour_posix == datetime & height == 0)
  # Filter only on datetime
  filtered_data_vertical <- subset(TOMST_hourly, datehour_posix == datetime)

  # Soil temperature as mean of all (20) sensors within the forest soil
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
  write.csv(TOMST_air_output, paste0("Data/TOMST_filtered_distance_temp_", format(datetime, "%Y%m%d_%H%M"), ".csv"), row.names = FALSE)
  write.csv(TOMST_air_output_vertical, paste0("Data/TOMST_filtered_height_temp_", format(datetime, "%Y%m%d_%H%M"), ".csv"), row.names = FALSE)

}


