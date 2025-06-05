library(dplyr)
library(lubridate)
library(ggplot2)
library(readxl)

#########
# INPUT #
#########

model_times <- lubridate::ymd_h(c("2023-07-08 00", "2023-07-08 04", "2023-07-08 08",
                                  "2023-07-08 12", "2023-07-08 16", "2023-07-08 20"))

res_list   <- list(
  readRDS("Data/data_for_animation/model_results_0h_vox.rds"),
  readRDS("Data/data_for_animation/model_results_4h_vox.rds"),
  readRDS("Data/data_for_animation/model_results_8h_vox.rds"),
  readRDS("Data/data_for_animation/model_results_12h_vox.rds"),
  readRDS("Data/data_for_animation/model_results_16h_vox.rds"),
  readRDS("Data/data_for_animation/model_results_20h_vox.rds")
)

tomst_files <- c(
  "Data/data_for_animation/TOMST_filtered_distance_temp_0h.csv",
  "Data/data_for_animation/TOMST_filtered_distance_temp_4h.csv",
  "Data/data_for_animation/TOMST_filtered_distance_temp_8h.csv",
  "Data/data_for_animation/TOMST_filtered_distance_temp_12h.csv",
  "Data/data_for_animation/TOMST_filtered_distance_temp_16h.csv",
  "Data/data_for_animation/TOMST_filtered_distance_temp_20h.csv"
)

tomst_times <- model_times

req_height = 1 # (m)
length_transect = 135 # (m)


####################################
# DATA FROM MODEL AND OBSERVATIONS #
####################################

# Model core data in long-form
model_core_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 15, x <= min(x)) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        id       = "Model-Core",
        Source       = "Model",
        Position = "Core",
        time        = model_times[i]
      )
  })
)

# Model edge data in long-form
model_edge_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 15, x == length_transect) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        id       = "Model-Edge",
        Source       = "Model",
        Position = "Edge",
        time        = model_times[i]
      )
  })
)

# Macro data in long-form
macro_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == max(z), y == 15, x == max(x)) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        id       = "Macro",
        Source       = "Model",
        Position = "Macro",
        time        = model_times[i]
      )
  })
)

# TOMST core data in long-form
tomst_core_df <- bind_rows(
  lapply(seq_along(tomst_files), function(i) {
    read.csv(tomst_files[i]) %>%
      transmute(
        x           = 135 - D_edge,
        temperature = Tair,
        id       = "TOMST-Core",
        Source       = "TOMST",
        Position = "Core",
        time        = tomst_times[i]
      ) %>%
      filter(x == 0)
  })
)

# TOMST core data in long-form
tomst_edge_df <- bind_rows(
  lapply(seq_along(tomst_files), function(i) {
    read.csv(tomst_files[i]) %>%
      transmute(
        x           = 135 - D_edge,
        temperature = Tair,
        id       = "TOMST-Edge",
        Source       = "TOMST",
        Position = "Edge",
        time        = tomst_times[i]
      ) %>%
      filter(x == 135)
  })
)


# combine dataframes
combined_df <- bind_rows(model_core_df, model_edge_df, tomst_core_df, tomst_edge_df, macro_df)
combined_df$time <- as.POSIXct(combined_df$time, format = "%Y-%m-%d %H:%M:%S")  # Pas dit format aan indien nodig
combined_df$Position <- factor(combined_df$Position, levels = c("Macro", "Edge", "Core"))



##############################################################
# PLOTTING TIMESERIES CORE VS EDGE AND MODEL VS OBSERVATIONS #
##############################################################

timeseries = ggplot(combined_df, aes(x = time, y = temperature, color = Position, linetype = Source)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +  # Loess curve without confidence interval
  labs(title = "Temperature time fluctuations",
       x = "Time",
       y = "Temperature (°C)",
       color = "Position",
       linetype = "Source") +
  scale_color_manual(
    values = c("Macro" = "red4",
                "Core" = "cornflowerblue",
                "Edge" = "orange"
               )) +
  scale_linetype_manual(values = c(
    "Model" = "solid",
    "TOMST" = "dotted")) +
  guides(
    linetype = guide_legend(override.aes = list(color = "black"))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
print(timeseries)

ggsave('Output/timeseries_core_vs_edge_AND_model_vs_observations.png', plot = timeseries, width = 9, height = 6, dpi = 300)
