library(dplyr)
library(ggplot2)
library(gganimate)
library(readxl)
library(lubridate)
library(ggimage)



#########
# INPUT #
#########

model_times <- lubridate::ymd_h(c("2023-07-08 00", "2023-07-08 04", "2023-07-08 08",
                                "2023-07-08 12", "2023-07-08 16", "2023-07-08 20"))

res_list   <- list(
  readRDS("Data/data_for_animation/model_results_0h.rds"),
  readRDS("Data/data_for_animation/model_results_4h.rds"),
  readRDS("Data/data_for_animation/model_results_8h.rds"),
  readRDS("Data/data_for_animation/model_results_12h.rds"),
  readRDS("Data/data_for_animation/model_results_16h.rds"),
  readRDS("Data/data_for_animation/model_results_20h.rds")
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

#########

# Model data in long-form
model_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 15, x <= length_transect) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        model       = "Modelled (every 1m)",
        time        = model_times[i]
      )
  })
)

# TOMST data in long-form
tomst_df <- bind_rows(
  lapply(seq_along(tomst_files), function(i) {
    read.csv(tomst_files[i]) %>%
      transmute(
        x           = 135 - D_edge,
        temperature = Tair,
        model       = "TOMST observations (every 15m)",
        time        = tomst_times[i]
      )
  })
)

# anim_df with columns x, temperature, model, time (POSIXct)
anim_df <- bind_rows(model_df, tomst_df)


# time steps at the exact hour
anim_df <- anim_df %>%
  mutate(time = floor_date(time, unit = "hour"),
         time_lbl = format(time, "%Y-%m-%d %Hh"))

# data frame with corner coordinates and labels
label_df <- anim_df %>%
  distinct(time_lbl) %>%
  mutate(
    x = -Inf,       # left edge panel
    y =  Inf        # upper edge panel
  )

# # data frame for position of sun and moon image
# pos_df <- anim_df %>%
#   distinct(time, time_lbl) %>%
#   mutate(
#     hour   = hour(time),
#     is_day = hour >= 6 & hour < 18,
#     frac   = case_when(
#       is_day     ~ (hour - 6) / 12,
#       hour >= 18 ~ (hour - 18) / 12,
#       TRUE       ~ (hour + 6) / 12
#     ),
#     x_sky = (1 - frac) * max(anim_df$x),
#     y_sky = 15 + sin(frac * pi) * (32 - 15),
#     image = ifelse(is_day, "Data/data_for_animation/sun.png", "Data/data_for_animation/moon.png")
#   )

# animation plot
p_anim <- ggplot(anim_df, aes(x = x, y = temperature, color = model, group = model)) +
  geom_line(data  = filter(anim_df, model == "Modelled (every 1m)")) +
  geom_line(data  = filter(anim_df, model == "TOMST observations (every 15m)")) +
  geom_point(data = filter(anim_df, model == "TOMST observations (every 15m)")) +

  # time label
  geom_text(
    data        = label_df,
    aes(x = x, y = y, label = time_lbl),
    inherit.aes = FALSE,
    hjust       = -0.05,
    vjust       = 1.5,
    size        = 10
  ) +

  # # dynamic sun and moon
  # geom_image(
  #   data = pos_df,
  #   aes(x = x_sky, y = y_sky, image = image),
  #   inherit.aes = FALSE,
  #   size = 0.1,
  #   by = "width"
  # ) +

  scale_color_manual(values = c(
    "Modelled (every 1m)"            = "cornflowerblue",
    "TOMST observations (every 15m)" = "darkgreen"
  )) +
  coord_cartesian(xlim = c(0, max(anim_df$x)), ylim = c(15, 32)) +
  labs(x = "Distance (m)", y = "Temperature (°C)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position      = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 12)
  ) +

  # animation states
  transition_states(
    states           = time_lbl,
    transition_length = 0,
    state_length      = 1,
    wrap              = FALSE
  ) +
  ease_aes('linear')

# render and save
n_states <- length(unique(anim_df$time_lbl))
anim <- animate(
  p_anim,
  nframes  = n_states,
  fps      = 0.5,
  width    = 1600,
  height   = 800,
  res      = 150,
  renderer = gifski_renderer()
)
anim_save("Output/temperature_sunmoon.gif", anim)
