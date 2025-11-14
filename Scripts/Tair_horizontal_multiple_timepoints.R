###############################################################################
#
# In this script multiple timepoints between model and observations are plot along the
# horizontal transect line.
#
###############################################################################

library(dplyr)
library(ggplot2)
library(gganimate)
library(readxl)
library(lubridate)
library(ggimage)
library(png)
library(grid)

#########
# INPUT #
#########

model_times <- lubridate::ymd_h(c("2023-07-08 01", "2023-07-08 08", "2023-07-08 12"))

res_list   <- list(
  readRDS("Output/2023-07-08_01/model_results.rds"),
  readRDS("Output/2023-07-08_08/model_results.rds"),
  readRDS("Output/2023-07-08_12/model_results.rds")
)

tomst_files <- c(
  "Output/2023-07-08_01/TOMST_filtered_distance_temp_20230708_0100.csv",
  "Output/2023-07-08_08/TOMST_filtered_distance_temp_20230708_0800.csv",
  "Output/2023-07-08_12/TOMST_filtered_distance_temp_20230708_1200.csv"
)

tomst_times <- model_times

req_height <- 1  # (m)
length_transect <- 135 # (m)

#########
# DATA  #
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

# Combine
anim_df <- bind_rows(model_df, tomst_df)

anim_df <- anim_df %>%
  mutate(time = floor_date(time, unit = "hour"),
         time_lbl = format(time, "%Y-%m-%d %Hh"))

label_df <- anim_df %>%
  distinct(time_lbl) %>%
  mutate(
    x = -Inf,
    y =  Inf
  )

anim_df$time_lbl <- factor(
  anim_df$time_lbl,
  levels = sort(unique(anim_df$time_lbl)),
  labels = c("01h (~night)", "08h (~morning)", "12h (~noon)")
)

####################
# BACKGROUND IMAGE #
####################

bg_img  <- png::readPNG("Output/horizontal_transect.png")
bg_grob <- grid::rasterGrob(bg_img, width = unit(1, "npc"), height = unit(1, "npc"))

############
#   PLOT   #
############

p_all <- ggplot(anim_df, aes(
  x = x, y = temperature,
  colour = time_lbl,
  linetype = model,
  shape = model,
  group = interaction(model, time_lbl)
)) +
  # background
  annotation_custom(
    grob = bg_grob,
    xmin = 0, xmax = length_transect,
    ymin = 15, ymax = 32
  ) +

  # Model-lines
  geom_line(
    data = dplyr::filter(anim_df, model == "Modelled (every 1m)"),
    linewidth = 1
  ) +

  # TOMST errorbars (±0.5 °C)
  geom_errorbar(
    data = dplyr::filter(anim_df, model == "TOMST observations (every 15m)"),
    inherit.aes = FALSE,
    aes(x = x,
        ymin = temperature - 0.5,
        ymax = temperature + 0.5,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    width = 3,
    linewidth = 0.6,
    alpha = 0.9,
    show.legend = FALSE
  ) +

  # TOMST points
  geom_point(
    data = dplyr::filter(anim_df, model == "TOMST observations (every 15m)"),
    size = 2.5
  ) +

  # scales and legends
  scale_colour_manual(
    name = "Time",
    values = c("01h (~night)"   = "blue",
               "08h (~morning)" = "orange",
               "12h (~noon)"    = "red")
  ) +
  scale_linetype_manual(
    name = "Source",
    values = c("Modelled (every 1m)" = "solid",
               "TOMST observations (every 15m)" = "blank")
  ) +
  scale_shape_manual(
    name = "Source",
    values = c("Modelled (every 1m)" = NA,
               "TOMST observations (every 15m)" = 16)
  ) +
  guides(
    colour   = guide_legend(order = 1, title = "Time", direction = "horizontal",override.aes = list(size = 3)),
    linetype = guide_legend(order = 2, title = "Source", direction = "horizontal", override.aes = list(linewidth  = 2)),
    shape    = guide_legend(order = 2, title = "Source", direction = "horizontal", override.aes = list(size = 3))
    ) +
  coord_cartesian(xlim = c(0, length_transect), ylim = c(15, 32)) +
  labs(x = "Distance from forest core (m)", y = "Temperature (°C)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position   = "top",
    legend.box        = "vertical",
    axis.title      = element_text(size = 35),
    axis.text       = element_text(size = 35),
    panel.background  = element_rect(fill = NA, colour = NA),
    legend.title      = element_text(size = 35),
    legend.text       = element_text(size = 35),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    plot.margin       = margin(10, 20, 10, 10)
  )

ggsave("Output/Tair_horizontal.png", p_all, width = 16, height = 9, dpi = 300)
