###############################################################################
#
# In this script multiple timepoints between model and observations are plot along the
# vertical transect line.
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
library(abind)

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
  "Output/2023-07-08_01/TOMST_filtered_height_temp_20230708_0100.csv",
  "Output/2023-07-08_08/TOMST_filtered_height_temp_20230708_0800.csv",
  "Output/2023-07-08_12/TOMST_filtered_height_temp_20230708_1200.csv"
)

tomst_times <- model_times

#########
# DATA  #
#########

# Model data in long-form (x == 75)
model_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z <= 35, y == 15, x == 75) %>%
      transmute(
        z           = z,
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
        z           = height,
        temperature = Tair,
        model       = "TOMST observations (every 7m)",
        time        = tomst_times[i]
      )
  })
)

anim_df <- bind_rows(model_df, tomst_df) %>%
  mutate(
    time     = lubridate::floor_date(time, unit = "hour"),
    hour_lbl = format(time, "%Hh"),
    time_lbl = factor(hour_lbl,
                      levels = c("01h", "08h", "12h"),
                      labels = c("01h (~night)", "08h (~morning)", "12h (~noon)"))
  )

temp_min <- floor(min(anim_df$temperature, na.rm = TRUE))
temp_max <- ceiling(max(anim_df$temperature, na.rm = TRUE))

####################
# BACKGROUND IMAGE #
####################

bg_img <- png::readPNG("Output/vertical_transect.png")

if (dim(bg_img)[3] == 4) {
  bg_img[,,4] <- bg_img[,,4] * 0.4  # transparent
} else {
  alpha_channel <- array(0.4, dim = c(dim(bg_img)[1], dim(bg_img)[2]))
  bg_img <- abind(bg_img, alpha_channel, along = 3)
}

bg_grob <- rasterGrob(bg_img, width = unit(1, "npc"), height = unit(1, "npc"))

temp_mid <- (temp_min + temp_max) / 2
width_range <- 8
xmin_bg <- temp_mid - width_range / 2
xmax_bg <- temp_mid + width_range / 2

############
#   PLOT   #
############

cap_h <- 0.35

p_vert <- ggplot() +
  # background
  annotation_custom(
    grob = bg_grob,
    xmin = xmin_bg,
    xmax = xmax_bg,
    ymin = 0,
    ymax = 35
  ) +

  # Model as path
  geom_path(
    data = dplyr::filter(anim_df, model == "Modelled (every 1m)"),
    aes(x = temperature, y = z, colour = time_lbl,
        linetype = model, shape = model, group = interaction(model, time_lbl)),
    linewidth = 1
  ) +

  # TOMST errorbars (±0.5 °C)
  geom_segment(
    data = dplyr::filter(anim_df, model == "TOMST observations (every 7m)"),
    inherit.aes = FALSE,
    aes(x = temperature - 0.5, xend = temperature + 0.5,
        y = z, yend = z,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    linewidth = 0.6, alpha = 0.9, show.legend = FALSE
  ) +
  geom_segment(
    data = dplyr::filter(anim_df, model == "TOMST observations (every 7m)"),
    inherit.aes = FALSE,
    aes(x = temperature - 0.5, xend = temperature - 0.5,
        y = z - cap_h, yend = z + cap_h,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    linewidth = 0.6, alpha = 0.9, show.legend = FALSE
  ) +
  geom_segment(
    data = dplyr::filter(anim_df, model == "TOMST observations (every 7m)"),
    inherit.aes = FALSE,
    aes(x = temperature + 0.5, xend = temperature + 0.5,
        y = z - cap_h, yend = z + cap_h,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    linewidth = 0.6, alpha = 0.9, show.legend = FALSE
  ) +

  # TOMST-points
  geom_point(
    data = dplyr::filter(anim_df, model == "TOMST observations (every 7m)"),
    aes(x = temperature, y = z, colour = time_lbl,
        linetype = model, shape = model, group = interaction(model, time_lbl)),
    size = 2.5
  ) +
  scale_colour_manual(
    name = "Time",
    values = c("01h (~night)" = "blue",
               "08h (~morning)" = "orange",
               "12h (~noon)"    = "red")
  ) +
  scale_linetype_manual(
    name = "Source",
    values = c("Modelled (every 1m)" = "solid",
               "TOMST observations (every 7m)"  = "blank")
  ) +
  scale_shape_manual(
    name = "Source",
    values = c("Modelled (every 1m)" = NA,
               "TOMST observations (every 7m)"  = 16)
  ) +
  guides(
    colour   = guide_legend(order = 1, title = "Time", direction = "horizontal"),
    linetype = guide_legend(order = 2, title = "Source", direction = "horizontal"),
    shape    = guide_legend(order = 2, title = "Source", direction = "horizontal")
  ) +
  coord_cartesian(ylim = c(0, 35), xlim = c(temp_min, temp_max)) +
  labs(x = "Temperature (°C)", y = "Height (m)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position   = "top",
    legend.box        = "vertical",
    axis.title      = element_text(size = 35),
    axis.text       = element_text(size = 35),
    legend.title      = element_text(size = 35),
    legend.text       = element_text(size = 35),
    panel.background  = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )

ggsave("Output/Tair_vertical.png", p_vert, width = 16, height = 18, dpi = 300)
