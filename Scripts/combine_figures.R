library(magick)

# Lees de afbeeldingen in

# validation plots
#-----------------
img_val_lr <- image_read("Output/validation/plots/model_vs_obs_over_time_and_position_winter_both.png")
img_val_res <- image_read("Output/validation/plots/time_series_residual_winter_both.png")

# Maak een witte tussenruimte van 200 pixels breed (zelfde hoogte als plots)
spacer_val <- image_blank(height = 500, width = image_info(img_val_lr)$width, color = "white")

combined_val <- image_append(c(img_val_lr, spacer_val, img_val_res), stack = TRUE)

# Opslaan
image_write(combined_val, path = "Output/combined_validation.png", format = "png")


# calibration plots
#------------------
img_cal_rmse <- image_read("Output/calibration/plots/RMSE_spring.png")
img_cal_pca <- image_read("Output/calibration/plots/PCA_space_spring.png")

spacer_cal <- image_blank(height = 500, width = image_info(img_cal_rmse)$width, color = "white")

combined_cal <- image_append(c(img_cal_rmse, spacer_cal, img_cal_pca), stack = TRUE)

image_write(combined_cal, path = "Output/combined_calibration.png", format = "png")


# sensitivity plots
#------------------
img_SA_Tair <- image_read("Output/sensitivity_analysis/Sobol_QoI/plots_output/Sobol_indices_airT_horizontal.png")
img_SA_Tsurf <- image_read("Output/sensitivity_analysis/Sobol_QoI/plots_output/Sobol_indices_forest_surfaceT_horizontal.png")

spacer_SA <- image_blank(height = image_info(img_SA_Tair)$height, width = 500, color = "white")

combined_SA <- image_append(c(img_SA_Tair, spacer_SA, img_SA_Tsurf))

image_write(combined_SA, path = "Output/combined_sensitivity.png", format = "png")


# flux/temperature plots
#-----------------------
img_SW <- image_read("Output/2023-07-08_12/SW_2D_F_b_down.png")
img_Tair <- image_read("Output/2023-07-08_12/temp_air.png")
img_Tsurf <- image_read("Output/2023-07-08_12/temp_surface.png")

spacer_tiles <- image_blank(height = 0, width = image_info(img_SW)$width, color = "white")

combined_tiles <- image_append(c(img_SW, spacer_tiles, img_Tsurf, spacer_tiles, img_Tair), stack = TRUE)

image_write(combined_tiles, path = "Output/combined_tiles.png", format = "png")
