# Title: Remotely Sensed Imagery for ENMs
# Author: Mariana Castaneda-Guzman
# Date created: 9/08/2023
# Date last updated: 3/10/2023

# Description: Remotely sensed imagery form 1950 to 2023 for sea surface
# temperature and salinity


# Packages ----------------------------------------------------------------

library(terra)
library(sf)
library(ggplot2)
library(ggpubr)

# Specify databases -------------------------------------------------------

# Database: MODIS 2003 - 2023
db <- "(2003-2023)"

# Study are extent --------------------------------------------------------

usa <- sf::read_sf("../../00 Global Data/USA_adm/USA_adm1.shp")
head(usa)

eastern_coastal_states <- c("Connecticut", "Delaware", "Florida",
                            "Georgia", "Maine", "Maryland", 
                            "Massachusetts", "New Hampshire", "New Jersey", 
                            "New York", "North Carolina", "Rhode Island", 
                            "South Carolina", "Virginia")

eastern_usa <- subset(usa, NAME_1 %in% eastern_coastal_states)
# plot(eastern_usa)

# Westernmost part of Florida: lat: 30.867218° long: -87.642057°
# Easternmost part of Maine: lat: 44.815109° long: -66.950473°

# Southernmost part of Florida: lat: 25.040514° long:-80.833716°
# Added a little more south to include some ocean 23.140514
# Northernmost part of Maine: lat: 47.459574° long: -69.224743°


study_area_extent <- st_bbox(c(xmin = -87.642057, xmax = -66.950473, 
                               ymin = 23.140514, ymax = 47.459574))

usa_crop <- st_crop(usa, study_area_extent)

plot(crop(sst_cropped[[1]], study_area_extent))
plot(usa_crop, add = T)
plot(eastern_usa, add = T, col = "lightgrey")


# Correlation Analysis ----------------------------------------------------

# Make composites of years
if(db == "(2003-2023)") var <- c("SST", "CHLO")

summary_stats <- c("mean", "min", "max", "sd", "range")
summary_stats_fun <- c(mean, min, max, sd, range)

for(i in var){
  for(s in 1:length(summary_stats)){
    ras_stat <- lapply(list.files(path = paste0("../Data/RS_", db,"/", i), 
                                  pattern = paste0(".*", summary_stats[s], "\\.tif$"),
                                  recursive = TRUE, 
                                  full.names = TRUE),
                       rast)
    
    ras_stat_stack <- terra::rast(ras_stat)
    
    if(summary_stats[s] == "range"){
      ras_stat_min <- terra::app(ras_stat_stack, fun = summary_stats_fun[[2]], na.rm = TRUE)
      ras_stat_max <- terra::app(ras_stat_stack, fun = summary_stats_fun[[3]], na.rm = TRUE)
      
      ras_stat <- ras_stat_max - ras_stat_min   
      
    }else{
      ras_stat <- terra::app(ras_stat_stack, 
                             fun = summary_stats_fun[[s]], na.rm = T)
    }
    
    stack_name <- paste0("../Data/RS_", db, "/", i, "/", 
                         tolower(i), "_", "cropped_", 
                         summary_stats[s], "_summary.tif")
    
    terra::writeRaster(x = ras_stat, 
                filename = stack_name, 
                overwrite = TRUE)
    
  }
}


# Correlation Plots -------------------------------------------------------


# For SST
# Read summary rasters
sst_files <- list.files(path = paste0("../Data/RS_", db, "/SST"), pattern = "(summary).*\\.tif$", full.names = TRUE)

# Stack the files
sst_rasters <- rast(sst_files)

# Rename raster stack layers so it shows better in plot
names(sst_rasters) <- c("SST_max", "SST_mean", "SST_min", "SST_range", "SST_sd")
plot(sst_rasters)

# correlation plot
tiff(filename = "../Figures/sst_cor.tiff", width = 7, height = 7, units = "in", compression = "lzw", res = 300)
pairs(sst_rasters,  main = paste0("SST ", db))
dev.off()


# For CHLO
# Read summary rasters
chlo_files <- list.files(path = paste0("../Data/RS_", db, "/CHLO"), pattern = "(summary).*\\.tif$", full.names = TRUE)

# Stack the files
chlo_rasters <- rast(chlo_files)

# Rename raster stack layers so it shows better in plot
names(chlo_rasters) <- c("CHLO_max", "CHLO_mean", "CHLO_min", "CHLO_range", "CHLO_sd")

plot(chlo_rasters)

# correlation plot
tiff(filename = "../Figures/chlo_cor.tiff", width = 7, height = 7, units = "in", compression = "lzw", res = 300)
pairs(chlo_rasters,  main = paste0("CHLO ", db))
dev.off()


# Stack all 

all <- c(sst_rasters, chlo_rasters)
plot(all)

# correlation plot
tiff(filename = paste0("../Figures/corr_pair_plot_all.tiff"), width = 7, height = 7,
     units = "in", compression = "lzw", res = 300)
pairs(all) 
dev.off()




# OLD ---------------------------------------------------------------------


# Locate and move databases 

# Remotely Sensed Imagery

# if(F){
#   
#   # 1854  to 2022 
#   rs_dir <- "//idstorage.cnre.vt.edu/IDStorage1/Data/Copernicus_SeaSurface/1950s/OceanDB/finalData/Historical_Reconstructed_Surface_Temp/SST_1850-2022"
#   
#   sst_cropped <- lapply(list.files(path = rs_dir,
#                                    pattern = "sst", 
#                                    recursive = TRUE, 
#                                    full.names = TRUE), 
#                         raster)
#   
#   
#   # 1993 to 2022- Unstable
#   rs_dir <- "//idstorage.cnre.vt.edu/IDStorage1/Data/Copernicus_SeaSurface/1950s/OceanDB/finalData/ST+SSS_1993-2022/"
#   
#   sss_cropped <- lapply(list.files(path = paste0(rs_dir, "SSS_1993_Cropped"),
#                                    pattern = "sss", 
#                                    recursive = TRUE, 
#                                    full.names = TRUE), 
#                         raster)
#   
#   
#   sst_cropped <- lapply(list.files(path = paste0(rs_dir, "SST_1993_Cropped"),
#                                    pattern = "sst", 
#                                    recursive = TRUE, 
#                                    full.names = TRUE), 
#                         raster)
#   
#   
#   sss_cropped <- lapply(sss_cropped,FUN = function(ras)
#     resample(ras, sst_cropped[[1]]))
#   
#   
#   # stack(sss_cropped, sst_cropped)
#   
#   # 1958 to 2022 - Stable
#   rs_dir <- "//idstorage.cnre.vt.edu/IDStorage1/Data/Copernicus_SeaSurface/1958-2022/ST+SSS/"
#   
#   
#   
#   sss_cropped <- lapply(list.files(path = paste0(rs_dir, "Salinity/SSS_1958-2022_Cropped"),
#                                    pattern = "sss", 
#                                    recursive = TRUE, 
#                                    full.names = TRUE), 
#                         raster)
#   
#   
#   sst_cropped <- lapply(list.files(path = paste0(rs_dir, "Surface_Temp/SST_1958-2022_Cropped"),
#                                    pattern = "sst", 
#                                    recursive = TRUE, 
#                                    full.names = TRUE), 
#                         raster)
#   
#   # 2003 to 2022 - Stable - Already clean
# }

# Clean data 
# if(T){
#   
#   if(db == "(1850-2022)"){
#     years <- 1854:2022
#     var <- c("SST")
#     
#   }else if(db == "(1958-2022)"){
#     years <- 1958:2022
#     var <- c("SST", "SSS")
#     
#   }else if(db == "(1993-2022)"){
#     years <- 1993:2022
#     var <- c("SST", "SSS")
#   }
#   
#   if(!dir.exists(paste0("../Data/RS_", db, "/"))) dir.create(paste0("../Data/RS_", db, "/"))
#   
#   for(i in var){
#     
#     if(db == "(1850-2022)"){
#       og_var_dir <- "//idstorage.cnre.vt.edu/IDStorage1/Data/Copernicus_SeaSurface/1950s/OceanDB/finalData/Historical_Reconstructed_Surface_Temp/SST_1850-2022/"
#     }else if(db == "(1958-2022)"){
#       if(i == "SST"){
#         og_var_dir <- "//idstorage.cnre.vt.edu/IDStorage1/Data/Copernicus_SeaSurface/1958-2022/ST+SSS/Surface_Temp/SST_1958-2022_Cropped/"
#       }else{
#         og_var_dir <- "//idstorage.cnre.vt.edu/IDStorage1/Data/Copernicus_SeaSurface/1958-2022/ST+SSS/Salinity/SSS_1958-2022_Cropped/"
#       }
#     }
#     
#     new_var_dir <- paste0("../Data/RS_", db, "/", i)
#     if(!dir.exists(new_var_dir)) dir.create(new_var_dir)
#     
#     for(y in years){
#       og_var_year_dir <- paste0(og_var_dir, y)
#       
#       og_var_year_ras <- lapply(list.files(path = og_var_year_dir, 
#                                            pattern="*.tif",                                         
#                                            recursive = TRUE,
#                                            full.names = TRUE), 
#                                 raster)
#       
#       new_var_year_dir <- paste0(new_var_dir, "/", y)
#       if(!dir.exists(new_var_year_dir)) dir.create(new_var_year_dir)
#       
#       for(ras in og_var_year_ras){
#         ras_cropped <- crop(ras, study_area_extent)
#         
#         # if(db == "(1993-2022)" && i == "SSS"){
#         #   ras_cropped <- resample(ras_cropped, crop(sst_cropped[[1]], study_area_extent))
#         # }
#         
#         ras_name <-  basename(ras@file@name)
#         
#         writeRaster(x = ras_cropped, 
#                     filename = paste0(new_var_year_dir, "/", ras_name), 
#                     format="GTiff", 
#                     overwrite=TRUE)
#         
#       }
#     }
#   }
# }
# 
# 




