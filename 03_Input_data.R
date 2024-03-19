# Title: Clean occurrences with Remotely Sensed Imagery
# Author: Mariana Castaneda-Guzman
# Date created: 9/12/2023
# Date last updated: 3/18/2023

# Description: Cleaning occurrences one per pixel per year to reduce sampling
# bias, or inflating the models

# Packages ----------------------------------------------------------------

if(!require(terra)) install.packages("terra")
if(!require(sf)) install.packages("sf")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")
if(!require(data.table)) install.packages("data.table")


# Specify Database --------------------------------------------------------

# Database
db <- "(2003-2023)"


# Create mask -------------------------------------------------------------
if(!file.exists(paste0("../Data/RS_", db, "/mask_", db, ".tif"))){
  # Select a random raster, could be any raster downloaded in the previous
  # scripts. This section will generate the calibration area
  mask <- rast(list.files(paste0("../Data/RS_", db), 
                          pattern = "*.tif", 
                          recursive = TRUE, 
                          full.names = TRUE)[1])
  
  # set value lengths to those of your mask (ensure no "hang over" of values)
  all_length <- length(values(mask)[!is.na(values(mask))])
  
  # Ignore value and assign to the a unique identifier (aka row or pixel index)
  values(mask)[!is.na(values(mask))] <- c(1:all_length)
  
  # overwrite your mask with the correct mask where value lengths match 
  writeRaster(mask, paste0("../Data/RS_", db, "/mask_" , db, ".tif"), overwrite = TRUE)
}else{
  
  #Read in your previously calibrated mask 
  mask <- rast(paste0("../Data/RS_", db, "/mask_" , db, ".tif"))
}

# With new filtering function ---------------------------------------------

# Occurrences
og_occ <- readRDS("../Data/GBIF/cv_occurrences_clean.RDS")
nrow(og_occ)

source("../../00 Global Code/move_points.R")

message(paste0("In database ", db))

# Create spatial points, some of which may be outside the extent
og_occ_sp <- st_as_sf(og_occ, 
                      coords = c("decimalLongitude", 
                                 "decimalLatitude"), 
                      crs = st_crs(mask))

if(!dir.exists(paste0("../Data/GBIF/RS_", db))) dir.create(paste0("../Data/GBIF/RS_", db))

if(!file.exists(paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_", db, ".RDS"))){
  og_occ_sp_corrected <- fix_to_extent(mask, og_occ_sp, spatial = TRUE)
  og_occ_sp_final <- find_nearest_non_na(r = mask, og_occ_sp_corrected, 
                                         spatial = TRUE)
  
  og_occ_sp_final_df <- as.data.frame(og_occ_sp_final)
  og_occ_sp_final_df$year <- og_occ$year
  
  saveRDS(og_occ_sp_final_df, paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_", db, ".RDS"))
  
  
  tiff(filename = paste0("../Figures/RS_", db, "/occ_corrected.tiff"), 
       width = 5, height = 6, units = "in", compression = "lzw", res = 300)
  
  pts <- st_coordinates(og_occ_sp)
  new.pts.df <- st_coordinates(og_occ_sp_final)
  nearest_values_df <- st_coordinates(og_occ_sp_corrected)
  
  mask_df <- as.data.frame(mask, xy = TRUE)
  
  lines_nearest <- data.frame(x = pts[, 1], 
                              y = pts[, 2],
                              xend = nearest_values_df[, 1],
                              yend = nearest_values_df[, 2])
  
  lines_new_pts <- data.frame(x = nearest_values_df[, 1], 
                              y = nearest_values_df[, 2],
                              xend = new.pts.df[, 1], 
                              yend = new.pts.df[, 2])
  
  
  ggplot() +
    geom_raster(data = mask_df , aes(x = x, y = y), fill = "lightblue") +
    geom_polygon(data = states %>% filter(!is.na(east_coast)), 
                 aes(x = long, y = lat, group = group), fill = "lightgrey",
                 color = "black") +
    geom_sf(data = og_occ_sp, shape = 21, col = "black")+
    geom_sf(data = og_occ_sp_corrected, shape = 1, col = "red") +
    geom_sf(data = og_occ_sp_final, shape = 17, col = "green") +
    xlab("Longitude") + ylab("Latitude") +
    geom_segment(data = lines_nearest,
                 mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 col='red', lwd = 1) +
    geom_segment(data = lines_new_pts, 
                 mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 col='green', lwd = 1) 
  
  # coord_fixed(1.3)
  
  dev.off()
  
  
  # Show only the pixels that had a difference
  diff <- data.frame(x_diff = pts[, 1] == new.pts.df[, 1],
                     y_diff = pts[, 2] == new.pts.df[, 2])
  
  diff$diff <- diff$x_diff + diff$y_diff
  
  count <- nrow(subset(diff, diff < 2))
  
  points_diff <- pts[diff$diff < 2, ]
  new.pts_diff <- new.pts.df[diff$diff < 2, ]
  
  tiff(filename = paste0("../Figures/RS_", db, "/occ_diff.tiff"), 
       width = 5, height = 6, units = "in", compression = "lzw", res = 300)
  # plot(mask, colNA = "lightgrey")
  # plot(as.polygons(mask, na.rm = FALSE, aggregate=FALSE), add=TRUE, border='black', lwd=1)
  # points(points_diff, pch = 21, col = "black")
  # points(new.pts_diff, pch = 17, col = "green")
  
  pts <- as.data.frame(points_diff)
  new.pts.df <- as.data.frame(new.pts_diff)
  
  mask_df <- as.data.frame(mask, xy = TRUE)
  
  lines_diff_new <- data.frame(x = pts[, 1], 
                              y = pts[, 2],
                              xend = new.pts.df[, 1], 
                              yend = new.pts.df[, 2])
  
  pts <- st_as_sf(pts, coords = c("X", "Y"))
  new.pts.df <- st_as_sf(new.pts.df, coords = c("X", "Y"))
  
  ggplot() +
    geom_raster(data = mask_df , aes(x = x, y = y), fill = "lightblue") +
    geom_polygon(data = states %>% filter(!is.na(east_coast)), 
                 aes(x = long, y = lat, group = group), fill = "lightgrey",
                 color = "black") +
    geom_sf(data = pts, shape = 21, col = "black") +
    geom_sf(data = new.pts.df, shape = 17, col = "green") +
    xlab("Longitude") + ylab("Latitude") +
    geom_segment(data = lines_diff_new,
                 mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 col='red', lwd = 1) +
    coord_sf()
  
  
  dev.off()
  
}

message(paste0("Done with database ", db))
  

# Filter occurrence points to one occurrence per pixel --------------------

message(paste0("In database ", db))

# convert the mask to points 
mask_p <- terra::as.data.frame(mask, xy = TRUE)

# use these points to create pixel index. This assigns each pixel it's own
# unique identifier. Same thing as to assign the row number as identifier
colnames(mask_p)[3] <- "pixel_index"

#read in occurrence points from previous code 
og_occ_c <- readRDS(paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_", db, ".RDS"))

# assign pixel index to occurrence data frame 
og_occ_c$pixel_index <- terra::extract(mask, 
                                       st_coordinates(st_as_sf(og_occ_c$geometry)))[,1]

og_occ_c <- og_occ_c[which(!is.na(og_occ_c$pixel_index)), ]

table(og_occ_c$year)
sum(table(og_occ_c$year))

# Clean to one per pixel
length(unique(og_occ_c$pixel_index)) # 871 

og_occ_unique <- unique(og_occ_c[, c("year", "pixel_index")])
nrow(og_occ_unique) # unique 2214 

table(og_occ_unique$year)
sum(table(og_occ_unique$year))

# remove NAs 
og_occ_unique <- og_occ_unique[!is.na(og_occ_unique$pixel_index),]

# Join the final occurrence set with the mask points, to get the latitude and
# longitude of each of the pixels
og_occ_unique <- left_join(og_occ_unique, mask_p, by="pixel_index")
nrow(og_occ_unique)

# Plot cleaned occurences for each db
tiff(filename = paste0("../Figures/RS_", db, "/occ_1p.tiff"), 
     width = 5, height = 6, units = "in", compression = "lzw", res = 300)

plot(mask, colNA = "lightgrey")
plot(st_as_sf(og_occ_unique, coords = c("x", "y")), pch = 21, col = "black", add = T)
title("One point per pixel, per year")
dev.off()

# save one per pixel occurrences as R object 
saveRDS(og_occ_unique, paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_one_per_pixel_", db, ".RDS"))


# Aggregate occurrence points per pixel --------------------------------

message(paste0("In database ", db))

# convert the mask to points 
mask_p <- terra::as.data.frame(mask, xy = TRUE)

# use these points to create pixel index. This assigns each pixel it's own
# unique identifier. Same thing as to assign the row number as identifier
colnames(mask_p)[3] <- "pixel_index"

#read in occurrence points from previous code 
og_occ_c <- readRDS(paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_", db, ".RDS"))

# assign pixel index to occurrence data frame 
og_occ_c$pixel_index <- terra::extract(mask, 
                                       st_coordinates(st_as_sf(og_occ_c$geometry)))[,1]

og_occ_c <- og_occ_c[which(!is.na(og_occ_c$pixel_index)), ]

table(og_occ_c$year)
sum(table(og_occ_c$year))


og_occ_agg <- og_occ_c %>% 
  group_by(pixel_index, year) %>% 
  mutate(obs = 1) %>% 
  summarise(obs = sum(obs))

og_occ_unique_agg <- left_join(og_occ_unique, og_occ_agg, by = c("pixel_index", "year"))

saveRDS(og_occ_unique_agg, paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_one_per_pixel_weighted_", db, ".RDS"))

# sum(og_occ_unique_agg$obs)

# # Buffers of Uncertainty --------------------------------------------------
# 
# # Occurrences
# og_occ <- readRDS("../Data/GBIF/cv_occurrences_clean.RDS")
# nrow(og_occ)
# 
# # Add and id column, this will be use to create the random points
# og_occ$id <-  1:nrow(og_occ)
# 
# 
# # Select columns to make a spatial object
# og_occ_sub <- og_occ[ , c("decimalLatitude","decimalLongitude")]
# 
# 
# # Force the file to be spatial
# # coordinates = Makes a new object spatial by making two values coordinates
# coordinates(og_occ_sub) <- ~ decimalLongitude + decimalLatitude
# length(og_occ_sub)
# 
# # Check the spatial og_occ
# # coords=coordinates
# head(og_occ_sub@coords)
# 
# # Check the Coordinates Reference System
# # non-capital fond crs=check CRS
# crs(og_occ_sub)
# 
# # Assign CRS
# myCRS1 <- CRS("+init=epsg:4326") # WGS 84, units=degrees, WGS=WGS84, ellipsoid=WGS84
# crs(og_occ_sub) <- myCRS1
# crs(og_occ_sub)
# 
# myCRS2 <- CRS("+init=epsg:3857") # Mercator, units=meters, WGS=WGS84, ellipsoid=WGS84
# og_occ_sub_projected <- spTransform(og_occ_sub, myCRS2)
# 
# # Create buffer zone:
# buffer_size <- ifelse(db == "(1958-2022)", 100000, 
#                       ifelse(db == "(1993-2022)", 10000, 10000))
# 
# og_occ_sub_buffer <- buffer(og_occ_sub_projected, width = buffer_size, dissolve = FALSE)
# 
# # Transform buffer back to WGS 84
# og_occ_sub_buffer <- spTransform(og_occ_sub_buffer, myCRS1)
# 
# 
# 
# # Create extended Occurrences -----------------------------------------------------
# 
# # Loop will iterate through each year and subset of polygons in that given year.
# # There is one polygon per line in the df
# 
# # new df to add all the new values for each point
# extended_df <- list()
# 
# for(y in min(og_occ$year, na.rm = T):max(og_occ$year, na.rm = T)){
#   
#   # The function which return a list of indexes of those observations with the
#   # same year. index = row number
#   polygon_index <- which(og_occ$year == y)
#   
#   if(length(polygon_index) > 0){
#     # This inner loop will iterate through the buffers, and create and extract the
#     # random point
#     extended_item <- list()
#     
#     for(i in polygon_index){
#       
#       # Create random point for unique buffer
#       random_points <- spsample(og_occ_sub_buffer@polygons[[i]], 100, type = "random")
#       crs(random_points) <- myCRS1
#       
#       # New df, this will only have the information of the points for this loop
#       # iteration.
#       new_rows <- as.data.frame(random_points)
#       
#       # Add the corresponding id to the random points. Id corresponds to the
#       # polygon number, that is the same as the row number in the original df
#       new_rows$id <- i
#       extended_item[[as.character(i)]] <- new_rows
#     }
#     
#     extended_item <- rbindlist(extended_item)
#     extended_item <- data.frame(extended_item)
#     
#     # Rename columns
#     names(extended_item) <- c("decimalLongitude", "decimalLatitude", "id")
#     
#     # Reorder columns (optional)
#     extended_item <- extended_item[ , c("id", "decimalLongitude", "decimalLatitude")]
#     
#     # Add new rows to extended df
#     extended_df[[as.character(y)]] <- extended_item
#   }
# }
# 
# extended_occ <- rbindlist(extended_df)
# 
# extended_occ_1 <- left_join(extended_occ, og_occ[ , c("id", "year")], by = "id")
# 
# extended_occ_1 <- data.frame(extended_occ_1)
# 
# saveRDS(extended_occ_1, paste0("../Data/GBIF/cv_occurrences_clean_extended_", db, ".RDS"))
# 
# 
# 
# # Filter occurrence points to one occurrence per pixel --------------------
# 
# # convert the mask to points 
# mask_p <- data.frame(rastToPoints(mask))
# 
# # use these points to create pixel index. This assigns each pixel it's own
# # unique identifier. Same thing as to assign the row number as identifier
# colnames(mask_p)[3] <- "pixel_index"
# 
# #read in occurrence points from previous code 
# og_occ_c <- readRDS(paste0("../Data/GBIF/cv_occurrences_clean_extended_", db, ".RDS"))
# 
# # assign pixel index to occurrence data frame 
# og_occ_c$pixel_index <- rast::extract(mask, og_occ_c[, c("decimalLongitude", "decimalLatitude")])
# 
# og_occ_c <- og_occ_c[which(!is.na(og_occ_c$pixel_index)), ]
# 
# # See how many observations per pixel there are
# sum(table(og_occ_c$pixel_index))
# 
# # Identify unique occurrences by pixel and year 
# og_occ_unique <- unique(og_occ_c[, c("year", "pixel_index")])
# nrow(og_occ_unique) # 5990 original occ, 444 unique occ
# 
# table(og_occ_unique$year[og_occ_unique$year > 1949])
# sum(table(og_occ_unique$year[og_occ_unique$year > 1949])) # 3711 unique observations since 1950
# 
# 
# # remove NAs 
# og_occ_unique <- og_occ_unique[!is.na(og_occ_unique$pixel_index),]
# 
# # Join the final occurrence set with the mask points, to get the latitude and
# # longitude of each of the pixels
# og_occ_unique <- left_join(og_occ_unique, mask_p, by="pixel_index")
# 
# # save one per pixel occurrences as R object 
# saveRDS(og_occ_unique, paste0("../Data/GBIF/cv_occurrences_clean_one_per_pixel_", db, ".RDS"))
# 
# 
