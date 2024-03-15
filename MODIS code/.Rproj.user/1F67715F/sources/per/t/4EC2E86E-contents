# Title: Remote Sensing and Cholera Code Supporting Function
# Authors: Mariana Castaneda-Guzman
# Date last updated: 06/17/2021

if(!require("terra")) install.packages("terra")
if(!require("sf")) install.packages("sf")
if(!require("parallel")) install.packages("parallel")
if(!require("RNetCDF")) install.packages("RNetCDF")

num_cores <- detectCores()
clust <- makeCluster(num_cores-2)

# Function to create sequence of specially formed URLS 
create_sequences <- function(name, format_1, format_2, years, months = 1:12){
  
  years_count <- length(years)
  months_count <- length(months)
  
  months <- sprintf("%.2d", months)
  sequence <- 1:(years_count*months_count)
  
  n = 1
  
  for(i in 1:years_count){
    for(j in 1:months_count){
      sequence[n] <- paste0("[(", years[i],"-", months[j],"-16):1:(", years[i],"-", months[j], "-16)]")
      n = n + 1
    }
  }
  
  
  full_sequence <- paste0(format_1, sequence, format_2)
  
  dates_1 <- data.frame(date = seq(as.Date(paste(years[1], months[1], "01", sep = "-")),
                                   as.Date(paste(years[years_count], months[months_count], "01", sep = "-")),
                                   by = "month"))
  
  # Check dates are in the subset of years
  dates <- subset(dates_1, format(dates_1$date, "%m") %in% months)
  dates <- subset(dates, format(dates$date, "%Y") %in% years)
  
  save <- paste0(name, "-", dates$date,'.nc')
  
  sequence <- list(full_sequence, save)
  
  return(sequence)
}

# Return the corresponding nc files (mostly useful in debug mode)
get_NetCDF_files <- function(folder, sequence){
  
  netcdf_files <- list.files(path = folder, 
                       pattern = "*.nc", 
                       full.names = TRUE,
                       recursive = TRUE)
  
  my_files <- netcdf_files[which(basename(netcdf_files) %in% sequence)]
  
  return(my_files)
}

# Writes rasters from the NetCDF files
get_raster <- function(name, dir, netCDF, min_lat, max_lat, max_lon, min_lon){
  
  raster_name <- gsub("\\..*", "", basename(netCDF))
  file_name <- paste0(dir, raster_name, ".tif")
  
  if(file.exists(file_name)){
    tryCatch({ 
      
      terra::rast(file_name)

    }, error = function(e){
      
      image_ncdf <- RNetCDF::read.nc(RNetCDF::open.nc(netCDF))
      
      
      if(name == "CHLO"){
        raster_extent <- terra::ext(min_lon, max_lon, min_lat, max_lat)
        data <- terra::rast(t(image_ncdf$chlor_a))
        
      }else if (name == "SST"){
        raster_extent <- terra::ext(min_lon, max_lon, min_lat, max_lat)
        data <- terra::rast(t(image_ncdf$sst))
      }
      
      terra::ext(data) <- raster_extent
      
      terra::crs(data) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      names(data) <- raster_name
      
      terra::writeRaster(x = data, 
                          filename = file_name, 
                          overwrite = TRUE)
      
    })
  }else{
    image_ncdf <- RNetCDF::read.nc(RNetCDF::open.nc(netCDF))
    
    if(name == "CHLO"){
      raster_extent <- terra::ext(min_lon, max_lon, min_lat, max_lat)
      # The column for chlo-a changes in the datasets, common error. change
      # here.
      # data <- terra::rast(t(image_ncdf$chlor_a))
      
      data <- terra::rast(t(image_ncdf$chlorophyll))
      
    }else if(name == "SST"){
      raster_extent <- terra::ext(min_lon, max_lon, min_lat, max_lat)
      data <- terra::rast(t(image_ncdf$sst))
    }
    
    terra::ext(data) <- raster_extent
    
    terra::crs(data) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    
    names(data) <- raster_name
    
    raster::writeRaster(x = data, 
                        filename = file_name,
                        overwrite = TRUE)
    
  }
  
  return(terra::rast(file_name))
  
}

# Crop raster to the correct extent, EEZ.
crop_raster <- function(dir, borders, raster){
  
  file_name <- paste0(dir, names(raster), ".tif")
  
  if(file.exists(file_name)){
    tryCatch({
      
      new_raster <- terra::rast(file_name)
      
    }, error = function(e){
      
      new_raster <- terra::mask(terra::crop(raster, borders), borders)
      
      terra::writeRaster(x = new_raster, 
                          filename = file_name, 
                          overwrite = TRUE)
    }
    )
  }else{
    
    new_raster <- terra::mask(terra::crop(raster, borders), borders)
    
    
    terra::writeRaster(x = new_raster, 
                       filename = file_name, 
                       overwrite = TRUE)
    
  }
  
  return(terra::rast(file_name))
  
}

# Group rasters
get_stacks <- function(name, pattern, dir, dest_dir, month = FALSE){
  
  if(month){
    
    month_pattern <- paste0(pattern, "[-\\.]01\\.tif$")
    
    
    files <-  list.files(path = dir, 
                         pattern = month_pattern,
                        full.names = TRUE, 
                        recursive = FALSE)
    
  }else{
    files <- list.files(path = dir, 
                        pattern = paste0("", pattern, ""),
                        full.names = TRUE, 
                        recursive = FALSE)
  }
  
  file_name <- paste0(dest_dir, name, "_stack_", pattern, ".tif")
  
  if(file.exists(file_name)){
    tryCatch({
      
      new_raster <- terra::rast(file_name)
      
    }, error = function(e){
      
      interval_raster_stack <- terra::rast(files)
      
      terra::writeRaster(x = interval_raster_stack, 
                          filename = file_name,
                          overwrite = TRUE)
    }
    )
  }else{
    
    interval_raster_stack <- terra::rast(files)
    
    terra::writeRaster(x = interval_raster_stack, 
                        filename = file_name,
                        overwrite = TRUE)
    
    
  }
  
  return(paste0(dest_dir, name, "_stack_", pattern, ".tif"))
}

# Function to calculate monthly and yearly statistics
calculate_statistics <- function(name, dir, file){
  
  rasters_stack <- rast(file)
  
  stats_list <- c(mean, max, min, range, sd)
  stats_list_name <- c("mean", "max", "min", "range", "sd")
  
  results <- lapply(1:length(stats_list), function(i){
    
    file_name <- paste0(name, "_", gsub("\\D", "", basename(file)), "_", stats_list_name[i], ".tif")
    
    if(file.exists(paste0(dir, file_name))){
      tryCatch({
        new_raster <- terra::rast(paste0(dir, file_name))
      }, error = function(e){
        if(stats_list_name[i] == "range"){
          stat_calc_min <- terra::app(rasters_stack, fun = stats_list[[3]], na.rm = TRUE)
          stat_calc_max <- terra::app(rasters_stack, fun = stats_list[[2]], na.rm = TRUE)
          
          stat_calc <- stat_calc_max - stat_calc_min
          names(stat_calc) <- "range"
        }else{
          stat_calc <- terra::app(rasters_stack, fun = stats_list[[i]], na.rm = TRUE)
        }
        
        raster::writeRaster(x = stat_calc,
                            filename = paste0(dir, file_name),
                            overwrite = TRUE)
      })
    }else{
      if(stats_list_name[i] == "range"){
        stat_calc_min <- terra::app(rasters_stack, fun = stats_list[[3]], na.rm = TRUE)
        stat_calc_max <- terra::app(rasters_stack, fun = stats_list[[2]], na.rm = TRUE)
        
        stat_calc <- stat_calc_max - stat_calc_min
        
      }else{
        stat_calc <- terra::app(rasters_stack, fun = stats_list[[i]], na.rm = TRUE)
      }
      
      raster::writeRaster(x = stat_calc,
                          filename = paste0(dir, file_name),
                          overwrite = TRUE)
    }
  })
}


