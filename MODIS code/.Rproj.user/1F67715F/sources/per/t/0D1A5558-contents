# Title: A Database of Global Coastal Conditions R Code
# Authors: Castaneda-Guzman, Mariana; Mantilla-Santos, Gabriel; Escobar, Luis E. 
# Date last updated: 3/7/2024

years <- 2023

months <- 1:12

setwd("../../Data/MODIS")

# Now, the code would automatically download both SST and CHLO. If you would
# like to download only information for an individual observation change the
# value after the assign operator '<-' to "SST", "CHLO". 
models <- "both"

# Specify below the coordinates that you wish to download. If you wish to
# download global data,leave next few lines as is. Note, this will take some
# time to run. Look at map for coordinates: https://www.satsig.net/lat_long.htm


# Latitude values should be in degrees north (current values = defaults for world)
max_lat <- 47.459574
min_lat <- 23.140514

# Longitude values should be in degrees east (current values = defaults for world)
min_lon <- -87.642057
max_lon <- -66.950473

# study_area_extent <- c(-87.642057, -66.950473, 23.140514, 47.459574)

source('../../Code/MODIS code/supporting_functions.r')

# (a) data procurement ----------------------------------------------------

# This section of code will create the downloading strings for both CHLO and SST
# but would only download the specific string specified above.

# Formats for SST, ERDP server
# erdMH1sstdmday_R2022SQMasked.nc?sstMask
# erdMH1sstdmdayR20190SQ.nc?sstMasked

format_SST <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1sstdmday_R2022SQMasked.nc?sstMask"
format_SST_2 <- paste0("[(", max_lat, "):1:(", min_lat, ")][(", min_lon, "):1:(", max_lon, ")]")

# Get sequences for download. function create_sequence will return a list with
# two lists. [[1]] the downloading string, [[2]] the name to give to the file
SST_sequence <- create_sequences(name = "SST", format_1 = format_SST, format_2 = format_SST_2, 
                                 years = years, months = months)


# SST_sequence <- SST_sequence[-c(226:228)]
# Formats for CHLO, ERDP server
# erdMH1chlamday.nc?chlorophyll
# erdMH1chlamday_R2022SQ.nc?chlor_a
# erdMH1chlamday_R2022NRT.nc?chlorophyll 

format_CHLO <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022NRT.nc?chlorophyll"
format_CHLO_2 <- paste0("[(", max_lat, "):1:(", min_lat, ")][(", min_lon, "):1:(", max_lon, ")]")

# Get sequences for download. function create_sequence will return a list with
# two lists. [[1]] the downloading string, [[2]] the name to give to the file
CHLO_sequence <- create_sequences(name = "CHLO", format_1 = format_CHLO, format_2 = format_CHLO_2, 
                                  years = years, months = months)


# The next few lines are to download NetCDF files from NOAA ERDAP Server. The
# file will download if and only if it does not exit in the given folder

# If models is only CHLO this if statement wont run, if is either both or SST
# this will run. 
if(models != "CHLO"){
  if(!dir.exists("SST")) dir.create("SST")
  if(!dir.exists("SST/SST_NetCDF")) dir.create("SST/SST_NetCDF")
  
  parLapply(cl = clust, 
            X = 1:length(SST_sequence[[1]]), 
            fun = function(i, SST_sequence){
              file_name = paste0("SST/SST_NetCDF/", SST_sequence[[2]][i])
              sequence = SST_sequence[[1]][i]
              
              if(!file.exists(file_name)){
                download.file(sequence, destfile = file_name, metohd = "auto", mode = "wb", quiet = TRUE)
              }
              
            }, SST_sequence)
}

# If models is only CHLO this if statement wont run, if is either both or SST
# this will run.

if(models != "SST"){
  if(!dir.exists("CHLO")) dir.create("CHLO")
  if(!dir.exists("CHLO/CHLO_NetCDF")) dir.create("CHLO/CHLO_NetCDF")
  
  parLapply(cl = clust, 
            X = 1:length(CHLO_sequence[[1]]), 
            fun = function(i, CHLO_sequence){
              file_name = paste0("CHLO/CHLO_NetCDF/", CHLO_sequence[[2]][i])
              sequence = CHLO_sequence[[1]][i]
              
              if(!file.exists(file_name)){
                download.file(sequence, destfile = file_name, metohd = "auto", mode = "wb", quiet = TRUE)
              }
              
            }, CHLO_sequence)
  
}



# (b) preparation ---------------------------------------------------------

if(models != "CHLO"){
  if(!dir.exists("SST/SST_raw_rasters")) dir.create("SST/SST_raw_rasters")
  
  parLapply(cl = clust, 
            X = 1:length(SST_sequence[[1]]),
            fun = function(i, SST_sequence, get_NetCDF_files, get_raster, min_lat, max_lat, max_lon, min_lon){
              get_raster(name = "SST",
                         dir = "SST/SST_raw_rasters/",
                         netCDF = get_NetCDF_files(folder = "SST/SST_NetCDF",sequence = SST_sequence[[2]][i]),
                         min_lat = min_lat, max_lat = max_lat, max_lon = max_lon, min_lon = min_lon)
            }, SST_sequence, get_NetCDF_files, get_raster, min_lat, max_lat, max_lon, min_lon)
}

if(models != "SST"){
  if(!dir.exists("CHLO/CHLO_raw_rasters")) dir.create("CHLO/CHLO_raw_rasters")
  
  # If this is throwing an error make sure the correct information if being red
  # in the netdf immage. look at get_raster second if statement.
  parLapply(cl = clust, 
            X = 1:length(CHLO_sequence[[1]]),
            fun = function(i, CHLO_sequence, get_NetCDF_files, get_raster, min_lat, max_lat, max_lon, min_lon){
              get_raster(name = "CHLO",
                         dir = "CHLO/CHLO_raw_rasters/",
                         netCDF = get_NetCDF_files(folder = "CHLO/CHLO_NetCDF", sequence = CHLO_sequence[[2]][i]),
                         min_lat = min_lat, max_lat = max_lat, max_lon = max_lon, min_lon = min_lon)
            }, CHLO_sequence, get_NetCDF_files, get_raster, min_lat, max_lat, max_lon, min_lon)
  
  
  
}


# (c) processing ----------------------------------------------------------

# Get Exclusive Economic Zones Shapefile (EZZ). The code below loads the file
# EZZ_borders.RDS, if it throws an error 'No such file or directory' is becuase
# file is not in your current directory make sure the file is in your directory,
# or specify the file path.

borders <- readRDS("../../Code/MODIS code/EEZ_borders.RDS")
borders <- sf::st_as_sf(borders)

borders_extent <- terra::ext(min_lon, max_lon, min_lat, max_lat)
borders <- st_crop(borders, borders_extent)

# This part might take some time, depending on the extent of your study area
if(models != "CHLO"){
  if(!dir.exists("SST/SST_cropped_rasters")) dir.create("SST/SST_cropped_rasters")
  
  
  if(!exists("SST_raw_rasters")){
    SST_raw_rasters <- lapply(list.files("SST/SST_raw_rasters/", 
                                         pattern = "2023",
                                         full.names = T),
                              rast)
    
  }else if(length(SST_raw_rasters) == 0){
    SST_raw_rasters <- lapply(list.files("SST/SST_raw_rasters/", 
                                         pattern = "2023",
                                         full.names = T),
                              rast)
  }
  
  
  # parLapply(cl = clust, SST_raw_rasters,
  #           fun = function(i, crop_raster, borders){
  #             crop_raster(dir = "SST/SST_cropped_rasters/",
  #                         borders = borders,
  #                         raster = i)
  #           }, crop_raster, borders)

  
  
  for(i in SST_raw_rasters){
    crop_raster(dir = "SST/SST_cropped_rasters/",
                borders = borders,
                raster = i)
  }


}

if(models != "SST"){
  if(!dir.exists("CHLO/CHLO_cropped_rasters")) dir.create("CHLO/CHLO_cropped_rasters")
  
  if(!exists("CHLO_raw_rasters")){
    CHLO_raw_rasters <- lapply(list.files("CHLO/CHLO_raw_rasters/", 
                                         pattern = "2023",
                                         full.names = T),
                              rast)
    
  }else if(length(CHLO_raw_rasters) == 0){
    CHLO_raw_rasters <- lapply(list.files("CHLO/CHLO_raw_rasters/", 
                                         pattern = "2023",
                                         full.names = T),
                              rast)
  }
  
  # parLapply(cl = clust, CHLO_raw_rasters, 
  #           fun = function(i, crop_raster, borders){
  #             crop_raster(dir = "CHLO/CHLO_cropped_rasters/", 
  #                         borders = borders, 
  #                         raster = i)
  #           }, crop_raster, borders)
  
  for(i in CHLO_raw_rasters){
    crop_raster(dir = "CHLO/CHLO_cropped_rasters/", 
                borders = borders, 
                raster = i)
  }
  
  
}


# (d) analysis ------------------------------------------------------------

years <- 2003:2023

# Create Stacks
if(models != "CHLO"){
  if(!dir.exists("SST/SST_statistical_results/")) dir.create("SST/SST_statistical_results/")
  if(!dir.exists("SST/SST_statistical_results/SST_yearly_statistics")) dir.create("SST/SST_statistical_results/SST_yearly_statistics")
  
  
  parLapply(cl = clust, X = years, fun = function(i, get_stacks){
    get_stacks(name = "SST",
               pattern = i, 
               dir = "SST/SST_cropped_rasters/",
               dest_dir = "SST/SST_statistical_results/SST_yearly_statistics/")
    
  }, get_stacks)

  
  if(!dir.exists("SST/SST_statistical_results/SST_monthly_statistics")) dir.create("SST/SST_statistical_results/SST_monthly_statistics")
  
  # Create monthly stacks
  months_pattern <- sprintf("%02d", months)

  parLapply(cl = clust, X = months_pattern, fun = function(i, get_stacks){
    get_stacks(name = "SST",
               pattern = i, 
               dir = "SST/SST_cropped_rasters/",
               dest_dir = "SST/SST_statistical_results/SST_monthly_statistics/",
               month = TRUE)
    
  }, get_stacks)

}

if(models != "SST"){
  if(!dir.exists("CHLO/CHLO_statistical_results/")) dir.create("CHLO/CHLO_statistical_results/")
  if(!dir.exists("CHLO/CHLO_statistical_results/CHLO_yearly_statistics")) dir.create("CHLO/CHLO_statistical_results/CHLO_yearly_statistics")
  
  
  parLapply(cl = clust, X = years, fun = function(i, get_stacks){
    get_stacks(name = "CHLO",
               pattern = i, 
               dir = "CHLO/CHLO_cropped_rasters/",
               dest_dir = "CHLO/CHLO_statistical_results/CHLO_yearly_statistics/")
    
  }, get_stacks)
  
  
  if(!dir.exists("CHLO/CHLO_statistical_results/CHLO_monthly_statistics")) dir.create("CHLO/CHLO_statistical_results/CHLO_monthly_statistics")
  
  # Create monthly stacks
  months_pattern <- sprintf("%02d", months)
  
  parLapply(cl = clust, X = months_pattern, fun = function(i, get_stacks){
    get_stacks(name = "CHLO",
               pattern = i, 
               dir = "CHLO/CHLO_cropped_rasters/",
               dest_dir = "CHLO/CHLO_statistical_results/CHLO_monthly_statistics/",
               month = TRUE)
    
  }, get_stacks)
}


# Calculate Statistics
if(models != "CHLO"){
  
  # Read in raster stacks
  list_files <- list.files(path = "SST/SST_statistical_results/SST_yearly_statistics/", 
                           pattern = "(stack).*\\.tif$",
                           full.names = TRUE,
                           recursive = FALSE)
  
  
  # Calculate statistics
  yearly_results <- lapply(list_files, FUN = function(file){
    calculate_statistics(name = "SST", 
                         dir = "SST/SST_statistical_results/SST_yearly_statistics/",
                         file = file)
  })
  
  
  # Read in raster stacks
  list_files <- list.files(path = "SST/SST_statistical_results/SST_monthly_statistics/", 
                           pattern = "(stack).*\\.tif$",
                           full.names = TRUE,
                           recursive = FALSE)
  
  
  # Calculate statistics
  monthly_results <- lapply(list_files, FUN = function(file){
    calculate_statistics(name = "SST", dir = "SST/SST_statistical_results/SST_monthly_statistics/", file = file)
  })
  
}

if(models != "SST"){
  
  # Read in raster stacks
  list_files <- list.files(path = "CHLO/CHLO_statistical_results/CHLO_yearly_statistics/", 
                           pattern = "(stack).*\\.tif$",
                           full.names = TRUE,
                           recursive = FALSE)
  
  
  # Calculate statistics
  yearly_results <- lapply(list_files, FUN = function(file){
    calculate_statistics(name = "CHLO", dir = "CHLO/CHLO_statistical_results/CHLO_yearly_statistics/", file = file)
  })
  
  
  # Read in raster stacks
  list_files <- list.files(path = "CHLO/CHLO_statistical_results/CHLO_monthly_statistics/", 
                           pattern = "(stack).*\\.tif$",
                           full.names = TRUE,
                           recursive = FALSE)
  
  
  # Calculate statistics
  monthly_results <- lapply(list_files, FUN = function(file){
    calculate_statistics(name = "CHLO", dir = "CHLO/CHLO_statistical_results/CHLO_monthly_statistics/", file = file)
  })
  
}

