# Title: ENMs parametrization, hypervolumes, and post-evaluation
# Author: Mariana Castaneda-Guzman
# Date created: 10/09/2023
# Date last updated: 3/28/2024

# Description: Running Hypervolumes for virginia oysters 


# Packages ----------------------------------------------------------------


library(hypervolume)
library(terra)
library(maps)
library(sf)
library(usdm)
library(ENMeval)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpubr)

# Specify Database --------------------------------------------------------

# Database
db <- "(2003-2023)"

# Input Data --------------------------------------------------------------

# Read mask to make sure same CRS
mask <- rast(paste0("../Data/RS_", db, "/mask_" , db, ".tif"))

# Load occurrences
occ_points <- readRDS(paste0("../Data/GBIF/RS_", db, "/cv_occurrences_corrected_one_per_pixel_weighted_", db, ".RDS"))

head(occ_points)
nrow(occ_points)

occ_points <- occ_points %>% filter(year >= 2003, year <= 2023)

length(unique(occ_points$pixel_index))
nrow(occ_points)

sum(occ_points$obs)

# This is where you could make unique
occ_coords <- st_as_sf(occ_points,
                       coords = c("x", "y"),
                       crs = st_crs(mask))
# Load layers
env_layers <- lapply(list.files(paste0("../Data/RS_", db, "/"), 
                                pattern = "summary", recursive = TRUE, 
                                full.names = TRUE), rast)

env_stack <- rast(env_layers)
names(env_stack) <- gsub("_summary\\.tif|cropped_", "", basename(sources(env_stack)))
# plot(env_stack)

# z-transform climate layers to make axes comparable
env_layers_scaled <- scale(env_stack, center = TRUE, scale = TRUE)
# plot(env_layers_scaled)

# Identify collinear variables that should be excluded
vif_cor <- usdm::vifcor(x = env_layers_scaled, th = 0.8)
vif_cor

# vifstep(occ_extract)

# Make the different set of variables

env_layers_vif <- exclude(env_layers_scaled, vif_cor)

if(length(which(names(env_layers_vif) == "chlo_sd")) > 0){
  env_layers_vif <- env_layers_vif[[-which(names(env_layers_vif) == "chlo_sd")]]
  env_layers_vif <- c(env_layers_vif,
                      env_layers_scaled[["chlo_range"]])
}

env_layers_cor1 <- env_layers_scaled[[c("sst_max",
                             "sst_min",
                             "sst_range",
                             "chlo_max",
                             "chlo_min")]]

env_layers_cor2 <- env_layers_scaled[[c("sst_mean",
                             "sst_sd",
                             "chlo_mean",
                             "chlo_range",
                             "chlo_min")]]

# usdm::vifcor(x = env_layers_cor2, th = 0.8)


# Final set of environmental layers + occurrence points
# Make them raster objects (IMPORTANT to make hypervolume work)
env_layers_vif <- raster::stack(env_layers_vif)
env_layers_cor1 <- raster::stack(env_layers_cor1)
env_layers_cor2 <- raster::stack(env_layers_cor2)


occ_extract_vif <- raster::extract(env_layers_vif, occ_coords)
occ_extract_cor1 <- raster::extract(env_layers_cor1, occ_coords)
occ_extract_cor2 <- raster::extract(env_layers_cor2, occ_coords)

head(occ_extract_cor2)
nrow(occ_extract_cor2)


# Test Combination of Env. Vars ------------------------------------------- This

# This is an optional step, just hepls visualize the general results you might
# get with all the variables Variables

env_combos <- list(env_layers_vif, env_layers_cor1, env_layers_cor2)
occ_combos <- list(occ_extract_vif, occ_extract_cor1, occ_extract_cor2)


eastern_coastal_states <- c("Connecticut", "Delaware", "Florida",
                            "Georgia", "Maine", "Maryland", 
                            "Massachusetts", "New Hampshire", "New Jersey", 
                            "New York", "North Carolina", "Rhode Island", 
                            "South Carolina", "Virginia")

states <- ggplot2::map_data("state")
states$east_coast <- ifelse(states$region %in%  tolower(eastern_coastal_states), 1, NA)

# To plot maps for each year
env_combos_plots <- list()

if(F){
  for(i in 1:length(env_combos)){
    # compute hyper volumes
    hv <- hypervolume_svm(occ_combos[[i]], name = 'v.oyster',
                          svm.nu = 0.01, svm.gamma = 0.5)
    
    # do species distribution modeling
    hv_map <- hypervolume_project(hv,
                                  env_combos[[i]],
                                  type = "inclusion", 
                                  fast.or.accurate = "fast")
    
    
    hv_map_df <- terra::as.data.frame(hv_map, xy = TRUE)
    names(hv_map_df)[3] <- "pred"
    
    hv_map_df <- hv_map_df %>% 
      dplyr::filter(!is.nan(pred))
    
    hv_map_df$pred <- ifelse(hv_map_df$pred >= 1, 1, 0)
    hv_map_df$pred <- as.factor(hv_map_df$pred)
    
    
    p <- ggplot( ) + 
      geom_tile(data = hv_map_df,
                aes(x = x, y = y, fill = pred)) +
      scale_fill_manual(values = c("lightblue", "green4")) +
      geom_polygon(data = states %>% dplyr::filter(!is.na(east_coast)), 
                   aes(x = long, y = lat, group = group),
                   color = "black", fill = "lightgrey") +
      xlab("Longitude") + ylab("Latitude") +
      scale_x_continuous(breaks = seq(-100, -60, 5)) +
      scale_y_continuous(breaks = seq(20, 50, 5)) +
      coord_fixed(1.3) +
      theme_bw()
    
    env_combos_plots[[i]] <- p
    
  }
  
  ggarrange(env_combos_plots[[1]], 
            env_combos_plots[[2]],
            env_combos_plots[[3]],
            common.legend = TRUE, 
            labels = c("A", "B", "C"),
            ncol = 3, 
            legend = "bottom")
}

# Hypervolume Parametrization ---------------------------------------------

set.seed(12345) # to be able to replicate your random draws

# Selecting calbration and evaluation sets

# Using the package ENMeval to select sets

# Enviroemtnal background
envs.bg <- env_layers_vif

# Make occurences a matrix with two columns (lat and long)
occs <- st_coordinates(occ_coords)

# Select a random sample as your background points
bg <- dismo::randomPoints(env_layers_vif[[1]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

# Get your 2-folds for validation and calibration 
folds <- get.checkerboard1(occs, envs.bg, bg, aggregation.factor=30)

# Plotting
evalplot.grps(pts = occs, pts.grp = folds$occs.grp, envs = envs.bg)
evalplot.grps(pts = bg, pts.grp = folds$bg.grp, envs = envs.bg)


# Set up parametrization loop 

# Environmental sets, with different sets of environmental layers
env_combos <- list(env_layers_vif, env_layers_cor1, env_layers_cor2)

# Parameters to test. Change based on preference. This runs 380 combinations,
# and includes default hypervolume parameters nu = 0.01 and gamma = 0.5
params <- expand.grid(nu = seq(0.005, 0.1, 0.005), 
                      gamma = seq(0.1, 1, 0.05))

params

model_results <- data.frame()

set.seed(1234)
message(paste0("Start time: ", format(Sys.time(), "%X")))

for(part in c(1,2)){
  
  if(part == 1){
    data_sf_train <- occ_coords[folds$occs.grp == 1, ] 
    data_sf_test <- occ_coords[folds$occs.grp == 2, ] 
  }else{
    data_sf_train <- occ_coords[folds$occs.grp == 2, ] 
    data_sf_test <- occ_coords[folds$occs.grp == 1, ] 
  }
  
  e_space_train_vif <- raster::extract(env_layers_vif, data_sf_train)
  e_space_test_vif <- raster::extract(env_layers_vif, data_sf_test)
  
  e_space_train_cor1 <- raster::extract(env_layers_cor1, data_sf_train)
  e_space_test_cor1 <- raster::extract(env_layers_cor1, data_sf_test)
  
  e_space_train_cor2 <- raster::extract(env_layers_cor2, data_sf_train)
  e_space_test_cor2 <- raster::extract(env_layers_cor2, data_sf_test)
  
  e_space_test <- list(e_space_test_vif, e_space_test_cor1, e_space_test_cor2)
  e_space_train <- list(e_space_train_vif, e_space_train_cor1, e_space_train_cor2)
  
  sub_model_results <- data.frame()

  for(set in 1:length(e_space_train)){
    message(paste0("In set: ", set, " partition:", part))

    for(i in 1:nrow(params)){
      message(paste0("In set: ", set, " partition:", part, "; combo: ", i, "/", nrow(params)))
      
      hv_train <- hypervolume_svm(e_space_train[[set]], # Points used for model training.
                                  name = "hv_svm", 
                                  svm.nu = params$nu[i], 
                                  svm.gamma = params$gamma[i], 
                                  samples.per.point = 10,
                                  verbose = FALSE) 
      
      # hv_projected_test <- hypervolume_estimate_probability(hv_train, 
      #                                                       e_space_test[[set]], 
      #                                                       verbose = FALSE) 
      
      # create hypervolume projection
      hv_projected_test <- hypervolume_project(hv_train,
                                               env_combos[[set]],
                                    type = "inclusion",
                                    fast.or.accurate = 'fast',
                                    verbose = FALSE)
      
  
      # hv_projected_test  <- raster::raster(paste0("../Outputs/hypervolumes/set", 1, "/hv_svm_map_",
      #                               sprintf("%d.tif", 2003)))
      
      hv_projected_p <- raster::rasterToPoints(hv_projected_test)[,3]
        
      hv_projected_test_p <- raster::extract(hv_projected_test, 
                                             data_sf_test)
      
      hv_projected_test_p[is.na(hv_projected_test_p)] <- 0
      
      p <- c() 
      
      # specify set, nu, gamma
      p$partition <- part
      p$set <- set
      p$nu <- params$nu[i]
      p$gamma <- params$gamma[i]

      # Get volume
      p$volume <- get_volume(hv_train) 
      
      # Calculate the true positives and false positive for the training and test
      # data sets
      
      p$TP_test <- length(hv_projected_test_p[hv_projected_test_p == 1]) 
      
      p$FN_test <- length(hv_projected_test_p[hv_projected_test_p == 0]) 
      

      # Calculate sensitivity and omission
      p$sensitivity_test <- p$TP_test/(p$TP_test + p$FN_test)
      p$omission_test <- 1 - p$sensitivity_test
      
      # Calculate CBP
      # num_successes <- p$TP_test
      p$total_area <- length(hv_projected_p)
      p$total_suitable_area <- length(hv_projected_p[hv_projected_p == 1])
      
      p$CBP <- dbinom(p$TP_test, 
                      size = p$TP_test + p$FN_test, 
                      prob = p$total_suitable_area/p$total_area)
      
      sub_model_results <- rbind(sub_model_results, p)
      
    }
    
  }
  
  model_results <- rbind(model_results, sub_model_results)
}

message(paste0("End time: ", format(Sys.time(), "%X")))

write.csv(model_results, 
          paste0("../Outputs/hypervolumes/model_params_hv_svm.csv"), 
          row.names = FALSE) 


# Figures for Parametrization ---------------------------------------------

# Set theme for plots 
theme_set(theme_bw() +
            theme(
              plot.title = element_text(family = "Times", size = 15),
              
              axis.title = element_text(family = "Times", size = 15),
              axis.text = element_text(family = "Times", size = 13),
              
              legend.title = element_text(family = "Times", size = 15),
              legend.text = element_text(family = "Times", size = 13),
              legend.position = "right",
              
              strip.text.x = element_text(family = "Times", size = 13)
            )
)

# Read in Parametrization
model_results <- read.csv("../Outputs/hypervolumes/model_params_hv_svm.csv")

model_results <- model_results %>%  
  mutate(CBP = dbinom(TP_test, 
                      TP_test + FN_test, 
                      total_suitable_area/total_area))

model_results <- model_results %>% 
  mutate(nu = as.factor(nu), 
         gamma = as.factor(gamma), 
         set = as.factor(set), partition = as.factor(partition),
         CBP_log = log(CBP)) %>% 
  # SALEN NEGATIVOS
  mutate(CBP_log = ifelse(is.infinite(CBP_log), 
                          min(model_results$CBP_log[!is.infinite(model_results$CBP_log)]),
                          CBP_log))


# Obtain parameters that minimize omission and volume, and are less than or
# equal to CBP = 0.001. Also aggregate data per eneviromatal layer set.
df <- model_results %>% 
  group_by(nu, gamma, set) %>% 
  filter(CBP <= 0.001) %>% 
  summarise(omission_test = mean(omission_test),
            volume = mean(volume))

# Models should have been reduce by at least a factor of 2
nrow(df)

# Plot scatter plot
ggplot(df) +
  geom_point(aes(x = omission_test, y = volume)) 

# # If you are running less that 6 different nu and 6 different gammas this plot
# # is prettier
# ggplot(df) + 
#   geom_point(aes(x = omission_test, y = volume, 
#                  color = nu, shape = gamma)) +
#   scale_color_brewer(palette = "Paired")

# Given that you want to minimize x (omission_test) and y (volume), you need to
# find a way to combine these two variables into a single objective function
# that you can minimize. Here's one possible approach:

# 1. Scale Variables (many many options):

# 1.1 Normalizing Transformation: You can use min-max scaling to bring both
# variables to the range [0, 1].

# 1.2 Rank Transformation: Instead of normalizing, you could rank the values of each
# variable. This transformation preserves the order of the data while making it
# less sensitive to extreme values. After ranking, you can then use the ranked
# values to create an objective function for minimization.

# 2. Combine the Variables: Once the variables are on the same scale, you can
# combine them into a single objective function. Since you want to minimize both
# variables, you can simply add them together. However, you might need to adjust
# the weights if one variable is more important than the other.


# Rank transformation
df_trans <- df 
df_trans$omission_t <- rank(df$omission_test, ties.method = "min")
df_trans$volume_t <- rank(df$volume, ties.method = "min")

# Plot scatter plot
ggplot(df_trans) +
  geom_point(aes(x = omission_t, y = volume_t)) 

# Calculate the objective function, your objective function can change depending
# on what you are looking for you might need to adjust the weights if one
# variable is more important than the other. Eg. multiple omission * 2

df_trans <- df_trans %>% 
  mutate(obj_fun = omission_t + volume_t) %>% 
  arrange(obj_fun)

# Save parameter combination with the lowest sum
min_op <- df_trans %>% 
  dplyr::filter(obj_fun == min(df_trans$obj_fun)) 

# table(min_op[, c("nu", "gamma")])

A <- ggplot(df_trans, aes(x = omission_t, y = volume_t)) +
  geom_point() +
  geom_point(data = min_op, aes(color = nu, shape = gamma),
             size = 3) +
  scale_color_brewer(palette = "Set1") +
  ylab("Volume (Rank Transformed)") + 
  xlab("Omission Rate (Rank Transformed)") 
  

B <- ggplot(df_sort,  aes(x = omission_test, y = volume)) +
  geom_point() +
  geom_point(data = min_op, aes(color = nu, shape = gamma),
             size = 3) +
  ylab("Volume") + xlab("Omission Rate") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05))

ggpubr::ggarrange(A, B, common.legend = TRUE)


# Save high resolution figure
tiff(filename = "../Figures/optimal_param.tiff", 
     width = 10, height = 5, units = "in", 
     compression = "lzw", res = 300)
ggpubr::ggarrange(A, B, labels = c("A", "B"),
                  common.legend = TRUE)

dev.off()

# Visualize all sets
in_all_sets <- df_trans %>% 
  dplyr::filter(nu == min_op$nu, gamma == min_op$gamma)


C <- ggplot(df_trans, aes(x = omission_t, y = volume_t)) +
  geom_point() +
  geom_point(data = in_all_sets, aes(color = set),
             size = 3) +
  scale_color_brewer(palette = "Paired") +
  ylab("Volume (Rank Transformed)") + 
  xlab("Omission Rate (Rank Transformed)") 


D <- ggplot(df_sort,  aes(x = omission_test, y = volume)) +
  geom_point() +
  geom_point(data = in_all_sets, aes(color = set),
             size = 3) +
  ylab("Volume") + xlab("Omission Rate") +
  scale_color_brewer(palette = "Paired") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05))

ggpubr::ggarrange(C, D, common.legend = TRUE)


# Save high resolution figure
tiff(filename = "../Figures/optimal_param_set.tiff", 
     width = 10, height = 5, units = "in", 
     compression = "lzw", res = 300)
ggpubr::ggarrange(C, D, labels = c("A", "B"),
                  common.legend = TRUE)

dev.off()



# Other plots for parametrization

# Plots to show distribution of volume and omission per parameter
g1 <- ggplot(model_results, 
       aes(x = nu, y = volume)) +
  geom_boxplot(aes(color = set)) +
  scale_color_brewer(palette = "Paired") +
  ylab("Volume") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
g1

g2 <- ggplot(model_results, 
             aes(x = nu, y = omission_test)) +
  geom_boxplot(aes(color = set)) +
  scale_color_brewer(palette = "Paired") +
  ylab("Omission Rate") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
g2

g3 <- ggplot(model_results, 
             aes(x = gamma, y = volume)) +
  geom_boxplot(aes(color = set)) +
  scale_color_brewer(palette = "Paired") +
  ylab("Volume") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
g3

g4 <- ggplot(model_results, 
             aes(x = gamma, y = omission_test)) +
  geom_boxplot(aes(color = set)) +
  scale_color_brewer(palette = "Paired") +
  ylab("Omission Rate") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
g4


ggpubr::ggarrange(g1,g2,g3,g4, 
                  ncol = 1, 
                  labels = c("A", "B", "C", "D"),
                  common.legend = TRUE)

tiff(filename = "../Figures/param_bh.tiff", width = 8, height = 11, 
     units = "in", compression = "lzw", res = 300)
ggpubr::ggarrange(g1,g2,g3,g4, 
                  ncol = 1, 
                  labels = c("A", "B", "C", "D"),
                  common.legend = TRUE)
dev.off()

# Anova and Tukey's test to check if there are significant differences
mod <- aov(volume ~ nu, model_results)
summary(mod) #no

mod <- aov(omission_test ~ nu, model_results)
summary(mod) #no

mod <- aov(volume ~ gamma, model_results)
summary(mod) #yes
TukeyHSD(mod)

mod <- aov(omission_test ~ gamma, model_results)
summary(mod) #yes
TukeyHSD(mod)

# Plot Tukey's test
plot(TukeyHSD(mod), las = 1)


# OTHER PLOTS
A <- ggplot(model_results, 
       aes(x = omission_test, y = volume, color = set)) +
  geom_point() +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.05, 0.5, 0.05)) +
  scale_color_brewer(palette = "Paired")
  # labs(color = "Env. Set")
A

tiff(filename = "../Figures/mod_res_set.tiff", width = 7, height = 4, 
     units = "in", compression = "lzw", res = 300)
A
dev.off()


# model_results <- model_results %>%  
#   mutate(NCBP = dnbinom(FN_test, 
#                         TP_test, 
#                         total_suitable_area/total_area))
# 
# ggplot(model_results, aes(x = NCBP, y = volume)) +
#   geom_point()


# Hypervolumes per year ---------------------------------------------------

# Run your hypervolume based on the chosen parameters. Default svm.nu = 0.01,
# svm.gamma = 0.5


# STEP 1: specify the parameters to use and method to use, and set of varibales

min_op

method <- "svm"
set <- 2
# if more than one just add in c() list, loop is set up to work with more than
# one combination
param_combo <- data.frame(nu = 0.065,
                          gamma = 0.5)


# Specify the number of variables to use
# STEP 2:

set.seed(1234)

env_combos <- list(env_layers_vif, env_layers_cor1, env_layers_cor2)
occ_combos <- list(occ_extract_vif, occ_extract_cor1, occ_extract_cor2)

for(combo in 1:length(param_combo)){
  
  # STEP 3: project to each year
  for(set in 1:length(occ_combos)){
    
    message(paste0("Param combo: ", combo, "; Set: ", set, "; Start time: ", format(Sys.time(), "%X")))
    
    # Create the main hypervolume
    hv <- hypervolume_svm(data = occ_combos[[set]],
                          svm.nu = param_combo$nu[combo],
                          svm.gamma = param_combo$gamma[combo],
                          verbose = FALSE)
    
    # Save hypervolume
    saveRDS(hv, paste0("../Outputs/hypervolumes/hv_svm_set_", set,
                       "_nu_", param_combo$nu[combo] , 
                       "_gamma_", param_combo$gamma[combo], 
                       "_", db, ".RDS"))
    
    for (year in c(2003:2023)){
      # to keep track of the year being projected
      message(paste0("Set: ", set, "; year: ", year, "; time:", format(Sys.time(), "%X")))
      
      
      # List all important variables
      if(set == 1){
        vars <- c("CHLO_%d_min.tif",
                  "CHLO_%d_range.tif",
                  "SST_%d_max.tif",
                  "SST_%d_min.tif",
                  "SST_%d_range.tif")
      }else if(set == 2){
        vars <- c("CHLO_%d_max.tif",
                  "CHLO_%d_min.tif",
                  "SST_%d_max.tif",
                  "SST_%d_min.tif",
                  "SST_%d_range.tif")
      }else{
        vars <- c("CHLO_%d_mean.tif",
                  "CHLO_%d_min.tif",
                  "CHLO_%d_range.tif",
                  "SST_%d_mean.tif",
                  "SST_%d_sd.tif")
      }
      
      # List all the files in the directory, this will help subset in the loop
      vars_files <- list.files(path = "../Data/RS_(2003-2023)/", 
                               patter = "\\.tif$", 
                               recursive = TRUE, 
                               full.names = TRUE)
      
      # Assign variable names, this are going to be used as column names in the line
      # below
      # List all important variables
      if(set == 1){
        vars_name <- c("CHLO_min",
                       "CHLO_range",
                       "SST_max",
                       "SST_min",
                       "SST_range")
      }else if(set == 2){
        vars_name <- c("CHLO_max",
                       "CHLO_min",
                       "SST_max",
                       "SST_min",
                       "SST_range")
      }else{
        vars_name <- c("CHLO_mean",
                       "CHLO_min",
                       "CHLO_range",
                       "SST_mean",
                       "SST_sd")
      }
      
      # stack all environmental layers for the given year
      env <- raster::stack(vars_files[which(basename(vars_files) %in% sprintf(vars, year))])
      names(env) <- vars_name
      
      env_layers_clean <- scale(env, center = TRUE, scale = TRUE)
      
      
      # create hypervolume projection
      result <- hypervolume_project(hv,
                                    env_layers_clean,
                                    type = "inclusion",
                                    fast.or.accurate='accurate',
                                    verbose = FALSE)
      
      
      # Write the resulting projection in G space
      writeRaster(x = result,
                  filename = paste0("../Outputs/hypervolumes/set", set, 
                  "/hv_svm_nu_", param_combo$nu[combo], 
                  "_gamma_",  param_combo$gamma[combo], 
                  "_map_", year, ".tif"),
                  overwrite = TRUE)
    }
    
    message(paste0("Param combo: ", combo, "; Set: ", set, "; End time: ", format(Sys.time(), "%X")))
    
  }
}


# Post model analysis -----------------------------------------------------

# Set theme for plots 
theme_set(theme_bw() +
            theme(
              plot.title = element_text(family = "Times", size = 15),
              
              axis.title = element_text(family = "Times", size = 15),
              axis.text = element_text(family = "Times", size = 13),
              
              legend.title = element_text(family = "Times", size = 15),
              legend.text = element_text(family = "Times", size = 13),
              legend.position = "right",
              
              strip.text.x = element_text(family = "Times", size = 13)
            )
)


eastern_coastal_states <- c("Connecticut", "Delaware", 
                            "Florida", "Georgia", "Maine",
                            "Maryland", "Massachusetts", 
                            "New Hampshire", "New Jersey",
                            "New York", "North Carolina", 
                            "Rhode Island", "South Carolina", "Virginia")

states <- map_data("state")
states$east_coast <- ifelse(states$region %in%  tolower(eastern_coastal_states), 1, NA)

# For levelplots
state_polygon <- states %>%
  st_as_sf(coords = c("long", "lat"), crs = st_crs(mask)) %>%
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
state_polygon <- sf:::as_Spatial(state_polygon$geometry)

east_state_polygon <- states %>%
  filter(east_coast == 1) %>% 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(mask)) %>%
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
east_state_polygon <- sf:::as_Spatial(east_state_polygon$geometry)



set <- "set1"

if(!dir.exists(paste0("../Figures/", set))) dir.create(paste0("../Figures/", set))

year_models <- lapply(list.files(paste0("../Outputs/hypervolumes/", set), 
                                 pattern = "hv_svm_", 
                                 full.names = TRUE), rast)

year_models <- lapply(year_models,
                         FUN = function(r){
                           m <- c(-Inf, 0, 0,
                                  0, Inf, 1)
                           rclmat <- matrix(m, ncol=3, byrow=TRUE)
                           r_rc <- classify(r, rclmat, include.lowest = TRUE)
                           return(r_rc)
                         })


# year_models <- year_models[-21]

models_summaries <- data.frame()
latitude_df <- data.frame()

for(ras in year_models){
  new_row <- data.frame(year = as.numeric(gsub("hv_svm_map_", "", names(ras))),
                        suitability = sum(values(ras, na.rm = TRUE) >= 1))
  
  ras_df <- as.data.frame(ras, xy = TRUE)
  
  year_lats <- data.frame(y = ras_df %>% 
                            filter(ras_df[,3] >= 1) %>% 
                            dplyr::select(y),
                          year = as.numeric(gsub("hv_svm_map_", "", names(ras))))
  
  latitude_df <- rbind(latitude_df, year_lats)
  
  
  long_lat_mean <- ras_df %>% 
    filter(ras_df[,3] >= 1) %>% 
    summarise(mean_long = mean(x), 
              mean_lat = mean(y),
              median_lat = median(y)) %>% 
    cbind(ras_df %>% 
            filter(ras_df[,3] >= 1, x >= -81.7975) %>% 
            summarise(mean_lat_x_81 = mean(y),
                      median_lat_x_81 = median(y))) %>%   
    as.data.frame()
  
  new_row <- cbind(new_row, long_lat_mean)
  
  models_summaries <- rbind(models_summaries, new_row)               
}



mod_A <- lm(suitability~year, models_summaries)
summary(mod_A)

fit_changepoint = cpt.mean(df$y)

# Return estimates
c(ints = param.est(fit_changepoint)$mean,
  cp = cpts(fit_changepoint))

percent <- format(abs(diff(param.est(fit_changepoint)$mean))/
  param.est(fit_changepoint)$mean[1] *100, digits = 3)


plot(fit_changepoint)

A <- ggplot(models_summaries, aes(x = year, y = suitability)) +
  geom_point() +
  geom_line() +
  xlab("Year") + ylab(expression(italic("C. virginica")~"Predicted Area (4"~km^2~")")) +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) +
  geom_segment(mapping = aes(x = 2003, y = param.est(fit_changepoint)$mean[1], 
                             xend = 2011, yend = param.est(fit_changepoint)$mean[1]),
               color = "red", linewidth = 1) +
  geom_segment(mapping = aes(x = 2011, y = param.est(fit_changepoint)$mean[2], 
                             xend = 2023, yend = param.est(fit_changepoint)$mean[2]),
               color = "red", linewidth = 1) +
  geom_segment(mapping = aes(x = 2011, y = param.est(fit_changepoint)$mean[1], 
                             xend = 2011, yend = param.est(fit_changepoint)$mean[2]),
               color = "red", linetype = 2) +
  annotate("text", x = 2010, y = 28500, 
           label = paste0(percent, "%"), color = "red") 
A


tiff(filename = paste0("../Figures/", set, "/figA_(pred_area).tiff"), 
     width = 7, height = 4, units = "in", res = 300, compression = "lzw")
A
dev.off()

mod_B <- lm(mean_lat~year, models_summaries)
summary(mod_B)

B <- ggplot(models_summaries, aes(x = year, y = mean_lat)) +
  geom_point() +
  geom_line() +
  geom_smooth(color = "red", method = "lm") +
  xlab("Year") + ylab("Mean Latitude") +
  annotate("text", x = 2008, y = 31.5,
           label = substitute(paste("Adj ", R^2, " = ", val1, ", ",
                                    beta, " = ", val2, " ", italic("p"), " = ",
                                   val3), 
                              list(val1 = signif(summary(mod_B)$adj.r.squared, 3), 
                                               val2 = signif(mod_B$coef[[2]], 3), 
                                               val3 = signif(summary(mod_B)$coef[2, 4], 3))),
           color = "red") +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) +
  scale_y_continuous(breaks = seq(30, 33, 0.35))
B


mod_B1 <- lm(median_lat~year, models_summaries)
summary(mod_B1)

B1 <- ggplot(models_summaries, aes(x = year, y = median_lat)) +
  geom_point() +
  geom_line() +
  geom_smooth(color = "red", method = "lm") +
  xlab("Year") + ylab("Median Latitude") +
  annotate("text", x = 2008, y = 29,
           label = substitute(paste("Adj ", R^2, " = ", val1, ", ",
                                    beta, " = ", val2, " ", italic("p"), " = ",
                                    val3), 
                              list(val1 = signif(summary(mod_B1)$adj.r.squared, 3), 
                                   val2 = signif(mod_B1$coef[[2]], 3), 
                                   val3 = signif(summary(mod_B1)$coef[2, 4], 3))),
           color = "red") +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) +
  scale_y_continuous(breaks = seq(28, 31.5, 0.5))

B1


mod_B2 <- lm(mean_lat_x_81~year, models_summaries)
summary(mod_B2)

B2 <- ggplot(models_summaries, aes(x = year, y = mean_lat_x_81)) +
  geom_point() +
  geom_line() +
  geom_smooth(color = "red", method = "lm") +
  xlab("Year") + ylab("Mean Latitude (lon > -81)") +
  labs(title = paste0("Adj R2=",signif(summary(mod_B2)$adj.r.squared, 3),
                      " B=",signif(mod_B2$coef[[2]], 3),
                      " P=",signif(summary(mod_B2)$coef[2,4], 3))) +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) 
B2


ggpubr::ggarrange(B, B1, labels = c("A", "B"), nrow = 2)


tiff(filename = paste0("../Figures/", set, "/figB_(pred_lat).tiff"), 
     width = 7, height = 6, units = "in", res = 300, compression = "lzw")
ggpubr::ggarrange(B, B1, labels = c("A", "B"), nrow = 2)
dev.off()



# LATITUDE DENSITIES

mod_B3 <- lm(y~year, latitude_df)

B3 <- ggplot(latitude_df %>% mutate(year = as.factor(year)),
             aes(x = year, y = y, fill = year)) +
  geom_violin() +
  scale_fill_viridis_d(direction = -1) +
  xlab("Year") + ylab("Latitude") +
  theme(legend.position = "none")

  
B3

B4 <- ggplot(latitude_df %>% 
         mutate(year = as.factor(year)), 
       aes(x = y, y = year, fill = year)) +
  geom_density_ridges() +
  scale_fill_viridis_d(direction = -1) +
  ylab("Year") + xlab("Latitude") +
  coord_flip() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(20, 50, 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

B4


ggpubr::ggarrange(ggpubr::ggarrange(B, B1, labels = c("A", "B"), nrow = 2),
                  B4, ncol = 2)


tiff(filename = paste0("../Figures/", set, "/figB2_(pred_lat_dens).tiff"), 
     width = 7, height = 4, units = "in", res = 300, compression = "lzw")
B4
dev.off()


B5 <- ggplot(latitude_df %>% 
               mutate(year = as.factor(year)) %>% 
               filter(year %in% c(2003, 2013, 2022)), 
             aes(x = y, group = year, fill = year)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(23, 47, 3)) +
  ylab("Density") + xlab("Latitude")

B5


tiff(filename = paste0("../Figures/", set, "/figB3_(pred_lat_dens_3y).tiff"), 
     width = 7, height = 4, units = "in", res = 300, compression = "lzw")
B5
dev.off()

mod_B6 <- aov(y~year, latitude_df %>% mutate(year = as.factor(year)))
summary(mod_B6)

mod_B6_Tukey <- TukeyHSD(mod_B6)
plot(mod_B6_Tukey, las = 1)


# MAPS
# Plots

result_df <- lapply(year_models, 
                    FUN = function(i) terra::as.data.frame(i, xy = TRUE))

# Year 2003

result_df_2003 <- result_df[[1]] %>% dplyr::filter(result_df[[1]][,3] <= 1)

result_df_2003[,3] <- as.factor(result_df_2003[,3])

C <- ggplot( ) + 
  geom_tile(data = result_df_2003,
            aes(x = x, y = y, fill = result_df_2003[,3])) +
  scale_fill_manual(values = c("lightblue", "green4")) +
  geom_polygon(data = states %>% dplyr::filter(!is.na(east_coast)), 
               aes(x = long, y = lat, group = group),
               color = "black", fill = "lightgrey") +
  geom_hline(models_summaries[c(20,1), ], 
             mapping = aes(yintercept = mean_lat,
                           color = as.factor(year), 
                           linetype = as.factor(year)),
             linewidth = 1) +
  scale_color_manual(values = c("orange", "tomato")) +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(breaks = seq(-100, -60, 5)) +
  scale_y_continuous(breaks = seq(23, 47, 3)) +
  coord_fixed(1.3) +
  ggtitle(sprintf("%d", 2003)) +
  labs(fill = "Area", color = "Mean\nLatitude", linetype = "Mean\nLatitude")
C


C_l <- levelplot(year_models[[1]], 
               xlab = list("Longitude", cex = 2),
               ylab = list("Latitude", cex = 2),
               main = list("2003", cex = 2),
               colorkey = FALSE,
               scales = list(cex = 1.5),
               margin = list(draw = TRUE, 
                             axis = list(col = 'black', fontsize = 13)), 
               col.regions = c("lightblue", "green4")) + 
  latticeExtra::layer(sp.polygons(state_polygon, lwd = 1, col='lightgrey')) +
  latticeExtra::layer(sp.polygons(east_state_polygon, lwd = 1, fill='lightgrey'))

C_l



# Year 2013
result_df_2013 <- result_df[[11]] %>% dplyr::filter(result_df[[11]][,3] <= 1)

result_df_2013[,3] <- as.factor(result_df_2013[,3])

C1 <- ggplot( ) + 
  geom_tile(data = result_df_2013,
            aes(x = x, y = y, fill = result_df_2013[,3])) +
  scale_fill_manual(values = c("lightblue", "green4")) +
  geom_polygon(data = states %>% dplyr::filter(!is.na(east_coast)), 
               aes(x = long, y = lat, group = group),
               color = "black", fill = "lightgrey") +
  geom_hline(models_summaries[c(20,1), ], 
             mapping = aes(yintercept = mean_lat,
                           color = as.factor(year), 
                           linetype = as.factor(year)),
             linewidth = 1) +
  scale_color_manual(values = c("orange","tomato")) +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(breaks = seq(-100, -60, 5)) +
  scale_y_continuous(breaks = seq(23, 47, 3)) +
  coord_fixed(1.3) +
  ggtitle(sprintf("%d", 2013)) +
  labs(fill = "Area", color = "Mean\nLatitude", linetype = "Mean\nLatitude")
C1

# Levelplot
C1_l <- levelplot(year_models[[11]], 
                 xlab = list("Longitude", cex = 2),
                 ylab = list("Latitude", cex = 2),
                 main = list("2013", cex = 2),
                 colorkey = FALSE,
                 scales = list(cex = 1.5),
                 margin = list(draw = TRUE, 
                               axis = list(col = 'black', fontsize = 13)), 
                 col.regions = c("lightblue", "green4")) + 
  latticeExtra::layer(sp.polygons(state_polygon, lwd = 1, col='lightgrey')) +
  latticeExtra::layer(sp.polygons(east_state_polygon, lwd = 1, fill='lightgrey'))

C1_l

# Year 2022
result_df_2023 <- result_df[[20]] %>% dplyr::filter(result_df[[20]][,3] <= 1)

result_df_2023[,3] <- as.factor(result_df_2023[,3])

C2 <- ggplot( ) + 
  geom_tile(data = result_df_2023,
            aes(x = x, y = y, fill = result_df_2022[,3])) +
  scale_fill_manual(values = c("lightblue", "green4")) +
  geom_polygon(data = states %>% dplyr::filter(!is.na(east_coast)), 
               aes(x = long, y = lat, group = group),
               color = "black", fill = "lightgrey") +
  geom_hline(models_summaries[c(21,1), ], 
             mapping = aes(yintercept = mean_lat,
                           color = as.factor(year), 
                           linetype = as.factor(year)),
             linewidth = 1) +
  scale_color_manual(values = c("orange", "tomato")) +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(breaks = seq(-100, -60, 5)) +
  scale_y_continuous(breaks = seq(23, 47, 3)) +
  coord_fixed(1.3) +
  ggtitle(sprintf("%d", 2023)) +
  labs(fill = "Area", color = "Mean\nLatitude", linetype = "Mean\nLatitude")
C2

C2_l <- levelplot(year_models[[21]], 
                  xlab = list("Longitude", cex = 2),
                  ylab = list("Latitude", cex = 2),
                  main = list("2023", cex = 2),
                  colorkey = FALSE,
                  scales = list(cex = 1.5),
                  margin = list(draw = TRUE, 
                                axis = list(col = 'black', fontsize = 13)), 
                  col.regions = c("lightblue", "green4")) + 
  latticeExtra::layer(sp.polygons(state_polygon, lwd = 1, col='lightgrey')) +
  latticeExtra::layer(sp.polygons(east_state_polygon, lwd = 1, fill='lightgrey'))

C2_l

C_maps_l <- ggpubr::ggarrange(C_l, C1_l, C2_l, 
                  nrow = 1)


tiff(filename = paste0("../Figures/", set, "/figC_l_(maps_3y).tiff"), 
     width = 10, height = 7, units = "in", res = 300, compression = "lzw")
C_maps_l
dev.off()


C_maps <- ggpubr::ggarrange(C, C1, C2, legend = "right", 
                            common.legend = TRUE,
                            nrow = 1)

tiff(filename = paste0("../Figures/", set, "/figC_(maps_3y).tiff"), 
     width = 10, height = 7, units = "in", res = 300, compression = "lzw")
C_maps
dev.off()


# Areas suitable across time
asat <- prod(rast(year_models))
plot(asat)

asat_df <- as.data.frame(asat, xy = TRUE)

asat_df <- asat_df %>% dplyr::filter(asat_df[,3] <= 1)

asat_df[,3] <- as.factor(asat_df[,3])

D <- ggplot( ) + 
  geom_tile(data = asat_df,
            aes(x = x, y = y, fill = asat_df[,3])) +
  scale_fill_manual(values = c("lightblue", "green4")) +
  geom_polygon(data = states %>% dplyr::filter(!is.na(east_coast)), 
               aes(x = long, y = lat, group = group),
               color = "black", fill = "lightgrey") +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(breaks = seq(-100, -60, 5)) +
  scale_y_continuous(breaks = seq(23, 47, 3)) +
  coord_fixed(1.3) +
  ggtitle("Areas suitable across time") +
  theme(legend.position = "none")
D

D_l <- levelplot(asat, 
                  xlab = list("Longitude", cex = 2),
                  ylab = list("Latitude", cex = 2),
                  main = list("Areas suitable across time", cex = 2),
                  colorkey = FALSE,
                  scales = list(cex = 1.5),
                  margin = list(draw = TRUE, 
                                axis = list(col = 'black', fontsize = 13)), 
                  col.regions = c("lightblue", "green4")) + 
  latticeExtra::layer(sp.polygons(state_polygon, lwd = 1, col='lightgrey')) +
  latticeExtra::layer(sp.polygons(east_state_polygon, lwd = 1, fill='lightgrey'))

D_l



# Luego sumar todos los raters,
# y hacer 1-la suma, para identificar zonas que se
# han perdido con el tiempo o shifts.

rshift <- 1 - sum(rast(year_models))
plot(rshift)

library(rasterVis)
library(viridisLite)

revMagma <- rasterTheme(region = rev(magma(10)))

D2_l <- levelplot(rshift, 
                 xlab = list("Longitude", cex = 2),
                 ylab = list("Latitude", cex = 2),
                 main = list("Shifts across time", cex = 2),
                 par.settings = revMagma,
                 scales = list(cex = 1.5),
                 margin = list(draw = TRUE, 
                               axis = list(col = 'black', fontsize = 13))) + 
  latticeExtra::layer(sp.polygons(state_polygon, lwd = 1, col='lightgrey')) +
  latticeExtra::layer(sp.polygons(east_state_polygon, lwd = 1, fill='lightgrey'))

D2_l


# More ENMeval function testing ------------------------------------------------


# We can increase the aggregation factor to give the groups bigger boxes.
# This can result in groups that are more environmentally different from each other.
envs.bg <- env_layers_vif

occs <- st_coordinates(occ_coords)

bg <- dismo::randomPoints(env_layers_vif[[1]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

plot(envs.bg[[1]])
points(bg, pch = 20, cex = 0.2)


block <- get.block(occs, bg)

evalplot.grps(pts = occs, pts.grp = block$occs.grp, envs = envs.bg) + 
  ggplot2::ggtitle("Spatial block partitions: occurrences")

evalplot.grps(pts = bg, pts.grp = block$bg.grp, envs = envs.bg) + 
  ggplot2::ggtitle("Spatial block partitions: background")


folds <- get.checkerboard1(occs, envs.bg, bg, aggregation.factor=30)

evalplot.grps(pts = occs, pts.grp = folds$occs.grp, envs = envs.bg)

evalplot.grps(pts = bg, pts.grp = folds$bg.grp, envs = envs.bg)






