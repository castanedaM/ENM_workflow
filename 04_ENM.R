# Title: ENMs parametrization, hypervolumes, and post-evaluation
# Author: Mariana Castaneda-Guzman
# Date created: 10/09/2023
# Date last updated: 3/15/2024

# Description: Running Hypervolumes for virginia oysters 

library(hypervolume)
library(terra)
library(maps)
library(sf)
library(usdm)
# library(stringr)
# install.packages("ENMeval")
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

occ_points <- occ_points %>% filter(year >=2003, year <= 2023)

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
plot(env_stack)

# z-transform climate layers to make axes comparable
env_layers_scaled <- scale(env_stack, center = TRUE, scale = TRUE)
plot(env_layers_scaled)

# Identify collinear variables that should be excluded
vif_cor <- usdm::vifcor(x = env_layers_scaled, th = 0.8)
vif_cor

# vifstep(occ_extract)

# Make the different set of variables

env_layers_vif <- exclude(env_layers_scaled, vif_cor)
env_layers_vif <- env_layers_vif[[-which(names(env_layers_vif) == "chlo_sd")]]
# env_layers_vif <- c(env_layers_vif, 
                    # env_layers_scaled[["chlo_cropped_range_summary"]])

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


# Test Combination of Env. Vars -------------------------------------------

# Variables
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
                                  type = "inclusion")
    
    
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

envs.bg <- env_layers_vif
occs <- st_coordinates(occ_coords)
bg <- dismo::randomPoints(env_layers_vif[[1]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(occs)

cb1 <- get.checkerboard1(occs, envs.bg, bg, aggregation.factor=30)

evalplot.grps(pts = occs, pts.grp = cb1$occs.grp, envs = envs.bg)
# evalplot.grps(pts = bg, pts.grp = cb1$bg.grp, envs = envs.bg)

env_combos <- list(env_layers_vif, env_layers_cor1, env_layers_cor2)

params <- expand.grid(nu = seq(0.005, 0.02, 0.005), 
                      gamma = seq(0.25, 0.75, 0.10))

params

model_results <- data.frame()

for(part in c(1,2)){
  
  if(part == 1){
    data_sf_train <- occ_coords[cb1$occs.grp == 1, ] 
    data_sf_test <- occ_coords[cb1$occs.grp == 2, ] 
  }else{
    data_sf_train <- occ_coords[cb1$occs.grp == 2, ] 
    data_sf_test <- occ_coords[cb1$occs.grp == 1, ] 
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
                                    fast.or.accurate = 'accurate',
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


write.csv(model_results, 
          paste0("../Outputs/hypervolumes/model_params_hv_svm.csv"), 
          row.names = FALSE) 


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
  mutate(CBP_log = ifelse(is.infinite(CBP_log), 
                          min(model_results$CBP_log[!is.infinite(model_results$CBP_log)]),
                          CBP_log))

df <- model_results %>% 
  group_by(nu, gamma) %>% 
  summarise(o_max = max(omission_test),
         o_mean = mean(omission_test),
         o_min = min(omission_test),
         v_max = max(volume), 
         v_mean = mean(volume),
         v_min = min(volume))


ggplot(df, 
       aes(x = o_mean, y = v_mean, color = nu, shape = gamma)) +
  geom_segment(aes(x = o_min, y = v_mean, xend = o_max, yend = v_mean)) +
  geom_segment(aes(x = o_mean, y = v_min, xend = o_mean, yend = v_max)) +
  geom_point(size = 3) +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05))


g1 <- ggplot(model_results, 
       aes(x = nu, y = volume)) +
  geom_boxplot() +
  geom_point(aes(shape = gamma, color = set), position = "jitter")
g2 <- ggplot(model_results, 
             aes(x = nu, y = omission_test)) +
  geom_boxplot() +
  geom_point(aes(shape = gamma, color = set), position = "jitter")
g3 <- ggplot(model_results, 
             aes(x = gamma, y = volume)) +
  geom_boxplot() +
  geom_point(aes(shape = nu, color = set), position = "jitter")
g4 <- ggplot(model_results, 
             aes(x = gamma, y = omission_test)) +
  geom_boxplot() +
  geom_point(aes(shape = nu, color = set), position = "jitter")


ggpubr::ggarrange(g1,g2,g3,g4)



mod <- aov(volume ~ nu, model_results)
summary(mod)
TukeyHSD(mod)
plot(TukeyHSD(mod), las = 1)


ggplot(model_results %>% filter(volume <= 2000, 
                                omission_test <= 0.2, 
                                CBP <= 0.001, set != 3),
            aes(x = omission_test, y = volume, color = nu, shape = gamma)) +
  geom_point() +
  ylab("Volume") + xlab("Omission Rate")


a <- ggplot(model_results, 
       aes(x = omission_test, y = volume, color = set)) +
  geom_point() +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05)) +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Env. Set")
a

b <- ggplot(model_results, 
            aes(x = omission_test, y = volume, color = nu, shape = gamma)) +
  geom_point() +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05))

b


c <- ggplot(model_results %>% 
              filter(CBP <= 0.001) %>% 
              mutate(nu = as.factor(nu), gamma = as.factor(gamma)), 
            aes(x = omission_test, y = volume, size = CBP_log, color = CBP_log)) +
  geom_point(alpha = 0.75) +
  scale_color_distiller(palette = "Greens") +
  ylab("Volume") + xlab("Omission Rate") +
  theme_bw()
c


c2 <- ggplot(model_results %>% 
              filter(CBP_log <= -600, volume <= 2000, omission_test <= 0.20) %>% 
              mutate(nu = as.factor(nu), gamma = as.factor(gamma)), 
            aes(x = CBP_log, y = volume, color = nu, shape = gamma)) +
  geom_point(alpha = 0.75) +
  # scale_color_distiller(palette = "Reds") +
  ylab("Volume") + xlab("CBP (log)") +
  theme_bw()
c2

e <- ggplot(model_results %>% 
              filter(CBP_log <= -650), 
            aes(x = omission_test, y = volume, size = CBP_log, color = CBP_log)) +
  geom_point() +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05))
e


d <- ggplot(model_results %>% 
              filter(volume <= 2000, CBP <= 0.001, 
                     set != 3, 
                     omission_test <= 0.20), 
            aes(x = omission_test, y = volume, color = nu, shape = gamma)) +
  geom_point() +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.1, 0.2, 0.05)) +
  facet_wrap(~partition)
d

ggpubr::ggarrange(ggpubr::ggarrange(a, b), c, nrow = 2)


# g1 <- model_results %>% filter(set == 1, nu == 0.02, gamma == 0.45)
# geom_point(data=g1, colour="blue") +



ggpubr::ggarrange(a, b, e, d, labels = c("A", "B", "C", "D"))

model_results <- model_results %>%  
  mutate(NCBP = dnbinom(FN_test, 
                        TP_test, 
                        total_suitable_area/total_area))

ggplot(model_results, aes(x = NCBP, y = volume)) +
  geom_point()


# Hypervolumes per year ---------------------------------------------------

# STEP 1: specify the parameters to use and method to use
# nu <- 0.005
# gamma <- 0.3
method <- "svm"

# svm.nu = 0.01, svm.gamma = 0.5
nu <- 0.01
gamma <- 0.45

# Specify the number of variables to use
# STEP 2:

set.seed(1234)


env_combos <- list(env_layers_vif, env_layers_cor1, env_layers_cor2)
occ_combos <- list(occ_extract_vif, occ_extract_cor1, occ_extract_cor2)


# STEP 3: project to each year
for(set in 1:length(occ_combos)){
  
  message(paste0("Set: ", set, "; start time: ", format(Sys.time(), "%X")))
  
  # Create the main hypervolume
  hv <- hypervolume_svm(data = occ_combos[[set]],
                        svm.nu = nu,
                        svm.gamma = gamma,
                        verbose = FALSE)
  
  # Save hypervolume
  saveRDS(hv, paste0("../Outputs/hypervolumes/set", set, "/hv_svm_", db, ".RDS"))
  
  for (year in c(2003:2022)){
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
                filename = paste0("../Outputs/hypervolumes/set", set, "/hv_svm_map_",
                                  sprintf("%d.tif", year)),
                overwrite = TRUE)
  }
  
  message(paste0("Set: ", set, "; end time: ", format(Sys.time(), "%X")))
  
}



# result_df <- terra::as.data.frame(result, xy = TRUE)
# names(result_df)[3] <- "pred"
# 
# table(result_df$pred)
# 
# result_df <- result_df %>% dplyr::filter(pred <= 1)
# 
# result_df$pred <- as.factor(result_df$pred)
# 
# f1 <- ggplot( ) + 
#   geom_tile(data = result_df,
#             aes(x = x, y = y, fill = pred)) +
#   scale_fill_manual(values = c("lightblue", "green4")) +
#   geom_polygon(data = states %>% dplyr::filter(!is.na(east_coast)), 
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = "lightgrey") +
#   xlab("Longitude") + ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-100, -60, 5)) +
#   scale_y_continuous(breaks = seq(22, 47, 5)) +
#   coord_fixed(1.3) +
#   ggtitle(sprintf("%d", year)) 
# f1
# 
# # ggpubr::ggarrange(a, f1)
# 
# 
# tiff(paste0("../Figures/RS_", db, "/hv_svm_map_",sprintf("%d.tiff", year)),
#      units = "in", width = 5, res = 300,
#      height = 6, compression = "lzw")
# f1
# dev.off()


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
              
              # panel.grid.major.x = element_blank(),
              # panel.grid.minor.x = element_blank(),
              # 
              # panel.grid.major.y = element_line(color = "grey"),
              # panel.grid.minor.y = element_line(color = "lightgrey"),
              
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

set <- "set1"


if(!dir.exists(paste0("../Figures/", set))) dir.create(paste0("../Figures/", set))

year_models <- lapply(list.files(paste0("../Outputs/hypervolumes/", set), 
                                 pattern = "hv_svm_map_", 
                                 full.names = TRUE), rast)


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

A <- ggplot(models_summaries, aes(x = year, y = suitability)) +
  geom_point() +
  geom_line() +
  geom_smooth(color = "red", method = "lm") +
  labs(title = paste0("Adj R2=",signif(summary(mod_A)$adj.r.squared, 3),
                      " B=",signif(mod_A$coef[[2]], 3),
                      " P=",signif(summary(mod_A)$coef[2,4], 3))) +
  xlab("Year") + ylab(expression(italic("C. virginica")~"Predicted Area (km^2)")) +
  scale_x_continuous(breaks = seq(2003, 2022, 2)) 
A


mod_B <- lm(mean_lat~year, models_summaries)
summary(mod_B)

B <- ggplot(models_summaries, aes(x = year, y = mean_lat)) +
  geom_point() +
  geom_line() +
  geom_smooth(color = "red", method = "lm") +
  xlab("Year") + ylab("Mean Latitude") +
  labs(title = paste0("Adj R2=",signif(summary(mod_B)$adj.r.squared, 3),
                      " B=",signif(mod_B$coef[[2]], 3),
                      " P=",signif(summary(mod_B)$coef[2,4], 3))) +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) 
B


mod_B1 <- lm(median_lat~year, models_summaries)
summary(mod_B1)

B1 <- ggplot(models_summaries, aes(x = year, y = median_lat)) +
  geom_point() +
  geom_line() +
  geom_smooth(color = "red", method = "lm") +
  xlab("Year") + ylab("Median Latitude") +
  labs(title = paste0("Adj R2=",signif(summary(mod_B1)$adj.r.squared, 3),
                      " B=",signif(mod_B1$coef[[2]], 3),
                      " P=",signif(summary(mod_B1)$coef[2,4], 3))) +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) 
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


ggpubr::ggarrange(B, B1, B2, labels = c("A", "B", "C"), nrow = 3)


# latitude_df <- latitude_df %>% mutate(year = as.factor(year)

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
  # theme_ridges() + 
  scale_fill_viridis_d(direction = -1) +
  ylab("Year") + xlab("Latitude") +
  theme(legend.position = "none")

B4


B5 <- ggplot(latitude_df %>% 
               mutate(year = as.factor(year)) %>% 
               filter(year %in% c(2003, 2013, 2022)), 
             aes(x = y, group = year, fill = year)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_continuous(breaks = seq(23, 47, 3)) +
  ylab("Year") + xlab("Latitude") 
  
B5

mod_B6 <- aov(y~year, latitude_df %>% mutate(year = as.factor(year)))
summary(mod_B6)

mod_B6_Tukey <- TukeyHSD(mod_B6)
plot(mod_B6_Tukey, las = 1)



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

# Year 2022
result_df_2022 <- result_df[[20]] %>% dplyr::filter(result_df[[20]][,3] <= 1)

result_df_2022[,3] <- as.factor(result_df_2022[,3])

C2 <- ggplot( ) + 
  geom_tile(data = result_df_2022,
            aes(x = x, y = y, fill = result_df_2022[,3])) +
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
  ggtitle(sprintf("%d", 2022)) +
  labs(fill = "Area", color = "Mean\nLatitude", linetype = "Mean\nLatitude")
C2


C_maps <- ggpubr::ggarrange(C, C1, C2, 
                  # labels = c("B", "B", "C"),
                  nrow = 1, 
                  common.legend = TRUE,
                  legend = "right")

ggpubr::ggarrange(A, B5, C_maps, nrow = 2) 


tiff(filename = paste0("../Figures/", set, "/lat_3y_panel.tiff"), 
     width = 10, height = 7, units = "in", res = 300, compression = "lzw")
ggpubr::ggarrange(B5, C_maps, nrow = 2) 
dev.off()


G <- ggplot(models_summaries) + 
  geom_hline(mapping = aes(yintercept = mean_lat,
                           color = year),
             linewidth = 2, alpha = 0.7) +
  scale_color_viridis_c() +
  ylab("Latitude") +
  labs(color = "Year")
G


tiff(filename = paste0("../Figures/", set, "/obs_hv_lat_lines.tiff"), 
     width = 4, height = 4, units = "in", res = 300, compression = "lzw")
G
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



# Luego sumar todos los raters, y hacer 1-la suma, para identificar zonas que se
# han perdido con el tiempo o shifts.
tmpfilter <- rast(year_models) > 1

filtered_image <- mask(inputimage, tmpfilter, maskvalue=1)

plot(rshift)

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


# Evaluation --------------------------------------------------------------


# Intersection
hypervolume_set()

# Centroid distance
hypervolume_distance(type = "centroid")

# Minimum distance
hypervolume_distance(type = "minumum")

# Minimum distance
hypervolume_overlap_statistics()

# 
# library(rgl)
# 
# env_points <- as.data.frame(as.points(env_layers_clean))
# 
# plot3d( 
#   x = env_points[,3], 
#   y = env_points[,4], 
#   z = env_points[,5], 
#   type = 's', 
#   radius = 1,
#   alpha = 0.5,
#   xlab = paste0(unlist(stringr::str_split(names(env_points)[3], "_"))[1], "_" , 
#                 unlist(stringr::str_split(names(env_points)[3], "_"))[3]),
#   ylab = paste0(unlist(stringr::str_split(names(env_points)[3], "_"))[1], "_" , 
#                 unlist(stringr::str_split(names(env_points)[4], "_"))[3]), 
#   zlab = paste0(unlist(stringr::str_split(names(env_points)[3], "_"))[1], "_" , 
#                 unlist(stringr::str_split(names(env_points)[5], "_"))[3])
#   )

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


cb1 <- get.checkerboard1(occs, envs.bg, bg, aggregation.factor=30)

evalplot.grps(pts = occs, pts.grp = cb1$occs.grp, envs = envs.bg)

evalplot.grps(pts = bg, pts.grp = cb1$bg.grp, envs = envs.bg)






