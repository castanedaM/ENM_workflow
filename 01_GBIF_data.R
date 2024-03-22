# Title: GBIF Occurrence data for Virginia oysters (Crassostrea virginica)
# Author: Mariana Castaneda Guzman
# Date created: 9/5/2023
# Date last updated: 3/08/2023

# Description: GBIF occurrence data for Virginia oysters (Crassostrea virginica)
# in Coastal states in Eastern United States


# Packages ----------------------------------------------------------------

library(rgbif) # for using GBIF data
library(maps) # for simple maps
library(ggplot2) # for plotting
library(RColorBrewer) # for editing plots
library(tidyverse) # for data handling

# Download GBIF data ------------------------------------------------------

va_oyster <- c("Crassostrea virginica")

# Download GBIF occurrence data for this species.
gbif_data <- occ_data(scientificName = va_oyster, 
                      hasCoordinate = TRUE, 
                      limit = 20000)

# check how the data are organized
View(gbif_data$data)
head(gbif_data$data)
nrow(gbif_data$data)
names(gbif_data)
names(gbif_data$meta)
names(gbif_data$data)

table(gbif_data$data$continent)
table(gbif_data$data$country)

# Visualize data ----------------------------------------------------------

# ALL OCCURRENCES

# cv = "Crassostrea virginica"

# get the columns that matter for mapping and cleaning the occurrence data
all_cv_coords <- gbif_data$data[ , c("decimalLongitude", "decimalLatitude", 
                                     "individualCount", "occurrenceStatus", 
                                     "coordinateUncertaintyInMeters", 
                                     "institutionCode", "references")]
head(all_cv_coords)


# map the occurrence data
maps::map("world", 
    xlim = range(all_cv_coords$decimalLongitude), 
    ylim = range(all_cv_coords$decimalLatitude)) 
points(all_cv_coords[ , c("decimalLongitude", "decimalLatitude")],
       pch = ".", cex = 3, col = "red")

# UNITED STATES OCCURRENCES

unique(gbif_data$data$country)

usa_cv <- subset(gbif_data$data, country == "United States of America")

usa_cv_coords <- usa_cv[ , c("decimalLongitude", "decimalLatitude", "individualCount", "occurrenceStatus", "coordinateUncertaintyInMeters", "institutionCode", "references")]
head(usa_cv_coords)

# map the occurrence data:
maps::map("state", 
    xlim = range(usa_cv_coords$decimalLongitude), 
    ylim = range(usa_cv_coords$decimalLatitude))  
points(usa_cv_coords[ , c("decimalLongitude", "decimalLatitude")], 
       pch = ".", cex = 3, col = "red")


# EASTERN UNITED STATES

# Westernmost part of Florida: lat: 30.867218° long: -87.642057°
# Easternmost part of Maine: lat: 44.815109° long: -66.950473°

# Southernmost part of Florida: lat: 25.040514° long:-80.833716°
# Northernmost part of Maine: lat: 47.459574° long: -69.224743°


unique(gbif_data$data$stateProvince)


east_usa_cv <- subset(gbif_data$data, decimalLongitude <= -66.950473 &
                        decimalLongitude >= -87.642057 &
                        decimalLatitude <= 47.459574 &
                        decimalLatitude >= 23.140514)
nrow(east_usa_cv)

                      
east_usa_cv_coords <- east_usa_cv[ , c("decimalLongitude", "decimalLatitude", 
                                       "individualCount", "occurrenceStatus", 
                                       "coordinateUncertaintyInMeters", 
                                       "institutionCode", "references")]


eastern_coastal_states <- c("Connecticut", "Delaware", 
                            "Florida", "Georgia", "Maine",
                            "Maryland", "Massachusetts", 
                            "New Hampshire", "New Jersey",
                            "New York", "North Carolina", 
                            "Rhode Island", "South Carolina", 
                            "Virginia")

# map the occurrence data:
maps::map(database = "state", 
    xlim = range(east_usa_cv_coords$decimalLongitude), 
    ylim = range(east_usa_cv_coords$decimalLatitude))
maps::map(database = "state",
    regions = east,col = "lightgrey",fill=T,add=TRUE)
points(east_usa_cv_coords[ , c("decimalLongitude", "decimalLatitude")],
       pch = 16, cex = log(east_usa_cv_coords$individualCount)/4, col = "red")
title("Virginia Oyster occurrences in Eastern United States")


# Clean Data --------------------------------------------------------------

# View(east_usa_cv)

east_usa_cv_clean <- east_usa_cv[, c(1,3,4,41,42,51,77,78,82,83)]
head(east_usa_cv_clean)


unique(east_usa_cv_clean$waterBody)
# View(east_usa_cv_clean)


east_usa_cv_clean$individualCount <- ifelse(is.na(east_usa_cv_clean$individualCount) |
                                              east_usa_cv_clean$individualCount == 0,
                                            1, east_usa_cv_clean$individualCount)

# View(east_usa_cv_clean)


mapview::mapview(east_usa_cv_clean, 
                 xcol = "decimalLongitude", ycol = "decimalLatitude")

# remove the ones that are way outside the area

rm_lat <- c(45.52583, 45.69736, 43.08194, 26.51667, 31.60127, 35.74117)
rm_long <- c(-79.15667, -73.51454, -75.14806, -77.36667, -86.65358, -79.0193)
rm_stateProvice <- c("Québec", "Bahama Islands", "Pennsylvania", "Alabama")
rm_key <- c("4006766611", "3947135221", "3966692565", "1698584671", "1234582853",
            "1456286540", "477730502", "477830558", "351320843", "476969788",
            "3947804987", "215837455", "1455937102", "477017595", "1698584257", 
            "4110413425", "1698584257", "1698750789")

rm_rows <- subset(east_usa_cv_clean, decimalLatitude %in% rm_lat &
                                decimalLongitude %in% rm_long |
                                stateProvince %in% rm_stateProvice | 
                                key %in% rm_key)$key

east_usa_cv_clean_1 <- east_usa_cv_clean %>% filter(!key %in% rm_rows, !is.na(year))


mapview::mapview(east_usa_cv_clean_1, 
                 xcol = "decimalLongitude", 
                 ycol = "decimalLatitude")


east_usa_cv_clean <- east_usa_cv_clean_1


east_usa_cv_clean <- east_usa_cv_clean %>% 
  mutate(stateProvince = ifelse(stateProvince == "Va", "Virginia", 
                                ifelse(stateProvince == "Md", "Maryland", 
                                       ifelse(stateProvince == "Florida (State)", "Florida",
                                              ifelse(stateProvince == "North carolina", "North Carolina", 
                                                     ifelse(stateProvince == "District of Columbia", "Virginia", 
                                                            ifelse(stateProvince == "South carolina", "South Carolina", stateProvince
                                                            )))))))




saveRDS(object = east_usa_cv_clean, file = "../Data/GBIF/cv_occurrences_clean.RDS")

nrow(east_usa_cv_clean)

# Figures -----------------------------------------------------------------

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


east_usa_cv_clean <- readRDS("../Data/GBIF/cv_occurrences_clean.RDS")

eastern_coastal_states <- c("Connecticut", "Delaware", 
                            "Florida", "Georgia", "Maine",
                            "Maryland", "Massachusetts", 
                            "New Hampshire", "New Jersey",
                            "New York", "North Carolina", 
                            "Rhode Island", "South Carolina", "Virginia")

og_occ_2 <- east_usa_cv_clean %>% 
  filter(year >= 2018, year <= 2023, !is.na(stateProvince)) %>% 
  mutate(year = as.factor(year)) 


p1 <- ggplot(og_occ_2 %>% 
         group_by(stateProvince, year) %>%
         mutate(count = 1) %>% 
         summarise(count = sum(count)), aes(y = reorder(stateProvince, count), x = count)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(aes(fill = year), 
             alpha = 0.75, size = 3, shape = 21, color = "black") +
  scale_fill_brewer(palette =  "RdYlBu") +
  ylab("State") + xlab(expression(italic("C. virginica")~"observations")) +
  scale_x_continuous(breaks = seq(0, 175, 25), limits = c(0, 175)) +
  labs(fill = "Year")
p1

# tiff(filename = paste0("../Figures/cv_obs_per_state.tiff"), 
#      width = 7, height = 5, units = "in", res = 300, compression = "lzw")
# p1
# dev.off()


p2 <- ggplot(data = east_usa_cv_clean %>% 
               filter(year >= 2003, year <= 2023) %>% 
               group_by(year) %>% mutate(count = 1) %>% 
               summarise(count=sum(count)),
             aes(x = year, y = count)) +
  geom_bar(stat= "identity", color = "black", fill = "lightgrey") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1100)) +
  scale_x_continuous(breaks = seq(2003, 2023, 4)) +
  xlab("Year") + ylab(expression("Eastern US"~italic("C. virginica")~"observations")) +
  geom_smooth(method = "gam", color = "firebrick1", se = F) 
p2

tiff(filename = paste0("../Figures/cv_obs_per_year.tiff"), 
     width = 7, height = 5, units = "in", res = 300, compression = "lzw")
p2
dev.off()


states <- map_data("state")
states$east_coast <- ifelse(states$region %in%  tolower(eastern_coastal_states), 1, NA)

p3 <- ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, fill = east_coast, group = group),
               color = "black") + 
  coord_fixed(1.3) +
  guides(fill = FALSE) +
  geom_point(data = east_usa_cv_clean, 
             aes(x = decimalLongitude, y = decimalLatitude), 
             color = "red", size = 0.5) +
  scale_x_continuous(breaks = seq(-150, -40, 10)) +
  xlab("Longitude") + ylab("Latitude")
p3


tiff(filename = paste0("../Figures/cv_obs_map.tiff"), 
     width = 7, height = 5, units = "in", res = 300, compression = "lzw")
p3
dev.off()


p2_3 <- p2 + annotation_custom(ggplotGrob(p3), xmin = 2000, xmax = 2017, 
                       ymin = 450, ymax = 1000)

p2_3

tiff(filename = paste0("../Figures/cv_obs_per_year_&_map.tiff"), 
     width = 7, height = 5, units = "in", res = 300, compression = "lzw")
p2_3
dev.off()


p4 <- ggarrange(ggarrange(p3, p2, labels = c("", "B"), nrow = 1),
         p1, nrow = 2, labels = c("A", "C"))

tiff(filename = paste0("../Figures/cv_obs_summary.tiff"), 
     width = 10, height = 9, units = "in", res = 300, compression = "lzw")
p4
dev.off()

# maps::map(database = "state", 
#           xlim = range(east_usa_cv_clean$decimalLongitude), 
#           ylim = range(east_usa_cv_clean$decimalLatitude))
# maps::map(database = "state",
#           regions = east,col = "lightgrey",fill=T,add=TRUE)
# points(east_usa_cv_clean[ , c("decimalLongitude", "decimalLatitude")],
#        pch = 16, col = east_usa_cv_clean$year,
#        # cex = log(east_usa_cv_clean$individualCount)/4,
#        # col = "red",
#        cex = 0.5)



p5 <- ggplot(data = east_usa_cv_clean %>% 
               group_by(year) %>% 
               filter(year >= 2003) %>% 
               summarise(individualCount_sum = sum(individualCount)),
             aes(x = year, y = individualCount_sum)) +
  geom_bar(stat = "identity") +
  ylab("individual count") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2003, 2023, 2)) +
  theme_classic()
p5


nrow(east_usa_cv_clean)
# table(east_usa_cv_clean$datasetName) 


p5 <- ggplot(og_occ_2 %>% 
               group_by(stateProvince, year) %>%
               summarise(individualCount = sum(individualCount)), 
             aes(y = reorder(stateProvince, individualCount), x = individualCount)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(aes(fill = year), 
             alpha = 0.75, size = 3, shape = 21, color = "black") +
  scale_fill_brewer(palette =  "RdYlBu") +
  ylab("State") + xlab(expression(italic("C. virginica")~"observations")) +
  scale_x_continuous(breaks = seq(0, 350, 50), limits = c(0, 350)) +
  labs(fill = "Year")
p5

tiff(filename = paste0("../Figures/cv_indv_obs_per_state.tiff"), 
     width = 7, height = 5, units = "in", res = 300, compression = "lzw")
p5
dev.off()







