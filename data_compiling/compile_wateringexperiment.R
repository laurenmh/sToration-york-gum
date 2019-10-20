library(tidyverse)

### Read in data from the watering experiment (Wainwright dataset)
water_abiotic <- read_csv(paste(datpath, "Watering_experiment/Abiotic_Wainwright_JEcol_180726.csv", sep="")) 
water_indiv <- read_csv(paste(datpath, "Watering_experiment/Individuals_Wainwright_JEcol_180725.csv", sep="")) 
water_neighbor <- read_csv(paste(datpath, "Watering_experiment/Neighbourhoods_Wainwright_JEcol_180725.csv", sep="")) 
water_spp <- read_csv(paste(datpath, "Watering_experiment/Species_attributes_Wainwright_JEcol_180725.csv", sep="")) 

