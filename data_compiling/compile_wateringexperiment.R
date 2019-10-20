library(tidyverse)

## Questions for Margie (or a closer read of the paper):
# 1) Can we confirm the nesting design (my understanding is: Reserve > Block > Plot.ID > Quadrant)
# 2) For water_indiv, what is the Plot column? (it's an integer, doesn't correspond to either Plot.ID or Quadrant)
# 3) What are the focal species? (A, H, T, W)

### Read in data from the watering experiment (Wainwright dataset)

# Plot level environmental data
water_abiotic <- read_csv(paste(datpath, "Watering_experiment/Abiotic_Wainwright_JEcol_180726.csv", sep=""))

# Quadrant-level focal individual fecundity data
water_indiv <- read_csv(paste(datpath, "Watering_experiment/Individuals_Wainwright_JEcol_180725.csv", sep="")) %>%
  mutate(Plot.ID = paste(Block, Trt.comp, Trt.water, sep = ""),
         Plot.ID.quadrant = paste(Plot.ID, Quadrant, sep=".")) %>%
  select(-Block.comp.water, -Block.water) %>%
  select(Reserve, Block, Plot.ID, Plot, Quadrant, Plot.ID.quadrant, Trt.comp, Trt.water, Focal.sp, Number.flowers.total, Seedcount.extrapolated.integer)

# Quadrant-level neighborhood species density
water_neighbor <- read_csv(paste(datpath, "Watering_experiment/Neighbourhoods_Wainwright_JEcol_180725.csv", sep="")) %>%
  mutate(Focal.sp = toupper(focal.sp), 
         Block = block, 
         Quadrant = quadrant,
         Reserve = reserve,
         Neighbor.sp = neighbour.id,
         Number = number) %>%
  select(Reserve, Block, Plot.ID, Quadrant, Plot.ID.quadrant, Trt.comp, Trt.water, Focal.sp, Neighbor.sp, Number)

# Species-level trait data
water_spp <- read_csv(paste(datpath, "Watering_experiment/Species_attributes_Wainwright_JEcol_180725.csv", sep="")) 

