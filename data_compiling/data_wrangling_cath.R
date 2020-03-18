# THE WATER EXPERIMENT has both shade and soil P data and focal plant and neighbour info
# but with a slightly different spatial design (25 by 25 quadrant in 50 by 50 plots, within block
# at both Perejnori (north - dry) and Bendering (south - wet))  and flowercount will need 
# to be extrapolated to seedcount - Cath 


# SO doing some data wrangling to match the watering datatsets together in a nice format for modelling - Cath
# TO RUN the following code, make sure compile_wateringexperiment.R has been run already 
# water_neighbour needs to be made into wide form 
water_spread <- water_neighbor
water_spread <- water_spread %>% 
  group_by_at(vars(-Number)) %>% # group by everything other than the number of neighbours
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index column. 
  spread(Neighbor.sp, Number, fill = 0) %>% 
  select(-row_id)  # drop the index

# # match the water_spread dataset with water_indiv (note the number of unique rows (individuals))
# dim(water_spread)
# dim(water_indiv)

# water_spread$Plot.ID[!(water_spread$Plot.ID %in% water_indiv$Plot.ID)]

water_full <- inner_join(water_spread, water_indiv, by = "Plot.ID.quadrant") 
dim(water_full)
NeededNames <- colnames(water_full)[c(1:69, 78:79)] 
water_full_mod <- subset(water_full, select = NeededNames)

# need to join with environmental variables shade (canopy cover) and Colwell P
# abiotic data is at plot level 
water_full_env <- water_full_mod
abiotic <- water_abiotic
abiotic <- rename(abiotic, Plot.ID.x = Plot.ID)
water_full_env <- full_join(water_full_env, abiotic, by = "Plot.ID.x")
water_full_env <- subset(water_full_env, water_full_env$Trt.water.x=="D") # only use dry (control) water treatment

# save as a csv file to be called by Run_fecundity_model_cath.R
write.csv(water_full_env, "water_full_env.csv")

# # make shade into a category with low and high 
# plot(density(water_full_env$Canopy, na.rm = T)) # cutting low at < 0.4 for now 
# water_full_env <- water_full_env %>% mutate(shade=cut(Canopy, breaks=c(0, 0.4, Inf), labels=c("low","high")))
# table(water_full_env$shade)
#








