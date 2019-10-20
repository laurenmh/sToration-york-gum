## Make sure to manually set your data pathway in set_pathway.R first!!

library(tidyverse)

## Read in data from the shade experiment 
load(paste(datpath, "Observational_data/complete.Rdata", sep="")) 


## Rename the file to match sToration naming conventions
obvs_fecund <- fecundity.data
rm(fecundity.data)

# Note: The observational data are in a list form, with a separate list for each of six focal species
# These are: Aira.caryophylla", "Hypochaeris.glabra", "Podotheca.gnaphalioides", "Trachymene.ornata", "Ursinia.anthemoides", "Waitzia.acuminata"  
# If you are more comfortable in dataframes, you can subset each species as its own dataframe. Just uncomment and run this code:
# aira_obvs_fecund <- obvs_fecund[[1]]
# hypochaeris_obvs_fecund <- obvs_fecund[[2]]
# podotheca_obvs_fecund <- obvs_fecund[[3]]
# trachymene_obvs_fecund <- obvs_fecund[[4]]
# ursinia_obvs_fecund <- obvs_fecund[[5]]
# waitzia_obvs_fecund <- obvs_fecund[[6]]
