## Make sure to manually set your data pathway in set_pathway.R first!!

library(tidyverse)

## Read in data from the shade experiment 
load(paste(datpath, "Shade_experiment/data.Rdata", sep="")) 

## Rename the file to match sToration naming conventions
shade_fecund <- fecundity.data
rm(fecundity.data)
