## Make sure to manually set your data pathway in set_pathway.R first!!

source("data_compiling/compile_observational_list.R")

## Note: The observational data are in a list form, with a separate list for each of six focal species
## These are: Aira.caryophylla", "Hypochaeris.glabra", "Podotheca.gnaphalioides", "Trachymene.ornata", "Ursinia.anthemoides", "Waitzia.acuminata"  
aira_obvs_fecund <- obvs_fecund[[1]]
hypochaeris_obvs_fecund <- obvs_fecund[[2]]
podotheca_obvs_fecund <- obvs_fecund[[3]]
trachymene_obvs_fecund <- obvs_fecund[[4]]
ursinia_obvs_fecund <- obvs_fecund[[5]]
waitzia_obvs_fecund <- obvs_fecund[[6]]

