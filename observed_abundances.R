dat <- read.csv("water_full_env.csv")

dat.sub <- subset(dat, select=c("Waitzia.acuminata", "Reserve.x",
                                "Block.x", "Plot.ID.x", "Quadrant.x",
                                "Plot.ID.quadrant"))
max(dat$Waitzia.acuminata)
 
                  