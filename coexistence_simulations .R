library(rstan)
library(coda)


# load in the rdata posteriors 
load("Topher model fits/ARCA_Phos_Finalfit.rdata")
arca.phos <- rstan::extract(FinalFit)
remove(FinalFit)

load("Topher model fits/ARCA_Shade_Finalfit.rdata")
arca.shade <- rstan::extract(FinalFit)
remove(FinalFit)

load("Topher model fits/Waitzia_Phos_Finalfit.rdata")
waac.phos <- rstan::extract(FinalFit)
remove(FinalFit)

load("Topher model fits/Waitzia_Shade_Finalfit.rdata")
waac.shade <- rstan::extract(FinalFit)
remove(FinalFit)


# load in s and g data 
load("SurvivalAndGermination/Germination.rdata")
germ <- rstan::extract(PrelimFit)
germ<-as.data.frame(germ)
remove(PrelimFit)

load("SurvivalAndGermination/Survival.rdata")
surv <- rstan::extract(PrelimFit)
surv<-as.data.frame(surv)
remove(PrelimFit)

germination <- data.frame(arca=germ$p.1, waac=germ$p.2)

survival <- data.frame(arca=surv$p.1, waac=surv$p.2)


# tryng to figure out the other values I need 
SpData <- read.csv("water_full_env.csv")
SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
SpData <- na.omit(SpData) 
SpDataFocal <- subset(SpData, Focal.sp.x == "W")

# From here we need to calculate and create the Intra vector, the
#   SpMatrix of abundances for all heterospecifics, and other necessary
#   objects like N, S, and Fecundity
shade <- as.vector(scale(SpDataFocal$Canopy))
phos <- as.vector(scale(SpDataFocal$Colwell.P))
Intra <- as.integer(SpDataFocal$Waitzia.acuminata)
N <- as.integer(dim(SpDataFocal)[1])
Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
Species <- names(SpDataFocal[10:69]) 
Species <- setdiff(Species, "Waitzia.acuminata")
TempSpMatrix <- subset(SpDataFocal, select = Species)
# Now discount any columns with 0 abundance
SpTotals <- colSums(TempSpMatrix)
SpToKeep <- SpTotals > 0
SpMatrix <- matrix(NA, nrow = nrow(TempSpMatrix), ncol = sum(SpToKeep)+1)
SpMatrix[,1] <- Intra
s <- 2
for(i in 1:ncol(TempSpMatrix)){
  if(SpToKeep[i]){
    SpMatrix[,s] <- TempSpMatrix[,i]
    s <- s + 1
  }
}
print(s)
S <- ncol(SpMatrix)
print(S)



# coexistence simulations - invasion growth rates -----------------------------------------

invader.abund <- 1

# things to set prior for hard coding 

germ=germination$waac
surv=survival$waac
lambdas=waac.phos$lambdas

dat <- read.csv("water_full_env.csv")

# get observed population sizes within rings (subset differently for community you want to invade into)
dat.sub <- subset(dat, select=c("Arctotheca.calendula", "Reserve.x",
                                "Block.x", "Plot.ID.x", "Quadrant.x",
                                "Plot.ID.quadrant"))
arca.abundance <- dat$Arctotheca.calendula[which(dat$Arctotheca.calendula>0)]
# invasion equation (maybe won't be useful because will need to hardcode the neighbour/community options)
# do.new.seeds <- function(lived, lambda, alpha_) {
#   N_tp1 <- lived*lambda/(1+alpha_intra*lived)
#   return(N_tp1)
# }

# run the simulation 
posteriors=100 # how many values from the posterior we want to sample?
runs=length(arca.abundance) # run through each option for resident community population size 
waac.into.generic <- matrix(NA, nrow=runs, ncol=posteriors) # make an empty matrix (array for env??)
#e = # to run through increments for each environmental condition
for (r in 1:runs) {
  for (p in 1:posteriors) {
    #for (e in 1:)
    x <- sample(seq(1, 9000),1)
    y <- sample(seq(1,4500),1)
    
    # calculate resident abundance
    #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
    Nj <- arca.abundance[r]
    germ_j <- germination$waac[y]
    lived.generic <- germ_j*Nj
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*invader.abund + 
      invader.abund*germ[y]*lambdas[x]/(1+waac.phos$alpha_hat_eij_tilde[x]*lived.generic)
    # calculate LDGR of WAAC
    waac.into.generic[r,p] <- log(waac_tp1/invader.abund)
    
  }
  print(r)
}

# QUESTIONS 
# how do I get alpha_eji for the specific species effects?
# environmental effects?
# how can i sort by the resident community size for the heat map 
hist(waac.into.generic)


# trying to figure out the environmental factor 
plot(density(waac.phos$alpha_hat_eij_tilde)) # use values above 1? 
plot(density(waac.phos$alpha_hat_ij_tilde))


