#### Running Fecundity model 
# This script comes after data_wrangling_cath.R and pairs with stan file fecundity_model_cath.stan 

# From Topher 
# Trial run of the "grand mean" model of species interactions. 
#    NOTE: We should come up with a better name at some point

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in the data and sort it
SpData <- read.csv("water_full_env.csv")
# We should think carefully about this step. na.omit() will omit any row where any of the data values are NA,
#   but we might now want that (e.g. if the NA is in a column we don't care about). 
SpData <- na.omit(SpData)     
SpDataFocal <- subset(SpData, Focal.sp.x == "W")
SpDataFocal <- SpDataFocal %>% mutate_at(c("Canopy", "Colwell.P"), ~(scale(.) %>% as.vector)) # not working?

# From here we need to calculate and create the Intra and Other vectors, the
#   SpMatrix of abundances for all species we are including, and other necessary
#   objects like N, S, shade, phos, and Fecundity
Intra <- as.integer(SpDataFocal$Waitzia.acuminata)
N <- as.integer(dim(SpDataFocal)[1])
Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
shade <- SpDataFocal$Canopy
phos <- SpDataFocal$Colwell.P

# This chunk of the code will implement your solution for tracking species identity
Species <- names(SpDataFocal[10:69])
Species.ID <- as.data.frame(Species)
Species.ID$Included <- rep(1,length(Species))
Species.ID$Column <- rep(NA, length(Species)) # The data frame is a great idea. We can also include a column for index value so we can easily

# Now, go through the interspecific abundances to initially filter out any that have too low of an abundance
Threshold <- 14
Other <- rep(0, N)
SpMatrixOriginal <- subset(SpDataFocal, select = Species)
SpTotals <- colSums(SpMatrixOriginal)
S <- sum(SpTotals > Threshold)
SpMatrix <- matrix(data = NA, nrow = N, ncol = S) 
k <- 1
for(s in 1:length(Species)){
  if(SpTotals[s] > Threshold){
    SpMatrix[,k] <- SpMatrixOriginal[,s] 
    Species.ID$Column[s] <- k # Track which column in SpMatrix corresponds to each original species
    k <- k + 1
  } else{
    Species.ID$Included[s] <- 0 # If the species does not meet the abundance threshold, then track that they are not included
    Other <- Other + SpMatrixOriginal[,s]
  }
}

# If we've processed the data correctly so far, all three of these should be the same
S
max(Species.ID$Column, na.rm = TRUE)
sum(Species.ID$Included)


# Do a preliminary fit to compile the stan model and check for convergence, mixing,
#    autocorrelation, etc.
PrelimFit <- stan(file = "fecundity_model_env_cath.stan", data = c("N", "S", "Fecundity", "Intra", "SpMatrix", "Other", "shade", "phos"),
                  iter = 3000, chains = 3, thin = 1, control = list(adapt_delta = 0.8, max_treedepth = 10))
# Things I'm trying 
# increasing iterations: 
# increasing thinning: 
# adapt_delta:



PrelimFit
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon", pars = c("alpha_sp"))
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon", pars = c("lambda", "b", "c"))

# Diagnostic plots
pairs(PrelimFit, pars = c("lambda", "b", "c"))
traceplot(PrelimFit, pars = c("lambda", "b", "c"))
# autocorrelation of the MCMC samples
quartz()
par(mfcol = c(2,2))
# Now extract the posteriors
posteriors <- extract(PrelimFit)
# Now plot each posterior
acf(posteriors$lambda)
acf(posteriors$alpha_intra)
acf(posteriors$alpha_mean)

# Once the chains have converged with no autocorrelation, sort through which
#    of the alpha_sp 95% credible intervals overlap 0 (i.e. no effect).
alpha_effects <- summary(PrelimFit, pars = "alpha_sp", probs = c(0.025, 0.975))$summary

for(i in 1:length(Species)){
  if(Species.ID$Included[i] == 1){ # If it was included in the preliminary model fit
    if(alpha_effects[Species.ID$Column[i],4] < 0 & alpha_effects[Species.ID$Column[i],5] > 0){ # Check if the 95% CI overlaps 0
      Species.ID$Included[i] <- 0 # If so, change this to 0 so it is not included in the final fit
      Other <- Other + SpMatrix[,Species.ID$Column[i]] # Update the Other vector with that species' abundance
    } # If the 95% CI doesn't overlap 0, then don't do anything and leave the Included entry as 1
  }
}

S <- sum(Species.ID$Included)
SpMatrixNew <- matrix(data = NA, nrow = N, ncol = S)
k <- 1
for(i in 1:length(Species)){
  if(Species.ID$Included == 1){
    SpMatrixNew[,k] <- SpMatrix[,Species.ID$Column[i]]
    Species.ID$Column[i] <- k
    k <- k + 1
  }
}
SpMatrix <- SpMatrixNew

FinalFit <- stan(file = "fecundity_model_env_cath.stan", data = c("N", "S", "Fecundity", "Intra", "SpMatrix", "Other", "shade", "phos"),
                 iter = 6000, chains = 3, thin = 1, control = list(adapt_delta = 0.8, max_treedepth = 10))

plot(FinalFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon", pars = c("alpha_intra", "alpha_mean", "alpha_sp"))
FinalPosteriors <- extract(FinalFit)

# Check on the realized alphas for the remaining 5 species
AlphaPosteriors <- vector(mode = "list", length = 5)
for(i in 1:5){
  AlphaPosteriors[[i]] <- FinalPosteriors$alpha_mean + FinalPosteriors$alpha_sp[,i]
}
quartz()
par(mfrow = c(2,3))
for(i in 1:5){
  hist(AlphaPosteriors[[i]], main = paste("Species ", i, sep = ""), col = "salmon")
}





