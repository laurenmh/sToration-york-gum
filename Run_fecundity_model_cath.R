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
SpDataFocal <- subset(SpData, SpData$Focal.sp.x=="A")
NeededNames <- colnames(SpDataFocal)[c(10:55,58)] 
ModData <- subset(SpData, select = NeededNames)
ModData <- na.omit(ModData)

# Determine the number of data points
N <- as.integer(dim(ModData)[1])
Fecundity <- as.integer(ModData$Number.flowers.total)
Intra <- as.integer(ModData$Arctotheca.calendula)
Other <- rep(0, N)
Threshold <- 13

# Create the model matrix
SpMatrixOriginal <- subset(ModData, select = setdiff(colnames(ModData), c("Number.flowers.total", "Arctotheca.calendula")))
SpMatrixOriginal <- as.matrix(SpMatrixOriginal)
SpTotals <- colSums(SpMatrixOriginal)
S <- sum(SpTotals > Threshold)
SpMatrix <- matrix(data = NA, nrow = N, ncol = S)
k <- 1
for(s in 1:ncol(SpMatrixOriginal)){
  if(SpTotals[s] > Threshold){
    SpMatrix[,k] <- SpMatrixOriginal[,s] 
    k <- k + 1
  } else{
    Other <- Other + SpMatrixOriginal[,s]
  }
}

# Do a preliminary fit to compile the stan model and check for convergence, mixing,
#    autocorrelation, etc.
PrelimFit <- stan(file = "fecundity_model_cath.stan", data = c("N", "S", "Fecundity", "Intra", "SpMatrix", "Other"),
                  iter = 6000, chains = 3, thin = 1, control = list(adapt_delta = 0.8, max_treedepth = 10))
PrelimFit
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon", pars = c("alpha_intra", "alpha_mean", "alpha_sp"))

# Diagnostic plots
pairs(PrelimFit, pars = c("lambda", "alpha_intra", "alpha_mean"))
traceplot(PrelimFit)
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
inclusion <- NULL
for(i in 1:S){
  if(alpha_effects[i,4] < 0 & alpha_effects[i,5] > 0){
    Other <- Other + SpMatrix[,i]
  }else{
    inclusion <- c(inclusion, i)
  }
}
SpMatrix <- SpMatrix[,inclusion]
S <- ncol(SpMatrix)

FinalFit <- stan(file = "GrandMeanInteraction.stan", data = c("N", "S", "Fecundity", "Intra", "SpMatrix", "Other"),
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



