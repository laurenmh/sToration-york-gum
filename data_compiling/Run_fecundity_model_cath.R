#### Running Fecundity model 
# This script comes after data_wrangling_cath.R and pairs with stan file fecundity_model_cath.stan 

library("rstan")
# This line of code detects the number of cores available to your machine so that
#    chains can be run in parallel. If you will be doing other things on your 
#    computer while the model is running, I suggest not running this line of code
#    as it can slow down other applications on your computer
options(mc.cores = parallel::detectCores())

# This line allows a compiled model to be saved directly, so that stan does
#    not need to recompile it if it hasn't changed which can save time
rstan_options(auto_write = TRUE)

# Load in the data

SpData <- read.csv("water_full_env.csv") 
# subset out the species that you want
SpData <- subset(SpData, focal.sp.x=="A") # something like that

# subset out the fecundity column and neighbour columns
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, exotic_graminoid, exotic_forb, native_forb, native_graminoid))
SpData <- na.omit(SpData)

# make the data you declared in fecundity_model_cath.stan 
N <- as.integer(dim(SpData)[1])
Plot <- as.integer(?)
Block <- as.integer(?)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic_graminoid, exotic_forb, native_forb, native_graminoid))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

# Do a preliminary model fit to compile the stan model and check for autocorrelation
#    and convergence
initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Arca/FecundityModel_random.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
#control = list(adapt_delta = 0.9
PrelimFit

# Take a look at the different convergence of the different chains
traceplot(PrelimFit)

# Check out the posteriors for the different parameters
stan_dens(PrelimFit, pars = "epsilon")
stan_dens(PrelimFit, pars = c("lambda" , "alpha_intra", "alpha_InvGraminoid", "alpha_NatForb", "alpha_NatGraminoid"))

# Plot them together with credible intervals
plot(PrelimFit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon")

# Finally check the autocorrelation of the MCMC samples
quartz()
par(mfcol = c(2,3))
# Now extract the posteriors
posteriors <- extract(PrelimFit)
# Now plot each posterior
acf(posteriors$lambda)
acf(posteriors$alpha_intra)
acf(posteriors$alpha_InvGrass)
acf(posteriors$alpha_InvForb)
acf(posteriors$alpha_NatForb)
acf(posteriors$alpha_Other)

# #acf(posteriors$sigma[,i]
# for(i in 1:6){
#         acf(posteriors$Beta[,i])
#         acf(posteriors$sigma[,i])
# }

# Since it looks like the chains have converged with little autocorrelation,
#    we can generate however many posterior samples we need for our analysis

save(PrelimFit, file = "Arca/Arca_posteriors.rdata")
fit <- extract(PrelimFit)

quartz(width=5, height=5)
plot(density(fit$alpha_intra), col="blue", xlim=c(-.3, .3), ylim=c(0,80), lwd=2, xlab="Alpha Range", ylab="Density", main="Momo")
lines(density(fit$alpha_InvForb), col="purple", lwd=2)
lines(density(fit$alpha_NatForb), col="orange", lwd=2)
lines(density(fit$alpha_InvGrass), col="red", lwd=2)
abline(h=0, col="black")
legend("topleft", c("Intraspecific", "Same Type", "Native Forb", "Invasive Grass"),
       bty="n", col=c("blue", "purple", "orange", "red"), lwd=2)

