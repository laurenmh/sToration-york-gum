# Load in the full data set and set working directories
library(here)
YorkGum <- read.csv(here("water_full_env.csv"))
FocalSp <- "W" # or H, T, A

#subset the data and create quatratic transformations of the environmental covariates
SpData <- subset(YorkGum, Focal.sp.x == FocalSp & Trt.comp.x == "S")
SpData$logLambda <- log(SpData$Number.flowers.total)
SpData$Phos2 <- SpData$Colwell.P ^ 2
SpData$Canopy2 <- SpData$Canopy ^ 2

# Fit the competing models with different quadratic options
LinearMod <- lm(logLambda ~ Colwell.P + Canopy + Colwell.P*Canopy, data = SpData)
QuadraticP <- lm(logLambda ~ Phos2 + Canopy + Phos2*Canopy, data = SpData)
QuadraticC <- lm(logLambda ~ Colwell.P + Canopy2 + Colwell.P*Canopy2, data = SpData)
QuadraticB <- lm(logLambda ~ Phos2 + Canopy2 + Phos2*Canopy2, data = SpData)

# Make predictions from each of the models
PredData <- subset(SpData, select = c("Colwell.P", "Canopy", "logLambda", "Phos2", "Canopy2"))
LinearPred <- cbind(PredData, predict(LinearMod, PredData, interval = "confidence"))
QuadP_pred <- cbind(PredData, predict(QuadraticP, PredData, interval = "confidence"))
QuadC_pred <- cbind(PredData, predict(QuadraticC, PredData, interval = "confidence"))
QuadB_pred <- cbind(PredData, predict(QuadraticB, PredData, interval = "confidence"))

# Make a plot with the residuals for each model across phosphorous and canopy cover
FigName <- paste("LambdaFunctionalForms/", FocalSp, "_residuals.pdf", sep = "")
pdf(file = FigName, width = 10, height = 6, onefile = FALSE, paper = "special")
  par(mfrow = c(2,4))
  # Phosphorous
  plot(x = LinearPred$Colwell.P, y = LinearPred$fit - LinearPred$logLambda, 
       xlab = "Phosphorous", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  mtext("Monotonic", side = 3, line = 0.5)
  
  plot(x = QuadP_pred$Colwell.P, y = QuadP_pred$fit - QuadP_pred$logLambda, 
       xlab = "Phosphorous", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  mtext("Phosphorous^2", side = 3, line = 0.5)
  
  plot(x = QuadC_pred$Colwell.P, y = QuadC_pred$fit - QuadC_pred$logLambda, 
       xlab = "Phosphorous", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  mtext("Canopy^2", side = 3, line = 0.5)
  
  plot(x = QuadB_pred$Colwell.P, y = QuadB_pred$fit - QuadB_pred$logLambda, 
       xlab = "Phosphorous", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  mtext("Both^2", side = 3, line = 0.5)
  
  # Canopy cover
  plot(x = LinearPred$Canopy, y = LinearPred$fit - LinearPred$logLambda, 
       xlab = "Canopy", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  
  plot(x = QuadP_pred$Canopy, y = QuadP_pred$fit - QuadP_pred$logLambda, 
       xlab = "Canopy", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  
  plot(x = QuadC_pred$Canopy, y = QuadC_pred$fit - QuadC_pred$logLambda, 
       xlab = "Canopy", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
  
  plot(x = QuadB_pred$Canopy, y = QuadB_pred$fit - QuadB_pred$logLambda, 
       xlab = "Canopy", ylab = "Residuals", main = "", ylim = c(-2.5, 2.5))
  abline(h = 0, lty = 2)
dev.off()




