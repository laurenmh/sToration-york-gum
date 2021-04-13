rm(list=ls())
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)
library("egg")

dataB <- readRDS(file = "Sim data/phos&comp_gradient_waac_bd.rds")
dataP <- readRDS(file = "Sim data/phos&comp_gradient_waac_pj.rds")


B <- ggplot(NULL) +
  stat_density_ridges(data=dataB, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("blue", "orange"), .3), name = "Phos", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Competitor manipulation") +
  scale_y_discrete(labels =c("Schoenus*3", "Schoenus*2", "WH.75", "WH.5", "observed", "double natives", "halve weeds", "remove weeds"))+
  theme_ridges(center = T) +
  ggtitle("Waitzia - Bendering")


P <- ggplot(NULL) +
  stat_density_ridges(data=dataP, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  
  scale_fill_manual(values = alpha(c("blue", "orange"), .3), name = "Phos", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Competitor manipulation") +
  ggtitle("Waitzia - Perenjori") +
  scale_y_discrete(labels =c("Pentaschistis*3", "Pentaschistis*2", "observed", "double natives", "halve weeds", "remove weeds"))+
  theme_ridges(center = T) 


plot <- ggarrange(B,P, ncol=2)
ggsave(filename = file.path("Figures","phos&comp_waac_deviation.pdf"), 
       plot, width = 10, height = 5)
