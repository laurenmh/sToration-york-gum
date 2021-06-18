rm(list=ls())
library(ggplot2)
library(ggridges)
library(cowplot)
library(rstan)
library(coda)

# Fig 2 ######

# arca Bendering 
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_arca_bd.rds")
dataBE <- dataBE[which(dataBE$comp==1),] #observed, double nat, remove all comp
dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]

dataB2E <- readRDS(file = "Sim data/shade&comp_gradient_arca_bd.rds")
dataB2E <- dataB2E[which(dataB2E$comp==1),] #observed, double nat, remove all comp
dataB2E <- dataB2E[which(dataB2E$e==3 | dataB2E$e==4 | dataB2E$e==5 | dataB2E$e==6 | dataB2E$e==7),]

# arca perenjori 
dataPE <- readRDS(file = "Sim data/phos&comp_gradient_arca_pj.rds")
dataPE <- dataPE[which(dataPE$comp==1),] #observed, double nat, remove all comp
dataPE <- dataPE[which(dataPE$e==3 | dataPE$e==4 | dataPE$e==5 | dataPE$e==6 | dataPE$e==7),]

dataP2E <- readRDS(file = "Sim data/shade&comp_gradient_arca_pj.rds")
dataP2E <- dataP2E[which(dataP2E$comp==1),] #observed, double nat, remove all comp
dataP2E <- dataP2E[which(dataP2E$e==3 | dataP2E$e==4 | dataP2E$e==5 | dataP2E$e==6 | dataP2E$e==7),]

arca.phos <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Bendering")) +
  stat_density_ridges(data=dataPE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Perenjori")) +
  scale_fill_manual(values = alpha(c("pink", "orange"), .5), name = "Site", labels = c("Bendering", "Perenjori")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

arca.shade <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Bendering")) +
  stat_density_ridges(data=dataP2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Perenjori")) +
  scale_fill_manual(values = alpha(c("pink", "orange"), .5), name = "Site", labels = c("Bendering", "Perenjori")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("LDGR")


# waac Bendering 
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_waac_bd.rds")
dataBE <- dataBE[which(dataBE$comp==1),]
#dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]

dataB2E <- readRDS(file = "Sim data/shade&comp_waac_bd_new.rds")
dataB2E <- dataB2E[which(dataB2E$comp==1),]
#dataB2E <- dataB2E[which(dataB2E$e==3 | dataB2E$e==4 | dataB2E$e==5 | dataB2E$e==6 | dataB2E$e==7),]

# waac perenjori 
dataPE <- readRDS(file = "Sim data/phos&comp_gradient_waac_pj.rds")
dataPE <- dataPE[which(dataPE$comp==1),]
#dataPE <- dataPE[which(dataPE$e==3 | dataPE$e==4 | dataPE$e==5 | dataPE$e==6 | dataPE$e==7),]

dataP2E <- readRDS(file = "Sim data/shade&comp_gradient_waac_pj.rds")
dataP2E <- dataP2E[which(dataP2E$comp==1),]
#dataP2E <- dataP2E[which(dataP2E$e==3 | dataP2E$e==4 | dataP2E$e==5 | dataP2E$e==6 | dataP2E$e==7),]

waac.phos <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Bendering")) +
  stat_density_ridges(data=dataPE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Perenjori")) +
  scale_fill_manual(values = alpha(c("pink", "orange"), .5), name = "Site", labels = c("Bendering", "Perenjori")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

waac.shade <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Bendering")) +
  stat_density_ridges(data=dataP2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "Perenjori")) +
  scale_fill_manual(values = alpha(c("pink", "orange"), .5), name = "Site", labels = c("Bendering", "Perenjori")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("LDGR")

legend <- get_legend(
  # create some space to the left of the legend
  waac.phos + theme(legend.box.margin = margin(0, 0, 0, 6))
)


legend <- get_legend(
  # create some space to the left of the legend
  waac.shade + theme(legend.box.margin = margin(0, 0, 0, 6))
)


a <- plot_grid(waac.phos + theme(legend.position="none"),
               arca.phos + theme(legend.position="none"),
               legend,
               waac.shade + theme(legend.position="none"),
               arca.shade + theme(legend.position="none"),
               labels = c("a) Native", "b) Exotic", "", "c) Native", "d) Exotic"),
               nrow = 2, ncol=3,  rel_widths = c(2,2,.7))


pdf("Figures/Fig2.pdf", width = 10, height = 7)
a
dev.off()









# Fig 3 ######

# arca Bendering 
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_arca_bd.rds")
dataBE <- dataBE[which(dataBE$comp==1 | dataBE$comp==2 |dataBE$comp==5),] #observed, double nat, remove all comp
dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]

# test <- dataBE[which(dataBE$comp==5 & dataBE$e==7),]$ldgr
# auc <- ecdf(test[which(test<HPDinterval((as.mcmc(as.numeric(test))))[,2])])
# 1-auc(0) # percent above zero (prob of coexistence)

dataB2E <- readRDS(file = "Sim data/shade&comp_gradient_arca_bd.rds")
dataB2E <- dataB2E[which(dataB2E$comp==1 | dataB2E$comp==2 |dataB2E$comp==5),] #observed, double nat, remove all comp
dataB2E <- dataB2E[which(dataB2E$e==3 | dataB2E$e==4 | dataB2E$e==5 | dataB2E$e==6 | dataB2E$e==7),]

arca.phos <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Neighbours", labels = c("observed", "double nats", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

arca.shade <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Neighbours", labels = c("observed", "double nats", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("LDGR")

# waac Bendering 
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_waac_bd.rds")
dataBE <- dataBE[which(dataBE$comp==1 | dataBE$comp==2 | dataBE$comp==4),]
#dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]

dataB2E <- readRDS(file = "Sim data/shade&comp_waac_bd_new.rds")
dataB2E <- dataB2E[which(dataB2E$comp==1 | dataB2E$comp==2 | dataB2E$comp==4),]
#dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]
# isolate a distribution to find AUC for 

waac.phos <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Neighbours", labels = c("observed", "double natives", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("sd-2", "sd-1", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

waac.shade <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Neighbours", labels = c("observed", "double natives", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  waac.phos + theme(legend.box.margin = margin(0, 0, 0, 6))
)


a <- plot_grid(waac.phos + theme(legend.position="none"),
               arca.phos + theme(legend.position="none"),
               legend,
               waac.shade + theme(legend.position="none"),
               arca.shade + theme(legend.position="none"),
               labels = c("a) Native", "b) Exotic", "", "c) Native", "d) Exotic"),
               nrow = 2, ncol=3,  rel_widths = c(2,2,.7))


pdf("Figures/Fig3BE.pdf", width = 10, height = 7)
a
dev.off()

# Fig 4 ######

# Perenjori plots 

# arca perenjori 
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_arca_pj.rds")
dataBE <- dataBE[which(dataBE$comp==1 | dataBE$comp==2 |dataBE$comp==5),] #observed, double nat, remove all comp
dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]

dataB2E <- readRDS(file = "Sim data/shade&comp_gradient_arca_pj.rds")
dataB2E <- dataB2E[which(dataB2E$comp==1 | dataB2E$comp==2 |dataB2E$comp==5),] #observed, double nat, remove all comp
dataB2E <- dataB2E[which(dataB2E$e==3 | dataB2E$e==4 | dataB2E$e==5 | dataB2E$e==6 | dataB2E$e==7),]

test <- dataB2E[which(dataB2E$comp==5 & dataB2E$e==7),]$ldgr
auc <- ecdf(test[which(test<HPDinterval((as.mcmc(as.numeric(test))))[,2])])
1-auc(0) # percent above zero (prob of coexistence)

arca.phos <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Neighbours", labels = c("observed", "double natives", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("sd-2", "sd-1", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

arca.shade <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Neighbours", labels = c("observed","double natives", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("LDGR")

# waac perenjori 
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_waac_pj.rds")
dataBE <- dataBE[which(dataBE$comp==1 | dataBE$comp==2 | dataBE$comp==5),]
#dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]

dataB2E <- readRDS(file = "Sim data/shade&comp_gradient_waac_pj.rds")
dataB2E <- dataB2E[which(dataB2E$comp==1 | dataB2E$comp==2 | dataB2E$comp==4),]
#dataBE <- dataBE[which(dataBE$e==3 | dataBE$e==4 | dataBE$e==5 | dataBE$e==6 | dataBE$e==7),]


waac.phos <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Comp", labels = c("observed", "double natives", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

waac.shade <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("gold", "skyblue", "pink"), .5), name = "Comp", labels = c("observed", "double natives", "none")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  waac.shade + theme(legend.box.margin = margin(0, 0, 0, 6))
)


a <- plot_grid(waac.phos + theme(legend.position="none"),
               arca.phos + theme(legend.position="none"),
               legend,
               waac.shade + theme(legend.position="none"),
               arca.shade + theme(legend.position="none"),
               labels = c("a) Native", "b) Exotic", "", "c) Native", "d) Exotic"),
               nrow = 2, ncol=3,  rel_widths = c(2,2,.7))


pdf("Figures/Fig4PJ.pdf", width = 10, height = 7)
a
dev.off()

