rm(list=ls())
library(ggplot2)
library(ggridges)
library(cowplot)

# ARCA --------------------------------------------------------
# Bendering
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_arca_bd.rds")
dataBE <- subset(dataBE, comp == "1")
dataBE <- dataBE[-which(dataBE$e=="3"),]

dataB2E <- readRDS(file = "Sim data/shade&comp_arca_bd.rds")
dataB2E <- subset(dataB2E, comp == "1")
dataB2E <- dataB2E[-which(dataB2E$e=="1"),]


BE <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("#F066EA"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phos") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

B2E <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("orange"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

dataB <- readRDS(file = "Sim data/phos&comp_gradient_arca_bd.rds")
dataB <- subset(dataB, e == "3")
#dataB <- dataB[-which(dataB$comp=="no.comp"),]
dataB2 <- readRDS(file = "Sim data/shade&comp_gradient_arca_bd.rds")
dataB2 <- subset(dataB2, e == "3")
#dataB2 <- dataB2[-which(dataB2$comp=="no.comp"),]


B <- ggplot(NULL) +
  stat_density_ridges(data=dataB, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("#39B600"), .3), name = "Phos", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")

B2 <- ggplot(NULL) +
  stat_density_ridges(data=dataB2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue"), .3), name = "Shade", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  BE + theme(legend.box.margin = margin(0, 0, 0, 6))
)
legend2 <- get_legend(
  # create some space to the left of the legend
  B2 + theme(legend.box.margin = margin(0, 0, 0, 6))
)

a <- plot_grid(BE + theme(legend.position="top"),
               B2E + theme(legend.position="top"),
               B + theme(legend.position="top"),
               B2 + theme(legend.position="top"),
               labels = c("a)", "b)", "", "c)", "d)", ""),
               nrow = 2, ncol=2,  rel_widths = c(1,1))


pdf("Figures/waac_new_bend.pdf", width = 10, height = 7)
a
dev.off()

# PJ
dataPE <- readRDS(file = "Sim data/phos&comp_gradient_arca_pj.rds")
dataPE <- subset(dataPE, comp == "1")
dataP2E <- readRDS(file = "Sim data/shade&comp_gradient_arca_pj.rds")
dataP2E <- subset(dataP2E, comp == "1")

PE <- ggplot(NULL) +
  stat_density_ridges(data=dataPE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("#F066EA"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phos") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

P2E <- ggplot(NULL) +
  stat_density_ridges(data=dataP2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("orange"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

dataP <- readRDS(file = "Sim data/phos&comp_gradient_arca_pj.rds")
dataP <- subset(dataP, e == "3")
dataP2 <- readRDS(file = "Sim data/shade&comp_gradient_arca_pj.rds")
dataP2 <- subset(dataP2, e == "3")


P <- ggplot(NULL) +
  stat_density_ridges(data=dataP, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("#39B600"), .3), name = "Phos", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")

P2 <- ggplot(NULL) +
  stat_density_ridges(data=dataP2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue"), .3), name = "Shade", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  PE + theme(legend.box.margin = margin(0, 0, 0, 6))
)
legend2 <- get_legend(
  # create some space to the left of the legend
  P2 + theme(legend.box.margin = margin(0, 0, 0, 6))
)

a <- plot_grid(PE + theme(legend.position="top"),
               P2E + theme(legend.position="top"),
               P + theme(legend.position="top"),
               P2 + theme(legend.position="top"),
               labels = c("a)", "b)", "", "c)", "d)", ""),
               nrow = 2, ncol=2,  rel_widths = c(1,1))


pdf("Figures/waac_new_pj.pdf", width = 10, height = 7)
a
dev.off()

# WAAC #--------------------------------------------------------------
# Bendering
dataBE <- readRDS(file = "Sim data/phos&comp_gradient_waac_bd.rds")
dataBE <- subset(dataBE, comp == "1")
dataBE <- dataBE[-which(dataBE$e=="3"),]

dataB2E <- readRDS(file = "Sim data/shade&comp_waac_bd_new.rds")
dataB2E <- subset(dataB2E, comp == "1")
dataB2E <- dataB2E[-which(dataB2E$e=="1"),]


BE <- ggplot(NULL) +
  stat_density_ridges(data=dataBE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("#F066EA"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phos") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

B2E <- ggplot(NULL) +
  stat_density_ridges(data=dataB2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("orange"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

dataB <- readRDS(file = "Sim data/phos&comp_gradient_waac_bd.rds")
dataB <- subset(dataB, e == "3")
#dataB <- dataB[-which(dataB$comp=="no.comp"),]
dataB2 <- readRDS(file = "Sim data/shade&comp_waac_bd_new.rds")
dataB2 <- subset(dataB2, e == "3")
#dataB2 <- dataB2[-which(dataB2$comp=="no.comp"),]


B <- ggplot(NULL) +
  stat_density_ridges(data=dataB, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("#39B600"), .3), name = "Phos", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")

B2 <- ggplot(NULL) +
  stat_density_ridges(data=dataB2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue"), .3), name = "Shade", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  BE + theme(legend.box.margin = margin(0, 0, 0, 6))
)
legend2 <- get_legend(
  # create some space to the left of the legend
  B2 + theme(legend.box.margin = margin(0, 0, 0, 6))
)

a <- plot_grid(BE + theme(legend.position="top"),
               B2E + theme(legend.position="top"),
               B + theme(legend.position="top"),
               B2 + theme(legend.position="top"),
               labels = c("a)", "b)", "", "c)", "d)", ""),
               nrow = 2, ncol=2,  rel_widths = c(1,1))


pdf("Figures/waac_new_bend.pdf", width = 10, height = 7)
a
dev.off()

# PJ
dataPE <- readRDS(file = "Sim data/phos&comp_gradient_waac_pj.rds")
dataPE <- subset(dataPE, comp == "1")
dataP2E <- readRDS(file = "Sim data/shade&comp_gradient_waac_pj.rds")
dataP2E <- subset(dataP2E, comp == "1")

PE <- ggplot(NULL) +
  stat_density_ridges(data=dataPE, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("#F066EA"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Phos") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

P2E <- ggplot(NULL) +
  stat_density_ridges(data=dataP2E, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = comp)) +
  scale_fill_manual(values = alpha(c("orange"), .3), name = "Comp", labels = c("observed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  theme_ridges(center = T) +
  xlab("")

dataP <- readRDS(file = "Sim data/phos&comp_gradient_waac_pj.rds")
dataP <- subset(dataP, e == "3")
dataP2 <- readRDS(file = "Sim data/shade&comp_gradient_waac_pj.rds")
dataP2 <- subset(dataP2, e == "3")


P <- ggplot(NULL) +
  stat_density_ridges(data=dataP, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("#39B600"), .3), name = "Phos", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")

P2 <- ggplot(NULL) +
  stat_density_ridges(data=dataP2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue"), .3), name = "Shade", labels = c("ambient")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  #scale_y_discrete(labels =c("sd-1", "sd-2", "sd", "sd+1", "sd+2"))+
  ylab("Comp") +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  PE + theme(legend.box.margin = margin(0, 0, 0, 6))
)
legend2 <- get_legend(
  # create some space to the left of the legend
  P2 + theme(legend.box.margin = margin(0, 0, 0, 6))
)

a <- plot_grid(PE + theme(legend.position="top"),
               P2E + theme(legend.position="top"),
               P + theme(legend.position="top"),
               P2 + theme(legend.position="top"),
               labels = c("a)", "b)", "", "c)", "d)", ""),
               nrow = 2, ncol=2,  rel_widths = c(1,1))


pdf("Figures/waac_new_pj.pdf", width = 10, height = 7)
a
dev.off()
