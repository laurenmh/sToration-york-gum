rm(list=ls())
library(ggplot2)
library(ggridges)
library(cowplot)

# ARCA PHOS --------------------------------------------------------

dataB <- readRDS(file = "Sim data/phos&comp_gradient_arca_bd.rds")
dataP <- readRDS(file = "Sim data/phos&comp_gradient_arca_pj.rds")
dataB2 <- readRDS(file = "Sim data/shade&comp_gradient_arca_bd.rds")
dataP2 <- readRDS(file = "Sim data/shade&comp_gradient_arca_pj.rds")

B <- ggplot(NULL) +
  stat_density_ridges(data=dataB, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue", "blue"), .3), name = "Phos", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Competitor manipulation") +
  #scale_y_discrete(labels =c("Schoenus*3", "Schoenus*2", "WH.75", "WH.5", "observed", "double natives", "halve weeds", "remove weeds"))+
  theme_ridges(center = T) +
  xlab("")

P <- ggplot(NULL) +
  stat_density_ridges(data=dataP, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue", "blue"), .3), name = "Phos", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("") +
  theme_ridges(center = T)+
  xlab("") 

B2 <- ggplot(NULL) +
  stat_density_ridges(data=dataB2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("pink", "red"), .3), name = "Shade", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  ylab("Competitor manipulation") +
  xlab("LDGR")
  
P2 <- ggplot(NULL) +
  stat_density_ridges(data=dataP2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("pink", "red"), .3), name = "Shade", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("") +
  theme_ridges(center = T) +
  xlab("LDGR")

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).

legend <- get_legend(
  # create some space to the left of the legend
  B + theme(legend.box.margin = margin(0, 0, 0, 3))
)
legend2 <- get_legend(
  # create some space to the left of the legend
  B2 + theme(legend.box.margin = margin(0, 0, 0, 3))
)

a <- plot_grid(B + theme(legend.position="none"),
          P + theme(legend.position="none"),
          legend, 
          B2 + theme(legend.position="none"),
          P2 + theme(legend.position="none"),
          legend2, 
          labels = c("A", "B", "", "C", "D", ""),
          nrow = 2, ncol=3,  rel_widths = c(1,1,.2))


pdf("Figures/arca_all.pdf", width = 10, height = 7)
a
dev.off()

# WAAC #--------------------------------------------------------------
dataB <- readRDS(file = "Sim data/phos&comp_gradient_waac_bd.rds")
dataP <- readRDS(file = "Sim data/phos&comp_gradient_waac_pj.rds")
dataB2 <- readRDS(file = "Sim data/shade&comp_waac_bd_new.rds")
dataP2 <- readRDS(file = "Sim data/shade&comp_gradient_waac_pj.rds")

B <- ggplot(NULL) +
  stat_density_ridges(data=dataB, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue", "blue"), .3), name = "Phos", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Competitor manipulation") +
  #scale_y_discrete(labels =c("Schoenus*3", "Schoenus*2", "WH.75", "WH.5", "observed", "double natives", "halve weeds", "remove weeds"))+
  theme_ridges(center = T) +
  xlab("")

P <- ggplot(NULL) +
  stat_density_ridges(data=dataP, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("skyblue", "blue"), .3), name = "Phos", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("") +
  theme_ridges(center = T)+
  xlab("") 

B2 <- ggplot(NULL) +
  stat_density_ridges(data=dataB2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("pink", "red"), .3), name = "Shade", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ridges(center = T) +
  ylab("Competitor manipulation") +
  xlab("LDGR")

P2 <- ggplot(NULL) +
  stat_density_ridges(data=dataP2, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = comp, fill = e)) +
  scale_fill_manual(values = alpha(c("pink", "red"), .2), name = "Shade", labels = c("low", "high")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("") +
  theme_ridges(center = T) +
  xlab("LDGR")


legend <- get_legend(
  # create some space to the left of the legend
  B + theme(legend.box.margin = margin(0, 0, 0, 6))
)
legend2 <- get_legend(
  # create some space to the left of the legend
  B2 + theme(legend.box.margin = margin(0, 0, 0, 6))
)

a <- plot_grid(B + theme(legend.position="none"),
               P + theme(legend.position="none"),
               legend, 
               B2 + theme(legend.position="none"),
               P2 + theme(legend.position="none"),
               legend2, 
               labels = c("A", "B", "", "C", "D", ""),
               nrow = 2, ncol=3,  rel_widths = c(1,1,.5))


pdf("Figures/waac_all.pdf", width = 10, height = 7)
a
dev.off()



