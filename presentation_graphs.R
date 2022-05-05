pdf(file = Obs, width = 5, height = 5, onefile = FALSE, paper = "special")

quartz(width=6.7, height=5.6)


ggplot(obs_waac_ldgr)+
  geom_density(aes(x=obs_ldgr, fill="red"))+
#  geom_density(aes(x=nocomp_ldgr, fill="pink"))+
#  geom_density(aes(x=lowp_ldgr, fill="green"))+
#  geom_density(aes(x=lowshade_ldgr, fill="yellow"))+
#  geom_density(aes(x=nocomp_lowp_ldgr, fill="blue")) +
  theme_classic() +
  scale_x_continuous(name="Waitzia Low Density Growth Rate", breaks=c(-0.5, 0, .5, 1, 1.5, 2),
                     labels=c("-0.5", "0", "0.5", "1", "1.5", "2"), lim=c(-0.75, 1.75)) +
  geom_vline(xintercept = 0, color="black", size=1.5) +
  scale_y_discrete(name="Probability Density", lim=c("0","1"), labels=c("", "")) +
  scale_fill_manual(guide='legend', name = 'Scenario',
                    values=alpha(c("red"=cbp2[2],"pink"=cbp2[8],
                                   "green"=cbp2[4], "yellow"=cbp2[1], "blue"=cbp2[6]),.5),
                    labels=c("Observed", "No Comp", "Low P",
                             "Low Shade", "No Comp, Low P")) +
  theme(text = element_text(size = 14), legend.position = c(.8,.85)) 

100*(length(which(obs_waac_ldgr$obs_ldgr < 0)) / length(obs_waac_ldgr$obs_ldgr))
