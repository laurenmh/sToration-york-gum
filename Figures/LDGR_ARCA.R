## Create plots for LDGR ##
load(here("Sparse_model_fits/WAAC_Shade_FinalFit.rdata"))
Post <- rstan::extract(FinalFit)
Shade <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
              reserve = reserve, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij,
              Inclusion_eij = Inclusion_eij, Intra = Intra)
Neighbors <- Shade$SpMatrix
reserve <- Shade$reserve

# load in s and g data 
load(here("SurvivalAndGermination/Germination.rdata"))
germ <- rstan::extract(PrelimFit)
germ<-as.data.frame(germ)
remove(PrelimFit)

load(here("SurvivalAndGermination/Survival.rdata"))
surv <- rstan::extract(PrelimFit)
surv<-as.data.frame(surv)
remove(PrelimFit)

germ <- germ$p.1
germ <- median(germ)
#germination <- data.frame(arca=germ$p.1, waac=germ$p.2)
surv <- surv$p.1
surv <- median(surv)
#survival <- data.frame(arca=surv$p.1, waac=surv$p.2)

# packages
library(tidyverse)
library(rstan)
library(here)
library(patchwork)
library(ggridges)

# user-defined functions
source(here("Sparse_model_fits/functions.R"))

# species and model definitions
foc_sp <- c("ARCA", "HYGL", "TRCY", "WAAC")
mod <- "PSI-lambda_PSI-alpha_re_plot_block" # see README.md for more info

FocalLetter <- "A" # "W", "A", "T", "H"
FocalPrefix <- "ARCA" # "WAAC", "ARCA", "HYGL", "TRCY"
FocalSpecies <- "Arctotheca.calendula" # "Waitzia.acuminata", "Arctotheca.calendula", "Trachymene.cyanopetala", "Hypochaeris.glabra"

#### ggridges plot for average competitive effect ####

alpha_post <- data.frame(NULL)

# load model fit and data
filename <- paste(FocalPrefix, mod, "FinalFit.rdata", sep = "_")

# update with where stan code output is saved
model_fit_path <- paste("~/Dropbox/Shared/sToration/YGW_subgroup/Stan_model_fits")
fp <- paste(model_fit_path, filename, sep = "/")
load(fp)

# get lambdas
dem_param <- "lambda"

# model form for lambda
mform <- formula("~ Reserve.x * Colwell.P_std * Canopy_std")

# Bayes estimators for log-linear model coefs on lambda
beta_lambda_be <- colMeans(
  rstan::extract(mfit_final, pars = "beta_lambda")[[1]]
)

# get lambda estimates for each reserve across canopy x phos
lambdas <- get_lambda(mform, param_estims = beta_lambda_be, sp_name = FocalPrefix, option = "cividis")

# get lambda estimates for each reserve across canopy x phos
beta_eta_be <- colMeans(
  rstan::extract(mfit_final, pars = "beta_eta")[[1]]
)

etas <- get_lambda(mform, param_estims = beta_eta_be, sp_name = FocalPrefix, option = "cividis")

mform_alpha <- formula("~ Colwell.P_std * Canopy_std")

beta_alpha_be <- colMeans(
  rstan::extract(mfit_final, pars = "beta_alpha")[[1]]
)

alphas <- get_lambda(mform_alpha, param_estims = beta_alpha_be, sp_name = FocalPrefix, option = "cividis")

plot_heatmap(mform_alpha, param_estims = beta_alpha_be, sp_name = FocalPrefix, option = "cividis", dem_param = "alpha")
plot_heatmap(mform, param_estims = beta_eta_be, sp_name = FocalPrefix, option = "cividis", dem_param = "eta")


Fecundity <- lambdas$pred / (1+(etas$pred*1)+(alphas$pred*49))
Nt_p_one <- germ*Fecundity + surv*(1-germ)*1
LDGR <- log(Nt_p_one/1)

ncells = 100
df_LDGR <- data.frame(
  Reserve.x = as.factor(rep(c("Bendering", "Perenjori"), each = ncells^2)),
  Colwell.P_std = rep(rep(seq(-1, 1, length.out = ncells), ncells), 2),
  Canopy_std = rep(rep(seq(-1, 1, length.out = ncells), each = ncells), 2),
  LDGR = LDGR
)

df_diff <- data.frame(
  Reserve.x = as.factor(rep(c("Bendering", "Perenjori"), each = ncells^2)),
  Colwell.P_std = rep(rep(seq(-1, 1, length.out = ncells), ncells), 2),
  Canopy_std = rep(rep(seq(-1, 1, length.out = ncells), each = ncells), 2),
  Diff = log(lambdas$pred)/LDGR
)

plotting_LDGR(sp_name = FocalPrefix, df_new = df_LDGR)

plotting_Diff(sp_name = FocalPrefix, df_new = df_diff)



# create plot
alpha_rp <- ggplot(data = alpha_post, aes(x = alpha, y = species, fill = prov)) +
  geom_density_ridges(alpha = 0.7, size = 0.1)+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
  theme(
    panel.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    axis.line = element_line(colour = "darkgrey"),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
    plot.title = element_text(size = 16)
  )+
  scale_fill_manual(values = c("pink", "orange"))+
  xlab("")+
  ylab("")+
  ggtitle(expression(bar(alpha)))



#### Non-generic effects for TRCY ####
rm(list = ls())

# user-defined functions
source(here("Sparse_model_fits/functions.R"))

# TRCY model fit and data
load(
  here("Sparse_model_fits/TRCY_PSI-lambda_PSI-alpha_re_plot_block_FinalFit.rdata")
)

# Get the species that had non-generic effects
ng_sp <- which(colSums(data_final$Q) > 0)
ng_sp_names <- colnames(data_final$A)[ng_sp]


# create new model matrix for prediction
n_new <- 100
df_new <- data.frame(
  Reserve.x = as.factor(rep(c("Bendering", "Perenjori"), each = ncells^2)),
  Colwell.P_std = rep(rep(seq(-1, 1, length.out = ncells), ncells), 2),
  Canopy_std = rep(rep(seq(-1, 1, length.out = ncells), each = ncells), 2)
)

#df_new_P <- data.frame(
#  Reserve.x = factor(
#    rep("perenjori", n_new),
#    levels = c("bendering","perenjori")
#  ),
#  Colwell.P_std = seq(-1, 1, length.out = n_new),
#  Canopy_std = rep(0, n_new)
#)

# construct model matrix
X_new <- model.matrix(
  m_formula,
  data = df_new
)

#X_alpha_new <- model.matrix(
#  ~ Colwell.P_std * Canopy_std,
#  data = df_new_P
#)
Z_alpha_new_P <- model.matrix(
  ~ Reserve.x * Colwell.P_std * Canopy_std,
  data = df_new_P
)

# extract draws from posteriors
posts <- rstan::extract(
  mfit_final,
  pars = c("beta_alpha", "B")
)
n_draws <- dim(posts$beta_alpha)[1]

# create posterior means
estims_mean <- apply(
  posts$beta_alpha, 
  1,
  function(x, mat){
    mat %*% x
  },
  mat = X_alpha_new_P
)

# add to plot df
plot_df_P <- data.frame(
  species = rep("generic", n_new),
  estim = exp(apply(estims_mean, 1, mean)),
  low = exp(apply(estims_mean, 1, quantile, probs = 0.025)),
  high = exp(apply(estims_mean, 1, quantile, probs = 0.975))
)

# create estims for non-generic species
estims_ng_P <- vector(mode = "list", length = 2)
# Outer loop is for each species
for(s in 1:2){
  estims_ng_i <- matrix(nrow = n_new, ncol = n_draws)
  for(i in 1:n_draws){
    estims_ng_i[,i] <- 
      X_alpha_new_P %*% posts$beta_alpha[i,] + Z_alpha_new_P %*% posts$B[i, , ng_sp[s]]
  }
  estims_ng_P[[s]] <- estims_ng_i
}

# create new dfs and add to existing
for(s in 1:length(estims_ng_P)){
  df_i <- data.frame(
    species = ng_sp_names[s],
    estim = exp(apply(estims_ng_P[[s]], 1, mean)),
    low = exp(apply(estims_ng_P[[s]],1, quantile, probs = 0.025)),
    high = exp(apply(estims_ng_P[[s]],1, quantile, probs = 0.975))
  )
  
  plot_df_P <- rbind(plot_df_P, df_i)
}

# add the x value
plot_df_P$phos <- rep(df_new_P$Colwell.P_std, 3)

# factor species in the order we want
plot_df_P$f.species <- factor(
  plot_df_P$species,
  levels = c(
    "generic",
    "Arctotheca.calendula",
    "Lawrencella.rosea"
  )
)







## Non-generic effects for shade
df_new_shade <-  data.frame(
  Reserve.x = factor(
    rep("perenjori", n_new),
    levels = c("bendering","perenjori")
  ),
  Colwell.P_std = rep(0, n_new),
  Canopy_std = seq(-1, 1, length.out = n_new)
)
X_alpha_new_shade <- model.matrix(
  ~ Colwell.P_std * Canopy_std,
  data = df_new_shade
)
Z_alpha_new_shade <- model.matrix(
  ~ Reserve.x * Colwell.P_std * Canopy_std,
  data = df_new_shade
)

# create posterior means
estims_shade <- apply(
  posts$beta_alpha, 
  1,
  function(x, mat){
    mat %*% x
  },
  mat = X_alpha_new_shade
)

# add to plot df
plot_df_shade <- data.frame(
  species = rep("generic", n_new),
  estim = exp(apply(estims_mean, 1, mean)),
  low = exp(apply(estims_mean, 1, quantile, probs = 0.025)),
  high = exp(apply(estims_mean, 1, quantile, probs = 0.975))
)


# create estims for non-generic species
estims_ng_shade <- matrix(nrow = n_new, ncol = n_draws)
for(i in 1:n_draws){
  estims_ng_shade[,i] <- 
    X_alpha_new_shade %*% posts$beta_alpha[i,] + Z_alpha_new_shade %*% posts$B[i, , ng_sp[3]]
}

# create new dfs and add to existing
df_ng_shade <- data.frame(
  species = ng_sp_names[3],
  estim = exp(apply(estims_ng_shade, 1, mean)),
  low = exp(apply(estims_ng_shade, 1, quantile, probs = 0.025)),
  high = exp(apply(estims_ng_shade, 1, quantile, probs = 0.975))
)
plot_df_shade <- rbind(plot_df_shade, df_ng_shade)

# add the x value
plot_df_shade$shade <- rep(df_new_shade$Canopy_std, 2)

# factor species in the order we want
plot_df_shade$f.species <- factor(
  plot_df_shade$species,
  levels = c(
    "generic",
    "Waitzia.acuminata"
  )
)





### Plotting ###

# define a theme
alpha_plot_theme <- theme(
  panel.background = element_blank(),
  plot.title = element_text(size = 12, face = "bold"),
  axis.line = element_line(color = "grey"),
  axis.title.y = element_text(size = 16),
  axis.title.x = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.key = element_blank(),
  legend.title = element_blank(),
  legend.text = ggtext::element_markdown(size = 12),
)

# create plot for phosphorous
plot_P <- ggplot(data = plot_df_P, aes(x = phos, y = estim, color = f.species))+
  geom_line(size = 1)+
  geom_line(aes(y = low), linetype = "dashed")+
  geom_line(aes(y = high), linetype = "dashed")+
  alpha_plot_theme+
  scale_color_manual(
    values = c("darkgrey", "green4", "orange"),
    labels = c("generic", "*A. calendula*", "*L. rosea*")
  )+
  ylab(expression(alpha[ij]))+
  xlab("Standardized Phosphorous")+
  ylim(c(0,0.15))+
  ggtitle("a)")


# create plot for shade
plot_shade <- ggplot(data = plot_df_shade, aes(x = shade, y = estim, color = f.species))+
  geom_line(size = 1)+
  geom_line(aes(y = low), linetype = "dashed")+
  geom_line(aes(y = high), linetype = "dashed")+
  alpha_plot_theme+
  scale_color_manual(
    values = c("darkgrey", "green4"),
    labels = c("generic", "*W. acuminata*")
  )+
  ylab(expression(alpha[ij]))+
  xlab("Standardized Canopy Cover")+
  ylim(c(0,0.15)) +
  ggtitle("b)")

# create total plot title
# plot_title <- ggplot()+
#   geom_text(
#     aes(x = 1, y = 1),
#     label = as.character(expression(paste(italic("T. cyanopetala"), " at Perenjori", sep=""))),
#     parse = T,
#     size = 6
#   )+
#   theme_void(
#   )

ggsave(
  filename = here("Figures/TRCY_alpha_ij_non-generic.svg"),
  plot = plot_P + plot_shade,
  width = 10,
  height = 4,
  units = "in"
)







#### generic effects reg. coefficients ####

# clear the workspace and start over
rm(list = ls())
source(here("Sparse_model_fits/functions.R"))

# load all model fits at once and store in a list
FocalPrefix <- "ARCA" # "WAAC", "ARCA", "HYGL", "TRCY"
FocalList <- c("ARCA", "WAAC", "HYGL", "TRCY")

#ARCA
fits <- vector("list")
model_fit_path <- paste("~/Dropbox/Shared/sToration/YGW_subgroup/Stan_model_fits")
filename <- paste(FocalList[1], mod, "FinalFit.rdata", sep = "_")
fp <- paste(model_fit_path, filename, sep = "/")
load(fp)
fits$ARCA <- mfit_final  

filename <- paste(FocalList[3], mod, "FinalFit.rdata", sep = "_")
fp <- paste(model_fit_path, filename, sep = "/")
#load(here("Sparse_model_fits/HYGL_PSI-lambda_PSI-alpha_re_plot_block_FinalFit.rdata"))
fits$HYGL <- mfit_final

filename <- paste(FocalList[4], mod, "FinalFit.rdata", sep = "_")
fp <- paste(model_fit_path, filename, sep = "/")
#load(here("Sparse_model_fits/TRCY_PSI-lambda_PSI-alpha_re_plot_block_FinalFit.rdata"))
fits$TRCY <- mfit_final

filename <- paste(FocalList[2], mod, "FinalFit.rdata", sep = "_")
fp <- paste(model_fit_path, filename, sep = "/")
#load(here("Sparse_model_fits/WAAC_PSI-lambda_PSI-alpha_re_plot_block_FinalFit.rdata"))
fits$WAAC <- mfit_final

rm(mfit_final, ppc_list, data_final, data_prelim)

# Get posterior means and credible intervals for each
beta_alpha_post <- map(
  fits,
  ~ rstan::extract(.x, pars = "beta_alpha")$beta_alpha
)

# define some useful objects
sp <- names(fits)
prov <- rep(c("exotic", "native"), each = 2)
params <- c("Intercept", "Phosphorous", "Canopy", "Phosphorous:Canopy")
os <- c(-0.2, -0.1, 0.1, 0.2)

# build dataframe
df_effs_alpha <- data.frame(
  species = rep(sp, each = ncol(beta_alpha_post[[1]])),
  provenance = rep(prov, each = ncol(beta_alpha_post[[1]])),
  effect = rep(params, ncol(beta_alpha_post[[1]])),
  # add posterior means
  estim = simplify(map(beta_alpha_post, ~ colMeans(.x))),
  # add bounds of cred. interval
  low = simplify(map(
    beta_alpha_post,
    ~ apply(.x, 2, quantile, probs = 0.025)
  )),
  high = simplify(map(
    beta_alpha_post,
    ~ apply(.x, 2, quantile, probs = 0.975)
  )),
  # add some offset for plotting
  y = simplify(map(
    os,
    ~ .x + c(4:1)
  ))
)

# theme for the plot
theme_effsplot <- theme(
  panel.background = element_rect(colour = "darkgrey", fill = "white"),
  panel.grid.major.x = element_line(colour = "grey", linetype = "dashed"),
  legend.title = element_blank(),
  legend.key = element_blank(),
  axis.text.y = element_text(color = "black")
)

# plot
effs_plot_alpha <- ggplot(data = df_effs_alpha, aes(x = estim, y = y))+
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0)+
  geom_point(
    aes(shape = species),
    size = 2
  )+
  geom_point(
    aes(shape = species, color = provenance),
    size = 1.5
  )+
  theme_effsplot+
  scale_color_manual(values = c("pink", "orange"))+
  scale_shape_manual(values = 15:18)+
  scale_y_continuous(breaks = 1:4, labels = params[4:1])+
  ylab("")+
  xlab(expression(hat(beta)[k]))+
  ggtitle(expression(alpha[ij]))

# save plot
ggsave(
  filename = here("Figures/effects_alpha_generic.svg"),
  plot = effs_plot_alpha,
  height = 4,
  width = 5,
  units = "in"
)





#### conspecific effects reg. coefficients ####

# Get posterior means and credible intervals for conspecific effects
beta_eta_post <- map(
  fits,
  ~ rstan::extract(.x, pars = "beta_eta")$beta_eta
)

# construct linear combos for four params for bendering and perenjori
bend <- rbind(
  c(1,0,0,0,0,0,0,0),
  c(0,0,1,0,0,0,0,0),
  c(0,0,0,1,0,0,0,0),
  c(0,0,0,0,0,0,1,0)
)

peren <- rbind(
  c(1,1,0,0,0,0,0,0),
  c(0,0,1,0,1,0,0,0),
  c(0,0,0,1,0,1,0,0),
  c(0,0,0,0,0,0,1,1)
)

# list of params for each site
b_beta_eta <- map(
  beta_eta_post,
  ~ .x %*% t(bend)
)
p_beta_eta <- map(
  beta_eta_post,
  ~ .x %*% t(peren)
)

# construct a dataframe for each

# build dataframe for bendering
df_effs_intra_b <- df_effs_alpha[,1:3]
df_effs_intra_b <- df_effs_intra_b %>% mutate(
  # add posterior means
  estim = simplify(map(b_beta_eta, ~ colMeans(.x))),
  # add bounds of cred. interval
  low = simplify(map(
    b_beta_eta,
    ~ apply(.x, 2, quantile, probs = 0.025)
  )),
  high = simplify(map(
    b_beta_eta,
    ~ apply(.x, 2, quantile, probs = 0.975)
  )),
  # add some offset for plotting
  y = simplify(map(
    os,
    ~ .x + c(4:1)
  ))
)

# build dataframe for perenjori
df_effs_intra_p <- df_effs_alpha[,1:3]
df_effs_intra_p <- df_effs_intra_p %>% mutate(
  # add posterior means
  estim = simplify(map(p_beta_eta, ~ colMeans(.x))),
  # add bounds of cred. interval
  low = simplify(map(
    p_beta_eta,
    ~ apply(.x, 2, quantile, probs = 0.025)
  )),
  high = simplify(map(
    p_beta_eta,
    ~ apply(.x, 2, quantile, probs = 0.975)
  )),
  # add some offset for plotting
  y = simplify(map(
    os,
    ~ .x + c(4:1)
  ))
)

# plots
intra_b <- ggplot(data = df_effs_intra_b, aes(x = estim, y = y))+
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0)+
  geom_point(
    aes(shape = species),
    size = 2
  )+
  geom_point(
    aes(shape = species, color = provenance),
    size = 1.5
  )+
  theme_effsplot+
  theme(legend.position = "none")+
  scale_color_manual(values = c("pink", "orange"))+
  scale_shape_manual(values = 15:18)+
  scale_y_continuous(breaks = 1:4, labels = params[4:1])+
  ylab("")+
  xlab(expression(hat(beta)[k]))+
  ggtitle(
    expression(paste(alpha[ii], ": Bendering"))
  )

intra_p <- ggplot(data = df_effs_intra_p, aes(x = estim, y = y))+
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0)+
  geom_point(
    aes(shape = species),
    size = 2
  )+
  geom_point(
    aes(shape = species, color = provenance),
    size = 1.5
  )+
  theme_effsplot+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values = c("pink", "orange"))+
  scale_shape_manual(values = 15:18)+
  ylab("")+
  xlab(expression(hat(beta)[k]))+
  ggtitle(
    expression(paste(alpha[ii], ": Perenjori"))
  )

# arrange figures into one
intra_comb <- intra_b + intra_p +
  plot_layout(widths = c(2,2))

# save plot  
ggsave(
  filename = here("Figures/effects_plot_intra.svg"),
  height = 4,
  width = 7,
  units = "in",
  device = "svg"
)





