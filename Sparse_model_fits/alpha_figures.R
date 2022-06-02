## Create plots for alpha values ##

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
  

#### ggridges plot for average competitive effect ####
  
  alpha_post <- data.frame(NULL)
  
# load model fit and data
  filename <- paste(foc_sp[1], mod, "FinalFit.rdata", sep = "_")
  fp <- paste(here("Sparse_model_fits"), filename, sep = "/")
  load(fp)
  
# extract posteriors
  beta_alpha <- as.data.frame(rstan::extract(mfit_final, pars = "beta_alpha"))
  
# add to dataframe
  alpha_post <- data.frame(
    alpha = exp(beta_alpha[,1]),
    species = foc_sp[1]
  )
  
# Remove model fits to avoid confusion
  rm(beta_alpha, data_final, data_prelim, mfit_final, ppc_list)
  
# add the remaining species
  for(s in 2:length(foc_sp)){
    
    # load model fit and data
    filename <- paste(foc_sp[s], mod, "FinalFit.rdata", sep = "_")
    fp <- paste(here("Sparse_model_fits"), filename, sep = "/")
    load(fp)
    
    # extract posteriors
    beta_alpha <- as.data.frame(rstan::extract(mfit_final, pars = "beta_alpha"))
    
    # create dataframe
    alpha_post_s <- data.frame(
      alpha = exp(beta_alpha[,1]),
      species = foc_sp[s]
    )
    
    # add to existing
    alpha_post <- rbind(alpha_post, alpha_post_s)
    
    # Remove model fits to avoid confusion
    rm(beta_alpha, data_final, data_prelim, mfit_final, ppc_list)
    
  }

# Add info about provenance
  prov <- data.frame(
    species = foc_sp,
    prov = rep(
      c("exotic", "native"),
      each = 2
    )
  )
  alpha_post <- left_join(
    alpha_post,
    prov,
    by = "species"
  )
  
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

 ggsave(
   filename = here("Figures/alpha_bar_ridgeplot.svg"),
   alpha_rp,
   width = 10,
   height = 8,
   units = "cm",
   device = "svg"
 )  

 

 
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
 df_new_P <- data.frame(
   Reserve.x = factor(
     rep("perenjori", n_new),
     levels = c("bendering","perenjori")
   ),
   Colwell.P_std = seq(-1, 1, length.out = n_new),
   Canopy_std = rep(0, n_new)
 )
 X_alpha_new_P <- model.matrix(
   ~ Colwell.P_std * Canopy_std,
   data = df_new_P
 )
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
  
  