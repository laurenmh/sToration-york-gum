// This Stan model is a fecundity model (Ricker version) for Arctotheca calendula
// seed count (already extrapolated)


data{
  // here we declare what parameters come directly from out data 
  int<lower = 1> N; // the number of observations (i.e. focal plants) 
  int Fecundity[N]; // the seed count value for each focal plant
  matrix[N,4] ModelMatrix; // matrix of neighbour abundances, with the no. cols = no. neighbours(/groups)
  // int<lower =1> C; where C = no. neighbour species 
  // matrix[N,C] ifm_alpha;
  int Plot;
  int Block;
  real<lower = 1> Phos[Block]; // soil Phosphorus (Colwell P) value (per Block) - being treated for now as continuous
  real<lower = 1> Shade[Plot]; // shade (canopy cover) percentage (per Plot) - being treated for now as continuous
}

parameters{
  // here we declare the parameters that we want the model to estimate 
  real Beta_1[5]; // Phos effect
  real Beta_2[5]; // Shade effect
  real Beta_3[5]; // interaction effect
  // should be the number of neighbours plus lambda in square brackets
}

model{
  // create a vector of predictions
  real F_hat[N];
  real lambda[N];

// something like the below code for an alternative to writing out all the alpha's but this is for multiple species 
    //vector[N] mu; // the linear predictor
     //for(n in 1:N) {
     //    mu[n] = exp(lambda[species_ID[n]] - dot_product(X[n], ifm_alpha[species_ID[n], ]));  
  //}
  
  real alpha_intra[N];
  real alpha_InvGrass[N];
  real alpha_InvForb[N];
  real alpha_NatForb[N];
  
  // set priors (super wide for now)
  Beta_1 ~ normal(0,1000);
  Beta_2 ~ normal(0,1000);
  Beta_3 ~ normal(0,1000);
  
  for(i in 1:N){
   lambda[i] = Beta_1[1]*Phos[Block[i]] + Beta_2[1]*Shade[Plot[i]] + Beta_3[1]*Phos[Block[i]]*Shade[Plot[i]]; // not sure about this indexing
   alpha_intra[i] = Beta_1[2]*Phos[i] + Beta_2[2]*Shade[i] + Beta_3[2]*Phos[i]*Shade[i]; 
   alpha_InvGrass[i] = Beta_1[3]*Phos[i] + Beta_2[3]*Shade[i] + Beta_3[3]*Phos[i]*Shade[i]; 
   alpha_InvForb[i] = Beta_1[4]*Phos[i] + Beta_2[4]*Shade[i] + Beta_3[4]*Phos[i]*Shade[i]; 
   alpha_NatForb[i] = Beta_1[5]*Phos[i] + Beta_2[5]*Shade[i] + Beta_3[5]*Phos[i]*Shade[i]; 
  }
  
}

  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda[i] * exp(alpha_intra[i] * ModelMatrix[i,1] + alpha_InvGrass[i] * ModelMatrix[i,2] + alpha_InvForb[i] * ModelMatrix[i,3] + alpha_NatForb[i] * ModelMatrix[i,4]); 
  }
  
  Fecundity ~ poisson(F_hat);
}







