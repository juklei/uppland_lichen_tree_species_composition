## lichen lpsac model
##
## First edit: 20190605
## Last edit: 201906118
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(p in 1:nplot){
    for(i in 1:ntree[p]){
      lambda_obs[i,p] <- (plot_richness[p]*i)/(sat_speed[p] + i)
      for(k in 1:nrep){
        obs[i,k,p] ~ dpois(lambda_obs[i,p])
  }}}
  
  ## Process model:
  for(p in 1:nplot){
    ## Richness:
    plot_richness[p] ~ dpois(lambda_rich[p])
    log(lambda_rich[p]) <- alpha_rich + 
                           beta_dec_rich*dec[p,1] +
                           #beta_spruce_rich*spruce[p,1] + 
                           #beta_pine_rich*pine[p,1] + 
                           #beta_nr_tsp_rich*nr_tsp[p,1] + 
                           # beta_dec1*dec1[p] + 
                           # beta_dec2*dec2[p] +
                           beta_dbh_rich*dbh[p,1]
    ## Saturation:
    sat_speed[p] ~ dpois(lambda_sat[p])
    log(lambda_sat[p]) <- alpha_sat + 
                          beta_dec_sat*dec[p,1] +
                          #beta_spruce_sat*spruce[p,1] + 
                          #beta_pine_sat*pine[p,1] +  
                          beta_dbh_sat*dbh[p,1]
  }
  
  ## Priors:
  
  alpha_rich ~ dnorm(0, 0.001) #dgamma(0.001, 0.001)
  alpha_sat ~ dnorm(0, 0.001) #dgamma(0.001, 0.001)
  beta_dec_rich ~ dnorm(0, 0.001)
  beta_dec_sat ~ dnorm(0, 0.001)
  # beta_spruce_rich ~ dnorm(0, 0.001)
  # beta_spruce_sat ~ dnorm(0, 0.001)
  # beta_pine_rich ~ dnorm(0, 0.001)
  # beta_pine_sat ~ dnorm(0, 0.001)
  # beta_nr_tsp_rich ~ dnorm(0, 0.001)
  # beta_dec1 ~ dnorm(0, 0.001)
  # beta_dec2 ~ dnorm(0, 0.001)
  beta_dbh_rich ~ dnorm(0, 0.001)
  beta_dbh_sat ~ dnorm(0, 0.001)
  
  ## Predictions:
  
  # ## For plotting the species accumulation curve:
  # for(p in 1:nplot){
  #     for(m in 1:51){
  #       obs_pred[m,p] ~ dpois((plot_richness[p]*(m-1))/(sat_speed[p] + (m-1)))
  # }}

  for(q in 1:length(dec_pred)){
    log(r_dec[q]) <- alpha_rich + beta_dec_rich*dec_pred[q]
  }
  
  # for(m in 1:length(spruce_pred)){
  #   log(r_spruce[m]) <- alpha_rich + beta_spruce_rich*spruce_pred[m]
  # }
  
  # for(m in 1:length(pine_pred)){
  #   log(r_pine[m]) <- alpha_rich + beta_pine_rich*pine_pred[m]
  # }

}

