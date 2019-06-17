## lichen lto model
##
## First edit: 20190605
## Last edit: 20190614
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
        obs_sim[i,k,p] ~ dpois(lambda_obs[i,p])
        
  }}}
  
  ## Process model:
  for(p in 1:nplot){
    
    plot_richness[p] ~ dpois(lambda_rich[p])
    log(lambda_rich[p]) <- alpha_rich + 
                           beta_dec_rich*dec[p,1] + 
                           beta_dbh_rich*dbh[p,1]
    
    sat_speed[p] ~ dpois(lambda_sat[p])
    log(lambda_sat[p]) <- alpha_sat + 
                          beta_dec_sat*dec[p,1] +
                          beta_dbh_sat*dbh[p,1]

  }
  
  ## Priors:
  
  alpha_rich ~ dnorm(0, 0.001)
  alpha_sat ~ dnorm(0, 0.001)
  beta_dec_rich ~ dnorm(0, 0.001)
  beta_dec_sat ~ dnorm(0, 0.001)
  beta_dbh_rich ~ dnorm(0, 0.001)
  beta_dbh_sat ~ dnorm(0, 0.001)
  # tau_plot[k] <- 1/sigma_plot^2
  # sigma_plot[k] ~ dgamma(0.001, 0.001)

  ## Predictions:
  
  ## Percent deciduous:
  for(m in 1:length(dec_pred)){
    
    log(r_dec[m]) <- alpha_rich + beta_dec_rich*dec_pred[m]

  }
  

  
  
}

