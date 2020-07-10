## lichen lpsac model without covariates
##
## First edit: 20190605
## Last edit: 20200710
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
    plot_richness[p] ~ dnorm(mu_rich, 1/sigma_rich^2)
    # pr_sim[p] ~ dnorm(mu_rich, 1/sigma_rich^2)
    ## Saturation:
    sat_speed[p] ~ dgamma(shape, rate)
    shape <- max(0.0001, mu_sat^2/sigma_sat^2)
    rate <- max(0.0001, mu_sat/sigma_sat^2)
  }
  
  ## Priors:
  
  mu_rich ~ dgamma(0.001, 0.001)
  sigma_rich ~ dgamma(0.001, 0.001)
  mu_sat ~ dunif(1, 10)
  sigma_sat ~ dgamma(0.001, 0.001)

  ## Predictions:
  
  ## For plotting the species accumulation curve:
  for(p in 1:nplot){
      for(m in 1:ntree[p]){
        obs_pred[m,p] ~ dpois((plot_richness[p]*m)/(sat_speed[p] + m))
  }}

}

