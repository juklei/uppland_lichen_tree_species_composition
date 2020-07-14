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
    for(k in 1:nrep){
      for(i in 1:ntree[p]){
        obs[i,k,p] ~ dpois(lambda_obs[i,k,p])
        lambda_obs[i,k,p] <- (plot_richness[k,p]*i)/(sat_speed[k,p] + i)
      }}}
  
  ## Process model:
  for(p in 1:nplot){
    for(k in 1:nrep){
      plot_richness[k,p] ~ dpois(lambda_rich[k])
      sat_speed[k,p] ~ dgamma(max(0.0001, mu_sat[k]^2/sigma_sat^2), 
                              max(0.0001, mu_sat[k]/sigma_sat^2))
    }}
  
  ## Priors:
  for(k in 1:nrep){
    lambda_rich[k] ~ dgamma(0.001, 0.001)
    mu_sat[k] ~ dunif(1, 10)
  }
  sigma_sat ~ dgamma(0.001, 0.001)

  ## Predictions:
  
  ## For plotting the species accumulation curve:
  for(p in 1:nplot){
    for(m in 1:ntree[p]){
      for(k in 1:nrep){
        obs_pred[m,k,p] <- (plot_richness[k,p]*m)/(sat_speed[k,p] + m)
      }
      mean_op[m,p] <- mean(obs_pred[m,,p])
      min_op[m,p] <- min(obs_pred[m,,p])
      max_op[m,p] <- max(obs_pred[m,,p])
  }}

}

