## lichen lpsac model % spruce as the covariates
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
    plot_richness[p] ~ dnorm(mu_rich[p], 1/sigma_rich^2)
    mu_rich[p] <- r_alpha +
                  r_beta_spruce*spruce[p,1] + r_beta2_spruce*spruce[p,1]^2 +
                  r_beta_dbh*dbh[p,1]
    ## Saturation:
    sat_speed[p] ~ dgamma(shape[p], rate[p])
    shape[p] <- max(0.0001, mu_sat[p]^2/sigma_sat^2)
    rate[p] <- max(0.0001, mu_sat[p]/sigma_sat^2)
    log(mu_sat[p]) <- s_alpha +
                      s_beta_spruce*spruce[p,1] + s_beta2_spruce*spruce[p,1]^2 +
                      s_beta_dbh*dbh[p,1]
  }
  
  ## Priors:
  
  ## Richness:
  sigma_rich ~ dgamma(0.001, 0.001)
  r_alpha ~ dgamma(0.001, 0.001)
  r_beta_spruce ~ dnorm(0, 0.001)
  r_beta2_spruce ~ dnorm(0, 0.001)
  r_beta_dbh ~ dnorm(0, 0.001)
  
  ## Saturation:
  sigma_sat ~ dgamma(0.001, 0.001)
  s_alpha ~ dnorm(0, 0.001)
  s_beta_spruce ~ dnorm(0, 0.001)
  s_beta2_spruce ~ dnorm(0, 0.001)
  s_beta_dbh ~ dnorm(0, 0.001)
  
  ## Plotting predictions and maxima for tree species percentages:
  
  for(m in 1:length(spruce_pred)){
    r_spruce[m] <- r_alpha +
                   r_beta_spruce*spruce_pred[m] +
                   r_beta2_spruce*spruce_pred[m]^2
    log(s_spruce[m]) <- s_alpha +
                        s_beta_spruce*spruce_pred[m] +
                        s_beta2_spruce*spruce_pred[m]^2
  }
  r_spruce_max <- - r_beta_spruce/(2*r_beta2_spruce)
  s_spruce_max <- - s_beta_spruce/(2*s_beta2_spruce)

}

