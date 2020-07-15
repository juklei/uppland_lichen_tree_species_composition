## lichen lpsac model % pine as the covariates
##
## First edit: 20190605
## Last edit: 20200714
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(p in 1:nplot){
    for(i in 1:ntree[p]){
      obs[i,p] ~ dnorm(mu_obs[i,p], 1/sd_obs[p]^2)T(0,)
      mu_obs[i,p] <- (plot_richness[p]*i)/(sat_speed[p] + i)
  }}
  
  ## Process model:
  for(p in 1:nplot){
    ## Richness:
    plot_richness[p] ~ dnorm(mu_rich[p], 1/sigma_rich^2)
    mu_rich[p] <- r_alpha +
                  r_beta_pine*pine[p,1] + r_beta2_pine*pine[p,1]^2 +
                  r_beta_dbh*dbh[p,1]
    ## Saturation:
    sat_speed[p] ~ dgamma(shape[p], rate[p])
    shape[p] <- max(0.0001, mu_sat[p]^2/sigma_sat^2)
    rate[p] <- max(0.0001, mu_sat[p]/sigma_sat^2)
    log(mu_sat[p]) <- s_alpha +
                      s_beta_pine*pine[p,1] + s_beta2_pine*pine[p,1]^2 +
                      s_beta_dbh*dbh[p,1]
  }
  
  ## Priors:
  
  for(p in 1:nplot){sd_obs[p] ~ dgamma(0.001, 0.001)}
  
  ## Richness:
  sigma_rich ~ dgamma(0.001, 0.001)
  r_alpha ~ dgamma(0.001, 0.001)
  r_beta_pine ~ dnorm(0, 0.001)
  r_beta2_pine ~ dnorm(0, 0.001)
  r_beta_dbh ~ dnorm(0, 0.001)
  
  ## Saturation:
  sigma_sat ~ dgamma(0.001, 0.001)
  s_alpha ~ dnorm(0, 0.001)
  s_beta_pine ~ dnorm(0, 0.001)
  s_beta2_pine ~ dnorm(0, 0.001)
  s_beta_dbh ~ dnorm(0, 0.001)
  
  ## Plotting predictions and maxima for tree species percentages:
  for(m in 1:length(pine_pred)){
    r_pine[m] <- r_alpha + r_beta_pine*pine_pred[m] + r_beta2_pine*pine_pred[m]^2
    log(s_pine[m]) <- s_alpha + s_beta_pine*pine_pred[m] + s_beta2_pine*pine_pred[m]^2
  }
  r_pine_max <- - r_beta_pine/(2*r_beta2_pine)
  s_pine_max <- - s_beta_pine/(2*s_beta2_pine)

}

