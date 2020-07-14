## lichen lpsac model % dec as the covariates
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
      for(k in 1:nrep){
        obs[i,k,p] ~ dpois(lambda_obs[i,k,p])
        lambda_obs[i,k,p] <- (plot_richness[k,p]*i)/(sat_speed[k,p] + i)
  }}}
  
  ## Process model:
  for(p in 1:nplot){
    for(k in 1:nrep){
      ## Richness:
      plot_richness[k,p] ~ dpois(lambda_rich[k,p])
      log(lambda_rich[k,p]) <- r_alpha + r_ek[k] +
                               r_beta_dec*dec[p,1] + 
                               r_beta2_dec*dec[p,1]^2 +
                               r_beta_dbh*dbh[p,1]
      ## Saturation:
      sat_speed[k,p] ~ dgamma(shape[k,p], rate[k,p])
      shape[k,p] <- max(0.0001, mu_sat[k,p]^2/sigma_sat^2)
      rate[k,p] <- max(0.0001, mu_sat[k,p]/sigma_sat^2)
      log(mu_sat[k,p]) <- s_alpha + s_ek[k] +
                          s_beta_dec*dec[p,1] + 
                          s_beta2_dec*dec[p,1]^2 +
                          s_beta_dbh*dbh[p,1]
  }}
  
  ## SAC effect:
  for(k in 1:nrep){
    r_ek[k] ~ dnorm(0, sd_r_ek)
    s_ek[k] ~ dnorm(0, sd_s_ek)
  }
  
  ## Priors:
  
  ## Richness:
  r_alpha ~ dnorm(0, 0.001)
  r_beta_dec ~ dnorm(0, 0.001)
  r_beta2_dec ~ dnorm(0, 0.001)
  r_beta_dbh ~ dnorm(0, 0.001)
  sd_r_ek ~ dgamma(0.001, 0.001)
  
  ## Saturation:
  sigma_sat ~ dgamma(0.001, 0.001)
  s_alpha ~ dnorm(0, 0.001)
  s_beta_dec ~ dnorm(0, 0.001)
  s_beta2_dec ~ dnorm(0, 0.001)
  s_beta_dbh ~ dnorm(0, 0.001)
  sd_s_ek ~ dgamma(0.001, 0.001)
  
  ## Plotting predictions and maxima for tree species percentages:
  
  for(m in 1:length(dec_pred)){
    log(r_dec[m]) <- r_alpha + 
                     r_beta_dec*dec_pred[m] + 
                     r_beta2_dec*dec_pred[m]^2
    log(s_dec[m]) <- s_alpha +
                     s_beta_dec*dec_pred[m] +
                     s_beta2_dec*dec_pred[m]^2
  }
  r_dec_max <- - r_beta_dec/(2*r_beta2_dec)
  s_dec_max <- - s_beta_dec/(2*s_beta2_dec)

}

