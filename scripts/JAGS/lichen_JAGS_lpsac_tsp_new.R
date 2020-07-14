## lichen lpsac model with the tsp number as the covariates
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
                               r_beta_2tsp*tsp_2[p] +
                               r_beta_3tsp*tsp_3[p] +
                               r_beta_4tsp*tsp_4[p] +
                               r_beta_dbh*dbh[p,1]
      ## Saturation:
      sat_speed[k,p] ~ dgamma(shape[k,p], rate[k,p])
      shape[k,p] <- max(0.0001, mu_sat[k,p]^2/sigma_sat^2)
      rate[k,p] <- max(0.0001, mu_sat[k,p]/sigma_sat^2)
      log(mu_sat[k,p]) <- s_alpha + s_ek[k] +
                          s_beta_2tsp*tsp_2[p] +
                          s_beta_3tsp*tsp_3[p] +
                          s_beta_4tsp*tsp_4[p] +
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
  r_beta_2tsp ~ dnorm(0, 0.001)
  r_beta_3tsp ~ dnorm(0, 0.001)
  r_beta_4tsp ~ dnorm(0, 0.001)
  r_beta_dbh ~ dnorm(0, 0.001)
  sd_r_ek ~ dgamma(0.001, 0.001)
  
  ## Saturation:
  sigma_sat ~ dgamma(0.001, 0.001)
  s_alpha ~ dnorm(0, 0.001)
  s_beta_2tsp ~ dnorm(0, 0.001)
  s_beta_3tsp ~ dnorm(0, 0.001)
  s_beta_4tsp ~ dnorm(0, 0.001)
  s_beta_dbh ~ dnorm(0, 0.001)
  sd_s_ek ~ dgamma(0.001, 0.001)
  
  ## Predictions:
  
  ## Richness:
  
  r_1tsp <- r_alpha
  r_2tsp <- r_alpha + r_beta_2tsp
  r_3tsp <- r_alpha + r_beta_3tsp
  r_4tsp <- r_alpha + r_beta_4tsp

  r_diff_21 <- r_2tsp - r_1tsp
  r_diff_31 <- r_3tsp - r_1tsp
  r_diff_41 <- r_4tsp - r_1tsp
  r_diff_32 <- r_3tsp - r_2tsp
  r_diff_42 <- r_4tsp - r_2tsp
  r_diff_43 <- r_4tsp - r_3tsp
  
  ## Saturation:
  
  log(s_1tsp) <- s_alpha
  log(s_2tsp) <- s_alpha + s_beta_2tsp
  log(s_3tsp) <- s_alpha + s_beta_3tsp
  log(s_4tsp) <- s_alpha + s_beta_4tsp
  
  s_diff_21 <- s_2tsp - s_1tsp
  s_diff_31 <- s_3tsp - s_1tsp
  s_diff_41 <- s_4tsp - s_1tsp
  s_diff_32 <- s_3tsp - s_2tsp
  s_diff_42 <- s_4tsp - s_2tsp
  s_diff_43 <- s_4tsp - s_3tsp
  
  
}

