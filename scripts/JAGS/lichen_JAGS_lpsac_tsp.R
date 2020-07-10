## lichen lpsac model with the tsp number as the covariates
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
                  r_beta_2tsp*tsp_2[p] +
                  r_beta_3tsp*tsp_3[p] +
                  r_beta_4tsp*tsp_4[p] +
                  r_beta_dbh*dbh[p,1]
    ## Saturation:
    sat_speed[p] ~ dgamma(shape[p], rate[p])
    shape[p] <- max(0.0001, mu_sat[p]^2/sigma_sat^2)
    rate[p] <- max(0.0001, mu_sat[p]/sigma_sat^2)
    log(mu_sat[p]) <- s_alpha +
                      s_beta_2tsp*tsp_2[p] +
                      s_beta_3tsp*tsp_3[p] +
                      s_beta_4tsp*tsp_4[p] +
                      s_beta_dbh*dbh[p,1]
  }
  
  ## Priors:
  
  ## Richness:
  sigma_rich ~ dgamma(0.001, 0.001)
  r_alpha ~ dgamma(0.001, 0.001)
  r_beta_2tsp ~ dnorm(0, 0.001)
  r_beta_3tsp ~ dnorm(0, 0.001)
  r_beta_4tsp ~ dnorm(0, 0.001)
  r_beta_dbh ~ dnorm(0, 0.001)
  
  ## Saturation:
  sigma_sat ~ dgamma(0.001, 0.001)
  s_alpha ~ dnorm(0, 0.001)
  s_beta_2tsp ~ dnorm(0, 0.001)
  s_beta_3tsp ~ dnorm(0, 0.001)
  s_beta_4tsp ~ dnorm(0, 0.001)
  s_beta_dbh ~ dnorm(0, 0.001)
  
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

