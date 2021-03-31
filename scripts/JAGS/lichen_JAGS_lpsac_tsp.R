## lichen lpsac model tsp as the covariates
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(p in 1:nplot){
    for(i in 1:ntree[p]){
      obs[i,p] ~ dnorm(mu_obs[i,p], 1/sd_obs[p]^2)T(0,)
      mu_obs[i,p] <- (gdiv[p]*i)/(bdiv[p]*u - u + i)
  }}
  
  ## Process model:
  for(p in 1:nplot){
    ## Gamma diversity:
    gdiv[p] ~ dnorm(mu_gdiv[p], 1/sigma_gdiv^2)
    mu_gdiv[p] <- g_icpt + 
                  g_2tsp*tsp_2[p] + 
                  g_3tsp*tsp_3[p] +
                  g_4tsp*tsp_4[p] +
                  g_dbh*dbh[p]
    ## Beta diversity:
    bdiv[p] ~ dgamma(mu_bdiv[p]^2/sigma_bdiv^2, mu_bdiv[p]/sigma_bdiv^2)
    log(mu_bdiv[p]) <- b_icpt +                  
                       b_2tsp*tsp_2[p] + 
                       b_3tsp*tsp_3[p] +
                       b_4tsp*tsp_4[p] + 
                       b_dbh*dbh[p]
    }
  
  ## Priors:
  
  ## Observation model:
  for(p in 1:nplot){sd_obs[p] ~ dgamma(0.001, 0.001)}
  
  ## Process model:
  ## Gamma diversity:
  sigma_gdiv ~ dt(0, pow(2.5,-2), 1)T(0,)
  g_icpt ~ dgamma(0.001, 0.001)
  g_2tsp ~ dnorm(0, 0.01)
  g_3tsp ~ dnorm(0, 0.01)
  g_4tsp ~ dnorm(0, 0.01)
  g_dbh ~ dnorm(0, 0.01)
  ## Beta diversity:
  sigma_bdiv ~ dt(0, pow(2.5,-2), 1)T(0,)
  b_icpt ~ dgamma(0.001, 0.001)
  b_2tsp ~ dnorm(0, 0.01)
  b_3tsp ~ dnorm(0, 0.01)
  b_4tsp ~ dnorm(0, 0.01)
  b_dbh ~ dnorm(0, 0.01)
  
  ## Predictions:
  
  ## Gamma diversity:
  gdiv_1tsp <- g_icpt
  gdiv_2tsp <- g_icpt + g_2tsp
  gdiv_3tsp <- g_icpt + g_3tsp
  gdiv_4tsp <- g_icpt + g_4tsp
  gdiv_diff_21 <- gdiv_2tsp - gdiv_1tsp
  gdiv_diff_31 <- gdiv_3tsp - gdiv_1tsp
  gdiv_diff_41 <- gdiv_4tsp - gdiv_1tsp
  gdiv_diff_32 <- gdiv_3tsp - gdiv_2tsp
  gdiv_diff_42 <- gdiv_4tsp - gdiv_2tsp
  gdiv_diff_43 <- gdiv_4tsp - gdiv_3tsp
  ## Beta diversity:
  log(bdiv_1tsp) <- b_icpt
  log(bdiv_2tsp) <- b_icpt + b_2tsp
  log(bdiv_3tsp) <- b_icpt + b_3tsp
  log(bdiv_4tsp) <- b_icpt + b_4tsp
  bdiv_diff_21 <- bdiv_2tsp - bdiv_1tsp
  bdiv_diff_31 <- bdiv_3tsp - bdiv_1tsp
  bdiv_diff_41 <- bdiv_4tsp - bdiv_1tsp
  bdiv_diff_32 <- bdiv_3tsp - bdiv_2tsp
  bdiv_diff_42 <- bdiv_4tsp - bdiv_2tsp
  bdiv_diff_43 <- bdiv_4tsp - bdiv_3tsp

  ## For plotting the species accumulation curve:
  for(p in 1:nplot){
    for(m in 1:ntree[p]){
      obs_pred[m,p] <- (gdiv[p]*m)/(bdiv[p]*u - u + m)
    }}
  
}

