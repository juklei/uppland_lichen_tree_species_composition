## lichen lpsac model % dec, spruce, or pine as the covariates
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
    mu_gdiv[p] <- g_icpt + g_perc*perc[p] + g_perc2*perc[p]^2 + g_dbh*dbh[p]
    ## Beta diversity:
    bdiv[p] ~ dgamma(mu_bdiv[p]^2/sigma_bdiv^2, mu_bdiv[p]/sigma_bdiv^2)
    log(mu_bdiv[p]) <- b_icpt + b_perc*perc[p] + b_perc2*perc[p]^2 + b_dbh*dbh[p]
  }

  ## Priors:
  
  ## Observation model:
  for(p in 1:nplot){sd_obs[p] ~ dgamma(0.001, 0.001)}
  
  ## Process model:
  ## Gamma diversity:
  sigma_gdiv ~ dt(0, pow(2.5,-2), 1)T(0,)
  g_icpt ~ dgamma(0.001, 0.001)
  g_perc ~ dnorm(0, 0.01)
  g_perc2 ~ dnorm(0, 0.01)
  g_dbh ~ dnorm(0, 0.01)
  ## Beta diversity:
  sigma_bdiv ~ dt(0, pow(2.5,-2), 1)T(0,)
  b_icpt ~ dgamma(0.001, 0.001)
  b_perc ~ dnorm(0, 0.01)
  b_perc2 ~ dnorm(0, 0.01)
  b_dbh ~ dnorm(0, 0.01) 
  
  ## Predictions:
  
  ## Diversity metrics predictions:
  for(m in 1:length(perc_pred)){
    ## Gamma diversity:
    gdiv_pred[m] <- g_icpt + g_perc*perc_pred[m] + g_perc2*perc_pred[m]^2
    ## Beta diversity:
    log(bdiv_pred[m]) <- b_icpt + b_perc*perc_pred[m] + b_perc2*perc_pred[m]^2
  }
  
  ## Maxima:
  g_perc_max <- - g_perc/(2*g_perc2)
  b_perc_max <- - b_perc/(2*b_perc2)
  
}

