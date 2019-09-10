## lichen lpsac model
##
## First edit: 20190605
## Last edit: 201906118
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
    # plot_richness[p] ~ dpois(lambda_rich[p])
    # log(lambda_rich[p]) <- alpha_rich + 
    plot_richness[p] ~ dnorm(mu_rich[p], 1/sigma_rich^2)
    # pr_sim[p] ~ dnorm(mu_rich[p], 1/sigma_rich^2)
    mu_rich[p] <- alpha_rich +
                  # beta_dec_rich*dec[p,1] +
                  # beta2_dec_rich*dec[p,1]^2 +
                  beta_spruce_rich*spruce[p,1] +
                  beta2_spruce_rich*spruce[p,1]^2 +
                  # beta_pine_rich*pine[p,1] +
                  # beta2_pine_rich*pine[p,1]^2 +
                  # beta_nr_tsp_rich*nr_tsp[p,1] +
                  beta_dbh_rich*dbh[p,1]
    ## Saturation:
    # sat_speed[p] ~ dpois(lambda_sat[p])
    # log(lambda_sat[p]) <- alpha_sat +
    sat_speed[p] ~ dnorm(mu_sat, 1/sigma_sat^2)
    # mu_sat[p] <- alpha_sat +
    #              beta_dec_sat*dec[p,1] +
    #              beta_spruce_sat*spruce[p,1] +
    #              beta_pine_sat*pine[p,1] +
    #              beta_dbh_sat*dbh[p,1]
  }
  
  ## Priors:
  
  # for(p in 1:nplot){sat_speed[p] ~ dunif(1, 10)}
  
  sigma_rich ~ dgamma(0.001, 0.001)
  sigma_sat ~ dgamma(0.001, 0.001)
  alpha_rich ~ dgamma(0.001, 0.001)
  mu_sat ~ dunif(1, 10)
  # beta_dec_rich ~ dnorm(0, 0.001)
  # beta2_dec_rich ~ dnorm(0, 0.001)
  # beta_dec_sat ~ dnorm(0, 0.001)
  beta_spruce_rich ~ dnorm(0, 0.001)
  beta2_spruce_rich ~ dnorm(0, 0.001)
  # beta_spruce_sat ~ dnorm(0, 0.001)
  # beta_pine_rich ~ dnorm(0, 0.001)
  # beta2_pine_rich ~ dnorm(0, 0.001)
  # beta_pine_sat ~ dnorm(0, 0.001)
  # beta_nr_tsp_rich ~ dnorm(0, 0.001)
  beta_dbh_rich ~ dnorm(0, 0.001)
  # beta_dbh_sat ~ dnorm(0, 0.001)
  
  ## Predictions:
  
  # ## For plotting the species accumulation curve:
  # for(p in 1:nplot){
  #     for(m in 1:51){
  #       obs_pred[m,p] ~ dpois((plot_richness[p]*(m-1))/(sat_speed[p] + (m-1)))
  # }}

  # for(q in 1:length(dec_pred)){
  #   r_dec[q] <- alpha_rich +
  #               beta_dec_rich*dec_pred[q] +
  #               beta2_dec_rich*dec_pred[q]^2
  # }
  
  for(m in 1:length(spruce_pred)){
    r_spruce[m] <- alpha_rich +
                   beta_spruce_rich*spruce_pred[m] +
                   beta2_spruce_rich*spruce_pred[m]^2
  }
  
  # for(m in 1:length(pine_pred)){
  #   r_pine[m] <- alpha_rich +
  #                beta_pine_rich*pine_pred[m] +
  #                beta2_pine_rich*pine_pred[m]^2
  # }

}

