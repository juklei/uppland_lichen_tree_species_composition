## lichen lpsac model
##
## First edit: 20190605
## Last edit: 20191014
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
    # pr_sim[p] ~ dnorm(mu_rich[p], 1/sigma_rich^2)
    # mu_rich[p] <- alpha +
    #               # beta_dec*dec[p,1] +
    #               # beta2_dec*dec[p,1]^2 +
    #               # beta_spruce*spruce[p,1] +
    #               # beta2_spruce*spruce[p,1]^2 +
    #               # beta_pine*pine[p,1] +
    #               # beta2_pine*pine[p,1]^2 +
    #               beta_2tsp*tsp_2[p] +
    #               beta_3tsp*tsp_3[p] +
    #               beta_4tsp*tsp_4[p] +
    #               beta_dbh*dbh[p,1]
    ## Saturation:
    sat_speed[p] ~ dnorm(mu_sat, 1/sigma_sat^2)
  }
  
  ## Priors:
  
  mu_rich ~ dgamma(0.001, 0.001)
  sigma_rich ~ dgamma(0.001, 0.001)
  # alpha ~ dgamma(0.001, 0.001)
  # beta_dec ~ dnorm(0, 0.001)
  # beta2_dec ~ dnorm(0, 0.001)
  # beta_spruce ~ dnorm(0, 0.001)
  # beta2_spruce ~ dnorm(0, 0.001)
  # beta_pine ~ dnorm(0, 0.001)
  # beta2_pine ~ dnorm(0, 0.001)
  # beta_2tsp ~ dnorm(0, 0.001)
  # beta_3tsp ~ dnorm(0, 0.001)
  # beta_4tsp ~ dnorm(0, 0.001)
  # beta_dbh ~ dnorm(0, 0.001)
  mu_sat ~ dunif(1, 10)
  sigma_sat ~ dgamma(0.001, 0.001)

  ## Predictions:
  
  # ## For plotting the species accumulation curve:
  # for(p in 1:nplot){
  #     for(m in 1:ntree[p]){
  #       obs_pred[m,p] ~ dpois((plot_richness[p]*m)/(sat_speed[p] + m))
  # }}
  #
  # # Plotting predictions for tree species percentages:
  # for(m in 1:length(dec_pred)){
  #   r_dec[m] <- alpha +
  #               beta_dec*dec_pred[m] +
  #               beta2_dec*dec_pred[m]^2
  # }
  # ## Maximum:
  # dec_max <- - beta_dec/(2*beta2_dec)
  # ## Non-presence vs. monoculture:
  # diff_100vs0_dec <- beta_dec*(dec_pred[length(dec_pred)] - dec_pred[1]) +
  #                    beta2_dec*(dec_pred[length(dec_pred)]^2 - dec_pred[1]^2)
  #
  # for(m in 1:length(spruce_pred)){
  #   r_spruce[m] <- alpha +
  #                  beta_spruce*spruce_pred[m] +
  #                  beta2_spruce*spruce_pred[m]^2
  # }
  # spruce_max <- - beta_spruce/(2*beta2_spruce)
  # diff_100vs0_spruce <- beta_spruce*(spruce_pred[length(spruce_pred)] -
  #                                    spruce_pred[1]) +
  #                       beta2_spruce*(spruce_pred[length(spruce_pred)]^2 -
  #                                     spruce_pred[1]^2)
  # 
  # for(m in 1:length(pine_pred)){
  #   r_pine[m] <- alpha +
  #                beta_pine*pine_pred[m] +
  #                beta2_pine*pine_pred[m]^2
  # }
  # pine_max <- - beta_pine/(2*beta2_pine)
  # diff_100vs0_pine <- beta_pine*(pine_pred[length(pine_pred)] - pine_pred[1]) +
  #                     beta2_pine*(pine_pred[length(pine_pred)]^2 -
  #                                 pine_pred[1]^2)
  #
  # ## Number of tree species:
  # r_1tsp <- alpha
  # r_2tsp <- alpha + beta_2tsp
  # r_3tsp <- alpha + beta_3tsp
  # r_4tsp <- alpha + beta_4tsp
  # 
  # diff_21 <- r_2tsp - r_1tsp
  # diff_31 <- r_3tsp - r_1tsp
  # diff_41 <- r_4tsp - r_1tsp
  # diff_32 <- r_3tsp - r_2tsp
  # diff_42 <- r_4tsp - r_2tsp
  # diff_43 <- r_4tsp - r_3tsp
  
}

