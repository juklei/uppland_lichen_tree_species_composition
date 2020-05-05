## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190612
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Tree level:
  for(i in 1:nobs){
  
    richness[i] ~ dbin(p, richness_true[i])
    richness_sim[i] ~ dbin(p, richness_true[i])
    
    richness_true[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha + plot_effect[plot[i]] + 
                      beta_pine*pine[i] + 
                      beta_spruce*spruce[i] +
                      beta_aspen*aspen[i] +
                      beta_oak*oak[i] +
                      beta_alder*alder[i] +
                      beta_dbh*dbh[i,1]
    
  }
  
  ## Plot level:
  for(j in 1:nplot){

    plot_effect[j] ~ dnorm(0, tau_plot)

  }
  
  ## Priors:
  
  p ~ dbeta(0.95, 0.05)
  
  alpha ~ dnorm(0, 0.01)
  beta_pine ~ dnorm(0, 0.01)
  beta_spruce ~ dnorm(0, 0.01)
  beta_aspen ~ dnorm(0, 0.01)
  beta_oak ~ dnorm(0, 0.01)
  beta_alder ~ dnorm(0, 0.01)
  beta_dbh ~ dnorm(0, 0.01)
  
  sigma_plot ~ dgamma(0.001, 0.001)
  tau_plot <- 1/sigma_plot^2

  # ## Model validation:
  # 
  # ## Bayesian p-value:
  # mean_richness <- mean(richness[])
  # mean_richness_sim <- mean(richness_sim[])
  # p_mean <- step(mean_richness_sim - mean_richness)
  # 
  # ## Coefficient of variation:
  # cv_richness <- sd(richness[])/mean_richness
  # cv_richness_sim <- sd(richness_sim[])/mean_richness_sim
  # p_cv <- step(cv_richness - cv_richness_sim)
  # 
  # ## Model fit:
  # for(m in 1:nobs){
  # 
  #   sq[m] <- (richness[m] - p*richness_true[m])^2
  #   sq_sim[m] <- (richness_sim[m] - p*richness_true[m])^2
  #   
  # }
  # 
  # fit <- sum(sq[])
  # fit_sim <- sum(sq_sim[])
  # p_fit <- step(fit_sim - fit)
  
  ## Predictions:
  
  ## Tree species:
  log(birch_mean) <- alpha
  log(pine_mean) <- alpha + beta_pine
  log(spruce_mean) <- alpha + beta_spruce
  log(aspen_mean) <- alpha + beta_aspen
  log(oak_mean) <- alpha + beta_oak
  log(alder_mean) <- alpha + beta_alder
  
  ## Calculate differences:
  aspen_birch <- birch_mean - aspen_mean
  aspen_oak <- oak_mean - aspen_mean
  aspen_alder <- alder_mean - aspen_mean
  aspen_pine <- pine_mean - aspen_mean
  aspen_spruce <- spruce_mean - aspen_mean
  birch_oak <- oak_mean - birch_mean
  birch_alder <- alder_mean - birch_mean
  birch_pine <- pine_mean - birch_mean
  birch_spruce <- spruce_mean - birch_mean
  oak_alder <- alder_mean - oak_mean
  oak_pine <- pine_mean - oak_mean
  oak_spruce <- spruce_mean - oak_mean
  alder_pine <- pine_mean - alder_mean
  alder_spruce <- spruce_mean - alder_mean
  pine_spruce <- spruce_mean - pine_mean
  
}


