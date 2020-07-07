## lichen lto model
##
## First edit: 20200706
## Last edit: 20200706
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Tree level:
  for(k in 1:nspecies){
    for(p in 1:nplot){
      for(i in 1:ntree[p]){
        obs[i,k,p] ~ dbern(p_occ[i,k,p])
        logit(p_occ[i,k,p]) <- alpha[k] + plot_effect[k,p] + 
                               beta_pine[k]*pine[i,p] + 
                               beta_spruce[k]*spruce[i,p] +
                               beta_aspen[k]*aspen[i,p] +
                               beta_oak[k]*oak[i,p] +
                               beta_alder[k]*alder[i,p] +
                               beta_dbh[k]*dbh[i,p]
  }}}
  
  ## Plot level:
  for(k in 1:nspecies){
    for(p in 1:nplot){
      plot_effect[k,p] ~ dnorm(0, 1/sigma_plot[k]^2)
  }}
  
  ## Priors:
  
  for(k in 1:nspecies){
    alpha[k]~ dnorm(mu_alpha, 1/sd_alpha^2)
    beta_pine[k] ~ dnorm(mu_b_pine, 1/sd_b_pine^2)
    beta_spruce[k] ~ dnorm(mu_b_spruce, 1/sd_b_spruce^2)
    beta_aspen[k] ~ dnorm(mu_b_aspen, 1/sd_b_aspen^2)
    beta_oak[k] ~ dnorm(mu_b_oak, 1/sd_b_oak^2)
    beta_alder[k] ~ dnorm(mu_b_alder, 1/sd_b_alder^2)
    beta_dbh[k] ~ dnorm(mu_b_dbh, 1/sd_b_dbh^2)
    sigma_plot[k] ~ dt(0, pow(u_sp,-2), 1)T(0,)
  }
  
  ## Hyperpriors:
  
  mu_alpha ~ dnorm(0, 0.1)
  sd_alpha ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pine ~ dnorm(0, 0.1)
  sd_b_pine ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_spruce ~ dnorm(0, 0.1)
  sd_b_spruce ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_aspen ~ dnorm(0, 0.1)
  sd_b_aspen ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_oak ~ dnorm(0, 0.1)
  sd_b_oak ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_alder ~ dnorm(0, 0.1)
  sd_b_alder ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_dbh ~ dnorm(0, 0.1)
  sd_b_dbh ~ dt(0, pow(2.5,-2), 1)T(0,)
  u_sp ~ dunif(0, 5)

  ## Model validation:

  ##...
  
  ## Posterior:
  
  ## Tree species specific and tree specific occurrence:
  for(k in 1:nspecies){
    logit(birch_out[k]) <- alpha[k]
    logit(pine_out[k]) <- alpha[k] + beta_pine[k]
    logit(spruce_out[k]) <- alpha[k] + beta_spruce[k]
    logit(aspen_out[k]) <- alpha[k] + beta_aspen[k]
    logit(oak_out[k]) <- alpha[k] + beta_oak[k]
    logit(alder_out[k]) <- alpha[k] + beta_alder[k]
    tree_out[k] <- sum(birch_out[k], pine_out[k], spruce_out[k], aspen_out[k],
                       oak_out[k], alder_out[k])/6
  }
  
  ## The tree species specific expected richness:
  birch_r <- sum(birch_out[])
  pine_r <- sum(pine_out[])
  spruce_r <- sum(spruce_out[])
  aspen_r <- sum(aspen_out[])
  oak_r <- sum(oak_out[])
  alder_r <- sum(alder_out[])
  tree_r <- sum(tree_out[])
  
  ## Calculate differences:
  birch_aspen <- birch_r - aspen_r
  oak_aspen <- oak_r - aspen_r
  alder_aspen <- alder_r - aspen_r
  pine_aspen <- pine_r - aspen_r
  spruce_aspen <- spruce_r - aspen_r
  oak_birch <- oak_r - birch_r
  alder_birch <- alder_r - birch_r
  pine_birch <- pine_r - birch_r
  spruce_birch <- spruce_r - birch_r
  alder_oak <- alder_r - oak_r
  pine_oak <- pine_r - oak_r
  spruce_oak <- spruce_r - oak_r
  pine_alder <- pine_r - alder_r
  spruce_alder <- spruce_r - alder_r
  spruce_pine <- spruce_r - pine_r
   
}


