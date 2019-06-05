## lichen lto model
##
## First edit: 20190319
## Last edit: 20190603
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(k in 1:nspecies){
    for(p in 1:nplot){
      for(i in 1:ntree[p]){
        obs[i,k,p] ~ dbern(p_occ[i,k,p])
        logit(p_occ[i,k,p]) <- alpha_plot_mean[k] + plot_effect[k,p] + 
                               beta_dbh[k]*stem_dbh[i,p] +
                               beta_td[k]*td[p,1] +
                               beta_td_quad[k]*td[p,1]^2
      }
    }
  }
  
  ## Plot level:
  for(k in 1:nspecies){
    for(p in 1:nplot){
      plot_effect[k,p] ~ dnorm(0, tau_plot[k])
      # mu[k,p] <- alpha_plot_mean[k] +
                 # beta_ud[k]*ud[p,1] +
                 # beta_ud_quad[k]*ud[p,1]^2 +
                 # beta_cd[k]*cd[p,1] +
                 # beta_cd_quad[k]*cd[p,1]^2 +
                 # beta_td[k]*td[p,1] +
                 # beta_td_quad[k]*td[p,1]^2
    }
  }
  
  ## Priors:
  
  for(k in 1:nspecies){
    beta_dbh[k] ~ dnorm(0, 0.001)
    alpha_plot_mean[k] ~ dnorm(0, 0.001)
    # beta_ud[k] ~ dnorm(0, 0.001)
    # beta_ud_quad[k] ~ dnorm(0, 0.001)
    # beta_cd[k] ~ dnorm(0, 0.001)
    # beta_cd_quad[k] ~ dnorm(0, 0.001)
    beta_td[k] ~ dnorm(0, 0.001)
    beta_td_quad[k] ~ dnorm(0, 0.001)
    tau_plot[k] <- 1/sigma_plot[k]^2
    sigma_plot[k] ~ dgamma(0.001, 0.001) T(0.0001, 50)
  }
  
  ## Model validation:

  ## Predictions:
  
  for(k in 1:nspecies){
    # ud_max[k] <- -beta_ud[k]/(2*beta_ud_quad[k])
    # cd_max[k] <- -beta_cd[k]/(2*beta_cd_quad[k])
    td_max[k] <- -beta_td[k]/(2*beta_td_quad[k])
  }
  
}


