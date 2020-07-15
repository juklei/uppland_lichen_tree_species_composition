## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We put explanatory variables on both the half saturation and the asymptote 
##
## First edit: 20190605
## Last edit: 20200710
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(runjags)
library(coda)
library(magrittr)
library(reshape2)
library(data.table)
library(parallel)
library(dclone)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

load("clean/species_accumulation_data.rda")
load("clean/sad_tree_part.rda")
str(sad)
str(sad_tree)

## 4. Prepare data and inits ---------------------------------------------------

## Calculate the number of trees by plot:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

obs <- apply(sad, c(1,3), mean)

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nplot = dim(sad)[3],
             ntree = ntree,
             obs = obs,
             dec = scale(sad_tree$dec),
             spruce = scale(sad_tree$spruce),
             pine = scale(sad_tree$pine),
             tsp_2 = ifelse(sad_tree$nr_tsp == 2, 1, 0),
             tsp_3 = ifelse(sad_tree$nr_tsp == 3, 1, 0),
             tsp_4 = ifelse(sad_tree$nr_tsp == 4, 1, 0),
             dbh = scale(sad_tree$dbh))

cor_list <- data[5:11]
cor_list$tsp1 <- ifelse(sad_tree$nr_tsp == 1, 1, 0)
cor_list$dec_quad <- data$dec^2
cor_list$pine_quad <- data$pine^2
cor_list$spruce_quad <- data$spruce^2

capture.output(cor(as.data.frame(cor_list))) %>% 
  write(., "results/cor_lpsac.txt")

## Add prediction data:
data$dec_pred <- seq(min(data$dec), max(data$dec), 0.05)
data$spruce_pred <- seq(min(data$spruce), max(data$spruce), 0.05)
data$pine_pred <- seq(min(data$pine), max(data$pine), 0.05)

str(data)

inits <- list(list(sd_obs = rep(1, data$nplot),
                   r_alpha = 10, r_beta_dbh = 1,
                   r_beta_dec = 1, r_beta2_dec = 1,
                   r_beta_spruce = 1, r_beta2_spruce = 1,
                   r_beta_pine = 1, r_beta2_pine = 1,
                   r_beta_2tsp = 1, r_beta_3tsp = 1, r_beta_4tsp = 1,
                   sigma_rich = 2, sigma_sat = 1, 
                   s_alpha = 3, s_beta_dbh = -0.2,
                   s_beta_dec = 0.1, s_beta2_dec = 0.1,
                   s_beta_spruce = 0.1, s_beta2_spruce = 0.1,
                   s_beta_pine = 0.1, s_beta2_pine = 0.1, 
                   s_beta_2tsp = 0, s_beta_3tsp = 0, s_beta_4tsp = 0),
              list(sd_obs = rep(1, data$nplot),
                   r_alpha = 18, r_beta_dbh = -1,
                   r_beta_dec = -1, r_beta2_dec = -1,
                   r_beta_spruce = -1, r_beta2_spruce = -1,
                   r_beta_pine = -1, r_beta2_pine = -1,
                   r_beta_2tsp = -1, r_beta_3tsp = -1, r_beta_4tsp = -1,
                   sigma_rich = 10, sigma_sat = 5, 
                   s_alpha = 1, s_beta_dbh = -0.2,
                   s_beta_dec = -1, s_beta2_dec = -1,
                   s_beta_spruce = -1, s_beta2_spruce = -1,
                   s_beta_pine = -1, s_beta2_pine = -1, 
                   s_beta_2tsp = 1, s_beta_3tsp = 1, s_beta_4tsp = 1),
              list(sd_obs = rep(1, data$nplot),
                   r_alpha = 40, r_beta_dbh = 0,
                   r_beta_dec = 0.01, r_beta2_dec = 0.1,
                   r_beta_spruce = 0.01, r_beta2_spruce = 0.1,
                   r_beta_pine = 0.01, r_beta2_pine = 0.1,
                   r_beta_2tsp = 0.01, r_beta_3tsp = 0.01, r_beta_4tsp = 0.01,
                   sigma_rich = 0.1, sigma_sat = 3, 
                   s_alpha = -1, s_beta_dbh = 0.01,
                   s_beta_dec = 0.01, s_beta2_dec = 0.01,
                   s_beta_spruce = 0.01, s_beta2_spruce = 0.01,
                   s_beta_pine = 0.01, s_beta2_pine = 0.01, 
                   s_beta_2tsp = 1, s_beta_3tsp = 1, s_beta_4tsp = 1))

## Load the models:
m.raw <- "scripts/JAGS/lichen_JAGS_lpsac_raw_new.R"
m.dec <- "scripts/JAGS/lichen_JAGS_lpsac_dec_new.R"
m.spruce <- "scripts/JAGS/lichen_JAGS_lpsac_spruce_new.R"
m.pine <- "scripts/JAGS/lichen_JAGS_lpsac_pine_new.R"
m.tsp <- "scripts/JAGS/lichen_JAGS_lpsac_tsp_new.R"

start <- Sys.time()

n.adapt <- 5000; n.iter <- 5000; samples <- 2500; n.thin <- 5

## 5. Run m.raw ----------------------------------------------------------------

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

## Load model:
parJagsModel(cl = cl, 
             name = "lpsac", 
             file = m.raw,
             data = data,
             inits = inits,
             n.chains = 3,
             n.adapt = n.adapt) 

## Update:
parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)

## Extract samples from estimated parameters:
zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("mu_rich", "sigma_rich",
                                        "mu_sat", "sigma_sat"),
                     n.iter = samples, 
                     thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_raw.txt")

## Test model convergence and mixing:
pdf("figures/lpsac_raw.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

## And with a table:
capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_raw.txt")

## Produce SAC from fitted model and compare with raw data:

zc_val <- parCodaSamples(cl = cl, model = "lpsac",
                         variable.names = "obs_pred",
                         n.iter = samples,
                         thin = n.thin)

pred <- summary(zc_val, quantiles = c(0.025, 0.5, 0.975))$quantiles
pred <- cbind(pred, "ntree" = unlist(lapply(data$ntree, function(x) 1:x)))
plot <- list() ; for(i in 1:data$nplot){plot[[i]] <- rep(i, ntree[i])}
pred <- as.data.frame(cbind(pred, "plot" = unlist(plot)))

dev.off()

pdf("figures/sim_vs_obs_new.pdf")
par(mfrow = c(3, 2))
for(i in 1:data$nplot) {
  x = c(0, pred[pred$plot == i, "ntree"])
  y = pred[pred$plot == i, ]
  plot(x, c(0, y$`97.5%`), 
       xlab = "tree nr", ylab = "richness",
       lty = "dashed", col = "blue", typ = "l")
  lines(x, c(0, y$`50%`), col = "blue")
  lines(x, c(0, y$`2.5%`), lty = "dashed", col = "blue")
  ## Real data:
  points(1:(length(x)-1), na.omit(as.vector(obs[,i])))
}
dev.off()

## Export estimated plot_richness and sat_speed for ternary plots:

zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("plot_richness", "sat_speed"),
                          n.iter = samples, 
                          thin = n.thin)

plot_estimates <- zj_pred
plot_estimates$dec <- sad_tree$dec
plot_estimates$spruce <- sad_tree$spruce
plot_estimates$pine <- sad_tree$pine
save(plot_estimates, file = "clean/plot_estimates.rda")

stopCluster(cl)

## 6. Run m.dec ----------------------------------------------------------------

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "lpsac", m.dec, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("sigma_rich", "r_alpha", "r_beta_dbh",
                                        "r_beta_dec", "r_beta2_dec",
                                        "sigma_sat", "s_alpha", "s_beta_dbh",
                                        "s_beta_dec", "s_beta2_dec"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_dec.txt")

pdf("figures/lpsac_dec.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_dec.txt")

zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("r_dec", "r_dec_max",
                                             "s_dec", "s_dec_max"),
                          n.iter = samples, thin = n.thin)

ld <- length(data$dec_pred)
export_dec <- list("r_dec" = zj_pred[, 1:ld],
                   "s_dec" = zj_pred[, (ld+2):(2*ld+1)],  ## TEST!
                   "r_dec_max" = backscale(unlist(zj_pred[, "r_dec_max"]), data$dec),
                   "s_dec_max" = backscale(unlist(zj_pred[, "s_dec_max"]), data$dec),
                   "dec_pred" = backscale(data$dec_pred, data$dec))
save(export_dec, file = "clean/export_dec.rda")

stopCluster(cl)

## 7. Run m.spruce -------------------------------------------------------------

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "lpsac", m.spruce, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("sigma_rich", "r_alpha", "r_beta_dbh",
                                        "r_beta_spruce", "r_beta2_spruce",
                                        "sigma_sat", "s_alpha", "s_beta_dbh",
                                        "s_beta_spruce", "s_beta2_spruce"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_spruce.txt")

pdf("figures/lpsac_spruce.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_spruce.txt")

zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("r_spruce", "r_spruce_max",
                                             "s_spruce", "s_spruce_max"),
                          n.iter = samples, thin = n.thin)

ls <- length(data$spruce_pred)
export_spruce <- list("r_spruce" = zj_pred[, 1:ls],
                      "s_spruce" = zj_pred[, (ls+2):(2*ls+1)],
                      "r_spruce_max" = backscale(unlist(zj_pred[, "r_spruce_max"]), data$spruce),
                      "s_spruce_max" = backscale(unlist(zj_pred[, "s_spruce_max"]), data$spruce),
                      "spruce_pred" = backscale(data$spruce_pred, data$spruce))
save(export_spruce, file = "clean/export_spruce.rda")

stopCluster(cl)

## 8. Run m.pine ---------------------------------------------------------------

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "lpsac", m.pine, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("sigma_rich", "r_alpha", "r_beta_dbh",
                                        "r_beta_pine", "r_beta2_pine",
                                        "sigma_sat", "s_alpha", "s_beta_dbh",
                                        "s_beta_pine", "s_beta2_pine"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_pine.txt")

pdf("figures/lpsac_pine.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_pine.txt")

zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("r_pine", "r_pine_max",
                                             "s_pine", "s_pine_max"),
                          n.iter = samples, thin = n.thin)

lp <- length(data$pine_pred)
export_pine <- list("r_pine" = zj_pred[, 1:lp],
                    "s_pine" = zj_pred[, (lp+2):(2*lp+1)],
                    "r_pine_max" = backscale(unlist(zj_pred[, "r_pine_max"]), data$pine),
                    "s_pine_max" = backscale(unlist(zj_pred[, "s_pine_max"]), data$pine),
                    "pine_pred" = backscale(data$pine_pred, data$pine))
save(export_pine, file = "clean/export_pine.rda")

stopCluster(cl)

## 9. Run m.tsp ----------------------------------------------------------------

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "lpsac", m.tsp, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("sigma_rich", "r_alpha", "r_beta_dbh",
                                        "r_beta_2tsp", "r_beta_3tsp", "r_beta_4tsp",
                                        "sigma_sat", "s_alpha", "s_beta_dbh",
                                        "s_beta_2tsp", "s_beta_3tsp", "s_beta_4tsp"),
                     n.iter = samples, thin = n.thin)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_tsp.txt")

pdf("figures/lpsac_tsp.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_tsp.txt")

zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("r_1tsp", "r_2tsp", "r_3tsp", "r_4tsp", 
                                             "r_diff_21", "r_diff_31", 
                                             "r_diff_41", "r_diff_32",
                                             "r_diff_42", "r_diff_43",
                                             "s_1tsp", "s_2tsp", "s_3tsp", "s_4tsp", 
                                             "s_diff_21", "s_diff_31", 
                                             "s_diff_41", "s_diff_32",
                                             "s_diff_42", "s_diff_43"),
                          n.iter = samples, thin = n.thin)

## Extract probability that difference between nr_tsp is bigger than 0:
tsp_diff <- combine.mcmc(zj_pred)[, c(5:10, 15:20)]
ANOVA_prob <- summary(tsp_diff)$quantiles[, c("50%", "2.5%", "97.5%")]
ANOVA_prob <- cbind(ANOVA_prob, 
                    "ecdf" = sapply(as.data.frame(tsp_diff), 
                                    function(x) 1-ecdf(x)(0))) 
write.csv(ANOVA_prob, "results/ANOVA_results_lpsac.csv")

export_tsp <- zj_pred[, c("r_1tsp", "r_2tsp", "r_3tsp", "r_4tsp",
                          "s_1tsp", "s_2tsp", "s_3tsp", "s_4tsp")]
save(export_tsp, file = "clean/export_tsp.rda")

stopCluster(cl)

end <- Sys.time()
end-start

## -------------------------------END-------------------------------------------
