## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We put an explanatory variable on both the half saturation and the asymptote 
##
## First edit: 20190605
## Last edit: 20190612
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)
library(reshape2)
library(data.table)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("results/ltsac")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Coefficient of variation:
CV <- function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)

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

capture.output(cor(sad_tree[, 2:6])) %>% write(., "results/cor_lpsac.txt")

## The variance decreases with the tree number looked at, e.g. at the asymptote
## the variance is 0. Here we look at sd, but in the JAGS module tau is used 
## instead. This is why the prior on beta_tau is strictly positive, even though 
## the trend in the graph below is negative.
T1 <- apply(sad, c(1, 3), sd, na.rm = TRUE)
T1 <- cbind(as.vector(T1), rep(1:dim(sad)[1], dim(sad)[3]))
T1 <- na.omit(T1)
plot(T1[, 2], T1[, 1], xlab = "nr.tree", ylab = "sd")
abline(lm(T1[, 1] ~ T1[, 2]))

## 4. The model ----------------------------------------------------------------

## Calculate the number of trees by plot:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nplot = dim(sad)[3],
             ntree = ntree,
             obs = sad,
             dec = scale(sad_tree$dec),
             nr_tsp = scale(sad_tree$nr_tsp),
             dbh = scale(sad_tree$dbh))

## Add prediction data:

## Percent deciduous:
data$dec_pred <- seq(min(data$dec), max(data$dec), 0.05)

## Nr. of tree species:
data$dec_pred <- seq(min(data$dec), max(data$dec), 0.05)

str(data)

## Prepare inits:

inits <- list(list(plot_richness = rep(20, data$nplot),
                   tau_obs = rep(1, dim(sad)[1]),
                   alpha_tau = 3,
                   beta_tau = 0.5,
                   sigma_tau = 0.5,
                   sat_speed = rep(5, data$nplot),
                   alpha_rich = 0,
                   alpha_sat = 0,
                   beta_dec_rich = 0.2,
                   beta_dec_sat = 0.1,
                   beta_dbh_rich = 0.2,
                   beta_dbh_sat = 0.1))

model <- "scripts/JAGS/lichen_JAGS_lpsac.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  5000

update(jm, n.iter = burn.in) 

samples <- 10000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("alpha_tau", 
                                      "beta_tau",
                                      "sigma_tau",
                                      "alpha_rich",
                                      "alpha_sat",
                                      "beta_dec_rich",
                                      "beta_dec_sat",
                                      "beta_dbh_rich",
                                      "beta_dbh_sat"),
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lpsac.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac.txt")

## Produce validation metrics:
zj_val <- jags.samples(jm,
                       variable.names = c("obs_sim"),
                       n.iter = 1000,
                       thin = 10)

## Plot the accumulation data per plot (plot in third dimesnion in array):

## Extract the mean of the simulated values:

## For the spread:
dsp <- zj_val$obs_sim[,1,,,1]

## For the mean:
dsp <- apply(zj_val$obs_sim, c(1,3,4), mean, na.rm = TRUE)


dev.off()

pdf("figures/sim_vs_obs.pdf")

par(mfrow = c(3, 2))

for(i in 1:data$nplot) {

## Sim data:
plot(rep(which(!is.na(dsp[,i,1])), dim(dsp)[3]), 
     na.omit(as.vector(dsp[,i,])), 
     col = "red", 
     xlab = "tree", 
     ylab = "richness")

## Real data:
points(rep(which(!is.na(sad[,1,i])), dim(sad)[2]), 
       na.omit(as.vector(sad[,,i])))

}

dev.off()

## 6. Produce and export figures -----------------------------------------------

zj_pred <- jags.samples(jm,
                        variable.names = c("r_dec"),
                        n.iter = samples,
                        thin = n.thin)

## -------------------------------END-------------------------------------------
