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

str(data)

## Prepare inits:

inits <- list(list(plot_richness = rep(20, data$nplot),
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
                 n.adapt = 500, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  1000

update(jm, n.iter = burn.in) 

samples <- 50000
n.thin <- 25

zc <- coda.samples(jm,
                   variable.names = c("alpha_rich",
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

## 6. Produce and export figures -----------------------------------------------


## -------------------------------END-------------------------------------------
