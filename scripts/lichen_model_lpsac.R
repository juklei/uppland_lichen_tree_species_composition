## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We put an explanatory variable on both the half saturation and the asymptote 
##
## First edit: 20190605
## Last edit: 20190904
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

## 4. The model ----------------------------------------------------------------

## Calculate the number of trees by plot:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nplot = dim(sad)[3],
             ntree = ntree,
             obs = sad,
             dec = scale(sad_tree$dec),
             spruce = scale(sad_tree$spruce),
             pine = scale(sad_tree$pine),
             tsp_2 = ifelse(sad_tree$nr_tsp == 2, 1, 0),
             tsp_3 = ifelse(sad_tree$nr_tsp == 3, 1, 0),
             tsp_4 = ifelse(sad_tree$nr_tsp == 4, 1, 0),
             dbh = scale(sad_tree$dbh))

capture.output(cor(as.data.frame(data[5:11]))) %>% 
                 write(., "results/cor_lpsac.txt")

## Add prediction data:
data$dec_pred <- seq(min(data$dec), max(data$dec), 0.05)
data$spruce_pred <- seq(min(data$spruce), max(data$spruce), 0.05)
data$pine_pred <- seq(min(data$pine), max(data$pine), 0.05)

str(data)

inits <- list(list(plot_richness = sample(1:50, data$nplot, replace = TRUE),
                   sat_speed = sample(1:10, data$nplot, replace = TRUE),
                   sigma_rich = 2,
                   sigma_sat = 2,
                   alpha = 29,
                   mu_sat = 5,
                   beta_dec = 1,
                   beta2_dec = 1,
                   beta_spruce = 0.5,
                   beta2_spruce = 0.5,
                   beta_pine = 1,
                   beta2_pine = 1,
                   beta_dbh = 0.2,
                   beta_2tsp = 0.1,
                   beta_3tsp = 0.1,
                   beta_4tsp = 0.1),
              list(plot_richness = sample(1:50, data$nplot, replace = TRUE),
                   sat_speed = sample(1:10, data$nplot, replace = TRUE),
                   sigma_rich = 5,
                   sigma_sat = 1,
                   alpha = 40,
                   mu_sat = 7,
                   beta_dec = -1,
                   beta2_dec = -1,
                   beta_spruce = -0.5,
                   beta2_spruce = -0.5,
                   beta_pine = -1,
                   beta2_pine = -1,
                   beta_dbh = -0.2,
                   beta_2tsp = -0.1,
                   beta_3tsp = -0.1,
                   beta_4tsp = -0.1),
              list(plot_richness = sample(1:50, data$nplot, replace = TRUE),
                   sat_speed = sample(1:10, data$nplot, replace = TRUE),
                   sigma_rich = 1,
                   sigma_sat = 1,
                   alpha = 15,
                   mu_sat = 1,
                   beta_dec = 0.1,
                   beta2_dec = 0.1,
                   beta_spruce = 0.1,
                   beta2_spruce = 0.1,
                   beta_pine = 0.1,
                   beta2_pine = 0.1,
                   beta_dbh = 0,
                   beta_2tsp = 0,
                   beta_3tsp = 0,
                   beta_4tsp = 0))

model <- "scripts/JAGS/lichen_JAGS_lpsac.R"

start <- Sys.time()

samples <- 1000; n.thin <- 10

## Model selection with one chain:

# jm <- jags.model(model, data, n.adapt = 2500, inits = inits[[1]], n.chains = 1)
# update(jm, n.iter = 2500)
# zj_modsel <- jags.samples(jm,
#                           variable.names = c("plot_richness", "pr_sim"),
#                           n.iter = samples,
#                           thin = n.thin)
# 
# PPL_quad_dec = sum((summary(zj_modsel$plot_richness, mean)$stat -
#                       summary(zj_modsel$pr_sim, mean)$stat)^2) +
#                sum((summary(zj_modsel$pr_sim, sd)$stat)^2)
# PPL_lin_dec = sum((summary(zj_modsel$plot_richness, mean)$stat - 
#                       summary(zj_modsel$pr_sim, mean)$stat)^2) + 
#               sum((summary(zj_modsel$pr_sim, sd)$stat)^2)

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

parJagsModel(cl = cl, 
             name = "lpsac", 
             file = model,
             data = data,
             n.adapt = 5000, 
             inits = inits, 
             n.chains = 3) 

parUpdate(cl = cl, object = "lpsac", n.iter = 1000)

zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("sigma_sat",
                                        "sigma_rich",
                                        "alpha",
                                        "mu_sat",
                                        "beta_dec",
                                        "beta2_dec",
                                        "dec_max",
                                        "diff_100vs0_dec",
                                        "beta_spruce",
                                        "beta2_spruce",
                                        "spruce_max",
                                        "diff_100vs0_spruce",
                                        "beta_pine",
                                        "beta2_pine",
                                        "pine_max",
                                        "diff_100vs0_pine",
                                        "beta_dbh",
                                        "beta_2tsp",
                                        "beta_3tsp",
                                        "beta_4tsp",
                                        "diff_21",
                                        "diff_31",
                                        "diff_41",
                                        "diff_32",
                                        "diff_42",
                                        "diff_43"),
                     n.iter = samples, 
                     thin = n.thin)

end <- Sys.time()
end-start

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_nr_tsp_normal_sat_model.txt")

## Extract probability that difference is bigger than 0:
ANOVA_prob <- data.frame("cat" = colnames(zc[[1]])[6:11],
                         "prob" = rep(NA, 6))
for(i in 1:6){
  ANOVA_prob$prob[i] <- 1-ecdf(unlist(zc[, paste(ANOVA_prob$cat[i])]))(0)
}

## Print and export:
capture.output(print("P(diff >= 0)"), as.matrix(ANOVA_prob)) %>%
  write(., "results/ANOVA_results_lpsac.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lpsac_nr_tsp_normal_sat_model.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_nr_tsp_normal_sat_model.txt")

# ## Produce validation metrics:
# zj_val <- jags.samples(jm,
#                        variable.names = c("obs_pred"),
#                        n.iter = samples,
#                        thin = n.thin)
# 
# pred <- summary(zj_val$obs_pred, quantile, c(.025,.5,.975))$stat
# x = 0:50
# 
# dev.off()
# 
# pdf("figures/sim_vs_obs_dec.pdf")
# 
# par(mfrow = c(3, 2))
# 
# for(i in 1:data$nplot) {
#   
#   y=pred[,,i]
#   plot(x,y[3,], lty="dashed", col="blue", xlab="tree nr", ylab="richness", typ="l")
#   lines(x,y[2,], col="blue")
#   lines(x,y[1,], lty="dashed", col="blue")
#   polygon(c(x,rev(x)), c(y[1,], rev(y[3,])), density=19, col="blue", angle=45)
#   
#   ## Real data:
#   points(rep(which(!is.na(sad[,1,i])), dim(sad)[2]),
#          na.omit(as.vector(sad[,,i])))
#   
# }
# 
# dev.off()

## 6. Produce and export figures -----------------------------------------------

zj_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("r_dec", 
                                             "r_spruce", 
                                             "r_pine", 
                                             "dec_max", 
                                             "spruce_max", 
                                             "pine_max", 
                                             "r_1tsp",
                                             "r_2tsp",
                                             "r_3tsp",
                                             "r_4tsp"),
                          n.iter = 2000, 
                          thin = 5)

# export_srd <- list("r_dec" = zj_pred,
#                    "dec_pred" = backscale(data$dec_pred, data$dec),
#                    "dec_max" = backscale(summary(zj_pred)$quantiles[length(data$dec_pred) + 1, ],
#                                          data$dec)*100)
# save(export_srd, file = "clean/sac_pred_r_dec.rda")

# export_srs <- list("r_spruce" = zj_pred,
#                    "spruce_pred" = backscale(data$spruce_pred, data$spruce),
#                    "spruce_max" = backscale(summary(zj_pred)$quantiles[length(data$spruce_pred) + 1, ], 
#                                             data$spruce)*100)
# save(export_srs, file = "clean/sac_pred_r_spruce.rda")

# export_srp <- list("r_pine" = zj_pred,
#                    "pine_pred" = backscale(data$pine_pred, data$pine),
#                    "pine_max" = backscale(summary(zj_pred)$quantiles[length(data$pine_pred) + 1, ],
#                                           data$pine)*100)
# save(export_srp, file = "clean/sac_pred_r_pine.rda")

export_srntsp <- zj_pred
save(export_srntsp, file = "clean/sac_pred_r_nr_tsp.rda")

stopCluster(cl)

## -------------------------------END-------------------------------------------
