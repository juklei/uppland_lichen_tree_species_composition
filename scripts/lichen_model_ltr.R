## model lichen richness in a hierarchical model with stem diameter and tree
## species at the tree level
##
## First edit: 20190125
## Last edit: 20190612
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)

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

load("clean/ltr_with_tree_data.rda")
head(ltr_tree)
str(ltr_tree)

ltr <- na.omit(ltr_tree)

## 4. The model ----------------------------------------------------------------

## Create model data set:
data <- list(nobs = nrow(ltr),
             plot = as.numeric(ltr$plot),
             nplot = length(unique(ltr$plot)),
             richness = ltr$richness,
             pine = ifelse(ltr$Tree.species == "Ps", 1, 0),
             spruce = ifelse(ltr$Tree.species == "Pa", 1, 0),
             aspen = ifelse(ltr$Tree.species == "Pt", 1, 0),
             oak = ifelse(ltr$Tree.species == "Qr", 1, 0),
             alder = ifelse(ltr$Tree.species == "Ag", 1, 0),
             dbh = scale(ltr$Tree.diameter.130.cm.above.ground))

str(data)

inits <-  list(list(p = 0.8,
                    richness_true = rep(8, data$nobs),
                    beta_pine = 0.5,
                    beta_spruce = 0.5,
                    beta_aspen = 0.5,
                    beta_oak = 0.5,
                    beta_alder = 0.5,
                    beta_dbh = 0.5,
                    alpha = 2,
                    sigma_plot = 2),
               list(p = 0.5,
                    richness_true = rep(5, data$nobs),
                    beta_pine = 0,
                    beta_spruce = 0,
                    beta_aspen = 0,
                    beta_oak = 0,
                    beta_alder = 0,
                    beta_dbh = 0,
                    alpha = 3,
                    sigma_plot = 0.1),
               list(p = 0.99,
                    richness_true = rep(11, data$nobs),
                    beta_pine = 1,
                    beta_spruce = 1,
                    beta_aspen = 1,
                    beta_oak = 1,
                    beta_alder = 1,
                    beta_dbh = 1,
                    alpha = 5,
                    sigma_plot = 5))

model <- "scripts/JAGS/lichen_JAGS_ltr.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 10000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("p",
                                      "alpha",
                                      "beta_pine",
                                      "beta_spruce",
                                      "beta_aspen",
                                      "beta_oak",
                                      "beta_alder",
                                      "beta_dbh",
                                      "sigma_plot"),
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_ltr.txt")

zj_diff <- jags.samples(jm, 
                        variable.names = c("aspen_birch",
                                           "aspen_oak",
                                           "aspen_alder",
                                           "aspen_pine",
                                           "aspen_spruce",
                                           "birch_oak", 
                                           "birch_alder", 
                                           "birch_pine",
                                           "birch_spruce",
                                           "oak_alder",
                                           "oak_pine",
                                           "oak_spruce",
                                           "alder_pine",
                                           "alder_spruce",
                                           "pine_spruce"), 
                        n.iter = samples, 
                        thin = n.thin)

## Extract 95% CIs of the difference:
ANOVA_diff <- lapply(zj_diff, 
                     function(x) summary(x, quantile, c(.025,.975))$stat)

## Extract probability that difference is bigger than 0:
ANOVA_prob <- lapply(zj_diff, function(x) 1-ecdf(x)(0)) 

## Print and export:
capture.output(print("95% Confidence Intervals"),
               as.data.frame(ANOVA_diff),
               print("P(diff >= 0)"),
               as.matrix(ANOVA_prob)) %>%
  write(., "results/ANOVA_results_ltr.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_ltr.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc), gelman.diag(zc)) %>% 
  write(., "results/diagnostics_ltr.txt")

# ## Produce validation metrics: 
# zj_val <- jags.samples(jm, 
#                        variable.names = c("mean_richness", 
#                                           "mean_richness_sim",
#                                           "p_mean", 
#                                           "cv_richness", 
#                                           "cv_richness_sim", 
#                                           "p_cv", 
#                                           "fit", 
#                                           "fit_sim",
#                                           "p_fit"), 
#                        n.iter = samples, 
#                        thin = n.thin)
# 
# ## Fit of mean:
# plot(zj_val$mean_richness, 
#      zj_val$mean_richness_sim, 
#      xlab = "mean real", 
#      ylab = "mean simulated", 
#      cex = .05)
# abline(0, 1)
# p <- summary(zj_val$p_mean, mean)
# text(x = 8, y = 10.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv_richness, 
#      zj_val$cv_richness_sim, 
#      xlab = "cv real", 
#      ylab = "cv simulated", 
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_cv, mean)
# text(x = .25, y = .335, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit, 
#      zj_val$fit_sim, 
#      xlab = "ssq real", 
#      ylab = "ssq simulated", 
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_fit, mean)
# text(x = 480, y = 650, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("spruce_mean",
                                           "pine_mean",
                                           "aspen_mean",
                                           "oak_mean",
                                           "alder_mean",
                                           "birch_mean"),
                        n.iter = samples, 
                        thin = n.thin)

## 7. Export data for fancy figures --------------------------------------------

export_ltr <- zj_pred

save(export_ltr, file = "clean/ltr_pred.rdata")

## -------------------------------END-------------------------------------------
