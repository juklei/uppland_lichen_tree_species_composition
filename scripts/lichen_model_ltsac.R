## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We put an explanatory variable on both the half saturation and the asymptote 
##
## First edit: 20190605
## Last edit: 20190605
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

lto <- read.csv("clean/lto_T_50.csv")
head(lto)
str(lto)
table(lto$species)

## 4. Create species accumulation data set:


## 4. The model ----------------------------------------------------------------

## Exclude NA in stem_dbh:
#lto <- lto[!is.na(lto$Tree.diameter.130.cm.above.ground), ]



## Maje a data array:
obs <- acast(lto[, c("species", "Tree.no", "plot", "observed")],
             formula = Tree.no ~ species ~ plot,
             value.var = "observed")
obs[is.na(obs)] <- 0


## Plot level explanatory variables need to be reduced to unique rows:
plu <- unique(lto[, c("plot", 
                      "nr_lov",
                      "nr_gran",
                      "nr_staende_dodved",
                      "nr_tall",
                      "average_dbh_lov",                    
                      "average_dbh_gran",
                      "average_dbh_staende_dodved",         
                      "average_dbh_tall",
                      "average_dbh_all_alive",
                      "PercentAbove5m",
                      "PercentBelow5m")])
nrow(plu)

## Create explanatory arrays:


stem_dbh <- acast(lto[, c("Tree.no", 
                          "plot", 
                          "Tree.diameter.130.cm.above.ground")],
                  formula = Tree.no ~ plot,
                  value.var = "Tree.diameter.130.cm.above.ground",
                  fun.aggregate = mean)

ntree <- unique(lto[, c("plot", "Tree.no")])
ntree <- as.data.table(ntree)
ntree <- ntree[, list("ntree" = max(Tree.no)), by = "plot"]

## Create model data set:
data <- list(nplot = length(unique(lto$plot)),
             nspecies = length(unique(lto$species)),
             ntree = ntree$ntree,
             obs = obs,
             stem_dbh = scale(stem_dbh),
             ud = scale(plu$PercentBelow5m),
             cd = scale(plu$PercentAbove5m),
             td = scale(plu$PercentBelow5m + plu$PercentAbove5m))

str(data)

## Prepare inits:

inits <-list(list(alpha_plot_mean = rep(0, data$nspecies),
                  beta_dbh = rep(0.5, data$nspecies),
                  beta_ud = rep(-0.5, data$nspecies),
                  beta_ud_quad = rep(-0.5, data$nspecies),
                  beta_cd = rep(-0.5, data$nspecies),
                  beta_cd_quad = rep(-0.5, data$nspecies),
                  beta_td = rep(-0.5, data$nspecies),
                  beta_td_quad = rep(-0.5, data$nspecies)
                  ))

model <- "scripts/JAGS/lichen_JAGS_lto.R"

T1 <- Sys.time()

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 50000
n.thin <- 25

zc <- coda.samples(jm,
                   variable.names = c("alpha_plot_mean",
                                      "tau_plot",
                                      "beta_dbh",
                                      "beta_ud",
                                      "beta_ud_quad",
                                      "beta_cd",
                                      "beta_cd_quad",
                                      "beta_td",
                                      "beta_td_quad"),
                   n.iter = samples, 
                   thin = n.thin)

T2 <- Sys.time()
T2-T1

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lichen.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lichen.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lichen.txt")

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("ud_max", "cd_max", "td_max"),
                        n.iter = samples, 
                        thin = n.thin)

ud_max <- summary(zj_pred$ud_max, quantile, c(.025,.5,.975))$stat
ud_max <- backscale(ud_max, data$ud)
cd_max <- summary(zj_pred$cd_max, quantile, c(.025,.5,.975))$stat
cd_max <- backscale(cd_max, data$cd)
td_max <- summary(zj_pred$td_max, quantile, c(.025,.5,.975))$stat
td_max <- backscale(td_max, data$td)

## Combine results:

result <- data.frame("species" = unique(lto$species),
                     "ud_max_median" = ud_max[2,],
                     "ud_max_2.5%" = ud_max[1,],
                     "ud_max_97.5%" = ud_max[3,],
                     "cd_max_median" = cd_max[2,],
                     "cd_max_2.5%" = cd_max[1,],
                     "cd_max_97.5%" = cd_max[3,])

result <- data.frame("species" = unique(lto$species),
                     "td_max_median" = td_max[2,],
                     "td_max_2.5%" = td_max[1,],
                     "td_max_97.5%" = td_max[3,])

## <0 and >100 replace by 0 resp. 100: 
result[result < 0] <- 0
result[result > 100] <- 100

write.csv(result, "results/optimal_vegetation densities.csv", row.names = FALSE)

## -------------------------------END-------------------------------------------
