## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We put explanatory variables on both the beta and gamma diversity 
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

## Calculate median and 95% CI for data.table:
quant_calc <- function(x) as.list(quantile(x, probs = c(0.025, 0.5, 0.975)))

## 3. Load and explore data ----------------------------------------------------

dir("clean")

load("clean/species_accumulation_data.rda")
load("clean/sad_tree_part.rda")
str(sad)
str(sad_tree)

## 4. Prepare data -------------------------------------------------------------

## Calculate the number of trees by plot:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

## Calculate mean accumulation curve across all 100 permutations:
obs <- apply(sad, c(1,3), mean)

## Create model data set:
data <- list(nplot = dim(sad)[3],
             ntree = ntree,
             obs = obs,
             u = 1, ## Evaluation unit of alpha diversity
             dec = scale(sad_tree$dec),
             spruce = scale(sad_tree$spruce),
             pine = scale(sad_tree$pine),
             tsp_2 = ifelse(sad_tree$nr_tsp == 2, 1, 0),
             tsp_3 = ifelse(sad_tree$nr_tsp == 3, 1, 0),
             tsp_4 = ifelse(sad_tree$nr_tsp == 4, 1, 0),
             dbh = scale(sad_tree$dbh)[, 1])

cor_list <- data[5:11]
cor_list$tsp1 <- ifelse(sad_tree$nr_tsp == 1, 1, 0)
cor_list$dec_quad <- data$dec^2
cor_list$pine_quad <- data$pine^2
cor_list$spruce_quad <- data$spruce^2

capture.output(cor(as.data.frame(cor_list))) %>% 
  write(., "results/cor_lpsac.txt")

str(data)

## Load the models:
m.perc <- "scripts/JAGS/lichen_JAGS_lpsac_perc.R"
m.tsp <- "scripts/JAGS/lichen_JAGS_lpsac_tsp.R"

start <- Sys.time()

n.adapt <- 1000; n.iter <- 1000; samples <- 5000; n.thin <- 5

## 5. Run m.perc ---------------------------------------------------------------

inits <- list(list(sd_obs = rep(1, data$nplot),
                   gdiv = rep(60, data$nplot), bdiv = rep(6, data$nplot),
                   sigma_gdiv = 1, 
                   g_icpt = 10, g_dbh = 1, g_perc = 1, g_perc2 = 1,
                   sigma_bdiv = 1, 
                   b_icpt = 1.1, b_dbh = 1, b_perc = 1, b_perc2 = 1),
              list(sd_obs = rep(0.5, data$nplot),
                   gdiv = rep(40, data$nplot), bdiv = rep(2, data$nplot),
                   sigma_gdiv = 5,
                   g_icpt = 33, g_dbh = -2, g_perc = 0, g_perc2 = -10,
                   sigma_bdiv = 0.5, 
                   b_icpt = 1.15, b_dbh = -0.1, b_perc = 0.1, b_perc2 = -0.1),
              list(sd_obs = rep(2, data$nplot),
                   gdiv = rep(20, data$nplot), bdiv = rep(3, data$nplot),
                   sigma_gdiv = 10,
                   g_icpt = 45, g_dbh = 5, g_perc = 10, g_perc2 = 0,
                   sigma_bdiv = 1.5, 
                   b_icpt = 1.5, b_dbh = 0.02, b_perc = 0.4, b_perc2 = 0))

perc_pred_export <- vector("list", 3)
names(perc_pred_export) <- c("dec", "spruce", "pine")
perc_max_export <- vector("list", 3)
names(perc_max_export) <- c("dec", "spruce", "pine")

for(i in c("dec", "spruce", "pine")){
  
  data$perc <- unlist(data[i])
  data$perc_pred <- seq(min(data$perc), max(data$perc), 0.05)
  
  cl <- makePSOCKcluster(3) 
  parJagsModel(cl, "lpsac", m.perc, data, inits, 3, n.adapt)
  parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)
  
  zc <- parCodaSamples(cl = cl, model = "lpsac",
                       variable.names = c("sd_obs",
                                          "sigma_gdiv",
                                          "g_icpt", "g_dbh",
                                          "g_perc", "g_perc2", 
                                          "sigma_bdiv", "b_icpt", "b_dbh",
                                          "b_perc", "b_perc2"),
                       n.iter = samples, thin = n.thin)
  
  capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
    write(., paste0("results/parameters_lpsac_", i, ".txt"))
  
  pdf(paste0("figures/lpsac_", i, ".pdf"))
  plot(zc); gelman.plot(zc) 
  dev.off()
  
  capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
    write(., paste0("results/diagnostics_lpsac_", i, ".txt"))
  
  ## Export predicted diversity metrics:
  zc_pred_perc <- parCodaSamples(cl = cl, model = "lpsac",
                                 variable.names = c("bdiv_pred", "gdiv_pred"),
                                 n.iter = samples, thin = n.thin)
  
  pred_perc <- as.data.frame(summary(zc_pred_perc)$quantile)
  pred_perc$perc_pred <- backscale(data$perc_pred, data[[i]])
  pred_perc$div_metric <- c("bdiv", "gdiv")
  pred_perc$div_metric <- sort(pred_perc$div_metric)
  
  perc_pred_export[[i]] <- pred_perc
  
  ## Export maximum percentage values for the quadratic curves of gdiv and bdiv:
  zc_pp_max <- parCodaSamples(cl = cl, model = "lpsac",
                              variable.names = c("b_perc_max", "g_perc_max"),
                              n.iter = samples, thin = n.thin)
  
  pp_max <- data.frame("b_max" = backscale(unlist(zc_pp_max[, "b_perc_max"]), 
                                           data[[i]]),
                       "g_max" = backscale(unlist(zc_pp_max[, "g_perc_max"]), 
                                           data[[i]]))
  
  perc_max_export[[i]] <- pp_max
  
  stopCluster(cl)
  
}

## Export for graphing:
pred_perc <- reshape2::melt(perc_pred_export,
                            measure.vars = colnames(perc_pred_export))
write.csv(pred_perc, paste0("clean/lpsac_collector_pred_perc.csv"))

## Export for graphing:
pp_max <- reshape2::melt(perc_max_export, 
                         measure.vars = colnames(perc_max_export))
write.csv(pp_max, "clean/lpsac_collector_max_perc.csv")

## 6. Run m.tsp ----------------------------------------------------------------

inits <- list(list(sd_obs = rep(1, data$nplot),
                   gdiv = rep(60, data$nplot), bdiv = rep(6, data$nplot),
                   sigma_gdiv = 1,
                   g_icpt = 30, g_dbh = 1, 
                   g_2tsp = 1, g_3tsp = 1, g_4tsp = 1,
                   sigma_bdiv = 1, 
                   b_icpt = 1, b_dbh = 1, 
                   b_2tsp = 1, b_3tsp = 1, b_4tsp = 1),
              list(sd_obs = rep(0.5, data$nplot),
                   gdiv = rep(40, data$nplot), bdiv = rep(2, data$nplot),
                   sigma_gdiv = 5,
                   g_icpt = 20, g_dbh = -2, 
                   g_2tsp = 0, g_3tsp = 0, g_4tsp = 5,
                   sigma_bdiv = 0.5, 
                   b_icpt = 0.6, b_dbh = -0.1, 
                   b_2tsp = -0.4, b_3tsp = 0, b_4tsp = 0.2),
              list(sd_obs = rep(2, data$nplot),
                   gdiv = rep(20, data$nplot), bdiv = rep(3, data$nplot),
                   sigma_gdiv = 10,
                   g_icpt = 35, g_dbh = 5, 
                   g_2tsp = 15, g_3tsp = 20, g_4tsp = 25,
                   sigma_bdiv = 1.2, 
                   b_icpt = 1.4, b_dbh = 0.5, 
                   b_2tsp = 0.5, b_3tsp = 0.8, b_4tsp = 0.8))

cl <- makePSOCKcluster(3) 
parJagsModel(cl, "lpsac", m.tsp, data, inits, 3, n.adapt)
parUpdate(cl = cl, object = "lpsac", n.iter = n.iter)

zc <- parCodaSamples(cl = cl, model = "lpsac",
                     variable.names = c("sd_obs",
                                        "sigma_gdiv", "g_icpt", "g_dbh",
                                        "g_2tsp", "g_3tsp", "g_4tsp", 
                                        "sigma_bdiv", "b_icpt", "b_dbh",
                                        "b_2tsp", "b_3tsp", "b_4tsp"),
                     n.iter = samples*2, thin = n.thin*2)

capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_tsp.txt")

pdf("figures/lpsac_tsp.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_tsp.txt")

zc_pred_tsp <- parCodaSamples(cl = cl, model = "lpsac",
                              variable.names = c("bdiv_1tsp", "bdiv_2tsp", 
                                                 "bdiv_3tsp", "bdiv_4tsp",
                                                 "gdiv_1tsp", "gdiv_2tsp", 
                                                 "gdiv_3tsp", "gdiv_4tsp"),
                              n.iter = samples*2, thin = n.thin*2)

pred_tsp <- as.data.frame(summary(zc_pred_tsp)$quantiles)
pred_tsp$div_metric <- c("bdiv", "gdiv")
pred_tsp$div_metric <- sort(pred_tsp$div_metric)
pred_tsp$nr_tsp <- rep(1:4, 2)

write.csv(pred_tsp, paste0("clean/lpsac_pred_tsp.csv"))

zc_tsp_diff <- parCodaSamples(cl = cl, model = "lpsac",
                              variable.names = c("bdiv_diff_21", "bdiv_diff_31", 
                                                 "bdiv_diff_41", "bdiv_diff_32",
                                                 "bdiv_diff_42", "bdiv_diff_43",
                                                 "gdiv_diff_21", "gdiv_diff_31", 
                                                 "gdiv_diff_41", "gdiv_diff_32",
                                                 "gdiv_diff_42", "gdiv_diff_43"),
                              n.iter = samples*2, thin = n.thin*2)

## Extract probability that difference between nr_tsp is bigger than 0:
ANOVA_prob <- summary(zc_tsp_diff)$quantiles
ANOVA_prob <- cbind(ANOVA_prob, 
                    "ecdf" = sapply(as.data.frame(combine.mcmc(zc_tsp_diff)), 
                                    function(x) 1-ecdf(x)(0))) 
write.csv(ANOVA_prob, "results/lpsac_tsp_ANOVA.csv")

## Use m.tsp to export diversity metrics for all sites to make ternary plots:

zc_pred <- parCodaSamples(cl = cl, model = "lpsac",
                          variable.names = c("bdiv", "gdiv"),
                          n.iter = samples,
                          thin = n.thin)
pred1 <- as.data.frame(summary(zc_pred)$quantiles)
pred1$mean <- summary(zc_pred)$statistics[, "Mean"]
pred1$div_metric <- gsub("[[[:digit:]]+]", "", rownames(pred1))
pred1$site <- 1:data$nplot
pred1 <- cbind(pred1, sad_tree[, c("dec", "spruce", "pine")])

## Export for graphing:
write.csv(pred1, "clean/lpsac_site_pred.csv")

## Produce SAC from fitted model and compare with raw data:

zc_val <- parCodaSamples(cl = cl, model = "lpsac",
                         variable.names = "obs_pred",
                         n.iter = samples/2,
                         thin = n.thin/2)

pred2 <- summary(zc_val, quantiles = c(0.025, 0.5, 0.975))$quantiles
pred2 <- cbind(pred2, "ntree" = unlist(lapply(data$ntree, function(x) 1:x)))
plot <- list() ; for(i in 1:data$nplot){plot[[i]] <- rep(i, ntree[i])}
pred2 <- as.data.frame(cbind(pred2, "plot" = unlist(plot)))

dev.off()

pdf("figures/sim_vs_obs.pdf")
par(mfrow = c(3, 2))
for(i in 1:data$nplot) {
  x = c(0, pred2[pred2$plot == i, "ntree"])
  y = pred2[pred2$plot == i, ]
  plot(x, c(0, y$`97.5%`), 
       xlab = "tree nr", ylab = "richness",
       lty = "dashed", col = "blue", typ = "l")
  lines(x, c(0, y$`50%`), col = "blue")
  lines(x, c(0, y$`2.5%`), lty = "dashed", col = "blue")
  ## Real data:
  points(1:(length(x)-1), na.omit(as.vector(obs[,i])))
}
dev.off()

stopCluster(cl)

end <- Sys.time()
end-start

## -------------------------------END-------------------------------------------
