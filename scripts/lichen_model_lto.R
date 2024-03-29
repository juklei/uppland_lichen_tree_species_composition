## model lichen occurence in a hierarchical model with stem diameter and tree
## species at the tree level. Takes 2.6 days with 2.6 GHz and 16GB RAM.
##
## First edit: 2020706
## Last edit: 20200731
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
library(gt)

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

lto <- read.csv("clean/lto.csv")
str(lto)
table(lto$species)

## 4. The model ----------------------------------------------------------------

## Exclude NA in dbh:
lto <- lto[!is.na(lto$DBH), ]

# ## Create table with nr.obswrved/non-observed per species:
# T1 <- table(lto[, c("species", "observed")])
# 
# ## Exclude species which where seen > x and not seen > x times:
# red_names <- names(which(T1[, "1"] > 50 & T1[, "0"] > 50))
# lto <- droplevels(lto[lto$species %in% red_names, ])

lto <- droplevels(lto)

## Create arrays:

obs <- acast(lto[, c("species", "Tree.no", "Plot.no.", "observed")],
             formula = Tree.no ~ species ~ Plot.no.,
             value.var = "observed")

## Create explanatory variables matrices:
dbh <- acast(lto,
             formula = Tree.no ~ Plot.no.,
             value.var = "DBH",
             fun.aggregate = mean)
tsp <- acast(unique(lto[, c("Tree.no", "Plot.no.", "tsp")]),
             formula = Tree.no ~ Plot.no.,
             value.var = "tsp")

## Create a vector showing how many trees there are per plot, to prevent JAGS
## from looping through too many trees per plot:
ntree <- unique(lto[, c("Plot.no.", "Tree.no")])
ntree <- as.data.table(ntree)
ntree <- ntree[, list("ntree" = max(Tree.no)), by = "Plot.no."]

## Create model data set:
data <- list(nplot = length(unique(lto$Plot.no.)),
             nspecies = length(unique(lto$species)),
             ntree = ntree$ntree,
             obs = obs,
             pine = matrix(ifelse(tsp == "Ps", 1, 0), nrow(tsp)),
             spruce = matrix(ifelse(tsp == "Pa", 1, 0), nrow(tsp)),
             aspen = matrix(ifelse(tsp == "Pt", 1, 0), nrow(tsp)),
             oak = matrix(ifelse(tsp == "Qr", 1, 0), nrow(tsp)),
             alder = matrix(ifelse(tsp == "Ag", 1, 0), nrow(tsp)),
             dbh = scale(dbh))

str(data)

## Prepare inits:

inits <-list(list(mu_alpha = 0, sd_alpha = 10,
                  mu_b_pine = 0, sd_b_pine = 1,
                  mu_b_spruce = 0, sd_b_spruce = 1,
                  mu_b_aspen = 0, sd_b_aspen = 1,
                  mu_b_oak = 0, sd_b_oak = 1,
                  mu_b_alder = 0, sd_b_alder = 1,
                  mu_b_dbh = 0, sd_b_dbh = 1,
                  u_sp = 1),
             list(mu_alpha = -8, sd_alpha = 2.5,
                  mu_b_pine = -3, sd_b_pine = 1,
                  mu_b_spruce = -1, sd_b_spruce = 2,
                  mu_b_aspen = -3, sd_b_aspen = 3,
                  mu_b_oak = -2, sd_b_oak = 1.5,
                  mu_b_alder = -1, sd_b_alder = 0.5,
                  mu_b_dbh = 0, sd_b_dbh = 0.05,
                  u_sp = 0.6),
             list(mu_alpha = -5, sd_alpha = 4,
                  mu_b_pine = 3, sd_b_pine = 3,
                  mu_b_spruce = 1.5, sd_b_spruce = 4,
                  mu_b_aspen = 3, sd_b_aspen = 7,
                  mu_b_oak = 2, sd_b_oak = 4.5,
                  mu_b_alder = 1.5, sd_b_alder = 2.5,
                  mu_b_dbh = 0.15, sd_b_dbh = 0.25,
                  u_sp = 1.5))

model <- "scripts/JAGS/lichen_JAGS_lto.R"

start <- Sys.time()

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

jm <- parJagsModel(cl = cl, 
                   name = "lto",
                   file = model,
                   data = data,
                   n.adapt = 5000, 
                   inits = inits,
                   n.chains = 3) 

parUpdate(cl = cl, object = "lto", n.iter = 5000)

end <- Sys.time()
end-start

samples <- 50000
n.thin <- 40

zc_1 <- parCodaSamples(cl = cl, model = "lto",
                       variable.names = c("mu_alpha", "sd_alpha",
                                          "mu_b_pine", "sd_b_pine",
                                          "mu_b_spruce", "sd_b_spruce",
                                          "mu_b_aspen", "sd_b_aspen",
                                          "mu_b_oak", "sd_b_oak",
                                          "mu_b_alder", "sd_b_alder",
                                          "mu_b_dbh", "sd_b_dbh",
                                          "u_sp"),
                       n.iter = samples, 
                       thin = n.thin)

end <- Sys.time()
end-start

## Export parameter estimates:
capture.output(summary(zc_1), HPDinterval(zc_1, prob = 0.95)) %>% 
  write(., "results/parameters_lto.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lto.pdf")
plot(zc_1); gelman.plot(zc_1) 
dev.off()

capture.output(raftery.diag(zc_1), heidel.diag(zc_1)) %>% 
  write(., "results/diagnostics_lto.txt")

# 6. Extract predictions from posterior ----------------------------------------

## Species specific differences in occuppancy per tree identity:

zc_2 <- parCodaSamples(cl = cl, model = "lto",
                       variable.names = c("birch_out", "birch_r", 
                                          "pine_out", "pine_r",
                                          "spruce_out", "spruce_r",
                                          "aspen_out", "aspen_r",
                                          "oak_out", "oak_r",
                                          "alder_out", "alder_r",
                                          "tree_out", "tree_r"),
                       n.iter = samples, 
                       thin = n.thin)

## Combine MCMC chains:
zc_2 <- combine.mcmc(zc_2)

## Extract P(occupancy) and species names:
out <- as.data.frame(summary(zc_2)$quantiles[, c("50%", "2.5%", "97.5%")])
out <- round(out, 2)
out$value <- paste0(out$`50%`, "(", out$`2.5%`, "-", out$`97.5%`, ")")
out$tsp <- unlist(tstrsplit(rownames(out), split = "_", keep = 1))
out$identity <- c(levels(lto$species), "Richness")

## Export the table:
write.csv(dcast(out, identity ~ tsp, value.var = "value"), 
          "results/lto_tsp.csv", 
          row.names = FALSE)

end <- Sys.time()
end-start

## DIfferences in Species richness per tree species:

zc_3 <- parCodaSamples(cl = cl, model = "lto",
                       variable.names = c("birch_aspen", "oak_aspen", 
                                          "alder_aspen", "pine_aspen",
                                          "spruce_aspen", "oak_birch", 
                                          "alder_birch", "pine_birch",
                                          "spruce_birch", "alder_oak",
                                          "pine_oak", "spruce_oak",
                                          "pine_alder", "spruce_alder",
                                          "spruce_pine"),
                       n.iter = samples, 
                       thin = n.thin)

## Combine MCMC chains:
zc_3 <- combine.mcmc(zc_3)

## Extract 95% CIs and P(diff >= 0) of the difference:
ANOVA <- as.data.frame(summary(zc_3)$quantiles[, c("50%", "2.5%", "97.5%")])
ANOVA$ecdf <- apply(zc_3, 2, function(x) 1-ecdf(x)(0)) 

write.csv(ANOVA, "results/ANOVA_results_lto.csv")

end <- Sys.time()
end-start

## -------------------------------END-------------------------------------------
