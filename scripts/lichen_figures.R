## Make figures for ltr and lpsac
##
## First edit: 20190612
## Last edit: 20190612
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")

## 2. Define or source functions used in this script ---------------------------

#...

## 3. Load and explore data ----------------------------------------------------

load("clean/ltr_pred.rdata")
load("clean/sac_pred.rdata")

## 4. Make graphs for ltr predictions ------------------------------------------

## Group mean plot for tree species:

## Make data set:
y_ls <- cbind(summary(export_ltr$pine_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_ltr$spruce_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_ltr$aspen_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_ltr$oak_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_ltr$alder_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_ltr$birch_mean, quantile, c(.025,.5,.975))$stat)
      
d_ls <- data.frame("richness" = t(y_ls)[,2],
                   "lower" = t(y_ls)[,1],
                   "upper" = t(y_ls)[,3],
                   "species" = c("pine", "spruce", "aspen", "oak", "alder", 
                                 "birch"))
d_ls$species <- factor(d_ls$species, 
                       levels =  c("aspen", "birch", "alder", "oak", "pine", 
                                   "spruce"))

g1 <- ggplot(d_ls, aes(x = species, y = richness))
g2 <- geom_point(size = 2)
g3 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)

png("figures/ltr_tsp.png", 10000/4, 7000/4, "px", res = 600/4)

g1 + g2 + g3 + theme_classic(40) + ylab("lichen richness per tree")

dev.off()

## 5. Make graphs for lpsac ----------------------------------------------------


## -------------------------------END-------------------------------------------
