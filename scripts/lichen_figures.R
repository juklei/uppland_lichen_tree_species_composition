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

load("clean/ltr_pred.rda")
load("clean/sac_pred_r_dec.rda")
load("clean/sac_pred_r_spruce.rda")
load("clean/sac_pred_r_pine.rda")

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

## Percentage deciduous:

y_dec <- summary(export_srd$r_dec, quantile, c(.025,.5,.975))$stat

d_ld <- data.frame("r" = t(y_dec)[,2],
                   "lower" = t(y_dec)[,1],
                   "upper" = t(y_dec)[,3],
                   "dec" = export_srd$dec_pred)

p1 <- ggplot(d_ld, aes(x = dec, y = r))
p2 <- geom_line(size = 2)
p3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/r_dec.png", 10000/4, 7000/4, "px", res = 600/4)

p1 + p2 + p3 +
  ylab("expected stand richness") + 
  xlab("percentage decidous trees") +
  theme_classic(40) 

dev.off()

## All percentages combined:

y_all <- cbind(summary(export_srd$r_dec, quantile, c(.025,.5,.975))$stat,
               summary(export_srs$r_spruce, quantile, c(.025,.5,.975))$stat,
               summary(export_srp$r_pine, quantile, c(.025,.5,.975))$stat)

d_all <- data.frame("r" = t(y_all)[,2],
                   "lower" = t(y_all)[,1],
                   "upper" = t(y_all)[,3],
                   "perc_tree" = c(export_srd$dec_pred,
                                   export_srs$spruce_pred,
                                   export_srp$pine_pred),
                   "tree_species" = c(rep("deciduous", 
                                          length(export_srd$dec_pred)),
                                      rep("spruce", 
                                          length(export_srs$spruce_pred)),
                                      rep("pine", 
                                          length(export_srp$pine_pred))))

q1 <- ggplot(d_all, 
             aes(x = perc_tree, y = r, fill = tree_species, color = tree_species))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/r_all.png", 10000/4, 7000/4, "px", res = 600/4)

q1 + q2 + q3 +
  ylab("expected stand richness") + 
  xlab("percentage of trees") +
  scale_color_manual(breaks = c("deciduous", "spruce", "pine"), 
                     values = c("blue", "black", "red")) + 
  scale_fill_manual(breaks = c("deciduous", "spruce", "pine"),
                    values = c("blue", "black", "red")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.15, 0.75), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## -------------------------------END-------------------------------------------
