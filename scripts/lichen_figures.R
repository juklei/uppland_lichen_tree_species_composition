## Make figures for ltr and lpsac
##
## First edit: 20190612
## Last edit: 20190612
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("ggplot2")
require("ggpubr")
require("cowplot")
require("rjags")
require("data.table")

## 2. Define or source functions used in this script ---------------------------

theme0 <- function(...) theme(legend.position = "none",
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.ticks = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.ticks.length = unit(0, "null"),
                              panel.border = element_blank(), ...)

## 3. Load and explore data ----------------------------------------------------

load("clean/ltr_pred.rdata")
load("clean/sac_pred_r_dec.rda")
load("clean/sac_pred_r_spruce.rda")
load("clean/sac_pred_r_pine.rda")
load("clean/sac_pred_r_nr_tsp.rda")

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
                       levels =  c("aspen", "birch", "oak", "alder", "pine", 
                                   "spruce"))

g1 <- ggplot(d_ls, aes(x = species, y = richness))
g2 <- geom_point(size = 2)
g3 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)
g4 <- annotate("text", 
               x = c("birch", "oak", "alder", "pine", "spruce"), 
               y = 5, 
               label = c("A", "A,B,C", "B", "B", "C"), 
               size = 10)

png("figures/ltr_tsp.png", 10000/4, 7000/4, "px", res = 600/4)

g1 + g2 + g3 + g4 + theme_classic(40) + ylab("lichen richness per tree") + 
  xlab("")

dev.off()

## 5. Make graphs for lpsac ----------------------------------------------------

## Nr. tree species:

y_nr_tsp <- rbind(summary(export_srntsp[, "r_1tsp"])$quantiles,
                  summary(export_srntsp[, "r_2tsp"])$quantiles,
                  summary(export_srntsp[, "r_3tsp"])$quantiles,
                  summary(export_srntsp[, "r_4tsp"])$quantiles)

d_lntsp <- data.frame("r" = y_nr_tsp[,3],
                      "lower" = y_nr_tsp[,1],
                      "upper" = y_nr_tsp[,5],
                      "nr_tsp" = 1:4)

p1 <- ggplot(d_lntsp, aes(x = nr_tsp, y = r))
p2 <- geom_line(size = 1, linetype = "dashed")
p3 <- geom_point(size = 2)
p4 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
p5 <- annotate("text", x = c(2, 3), y = 20, label = "A", size = 10)

png("figures/r_nr_tsp.png", 9000/4, 7000/4, "px", res = 600/4)

p1 + p2 + p3 + p4 + p5 +
  ylab("expected stand richness") + 
  xlab("number of tree species") +
  theme_classic(40) 

dev.off()

## All percentages combined:

y_all <- rbind(summary(export_srd$r_dec)$quantiles[
               1:length(export_srd$dec_pred), ],
               summary(export_srs$r_spruce)$quantiles[
               1:length(export_srs$spruce_pred), ],
               summary(export_srp$r_pine)$quantiles[
               1:length(export_srp$pine_pred), ])

d_all <- data.frame("r" = y_all[,3],
                    "lower" = y_all[,1],
                    "upper" = y_all[,5],
                    "perc_tree" = c(export_srd$dec_pred,
                                    export_srs$spruce_pred,
                                    export_srp$pine_pred),
                    "tree_species" = c(rep("deciduous", 
                                           length(export_srd$dec_pred)),
                                       rep("spruce", 
                                           length(export_srs$spruce_pred)),
                                       rep("pine", 
                                           length(export_srp$pine_pred))))
## Maxima:
maxima <- data.frame("value" = c(export_srd$dec_max_all,
                                 export_srs$spruce_max_all,
                                 export_srp$pine_max_all),
                     "tree_species" = c(rep("deciduous", length(export_srd$dec_max_all)),
                                        rep("spruce", length(export_srs$spruce_max_all)),
                                        rep("pine", length(export_srp$pine_max_all))))

## Style for both:
s1 <- scale_color_manual(breaks = c("deciduous", "pine", "spruce"), 
                         values = c("#ffc425", "#d11141", "black"))
s2 <- scale_fill_manual(breaks = c("deciduous", "pine", "spruce"),
                        values = c("#ffc425", "#d11141", "black"))

## Predictions:
q1 <- ggplot(d_all, 
             aes(x = perc_tree*100, 
                 y = r, 
                 fill = tree_species, 
                 color = tree_species,
                 lty = tree_species))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, colour = NA)
Q <- q1 + q2 + q3 +
     ylab("expected stand richness") + 
     xlab("percentage of trees") +
     s1 + s2 +
     theme_classic(40) +                  
     theme(legend.position = c(0.29, 0.10), 
           legend.title = element_blank(),
           legend.key.size = unit(3, 'lines'),
           legend.direction = "horizontal")

## Maxima:
p1 <- ggplot(data = maxima_all, 
             aes(x = value*100, 
                 color = tree_species, 
                 fill = tree_species,
                 lty = tree_species))
p2 <- geom_density(alpha = .2)
P <- p1 + p2 + xlim(0, 100) + s1 + s2 + theme0()
  
png("figures/r_all_quad.png", 10000/4, 9000/4, "px", res = 600/4)

plot_grid(P, Q, align = "v", nrow = 2, rel_heights = c(1/5, 4/5))

dev.off()

## -------------------------------END-------------------------------------------
