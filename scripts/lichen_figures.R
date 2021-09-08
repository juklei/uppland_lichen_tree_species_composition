## Make figures for ltr and lpsac
##
## First edit: 20201120
## Last edit: 20201120
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require(ggtern)
require(ggplot2)
require(ggtern)
require(cowplot)
require(data.table)
require(scales)

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

dir("clean")

site_pred <- read.csv("clean/lpsac_site_pred.csv")
pred_perc <- read.csv("clean/lpsac_pred_perc.csv")            
max_perc <- read.csv("clean/lpsac_max_perc.csv")             
pred_tsp <- read.csv("clean/lpsac_pred_tsp.csv") 
lto <- read.csv("results/lto_tsp.csv")
lto <- lto[lto$identity == "Richness", - 1]

## 4. Make graphs for lto predictions ------------------------------------------

## Group mean plot for tree species:

## Make data set:
lto <- t(lto)

lto_r <- sapply(strsplit(lto, split = "(", fixed = TRUE), FUN = "[[", 1)
lto_lu <- sapply(strsplit(lto, split = "(", fixed = TRUE), FUN = "[[", 2)
lto_l <- sapply(strsplit(lto_lu, split = "-"), FUN = "[[", 1)
lto_u <- sapply(strsplit(lto_lu, split = "-"), FUN = "[[", 2)
lto_u <- unlist(strsplit(lto_u, ")"))

d_ls <- data.frame("richness" = as.numeric(lto_r), 
                   "lower" = as.numeric(lto_l),
                   "upper" = as.numeric(lto_u),
                   "species" = rownames(lto))
d_ls$species <- factor(d_ls$species, 
                       levels =  c("tree", "aspen", "birch", "oak", "alder", 
                                   "pine", "spruce"))
levels(d_ls$species) <- c("tree", "Aspen", "Birch spp.", "Oak", "Alder", "Pine", "Spruce")

g1 <- ggplot(d_ls[d_ls$species != "tree", ], aes(x = species, y = richness))
g2 <- geom_point(size = 4)
g3 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)
g4 <- geom_hline(yintercept = unlist(d_ls[d_ls$species == "tree", 1:3]),
                 linetype = c("solid", "dashed", "dashed"),
                 alpha = 0.4)
# g5 <- annotate("text", 
#                x = c("Birch spp.", "Oak", "Alder", "Pine", "Spruce"), 
#                y = 4, 
#                label = c("A", "A,B,C", "B", "C,D", "D"), 
#                size = 12)

png("figures/Fig1.png", 10000, 7000, "px", res = 600)
ggplot2:::print.ggplot(g1 + g2 + g3 + g4 + #g5 + 
                       theme_classic(45) + 
                       ylab("tree level alpha diversity") + xlab(""))
dev.off()

## 5. Make graphs for lpsac with nr. of tree species ---------------------------

head(pred_tsp)

levels(pred_tsp$div_metric) <- c("beta diversity", "gamma diversity")

p1 <- ggplot(pred_tsp, aes(x = nr_tsp, y = X50.))
p2 <- geom_line(size = 1, linetype = "dashed", colour = "grey")
p3 <- geom_point(size = 3)
p4 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), width = 0.2)
p5 <- facet_grid(div_metric ~ ., scales = "free_y")

png("figures/Fig2.png", 6000, 8000, "px", res = 600)
ggplot2:::print.ggplot(p1 + p2 + p3 + p4 + p5 +
                       ylab("number/fraction of epiphytic lichen species") + 
                       xlab("number of tree species") +
                       theme_classic(40))
dev.off()

## 6. All tree species group percentages combined ------------------------------

## All percentages combined:

## Chose diversity index here:
# div <- "bdiv"
div <- "gdiv"

head(pred_perc)
head(max_perc)

## Style for both:
s1 <- scale_color_manual(breaks = c("dec", "pine", "spruce"), 
                         values = c("#ffc425", "#d11141", "black"))
s2 <- scale_fill_manual(breaks = c("dec", "pine", "spruce"),
                        values = c("#ffc425", "#d11141", "black"))

## Ylab for ggplot & ggtern:
if(div == "bdiv"){ylab <- "beta diversity"} 
if(div == "gdiv"){ylab <- "gamma diversity"} 

## Predictions:
q1 <- ggplot(droplevels(pred_perc[pred_perc$div_metric == div, ]), 
             aes(x = perc_pred*100, y = X50., 
                 fill = L1, color = L1, lty = L1))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = .2, colour = NA)

Q <- q1 + q2 + q3 +
     ylab(ylab) + xlab("percentage of trees") +
     s1 + s2 +
     theme_classic(40) +       
     theme(legend.position = c(0.4, 0.1), 
           legend.title = element_blank(),
           legend.key.size = unit(3, 'lines'),
           legend.background =  element_blank(),
           legend.direction = "horizontal")
if(div == "adiv"){Q <- Q + coord_trans(y = "log10")}

## Maxima:
if(div == "bdiv"){
  max_perc$b_max[max_perc$b_max < 0 | max_perc$b_max >1] <- NA
  p1 <- ggplot(max_perc, aes(b_max*100, color = L1, fill = L1, lty = L1)) 
}
if(div == "gdiv"){
  max_perc$g_max[max_perc$g_max < 0 | max_perc$g_max >1] <- NA
  p1 <- ggplot(max_perc, aes(g_max*100, color = L1, fill = L1, lty = L1))
}

p2 <- geom_density(alpha = .2)
P <- p1 + p2 + xlim(0, 100) + s1 + s2 + theme0()
Q <- plot_grid(P, Q, align = "v", nrow = 2, rel_heights = c(1/5, 4/5))

png(paste0("figures/lpsac_perc_quad_", div, ".png"), 4000, 4000, "px", res = 300)
ggplot2:::print.ggplot(Q)
dev.off()
  
## Ternary plot:

head(site_pred)

R <- ggtern(droplevels(site_pred[site_pred$div_metric == div, ]), 
            aes(x = dec, y = spruce, z = pine)) +
  geom_mask() +
  geom_point(aes(colour = X50.), size = 10) +
  scale_colour_gradient(low = "yellow", high = "blue", breaks = c(3, seq(5, 60, 10))) + 
  xlab("") + ylab("") + ggtern::zlab("") +
  labs(colour = ylab, size = 10) +
  theme_classic(40) + 
  theme(legend.position = "bottom",
        legend.key.size = unit(3, 'lines'),
        legend.direction = "horizontal") +
  Tarrowlab("% spruce") + 
  Larrowlab("% deciduous") + 
  Rarrowlab("% pine") +
  theme_showarrows()

png(paste0("figures/lpsac_perc_tern_", div, ".png"), 4000, 4000, "px", res = 300)
R
dev.off()

## -------------------------------END-------------------------------------------
