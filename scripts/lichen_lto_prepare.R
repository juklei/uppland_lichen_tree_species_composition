## In this script lichen data is prepared and merged with forest data. 
##
## First edit: 20190318
## Last edit: 20191016
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(dplyr)
library(corrgram)

## 2. Load and explore data ----------------------------------------------------

dir("data")
from_GT <- read.csv("data/from_GT_observable_lichen.csv")
l_obs <- read.csv("data/Data_lavar_Almunge15_March_2019.csv")

head(l_obs)

## 3. Rearrange lichen data set to merge with forest data for analysis. -------- 

## Exclude data rows with unknown tree species:
l_obs <- l_obs[!(l_obs$Tree.species == "check" | is.na(l_obs$Tree.species)), ] 

## Reduce l_obs to only uncut trees:
l_obs <- l_obs[!is.na(l_obs$Tree.no), ]

## Merge Physcia adscendens, P. tenella och P. adscendens/P. tenella to
## P. adscendens/P. tenella"
T1 <- l_obs[, c("Physcia.adscendens", 
                "Physcia.tenella",
                "P..adscendens.tenella")]
l_obs$P..adscendens.tenella <- ifelse(rowSums(T1, na.rm = TRUE) > 0, 1, NA)

## Reduce data set to desired columns:
lo_red <- as.data.table(l_obs[, c(2, 10:12, 14, 17:29, 31:110, 112, 114:132, 134:137)])

## Make data frame based on occupancy of all species per tree:
l_occ <- melt(lo_red, 
              id.vars = colnames(lo_red)[c(1:5)],
              measure.vars = colnames(lo_red)[c(6:length(lo_red))],
              value.name = "observed",
              variable.name = "species")

## Replace NA with 0 in species list to get presence absence per tree:
l_occ$observed[is.na(l_occ$observed)] <- 0

## Calculate total observed per tree:

l_tot <- l_occ[, list("observed" = ifelse(sum(observed) > 0, 1, 0),
                      "tsp" = unique(Tree.species),                     
                      "DBH" = mean(Tree.diameter.130.cm.above.ground)), 
                      by = c("Plot.no.", "Tree.no", "species")]

## 5. Export clean data sets used in the jags models here: --------------------- 

dir.create("clean")

## Export:
write.csv(l_tot, "clean/lto.csv", row.names = FALSE)

## -------------------------------END-------------------------------------------
