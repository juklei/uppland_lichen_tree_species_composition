## We want to shuffle the order of the trees on each plot n times.
## For each shuffle we virtually put the species seen on the first tree
## on all the following trees for that shuffle round and by plot
## For that I create a function which does n shuffelings per plot and apply it 
## then to a list.
##
## First edit: 20190605
## Last edit: 20190605
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)

## 2. Load and explore data ----------------------------------------------------

l_obs <- as.data.table(read.csv("data/Data_lavar_Almunge15_March_2019.csv"))
f_subplot <- read.csv("data/forest_data_uppland_subplot.csv")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")

## 3. Create list of matrices per plot that goes into the shuffeling ----------- 

## Exclude data rows with unknown tree species:
l_obs <- l_obs[l_obs$Tree.species != "check", ]

## Change names of plot and circles for consistency with forest data:
colnames(l_obs)[2] <- "plot"
l_obs$plot <- as.factor(paste0("plot_", l_obs$plot))
colnames(l_obs)[3] <- "circle_10m"
levels(l_obs$circle_10m) <- c("middle", "east", "west")

## Reduce to trees of known species:
l_obs <- droplevels(l_obs[!is.na(l_obs$Tree.species), ])

## Reduce l_obs to only uncut trees:
l_obs <- l_obs[!is.na(l_obs$Tree.no), c(2:7, 10:14, 17:137)]

## Merge Physcia adscendens, P. tenella och P. adscendens/P. tenella to
## P. adscendens/P. tenella"
T1 <- l_obs[, c("Physcia.adscendens", 
                "Physcia.tenella",
                "P..adscendens.tenella")]
l_obs$P..adscendens.tenella <- ifelse(rowSums(T1, na.rm = TRUE) > 0, 1, NA)

## Reduce data set to desired columns:
lo_red <- l_obs[, c(1, 7, 12:24, 26:105, 107, 109:127, 129:132)]

## Sum species obs to total:
lo_red <- lo_red[, lapply(.SD, sum, na.rm = TRUE), by = c("plot", "Tree.no")]
lo_red[lo_red > 1] <- 1

## Split into a list of data frames per plot:
lo_list <- split(lo_red, by = "plot", keep.by = FALSE)

## Turn the elements into matrices without the tree number:
lo_list <- lapply(lo_list, function(x) as.matrix(x[, -1]))

## Delete all species columns which have never been seen on a plot:
lo_list <- lapply(lo_list, function(x) x[, colSums(x) != 0])

## 4. Create the data set which is needed to fit the accumultation curve:

## How many shuffelings?
S <- 100

## Define the function:
sac_create <- function(x) {
  
  out <- matrix(NA, max(l_obs$Tree.no), S)
  for(i in 1:S) {
    ## Select random rows in x:
    D <- x[sample(nrow(x)), ]
    for(j in 1:ncol(D)) {D[min(which(D[, j] == 1)):nrow(D), j] <- 1}
    out[1:nrow(D), i] <- rowSums(D)
  }
  return(out)

}

## Apply the function:
l_sac <- lapply(lo_list, sac_create)

## Store order of the list elements for ordering forest data:
plot_order <- names(l_sac)

## Turn l_sac into an array:
sad <- array(unlist(l_sac), dim = c(max(lengths(l_sac)/S), S, length(l_sac)))

## 4. Make a data frame with the tree data with the same plot order as the -----
##    list above

lof_tree <- as.data.table(l_obs[, c(1:2, 8:9)])
lt_red <- lof_tree[, list("circle_10m" = unique(circle_10m),
                          "pine" = sum(Tree.species == "Ps")/nrow(.SD),
                          "spruce" = sum(Tree.species == "Pa")/nrow(.SD),
                          "dbh" = mean(Tree.diameter.130.cm.above.ground,
                                       na.rm = TRUE)),
                   by = "plot"]
lt_red$dec <- 1-(lt_red$pine + lt_red$spruce)

## From plot level We want to have the average dbh per plot as a proxy for age:
lof_tree <- merge(lt_red, f_plot[, c(1, 20)], all.x = TRUE, by = "plot")

## Then we want to add from the subplot level, the lidar measurments:
lof_tree <- merge(lof_tree,
                  f_subplot[f_subplot$buffer == 10, c(1:2, 6:8)],
                  all.x = TRUE,
                  by = c("plot", "circle_10m"),
                  allow.cartesian = TRUE)

## Change the order of lof_tree according to l_sac above:
lof_tree <- lof_tree[order(match(plot, plot_order))]

## Check if the order of l_sac and lof_tree are the same:
all(lof_tree$plot == names(l_sac))

## 5. Store the data sets used in the jags models here: ------------------------ 

dir.create("clean")

## Export:
save(sad, file = "clean/species_accumulation_data.rda")
save(lof_tree, file = "clean/sad_forest_part.rda")

## -------------------------------END-------------------------------------------
