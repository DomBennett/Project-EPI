## 21/11/2013
## D.J. Bennett
## Using the body mass data to determine morphological stasis
## TODO: Expand to include clades, not just species.

## Libraries
source(file.path("Functions", "EcoDataTools.R"))
require("reshape2")

## Dirs
input.dir <- "0_data"
output.dir <- "2_mmetrics"

## Data
data <- read.delim(file.path(input.dir, "qbodydata.txt"), na.strings = -999)
data <- na.omit(data)
tree <- read.tree(file.path(input.dir, "bininda.txt"))
tree$tip.label <- sub("_", " ", tree$tip.label)

## Data manip
data$Binomial <- paste(data$Genus, data$Species)
molten <- melt(data, id.vars = c("Continent", "Status", "Order", "Family", "Genus",
                                 "Species", "References", "Binomial"))
data <- dcast(molten, formula = Binomial ~ variable, mean)
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data$Binomial])
data <- data[data$Binomial %in% tree$tip.label, ]

## Remove polytomies
sample.names <- tree$tip.label[split0(1:3237, 10)] # take uniform sample to avoid polys
sampled.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% sample.names])
sampled.data <- data[data$Binomial %in% sample.names,]
plot(sampled.tree)

## Ancestral Character Estimation
ace.res <- ace(sampled.data$Combined_mass, sampled.tree)

## Plotting
plotting.tree <- sampled.tree
plotting.tree$tip.label <- sampled.data$Combined_mass
plot(plotting.tree)
nodelabels(text = signif(ace.res$ace, 2))
signif(ace.res$ace, 2)

ace.res$ace[1] - 12.5