## 21/11/2013
## D.J. Bennett
## Using the body mass data to determine ecological niche stasis
## TODO: Use tree distribution to avoid polytomies
## TODO: consider the effect changing the number of species represented in the phylogeny affects
##  ancestral characters estimates
## TODO: ML doesn't work -- why?
## TODO: consider toehr ancestral character options. See ouch and phytools

## Libraries
source(file.path("functions", "EcoDataTools.R"))
source(file.path("functions", "taxaResolve.R"))
require(geiger)
require("reshape2")
precedingNode <- function(phylo, node) {
  return (phylo$edge[phylo$edge[,2] == node, 1])
}
precedingEdge <- function(phylo, node) {
  return (phylo$edge.length[phylo$edge[,2] == node])
}
findCladeName <- function(names) { # this really slows down the code...
  ## takes a taxaResolve object and returns the name of the lowest shared taxonomic group
  ## TODO: not tested
  if (length(names) > 1) {
    split.names <- strsplit(names, " ")
    if (all(lapply(split.names, function(x) x[1]) %in% split.names[[1]][1])) {
      return(split.names[[1]][1])
    } else {
      taxa.resolved <- taxaResolve(names)
      if (any(is.na(taxa.resolved))) {
        return ("Failed taxon resolution")
      } else {
        lineages <- strsplit(as.vector(taxa.resolved$lineage), "\\|")
        previous <- lineages[[1]][2] # start at 2, frist element is empty
        success <- FALSE
        for (i in 3:length(lineages[[1]])) {
          slice <- unlist(lapply(lineages, function(x) x[i]))
          if (all(slice %in% slice[1])) {
            previous <- slice[1]
            next
          } else {
            return (previous)
          }
        }
        return ("findCladeNames failed")
      }
    }
  } else {
    return (names)
  }
}

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
# manageable size
# sample.names <- tree$tip.label[sample(1:3237, 500)] # take uniform sample to avoid polys
# tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% sample.names])
#plot(sampled.tree)
binary.tree <- multi2di(tree)
#plot(binary.tree)
binary.data <- data[match(binary.tree$tip.label, data$Binomial),]

## Ancestral Character Estimation
ace.res <- ace(binary.data$Combined_mass, binary.tree, type="continuous", method="pic")
masses <- c(binary.data$Combined_mass, as.vector(ace.res$ace))

## Plotting
#plotting.tree <- binary.tree
#plotting.tree$tip.label <- binary.data$Combined_mass
#plot(plotting.tree)
#nodelabels(text = signif(ace.res$ace, 2))

## Calculating change
nodes <- c(1:length(binary.tree$tip.label), 2:binary.tree$Nnode + length(binary.tree$tip.label))
# all the terminal nodes and all the internal excluding root
mass.changes <- edge.lengths <- preceding.nodes <- mass.nodes <- mass.preceding.nodes <-
  taxa <- clades <- rep(NA, length(nodes))
for (i in 1:length(nodes)) {
  print (i)
  temp.taxa <- tips(binary.tree, nodes[i])
  clades[i] <- findCladeName(temp.taxa)
  taxa[i] <- paste(tips(binary.tree, nodes[i]), collapse = "|")
  preceding.nodes[i] <- precedingNode(binary.tree, nodes[i])
  edge.lengths[i] <- precedingEdge(binary.tree, nodes[i])
  mass.nodes[i] <- masses[nodes[i]]
  mass.preceding.nodes[i] <- masses[preceding.nodes[i]]
}
mass.change <- mass.preceding.nodes - mass.nodes
res <- data.frame(taxon = taxa, clade = clades, node = nodes, preceding.node = preceding.nodes,
                  mass.node = mass.nodes, mass.preceding.node = mass.preceding.nodes,
                  mass.change = mass.change, edge.length = edge.lengths,
                  mass.by.edge = abs(mass.change/edge.lengths))
res <- res[res$mass.by.edge != Inf,] # remove inf
res <- res[order(res$mass.by.edge),]
write.csv(x = res, file = file.path(output.dir, "LFI_mammal_masses.csv"), row.names = FALSE)