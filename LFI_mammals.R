## 07/01/2014
## D.J. Bennett
## Developing the LFI for real data

## Libraries
require(ape)
require(geiger)

## Functions
source("readNexusData.R")
precedingEdge <- function(phylo, node) {
  return (phylo$edge.length[phylo$edge[,2] == node])
}
addOutgroup <- function(phylo, outgroup.factor = 100) {
  # https://stat.ethz.ch/pipermail/r-sig-phylo/2012-November/002417.html
  rtt.dist <- mean(diag(vcv.phylo(phylo)))
  dist <- rtt.dist/outgroup.factor
  edge.matrix <- matrix(c(3,2,3,1), 2, 2, byrow = TRUE)
  tip <- list(edge = edge.matrix, tip.label = c("outgroup", "clade.to.be"),
              edge.length = c(dist, rtt.dist + dist), Nnode = 1)
  class(tip)<-"phylo"
  return (bind.tree(tip, phylo, where = 2))
}

## Dirs
input.dir <- "0_data"

## Input
chars <- readNexusData(file.path(input.dir, "morpho_matrix_forR.nex"))
lapply(chars)
chars[[1]]

# for each character i need to know if its ordered