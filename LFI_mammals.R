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
file <- "mbank_X1766_1-8-2014_922_no_notes.nex"
chars <- readNexusData(file.path(input.dir, file))
phylo <- read.tree(file.path(input.dir, "bininda.txt"))
phylo$tip.label <- sub("_", " ", phylo$tip.label)
# only using names in phylo
chars <- chars[(rownames(chars) %in% phylo$tip.label),]
phylo <- drop.tip(phylo, phylo$tip.label[!phylo$tip.label %in% rownames(chars)])

# for each character i need to know if it's ordered
for (i in ncol(chars)) {
  # extract phylogeny for species with character
  char <- chars[ ,i]
  if (any(is.na(char))) {
    char <- char[!is.na(char)]
    temp.phylo <- drop.tip(phylo, phylo$tip.label[!phylo$tip.label %in% names(char)])
  } else {
    temp.phylo <- phylo
  }
  # assume lowest figure is most primitive state
  as.numeric(unique(char))
  as.numeric(c("a"))
}