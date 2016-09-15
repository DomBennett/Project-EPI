
library(ape)

addOutgroup <- function(phylo, outgroup.factor = 100) {
  # Add an outgroup to phylogeny
  # (https://stat.ethz.ch/pipermail/r-sig-phylo/2012-November/002417.html)
  #
  # Args:
  #  phylo: phylogeny(ape class)
  #  outgroup.factor: determines length of branch connecting outgroup to the rest of phylo,
  #   large values produce smaller branches(default 100)
  #
  # Return:
  #  phylogeny(ape class)
  rtt.dist <- mean(diag(vcv.phylo(phylo)))
  dist <- rtt.dist/outgroup.factor
  edge.matrix <- matrix(c(3,2,3,1), 2, 2, byrow = TRUE)
  tip <- list(edge = edge.matrix, tip.label = c("outgroup", "clade.to.be"),
              edge.length = c(dist, rtt.dist + dist), Nnode = 1)
  class(tip) <- "phylo"
  res <- bind.tree(tip, phylo, where = 2)
  res <- reorder.phylo(res, order = "cladewise")
  return(res)
}

tree <- rtree(4)
chars <- sample(0:2, 4, replace=TRUE)
chars <- c(0,1,0,2)
chars <- c('a', 'b', 'c', 'd')
reduced.tree <- addOutgroup(tree)
reduced.tree <- unroot(reduced.tree)
chars <- c(chars, 0)
plot(reduced.tree, show.tip.label = F);nodelabels();
tiplabels(text=chars, tip=1:length(chars), frame = "none", adj=-1)
MPR(as.character(chars), reduced.tree, "outgroup")
MPR(chars, reduced.tree, "outgroup")

abs(sum(c(0,0)) - sum(c(0,1)))/2
abs(sum(c(0,0)) - sum(c(1,1)))/2
