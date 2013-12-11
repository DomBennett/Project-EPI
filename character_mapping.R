## 05/12/2013
## D.J. Bennett
## Exploring how to map characters to a tree

## Libraries
require(ape)
require(phytools)
require(geiger)
require(diversitree)

## Functions
nexusEquate <- function(character.list) {
  lapplyFun <- function(x) {
    for (i in 1:length(x)) {
      if (x[i] == "r") {
        x[i] <- "01"
      } else if (x[i] == "s") {
        x[i] <- "12"
      } else if (x[i] == "t") {
        x[i] <- "123"
      } else if (x[i] == "u") {
        x[i] <- "23"
      } else {
        # stop(paste0("Unknown character: ", x[i]))
      }
    }
    return(x)
  }
  return(lapply(character.list, lapplyFun))
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

## Data
characters <- nexusEquate(read.nexus.data(file.path(input.dir, "SimpleNexus.nex")))
phylo <- sim.bdtree(b = 1, d = 0, n = length(characters) - 1)
phylo <- addOutgroup(phylo)
phylo$tip.label <- sample(names(characters))
plot(phylo)
# NB each character transistion must sum to 0
#q <- list(rbind(c(-1, 1/3, 1/3, 1/3), c(1/3, -1, 1/3, 1/3), c(1/3, 1/3, -1, 1/3), 
#                c(1/3, 1/3, 1/3, -1)))
q <- list(rbind(c(-0.1, 0.1), c(0.01, -0.01)))
x <-sim.char(phylo, q, model="discrete", n=1)[,,1]

## Using parsimony
outgroup <- phylo$tip.label[1]
uphylo <- unroot(phylo)
res <- MPR(x[uphylo$tip.label], uphylo, outgroup)
plot(uphylo, no.margin = TRUE, show.tip.label = FALSE)
nodetexts <- apply(res, 1, function(t) paste0(t, collapse = ","))
nodelabels(text = nodetexts, bg = "white", cex = .75)
tiplabels(x[phylo$tip.label], adj = -2, frame = "none", cex = .75)

## Using ML
# https://stat.ethz.ch/pipermail/r-sig-phylo/2013-May/002724.html
#http://www.springer.com/life+sciences/evolutionary+%26+developmental+biology/book/978-1-4614-1742-2
#x <- unlist(lapply(characters, function(t) t[[1]]))
mphylo <- make.simmap(phylo, x[phylo$tip.label], model = "SYM")
plotSimmap(mphylo)
names(mphylo)
mphylo$mapped.edge
mphylo$logL


## MPR
#res <- list()
#for (i in 1:length(characters[[1]])) {
#  x <- unlist(lapply(characters, function(t) t[[i]]))
#  mpr.res <- MPR(x[uphylo$tip.label], uphylo, outgroup)
#  plot(uphylo, no.margin = TRUE, show.tip.label = FALSE)
#  nodetexts <- apply(mpr.res, 1, function(t) paste0(t, collapse = ","))
#  nodelabels(text = nodetexts, frame = "none")
#  tiplabels(x[phylo$tip.label], adj = -2, frame = "none")
#  res <- c(res, list(mpr.res))
#}