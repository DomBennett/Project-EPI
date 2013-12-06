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

## Using parsimony
outgroup <- phylo$tip.label[1]
uphylo <- unroot(phylo)
res <- list()
for (i in 1:length(characters[[1]])) {
  x <- unlist(lapply(characters, function(t) t[[i]]))
  mpr.res <- MPR(x[uphylo$tip.label], uphylo, outgroup)
  plot(uphylo, no.margin = TRUE, show.tip.label = FALSE)
  nodetexts <- apply(mpr.res, 1, function(t) paste0(t, collapse = ","))
  nodelabels(text = nodetexts, frame = "none")
  tiplabels(x[phylo$tip.label], adj = -2, frame = "none")
  res <- c(res, list(mpr.res))
}

## Using ML

mphylo <- make.simmap(phylo, x)

plotSimmap(mphylo)

?vcv

ace(x, phylo, type = "discrete")

plot(phylo)



plot(rtree(30, rooted = FALSE))

tr <- read.tree(text = "(((i,j)c,(k,l)b)a,(h,g)e,f)d;")
x <- c(1, 3, 0, 6, 5, 2, 4)
names(x) <- letters[6:12]
(o <- MPR(x, tr, "f"))
plot(tr)
nodelabels(paste("[", o[, 1], ",", o[, 2], "]", sep = ""))
tiplabels(x[tr$tip.label], adj = -2)

x <- rpois(30, 1)
tr <- rtree(30, rooted = FALSE)
plot(tr, show.tip.label = FALSE)
o <- MPR(x, tr, "t1")
nodelabels(paste("[", o[, 1], ",", o[, 2], "]", sep = ""))
tiplabels(x[tr$tip.label])#, adj = -2)