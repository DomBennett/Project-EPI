## 05/12/2013
## D.J. Bennett
## Exploring how to map characters to a tree

## Libraries
require(ape)
require(phytools)
require(geiger)
require(diversitree)

## Generating data
phylo <- sim.bdtree(b = 1, d = 0, n = 100)
x <- sim.char(phy = phylo, par = list(matrix(c(-0.2,0.2,0.2,-0.2),2,2)), model = "discrete")[,,1]
mphylo <- make.simmap(phylo, x)

plotSimmap(mphylo)

?vcv