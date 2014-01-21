## 19/12/2013
## D.J. Bennett
## Developing the LFI

## Libraries
require(ape)
require(geiger)
require(phytools)

## Functions
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

## Parameters
n <- 20
n.chars <- 100
use.parsimony <- TRUE

## Generating data
## Generating tree
# simulated
#phylo <- stree(n, "left")
#phylo <- compute.brlen(phylo, method="Grafen")
# phytools methods
phylo <- sim.bdtree(b = 1, d = 0, n = n)
# ape methods
#plot(rlineage(0.1, 0)) # Yule process with lambda = 0.1
#plot(rlineage(0.1, 0.05)) # simple birth-death process
#b <- function(t) 1/(1 + exp(0.2*t - 1)) # logistic
#layout(matrix(0:3, 2, byrow = TRUE))
#curve(b, 0, 50, xlab = "Time", ylab = "")
#mu <- 0.07
#phylo <- rbdtree(b, mu)
#ltt.plot(phylo)
#plot(drop.fossil(phylo))
if (use.parsimony) {
  phylo <- addOutgroup(phylo) 
}
plot(phylo)
ltt.plot(phylo)

## Generating characters
# phytools
#q <- list()
#for (i in 1:n.chars) {
#  q <- c(q, list(rbind(c(-0.1, 0.1), c(0.1, -0.1))))
#}
#chars <- sim.char(phylo, q, model="discrete")[,,1]
# ape
chars <- matrix(ncol = n.chars, nrow = length(phylo$tip.label))
rownames(chars) <- phylo$tip.label
for (i in 1:n.chars) {
  chars[ ,i] <- rTraitDisc(phylo, k = 2, rate = 0.1)
}

## Calculate change along branches
if (use.parsimony) {
  phylo <- unroot(phylo)
  ancestral.states <- list()
  for (i in 1:ncol(chars)) {
    res <- MPR(chars[phylo$tip.label,i], phylo, "outgroup")
    # rownames(res) <- as.character(as.numeric(rownames(res)) - 1)
    res <- rbind(matrix(rep(chars[ ,i], each = 2), ncol = 2, byrow = TRUE), res)
    ancestral.states <- c(ancestral.states, list(res))
  }
  ## Calculate number of changes for all clades
  morpho.changes <- list()
  for (i in 1:length(ancestral.states)) {
    character <- ancestral.states[[i]]
    nodes <- 1:(phylo$Nnode + length(phylo$tip.label))
    nodes <- nodes[-(length(phylo$tip.label) + 1)] # don't need root node
    temp.res <- temp.names <- rep(NA, length(nodes))
    for (j in 1:length(nodes)) {
      node <- nodes[j]
      descendents <- tips(phylo, node)
      edges <- which.edge(phylo, descendents)
      tot.changes <- 0
      for (edge in edges) {
        edge <- phylo$edge[edge, ]
        start <- character[edge[1], ]
        end <- character[edge[2], ]
        tot.changes <- (sum(start != end)/2) + tot.changes
      }
      temp.res[j] <- tot.changes
      temp.names[j] <- paste(tips(phylo, node), collapse = "|")
    }
    names(temp.res) <- temp.names
    morpho.changes <- c(morpho.changes, list(temp.res))
  }
  phylo <- drop.tip(root(phylo, "outgroup"), "outgroup")
} else {
  ## TODO: Need to re-work this section to conform to the MPR format above
  ## Current thinking: sum the proportion of changed branch length for each character
  morpho.changes3 <- list()
  for (i in 1:ncol(chars)) {
    # http://blog.phytools.org/2011/06/stochastic-character-mapping-on-tree.html
    # I'm not sure how this works.... some sort of MCMC that uses the probabilities caculated
    # by the ML model....
    temp.chars <- chars[phylo$tip.label, i]
    if (length(unique(temp.chars)) > 1) {
      temp.res <- make.simmap(phylo, temp.chars, model = "SYM")
      # amount of change for each edge
      temp.res <- rowSums(temp.res$mapped.edge) -
        abs(temp.res$mapped.edge[,1] - temp.res$mapped.edge[,2])
      morpho.changes3 <- c(morpho.changes3, list(temp.res))
    }
  }
}

## Calc LFI for each node
nodes <- 1:(phylo$Nnode + length(phylo$tip.label))
nodes <- nodes[-(length(phylo$tip.label) + 1)] # don't need root node
morpho <- time <- n.taxa <- n.changes <- pd <- rep(NA, length(nodes))
ntips <- length(phylo$tip.label)
for (i in 1:length(nodes)) {
  node <- nodes[i]
  taxa <- paste(tips(phylo, node), collapse = "|")
  n.changes[i] <- 0
  for (chars in morpho.changes) {
    n.changes[i] <- chars[names(chars) == taxa] + n.changes[i]
  }
  e1 <- precedingEdge(phylo, node) # alternatively use branching.times
  if (node <= ntips) {
    time[i] <- e1
    pd[i] <- e1
    n.taxa[i] <- 1
  } else {
    lf.clade <- extract.clade(phylo, node)
    e2 <- mean(diag(vcv.phylo(lf.clade)))
    time[i] <- e1 + e2
    pd[i]<- sum(lf.clade$edge.length) + e1
    taxa <- lf.clade$tip.label
    n.taxa[i] <- length(taxa)
  }
  morpho[i] <- n.changes[i]/pd[i]
}
res <- data.frame(nodes, n.taxa, pd, time, n.changes, morpho)
# res <- res[n.taxa > 2,] # bias of few changes for small taxonomic groups (Parsimony error?)

## Plotting
plot(res$pd, res$n.changes, xlab = "Phylogenetic diversity of clade",
     ylab = "Number of changes along branch")
plot(res$n.taxa, res$n.changes, xlab = "Number of taxa", ylab = "Number of changes along a branch")
hist(log(res$time), main = 'Log time', xlab = NULL)
hist(log(res$n.taxa), main = 'Log ntaxa', xlab = NULL)
hist(log(res$morpho + 1), main = 'Log Morphological Change', xlab = NULL)
plot(log(res$n.taxa), log(res$morpho + 1), xlab = "Log ntaxa", ylab = "Log Morphological Change")
abline(h = mean(log(res$morpho + 1)), col = "red", lty = 2)
plot(log(res$time), log(res$morpho + 1), xlab = "Log time", ylab = "Log Morphological Change")
abline(h = mean(log(res$morpho + 1)), col = "red", lty = 2)
plot(res$n.taxa, res$time, xlab= "Number of taxa", ylab = "Age")
# variance among lfi stats
sd(res$time/sum(res$time))
sd(res$morpho/sum(res$morpho)) # morphological change has greatest variance
sd(res$n.taxa/sum(res$n.taxa))
# I get a dispersion