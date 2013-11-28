## 18/11/2013
## D.J. Bennett
## Living Fossil Index Development: Calculating evolutionary metrics from mammal tree
## TODO: nearest-neighbour distance doesn't seem to work.

## Libraries
source(file.path("Functions", "EcoDataTools.R"))

## Dirs
input.dir <- "0_data"
output.dir <- "1_emetrics"

## Import
tree <- read.tree(file.path(input.dir, "bininda.txt"))
tree$tip.label <- sub("_", " ", tree$tip.label)

## ED
relative.tree <- tree # make branc lengths relative to total tree size
relative.tree$edge.length <- tree$edge.length*100/sum(tree$edge.length)
group <- n.spp <- s.edge <- nnnd <- d.edges <- rep(NA, nrow(relative.tree$edge))
for (i in 1:nrow(relative.tree$edge)) {
  node <- relative.tree$edge[i,2]
  temp.nnnd <- nearestNodeLength(relative.tree, node) # this is not right for large trees
  descendents <- nodeDescendents(relative.tree, node)
  temp.edges <- extractEdges(relative.tree, descendents, type = 3)
  group[i] <- paste(descendents, collapse = "|")
  n.spp[i] <- length(descendents)
  s.edge[i] <- relative.tree$edge.length[i]
  nnnd[i] <- relative.tree$edge.length[i] + temp.nnnd
  if (n.spp[i] < 2) {
    d.edges[i] <- 0
  } else {
    d.edges[i] <- sum(relative.tree$edge.length[temp.edges])
  }
}
res <- data.frame(group, Nspp = n.spp, Supporting_edge_length = s.edge,
                  Nearest_neighbour_node_distance = nnnd, Descendent_edge_length = d.edges)

## Output
write.csv(x = res, file = file.path(output.dir, "emetrics.csv"))