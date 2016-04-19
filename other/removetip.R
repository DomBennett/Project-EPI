removeTip <- function (tree, tip.name) {
  edges.to.drop <- c ()
  # find connecting nodes and edges
  tip.node <- which (tree$tip.label == tip.name)
  tip.edge <- which (tree$edge[ ,2] == tip.node)
  edges.to.drop <- c (edges.to.drop, tip.edge)
  internal.node <- tree$edge[tip.edge, 1]
  # if internal node is root, there is no internal edge
  if (internal.node != length (tree$tip.label) + 1) {
    internal.edge <- which (tree$edge[ ,2] == internal.node)
    edges.to.drop <- c (edges.to.drop, internal.edge)
    # if internal node has more than 1 child edge, must
    #  add this length to the already existing edge
    if (sum (tree$edge[, 1] == internal.node) > 1) {
      corres.edge <- which (tree$edge[, 1] == internal.node)
      corres.edge <- corres.edge[corres.edge != tip.edge]
      corres.node <- tree$edge[corres.edge, 2]
      internal.edge.length <- tree$edge.length[internal.edge]
      tree$edge.length[corres.edge] <-
        tree$edge.length[corres.edge] + internal.edge.length
    }
  }
  # remove edges to drop
  new.edges <- tree$edge[-edges.to.drop, ]
  # re-number nodes
  # remove 1 from all nodes greater than or equal to
  #  the internal node
  new.edges[new.edges[ ,1] >= internal.node,1] <-
    new.edges[new.edges[ ,1] >= internal.node,1] - 1
  new.edges[new.edges[ ,2] >= internal.node,2] <-
    new.edges[new.edges[ ,2] >= internal.node,2] - 1
  # remove 1 from all nodes above tip node
  new.edges[new.edges[ ,1] > tip.node,1] <-
    new.edges[new.edges[ ,1] > tip.node,1] - 1
  new.edges[new.edges[ ,2] > tip.node,2] <-
    new.edges[new.edges[ ,2] > tip.node,2] - 1
  # if there is a tripartition ...
  # ordering has changed as a result of dropping
  trip.bool <- sum (new.edges[ ,1] == internal.node - 2) > 2
  if (trip.bool) {
    
  }
  # replace old with new
  tree$edge <- new.edges
  tree$edge.length <- tree$edge.length[-edges.to.drop]
  tree$tip.label <- tree$tip.label[-tip.node]
  # update Nnode
  tree$Nnode <- tree$Nnode - 1
  # tidy up new tree
  if (!is.null (attr (tree, "order"))) 
    attr(tree, "order") <- NULL
  if (!is.null (tree$node.label)) {
    tree$node.label <- tree$node.label[c (-tip.node,
                                          -internal.node)]
  }
  if (!is.null (tree$node.ages)) {
    tree.node.ages <- tree$node.ages[c (-tip.node,
                                        -internal.node)]
  }
  tree
}