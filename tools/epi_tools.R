
addNodeAges <- function (phylo) {
  # Add node.ages to phylo object
  phylo.age <- max (diag (vcv.phylo (phylo)))
  node.ages <-
    mdply (.data = data.frame (node = 1:(length (phylo$tip.label) + phylo$Nnode)),
           .progress = create_progress_bar (name = "time"), .fun = getNodeAge, phylo,
           phylo.age)[ ,2]
  phylo$node.age <- node.ages
  phylo
}

addEdgeLabels <- function (phylo) {
  findLabel <- function (i) {
    node <- phylo$edge[i,1]
    phylo$clade_labels[i]
  }
  phylo$edge.label <- mdply (.data = data.frame (i = 1:nrow (phylo$edge)),
                             .fun = findLabel)[ ,2]
  phylo
}


plotEdgeChanges <- function (phylo, by.char = FALSE) {
  # Plot phylogeny with character labels on tips, edges and nodes (for sanity checking)
  #
  # Args:
  #  phylo: phylogeny (ape class) with edge.changes and edge.change.obj
  #
  # Returns:
  #  None
  if (by.char) {
    for (i in 1:length (phylo$edge.change.obj)) {
      node.states <- phylo$edge.change.obj[[i]][['node.states']]
      reduced.tree <- phylo$edge.change.obj[[i]][['reduced.tree']]
      edge.changes <- phylo$edge.change.obj[[i]][['changes']]
      character.name <- phylo$edge.change.obj[[i]][['character.name']]
      edge.scores <- edge.changes[2, edge.changes[1, ] != 0]
      corres.edges <- edge.changes[1, edge.changes[1, ] != 0]
      plot (reduced.tree, show.tip.label = FALSE, main = character.name)
      nodelabels (paste("[", node.states[-(1:length (reduced.tree$tip.label)), 1], ",",
                       node.states[-(1:length (reduced.tree$tip.label)), 2], "]",
                       sep = ""))
      tiplabels (node.states[1:length(reduced.tree$tip.label),1], adj = -2)
      edgelabels (text = edge.scores, edge = corres.edges, adj = c(0, 1))
      x <- readline (paste0 ("Character [", i,
                            "]. Press return for next character or Esc to exit."))
    }
  } else {
    plot (phylo, main = "Mean number of changes by branch")
    edgelabels (text = round (phylo$edge.changes, 2))
  }
}
  

calcMetrics <- function (phylo) {
  # Calculate time, success, change and other metrics per node
  #
  # Args:
  #  phylo: phylogeny with edge.changes, eds and n.desc
  #
  # Return:
  #  data.frame
  calcEachNode <- function (node) {
    #node <- phylo$edge[i,2]
    edge <- which (phylo$edge[ ,2] == node)
    clade_label <- phylo$clade_labels[node]
    temp.edges <- MoreTreeTools::getEdges(phylo, node=node, type=3)
    descs <- phylo$desc[[node]]
    n <- length(descs)
    #n <- sum (phylo$scope[phylo$tip.label %in% descs])
    s.edge.length <- phylo$edge.length[edge]
    s.edge.change <- phylo$edge.changes[edge]
    if (n < 2) {
      d.edge.length <- 0
      d.edge.change <- 0
      time.split <- phylo$edge.length[edge]
    } else {
      # make it 0 if NA
      if(is.na(s.edge.change)) {
        s.edge.change <- 0
      }
      lf.clade <- extract.clade(phylo, node)
      rtt <- mean(diag(vcv.phylo(lf.clade)))
      time.split <- phylo$edge.length[edge] + rtt
      # Only sum comparable edges -- i.e. those with change
      temp.edge.changes <- phylo$edge.changes[temp.edges]
      comp.edges <- temp.edges[which (!is.na (temp.edge.changes))]
      d.edge.change <- sum (phylo$edge.changes[comp.edges])
      d.edge.length <- sum (phylo$edge.length[comp.edges])
    }
    mean.change <- (s.edge.change + d.edge.change)/(s.edge.length + d.edge.length)
    pd <- s.edge.length + d.edge.length
    temp.eds <- phylo$eds[phylo$eds$sp %in% descs,2]
    mean.ed <- mean (temp.eds)
    sd.ed <- sd (temp.eds)
    parent <- MoreTreeTools::getParent(phylo, node=node)
    if(sum(phylo$edge[ ,1] == parent) > 2) {
      # polytomous parent
      sister.node <- NA
    } else {
      sister.node <- MoreTreeTools::getSister(phylo, node=node)
    }
    data.frame (node, clade_label, sister.node, n, s.edge.length, s.edge.change,
                d.edge.length, d.edge.change, pd, time.split, mean.change, mean.ed,
                sd.ed, stringsAsFactors=FALSE)
  }
  addSisterContrasts <- function (i) {
    sister.node <- res[i,'sister.node']
    if(is.na(sister.node)) {
      contrast.change <- contrast.n <- contrast.ed <-
        contrast.pd <- NA
    } else {
      sister.i <- which (res$node == sister.node)
      # add 1 and remove 1 to avoid Inf
      contrast.change <- (res$mean.change[i] + 1)/
        (res$mean.change[sister.i] + 1)
      contrast.change <- contrast.change - 1
      contrast.n <- res$n[i]/res$n[sister.i]
      contrast.ed <- res$mean.ed[i]/res$mean.ed[sister.i]
      contrast.pd <- res$pd[i]/res$pd[sister.i]
    }
    data.frame (contrast.change, contrast.n, contrast.pd,
                contrast.ed, stringsAsFactors=FALSE)
  }
  # all nodes apart from root
  nodes <- c (1:length (phylo$tip.label),
              (length (phylo$tip.label) + 2): (length (phylo$tip.label) + phylo$Nnode))
  res <- mdply (.data = data.frame (node = nodes, stringsAsFactors=FALSE),
                .fun = calcEachNode, .progress = create_progress_bar (name = 'time'))
  contrast.res <- mdply (.data = data.frame (i = 1:nrow (res), stringsAsFactors=FALSE),
                         .fun = addSisterContrasts)
  cbind (res, contrast.res[ ,-1])
}

EPIChecker <- function (metrics, cut) {
  hist(metrics$time)
  hist(metrics$change)
  hist(metrics$success)
  plot(time ~ epi, data=metrics)
  abline (lm (time ~ epi, data=metrics), col = "red")
  plot(change ~ epi, data=metrics)
  abline (lm (change ~ epi, data=metrics), col = "red")
  plot(success ~ epi, data=metrics)
  abline (lm (success ~ epi, data=metrics), col = "red")
  plot(epi_nc ~ epi, data=metrics)
  abline (lm (epi_nc ~ epi, data=metrics), col = "red")
  hist (metrics$epi, xlab = 'EPI', ylab = NULL, main = NULL,
        col = 'cornflowerblue')
  cutoff <- quantile (metrics$epi, probs = cut, na.rm=TRUE)
  abline (v = cutoff, col = "red", lwd = 2)
  text (labels = '<-- Living fossils', x=cutoff, y=nrow(metrics)/100, cex = 0.8)
  hist (metrics$epi_nc, xlab = 'EPI (No change)', ylab = NULL, main = NULL,
        col = 'cornflowerblue')
  cutoff <- quantile (metrics$epi_nc, probs = cut, na.rm=TRUE)
  abline (v = cutoff, col = "red", lwd = 2)
  text (labels = '<-- Living fossils', x=cutoff, y=nrow(metrics)/100, cex = 0.8)
}
