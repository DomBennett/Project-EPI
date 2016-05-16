
normalise <- function(xs) {
  (xs - min(xs, na.rm=TRUE)) /
    ( max(xs, na.rm=TRUE) - min(xs, na.rm=TRUE) )
}

getSisters <- function(phylo, node) {
  parent <- MoreTreeTools::getParent(phylo, node=node)
  sister.edges <- phylo$edge[ ,1] == parent
  sister.nodes <- phylo$edge[sister.edges, 2]
  sister.nodes <- sister.nodes[sister.nodes != node]
  sister.nodes
}

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

calcMetrics <- function(phylos) {
  p <- phylos[[1]]
  res <- calcMetricsPhylo(p)
  # TODO: how to handle node differences between trees?
  # if(length(phylos) > 1) {
  #   clbool <- !colnames(res) %in% c('node', 'clade_label')
  #   for(p in phylos[2:length(phylos)]) {
  #     tmp <- calcMetricsPhylo(p)
  #     rwtmp <- match(tmp[['node']], res[['node']])
  #     rwtmp <- rwtmp[!is.na(rwtmp)]
  #     rwres <- match(res[['node']], tmp[['node']])
  #     rwres <- rwres[!is.na(rwres)]
  #     res[rwres, clbool] <- (res[rwres, clbool] + tmp[rwtmp, clbool])/2
  #   }
  # }
  res
}

calcMetricsPhylo <- function (phylo) {
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
      pd <- s.edge.length
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
      pd <- sum(phylo$edge.length[temp.edges])
    }
    mean.change <- (s.edge.change + d.edge.change)/(s.edge.length + d.edge.length)
    temp.eds <- phylo$eds[phylo$eds$sp %in% descs,2]
    mean.ed <- mean (temp.eds)
    sd.ed <- sd (temp.eds)
    data.frame(node, clade_label, n, s.edge.length, s.edge.change,
         d.edge.length, d.edge.change, pd, time.split, mean.change, mean.ed,
         sd.ed, stringsAsFactors=FALSE)
  }
  addSisterContrasts <- function (i) {
    sister.node <- getSisters(phylo, res$node[i])
    sister.i <- which(res$node %in% sister.node)
    # add 1 and remove 1 to avoid Inf
    contrast.change <- (res$mean.change[i] + 1)/
      (mean(res$mean.change[sister.i], na.rm=TRUE) + 1)
    contrast.change <- contrast.change - 1
    contrast.n <- res$n[i]/mean(res$n[sister.i], na.rm=TRUE)
    contrast.ed1 <- res$mean.ed[i]/mean(res$mean.ed[sister.i], na.rm=TRUE)
    contrast.ed2 <- (res$pd[i]/res$n[i])/(mean(res$pd[sister.i], na.rm=TRUE)/
                            mean(res$n[sister.i], na.rm=TRUE))
    contrast.pd1 <- res$pd[i] / mean(res$pd[sister.i], na.rm=TRUE)
    contrast.pd2 <- res$pd[i] - mean(res$pd[sister.i], na.rm=TRUE)
    contrast.s.edge.length <- res$s.edge.length[i]/
      mean(res$s.edge.length[sister.i], na.rm=TRUE)
    data.frame (contrast.change, contrast.n, contrast.pd1, contrast.pd2,
                contrast.ed1, contrast.ed2, contrast.s.edge.length, stringsAsFactors=FALSE)
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
