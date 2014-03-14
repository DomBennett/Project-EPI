## No copyright, no warranty
## Dominic John Bennett
## Functions for calculating LFI
## 07/02/2014

## Dependencies
require (ape)
require (geiger)
require (plyr)

## Functions
nearestNodeDistance <- function(phylo, node, display = FALSE) {
  # Return the length of the to the next nearest node
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  node: the node number in phylo
  #
  # Return:
  #  numeric
  # find next node one step in
  inner.node <- phylo$edge[match(node, phylo$edge[,2]),1]
  # find all nodes and edges descending from inner node
  d.edges <- which(phylo$edge[,1] %in% inner.node)
  d.nodes <- phylo$edge[d.edges, 2]
  # remove starting node
  nearest.node <- d.nodes[!d.nodes %in% node][1]
  nearest.node.edge <- d.edges[!d.nodes %in% node][1] #in case of polytomies
  if (display) {
    edge.lties <- ifelse(1:nrow(phylo$edge) %in% nearest.node.edge, 1, 3)
    plot.phylo(phylo, edge.lty = edge.lties, show.tip.label = TRUE)
    nodelabels("node", node)
    nodelabels("innernode", inner.node)
    nodelabels("nextnode", nearest.node)
  }
  return(phylo$edge.length[nearest.node.edge])
}

nodeDescendants <- function(phylo, node, display = FALSE) {
  # Return the descendant species from a node
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  node: the node number in phylo
  #
  # Return:
  #  vector of tip labels
  if (!is.numeric(node)) {
    stop("node is not numeric!")
  }
  if (node > phylo$Nnode + length(phylo$tip.label)) {
    stop("node is greater than the number of nodes in phylo!")
  }
  if (node <= length(phylo$tip.label)) {
    term.nodes <- node
  } else {
    term.nodes <- vector()
    temp.nodes <- node
    while (length(temp.nodes) > 0) {
      connecting.nodes <- phylo$edge[phylo$edge[,1] %in% temp.nodes, 2]
      term.nodes <- c(term.nodes, connecting.nodes[connecting.nodes <= length(phylo$tip.label)])
      temp.nodes <- connecting.nodes[connecting.nodes > length(phylo$tip.label)]
    }
  }
  descendants <- phylo$tip.label[term.nodes]
  if (display) {
    tip.cols <- ifelse(phylo$tip.label %in% descendants, "black", "grey")
    plot.phylo(phylo, tip.color = tip.cols, show.tip.label = TRUE)
    nodelabels("node", node)
  }
  return (descendants)  
}

listClades <- function (phylo) {
  # List all descendants of each node in phylogeny.
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #
  # Returns:
  #  list of list of clades and node number
  nodes <- 1:(phylo$Nnode + length (phylo$tip.label))
  nodes <- nodes[-(length (phylo$tip.label) + 1)]
  clades <- list()
  for (node in nodes) {
    clades <- c (clades, list (tips(phylo, node)))
  }
  return (list (clades, nodes))
}

matchClades <- function (qclades, sclades, exact = FALSE) {
  # Phylogenies may have different numbers of species, as such node numbers will not be
  #  directly comparable. This function matchs nodes between phylogenies using clade names.
  #
  # Args:
  #  qclades: query clades, a list of clades, the clades to be mapped to the subject clades
  #  sclades: subject clades
  #  exact: if true, nodes will be mapped to one another if they have all of the same
  #   descendants. If false, nodes will be mapped if descedants are missing to allow for 
  #   missing taxa (i.e. smallest node with all the same descendants).
  #
  # Returns:
  #  vector of query node numbers
  if (exact) {
    res <- match (qclades, sclades)
    res <- res[!is.na (res)]
    res
  } else {
    partMatch <- function (clade, clades) {
      scores <- unlist (lapply (clades, function (x) sum (clade %in% x)))
      if (sum (scores) == 0) {
        NA
      } else {
        matching.clades <- clades[scores == max (scores)]
        matching.pos <- which(scores == max (scores))
        matching.sizes <- sapply (matching.clades, length)
        # smallest clade with all members is best match
        best.match <- matching.pos[matching.sizes == min (matching.sizes)]
        best.match
      }
    }
    # Before partial matching, some time saving steps:
    # 1. find identical matches
    res <- match (qclades, sclades)
    unresolved.qclades <- qclades[which(is.na (res))]
    # 2. Replace resolved clades with NA
    sclades[res[!is.na (res)]] <- NA
    part.res <- unlist (lapply (unresolved.qclades, partMatch, clades = sclades))
    res[is.na (res)] <- part.res
    res
  }
}

addOutgroup <- function (phylo, outgroup.factor = 100) {
  # Add an outgroup to phylogeny
  #  (https://stat.ethz.ch/pipermail/r-sig-phylo/2012-November/002417.html)
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  outgroup.factor: determines length of branch connecting outgroup to the rest of phylo,
  #   large values produce smaller branches (default 100)
  #
  # Return:
  #  phylogeny (ape class)
  rtt.dist <- mean (diag (vcv.phylo (phylo)))
  dist <- rtt.dist/outgroup.factor
  edge.matrix <- matrix (c (3,2,3,1), 2, 2, byrow = TRUE)
  tip <- list (edge = edge.matrix, tip.label = c ("outgroup", "clade.to.be"),
              edge.length = c (dist, rtt.dist + dist), Nnode = 1)
  class (tip) <- "phylo"
  res <- bind.tree (tip, phylo, where = 2)
  res <- reorder.phylo (res, order = "cladewise")
  return (res)
}

parsimonyReconstruction <- function (chars, phylo, order.numeric = TRUE,
                                     min.n = 3) {
  # Return upper and lower character states for all nodes in phylogeny for a list of
  #  characters using parsimony. Where missing data occur, the phylogeny is reduced by
  #  dropping tips. The ancestral character is determined as the most primitive state
  #  for ordered characters and randomly for unordered characters.
  #
  # Args:
  #  chars: matrix with species names as rows and character states as columns
  #  phylo: phylogeny (ape class)
  #  order.numeric: if true, assumes all numeric characters are ordered
  #  min.n: the minimum number of species with character data
  #
  # Returns:
  #  list of node states and reduced phylogeny for each character (reconstruction.obj)
  #if (!is.ultrametric (phylo)) {
  #  stop("Phylogeny must be ultrametric (i.e. w/o polytomies)")
  #}
  if (!all (rownames (chars) %in% phylo$tip.label)) {
    stop("Chars must be a matrix w/ its rownames matching tips in phylogeny")
  }
  reconstruction.obj <- list ()
  for (i in 1:ncol (chars)) {
    temp.chars <- chars[phylo$tip.label,i]
    names (temp.chars) <- phylo$tip.label
    tips.to.drop <- names (temp.chars)[is.na (temp.chars)]
    temp.chars <- temp.chars[!is.na (temp.chars)]
    if ((length (phylo$tip.label) - length (tips.to.drop)) < min.n) {
      cat ( paste0 ("Dropping [", i, "th] character -- too many missing values\n"))
      next
    } else if (length (tips.to.drop) > 0) {
      reduced.tree <- drop.tip (phylo, tips.to.drop)
    } else {
      reduced.tree <- phylo
    }
    # Convert chars to numeric and use primitive state for outgroup
    non.numeric.chars <- FALSE
    if (order.numeric) {
      temp.chars.names <- names (temp.chars)
      temp.chars <- tryCatch (as.numeric (temp.chars),
                             warning = function (c) temp.chars)
      if(is.numeric (temp.chars)) {
        names (temp.chars) <- temp.chars.names
        reduced.tree <- addOutgroup (reduced.tree)
        temp.chars <- c (temp.chars, min (temp.chars))
        # min is the most primitive state
        names (temp.chars)[length (temp.chars)] <- "outgroup"
      } else {
        non.numeric.chars <- TRUE
      }
    }
    # Else, use random character for outgroup
    if (!order.numeric | non.numeric.chars) {
      reduced.tree <- addOutgroup (reduced.tree)
      temp.chars <- c (temp.chars, sample (temp.chars, 1))
      names (temp.chars)[length (temp.chars)] <- "outgroup"
    }
    if (is.rooted (reduced.tree)) {
      reduced.tree <- unroot (reduced.tree)
    }
    # calculate for nodes using parsimony
    res <- MPR (temp.chars[reduced.tree$tip.label], reduced.tree, "outgroup")
    reduced.tree <- root(reduced.tree, outgroup = "outgroup")
    reduced.tree <- drop.tip(reduced.tree, tip = 1)
    # add tip node states
    res <- rbind (matrix (rep (temp.chars[reduced.tree$tip.label], each = 2), ncol = 2,
                        byrow = TRUE), res)
    # node numbers are same as tree with outgroup, but 1 fewer node
    rownames (res) <- 1:(length (reduced.tree$tip.label) + reduced.tree$Nnode)
    if (!is.null (colnames (chars))) {
      res <- list (node.states = res, reduced.tree = reduced.tree,
                   character.name = names (chars)[i])
    } else {
      res <- list (node.states = res, reduced.tree = reduced.tree,
                   character.name = paste0 ("Character [", i,"]"))
    }
    reconstruction.obj <- c (reconstruction.obj, list (res))
  }
  return (reconstruction.obj)
}

calcEdgeChanges <- function (f.phylo, reconstruction.obj, weight.by.edge = TRUE) {
  # Count the number of changes that have occured along each branch in phylogeny given
  #  a reconstruction object using patristic distance.
  #
  # Args:
  #  f.phylo: full phylogeny (ape class)
  #  reconstruction.obj: a list for multiple characters containing a matrix of upper and lower
  #   estimates of node states and a reduced phylogeny (i.e. return from
  #   parsominyReconstruction)
  #
  # Returns:
  #  phylo with new factor: edge.changes (mean changes per edge)
  # Internal functions
  calcEachPhylo <- function (part.reconstruction.obj) {
    calcEachEdge <- function (f.node, r.node) {
      for (each in f.node) { # for multiple possible matching f.nodes
        connecting.node <- which (r.phylo$edge[ ,2] == r.node)
        edge <- r.phylo$edge[connecting.node, ]
        start <- r.node.states[edge[1], ]
        end <- r.node.states[edge[2], ]
        if (is.numeric (start)) { # for ordered states
          change <- abs (sum (start) - sum (end))/2 
        } else { # for unordered states
          change <- sum (start != end)/2
        }
        if (weight.by.edge) {
          edge.length <- f.phylo$edge.length[which (f.phylo$edge[ ,2] == each)]
          change <- change/edge.length
        } else {
          change <- change/length(f.node)
        }
      }
      change
    }
    r.node.states <- part.reconstruction.obj[['node.states']]
    r.phylo <- part.reconstruction.obj[['reduced.tree']]
    r.clade.node <- listClades (r.phylo)
    r.clades <- r.clade.node[[1]]
    r.nodes <- r.clade.node[[2]]
    # find shared nodes between full and reduced trees
    rf.match <- matchClades (r.clades, f.clades)
    res <- mdply (data.frame (r.node = r.nodes,
                              f.node = rf.match), calcEachEdge)
    names (res)[3] <- "change"
    res
  }
  f.clade.node <- listClades (f.phylo)
  f.clades <- f.clade.node[[1]]
  f.nodes <- f.clade.node[[2]]
  res <- ldply (.data = reconstruction.obj, .fun = calcEachPhylo,
                .progress = create_progress_bar (name = "time"))
  # Calculate mean change for each edge across all chars
  edge.changes <- ddply (.data = res, .variables = .(f.node), .fun = summarize,
                         mean.change = mean (change), n = length (change))
  mean.changes <- edge.changes$mean.change[match (phylo$edge[ ,2], edge.changes$f.node)]
  f.phylo$edge.changes <- mean.changes
  f.phylo
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
  

calcLFIMeasures <- function (phylo) {
  # Calculate a living fossil measures based on branch changes
  #
  # Args:
  #  phylo: phylogeny with edge.changes
  #
  # Return:
  #  data.frame with LFI measures
  if (is.null (phylo$edge.changes)) {
    stop ("Phylo class has no edge.changes.")
  }
  calcEachNode <- function (i) {
    node <- phylo$edge[i,2]
    descendants <- nodeDescendants(phylo, node)
    temp.edges <- extractEdges(phylo, descendants, type = 3)
    clade <- paste(descendants, collapse = "|")
    n <- length(descendants)
    s.edge.length <- phylo$edge.length[i]
    s.edge.change <- phylo$edge.changes[i]
    nnnd <- nearestNodeDistance(phylo, node)
    if (n < 2) {
      d.edge.length <- 0
      d.edge.change <- 0
      time <- phylo$edge.length[i]
    } else {
      lf.clade <- extract.clade(phylo, node)
      rtt <- mean(diag(vcv.phylo(lf.clade)))
      time <- phylo$edge.length[i] + rtt
      # Only sum comparable edges -- i.e. those with change
      temp.edge.changes <- phylo$edge.changes[temp.edges]
      comp.edges <- temp.edges[which (!is.na (temp.edge.changes))]
      d.edge.change <- sum (phylo$edge.changes[comp.edges])
      d.edge.length <- sum (phylo$edge.length[comp.edges])
    }
    data.frame (node, clade, n, nnnd, s.edge.length, s.edge.change, d.edge.length,
                d.edge.change, time)
  }
  res <- mdply (.data = data.frame (i = 1:nrow (phylo$edge)), .fun = calcEachNode)
  res[ ,-1]
}

lfiChecker <- function (time, change, performance, cut) {
  hist(time)
  hist(change)
  hist(performance)
  lfi <- time - (change + performance)/2
  plot(time ~ lfi)
  abline (lm (time ~ lfi), col = "red")
  plot(change ~ lfi)
  abline (lm (change ~ lfi), col = "red")
  plot(performance ~ lfi)
  abline (lm (performance ~ lfi), col = "red")
  hist (lfi)
  cutoff <- quantile (lfi, probs = cut)
  abline (v = cutoff, col = "red")
  lfi
}

plotLFI <- function (phylo, lfi.res) {
  # Plot the results of LFI
  return (NA)
}
