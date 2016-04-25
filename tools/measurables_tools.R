# In-house functions for running pipeline

calcSuccess <- function(phylo) {
  nodes <- 1:(length(phylo$tip.label) + phylo$Nnode)
  res <- mlply(.data = data.frame(node = nodes),
                .fun = MoreTreeTools::getChildren,
                .progress = create_progress_bar(name = "time"),
               tree=phylo)
  phylo$desc <- res
  phylo
}

calcTime <- function(phylo) {
  # Current: see calcED() in MoreTreeTools
  countDescendants <- function(node) {
    length(MoreTreeTools::getChildren(phylo, node))
  }
  calcSpecies <- function(sp) {
    edges <- MoreTreeTools::getEdges(phylo, tips=as.character(sp), type=2)
    n.descs = mdply(.data = data.frame(node = phylo$edge[edges, 2]),
                     .fun = countDescendants)[ ,2]
    ed <- sum(phylo$edge.length[edges]/n.descs)
    ed
  }
  res <- mdply(.data = data.frame(sp = phylo$tip.label),
                .fun = calcSpecies, .progress = create_progress_bar(name = "time"))
  phylo$eds <- res
  phylo
}

matchClades <- function(q.clade.node, s.clade.node) {
  # Phylogenies may have different numbers of species, as such node numbers will not be
  #  directly comparable. This function matchs nodes between phylogenies using clade names.
  #
  # Args:
  #  q.clade.node: result from getClades for query phylogeny
  #  s.clade.node: subject
  #
  # Returns:
  #  list of subject node numbers in order of query nodes
  s.clade.env <- new.env()
  local(s.clades <- s.clade.node[[1]], env = s.clade.env)
  local(s.nodes <- s.clade.node[[2]],  env = s.clade.env)
  matchQCladeInSClades <- local(function(q.clade) {
    compClades <- function(s.clade) {
      sum(q.clade %in% s.clade) / length(q.clade)
    }
    match.scores <- ldply(.data = s.clades, .fun = compClades)[ ,1]
    match.bool <- match.scores == max(match.scores) & match.scores > 0
    res <- s.nodes[match.bool]
    s.clades <<- s.clades[!match.bool]
    s.nodes <<- s.nodes[!match.bool]
    if(length(res) == 0) {
      NA
    } else {
      res
    }
  }, env = s.clade.env)
  res <- llply(.data = q.clade.node[[1]], .fun = matchQCladeInSClades)
  res
}

calcChange <- function(f.phylo, reconstruction.obj, weight.by.edge = TRUE,
                       parallel=FALSE) {
  # Count the number of changes that have occured along each branch in phylogeny given
  #  a reconstruction object using patristic distance.
  #
  # Args:
  #  f.phylo: full phylogeny(ape class)
  #  reconstruction.obj: a list for multiple characters containing a matrix of upper and lower
  #   estimates of node states and a reduced phylogeny(i.e. return from
  #   parsominyReconstruction)
  #  weight.by.edge: if true, changes that have occurred for multiple possible branches in the full
  #   tree will be mapped on to the full tree in proportion to their branch length. Else, the change
  #   is equally divided between all possible branches. Default true.
  #
  # Returns:
  #  phylo with new factor: edge.changes(mean changes per edge)
  # Internal functions
  mapOnFPhylo <- function(r.res, matching.f.nodes, weight.by.edge) {
    calcEdgeLengthProps <- function(f.nodes) {
      if(length(f.nodes) == 1) {
        return(1)
      }
      res <- rep(NA, length(f.nodes))
      for(i in 1:length(f.nodes)) {
        res[i] <- f.phylo$edge.length[which(f.phylo$edge[ ,2] == f.nodes[i])]
      }
      return(res / sum(res))
    }
    eachRNode <- function(i) {
      f.node <- matching.f.nodes[[i]]
      change <- rep(r.res[['change']][i], length(f.node))
      r.node <- rep(r.res[['r.node']][i], length(f.node))
      if(weight.by.edge) {
        change <- change*calcEdgeLengthProps(f.node)
      } else {
        change <- change/length(f.node)
      }
      data.frame(f.node, r.node, change)
    }
    res <- mdply(.data = data.frame(i = 1:nrow(r.res)), .fun = eachRNode)
    res
  }
  calcEachRPhylo <- function(part.reconstruction.obj) {
    calcEachRPhyloEdge <- function(r.node) {
      connecting.node <- which(r.phylo$edge[ ,2] == r.node)
      edge <- r.phylo$edge[connecting.node, ]
      start <- r.node.states[edge[1], ]
      end <- r.node.states[edge[2], ]
      if(is.numeric(start)) { # for ordered states
        change <- abs(sum(start) - sum(end))/2
      } else { # for unordered states
        change <- sum(start != end)/2
      }
      change
    }
    r.node.states <- part.reconstruction.obj[['node.states']]
    r.phylo <- part.reconstruction.obj[['reduced.tree']]
    r.clade.node <- MoreTreeTools::getClades(r.phylo)
    r.nodes <- r.clade.node[[2]]
    r.res <- mdply(data.frame(r.node = r.nodes), calcEachRPhyloEdge)
    names(r.res)[2] <- "change"
    matching.f.nodes <- matchClades(r.clade.node, f.clade.node)
    res <- mapOnFPhylo(r.res, matching.f.nodes, weight.by.edge = TRUE)
    res
  }
  f.clade.node <- MoreTreeTools::getClades(f.phylo)
  f.clades <- f.clade.node[[1]]
  f.nodes <- f.clade.node[[2]]
  res <- ldply(.data = reconstruction.obj, .fun = calcEachRPhylo,
               .progress = create_progress_bar(name = "time"),
               .parallel=parallel)
  # Calculate mean change for each edge across all chars
  edge.changes <- ddply(.data = res, .variables = .(f.node), .fun = summarize,
                        mean.change = mean(change), n = length(change),
                        .parallel=parallel)
  mean.changes <- edge.changes$mean.change[match(phylo$edge[ ,2], edge.changes$f.node)]
  f.phylo$edge.changes <- mean.changes
  f.phylo
}

parsimonyReconstruction <- function(chars, phylo, order.numeric = TRUE,
                                     min.n = 3) {
  # Return upper and lower character states for all nodes in phylogeny for a list of
  #  characters using parsimony. Where missing data occur, the phylogeny is reduced by
  #  dropping tips. The ancestral character is determined as the most primitive state
  #  for ordered characters and randomly for unordered characters.
  #
  # Args:
  #  chars: matrix with species names as rows and character states as columns
  #  phylo: phylogeny(ape class)
  #  order.numeric: if true, assumes all numeric characters are ordered
  #  min.n: the minimum number of species with character data
  #
  # Returns:
  #  list of node states and reduced phylogeny for each character(reconstruction.obj)
  #if(!is.ultrametric(phylo)) {
  #  stop("Phylogeny must be ultrametric(i.e. w/o polytomies)")
  #}
  if(!all(rownames(chars) %in% phylo$tip.label)) {
    stop("Chars must be a matrix w/ its rownames matching tips in phylogeny")
  }
  reconstruction.obj <- list()
  for(i in 1:ncol(chars)) {
    temp.chars <- chars[phylo$tip.label,i]
    names(temp.chars) <- phylo$tip.label
    tips.to.drop <- names(temp.chars)[is.na(temp.chars)]
    temp.chars <- temp.chars[!is.na(temp.chars)]
    if((length(phylo$tip.label) - length(tips.to.drop)) < min.n) {
      cat( paste0("Dropping [", i, "th] character -- too many missing values\n"))
      next
    } else if(length(tips.to.drop) > 0) {
      reduced.tree <- drop.tip(phylo, tips.to.drop)
    } else {
      reduced.tree <- phylo
    }
    # Convert chars to numeric and use primitive state for outgroup
    non.numeric.chars <- FALSE
    if(order.numeric) {
      temp.chars.names <- names(temp.chars)
      temp.chars <- tryCatch(as.numeric(temp.chars),
                              warning = function(c) temp.chars)
      if(is.numeric(temp.chars)) {
        names(temp.chars) <- temp.chars.names
        reduced.tree <- addOutgroup(reduced.tree)
        temp.chars <- c(temp.chars, min(temp.chars))
        # min is the most primitive state
        names(temp.chars)[length(temp.chars)] <- "outgroup"
      } else {
        non.numeric.chars <- TRUE
      }
    }
    # Else, use random character for outgroup
    if(!order.numeric | non.numeric.chars) {
      reduced.tree <- addOutgroup(reduced.tree)
      temp.chars <- c(temp.chars, sample(temp.chars, 1))
      names(temp.chars)[length(temp.chars)] <- "outgroup"
      uchars <- sort(unique(temp.chars))
      temp.chars.char <- temp.chars
      temp.chars <- sapply(temp.chars, FUN = function(x) match(x, uchars))
    }
    if(is.rooted(reduced.tree)) {
      reduced.tree <- unroot(reduced.tree)
    }
    # calculate for nodes using parsimony
    res <- MPR(temp.chars[reduced.tree$tip.label], reduced.tree, "outgroup")
    if(!order.numeric | non.numeric.chars) {
      res <- matrix(sapply(res, function(x) uchars[x]), ncol = 2)
      colnames(res) <- c('upper', 'lower')
    }
    reduced.tree <- root(reduced.tree, outgroup = "outgroup")
    reduced.tree <- drop.tip(reduced.tree, tip = 1)
    # add tip node states
    res <- rbind(matrix(rep(temp.chars[reduced.tree$tip.label], each = 2), ncol = 2,
                          byrow = TRUE), res)
    # node numbers are same as tree with outgroup, but 1 fewer node
    rownames(res) <- 1:(length(reduced.tree$tip.label) + reduced.tree$Nnode)
    if(!is.null(colnames(chars))) {
      res <- list(node.states = res, reduced.tree = reduced.tree,
                   character.name = names(chars)[i])
    } else {
      res <- list(node.states = res, reduced.tree = reduced.tree,
                   character.name = paste0("Character [", i,"]"))
    }
    reconstruction.obj <- c(reconstruction.obj, list(res))
  }
  return(reconstruction.obj)
}

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