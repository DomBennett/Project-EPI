
calcRsqs <- function(chars, parallel=FALSE) {
  # Calculate rsqs for every possible combination of character
  # Return list of combinations + rsq
  .calc <- function(i) {
    c1 <- chars[,cmbs[1,i]]
    c2 <- chars[,cmbs[2,i]]
    pull <- !is.na(c1) & !is.na(c2)
    if(sum(pull) > 20) {
      tst_res <- suppressWarnings(cor.test(c1[pull], c2[pull]))
      return(abs(tst_res$estimate[[1]]))
    }
    NA
  }
  cmbs <- combn(x=colnames(chars), m=2)
  rsqs <- plyr::maply(1:ncol(cmbs), .fun=.calc, .parallel=parallel)
  list('cmbs'=cmbs, 'rsq'=rsqs)
}

matchClades <- function(q_clade_obj, s_clade_obj) {
  # Phylogenies may have different numbers of species, as such node numbers will not be
  #  directly comparable. This function matchs nodes between phylogenies using clade names.
  #
  # Args:
  #  q_clade_obj: result from getClades for query phylogeny
  #  s_clade_obj: subject
  #
  # Returns:
  #  list of subject node numbers in order of query nodes
  .mtch <- function(q) {
    .mtchScr <- function(s) {
      (sum(q %in% s) / length(q)) +
        (sum(s %in% q) / length(s))
    }
    sbst <- which(len_s >= length(q))
    mtch_scrs <- sapply(sclades[sbst], .mtchScr)
    if(all(mtch_scrs < 1)) return(NA)
    snodes[sbst[which.max(mtch_scrs)]]
  }
  sclades <- s_clade_obj[[1]]
  snodes <- s_clade_obj[[2]]
  qclades <- q_clade_obj[[1]]
  len_s <- sapply(sclades, length)
  sapply(qclades, .mtch)
}

getClades <- function(tree) {
  nodes <- 1:(tree$Nnode + length(tree$tip.label))
  clades <- plyr::mlply(.data = data.frame(node = nodes), .fun = getChildren, 
                        tree)
  if (!is.null(tree$all.node.label)) {
    names(clades) <- tree$all.node.label
  }
  sizes <- plyr::ldply(.data = clades, .fun = length)[, 2]
  clades <- clades[order(sizes, decreasing = TRUE)]
  nodes <- nodes[order(sizes, decreasing = TRUE)]
  list(clade.children = clades, clade.node = nodes)
}

getChildren <- function(tree, node) {
  if (!is.numeric(node)) {
    stop("Node is not numeric!")
  }
  if (node > tree$Nnode + length(tree$tip.label)) {
    stop("Node is greater than the number of nodes in tree!")
  }
  if (node <= length(tree$tip.label)) {
    term.nodes <- node
  }
  else {
    term.nodes <- vector()
    temp.nodes <- node
    while (length(temp.nodes) > 0) {
      connecting.nodes <- tree$edge[tree$edge[, 1] %in% 
                                      temp.nodes, 2]
      term.nodes <- c(term.nodes, connecting.nodes[connecting.nodes <= 
                                                     length(tree$tip.label)])
      temp.nodes <- connecting.nodes[connecting.nodes > 
                                       length(tree$tip.label)]
    }
  }
  term.nodes
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
    res <- plyr::mdply(.data = data.frame(i = 1:nrow(r.res)), .fun = eachRNode)
    res
  }
  calcEachRPhylo <- function(part.reconstruction.obj) {
    calcEachRPhyloEdge <- function(r.node) {
      # change is calculated for the PRECEDING edge of a given node
      connecting.edge <- which(r.phylo$edge[ ,2] == r.node)
      nds <- r.phylo$edge[connecting.edge, ]
      start <- r.node.states[nds[1], ]
      end <- r.node.states[nds[2], ]
      score <- (1 + abs(sum(start) - sum(end))/2)/(nstates+1)
      score
    }
    nstates <- max(part.reconstruction.obj[['node.states']])
    r.node.states <- part.reconstruction.obj[['node.states']]
    r.phylo <- part.reconstruction.obj[['reduced.tree']]
    r.clade.node <- MoreTreeTools::getClades(r.phylo)
    r.nodes <- r.clade.node[[2]]
    r.res <- plyr::mdply(data.frame(r.node = r.nodes), calcEachRPhyloEdge)
    # NAs introduced by root
    r.res <- r.res[!is.na(r.res[ ,2]), ]
    names(r.res)[2] <- "change"
    matching.f.nodes <- matchClades(r.clade.node, f.clade.node)
    res <- mapOnFPhylo(r.res, matching.f.nodes, weight.by.edge = TRUE)
    res
  }
  f.clade.node <- MoreTreeTools::getClades(f.phylo)
  f.clades <- f.clade.node[[1]]
  f.nodes <- f.clade.node[[2]]
  res <- plyr::llply(.data = reconstruction.obj, .fun = calcEachRPhylo,
               .parallel=parallel)
  res
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
    # make sure bifurcating
    if (!is.binary.tree (reduced.tree)) {
      reduced.tree <- multi2di (reduced.tree)
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
                   character.name = colnames(chars)[i])
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

changesByClade <- function(nds, changes, tree) {
  # loop through nd indexes and pull out change
  # by character from "changes"
  # Args:
  #  nds: indexes of nodes in tree, derived from clades_phylo
  #  changes: changes by node for each character, derived from calcChange
  #  tree: ape class tree
  # Returns
  #  list of change by character for each node, NA if no info for character
  #   order is same as the character order in changes.
  getChangesByNode <- function(nd) {
    getChange <- function(char) {
      if(nd %in% char[['f.node']]) {
        # take mean, might be more than one possible match in f.phylo
        cs <- char[['change']][char[['f.node']] == nd]
        return(mean(cs, na.rm=TRUE))
      }
      NA
    }
    sapply(changes, getChange)
  }
  getChangesByClade <- function(nd) {
    getChange <- function(i) {
      char_states <- sapply(tmp_nds, function(x) changes_by_node[[x]][[i]])
      if(any(!is.na(char_states))) {
        return(mean(char_states, na.rm=TRUE))
      }
      NA
    }
    edgs <- MoreTreeTools::getEdges(tree, node=nd)
    # change is calculated for the PRECEDING edge of a given node
    tmp_nds <- c(nd, tree$edge[edgs, 2])  # all descending nodes + node itself
    sapply(1:length(changes), getChange)
  }
  nds <- clades_phylo[['clade.node']]
  changes_by_node <- lapply(nds, getChangesByNode)
  res <- lapply(nds, getChangesByClade)
  res
}