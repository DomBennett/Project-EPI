## No copyright, no warranty
## Dominic John Bennett
## Functions for calculating LFI
## 07/02/2014

## Dependencies
require (ape)
require (geiger)
require (plyr)

## Functions
safeFromJSON <- function (url, max.trys = 10) {
  # Wrapper for fromJSON
  trys <- 0
  while (trys < max.trys) {
    json.obj <- try (fromJSON(url)[[1]], silent = TRUE)
    if (class(json.obj) == 'try-error') {
      cat ('---- Connection failed: trying again ----\n')
      trys <- trys + 1
      Sys.sleep (10)
    } else {
      return (json.obj)
    }
  }
  stop ("Failed to connect, server may be down.")
}

findParent <- function (name) {
  # Pull lineage for a fossil record
  extractRankNumber <- function (rec) {
    if (rec$sta == "belongs to") {
      return (rec$rnk)
    }
  }
  extractNames <- function (rec) {
    if (rec$sta == "belongs to") {
      return (rec$nam)
    }
  }
  url <- "http://paleobiodb.org/data1.1/taxa/"
  url <- paste0(url, "list.json?name=", name, "&rel=all_parents")
  json.obj <- safeFromJSON (url)
  if (length (json.obj) > 0) {
    return (rev (ldply (.data = json.obj, .fun = extractNames)[ ,1])) # returned in order
  } else {
    return (c ())
  }
}

findClade <- function (lineages) {
  for (i in length (lineages[[1]]):1) {
    subj <- lineages[[1]][i]
    j <- 2
    while (TRUE) {
      if (j > length (lineages)) {
        success <- TRUE
        break
      }
      query <- lineages[[j]]
      if (subj %in% query) {
        j <- j + 1
      } else {
        success <- FALSE
        break
      }
    }
    if (success) {
      return (subj)
    }
  }
  NA
}

labelNodes <- function (phylo) {
  # Use GNR to label all nodes in a phylogeny
  taxa.res <- taxaResolve (phylo$tip.label)
  nodes <- 1:(length (phylo$tip.label) + phylo$Nnode)
  node.label <- rep (FALSE, length (nodes))
  node.label[1:length (phylo$tip.label)] <-
    unlist(lapply (strsplit(phylo$tip.label, " "), function (x) x[1]))
  for (i in (length (phylo$tip.label) + 1):length (nodes)) {
    descendants <- nodeDescendants(phylo, node = i)
    genus.names <- unlist(lapply (strsplit(descendants, " "), function (x) x[1]))
    if (all (genus.names == genus.names[1])) {
      node.label[i] <- genus.names[1]
    } else {
      lineages <- as.character (taxa.res[taxa.res$search.name %in% descendants, "lineage"])
      lineages <- strsplit (lineages, "\\|")
      lineages <- lineages[!is.na (lineages)]
      if (length (lineages) > 0) {
        node.label[i] <- findClade (lineages)
      }
    }
  }
  phylo$node.label <- node.label
  phylo
}

palaeoPull <- function (clade, limit = 'all') {
  # Pull all records for a clade name PBDB
  url <- paste0 ("http://paleobiodb.org/data1.1/occs/list.json?base_name=",
                 clade,"&limit=", limit)
  json.obj <- safeFromJSON(url)
  extractAgeName <- function (rec) {
    data.frame (name = rec$tna, max.age = rec$eag, min.age = rec$lag)
  }
  downloaded <- ldply (.data = json.obj, .fun = extractAgeName)
  downloaded <- downloaded[!duplicated (downloaded), ]
  downloaded
}

getNodeAge <- function (phylo, node, phylo.age = NA) {
  # Get age of node from root
  term.node <- length (phylo$tip.label) + 1
  if (is.na (phylo.age)) {
    phylo.age <- max (diag (vcv.phylo (phylo)))
  }
  if (term.node == node) {
    return (phylo.age)
  }
  if (node < length (phylo$tip.label) & is.ultrametric (phylo)) {
    return (0)
  }
  edges <- c ()
  while (node != term.node) {
    edges <- c (edges, which (phylo$edge[ ,2] == node))
    node <- phylo$edge[phylo$edge[ ,2] == node, 1]
  }
  return (phylo.age - sum (phylo$edge.length[edges]))
}

addTip <- function (phylo, phylo.edge, tip.name, tip.age, node.age,
                    node.label) {
  # Add tip to phylogeny (modified from multi2di)
  #
  # Args
  #  phylo: phylogeny on which the new tip will be added
  #  phylo.edge: edge of phylogeny where tip will be added
  #  tip.name: name of new tip
  #  tip.age: age of new tip (time from root)
  #  node.age: age of newly created node
  #
  # Return
  #   phylo
  insert <- function (target.vector, source.vector, index) {
    #http://stackoverflow.com/questions/1493969/how-to-insert-elements-into-a-vector
    index <- index - 1
    positions <- c (seq_along (target.vector), index + 0.5)
    res <- c (target.vector, source.vector)
    res[order (positions)]
  }
  node.1 <- phylo$edge[phylo.edge, 1]
  node.2 <- phylo$edge[phylo.edge, 2]
  new.tip.edge.length <- node.age - tip.age
  if (new.tip.edge.length < 0) {
    stop ("Node age must be greater than tip age.")
  }
  new.node.edge.length <- phylo$node.age[node.1] - node.age
  if (new.node.edge.length < 0) {
    stop ("Node age must be greater than incipient node age")
  }
  edges.to.replace <- c(phylo.edge, getEdges (phylo, node.2))
  new.edge.lengths <- c (new.node.edge.length, new.tip.edge.length,
                         phylo$edge.length[phylo.edge] - new.node.edge.length,
                         phylo$edge.length[edges.to.replace][-1])
  target <- node.1 + 1
  phylo$edge[phylo$edge > length (phylo$tip.label)] <- 
    phylo$edge[phylo$edge > length (phylo$tip.label)] + 1
  phylo$Nnode <- phylo$Nnode + 1
  new.tip.labels <- c (phylo$tip.label, tip.name)
  new.tip <- length (new.tip.labels)
  new.node <- new.tip + phylo$Nnode
  new.edges <- matrix (nrow = length (edges.to.replace) + 2, ncol = 2)
  new.edges[1, ] <- c (target, new.node)
  new.edges[2, ] <- c (new.node, new.tip)
  for (i in 1:length (edges.to.replace)) {
    new.edge <- phylo$edge[edges.to.replace[i],]
    new.edge[new.edge == target] <- new.node
    new.edges[i+2, ] <- new.edge
  }
  phylo$edge <- rbind (phylo$edge[-edges.to.replace, ], new.edges)
  phylo$edge.length <- c (phylo$edge.length[-edges.to.replace],
                          new.edge.lengths)
  phylo$tip.label <- new.tip.labels
  if (!is.null(attr(phylo, "order"))) 
    attr(phylo, "order") <- NULL
  phylo$node.label <- insert (phylo$node.label, rep(node.label, 2),
                              c (new.tip, new.node))
  phylo$node.age <- insert (phylo$node.age, c (tip.age, node.age),
                            c (new.tip, new.node))
  phylo <- reorder(phylo)
  newNb <- integer(phylo$Nnode)
  n <- length (phylo$tip.label)
  newNb[1] <- n + 1L
  sndcol <- phylo$edge[, 2] > n
  o <- 1 + rank (phylo$edge[sndcol, 2])
  int.node.i <- (length (phylo$tip.label) + 1):
    (length (phylo$tip.label) + phylo$Nnode)
  int.node.label <- phylo$node.label[int.node.i]
  int.node.age <- phylo$node.age[int.node.i]
  int.node.label <- int.node.label[c(1, o)]
  int.node.age <- int.node.age[c(1, o)]
  phylo$node.label[int.node.i] <- int.node.label
  phylo$node.age[int.node.i] <- int.node.age
  phylo$edge[sndcol, 2] <- newNb[phylo$edge[sndcol, 2] - n] <- n + 
    2:phylo$Nnode
  phylo$edge[, 1] <- newNb[phylo$edge[, 1] - n]
  phylo
}

addFossilsToPhylogeny <- function (phylo, records, ex.age = 'min.age') {
  # Add PBDB fossil records to a phylogeny by name and age matching
  #
  # Args
  #  phylo: phylogeny on which fossils are to be added
  #  records: fossil records from PBDB (from palaeoPull)
  #  ex.age: estimate of extinction, either min.age or max.age
  #
  # Return
  #  phylo
  phylo.env <- new.env ()
  local (findEdge <- function (name, age) {
    # Finding matching node in phylogeny by name, then by age
    # Add to oldest matched named node
    matching.nodes <- which (phylo$node.label == name)
    matching.node.ages <- phylo$node.age[matching.nodes]
    node.1 <- matching.nodes[matching.node.ages == max (matching.node.ages)]
    if (length (node.1) > 1) {
      # if more than one matching named node of same age, choose at random
      node.1 <- sample (node.1, 1)
    }
    if (sum (matching.node.ages > age) > 0) {
      node.2 <- phylo$edge[phylo$edge[ ,1] == node.1, 2]
      if (sum (node.2 %in% matching.nodes) != 1) {
        # if all or no decendants of node.1 match name, choose at random
        node.2 <- sample (node.2, 1)
      } else {
        node.2 <- node.2[node.2 %in% matching.nodes]
      }
    } else {
      while (TRUE) {
        node.age <- phylo$node.age[node.1]
        if (node.age > age) {
          break
        }
        node.2 <- node.1
        node.1 <- phylo$edge[phylo$edge[ ,2] == node.1, 1]
      }
    }
    which (phylo$edge[, 1] == node.1 & phylo$edge[ ,2] == node.2)
  }, env = phylo.env)
  local (matchFossilToPhylogeny <- function (record) {
    name <- as.character(record["name"])
    age <- as.numeric (record[ex.age])
    if (name %in% phylo$tip.label) {
      return (NA)
    }
    genus.name <- strsplit(as.character (name), " ")[[1]][1]
    if (genus.name %in% phylo$node.label) {
      lineage <- genus.name
      edge <- findEdge (genus.name, age)
      return (data.frame (name, age, lineage, edge))
    }
    parent <- findParent(genus.name)
    if (length (parent) == 0) {
      print (paste0 ("No lineage record for [", name, "] ..."))
      return (NA)
    }
    matching.clade <- parent[match (TRUE, parent %in% phylo$node.label)]
    if (is.na (matching.clade)) {
      return (NA)
    }
    lineage <- paste (parent[1:which (parent == matching.clade)], collapse = "|")
    edge <- findEdge (matching.clade, age)
    data.frame (name, age, lineage, edge)
  }, env = phylo.env)
  addFossil <- local (function (i) {
    record <- records[i, ]
    record <- matchFossilToPhylogeny (record)
    #print (record)
    if (!is.na (record[[1]])) {
      phylo.edge <- as.integer (record['edge'])
      tip.name <- as.character (record[['name']])
      tip.age <- as.numeric (record['age'])
      if (tip.age == 0) {
        tip.age <- 0.01
      }
      node.age <- phylo$node.age[phylo$edge[phylo.edge, 1]]
      node.age <- node.age - 0.01
      node.label <- as.character (record[['lineage']])
      phylo <<- addTip (phylo, phylo.edge, tip.name, tip.age, node.age, node.label)
    }
  }, env = phylo.env)
  dups <- records$name[duplicated (records$name)]
  for (dup in dups) {
    max.age <- max (records[records$name %in% dup, 'max.age'])
    min.age <- min (records[records$name %in% dup, 'min.age'])
    record <- data.frame (name = dup, min.age, max.age)
    records <- rbind (records[!records$name %in% dup, ], record)
  }
  phylo.age <- max (phylo$node.age)
  records <- records[records[ex.age] < phylo.age, ]
  records <- records [order (records[ex.age], decreasing = TRUE), ]
  local (records <- records, env = phylo.env)
  local (phylo <- phylo, env = phylo.env)
  m_ply (.data = data.frame (i = 1:nrow (records)),
         .fun = addFossil, .expand = FALSE,
         .progress = create_progress_bar (name = "time"))
  get (x = 'phylo', envir = phylo.env)
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

calcFairProportion <- function (phylo) {
  countDescendants <- function (node) {
    length (nodeDescendants (phylo, node))
  }
  calcSpecies <- function (sp) {
    edges <- extractEdges (phylo, as.character (sp), type = 2)
    n.descs = mdply (.data = data.frame (node = phylo$edge[edges, 2]),
                     .fun = countDescendants)[ ,2]
    sum (phylo$edge.length[edges]/n.descs)
  }
  EDs <- mdply (.data = data.frame (sp = phylo$tip.label),
                .fun = calcSpecies, .progress = create_progress_bar (name = "time"))
  EDs
}

findSisterNode <- function (phylo, node) {
  inner.node <- phylo$edge[match(node, phylo$edge[,2]),1]
  # find all nodes and edges descending from inner node
  d.edges <- which(phylo$edge[,1] %in% inner.node)
  d.nodes <- phylo$edge[d.edges, 2]
  # remove starting node
  nearest.node <- d.nodes[!d.nodes %in% node][1]
  nearest.node
}

getEdges <- function (phylo, node) {
  ## Find all edges from given node to tips
  edges <- c ()
  while (TRUE) {
    bool <- phylo$edge[ ,1] %in% node
    if (sum (bool) > 1) {
      node <- phylo$edge[bool, 2]
      edges <- c (edges, which (bool))
    } else {
      break
    }
  }
  edges
}

getNodes <- function (phylo, node) {
  ## Find all nodes from given node to root
  base.node <- length (phylo$tip.label) + 1
  nodes <- c ()
  while (node != base.node) {
    node <- phylo$edge[phylo$edge[ ,2] == node,1]
    nodes <- c (nodes, node)
  }
  nodes
}

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

listClades <- function (phylo) {
  # List all descendants of each node in phylogeny.
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #
  # Returns:
  #  list of list of clades and node number ordered by clade size
  nodes <- 1:(phylo$Nnode + length (phylo$tip.label))
  nodes <- nodes[-(length (phylo$tip.label) + 1)]
  clades <- list()
  for (node in nodes) {
    clades <- c (clades, list (tips(phylo, node)))
  }
  sizes <- ldply (.data = clades, .fun = length)[ ,1]
  clades <- clades[order (sizes, decreasing = TRUE)]
  nodes <- nodes[order (sizes, decreasing = TRUE)]
  list (clades, nodes)
}
  
matchClades <- function (q.clade.node, s.clade.node) {
  # Phylogenies may have different numbers of species, as such node numbers will not be
  #  directly comparable. This function matchs nodes between phylogenies using clade names.
  #
  # Args:
  #  q.clade.node: result from listClades for query phylogeny
  #  s.clade.node: subject
  #
  # Returns:
  #  list of subject node numbers in order of query nodes
  s.clade.env <- new.env ()
  local (s.clades <- s.clade.node[[1]], env = s.clade.env)
  local (s.nodes <- s.clade.node[[2]],  env = s.clade.env)
  matchQCladeInSClades <- local (function (q.clade) {
    compClades <- function (s.clade) {
      sum (q.clade %in% s.clade) / length (q.clade)
    }
    match.scores <- ldply (.data = s.clades, .fun = compClades)[ ,1]
    match.bool <- match.scores == max(match.scores) & match.scores > 0
    res <- s.nodes[match.bool]
    s.clades <<- s.clades[!match.bool]
    s.nodes <<- s.nodes[!match.bool]
    if (length (res) == 0) {
      NA
    } else {
      res
    }
  }, env = s.clade.env)
  res <- llply (.data = q.clade.node[[1]], .fun = matchQCladeInSClades)
  res
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
      uchars <- sort (unique (temp.chars))
      temp.chars.char <- temp.chars
      temp.chars <- sapply (temp.chars, FUN = function (x) match (x, uchars))
    }
    if (is.rooted (reduced.tree)) {
      reduced.tree <- unroot (reduced.tree)
    }
    # calculate for nodes using parsimony
    res <- MPR (temp.chars[reduced.tree$tip.label], reduced.tree, "outgroup")
    if (!order.numeric | non.numeric.chars) {
      res <- matrix (sapply (res, function (x) uchars[x]), ncol = 2)
      colnames (res) <- c('upper', 'lower')
    }
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
  #  weight.by.edge: if true, changes that have occurred for multiple possible branches in the full
  #   tree will be mapped on to the full tree in proportion to their branch length. Else, the change
  #   is equally divided between all possible branches. Default true.
  #
  # Returns:
  #  phylo with new factor: edge.changes (mean changes per edge)
  # Internal functions
  mapOnFPhylo <- function (r.res, matching.f.nodes, weight.by.edge) {
    calcEdgeLengthProps <- function (f.nodes) {
      if (length (f.nodes) == 1) {
        return (1)
      }
      res <- rep (NA, length (f.nodes))
      for (i in 1:length (f.nodes)) {
        res[i] <- f.phylo$edge.length[which (f.phylo$edge[ ,2] == f.nodes[i])]
      }
      return (res / sum (res))
    }
    eachRNode <- function (i) {
      f.node <- matching.f.nodes[[i]]
      change <- rep (r.res[['change']][i], length (f.node))
      r.node <- rep (r.res[['r.node']][i], length (f.node))
      if (weight.by.edge) {
        change <- change*calcEdgeLengthProps (f.node)
      } else {
        change <- change/length (f.node)
      }
      data.frame (f.node, r.node, change)
    }
    res <- mdply (.data = data.frame (i = 1:nrow (r.res)), .fun = eachRNode)
    res
  }
  calcEachRPhylo <- function (part.reconstruction.obj) {
    calcEachRPhyloEdge <- function (r.node) {
      connecting.node <- which (r.phylo$edge[ ,2] == r.node)
      edge <- r.phylo$edge[connecting.node, ]
      start <- r.node.states[edge[1], ]
      end <- r.node.states[edge[2], ]
      if (is.numeric (start)) { # for ordered states
        change <- abs (sum (start) - sum (end))/2 
      } else { # for unordered states
        change <- sum (start != end)/2
      }
      change
    }
    r.node.states <- part.reconstruction.obj[['node.states']]
    r.phylo <- part.reconstruction.obj[['reduced.tree']]
    r.clade.node <- listClades (r.phylo)
    r.nodes <- r.clade.node[[2]]
    r.res <- mdply (data.frame (r.node = r.nodes), calcEachRPhyloEdge)
    names (r.res)[2] <- "change"
    matching.f.nodes <- matchClades (r.clade.node, f.clade.node)
    res <- mapOnFPhylo (r.res, matching.f.nodes, weight.by.edge = TRUE)
    res
  }
  f.clade.node <- listClades (f.phylo)
  f.clades <- f.clade.node[[1]]
  f.nodes <- f.clade.node[[2]]
  res <- ldply (.data = reconstruction.obj, .fun = calcEachRPhylo,
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
  EDs <- evol.distinct (phylo, type = "fair.proportion")
  calcEachNode <- function (i) {
    node <- phylo$edge[i,2]
    descendants <- nodeDescendants(phylo, node)
    temp.edges <- extractEdges(phylo, descendants, type = 3)
    clade <- paste(descendants, collapse = "|")
    n <- length (descendants)
    #n <- sum (phylo$scope[phylo$tip.label %in% descendants])
    s.edge.length <- phylo$edge.length[i]
    s.edge.change <- phylo$edge.changes[i]
    nnnd <- nearestNodeDistance(phylo, node)
    if (length (descendants) < 2) {
      d.edge.length <- 0
      d.edge.change <- 0
      time.split <- phylo$edge.length[i]
    } else {
      lf.clade <- extract.clade(phylo, node)
      rtt <- mean(diag(vcv.phylo(lf.clade)))
      time.split <- phylo$edge.length[i] + rtt
      # Only sum comparable edges -- i.e. those with change
      temp.edge.changes <- phylo$edge.changes[temp.edges]
      comp.edges <- temp.edges[which (!is.na (temp.edge.changes))]
      d.edge.change <- sum (phylo$edge.changes[comp.edges])
      d.edge.length <- sum (phylo$edge.length[comp.edges])
    }
    mean.change <- (s.edge.change + d.edge.change)/s.edge.length + d.edge.length
    mean.ED <- mean (EDs[descendants])
    sd.ED <- sd (EDs[descendants])
    sister.node <- findSisterNode (phylo, node)
    data.frame (node, clade, sister.node, n, nnnd, s.edge.length, s.edge.change, d.edge.length,
                d.edge.change, time.split, mean.change, mean.ED, sd.ED)
  }
  addSisterContrasts <- function (i) {
    sister.node <- res[i,'sister.node']
    sister.i <- which (res$node == sister.node)
    contrast.change <- res$change[i]/res$change[sister.i]
    contrast.n <- res$n[i]/res$n[sister.i]
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
